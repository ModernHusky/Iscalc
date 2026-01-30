"""Solver Loop 模块

协调LLM和Iscalc的迭代求解过程。
支持详细的日志调试输出，显示完整调用链和工作流程。
"""

import asyncio
import os
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any, Callable, Awaitable
from enum import Enum

from .logger_config import get_phase_logger, print_separator, format_command_result
from .llm_engine import LLMEngine, LLMResponse
from .command_executor import CommandExecutor, ExecutionResult
from .error_handler import ErrorHandler
from .config import SolverConfig


class EventType(Enum):
    """事件类型"""
    THINKING = "thinking"
    COMMAND = "command"
    RESULT = "result"
    ERROR = "error"
    COMPLETE = "complete"
    LOOP_DETECTED = "loop_detected"
    LOOP_RETRY = "loop_retry"  # 循环重试事件
    MAX_ITERATIONS = "max_iterations"


@dataclass
class SolveEvent:
    """求解事件"""
    type: EventType
    content: str
    latex: Optional[str] = None
    step: int = 0
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class SolveResult:
    """求解结果"""
    success: bool
    final_expression: str
    latex_solution: str
    steps: List[Dict[str, Any]]
    error: Optional[str] = None
    iterations: int = 0


class SolverLoop:
    """求解循环"""
    
    def __init__(
        self,
        llm_engine: LLMEngine,
        executor: CommandExecutor,
        config: Optional[SolverConfig] = None
    ):
        self.llm = llm_engine
        self.executor = executor
        self.config = config or SolverConfig()
        self.error_handler = ErrorHandler(self.config.max_consecutive_errors)
        self.logger = get_phase_logger(__name__)
    
    async def solve(
        self,
        expression: str,
        callback: Optional[Callable[[SolveEvent], Awaitable[None]]] = None,
        conditions: Optional[List[str]] = None,
        user_instruction: Optional[str] = None
    ) -> SolveResult:
        """执行完整的求解循环"""
        import sys
        
        # 启动阶段日志
        print_separator()
        self.logger.solver_start("表达式: '%s'", expression)
        sys.stdout.flush()  # 强制刷新
        
        if conditions:
            self.logger.info("   条件: %s", conditions)
        if user_instruction:
            self.logger.info("   用户指导: %s", user_instruction)
        print_separator("─", 40)
        sys.stdout.flush()  # 强制刷新
        
        # 初始化
        init_result = self.executor.initialize(expression, conditions)
        if not init_result.success:
            self.logger.error(f"❌ Initialization failed: {init_result.error}")
            if callback:
                await callback(SolveEvent(
                    type=EventType.ERROR,
                    content=init_result.error or "初始化失败"
                ))
            return SolveResult(
                success=False,
                final_expression=expression,
                latex_solution="",
                steps=[],
                error=init_result.error
            )
        
        self.logger.info(f"✅ Initialized. Result: {init_result.result}")
        if callback:
            await callback(SolveEvent(
                type=EventType.RESULT,
                content=f"初始表达式: {init_result.result}",
                latex=init_result.latex_result,
                step=0
            ))
        
        steps = []
        iteration = 0
        last_error = None
        expression_history = [init_result.result]
        loop_retry_count = 0  # 循环重试计数
        max_loop_retries = 3  # 最大重试次数
        
        # 动态技能状态
        active_skills: List[str] = []
        loaded_skill_names: List[str] = []
        
        from .skills import get_skill_loader
        
        while iteration < self.config.max_iterations:
            iteration += 1
            current_expr = self.executor.get_current_expression()
            current_state = self.executor.get_current_state_name()
            
            # 迭代阶段日志
            self.logger.iteration("迭代 %d | 状态: %s", iteration, current_state)
            self.logger.info("   当前表达式: %s", current_expr[:100] if current_expr else "None")
            
            # 检查是否完成
            if self.executor.is_finished():
                self.logger.complete("求解完成 (executor 检测到)")
                if callback:
                    await callback(SolveEvent(
                        type=EventType.COMPLETE,
                        content="求解完成！",
                        latex=self.executor.get_current_latex()
                    ))
                break
            
            # 检测循环
            loop_expr = self._detect_loop_expr(expression_history)
            if loop_expr:
                self.logger.warning(f"⚠️ Loop detected: {loop_expr}")
                loop_retry_count += 1
                if loop_retry_count > max_loop_retries:
                    self.logger.error("❌ Max loop retries exceeded. Aborting.")
                    # 超过重试次数，停止求解
                    if callback:
                        await callback(SolveEvent(
                            type=EventType.LOOP_DETECTED,
                            content="检测到循环，已重试3次仍无法解决，停止求解"
                        ))
                    return SolveResult(
                        success=False,
                        final_expression=current_expr or expression,
                        latex_solution=self._format_solution_path(steps),
                        steps=steps,
                        error="检测到表达式循环，重试失败",
                        iterations=iteration
                    )
                
                # 通知LLM循环并要求重新思考
                if callback:
                    await callback(SolveEvent(
                        type=EventType.LOOP_RETRY,
                        content=f"检测到循环 (重试 {loop_retry_count}/{max_loop_retries})，尝试其他方法...",
                        metadata={"retry_count": loop_retry_count, "loop_expr": loop_expr}
                    ))
                
                # 移除最后一个重复的表达式，回退状态
                if len(expression_history) > 1:
                    expression_history.pop()
                    # 同时移除对应的step记录
                    if steps:
                        steps.pop()
                
                # 设置错误信息，让LLM知道需要换方法
                last_error = f"检测到循环！表达式 '{loop_expr}' 重复出现。请使用完全不同的方法或命令来简化表达式，不要重复之前的操作。"
                continue
            
            # 生成命令
            if callback:
                await callback(SolveEvent(
                    type=EventType.THINKING,
                    content="正在思考下一步...",
                    step=iteration
                ))
            
            try:
                # 收集流式响应
                full_response = ""
                
                # LLM 思考阶段日志
                self.logger.llm_thinking("生成命令中... (状态: %s, 已装备技能: %d 个)", 
                                         current_state, len(active_skills))
                
                # 根据配置选择生成模式
                if self.config.use_tool_calling:
                    # 使用 Function Calling 模式：LLM 可以实时调用工具读取技能
                    def on_tool_call(tool_name, result):
                        self.logger.info("   工具 %s 返回: %s", tool_name, result[:100] + "..." if len(result) > 100 else result)
                    
                    full_response = await self.llm.generate_with_tools(
                        current_expr or expression,
                        self.executor.get_history(),
                        last_error,
                        current_state,
                        conditions=conditions,
                        user_instruction=user_instruction,
                        active_skills=active_skills,
                        on_tool_call=on_tool_call
                    )
                    
                    if callback:
                        await callback(SolveEvent(
                            type=EventType.THINKING,
                            content=full_response,
                            step=iteration,
                            metadata={"streaming": False, "tool_calling": True}
                        ))
                else:
                    # 使用流式生成模式
                    async for chunk in self.llm.generate_command(
                        current_expr or expression,
                        self.executor.get_history(),
                        last_error,
                        current_state,
                        conditions=conditions,
                        user_instruction=user_instruction,
                        active_skills=active_skills  # 注入已加载的技能
                    ):
                        full_response += chunk
                        if callback:
                            await callback(SolveEvent(
                                type=EventType.THINKING,
                                content=full_response,
                                step=iteration,
                            metadata={"streaming": True}
                        ))
                
                llm_response = self.llm.parse_response(full_response)
                
                # LLM 响应日志
                self.logger.llm_response("命令: %s", llm_response.command)
                self.logger.info("   思考: %s", llm_response.thinking[:150] + "..." if len(llm_response.thinking) > 150 else llm_response.thinking)
                self.logger.info("   解释: %s", llm_response.explanation)
                
                # ============ 自动技能检测 ============
                # 如果 LLM 在 thinking 中提到某个技能但没有执行 cat 命令加载它，
                # 系统会自动检测并加载该技能
                from .skills import detect_skill_mentions
                
                mentioned_skills = detect_skill_mentions(
                    llm_response.thinking + " " + llm_response.explanation,
                    loaded_skill_names
                )
                
                if mentioned_skills:
                    for skill_info in mentioned_skills:
                        skill_path = skill_info["full_path"]
                        skill_rel_path = skill_info["path"]
                        skill_name = skill_info["name"]
                        
                        # 自动加载该技能
                        try:
                            with open(skill_path, 'r', encoding='utf-8') as f:
                                skill_content = f.read()
                            
                            # 醒目的自动加载日志
                            print("")
                            print_separator("═", 60)
                            self.logger.skill_load("🔮 自动检测并装备技能: %s", skill_name)
                            print_separator("─", 60)
                            self.logger.info("   📁 路径: %s", skill_path)
                            self.logger.info("   📄 大小: %d 字节", len(skill_content))
                            self.logger.info("   💡 原因: LLM 在 thinking 中提到了此技能")
                            print_separator("═", 60)
                            print("")
                            
                            # 添加到上下文
                            skill_dir = os.path.dirname(skill_path)
                            skill_block = f"\n## 技能文件 (自动加载): {skill_rel_path}\n> Directory: {skill_dir}\n\n{skill_content}\n\n---\n"
                            active_skills.append(skill_block)
                            loaded_skill_names.append(skill_path)
                        except Exception as e:
                            self.logger.warning("自动加载技能失败: %s - %s", skill_name, str(e))
                
            except Exception as e:
                self.logger.error_phase("LLM 调用失败: %s", str(e))
                error_info = self.error_handler.handle(e, {"expression": current_expr})
                if callback:
                    await callback(SolveEvent(
                        type=EventType.ERROR,
                        content=f"LLM调用失败: {str(e)}",
                        step=iteration
                    ))
                last_error = str(e)
                continue
            
            # 检查是否完成
            if llm_response.is_final:
                self.logger.complete("LLM 判断已完成: %s", llm_response.explanation)
                if callback:
                    await callback(SolveEvent(
                        type=EventType.COMPLETE,
                        content=f"LLM判断已完成: {llm_response.explanation}",
                        latex=self.executor.get_current_latex()
                    ))
                break
            
            # 处理命令
            command = llm_response.command.strip()
            if not command:
                self.logger.warning("Empty command generated.")
                if callback:
                    await callback(SolveEvent(
                        type=EventType.ERROR,
                        content="LLM未生成有效命令",
                        step=iteration
                    ))
                last_error = "未生成有效命令"
                continue
            
            if callback:
                await callback(SolveEvent(
                    type=EventType.COMMAND,
                    content=command,
                    step=iteration,
                    metadata={"thinking": llm_response.thinking, "explanation": llm_response.explanation}
                ))
            
            # Meta-Command Interception: File System Interaction
            if command.startswith("cat ") or command.startswith("read_skill "):
                # Load File / Skill
                target = command.split(" ", 1)[1].strip()
                self.logger.loading("读取文件: %s", target)
                
                # Resolve path (target is relative to llm_iscalc usually, but let's try to be robust)
                # LLM sees paths like "skills/commands/..." which are relative to llm_iscalc
                # We need to map this to absolute path.
                try:
                    base_dir = os.path.dirname(os.path.abspath(__file__))
                    # Check if absolute or relative
                    if os.path.isabs(target):
                        abs_path = target
                    else:
                        abs_path = os.path.join(base_dir, target)
                    
                    if not os.path.exists(abs_path):
                        # Try relative to CWD just in case
                        if os.path.exists(target):
                            abs_path = os.path.abspath(target)
                        else:
                            raise FileNotFoundError(f"File not found: {target}")
                    
                    # Read content
                    with open(abs_path, 'r', encoding='utf-8') as f:
                        file_content = f.read()

                    # Check if it's a SKILL.md -> Equip it
                    is_skill = target.endswith("SKILL.md") or "SKILL.md" in target
                    
                    if is_skill:
                        # Extract skill name from path or content? 
                        # We just need the content with header.
                        # We use the existing logic to format it nicely with directory info
                        # But we already have the content.
                        
                        skill_dir = os.path.dirname(abs_path)
                        
                        # Add to active_skills if not present
                        # Use a hash or just text check to avoid duplicates?
                        # Simplest: Check if path is in loaded_skill_names (we can store path there now)
                        if abs_path in loaded_skill_names:
                             self.logger.info("   技能已加载过: %s", abs_path)
                             last_error = f"File '{target}' 已经加载过了。"
                        else:
                             # 计算文件大小
                             file_size = len(file_content)
                             
                             # 醒目的技能加载日志
                             print("")
                             print_separator("═", 60)
                             self.logger.skill_load("📚 装备技能: %s", target)
                             print_separator("─", 60)
                             self.logger.info("   📁 完整路径: %s", abs_path)
                             self.logger.info("   📂 所在目录: %s", skill_dir)
                             self.logger.info("   📄 文件大小: %d 字节 (%d 行)", file_size, file_content.count('\n') + 1)
                             print_separator("═", 60)
                             print("")
                             
                             # Format similar to get_skill_instructions
                             skill_block = f"\n## 技能文件: {target}\n> Directory: {skill_dir}\n\n{file_content}\n\n---\n"
                             active_skills.append(skill_block)
                             loaded_skill_names.append(abs_path)
                             
                             if callback:
                                await callback(SolveEvent(
                                    type=EventType.RESULT,
                                    content=f"已加载技能文件: {target}",
                                    step=iteration
                                ))
                             
                             # Success record
                             steps.append({
                                "step": iteration, "command": command, "thinking": llm_response.thinking,
                                "explanation": llm_response.explanation, "success": True, "expr_before": current_expr,
                                "expr_after": current_expr, "latex_after": self.executor.get_current_latex(),
                                "error": None, "changed": False
                             })
                             last_error = None
                        continue
                        
                    else:
                        self.logger.loading("读取资源文件: %s (%d 字节)", target, len(file_content))
                        # Normal file read (Layer 3 resource)
                        # Just return content to history/last_error for immediate consumption
                        result_msg = f"File Content ({target}):\n\n{file_content[:2000]}..." # Truncate if too long?
                        if len(file_content) > 2000:
                            result_msg += f"\n(Truncated, total {len(file_content)} chars)"
                            
                        # Pass via last_error mechanism to show in next prompt
                        last_error = f"Command '{command}' Output:\n{result_msg}"
                        
                        if callback:
                             await callback(SolveEvent(
                                type=EventType.RESULT,
                                content=f"读取文件: {target}",
                                step=iteration
                            ))
                        
                        steps.append({
                            "step": iteration, "command": command, "thinking": llm_response.thinking,
                            "explanation": llm_response.explanation, "success": True, "expr_before": current_expr,
                            "expr_after": current_expr, "latex_after": self.executor.get_current_latex(),
                            "error": result_msg, "changed": False
                        })
                        continue

                except Exception as e:
                    self.logger.error_phase("读取文件失败: %s - %s", target, str(e))
                    last_error = f"Read File Failed: {e}"
                    self.executor.history.append({
                             "command": command, "expr_before": current_expr, "success": False, "error": last_error
                    })
                    continue

                
            if command.startswith("run_script ") or command.startswith("python "):
                self.logger.execution("执行脚本: %s", command)
                # Execute Script
                import subprocess
                cmd_parts = command.split()
                
                # Safety check? We trust the agent for now (local execution).
                try:
                    # Capture stdout/stderr
                    # If command is "run_script foo.py", run "python foo.py"
                    if command.startswith("run_script "):
                        script_path = command.split(" ", 1)[1].strip()
                        run_cmd = ["python", script_path]
                    else:
                        run_cmd = cmd_parts # "python foo.py ..."
                        
                    proc = await asyncio.create_subprocess_exec(
                        *run_cmd,
                        stdout=asyncio.subprocess.PIPE,
                        stderr=asyncio.subprocess.PIPE,
                        cwd="." # Or project root
                    )
                    stdout, stderr = await proc.communicate()
                    
                    output = stdout.decode().strip()
                    err_out = stderr.decode().strip()
                    
                    result_msg = f"Stdout:\n{output}\nStderr:\n{err_out}"
                    self.logger.info(f"Script Result: ReturnCode={proc.returncode}, Output len={len(output)}")
                    
                    self.executor.history.append({
                        "command": command,
                        "expr_before": current_expr,
                        "expr_after": current_expr,
                        "success": proc.returncode == 0,
                        "changed": False,
                        "intermediate_steps": [], # No latex steps
                        "error": result_msg if proc.returncode != 0 else None
                        # "output": result_msg # Maybe add output field to executor history schema later?
                        # For now, put output in "error" if failed, or "expr_after" (which is wrong).
                        # Actually executor history schema is strict.
                        # Let's append to "error" if we want LLM to see it easily in history template,
                        # or we modify HISTORY_TEMPLATE to show "Result/Output".
                        # Current HISTORY_TEMPLATE shows "结果: {expr_after}" or "错误: {error}".
                        # We can put result in expr_after for now even if unchanged? No, that's confusing.
                        # Let's put it in `error` field with a prefix "Output:" so it shows up?
                        # Or better: Just rely on it being the last action.
                    })
                    
                    # Update: Modify history template logic in llm_engine to show "Output" field if present.
                    # For now, append to steps.
                    
                    if callback:
                         await callback(SolveEvent(
                            type=EventType.RESULT,
                            content=f"脚本输出:\n{result_msg}",
                            step=iteration
                        ))
                    
                    # We inject the output into the history via a hack: put it in error field but mark success=True?
                    # No, let's put it in a new way.
                    # Actually, for the LLM next turn, we want it to see the output.
                    # We can append the output to the prompt manually next time?
                    # Or modify `executor.history` to accept `output` field and update `llm_engine._build_messages`.
                    # Given the constraints, let's put it in `error` field but prefix with "Output:"
                    # and pass success=True? No, existing logic might ignore error if success=True.
                    
                    # Hack: Treat it as an error (so it shows in history) but with success=False logical flow?
                    # No, loop continues.
                    # Let's put Output in `error` field but call it "Execution Result"
                    
                    # Better: Just put it in active_skills as a temporary "Last Output" block?
                    # Or just rely on `last_error` which is fed into prompt?
                    # The prompt uses `last_error` for "Last command failed...".
                    
                    # Let's use `last_error` variable to pass output to next turn, 
                    # but format it as "Last command output: ...".
                    # And set successful=True in history so it doesn't look like a failure in list.
                    
                    last_error = f"Command '{command}' Output:\n{result_msg}"
                    
                    steps.append({
                        "step": iteration, "command": command, "thinking": llm_response.thinking,
                        "explanation": llm_response.explanation, "success": proc.returncode == 0,
                        "expr_before": current_expr, "expr_after": current_expr, "latex_after": self.executor.get_current_latex(),
                        "error": result_msg, "changed": False
                    })
                    continue
                    
                except Exception as e:
                    self.logger.error_phase("脚本执行失败: %s", str(e))
                    last_error = f"Script Execution Failed: {e}"
                     # ... record failure
                    continue

            # 标准 iscalc 命令执行
            self.logger.execution("执行命令: %s", command)
            exec_result = self.executor.execute(command)
            
            step_record = {
                "step": iteration,
                "command": command,
                "thinking": llm_response.thinking,
                "explanation": llm_response.explanation,
                "success": exec_result.success,
                "expr_before": current_expr,
                "expr_after": exec_result.result if exec_result.success else None,
                "latex_after": exec_result.latex_result if exec_result.success else None,
                "error": exec_result.error if not exec_result.success else None,
                "changed": exec_result.changed if exec_result.success else False
            }
            steps.append(step_record)
            
            if exec_result.success:
                self.logger.result("成功 | 变更=%s", exec_result.changed)
                if exec_result.changed:
                    self.logger.info("   新表达式: %s", exec_result.result[:100] if exec_result.result else "None")
                self.error_handler.reset_consecutive_errors()
                last_error = None
                expression_history.append(exec_result.result)
                
                if callback:
                    await callback(SolveEvent(
                        type=EventType.RESULT,
                        content=exec_result.result or "",
                        latex=exec_result.latex_result,
                        step=iteration,
                        metadata={
                            "changed": exec_result.changed,
                            "intermediate_steps": exec_result.intermediate_steps
                        }
                    ))
            else:
                self.logger.error_phase("执行失败: %s", exec_result.error)
                error_info = self.error_handler.handle(
                    Exception(exec_result.error or "未知错误"),
                    {"command": command, "expression": current_expr}
                )
                last_error = self.error_handler.format_for_llm(error_info)
                
                if callback:
                    await callback(SolveEvent(
                        type=EventType.ERROR,
                        content=exec_result.error or "执行失败",
                        step=iteration
                    ))
                
                if self.error_handler.should_suggest_check_input():
                    self.logger.warning("Suggesting input check due to consecutive errors.")
                    if callback:
                        await callback(SolveEvent(
                            type=EventType.ERROR,
                            content="连续错误过多，建议检查输入表达式",
                            step=iteration
                        ))
        
        # 检查是否达到最大迭代次数
        if iteration >= self.config.max_iterations:
            self.logger.warning("Reached max iterations.")
            if callback:
                await callback(SolveEvent(
                    type=EventType.MAX_ITERATIONS,
                    content=f"达到最大迭代次数 ({self.config.max_iterations})"
                ))
        
        final_expr = self.executor.get_current_expression() or expression
        print_separator()
        self.logger.complete("最终结果: %s", final_expr[:100] if final_expr else "None")
        self.logger.info("   总迭代次数: %d", iteration)
        print_separator()
        return SolveResult(
            success=self.executor.is_finished() or iteration < self.config.max_iterations,
            final_expression=final_expr,
            latex_solution=self._format_solution_path(steps),
            steps=steps,
            iterations=iteration
        )
    
    def _detect_loop(self, expression_history: List[str]) -> bool:
        """检测是否陷入循环"""
        if len(expression_history) < self.config.loop_detection_window:
            return False
        
        recent = expression_history[-self.config.loop_detection_window:]
        # 检查是否有重复表达式
        for expr in recent:
            if recent.count(expr) >= 2:
                return True
        return False
    
    def _detect_loop_expr(self, expression_history: List[str]) -> Optional[str]:
        """检测循环并返回重复的表达式，如果没有循环返回None"""
        if len(expression_history) < self.config.loop_detection_window:
            return None
        
        recent = expression_history[-self.config.loop_detection_window:]
        # 检查是否有重复表达式
        for expr in recent:
            if recent.count(expr) >= 2:
                return expr
        return None
    
    def _format_solution_path(self, steps: List[Dict[str, Any]]) -> str:
        """格式化完整求解路径为LaTeX"""
        if not steps:
            return ""
        
        lines = ["\\begin{align*}"]
        
        for step in steps:
            if step["success"] and step.get("latex_after"):
                cmd = step["command"].replace("_", "\\_")
                lines.append(f"&\\text{{{cmd}}} \\\\")
                lines.append(f"&= {step['latex_after']} \\\\")
        
        lines.append("\\end{align*}")
        return "\n".join(lines)
