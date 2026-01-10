"""Solver Loop 模块

协调LLM和Iscalc的迭代求解过程。
"""

import asyncio
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any, Callable, Awaitable
from enum import Enum

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
    
    async def solve(
        self,
        expression: str,
        callback: Optional[Callable[[SolveEvent], Awaitable[None]]] = None,
        conditions: Optional[List[str]] = None
    ) -> SolveResult:
        """执行完整的求解循环"""
        
        # 初始化
        init_result = self.executor.initialize(expression, conditions)
        if not init_result.success:
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
        
        while iteration < self.config.max_iterations:
            iteration += 1
            current_expr = self.executor.get_current_expression()
            
            # 检查是否完成
            if self.executor.is_finished():
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
                loop_retry_count += 1
                if loop_retry_count > max_loop_retries:
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
                async for chunk in self.llm.generate_command(
                    current_expr or expression,
                    self.executor.get_history(),
                    last_error
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
                
            except Exception as e:
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
                if callback:
                    await callback(SolveEvent(
                        type=EventType.COMPLETE,
                        content=f"LLM判断已完成: {llm_response.explanation}",
                        latex=self.executor.get_current_latex()
                    ))
                break
            
            # 执行命令
            command = llm_response.command.strip()
            if not command:
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
                self.error_handler.reset_consecutive_errors()
                last_error = None
                expression_history.append(exec_result.result)
                
                if callback:
                    await callback(SolveEvent(
                        type=EventType.RESULT,
                        content=exec_result.result or "",
                        latex=exec_result.latex_result,
                        step=iteration,
                        metadata={"changed": exec_result.changed}
                    ))
            else:
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
                    if callback:
                        await callback(SolveEvent(
                            type=EventType.ERROR,
                            content="连续错误过多，建议检查输入表达式",
                            step=iteration
                        ))
        
        # 检查是否达到最大迭代次数
        if iteration >= self.config.max_iterations:
            if callback:
                await callback(SolveEvent(
                    type=EventType.MAX_ITERATIONS,
                    content=f"达到最大迭代次数 ({self.config.max_iterations})"
                ))
        
        final_expr = self.executor.get_current_expression() or expression
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
