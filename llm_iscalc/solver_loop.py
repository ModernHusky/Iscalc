"""Solver Loop 模块

协调LLM和Iscalc的迭代求解过程。
支持详细的日志调试输出，显示完整调用链和工作流程。
"""

import asyncio
import os
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any, Callable, Awaitable, Set
from enum import Enum

from .logger_config import get_phase_logger, print_separator, format_command_result
from .llm_engine import LLMEngine, LLMResponse
from .command_executor import CommandExecutor, ExecutionResult
from .error_handler import ErrorHandler
from .config import SolverConfig
from .models import TokenUsage, PromptComponent, RoundLog


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
    SKILL_LOAD = "skill_load"  # 技能加载事件（Search-o1 风格）
    COMMAND_SUCCESS = "command_success"  # 命令执行成功（2026-02-01 架构重构）
    COMMAND_FAILURE = "command_failure"  # 命令执行失败（2026-02-01 架构重构）
    PLAN_GENERATED = "plan_generated"    # 计划生成
    PLAN_STEP_START = "plan_step_start"  # 计划步骤开始
    PLAN_STEP_COMPLETE = "plan_step_complete" # 计划步骤完成
    PLAN_UPDATED = "plan_updated"        # 计划被动态修改

    SUB_PLAN_GENERATED = "sub_plan_generated" # 子计划生成或更新
    SUB_PLAN_STEP_START = "sub_plan_step_start" # 子计划步骤开始
    SUB_PLAN_STEP_COMPLETE = "sub_plan_step_complete" # 子计划步骤完成


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
    total_token_usage: TokenUsage = field(default_factory=TokenUsage)
    round_logs: List[RoundLog] = field(default_factory=list)


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

        # Search-o1 风格：技能加载次数限制
        self.max_skill_loads_per_solve = 10  # 单次求解最多加载技能数

        # 记忆机制：累积每轮迭代的经验教训
        self.memory = ""  # 记忆内容，会在每轮迭代后更新

    def _detect_and_load_skills(
        self,
        text: str,
        active_skills: List[str],
        loaded_skill_names: Set[str],
        callback = None
    ) -> tuple:
        """
        Search-o1 风格：检测文本中的技能加载标记并加载技能

        借鉴 Search-o1 的核心算法：
        1. 检测 <|load_skill|>技能名<|end_load_skill|> 标记
        2. 加载技能内容
        3. 将标记替换为加载确认信息

        Args:
            text: 待检测的文本（LLM 的 thinking 输出）
            active_skills: 当前活跃的技能内容列表
            loaded_skill_names: 已加载的技能名列表
            callback: 事件回调函数（可选）

        Returns:
            tuple: (处理后的文本, 是否检测到标记, 新加载的技能名列表)
        """
        import re
        from .prompts import SKILL_LOAD_BEGIN, SKILL_LOAD_END
        from .skills import get_skill_loader

        # 🔧 关键修复：规范化 LLM 可能输出的错误格式标记
        # LLM 有时会输出 <load_skill>...</load_skill> 而不是正确的 <|load_skill|>...<|end_load_skill|>
        # 在检测前先转换为正确格式
        text = re.sub(r'<load_skill>([^<]+)</load_skill>',
                      f'{SKILL_LOAD_BEGIN}\\1{SKILL_LOAD_END}', text)
        text = re.sub(r'<<\s*([a-zA-Z_-]+)\s*>>',
                      f'{SKILL_LOAD_BEGIN}\\1{SKILL_LOAD_END}', text)

        # 构建正则匹配模式
        pattern = re.escape(SKILL_LOAD_BEGIN) + r'(.+?)' + re.escape(SKILL_LOAD_END)
        matches = list(re.finditer(pattern, text, re.DOTALL))

        if not matches:
            return text, False, []

        newly_loaded: List[str] = []
        result_text = text

        loader = get_skill_loader()
        try:
            path_to_name = {s.path: s.name for s in loader.discover_skills() if s.path and s.name}
        except Exception:
            path_to_name = {}

        for match in matches:
            requested = match.group(1).strip()

            # 检查加载次数限制
            if len(loaded_skill_names) >= self.max_skill_loads_per_solve:
                result_text = result_text.replace(
                    match.group(0),
                    f"[⚠ 已达到单次求解技能加载上限 ({self.max_skill_loads_per_solve})]"
                )
                continue

            # 定位技能文件（支持 skill 名 / 路径 / 关键词查询）
            skill_path = self._resolve_skill_path(requested)
            if skill_path is None:
                result_text = result_text.replace(
                    match.group(0),
                    f"[✗ 未找到技能: {requested}]"
                )
                continue

            canonical_name = path_to_name.get(skill_path) or requested

            # 检查是否已加载（以 canonical name 为准）
            if canonical_name in loaded_skill_names:
                result_text = result_text.replace(
                    match.group(0),
                    f"[✓ 技能 {canonical_name} 已在上下文中]"
                )
                continue

            # 加载技能内容（优先走 SkillLoader 以跳过 frontmatter + 复用缓存）
            try:
                skill_content_obj = loader.load_skill_content(canonical_name)
                if skill_content_obj:
                    skill_body = skill_content_obj.full_content
                    skill_dir = os.path.dirname(skill_content_obj.metadata.path)
                else:
                    # Fallback: raw read (should be rare)
                    with open(skill_path, 'r', encoding='utf-8-sig') as f:
                        skill_body = f.read()
                    skill_dir = os.path.dirname(skill_path)

                skill_block = (
                    f"\n## 技能文件 (主动加载): {canonical_name}\n"
                    f"> Path: {skill_path}\n"
                    f"> Directory: {skill_dir}\n\n"
                    f"{skill_body}\n\n---\n"
                )
                active_skills.append(skill_block)
                loaded_skill_names.add(canonical_name)
                newly_loaded.append(canonical_name)

                # 替换标记为加载确认
                result_text = result_text.replace(
                    match.group(0),
                    f"[✓ 已加载技能: {canonical_name}]"
                )

                # 打印醒目的日志
                print("")
                print_separator("═", 60)
                self.logger.skill_load("🔖 主动加载技能 (Search-o1 风格): %s", canonical_name)
                print_separator("─", 60)
                self.logger.info("   📁 路径: %s", skill_path)
                self.logger.info("   📄 大小: %d 字节", len(skill_body))
                print_separator("═", 60)
                print("")

            except Exception as e:
                self.logger.warning("技能加载失败: %s - %s", canonical_name, str(e))
                result_text = result_text.replace(
                    match.group(0),
                    f"[✗ 加载失败: {canonical_name} - {str(e)}]"
                )

        return result_text, len(newly_loaded) > 0, newly_loaded

    def _resolve_skill_path(self, skill_name: str) -> Optional[str]:
        """解析技能名称到实际文件路径。

        支持多种格式：
        - 短名称: "rewrite", "strategy-integral"
        - 完整路径: "skills/commands/rewrite/SKILL.md"
        - 关键词查询（Search-o1 风格）:
          - "search: partial fraction decomposition"
          - "partial fraction decomposition" (系统将做模糊匹配/关键词匹配)
        """
        base_dir = os.path.dirname(os.path.abspath(__file__))

        raw = (skill_name or "").strip()
        if not raw:
            return None

        # Optional: allow `search:`-style prefixes
        raw_lower = raw.lower()
        for prefix in ("search:", "query:", "kw:", "关键词:", "关键字:"):
            if raw_lower.startswith(prefix):
                raw = raw[len(prefix):].strip()
                raw_lower = raw.lower()
                break

        # 1) 如果是完整路径
        if raw_lower.endswith('.md') or '/' in raw or '\\' in raw:
            full_path = os.path.join(base_dir, raw.replace('/', os.sep))
            if os.path.exists(full_path):
                return full_path

        # 2) 搜索技能目录（仅对“看起来像技能名”的输入做直拼路径，避免空格导致错误路径）
        if ' ' not in raw:
            skill_dirs = ['commands', 'strategies', 'states', 'general']
            for category in skill_dirs:
                # 尝试直接匹配
                skill_path = os.path.join(base_dir, 'skills', category, raw, 'SKILL.md')
                if os.path.exists(skill_path):
                    return skill_path

                # 尝试带 state- 前缀
                if category == 'states' and not raw.startswith('state-'):
                    skill_path = os.path.join(base_dir, 'skills', category, f'state-{raw}', 'SKILL.md')
                    if os.path.exists(skill_path):
                        return skill_path

        # 3) 使用 SkillLoader 搜索（精确 + 关键词/模糊）
        try:
            from .skills import get_skill_loader

            loader = get_skill_loader()

            # Exact name match first
            for s in loader.discover_skills():
                if (s.name or "").lower() == raw_lower:
                    return s.path

            # Fuzzy/keyword search
            hits = loader.search_skills(raw, limit=1)
            if hits:
                return hits[0].path
        except Exception:
            pass

        return None

    async def _inline_continue_on_read_skill(
        self,
        full_response: str,
        llm_response: LLMResponse,
        current_expr: Optional[str],
        current_state: str,
        conditions: Optional[List[str]],
        user_instruction: Optional[str],
        active_skills: List[str],
        loaded_skill_names: Set[str],
        loaded_skill_paths: Set[str],
        callback: Optional[Callable[[SolveEvent], Awaitable[None]]],
        iteration: int,
        max_rounds: int = 3,
    ) -> tuple:
        """如果模型输出 read_skill/cat，则在同一步内加载技能并继续生成最终命令。

        目的：避免“read_skill”成为单独步骤，保证在同一迭代内给出可执行命令。
        """
        from .skills import get_skill_loader

        rounds = 0
        while rounds < max_rounds:
            command = (llm_response.command or "").strip()
            if not (command.startswith("read_skill ") or command.startswith("cat ")):
                break

            target = command.split(" ", 1)[1].strip() if " " in command else ""
            if not target:
                break

            # 仅处理 SKILL.md
            if "SKILL.md" not in target and not target.endswith("SKILL.md"):
                break

            # Resolve absolute path
            try:
                base_dir = os.path.dirname(os.path.abspath(__file__))
                if os.path.isabs(target):
                    abs_path = target
                else:
                    abs_path = os.path.join(base_dir, target)
                if not os.path.exists(abs_path):
                    if os.path.exists(target):
                        abs_path = os.path.abspath(target)
                    else:
                        break
            except Exception:
                break

            abs_norm = os.path.abspath(abs_path)

            # Resolve canonical skill name + content via SkillLoader
            loader = get_skill_loader()
            meta = None
            try:
                for s in loader.discover_skills():
                    if not getattr(s, "path", None):
                        continue
                    if os.path.normcase(os.path.abspath(s.path)) == os.path.normcase(abs_norm):
                        meta = s
                        break
            except Exception:
                meta = None

            canonical_name = meta.name if meta and getattr(meta, "name", None) else os.path.basename(os.path.dirname(abs_norm))

            if canonical_name not in loaded_skill_names or abs_norm not in loaded_skill_paths:
                content_obj = loader.load_skill_content(canonical_name) if meta else None
                if content_obj:
                    skill_body = content_obj.full_content
                    skill_dir = os.path.dirname(content_obj.metadata.path)
                else:
                    with open(abs_norm, 'r', encoding='utf-8-sig') as f:
                        raw = f.read()
                    # strip frontmatter to reduce tokens
                    lines = raw.splitlines()
                    if lines and lines[0].strip() == "---":
                        for i in range(1, len(lines)):
                            if lines[i].strip() == "---":
                                raw = "\n".join(lines[i+1:]).strip()
                                break
                    skill_body = raw.strip()
                    skill_dir = os.path.dirname(abs_norm)

                active_skills.append(
                    f"\n## 技能: {canonical_name}\n"
                    f"> Path: {abs_norm}\n"
                    f"> Directory: {skill_dir}\n\n"
                    f"{skill_body}\n\n---\n"
                )
                loaded_skill_names.add(canonical_name)
                loaded_skill_paths.add(abs_norm)

                if callback:
                    await callback(SolveEvent(
                        type=EventType.SKILL_LOAD,
                        content=f"🧠 正在查阅本地技能库: {canonical_name}\n   (路径: {abs_norm})",
                        step=iteration,
                        metadata={"skill_name": canonical_name, "path": abs_norm, "round": rounds + 1, "trigger": "read_skill"}
                    ))

            # Rebuild messages with newly loaded skills and continue generation
            base_messages = self.llm._build_messages(
                current_expr or "",
                self.executor.get_history(),
                None,
                current_state=current_state,
                conditions=conditions,
                user_instruction=user_instruction,
                active_skills=active_skills
            )

            async for chunk in self.llm.continue_generation_v2(
                base_messages,
                full_response,
                canonical_name
            ):
                full_response += chunk
                if callback:
                    await callback(SolveEvent(
                        type=EventType.THINKING,
                        content=full_response,
                        step=iteration,
                        metadata={"streaming": True, "inline_read_skill": True, "round": rounds + 1}
                    ))

            llm_response = self.llm.parse_response(full_response)
            rounds += 1

        return full_response, llm_response

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
                error=init_result.error,
                total_token_usage=TokenUsage(),
                round_logs=[]
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
        max_loop_retries = max(0, int(getattr(self.config, 'max_loop_retries', 0)))  # 默认禁用循环自动重试

        # Token 统计
        total_token_usage = TokenUsage()
        round_logs: List[RoundLog] = []

        # 重置记忆（每次新的求解任务开始时清空）
        self.memory = ""

        # 动态技能状态
        active_skills: List[str] = []
        loaded_skill_names: Set[str] = set()
        loaded_skill_paths: Set[str] = set()
        loaded_skill_refs: Set[str] = set()  # raw <|load_skill|> payloads (may be keywords)

        # ============ Plan-and-Execute 状态（跨迭代持久化） ============
        # 移至外层以支持在内层循环被打断后（或进入下一轮 iteration），LLM能延续原计划
        current_plan = None
        last_plan_str = ""
        current_step_number = 0

        current_sub_plan = None
        last_sub_plan_str = ""
        current_sub_step_number = 0

        from .skills import get_skill_loader

        # ============ 新策略：元数据清单 + 按需加载 ============
        # 不再预加载完整技能内容，而是让 LLM 根据技能清单主动选择并加载
        # 技能清单（元数据）已注入 System Prompt 中
        # LLM 通过输出 <|load_skill|>技能名<|end_load_skill|> 标记动态加载
        loader = get_skill_loader()
        self.logger.info("📋 技能清单已注入 System Prompt（共 %d 个技能）",
                        len(loader.discover_skills()))

        # 思考分组：只在外层 plan step 变化时递增，使同一 TODO 内的所有 iteration 共享一个思考 UI 区块
        thinking_group = 0
        prev_step_number = 0
        accumulated_display = ""

        while iteration < self.config.max_iterations:

            iteration += 1
            current_expr = self.executor.get_current_expression()
            current_state = self.executor.get_current_state_name()

            # 迭代阶段日志
            self.logger.iteration("迭代 %d | 状态: %s", iteration, current_state)
            self.logger.info("   当前表达式: %s", current_expr[:100] if current_expr else "None")

            # 检查是否完成
            if self.executor.is_finished():
                if current_sub_step_number > 0 and callback:
                    await callback(SolveEvent(
                        type=EventType.SUB_PLAN_STEP_COMPLETE,
                        content=str(current_sub_step_number),
                        step=iteration,
                        metadata={"sub_step_number": current_sub_step_number}
                    ))
                if current_step_number > 0 and callback:
                    await callback(SolveEvent(
                        type=EventType.PLAN_STEP_COMPLETE,
                        content=str(current_step_number),
                        step=iteration,
                        metadata={"step_number": current_step_number}
                    ))
                current_sub_plan = None
                last_sub_plan_str = ""
                current_sub_step_number = 0
                current_step_number = 0
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

                if max_loop_retries <= 0:
                    self.logger.error("❌ Loop retry disabled. Aborting on first loop detection.")
                    if callback:
                        await callback(SolveEvent(
                            type=EventType.LOOP_DETECTED,
                            content="检测到表达式循环，已关闭自动重试并停止求解",
                            step=iteration,
                            metadata={"loop_expr": loop_expr, "retry_enabled": False}
                        ))
                    return SolveResult(
                        success=False,
                        final_expression=current_expr or expression,
                        latex_solution=self._format_solution_path(steps),
                        steps=steps,
                        error="检测到表达式循环（自动重试已禁用）",
                        iterations=iteration,
                        total_token_usage=total_token_usage,
                        round_logs=round_logs
                    )

                loop_retry_count += 1
                if loop_retry_count > max_loop_retries:
                    self.logger.error("❌ Max loop retries exceeded. Aborting.")
                    # 超过重试次数，停止求解
                    if callback:
                        await callback(SolveEvent(
                            type=EventType.LOOP_DETECTED,
                            content=f"检测到循环，已重试{max_loop_retries}次仍无法解决，停止求解",
                            step=iteration,
                            metadata={"loop_expr": loop_expr, "retry_enabled": True, "max_loop_retries": max_loop_retries}
                        ))
                    return SolveResult(
                        success=False,
                        final_expression=current_expr or expression,
                        latex_solution=self._format_solution_path(steps),
                        steps=steps,
                        error="检测到表达式循环，重试失败",
                        iterations=iteration,
                        total_token_usage=total_token_usage,
                        round_logs=round_logs
                    )

                # 通知LLM循环并要求重新思考
                if callback:
                    await callback(SolveEvent(
                        type=EventType.LOOP_RETRY,
                        content=f"检测到循环 (重试 {loop_retry_count}/{max_loop_retries})，尝试其他方法...",
                        step=iteration,
                        metadata={"retry_count": loop_retry_count, "loop_expr": loop_expr, "retry_enabled": True}
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
                    # 先发送初始"正在思考"状态
                    if callback:
                        await callback(SolveEvent(
                            type=EventType.THINKING,
                            content='{"thinking": "正在分析问题并决定策略...", "command": "", "explanation": "", "is_final": false}',
                            step=iteration,
                            metadata={"streaming": True, "tool_calling": True, "phase": "init"}
                        ))

                    # 工具调用回调：发送实时状态更新
                    async def on_tool_call_async(tool_name, result):
                        self.logger.info("   工具 %s 返回: %s", tool_name, result[:100] + "..." if len(result) > 100 else result)
                        # 发送技能加载状态更新
                        if callback:
                            skill_name = result.split('\n')[0] if '\n' in result else tool_name
                            skill_name = skill_name[:50] + "..." if len(skill_name) > 50 else skill_name
                            await callback(SolveEvent(
                                type=EventType.THINKING,
                                content=f'{{"thinking": "正在查阅技能文档: {tool_name}...\\n已加载并学习相关内容。", "command": "", "explanation": "", "is_final": false}}',
                                step=iteration,
                                metadata={"streaming": True, "tool_calling": True, "phase": "tool_call", "tool": tool_name}
                            ))

                    # 使用同步回调包装异步回调
                    tool_call_events = []
                    def on_tool_call(tool_name, result):
                        self.logger.info("   工具 %s 返回: %s", tool_name, result[:100] + "..." if len(result) > 100 else result)
                        tool_call_events.append((tool_name, result))

                    # 将记忆附加到 last_error 中
                    error_with_memory = last_error
                    if self.memory:
                        memory_section = f"\n\n## 历史经验教训\n\n{self.memory}"
                        error_with_memory = (last_error or "") + memory_section

                    full_response = await self.llm.generate_with_tools(
                        current_expr or expression,
                        self.executor.get_history(),
                        error_with_memory,
                        current_state,
                        conditions=conditions,
                        user_instruction=user_instruction,
                        active_skills=active_skills,
                        on_tool_call=on_tool_call
                    )

                    # 发送工具调用事件（如果有）
                    for tool_name, result in tool_call_events:
                        if callback:
                            await callback(SolveEvent(
                                type=EventType.THINKING,
                                content=f'{{"thinking": "已查阅技能: {tool_name}\\n正在根据技能内容决定下一步...", "command": "", "explanation": "", "is_final": false}}',
                                step=iteration,
                                metadata={"streaming": True, "tool_calling": True, "phase": "tool_result", "tool": tool_name}
                            ))

                    if callback:
                        await callback(SolveEvent(
                            type=EventType.THINKING,
                            content=full_response,
                            step=iteration,
                            metadata={"streaming": False, "tool_calling": True, "phase": "complete"}
                        ))

                    # Search-o1 风格：检测 Function Calling 响应中的技能加载标记
                    full_response, skills_loaded, loaded_names = self._detect_and_load_skills(
                        full_response, active_skills, loaded_skill_names, callback
                    )
                    if skills_loaded and callback:
                        for name in loaded_names:
                            await callback(SolveEvent(
                                type=EventType.SKILL_LOAD,
                                content=f"已加载技能: {name}",
                                step=iteration,
                                metadata={"skill_name": name, "trigger": "marker"}
                            ))
                else:
                    # ============ Search-o1 风格：流式生成 + 边思考边查询 ============
                    from .prompts import SKILL_LOAD_BEGIN, SKILL_LOAD_END, EXECUTE_BEGIN, EXECUTE_END
                    from .prompts import PLAN_BEGIN, PLAN_END, PLAN_STEP_BEGIN, PLAN_STEP_END
                    from .prompts import SUB_PLAN_BEGIN, SUB_PLAN_END, SUB_STEP_BEGIN, SUB_STEP_END
                    import re
                    import json

                    # 保存原始消息用于继续生成
                    base_messages = self.llm._build_messages(
                        current_expr or expression,
                        self.executor.get_history(),
                        last_error,
                        current_state=current_state,
                        conditions=conditions,
                        user_instruction=user_instruction,
                        active_skills=active_skills
                    )

                    max_action_rounds = 30  # 最多加载30次技能或执行命令
                    action_round = 0
                    executed_commands_count = 0  # 追踪已经执行了几个 <|execute|> 块

                    # 思考分组：有计划时按外层步骤分组，无计划时每次 iteration 独立
                    if current_plan is None or current_step_number != prev_step_number:
                        thinking_group += 1
                        prev_step_number = current_step_number
                        display_text = ""
                        accumulated_display = ""
                    else:
                        # 同一外层步骤内继续累积
                        display_text = accumulated_display

                    raw_text = ""  # 维护原始文本用于 LLM 上下文（保留标记）
                    is_first_generation = True
                    last_continuation_text = ""  # 追踪最后一轮续写输出，供 parse_response 使用
                    pre_skill_display = ""  # 追踪续写前的 display_text 快照
                    continuation_nudge_count = 0
                    max_continuation_nudges = 6

                    def _get_next_declared_step(plan_items, current_step):
                        declared_steps = sorted({
                            int(item.get('step', 0))
                            for item in (plan_items or [])
                            if int(item.get('step', 0)) > 0
                        })
                        for step_number in declared_steps:
                            if step_number > current_step:
                                return step_number
                        return None

                    def _get_step_description(plan_items, step_number):
                        for item in (plan_items or []):
                            try:
                                if int(item.get('step', 0)) == step_number:
                                    return item.get('description', '')
                            except Exception:
                                continue
                        return ""

                    def _get_sub_plan_progress():
                        declared_steps = sorted({
                            int(item.get('step', 0))
                            for item in (current_sub_plan or [])
                            if int(item.get('step', 0)) > 0
                        })
                        if not declared_steps:
                            return {
                                "has_plan": False,
                                "next_step": None,
                                "max_step": 0,
                                "has_started": False,
                                "needs_more_steps": False,
                                "is_complete": False,
                            }

                        next_step = None
                        if current_sub_step_number <= 0:
                            next_step = declared_steps[0]
                        else:
                            for step_number in declared_steps:
                                if step_number > current_sub_step_number:
                                    next_step = step_number
                                    break

                        max_step = declared_steps[-1]
                        has_started = current_sub_step_number > 0
                        is_complete = has_started and current_sub_step_number >= max_step

                        return {
                            "has_plan": True,
                            "next_step": next_step,
                            "max_step": max_step,
                            "has_started": has_started,
                            "needs_more_steps": next_step is not None,
                            "is_complete": is_complete,
                        }

                    def _has_unprocessed_action(current_raw_text: str) -> bool:
                        pattern_exec = re.escape(EXECUTE_BEGIN) + r'(.+?)' + re.escape(EXECUTE_END)
                        matches_exec = list(re.finditer(pattern_exec, current_raw_text, re.DOTALL))
                        if len(matches_exec) > executed_commands_count:
                            return True

                        pattern_skill = re.escape(SKILL_LOAD_BEGIN) + r'(.+?)' + re.escape(SKILL_LOAD_END)
                        for match in re.finditer(pattern_skill, current_raw_text, re.DOTALL):
                            if match.group(1).strip() not in loaded_skill_refs:
                                return True
                        return False

                    async def _complete_active_plan_scope():
                        nonlocal current_step_number
                        nonlocal current_sub_plan, last_sub_plan_str, current_sub_step_number

                        completed_step_number = current_step_number
                        completed_sub_step_number = current_sub_step_number

                        if completed_sub_step_number > 0 and callback:
                            await callback(SolveEvent(
                                type=EventType.SUB_PLAN_STEP_COMPLETE,
                                content=str(completed_sub_step_number),
                                step=iteration,
                                metadata={"sub_step_number": completed_sub_step_number}
                            ))

                        current_sub_plan = None
                        last_sub_plan_str = ""
                        current_sub_step_number = 0

                        if completed_step_number > 0 and callback:
                            await callback(SolveEvent(
                                type=EventType.PLAN_STEP_COMPLETE,
                                content=str(completed_step_number),
                                step=iteration,
                                metadata={"step_number": completed_step_number}
                            ))

                        current_step_number = 0

                    async def _continue_generation_with_nudge(
                        nudge_msg: str,
                        continue_key: str,
                        ui_message: str,
                        log_label: str,
                    ) -> bool:
                        nonlocal raw_text, display_text, pre_skill_display
                        nonlocal last_continuation_text, continuation_nudge_count

                        if continuation_nudge_count >= max_continuation_nudges:
                            self.logger.warning("⚠ 连续续写提示过多，停止继续追问当前计划节点")
                            return False

                        continuation_nudge_count += 1
                        nudge = f"\n<|observation|>\n[系统提示: {nudge_msg}]\n<|end_observation|>\n\n"
                        raw_text += nudge
                        if ui_message:
                            display_text += f"\n[{ui_message}]\n"

                        self.logger.info("   🔄 %s", log_label)

                        base_messages = self.llm._build_messages(
                            current_expr or expression,
                            self.executor.get_history(),
                            None,
                            current_state=current_state,
                            conditions=conditions,
                            user_instruction=user_instruction,
                            active_skills=active_skills
                        )

                        pre_skill_display = display_text
                        continuation_text = ""
                        last_continuation_text = ""

                        try:
                            async for chunk in self.llm.continue_generation_v2(
                                base_messages,
                                raw_text,
                                continue_key
                            ):
                                continuation_text += chunk
                                display_text += chunk
                                raw_text += chunk
                                await _check_plan_events(raw_text)

                                if callback:
                                    await callback(SolveEvent(
                                        type=EventType.THINKING,
                                        content=display_text,
                                        step=thinking_group,
                                        metadata={
                                            "streaming": True,
                                            "action_round": action_round,
                                            "pre_skill_display": ""
                                        }
                                    ))

                                if _has_unprocessed_action(raw_text):
                                    self.logger.info("   ⏸ %s中检测到新动作", log_label)
                                    break

                            last_continuation_text = continuation_text
                            self.logger.info("   ✅ %s: %d 字符", log_label, len(continuation_text))
                            return True
                        except Exception as e:
                            self.logger.warning("%s失败: %s", log_label, str(e))
                            return False

                    async def _check_plan_events(current_raw_text: str):
                        nonlocal last_plan_str, current_plan, current_step_number
                        nonlocal last_sub_plan_str, current_sub_plan, current_sub_step_number

                        pattern_plan = re.escape(PLAN_BEGIN) + r'(.*?)' + re.escape(PLAN_END)
                        matches_plan = list(re.finditer(pattern_plan, current_raw_text, re.DOTALL))
                        if matches_plan:
                            latest_plan_str = matches_plan[-1].group(1).strip()
                            if latest_plan_str != last_plan_str:
                                last_plan_str = latest_plan_str
                                try:
                                    plan_data = json.loads(latest_plan_str)
                                    is_update = current_plan is not None
                                    current_plan = plan_data
                                    if callback:
                                        await callback(SolveEvent(
                                            type=EventType.PLAN_UPDATED if is_update else EventType.PLAN_GENERATED,
                                            content=latest_plan_str,
                                            step=iteration,
                                            metadata={"active_step_number": current_step_number}
                                        ))
                                except Exception as e:
                                    self.logger.warning(f"解析计划 JSON 失败: {e}")

                        pattern_step = re.escape(PLAN_STEP_BEGIN) + r'(\d+)' + re.escape(PLAN_STEP_END)
                        matches_step = list(re.finditer(pattern_step, current_raw_text))
                        if matches_step:
                            latest_step_str = matches_step[-1].group(1)
                            try:
                                step_num = int(latest_step_str)
                                if step_num != current_step_number:
                                    await _complete_active_plan_scope()
                                    current_step_number = step_num
                                    if callback:
                                        await callback(SolveEvent(
                                            type=EventType.PLAN_STEP_START,
                                            content=str(current_step_number),
                                            step=iteration,
                                            metadata={"step_number": current_step_number}
                                        ))
                            except ValueError:
                                pass

                        # ========= 检测内层子计划 (sub_plan) =========

                        pattern_sub_plan = re.escape(SUB_PLAN_BEGIN) + r'([\s\S]*?)' + re.escape(SUB_PLAN_END)
                        matches_sub_plan = list(re.finditer(pattern_sub_plan, current_raw_text))
                        if matches_sub_plan:
                            latest_sub_plan_str = matches_sub_plan[-1].group(1).strip()
                            if latest_sub_plan_str != last_sub_plan_str:
                                last_sub_plan_str = latest_sub_plan_str
                                try:
                                    sub_plan_data = json.loads(latest_sub_plan_str)
                                    current_sub_plan = sub_plan_data
                                    if callback:
                                        await callback(SolveEvent(
                                            type=EventType.SUB_PLAN_GENERATED,
                                            content=latest_sub_plan_str,
                                            step=iteration,
                                            metadata={
                                                "parent_step_number": current_step_number,
                                                "active_sub_step_number": current_sub_step_number
                                            }
                                        ))
                                except Exception as e:
                                    self.logger.warning(f"解析子计划 JSON 失败: {e}")

                        pattern_sub_step = re.escape(SUB_STEP_BEGIN) + r'(\d+)' + re.escape(SUB_STEP_END)
                        matches_sub_step = list(re.finditer(pattern_sub_step, current_raw_text))
                        if matches_sub_step:
                            latest_sub_step_str = matches_sub_step[-1].group(1)
                            try:
                                sub_step_num = int(latest_sub_step_str)
                                if sub_step_num != current_sub_step_number:
                                    if current_sub_step_number > 0 and callback:
                                        await callback(SolveEvent(
                                            type=EventType.SUB_PLAN_STEP_COMPLETE,
                                            content=str(current_sub_step_number),
                                            step=iteration,
                                            metadata={"sub_step_number": current_sub_step_number}
                                        ))
                                    current_sub_step_number = sub_step_num
                                    if callback:
                                        await callback(SolveEvent(
                                            type=EventType.SUB_PLAN_STEP_START,
                                            content=str(current_sub_step_number),
                                            step=iteration,
                                            metadata={"sub_step_number": current_sub_step_number}
                                        ))
                            except ValueError:
                                pass

                    while action_round < max_action_rounds:
                        # 开始一轮生成（第一轮或继续生成）
                        self.logger.info("🔄 Search-o1/ReAct 循环: 第 %d 轮 (is_first=%s)", action_round, is_first_generation)
                        round_text = ""
                        found_new_action = False
                        action_type = None  # 'skill' or 'execute'
                        action_payload = None

                        if is_first_generation:
                            # 第一轮：正常流式生成
                            # 将记忆附加到 last_error 中
                            error_with_memory = last_error
                            if self.memory:
                                memory_section = f"\n\n## 历史经验教训\n\n{self.memory}"
                                error_with_memory = (last_error or "") + memory_section

                            # 将遗留计划注入当前状态（避免重做）
                            if last_plan_str:
                                plan_reminder = f"\n\n[系统提示] 你目前正在执行此前未完成的计划，当前计划内容如下：\n<|plan|>\n{last_plan_str}\n<|end_plan|>\n请直接从当前最新的 `<|step|>` 开始继续执行，无需再重新输出原有计划。"
                                error_with_memory = (error_with_memory or "") + plan_reminder

                            async for chunk in self.llm.generate_command(
                                current_expr or expression,
                                self.executor.get_history(),
                                error_with_memory,
                                current_state,
                                conditions=conditions,
                                user_instruction=user_instruction,
                                active_skills=active_skills
                            ):
                                round_text += chunk
                                display_text += chunk
                                raw_text += chunk
                                await _check_plan_events(raw_text)

                                # 实时发送流式更新
                                if callback:
                                    await callback(SolveEvent(
                                        type=EventType.THINKING,
                                        content=display_text,
                                        step=thinking_group,
                                        metadata={"streaming": True, "action_round": action_round}
                                    ))

                                # 实时检测技能标记和命令标记
                                # 先规范化可能的错误格式标记
                                raw_text_normalized = re.sub(r'<load_skill>([^<]+)</load_skill>',
                                                              f'{SKILL_LOAD_BEGIN}\\1{SKILL_LOAD_END}', raw_text)
                                raw_text_normalized = re.sub(r'<<\s*([a-zA-Z_-]+)\s*>>',
                                                              f'{SKILL_LOAD_BEGIN}\\1{SKILL_LOAD_END}', raw_text_normalized)
                                if raw_text_normalized != raw_text:
                                    raw_text = raw_text_normalized
                                    display_text = re.sub(r'<load_skill>([^<]+)</load_skill>',
                                                          f'{SKILL_LOAD_BEGIN}\\1{SKILL_LOAD_END}', display_text)
                                    display_text = re.sub(r'<<\s*([a-zA-Z_-]+)\s*>>',
                                                          f'{SKILL_LOAD_BEGIN}\\1{SKILL_LOAD_END}', display_text)

                                # 检查 <|load_skill|>
                                pattern_skill = re.escape(SKILL_LOAD_BEGIN) + r'(.+?)' + re.escape(SKILL_LOAD_END)
                                matches_skill = list(re.finditer(pattern_skill, raw_text, re.DOTALL))
                                for match in matches_skill:
                                    skill_name = match.group(1).strip()
                                    if skill_name not in loaded_skill_refs:
                                        found_new_action = True
                                        action_type = "skill"
                                        action_payload = skill_name

                                        marker_end_pos = match.end()
                                        raw_text = raw_text[:marker_end_pos]
                                        display_text = display_text[:marker_end_pos]
                                        self.logger.info("   ✂️ 截断到技能标记位置: %d 字符", marker_end_pos)
                                        break
                                    else:
                                        # 技能已经加载过，替换掉原始标记为加载确认，打破LLM在后续续写中重复输出标签的死循环
                                        marker_str = match.group(0)
                                        if f"[✓ 已加载技能: {skill_name}]" not in display_text:
                                            # 使用不会出错的正则替换（注意转义）
                                            # 注意：由于流式输出中可能存在转义，使用安全的替换方式
                                            # 我们直接使用简单的 replace 即可，因为 marker_str 是从 raw_text 提取的精确字符串
                                            raw_text = raw_text.replace(marker_str, f"\n[✓ 已加载技能: {skill_name}]\n", 1)
                                            display_text = display_text.replace(marker_str, f"\n[✓ 已加载技能: {skill_name}]\n", 1)
                                if found_new_skill:
                                    break

                            is_first_generation = False

                        # 检测是否有未加载的技能或未执行的命令
                        if not found_new_action:
                            pattern_skill = re.escape(SKILL_LOAD_BEGIN) + r'(.+?)' + re.escape(SKILL_LOAD_END)
                            for match in re.finditer(pattern_skill, raw_text, re.DOTALL):
                                if match.group(1).strip() not in loaded_skill_refs:
                                    found_new_action = True
                                    action_type = "skill"
                                    action_payload = match.group(1).strip()
                                    break

                            if not found_new_action:
                                pattern_exec = re.escape(EXECUTE_BEGIN) + r'(.+?)' + re.escape(EXECUTE_END)
                                matches_exec = list(re.finditer(pattern_exec, raw_text, re.DOTALL))
                                if len(matches_exec) > executed_commands_count:
                                    match = matches_exec[executed_commands_count]
                                    found_new_action = True
                                    action_type = "execute"
                                    action_payload = match.group(1).strip()

                        # 如果没有任何新动作，检查是否需要继续推进当前计划节点
                        if not found_new_action:
                            sub_plan_progress = _get_sub_plan_progress()
                            next_outer_step = _get_next_declared_step(current_plan, current_step_number)
                            has_is_final = '"is_final": true' in raw_text or '"is_final":true' in raw_text

                            if current_step_number > 0 and sub_plan_progress["has_plan"] and sub_plan_progress["needs_more_steps"]:
                                next_sub_step = sub_plan_progress["next_step"]
                                nudge_msg = (
                                    f"你仍在执行外层步骤{current_step_number}。当前子计划还没有完成，不要跳到新的 <|step|>。"
                                )
                                if next_sub_step is not None:
                                    nudge_msg += (
                                        f" 如果当前子步骤已经完成，请立刻开始 <|sub_step|>{next_sub_step}<|end_sub_step|>；"
                                        "如果需要调整方法，请重新输出 <|sub_plan|>，但只保留剩余子步骤，然后继续执行。"
                                    )
                                    ui_message = f"当前 todo 尚未完成，继续推进子步骤 {next_sub_step}..."
                                else:
                                    nudge_msg += " 请先开始当前子计划的第一个 <|sub_step|>，或重写剩余子计划后继续执行。"
                                    ui_message = "当前 todo 尚未完成，继续推进当前子计划..."

                                if await _continue_generation_with_nudge(
                                    nudge_msg,
                                    f"subplan_step_{current_step_number}_{next_sub_step or 'resume'}",
                                    ui_message,
                                    "子计划续写"
                                ):
                                    continue

                            elif current_step_number > 0 and sub_plan_progress["has_plan"] and sub_plan_progress["is_complete"] and next_outer_step is not None:
                                next_desc = _get_step_description(current_plan, next_outer_step)
                                nudge_msg = f"外层步骤{current_step_number}的当前子计划已经执行到最后一个子步骤。"
                                if next_desc:
                                    nudge_msg += (
                                        f" 如果该外层步骤已经完成，请继续 <|step|>{next_outer_step}<|end_step|> 并执行步骤{next_outer_step}: {next_desc}。"
                                    )
                                else:
                                    nudge_msg += f" 如果该外层步骤已经完成，请继续 <|step|>{next_outer_step}<|end_step|>。"
                                nudge_msg += " 如果该外层步骤其实还没完成，请在当前步骤内重新输出只包含剩余任务的 <|sub_plan|>，然后继续执行。不要直接结束。"

                                if await _continue_generation_with_nudge(
                                    nudge_msg,
                                    f"plan_step_{next_outer_step}",
                                    f"外层步骤 {current_step_number} 子计划已收束，准备切换到步骤 {next_outer_step}...",
                                    "外层步骤续写"
                                ):
                                    continue

                            elif current_plan and current_step_number > 0 and next_outer_step is not None:
                                next_desc = _get_step_description(current_plan, next_outer_step)
                                nudge_msg = f"步骤{current_step_number}已完成。请继续执行计划的步骤{next_outer_step}"
                                if next_desc:
                                    nudge_msg += f": {next_desc}"
                                nudge_msg += "。不要重新制定计划，直接用 <|step|> 标记并执行。"

                                if await _continue_generation_with_nudge(
                                    nudge_msg,
                                    f"plan_step_{next_outer_step}",
                                    f"步骤 {current_step_number} 完成，继续执行步骤 {next_outer_step}...",
                                    "外层步骤续写"
                                ):
                                    continue

                            elif current_plan and current_step_number > 0 and not has_is_final:
                                if sub_plan_progress["has_plan"] and not sub_plan_progress["is_complete"]:
                                    nudge_msg = "宏观计划已经走到最后一个外层步骤，但当前子计划还没有闭合。请先补齐剩余 <|sub_step|>，或重写当前 <|sub_plan|> 后继续执行，不要直接结束。"
                                    ui_message = "最后一个 todo 还没有闭合，继续补齐子计划..."
                                    continue_key = f"subplan_finalize_{current_step_number}"
                                else:
                                    nudge_msg = "你的所有计划步骤都已执行完毕。请评估当前表达式：如果已经化简到最终结果，请直接输出 is_final: true 的 JSON 结束；如果还需要更多操作，请用 <|plan|> 输出补充计划并继续执行。"
                                    ui_message = "计划全部完成，等待确认或补充计划..."
                                    continue_key = "plan_finalize"

                                if await _continue_generation_with_nudge(
                                    nudge_msg,
                                    continue_key,
                                    ui_message,
                                    "计划终态续写"
                                ):
                                    continue

                            # 无计划或计划已完成且已确认，正常退出内层循环
                            self.logger.info("   ⏳ 没有检测到新的动作指令(<|execute|> / <|load_skill|>)，结束内层循环等待返回...")
                            break

                        action_round += 1
                        continuation_nudge_count = 0

                        # 处理技能加载
                        if action_type == "skill":
                            requested_skill = action_payload
                            self.logger.info("🔖 Search-o1 第%d轮: 加载技能 %s", action_round, requested_skill)


                            try:
                                skill_path = self._resolve_skill_path(requested_skill)
                                if not skill_path:
                                    self.logger.warning("技能未找到: %s", requested_skill)
                                    # 替换标记为错误提示 (UI)
                                    pattern = re.escape(SKILL_LOAD_BEGIN) + re.escape(requested_skill) + re.escape(SKILL_LOAD_END)
                                    display_text = re.sub(pattern, f"[技能 {requested_skill} 未找到]", display_text)
                                    # 记录该 marker，防止死循环
                                    loaded_skill_refs.add(requested_skill)
                                    continue

                                # Resolve canonical name (for de-dup) + load body (strip frontmatter) via SkillLoader
                                loader = get_skill_loader()
                                try:
                                    meta = next((s for s in loader.discover_skills() if s.path == skill_path), None)
                                except Exception:
                                    meta = None

                                canonical_name = meta.name if meta and meta.name else requested_skill

                                # Mark the marker payload as seen
                                loaded_skill_refs.add(requested_skill)
                                # Also mark canonical skill name so future markers using the name don't interrupt again
                                if canonical_name:
                                    loaded_skill_refs.add(canonical_name)

                                # Only load content if this canonical skill is new
                                if canonical_name not in loaded_skill_names:
                                    content_obj = loader.load_skill_content(canonical_name)
                                    if content_obj:
                                        skill_body = content_obj.full_content
                                        skill_dir = os.path.dirname(content_obj.metadata.path)
                                    else:
                                        with open(skill_path, 'r', encoding='utf-8-sig') as f:
                                            skill_body = f.read()
                                        skill_dir = os.path.dirname(skill_path)

                                    loaded_skill_names.add(canonical_name)
                                    loaded_skill_paths.add(skill_path)
                                    active_skills.append(
                                        f"\n## 技能: {canonical_name}\n"
                                        f"> Path: {skill_path}\n"
                                        f"> Directory: {skill_dir}\n\n"
                                        f"{skill_body}\n\n---\n"
                                    )

                                    # 发送技能加载通知
                                    if callback:
                                        await callback(SolveEvent(
                                            type=EventType.SKILL_LOAD,
                                            content=f"🧠 正在查阅本地技能库: {canonical_name}\n   (路径: {skill_path})",
                                            step=iteration,
                                            metadata={
                                                "skill_name": canonical_name,
                                                "path": skill_path,
                                                "round": action_round,
                                                "requested": requested_skill,
                                            }
                                        ))

                                # Rebuild base_messages with (potentially) new loaded skills
                                base_messages = self.llm._build_messages(
                                    current_expr or expression,
                                    self.executor.get_history(),
                                    last_error,
                                    current_state=current_state,
                                    conditions=conditions,
                                    user_instruction=user_instruction,
                                    active_skills=active_skills
                                )

                                # 打印日志
                                print("")
                                print_separator("═", 60)
                                self.logger.skill_load(f"🔮 Search-o1 第{action_round}轮: {canonical_name}")
                                self.logger.info("   📁 路径: %s", skill_path)
                                print_separator("═", 60)

                                # 替换标记为简洁提示 (仅在 UI display_text 中替换)
                                pattern_marker = re.escape(SKILL_LOAD_BEGIN) + re.escape(requested_skill) + re.escape(SKILL_LOAD_END)
                                pattern_str = r'\\?\s*' + pattern_marker
                                display_text = re.sub(pattern_str, f"\n[✓ 已加载技能: {canonical_name}]\n", display_text, count=1)

                                # 发送技能加载完成状态 (更新 Thinking 内容)
                                if callback:
                                    await callback(SolveEvent(
                                        type=EventType.THINKING,
                                        content=display_text,
                                        step=thinking_group,
                                        metadata={"streaming": True, "action_round": action_round, "skill_loaded": canonical_name}
                                    ))

                                self.logger.info("   🔄 继续生成...")

                                # 保存续写前的 display_text 快照（旧思考 + 技能标记）
                                pre_skill_display = display_text

                                # 继续生成（使用 assistant prefix 模式，技能内容已注入 system prompt）
                                continuation_text = ""
                                last_continuation_text = ""  # 重置，准备追踪本轮续写
                                async for chunk in self.llm.continue_generation_v2(
                                    base_messages,
                                    raw_text,
                                    canonical_name
                                ):
                                    continuation_text += chunk
                                    display_text += chunk
                                    raw_text += chunk

                                    # 实时发送流式更新
                                    # 核心：发送 continuation_text（单 JSON），旧思考通过 metadata 传递
                                    # 这样 app_flask.py 的 parser 始终只处理单 JSON，无需多 JSON 拼接解析
                                    if callback:
                                        await callback(SolveEvent(
                                            type=EventType.THINKING,
                                            content=continuation_text,
                                            step=thinking_group,
                                            metadata={
                                                "streaming": True,
                                                "action_round": action_round,
                                                "pre_skill_display": pre_skill_display
                                            }
                                        ))

                                    # 实时检测新的技能标记
                                    # 🔧 先规范化可能的错误格式标记
                                    raw_text = re.sub(r'<load_skill>([^<]+)</load_skill>',
                                                      f'{SKILL_LOAD_BEGIN}\\1{SKILL_LOAD_END}', raw_text)
                                    raw_text = re.sub(r'<<\s*([a-zA-Z_-]+)\s*>>',
                                                      f'{SKILL_LOAD_BEGIN}\\1{SKILL_LOAD_END}', raw_text)
                                    display_text = re.sub(r'<load_skill>([^<]+)</load_skill>',
                                                          f'{SKILL_LOAD_BEGIN}\\1{SKILL_LOAD_END}', display_text)
                                    display_text = re.sub(r'<<\s*([a-zA-Z_-]+)\s*>>',
                                                          f'{SKILL_LOAD_BEGIN}\\1{SKILL_LOAD_END}', display_text)

                                    found_new_action_in_chunk = False
                                    pattern_exec = re.escape(EXECUTE_BEGIN) + r'(.+?)' + re.escape(EXECUTE_END)
                                    matches_exec = list(re.finditer(pattern_exec, raw_text, re.DOTALL))
                                    if len(matches_exec) > executed_commands_count:
                                        found_new_action_in_chunk = True

                                    if not found_new_action_in_chunk:
                                        pattern_skill = re.escape(SKILL_LOAD_BEGIN) + r'(.+?)' + re.escape(SKILL_LOAD_END)
                                        matches_skill = list(re.finditer(pattern_skill, raw_text, re.DOTALL))
                                        for m in matches_skill:
                                            if m.group(1).strip() not in loaded_skill_refs:
                                                found_new_action_in_chunk = True
                                                break

                                    if found_new_action_in_chunk:
                                        self.logger.info("   ⏸ 检测到新动作请求，准备下一轮")
                                        break

                                # 日志：记录续写结果
                                self.logger.info("   ✅ 续写完成: %d 字符", len(continuation_text))
                                last_continuation_text = continuation_text  # 保存最后一轮续写结果

                            except Exception as e:
                                self.logger.warning("技能加载失败: %s - %s", requested_skill, str(e))
                                display_text += f"\n[技能 {requested_skill} 加载失败: {str(e)}]\n"
                                loaded_skill_refs.add(requested_skill)

                        elif action_type == "execute":
                            command = action_payload
                            executed_commands_count += 1
                            self.logger.execution("⚡ ReAct 第%d轮: 立即执行命令 %s", action_round, command)

                            # 执行命令
                            exec_result = self.executor.execute(command)

                            # 收集结果
                            if exec_result.success:
                                self.logger.result("成功 | 变更=%s", exec_result.changed)
                                observation = exec_result.result if exec_result.result else "成功 (无新结果返回)"
                                observation_str = f"成功 | 结果:\n{observation}"
                                last_error = None
                                self.error_handler.reset_consecutive_errors()

                                # 将步骤添加到列表
                                steps.append({
                                    "step": iteration,
                                    "command": command,
                                    "thinking": "",  # 流式过程中的步骤
                                    "explanation": "系统自动提取的执行结果",
                                    "success": True,
                                    "expr_before": current_expr,
                                    "expr_after": exec_result.result,
                                    "latex_after": exec_result.latex_result,
                                    "error": None,
                                    "changed": exec_result.changed
                                })
                                current_expr = exec_result.result

                                # 🎯 关键修复：执行命令后立即检测 executor 是否已完成
                                # 当 done 命令使状态回到 InitialState 时，不再续写，直接跳出内层循环
                                if self.executor.is_finished():
                                    self.logger.complete("🎯 executor 检测到证明完成 (执行 '%s' 后状态回到 initial)", command)
                                    break

                                if callback:
                                    import asyncio
                                    await callback(SolveEvent(
                                        type=EventType.COMMAND_SUCCESS,
                                        content=command,
                                        step=iteration,
                                        metadata={"changed": exec_result.changed}
                                    ))
                                    await callback(SolveEvent(
                                        type=EventType.RESULT,
                                        content=exec_result.result or "",
                                        latex=exec_result.latex_result,
                                        step=iteration,
                                        metadata={"changed": exec_result.changed}
                                    ))
                            else:
                                self.logger.error_phase("执行失败: %s", exec_result.error)
                                observation_str = f"失败 | 错误:\n{exec_result.error}"
                                error_info = self.error_handler.handle(
                                    Exception(exec_result.error or "未知错误"),
                                    {"command": command, "expression": current_expr}
                                )
                                last_error = self.error_handler.format_for_llm(error_info)

                                steps.append({
                                    "step": iteration,
                                    "command": command,
                                    "thinking": "",
                                    "explanation": "",
                                    "success": False,
                                    "expr_before": current_expr,
                                    "expr_after": None,
                                    "latex_after": None,
                                    "error": exec_result.error,
                                    "changed": False
                                })

                                if callback:
                                    import asyncio
                                    await callback(SolveEvent(
                                        type=EventType.COMMAND_FAILURE,
                                        content=command,
                                        step=iteration,
                                        metadata={"error": exec_result.error}
                                    ))
                                    await callback(SolveEvent(
                                        type=EventType.ERROR,
                                        content=exec_result.error or "执行失败",
                                        step=iteration
                                    ))

                            # 更新上下文
                            replacement = f"\n<|observation|>\n{observation_str}\n<|end_observation|>\n"

                            # 记录快照
                            pre_skill_display = display_text + replacement

                            # 继续生成
                            self.logger.info("   🔄 观察结果后继续生成...")
                            continuation_text = ""
                            last_continuation_text = ""

                            raw_text += replacement
                            display_text += f"\n[✓ 执行结果: {observation_str}]\n"

                            base_messages = self.llm._build_messages(
                                current_expr or expression,
                                self.executor.get_history(),
                                None,
                                current_state=current_state,
                                conditions=conditions,
                                user_instruction=user_instruction,
                                active_skills=active_skills
                            )

                            try:
                                async for chunk in self.llm.continue_generation_v2(
                                    base_messages,
                                    raw_text,
                                    command
                                ):
                                    continuation_text += chunk
                                    display_text += chunk
                                    raw_text += chunk
                                    await _check_plan_events(raw_text)

                                    if callback:
                                        import asyncio
                                        await callback(SolveEvent(
                                            type=EventType.THINKING,
                                            content=continuation_text,
                                            step=thinking_group,
                                            metadata={
                                                "streaming": True,
                                                "action_round": action_round,
                                                "pre_skill_display": pre_skill_display
                                            }
                                        ))

                                    # 检测新的动作
                                    found_new_action_in_chunk = False
                                    pattern_exec = re.escape(EXECUTE_BEGIN) + r'(.+?)' + re.escape(EXECUTE_END)
                                    matches_exec = list(re.finditer(pattern_exec, raw_text, re.DOTALL))
                                    if len(matches_exec) > executed_commands_count:
                                        found_new_action_in_chunk = True

                                    if not found_new_action_in_chunk:
                                        pattern_skill = re.escape(SKILL_LOAD_BEGIN) + r'(.+?)' + re.escape(SKILL_LOAD_END)
                                        matches_skill = list(re.finditer(pattern_skill, raw_text, re.DOTALL))
                                        for m in matches_skill:
                                            if m.group(1).strip() not in loaded_skill_refs:
                                                found_new_action_in_chunk = True
                                                break

                                    if found_new_action_in_chunk:
                                        self.logger.info("   ⏸ 检测到新动作请求，准备下一轮")
                                        break
                                    else:
                                        marker_str = match.group(0)
                                        if f"[✓ 已加载技能: {next_skill}]" not in display_text:
                                            raw_text = raw_text.replace(marker_str, f"\n[✓ 已加载技能: {next_skill}]\n", 1)
                                            display_text = display_text.replace(marker_str, f"\n[✓ 已加载技能: {next_skill}]\n", 1)
                                            continuation_text = continuation_text.replace(marker_str, f"\n[✓ 已加载技能: {next_skill}]\n", 1)
                                if found_new_skill:
                                    break
                            
                            # 日志：记录续写结果
                            self.logger.info("   ✅ 续写完成: %d 字符", len(continuation_text))
                            last_continuation_text = continuation_text  # 保存最后一轮续写结果

                                self.logger.info("   ✅ 续写完成: %d 字符", len(continuation_text))
                                last_continuation_text = continuation_text

                            except Exception as e:
                                self.logger.warning("执行后继续生成失败: %s", str(e))
                                display_text += f"\n[执行后继续生成失败: {str(e)}]\n"

                    # 循环结束：可能是达到最大轮次或没有更多技能/动作
                    # 保存累积的思考文本，供下一个 iteration（同一外层步骤）恢复使用
                    accumulated_display = display_text

                    # 核心修复：如果有续写输出，只把续写部分传给 parse_response
                    # 这避免了 parse_response 扫描整个多 JSON 拼接的 display_text，
                    # 从而防止策略3从旧 thinking 中误匹配 rewrite/substitute 等词作为命令
                    full_response = last_continuation_text if last_continuation_text else display_text

                    if action_round >= max_action_rounds:
                        self.logger.warning("⚠ 达到最大动作轮次 (%d)", max_action_rounds)
                        full_response += f"\n\n[系统提示: 已达到最大动作轮次 {max_action_rounds}]"

                    # 发送最终响应（含 pre_skill_display 供前端正确显示旧思考）
                    if callback:
                        await callback(SolveEvent(
                            type=EventType.THINKING,
                            content=full_response,
                            step=thinking_group,
                            metadata={
                                "streaming": False,
                                "action_rounds": action_round,
                                "pre_skill_display": pre_skill_display if last_continuation_text else ""
                            }
                        ))

                llm_response = self.llm.parse_response(full_response)

                # ReAct completion fix: if command is empty but is_final is true
                if executed_commands_count > 0 and llm_response.is_final:
                    llm_response.command = ""

                # 🔧 关键修复：如果技能加载后没有有效命令，强制继续生成
                # 这确保模型在同一步骤内完成技能加载和命令输出
                if not llm_response.command and loaded_skill_names:
                    self.logger.warning("⚠ 技能已加载但未生成有效命令，强制继续生成...")

                    # 使用最后加载的技能名称
                    last_skill_name = list(loaded_skill_names)[-1]

                    # 重建 messages 并调用回退模式
                    # 🔧 用 display_text（完整上下文）而非 full_response（可能只有 "{"）
                    base_messages = self.llm._build_messages(
                        expression=current_expr,
                        history=self.executor.get_history(),
                        last_error=None,
                        current_state=current_state,
                        conditions=conditions,
                        user_instruction=user_instruction,
                        active_skills=active_skills
                    )

                    # 添加当前输出作为助手消息（用 display_text 保留完整上下文）
                    base_messages.append({
                        "role": "assistant",
                        "content": display_text
                    })

                    # 调用回退模式继续生成
                    fallback_text = ""
                    # 确定旧思考前缀（用于前端显示）
                    fallback_pre_skill = pre_skill_display if last_continuation_text else display_text
                    async for chunk in self.llm._continue_fallback(base_messages, last_skill_name):
                        fallback_text += chunk
                        if callback:
                            await callback(SolveEvent(
                                type=EventType.THINKING,
                                content=fallback_text,
                                step=thinking_group,
                                metadata={
                                    "streaming": True,
                                    "continuation": True,
                                    "pre_skill_display": fallback_pre_skill
                                }
                            ))

                    # parse_response 只传回退文本（干净单 JSON）
                    full_response = fallback_text
                    llm_response = self.llm.parse_response(full_response)
                    self.logger.info("✅ 续写完成，命令: %s", llm_response.command[:50] if llm_response.command else "无")

                # 如果模型输出 read_skill/cat，则在同一步内继续生成最终命令
                full_response, llm_response = await self._inline_continue_on_read_skill(
                    full_response,
                    llm_response,
                    current_expr,
                    current_state,
                    conditions,
                    user_instruction,
                    active_skills,
                    loaded_skill_names,
                    loaded_skill_paths,
                    callback,
                    iteration,
                )

                # LLM 响应日志
                self.logger.llm_response("命令: %s", llm_response.command)
                self.logger.info("   思考: %s", llm_response.thinking[:150] + "..." if len(llm_response.thinking) > 150 else llm_response.thinking)
                self.logger.info("   解释: %s", llm_response.explanation)

                # ============ 自动技能检测（可选） ============
                # 说明：默认关闭，避免无意中加载大量技能（并降低每轮开销）。
                # 如需开启，请在 config.SolverConfig.auto_load_mentioned_skills = True。
                if getattr(self.config, "auto_load_mentioned_skills", False):
                    from .skills import detect_skill_mentions, get_skill_loader

                    mentioned_skills = detect_skill_mentions(
                        llm_response.thinking + " " + llm_response.explanation,
                        list(loaded_skill_paths)
                    )

                    if mentioned_skills:
                        loader = get_skill_loader()
                        for skill_info in mentioned_skills:
                            skill_path = skill_info["full_path"]
                            skill_rel_path = skill_info["path"]
                            skill_name = skill_info["name"]

                            # Skip if already loaded
                            if skill_name in loaded_skill_names or skill_path in loaded_skill_paths:
                                continue

                            try:
                                content_obj = loader.load_skill_content(skill_name)
                                if content_obj:
                                    skill_body = content_obj.full_content
                                else:
                                    with open(skill_path, 'r', encoding='utf-8-sig') as f:
                                        skill_body = f.read()

                                # 醒目的自动加载日志
                                print("")
                                print_separator("═", 60)
                                self.logger.skill_load("🔮 自动检测并装备技能: %s", skill_name)
                                print_separator("─", 60)
                                self.logger.info("   📁 路径: %s", skill_path)
                                self.logger.info("   📄 大小: %d 字节", len(skill_body))
                                self.logger.info("   💡 原因: LLM 在 thinking 中提到了此技能")
                                print_separator("═", 60)
                                print("")

                                # 添加到上下文
                                skill_dir = os.path.dirname(skill_path)
                                skill_block = (
                                    f"\n## 技能文件 (自动加载): {skill_rel_path}\n"
                                    f"> Directory: {skill_dir}\n\n"
                                    f"{skill_body}\n\n---\n"
                                )
                                active_skills.append(skill_block)
                                loaded_skill_names.add(skill_name)
                                loaded_skill_paths.add(skill_path)

                                # If no command was generated, treat this as a successful read_skill step.
                                if not llm_response.command.strip():
                                    self.logger.info("   🔄 Auto-converting skill load of '%s' to explicit command step", skill_name)
                                    llm_response.command = f"read_skill {skill_rel_path}"
                                    llm_response.explanation = f"系统自动检测并加载技能: {skill_name}"
                                    break

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
                if current_sub_step_number > 0 and callback:
                    await callback(SolveEvent(
                        type=EventType.SUB_PLAN_STEP_COMPLETE,
                        content=str(current_sub_step_number),
                        step=iteration,
                        metadata={"sub_step_number": current_sub_step_number}
                    ))
                if current_step_number > 0 and callback:
                    await callback(SolveEvent(
                        type=EventType.PLAN_STEP_COMPLETE,
                        content=str(current_step_number),
                        step=iteration,
                        metadata={"step_number": current_step_number}
                    ))
                current_sub_plan = None
                last_sub_plan_str = ""
                current_sub_step_number = 0
                current_step_number = 0
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
                        from .skills import get_skill_loader

                        skill_dir = os.path.dirname(abs_path)
                        abs_norm = os.path.abspath(abs_path)

                        loader = get_skill_loader()

                        # Best-effort: map file path -> canonical skill name
                        meta = None
                        try:
                            for s in loader.discover_skills():
                                if not getattr(s, "path", None):
                                    continue
                                if os.path.normcase(os.path.abspath(s.path)) == os.path.normcase(abs_norm):
                                    meta = s
                                    break
                        except Exception:
                            meta = None

                        canonical_name = meta.name if meta and getattr(meta, "name", None) else os.path.basename(skill_dir)

                        already_loaded = canonical_name in loaded_skill_names or abs_norm in loaded_skill_paths

                        if already_loaded:
                            self.logger.info("   技能已加载过: %s (%s)", canonical_name, abs_norm)
                            self.logger.info("   (Treating as success since it is available)")
                            last_error = None
                        else:
                            # Strip YAML frontmatter to reduce token waste
                            def _strip_frontmatter(md: str) -> str:
                                lines = (md or "").splitlines()
                                if lines and lines[0].strip() == "---":
                                    for i in range(1, len(lines)):
                                        if lines[i].strip() == "---":
                                            return "\n".join(lines[i+1:]).strip()
                                return (md or "").strip()

                            # Prefer cached/parsed body via SkillLoader
                            content_obj = loader.load_skill_content(canonical_name) if meta else None
                            skill_body = content_obj.full_content if content_obj else _strip_frontmatter(file_content)

                            file_size = len(skill_body)

                            # 醒目的技能加载日志
                            print("")
                            print_separator("═", 60)
                            self.logger.skill_load("📚 装备技能: %s", canonical_name)
                            print_separator("─", 60)
                            self.logger.info("   📁 完整路径: %s", abs_norm)
                            self.logger.info("   📂 所在目录: %s", skill_dir)
                            self.logger.info("   📄 文件大小: %d 字节", file_size)
                            print_separator("═", 60)
                            print("")

                            skill_block = (
                                f"\n## 技能文件: {canonical_name}\n"
                                f"> Path: {abs_norm}\n"
                                f"> Directory: {skill_dir}\n\n"
                                f"{skill_body}\n\n---\n"
                            )
                            active_skills.append(skill_block)
                            loaded_skill_names.add(canonical_name)
                            loaded_skill_paths.add(abs_norm)

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

            # 收集 Token 信息
            token_usage = self.llm.last_token_usage or TokenUsage()
            total_token_usage.add(token_usage)

            # 创建轮次日志
            round_log = RoundLog(
                round=iteration,
                prompt_text=self.llm.last_prompt_text,
                response_text=full_response,
                state_before=current_expr or "",
                state_after=exec_result.result if exec_result.success else current_expr or "",
                token_usage=token_usage,
                prompt_components=self.llm.last_prompt_components.copy()
            )
            round_logs.append(round_log)

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

            # 生成本轮记忆
            try:
                memory_entry = await self._generate_memory(
                    response=full_response,
                    result=exec_result.result if exec_result.success else "",
                    error=exec_result.error if not exec_result.success else None
                )
                if memory_entry:
                    self.memory += f"\n### Round {iteration}\n\n{memory_entry}\n"
                    self.logger.info(f"记忆已更新: {memory_entry[:50]}...")
            except Exception as e:
                self.logger.warning(f"生成记忆失败: {e}")

            if exec_result.success:
                self.logger.result("成功 | 变更=%s", exec_result.changed)
                if exec_result.changed:
                    self.logger.info("   新表达式: %s", exec_result.result[:100] if exec_result.result else "None")
                self.error_handler.reset_consecutive_errors()
                last_error = None
                expression_history.append(exec_result.result)

                if callback:
                    # 明确告诉前端命令执行成功
                    # 记录调试日志：追踪为什么会有多余的 simplify
                    self.logger.info(f"发送 COMMAND_SUCCESS: {command} (changed={exec_result.changed}, step={iteration})")

                    await callback(SolveEvent(
                        type=EventType.COMMAND_SUCCESS,
                        content=command,
                        step=iteration,
                        metadata={
                            "explanation": llm_response.explanation,
                            "changed": exec_result.changed
                        }
                    ))
                    # 继续发送 RESULT 事件（包含结果数据）
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
                    # 2026-02-01 架构重构：明确告诉前端命令执行失败
                    await callback(SolveEvent(
                        type=EventType.COMMAND_FAILURE,
                        content=command,
                        step=iteration,
                        metadata={
                            "error": exec_result.error,
                            "explanation": llm_response.explanation
                        }
                    ))
                    # 继续发送 ERROR 事件（包含错误详情）
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
            iterations=iteration,
            total_token_usage=total_token_usage,
            round_logs=round_logs
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

    async def _generate_memory(self, response: str, result: str, error: Optional[str] = None) -> str:
        """生成本轮迭代的记忆摘要

        Args:
            response: LLM 的原始响应
            result: 执行结果
            error: 错误信息（如果有）

        Returns:
            记忆摘要字符串
        """
        # 构建记忆生成的提示词
        memory_prompt = f"""## 记忆更新指令

在这一步，我们需要为当前轮次的交互生成一个简洁的记忆摘要。

## 本轮交互信息

### LLM 的响应
{response}

### 执行结果
{result}

### 错误信息（如果有）
{error or '无'}

## 记忆生成要求

请生成一个简洁的记忆摘要（不超过100个token），包含以下内容：

1. **本轮尝试的方法**：简要描述使用了什么策略或命令
2. **执行结果**：成功还是失败
3. **失败原因**（如果失败）：记录错误信息和可能的原因
4. **经验教训**：避免在后续轮次中重复相同的错误
5. **替代方案**（如果有）：如果当前方法不可行，记录其他可能的方法

**注意**：
- 保持简洁，不超过100个token
- 重点记录**失败的尝试**，避免重复错误
- 不要使用 Markdown 标题（如 `## Updated Memory`）
- 直接输出记忆内容，不需要额外的格式
"""

        try:
            # 调用 API 生成记忆
            client = await self.llm._get_client()
            response_api = await client.post(
                f"{self.llm.config.api_base}/chat/completions",
                headers={
                    "Authorization": f"Bearer {self.llm.config.api_key}",
                    "Content-Type": "application/json"
                },
                json={
                    "model": self.llm.config.model,
                    "messages": [
                        {"role": "user", "content": memory_prompt}
                    ],
                    "max_tokens": 150,
                    "temperature": 0.7
                },
                timeout=30.0
            )

            if response_api.status_code != 200:
                raise Exception(f"API错误 {response_api.status_code}: {response_api.text}")

            response_data = response_api.json()
            memory_text = response_data["choices"][0]["message"]["content"]
            return memory_text.strip()
        except Exception as e:
            self.logger.warning(f"生成记忆失败: {e}")
            # 如果生成失败，返回简单的默认记忆
            if error:
                return f"本轮尝试失败，错误：{error[:50]}"
            else:
                return f"本轮执行成功，结果：{result[:50]}"

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
