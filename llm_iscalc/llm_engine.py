"""LLM Engine 模块

负责与LLM API交互，生成求解策略和命令。
支持详细的日志调试输出。
支持 Function Calling 让 LLM 实时读取技能文件。
"""

import json
import asyncio
import os
from typing import Optional, List, Dict, Any, AsyncGenerator, Callable
from dataclasses import dataclass

import httpx

from .logger_config import get_phase_logger, print_separator
from .prompts import SYSTEM_PROMPT, USER_MESSAGE_TEMPLATE, HISTORY_TEMPLATE, ERROR_FEEDBACK_TEMPLATE, build_dynamic_system_prompt
from .config import LLMConfig


# ============ Function Calling 工具定义 ============

SKILL_TOOLS = [
    {
        "type": "function",
        "function": {
            "name": "read_skill",
            "description": "读取技能文件(SKILL.md)以学习如何处理特定类型的问题。当你不确定如何使用某个命令或处理某类表达式时，使用此工具加载相关技能。",
            "parameters": {
                "type": "object",
                "properties": {
                    "skill_path": {
                        "type": "string",
                        "description": "技能文件的相对路径，如 'skills/strategies/strategy-integral/SKILL.md' 或 'skills/commands/rewrite/SKILL.md'"
                    }
                },
                "required": ["skill_path"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "list_skill_resources",
            "description": "列出技能的扩展资源。某些技能有 references/ 目录（包含详细文档）或 scripts/ 目录（包含辅助脚本）。使用此工具查看有哪些扩展资源可用。",
            "parameters": {
                "type": "object",
                "properties": {
                    "skill_name": {
                        "type": "string",
                        "description": "技能名称，如 'rewrite' 或 'strategy-integral'"
                    }
                },
                "required": ["skill_name"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "read_skill_resource",
            "description": "读取技能的扩展资源内容。当 SKILL.md 中提到需要参考更详细的文档或脚本时，使用此工具加载这些资源。",
            "parameters": {
                "type": "object",
                "properties": {
                    "skill_name": {
                        "type": "string",
                        "description": "技能名称"
                    },
                    "resource_name": {
                        "type": "string",
                        "description": "资源文件名，如 'ADVANCED.md' 或 'helper.py'"
                    }
                },
                "required": ["skill_name", "resource_name"]
            }
        }
    }
]


@dataclass
class LLMResponse:
    """LLM响应结构"""
    thinking: str
    command: str
    explanation: str
    is_final: bool
    raw_response: str


class LLMEngine:
    """LLM引擎"""
    
    def __init__(self, config: LLMConfig):
        self.config = config
        self.logger = get_phase_logger(__name__)
        self._client: Optional[httpx.AsyncClient] = None
    
    async def _get_client(self) -> httpx.AsyncClient:
        """获取HTTP客户端"""
        if self._client is None or self._client.is_closed:
            self._client = httpx.AsyncClient(timeout=self.config.timeout)
        return self._client
    
    async def close(self):
        """关闭客户端"""
        if self._client and not self._client.is_closed:
            await self._client.aclose()
    
    def _build_messages(
        self,
        expression: str,
        history: List[Dict[str, Any]],
        last_error: Optional[str] = None,
        use_dynamic_prompt: bool = True,
        current_state: str = "CALCULATE",
        conditions: Optional[List[str]] = None,
        user_instruction: Optional[str] = None,
        active_skills: Optional[List[str]] = None  # New: Loaded skills content (Layer 2)
    ) -> List[Dict[str, str]]:
        """构建消息列表
        
        Args:
            expression: 当前表达式
            history: 历史记录
            last_error: 上次错误
            use_dynamic_prompt: 是否使用渐进式披露的动态提示词
            current_state: 当前状态名称（如 'CALCULATE', 'PROVE', 'INDUCTION'）
            active_skills: 已加载的技能内容列表
        """
        # 根据配置选择提示词构建方式
        if use_dynamic_prompt:
            # 基础提示词构建
            base_system_prompt = build_dynamic_system_prompt(
                expression,
                include_examples=len(history) < 3,  # 前几步包含示例
                include_all_commands=False,
                current_state=current_state
            )
            
            # 注入已激活的技能 (Layer 2 Content)
            skill_section = ""
            if active_skills:
                skill_section = "\n\n# 已加载的技能详情 (Loaded Skills)\n" + "\n".join(active_skills)
            
            system_prompt = base_system_prompt + skill_section
            
        else:
            # 完整模式：加载所有命令
            system_prompt = SYSTEM_PROMPT
        
        messages = [{"role": "system", "content": system_prompt}]
        
        # 构建历史部分
        history_section = ""
        if history:
            steps = []
            for i, h in enumerate(history[-10:], 1):  # 只保留最近10步
                if h.get("success"):
                    steps.append(f"{i}. 命令: {h['command']}\n   结果: {h['expr_after']}")
                else:
                    steps.append(f"{i}. 命令: {h['command']}\n   错误: {h.get('error', '未知错误')}")
            history_section = HISTORY_TEMPLATE.format(steps="\n".join(steps))
        
        # 添加错误反馈
        if last_error:
            history_section += "\n" + ERROR_FEEDBACK_TEMPLATE.format(
                command=history[-1]["command"] if history else "unknown",
                error=last_error
            )
        
        user_message = USER_MESSAGE_TEMPLATE.format(
            expression=expression,
            history_section=history_section,
            current_state=current_state,

            conditions=", ".join(conditions) if conditions else "无",
            user_instruction=user_instruction or "无"
        )
        
        messages.append({"role": "user", "content": user_message})
        return messages
    
    async def generate_command(
        self,
        expression: str,
        history: List[Dict[str, Any]],
        last_error: Optional[str] = None,
        current_state: str = "CALCULATE",
        conditions: Optional[List[str]] = None,
        user_instruction: Optional[str] = None,
        active_skills: Optional[List[str]] = None
    ) -> AsyncGenerator[str, None]:
        """流式生成命令"""
        messages = self._build_messages(
            expression, 
            history, 
            last_error, 
            current_state=current_state, 
            conditions=conditions, 
            user_instruction=user_instruction,
            active_skills=active_skills
        )
        
        # 请求日志
        system_len = len(messages[0]["content"])
        user_len = len(messages[1]["content"]) if len(messages) > 1 else 0
        
        self.logger.llm_thinking("发送请求到 LLM API")
        self.logger.info("   模型: %s", self.config.model)
        self.logger.info("   系统提示词: %d 字符", system_len)
        self.logger.info("   用户消息: %d 字符", user_len)
        self.logger.info("   已装备技能: %d 个", len(active_skills or []))
                
        for attempt in range(self.config.max_retries):
            try:
                client = await self._get_client()
                
                async with client.stream(
                    "POST",
                    f"{self.config.api_base}/chat/completions",
                    headers={
                        "Authorization": f"Bearer {self.config.api_key}",
                        "Content-Type": "application/json"
                    },
                    json={
                        "model": self.config.model,
                        "messages": messages,
                        "temperature": self.config.temperature,
                        "max_tokens": self.config.max_tokens,
                        "stream": True
                    }
                ) as response:
                    if response.status_code != 200:
                        error_text = await response.aread()
                        raise Exception(f"API错误 {response.status_code}: {error_text.decode()}")
                    
                    async for line in response.aiter_lines():
                        if line.startswith("data: "):
                            data = line[6:]
                            if data == "[DONE]":
                                break
                            try:
                                chunk = json.loads(data)
                                if "choices" in chunk and chunk["choices"]:
                                    delta = chunk["choices"][0].get("delta", {})
                                    content = delta.get("content", "")
                                    if content:
                                        yield content
                            except json.JSONDecodeError:
                                continue
                return
                
            except Exception as e:
                self.logger.warning(f"API调用失败 (尝试 {attempt + 1}/{self.config.max_retries}): {e}")
                if attempt < self.config.max_retries - 1:
                    await asyncio.sleep(self.config.retry_delay * (2 ** attempt))
                else:
                    raise
    
    async def continue_generation(
        self,
        base_messages: List[Dict[str, str]],
        previous_output: str,
        skill_result: str,
        skill_name: str
    ) -> AsyncGenerator[str, None]:
        """Search-o1 风格：继续生成
        
        当 LLM 输出中包含技能加载标记时，中断生成，加载技能后调用此方法继续。
        
        Args:
            base_messages: 原始的 system + user 消息
            previous_output: 之前生成的部分输出（到技能标记处）
            skill_result: 加载的技能内容
            skill_name: 技能名称
        
        Yields:
            继续生成的文本块
        """
        # 构建消息：原始消息 + assistant的部分输出 + 技能注入 + 继续提示
        messages = base_messages.copy()
        
        # 添加之前的输出作为 assistant 消息（表示 LLM 已经输出了这些内容）
        messages.append({
            "role": "assistant",
            "content": previous_output
        })
        
        # 添加技能结果作为 user 消息（模拟系统注入）
        continue_prompt = f"""[系统消息] 技能文档已加载: {skill_name}

技能内容：
{skill_result}

现在你已经拥有了这个技能的完整文档。请**立即恢复你的思考**，并结合这些新知识解决问题：

1. **继续分析**：
   - 既然已经知道了 `{skill_name}` 的具体用法，现在应该怎么做？
   - 验证之前的假设，并确定具体的参数。

2. **必须输出行动**：
   - **不要停止！** 请继续输出你的思考过程。
   - 最终必须以 **JSON 格式命令** 结束本次回答。
     格式：`{{"thinking": "...", "command": "...", "explanation": "...", "is_final": ...}}`

请紧接着上文继续输出（不要重复上文，直接继续思考）："""
        
        messages.append({
            "role": "user", 
            "content": continue_prompt
        })
        
        self.logger.info("🔄 继续生成 (技能: %s)", skill_name)
        
        # 流式生成続きの部分
        for attempt in range(self.config.max_retries):
            try:
                async with httpx.AsyncClient(timeout=60.0) as client:
                    async with client.stream(
                        "POST",
                        f"{self.config.api_base}/chat/completions",
                        headers={"Authorization": f"Bearer {self.config.api_key}"},
                        json={
                            "model": self.config.model,
                            "messages": messages,
                            "temperature": self.config.temperature,
                            "stream": True
                        }
                    ) as response:
                        async for line in response.aiter_lines():
                            if line.startswith("data: "):
                                data = line[6:]
                                if data == "[DONE]":
                                    break
                                try:
                                    chunk = json.loads(data)
                                    if chunk["choices"][0]["delta"].get("content"):
                                        content = chunk["choices"][0]["delta"]["content"]
                                        yield content
                                except json.JSONDecodeError:
                                    continue
                return
                
            except Exception as e:
                self.logger.warning(f"继续生成失败 (尝试 {attempt + 1}/{self.config.max_retries}): {e}")
                if attempt < self.config.max_retries - 1:
                    await asyncio.sleep(self.config.retry_delay * (2 ** attempt))
                else:
                    raise
    
    async def generate_command_sync(
        self,
        expression: str,
        history: List[Dict[str, Any]],
        last_error: Optional[str] = None,
        current_state: str = "CALCULATE",
        conditions: Optional[List[str]] = None,
        user_instruction: Optional[str] = None
    ) -> LLMResponse:
        """非流式生成命令（收集完整响应）"""
        full_response = ""
        async for chunk in self.generate_command(expression, history, last_error, current_state, conditions, user_instruction):
            full_response += chunk
        
        return self.parse_response(full_response)
    
    def parse_response(self, response: str) -> LLMResponse:
        """解析LLM响应"""
        import re
        
        # 策略1：尝试提取 JSON
        try:
            # 优先 1：找 Markdown 代码块
            json_block_match = re.search(r'```(?:json)?\s*\n?(\{.*?\})\s*\n?```', response, re.DOTALL)
            if json_block_match:
                json_str = json_block_match.group(1)
                data = json.loads(json_str)
                return LLMResponse(
                    thinking=data.get("thinking", ""),
                    command=data.get("command", ""),
                    explanation=data.get("explanation", ""),
                    is_final=data.get("is_final", False),
                    raw_response=response
                )
            
            # 优先 2：稳健的 JSON 提取（从后向前搜索）
            # 因为普通文本中可能包含 '{'（如数学公式），简单 find('{') 容易出错
            # 我们假设有效的 JSON 响应通常在输出的末尾
            
            # 找到最后一个 '}'
            end_idx = response.rfind("}")
            if end_idx != -1:
                # 从这一点向前扫描，寻找匹配的 '{'，使得 parse 成功
                # 为了效率，我们限制向前扫描的范围（例如最后 4000 个字符）
                scan_start = max(0, len(response) - 4000)
                subset = response[scan_start:end_idx+1]
                
                # 在 subset 中寻找所有 '{' 的位置
                start_indices = [m.start() for m in re.finditer(r'\{', subset)]
                
                # 从最靠后的 '{' 开始尝试解析
                for rel_start in reversed(start_indices):
                    candidate = subset[rel_start:]
                    try:
                        data = json.loads(candidate)
                        # 验证关键字段是否存在，避免解析到无关的 JSON
                        if "command" in data or "thinking" in data:
                            return LLMResponse(
                                thinking=data.get("thinking", ""),
                                command=data.get("command", ""),
                                explanation=data.get("explanation", ""),
                                is_final=data.get("is_final", False),
                                raw_response=response
                            )
                    except json.JSONDecodeError:
                        continue
            
        except Exception as e:
            self.logger.warning(f"JSON 解析失败: {e}")

            pass
        
        # 策略2：从 Markdown 格式中提取（如 ## 命令\n```\nxxx\n```）
        command_block = re.search(r'(?:命令|command)[：:]\s*[`"]?([^\n`"]+)[`"]?', response, re.IGNORECASE)
        if command_block:
            cmd = command_block.group(1).strip()
            thinking_match = re.search(r'(?:思考|thinking)[：:]\s*(.+?)(?=(?:命令|command|解释|explanation|$))', response, re.IGNORECASE | re.DOTALL)
            explanation_match = re.search(r'(?:解释|explanation)[：:]\s*(.+?)(?=$|\n\n)', response, re.IGNORECASE | re.DOTALL)
            
            return LLMResponse(
                thinking=thinking_match.group(1).strip() if thinking_match else "",
                command=cmd,
                explanation=explanation_match.group(1).strip() if explanation_match else "",
                is_final=False,
                raw_response=response
            )
        
        # 策略3：尝试识别常见命令模式
        common_commands = [
            r'(rewrite\s+.+)',
            r'(simplify(?:\s+.+)?)',
            r'(substitute\s+.+)',
            r'(apply\s+.+)',
            r'(integrate\s+.+)',
            r'(partial\s+fraction\s+decomposition)',
            r'(cat\s+.+)',
        ]
        
        for pattern in common_commands:
            match = re.search(pattern, response, re.IGNORECASE)
            if match:
                cmd = match.group(1).strip()
                self.logger.warning("从文本中提取命令: %s", cmd[:50])
                return LLMResponse(
                    thinking=response[:200] if len(response) > 200 else response,
                    command=cmd,
                    explanation="命令从响应文本中提取",
                    is_final=False,
                    raw_response=response
                )
        
        # 所有策略都失败
        self.logger.warning("无法从响应中解析命令，原始响应: %s", response[:300])
        return LLMResponse(
            thinking=response,
            command="",
            explanation="无法解析响应",
            is_final=False,
            raw_response=response
        )
    
    def _execute_read_skill(self, skill_path: str) -> str:
        """执行 read_skill 工具调用
        
        Args:
            skill_path: 技能文件的相对路径
            
        Returns:
            技能文件内容或错误信息
        """
        base_dir = os.path.dirname(os.path.abspath(__file__))
        
        # 标准化路径
        if os.path.isabs(skill_path):
            abs_path = skill_path
        else:
            abs_path = os.path.join(base_dir, skill_path)
        
        if not os.path.exists(abs_path):
            # 尝试从当前目录解析
            if os.path.exists(skill_path):
                abs_path = os.path.abspath(skill_path)
            else:
                return f"错误: 找不到技能文件 '{skill_path}'"
        
        try:
            with open(abs_path, 'r', encoding='utf-8') as f:
                content = f.read()
            
            # 醒目日志
            print("")
            print_separator("═", 60)
            self.logger.skill_load("🔮 LLM 实时调用读取技能: %s", skill_path)
            print_separator("─", 60)
            self.logger.info("   📁 路径: %s", abs_path)
            self.logger.info("   📄 大小: %d 字节", len(content))
            print_separator("═", 60)
            print("")
            
            return content
        except Exception as e:
            return f"错误: 读取技能文件失败 - {str(e)}"
    
    def _execute_list_skill_resources(self, skill_name: str) -> str:
        """执行 list_skill_resources 工具调用
        
        列出技能的扩展资源（第三层）。
        """
        from .skills import get_skill_loader
        
        loader = get_skill_loader()
        resources = loader.list_skill_resources(skill_name)
        
        if not resources:
            return f"技能 '{skill_name}' 没有扩展资源（references/ 或 scripts/ 目录为空）"
        
        result_lines = [f"技能 '{skill_name}' 的扩展资源:\n"]
        for res in resources:
            result_lines.append(f"- [{res.resource_type}] {res.name}")
        
        self.logger.info("   📂 列出扩展资源: %s 有 %d 个资源", skill_name, len(resources))
        return "\n".join(result_lines)
    
    def _execute_read_skill_resource(self, skill_name: str, resource_name: str) -> str:
        """执行 read_skill_resource 工具调用
        
        读取技能的扩展资源内容（第三层）。
        """
        from .skills import get_skill_loader
        
        loader = get_skill_loader()
        resource = loader.load_skill_resource(skill_name, resource_name)
        
        if not resource or not resource.content:
            return f"错误: 找不到扩展资源 '{skill_name}/{resource_name}'"
        
        # 醒目日志
        print("")
        print_separator("═", 60)
        self.logger.skill_load("📚 LLM 加载扩展资源: %s/%s", skill_name, resource_name)
        print_separator("─", 60)
        self.logger.info("   📁 路径: %s", resource.path)
        self.logger.info("   📄 类型: %s", resource.resource_type)
        self.logger.info("   📦 大小: %d 字节", len(resource.content))
        print_separator("═", 60)
        print("")
        
        return resource.content
    
    async def generate_with_tools(
        self,
        expression: str,
        history: List[Dict[str, Any]],
        last_error: Optional[str] = None,
        current_state: str = "CALCULATE",
        conditions: Optional[List[str]] = None,
        user_instruction: Optional[str] = None,
        active_skills: Optional[List[str]] = None,
        on_tool_call: Optional[Callable[[str, str], None]] = None,
        max_tool_calls: int = 5
    ) -> str:
        """支持 Function Calling 的非流式生成
        
        LLM 可以在思考过程中调用 read_skill 工具获取技能内容。
        
        Args:
            expression: 当前表达式
            history: 历史记录
            last_error: 上次错误
            current_state: 当前状态
            conditions: 条件列表
            user_instruction: 用户指导
            active_skills: 已加载的技能
            on_tool_call: 工具调用回调函数 (tool_name, result)
            max_tool_calls: 最大工具调用次数
            
        Returns:
            最终响应内容
        """
        messages = self._build_messages(
            expression,
            history,
            last_error,
            current_state=current_state,
            conditions=conditions,
            user_instruction=user_instruction,
            active_skills=active_skills
        )
        
        self.logger.llm_thinking("发送请求到 LLM API (支持工具调用)")
        self.logger.info("   模型: %s", self.config.model)
        self.logger.info("   工具: read_skill")
        
        tool_call_count = 0
        
        while tool_call_count < max_tool_calls:
            client = await self._get_client()
            
            # 发送请求 (非流式，以支持工具调用)
            response = await client.post(
                f"{self.config.api_base}/chat/completions",
                headers={
                    "Authorization": f"Bearer {self.config.api_key}",
                    "Content-Type": "application/json"
                },
                json={
                    "model": self.config.model,
                    "messages": messages,
                    "temperature": self.config.temperature,
                    "max_tokens": self.config.max_tokens,
                    "tools": SKILL_TOOLS,
                    "tool_choice": "auto"  # 让模型自己决定是否调用工具
                },
                timeout=self.config.timeout
            )
            
            if response.status_code != 200:
                raise Exception(f"API错误 {response.status_code}: {response.text}")
            
            result = response.json()
            choice = result["choices"][0]
            message = choice["message"]
            
            # 检查是否有工具调用
            if message.get("tool_calls"):
                tool_calls = message["tool_calls"]
                
                # 将 assistant 消息添加到历史
                messages.append(message)
                
                # 处理每个工具调用
                for tool_call in tool_calls:
                    tool_call_count += 1
                    func_name = tool_call["function"]["name"]
                    func_args = json.loads(tool_call["function"]["arguments"])
                    
                    self.logger.info("   🔧 工具调用 #%d: %s(%s)", 
                                    tool_call_count, func_name, func_args)
                    
                    # 执行工具
                    if func_name == "read_skill":
                        tool_result = self._execute_read_skill(func_args.get("skill_path", ""))
                    elif func_name == "list_skill_resources":
                        tool_result = self._execute_list_skill_resources(func_args.get("skill_name", ""))
                    elif func_name == "read_skill_resource":
                        tool_result = self._execute_read_skill_resource(
                            func_args.get("skill_name", ""),
                            func_args.get("resource_name", "")
                        )
                    else:
                        tool_result = f"未知工具: {func_name}"
                    
                    # 回调通知
                    if on_tool_call:
                        on_tool_call(func_name, tool_result[:200] + "..." if len(tool_result) > 200 else tool_result)
                    
                    # 将工具结果添加到消息
                    messages.append({
                        "role": "tool",
                        "tool_call_id": tool_call["id"],
                        "content": tool_result
                    })
                
                # 继续循环，让模型处理工具结果
                continue
            
            else:
                # 没有工具调用，返回最终内容
                return message.get("content", "")
        
        # 达到最大工具调用次数，返回最后的内容
        self.logger.warning("达到最大工具调用次数 (%d)", max_tool_calls)
        return message.get("content", "")
