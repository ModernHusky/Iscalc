"""LLM Engine 模块

负责与LLM API交互，生成求解策略和命令。
"""

import json
import asyncio
import logging
from typing import Optional, List, Dict, Any, AsyncGenerator
from dataclasses import dataclass

import httpx

from .prompts import SYSTEM_PROMPT, USER_MESSAGE_TEMPLATE, HISTORY_TEMPLATE, ERROR_FEEDBACK_TEMPLATE
from .config import LLMConfig


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
        self.logger = logging.getLogger(__name__)
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
        last_error: Optional[str] = None
    ) -> List[Dict[str, str]]:
        """构建消息列表"""
        messages = [{"role": "system", "content": SYSTEM_PROMPT}]
        
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
            history_section=history_section
        )
        
        messages.append({"role": "user", "content": user_message})
        return messages
    
    async def generate_command(
        self,
        expression: str,
        history: List[Dict[str, Any]],
        last_error: Optional[str] = None
    ) -> AsyncGenerator[str, None]:
        """流式生成命令"""
        messages = self._build_messages(expression, history, last_error)
        
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
    
    async def generate_command_sync(
        self,
        expression: str,
        history: List[Dict[str, Any]],
        last_error: Optional[str] = None
    ) -> LLMResponse:
        """非流式生成命令（收集完整响应）"""
        full_response = ""
        async for chunk in self.generate_command(expression, history, last_error):
            full_response += chunk
        
        return self.parse_response(full_response)
    
    def parse_response(self, response: str) -> LLMResponse:
        """解析LLM响应"""
        try:
            # 尝试提取JSON
            json_start = response.find("{")
            json_end = response.rfind("}") + 1
            
            if json_start != -1 and json_end > json_start:
                json_str = response[json_start:json_end]
                data = json.loads(json_str)
                
                return LLMResponse(
                    thinking=data.get("thinking", ""),
                    command=data.get("command", ""),
                    explanation=data.get("explanation", ""),
                    is_final=data.get("is_final", False),
                    raw_response=response
                )
        except json.JSONDecodeError:
            pass
        
        # 如果JSON解析失败，尝试从文本中提取命令
        return LLMResponse(
            thinking=response,
            command="",
            explanation="无法解析响应",
            is_final=False,
            raw_response=response
        )
