"""数据模型定义"""
from dataclasses import dataclass, field
from typing import List, Optional
from datetime import datetime


@dataclass
class TokenUsage:
    """Token 使用统计"""
    prompt_tokens: int = 0
    completion_tokens: int = 0
    total_tokens: int = 0
    prompt_cache_hit_tokens: int = 0
    prompt_cache_miss_tokens: int = 0
    cost: float = 0.0
    
    def add(self, other: 'TokenUsage'):
        """累加另一个 TokenUsage"""
        self.prompt_tokens += other.prompt_tokens
        self.completion_tokens += other.completion_tokens
        self.total_tokens += other.total_tokens
        self.prompt_cache_hit_tokens += other.prompt_cache_hit_tokens
        self.prompt_cache_miss_tokens += other.prompt_cache_miss_tokens
        self.cost += other.cost


@dataclass
class PromptComponent:
    """提示词成分"""
    name: str
    token_count: int


@dataclass
class RoundLog:
    """单轮交互日志"""
    round: int
    prompt_text: str
    response_text: str
    state_before: str
    state_after: str
    token_usage: TokenUsage
    prompt_components: List[PromptComponent]
    timestamp: str = field(default_factory=lambda: datetime.now().isoformat())
