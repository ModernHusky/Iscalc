"""配置文件

定义系统的各种配置参数。
"""

import os
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class LLMConfig:
    """LLM配置"""
    # DeepSeek API 配置
    api_key: str = "sk-58fe61d83d944b91ba997b3d397289fd"
    api_base: str = "https://api.deepseek.com/v1"
    model: str = "deepseek-chat"
    temperature: float = 0.1
    max_tokens: int = 2000
    timeout: float = 60.0
    max_retries: int = 3
    retry_delay: float = 1.0


@dataclass
class SolverConfig:
    """求解器配置"""
    max_iterations: int = 50
    max_consecutive_errors: int = 3
    timeout_seconds: float = 300.0
    loop_detection_window: int = 5


@dataclass
class IscalcConfig:
    """Iscalc配置"""
    base_theory: str = "standard"
    default_conditions: list = field(default_factory=list)


@dataclass
class Config:
    """系统总配置"""
    llm: LLMConfig = field(default_factory=LLMConfig)
    solver: SolverConfig = field(default_factory=SolverConfig)
    iscalc: IscalcConfig = field(default_factory=IscalcConfig)
    debug: bool = False


default_config = Config()


def create_config(
    api_key: Optional[str] = None,
    api_base: Optional[str] = None,
    model: Optional[str] = None,
    max_iterations: Optional[int] = None,
    debug: bool = False
) -> Config:
    """创建配置实例"""
    config = Config(debug=debug)
    
    if api_key:
        config.llm.api_key = api_key
    if api_base:
        config.llm.api_base = api_base
    if model:
        config.llm.model = model
    if max_iterations:
        config.solver.max_iterations = max_iterations
    
    return config
