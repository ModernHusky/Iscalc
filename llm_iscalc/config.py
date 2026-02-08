"""配置文件

定义系统的各种配置参数。
"""

import os
from dataclasses import dataclass, field
from typing import Optional
from pathlib import Path

# Import the .env file from project root
try:
    from dotenv import load_dotenv
    project_root = Path(__file__).parent.parent
    root_env_path = project_root / '.env'
    if root_env_path.exists():
        load_dotenv(dotenv_path=root_env_path)
    else:
        load_dotenv()  # 尝试从当前工作目录加载
except ImportError:
    pass


@dataclass
class LLMConfig:
    """LLM配置"""
    # # 智谱 GLM API 配置
    # api_key: str = "a4810096a77743b3883be03959158863.oLZbloLwun4SRHO2"
    # api_base: str = "https://open.bigmodel.cn/api/paas/v4/"  # 注意：只需要 base URL，不包括具体 endpoint
    # model: str = "GLM-4.7"  # 使用 Plus 版本以获得更强的数学推理能力
    # DeepSeek API 配置（备用）
    api_key: str = field(default_factory=lambda: os.getenv("API_KEY", ""))
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
    max_iterations: int = 50    # 最大迭代次数
    max_consecutive_errors: int = 3    # 最大连续错误次数
    timeout_seconds: float = 300.0    # 超时时间（秒）
    loop_detection_window: int = 5    # 循环检测窗口大小
    use_tool_calling: bool = False  # Search-o1 风格：流式生成 + 动态技能加载
    auto_load_mentioned_skills: bool = False  # 默认关闭：只在显式 <|load_skill|> / read_skill 时加载


@dataclass
class IscalcConfig:
    """Iscalc配置"""
    base_theory: str = "base"
    default_conditions: list = field(default_factory=list)


@dataclass
class Config:
    """系统总配置"""
    llm: LLMConfig = field(default_factory=LLMConfig)  # LLM 引擎配置（API、模型等）
    solver: SolverConfig = field(default_factory=SolverConfig)  # 求解器配置（迭代、超时等）
    iscalc: IscalcConfig = field(default_factory=IscalcConfig)  # IsCalc 系统配置（理论文件等）
    debug: bool = False  # 调试模式开关


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
    
    # 更新 LLM 配置
    if api_key:
        config.llm.api_key = api_key
    if api_base:
        config.llm.api_base = api_base
    if model:
        config.llm.model = model
    
    # 更新求解器配置
    if max_iterations:
        config.solver.max_iterations = max_iterations
    
    return config
