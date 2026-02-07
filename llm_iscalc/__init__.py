"""LLM-Iscalc Integration System

一个将大语言模型(LLM)与Iscalc数学表达式求解引擎相结合的系统。
用户输入数学表达式后，LLM分析表达式并生成Iscalc命令，
系统执行命令并将结果反馈给LLM进行下一步规划，直到表达式化为最简形式。
"""

__version__ = "0.1.0"

from .command_executor import CommandExecutor, ExecutionResult
from .error_handler import ErrorHandler, ErrorInfo
from .llm_engine import LLMEngine
from .session_manager import SessionManager, Session
from .solver_loop import SolverLoop, SolveEvent, SolveResult

__all__ = [
    "CommandExecutor",
    "ExecutionResult",
    "ErrorHandler",
    "ErrorInfo",
    "LLMEngine",
    "SessionManager",
    "Session",
    "SolverLoop",
    "SolveEvent",
    "SolveResult",
]
