"""Error Handler 模块

负责处理和分类错误，提供错误恢复建议。
"""

import logging
import traceback
from dataclasses import dataclass, field
from datetime import datetime
from typing import Optional, List, Dict, Any
from enum import Enum


class ErrorCategory(Enum):
    """错误分类"""
    SYNTAX = "syntax"
    RULE = "rule"
    SYSTEM = "system"
    API = "api"
    TIMEOUT = "timeout"
    UNKNOWN = "unknown"


@dataclass
class ErrorInfo:
    """错误信息结构"""
    category: ErrorCategory
    message: str
    suggestion: str
    recoverable: bool
    original_exception: Optional[Exception] = None
    context: Dict[str, Any] = field(default_factory=dict)
    timestamp: datetime = field(default_factory=datetime.now)


class ErrorHandler:
    """错误处理器"""
    
    ERROR_PATTERNS = {
        ErrorCategory.SYNTAX: [
            "ParseException", "SyntaxError", "解析错误", "语法错误", "命令语法错误",
        ],
        ErrorCategory.RULE: [
            "RuleException", "AssertionError", "规则应用错误", "cannot apply", "rule failed",
        ],
        ErrorCategory.API: [
            "APIError", "RateLimitError", "AuthenticationError", "ConnectionError", "HTTPError",
        ],
        ErrorCategory.TIMEOUT: [
            "TimeoutError", "Timeout", "timed out",
        ],
    }
    
    RECOVERY_SUGGESTIONS = {
        ErrorCategory.SYNTAX: "命令语法错误。不要猜测命令格式：请先用 <|load_skill|>search: 关键词<|end_load_skill|> 加载对应技能（或 read_skill 加载具体 SKILL.md），再严格按技能文档给出正确的Iscalc可执行命令。",
        ErrorCategory.RULE: "该规则无法应用于当前表达式。请尝试其他化简方法。",
        ErrorCategory.API: "API调用失败。系统将自动重试。",
        ErrorCategory.TIMEOUT: "操作超时。请尝试简化表达式或分步执行。",
        ErrorCategory.SYSTEM: "系统内部错误。请检查输入表达式是否有效。",
        ErrorCategory.UNKNOWN: "发生未知错误。请检查输入并重试。",
    }
    
    def __init__(self, max_consecutive_errors: int = 3):
        self.error_log: List[ErrorInfo] = []
        self.consecutive_errors: int = 0
        self.max_consecutive_errors = max_consecutive_errors
        self.logger = logging.getLogger(__name__)
    
    def handle(self, exception: Exception, context: Optional[Dict[str, Any]] = None) -> ErrorInfo:
        """处理异常"""
        context = context or {}
        category = self.categorize(exception)
        message = self._format_error_message(exception)
        suggestion = self._get_suggestion(category)
        recoverable = self._is_recoverable(category)
        
        error_info = ErrorInfo(
            category=category,
            message=message,
            suggestion=suggestion,
            recoverable=recoverable,
            original_exception=exception,
            context=context
        )
        
        self._log_error(error_info)
        self.consecutive_errors += 1
        
        return error_info
    
    def categorize(self, exception: Exception) -> ErrorCategory:
        """分类错误类型"""
        exception_name = type(exception).__name__
        exception_str = str(exception)
        
        for category, patterns in self.ERROR_PATTERNS.items():
            for pattern in patterns:
                if pattern in exception_name or pattern in exception_str:
                    return category
        
        return ErrorCategory.UNKNOWN
    
    def format_for_llm(self, error_info: ErrorInfo) -> str:
        """格式化错误信息供LLM理解"""
        parts = [
            f"错误类型: {error_info.category.value}",
            f"错误信息: {error_info.message}",
            f"建议: {error_info.suggestion}",
        ]
        
        if error_info.context:
            if "command" in error_info.context:
                parts.append(f"失败的命令: {error_info.context['command']}")
            if "expression" in error_info.context:
                parts.append(f"当前表达式: {error_info.context['expression']}")
        
        if self.consecutive_errors >= self.max_consecutive_errors:
            parts.append(f"\n警告: 已连续发生{self.consecutive_errors}次错误，请考虑尝试完全不同的方法。")
        
        return "\n".join(parts)
    
    def reset_consecutive_errors(self) -> None:
        """重置连续错误计数"""
        self.consecutive_errors = 0
    
    def should_suggest_check_input(self) -> bool:
        """是否应该建议用户检查输入"""
        return self.consecutive_errors >= self.max_consecutive_errors
    
    def get_error_log(self) -> List[ErrorInfo]:
        return self.error_log.copy()
    
    def clear_error_log(self) -> None:
        self.error_log = []
        self.consecutive_errors = 0
    
    def _format_error_message(self, exception: Exception) -> str:
        return f"{type(exception).__name__}: {str(exception)}"
    
    def _get_suggestion(self, category: ErrorCategory) -> str:
        return self.RECOVERY_SUGGESTIONS.get(category, self.RECOVERY_SUGGESTIONS[ErrorCategory.UNKNOWN])
    
    def _is_recoverable(self, category: ErrorCategory) -> bool:
        return category in [ErrorCategory.SYNTAX, ErrorCategory.RULE, ErrorCategory.API]
    
    def _log_error(self, error_info: ErrorInfo) -> None:
        self.error_log.append(error_info)
        self.logger.error(f"[{error_info.category.value}] {error_info.message}")
