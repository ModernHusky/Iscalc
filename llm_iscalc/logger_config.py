"""统一日志配置模块

为 llm_iscalc 系统提供结构化日志支持。
支持阶段标记、彩色输出、工作流追踪。
"""

import logging
import sys
from typing import Optional
from datetime import datetime


# ============ 阶段标记常量 ============

class Phase:
    """工作流阶段标记"""
    SOLVER_START = "SOLVER START"
    DISCOVERY = "DISCOVERY"
    MATCHING = "MATCHING"
    LOADING = "LOADING"
    ITERATION = "ITERATION"
    LLM_THINKING = "LLM THINKING"
    LLM_RESPONSE = "LLM RESPONSE"
    EXECUTION = "EXECUTION"
    RESULT = "RESULT"
    SKILL_LOAD = "SKILL LOAD"
    COMPLETE = "COMPLETE"
    ERROR = "ERROR"


# ============ ANSI 颜色代码 ============

class Colors:
    """ANSI 颜色代码（支持 Windows 10+）"""
    RESET = "\033[0m"
    BOLD = "\033[1m"
    
    # 前景色
    RED = "\033[31m"
    GREEN = "\033[32m"
    YELLOW = "\033[33m"
    BLUE = "\033[34m"
    MAGENTA = "\033[35m"
    CYAN = "\033[36m"
    WHITE = "\033[37m"
    
    # 亮色
    BRIGHT_RED = "\033[91m"
    BRIGHT_GREEN = "\033[92m"
    BRIGHT_YELLOW = "\033[93m"
    BRIGHT_BLUE = "\033[94m"
    BRIGHT_MAGENTA = "\033[95m"
    BRIGHT_CYAN = "\033[96m"
    
    # 背景色
    BG_BLACK = "\033[40m"
    BG_BLUE = "\033[44m"


# 阶段颜色映射
PHASE_COLORS = {
    Phase.SOLVER_START: Colors.BRIGHT_GREEN + Colors.BOLD,
    Phase.DISCOVERY: Colors.BRIGHT_CYAN,
    Phase.MATCHING: Colors.BRIGHT_MAGENTA,
    Phase.LOADING: Colors.BRIGHT_BLUE,
    Phase.ITERATION: Colors.BRIGHT_YELLOW,
    Phase.LLM_THINKING: Colors.CYAN,
    Phase.LLM_RESPONSE: Colors.BRIGHT_CYAN,
    Phase.EXECUTION: Colors.BRIGHT_YELLOW + Colors.BOLD,
    Phase.RESULT: Colors.GREEN,
    Phase.SKILL_LOAD: Colors.BRIGHT_BLUE,
    Phase.COMPLETE: Colors.BRIGHT_GREEN + Colors.BOLD,
    Phase.ERROR: Colors.BRIGHT_RED + Colors.BOLD,
}

# 阶段图标映射
PHASE_ICONS = {
    Phase.SOLVER_START: "🚀",
    Phase.DISCOVERY: "📂",
    Phase.MATCHING: "🔍",
    Phase.LOADING: "📖",
    Phase.ITERATION: "🔄",
    Phase.LLM_THINKING: "🤖",
    Phase.LLM_RESPONSE: "💡",
    Phase.EXECUTION: "⚡",
    Phase.RESULT: "📊",
    Phase.SKILL_LOAD: "📦",
    Phase.COMPLETE: "🎉",
    Phase.ERROR: "❌",
}


# ============ 自定义日志格式化器 ============

class ColoredFormatter(logging.Formatter):
    """支持彩色输出和阶段标记的日志格式化器"""
    
    def __init__(self, use_colors: bool = True):
        super().__init__()
        self.use_colors = use_colors
        # 在 Windows 上启用 ANSI 支持
        if sys.platform == 'win32':
            self._enable_windows_ansi()
    
    def _enable_windows_ansi(self):
        """在 Windows 上启用 ANSI 转义序列支持"""
        try:
            import ctypes
            kernel32 = ctypes.windll.kernel32
            # 启用 ENABLE_VIRTUAL_TERMINAL_PROCESSING
            kernel32.SetConsoleMode(kernel32.GetStdHandle(-11), 7)
        except Exception:
            self.use_colors = False
    
    def format(self, record: logging.LogRecord) -> str:
        # 提取阶段信息（如果有）
        phase = getattr(record, 'phase', None)
        
        # 构建时间戳
        timestamp = datetime.now().strftime("%H:%M:%S.%f")[:-3]
        
        # 构建日志消息
        if phase:
            icon = PHASE_ICONS.get(phase, "")
            color = PHASE_COLORS.get(phase, "") if self.use_colors else ""
            reset = Colors.RESET if self.use_colors else ""
            
            # 格式: [时间] 图标 [阶段] 消息
            formatted = f"[{timestamp}] {icon} {color}[{phase}]{reset} {record.getMessage()}"
        else:
            # 普通日志格式
            level_color = ""
            if self.use_colors:
                if record.levelno >= logging.ERROR:
                    level_color = Colors.BRIGHT_RED
                elif record.levelno >= logging.WARNING:
                    level_color = Colors.BRIGHT_YELLOW
                elif record.levelno >= logging.INFO:
                    level_color = Colors.WHITE
                else:
                    level_color = Colors.CYAN
            
            reset = Colors.RESET if self.use_colors else ""
            formatted = f"[{timestamp}] {level_color}{record.levelname:8}{reset} {record.getMessage()}"
        
        # 添加异常信息
        if record.exc_info:
            formatted += "\n" + self.formatException(record.exc_info)
        
        return formatted


# ============ 阶段日志适配器 ============

class PhaseLoggerAdapter(logging.LoggerAdapter):
    """支持阶段标记的日志适配器"""
    
    def process(self, msg, kwargs):
        # 从 extra 中提取 phase
        extra = kwargs.get('extra', {})
        if 'phase' not in extra and hasattr(self, 'current_phase'):
            extra['phase'] = self.current_phase
            kwargs['extra'] = extra
        return msg, kwargs
    
    def phase(self, phase: str, msg: str, *args, **kwargs):
        """使用指定阶段记录日志"""
        kwargs.setdefault('extra', {})['phase'] = phase
        self.info(msg, *args, **kwargs)
    
    def solver_start(self, msg: str, *args, **kwargs):
        self.phase(Phase.SOLVER_START, msg, *args, **kwargs)
    
    def discovery(self, msg: str, *args, **kwargs):
        self.phase(Phase.DISCOVERY, msg, *args, **kwargs)
    
    def matching(self, msg: str, *args, **kwargs):
        self.phase(Phase.MATCHING, msg, *args, **kwargs)
    
    def loading(self, msg: str, *args, **kwargs):
        self.phase(Phase.LOADING, msg, *args, **kwargs)
    
    def iteration(self, msg: str, *args, **kwargs):
        self.phase(Phase.ITERATION, msg, *args, **kwargs)
    
    def llm_thinking(self, msg: str, *args, **kwargs):
        self.phase(Phase.LLM_THINKING, msg, *args, **kwargs)
    
    def llm_response(self, msg: str, *args, **kwargs):
        self.phase(Phase.LLM_RESPONSE, msg, *args, **kwargs)
    
    def execution(self, msg: str, *args, **kwargs):
        self.phase(Phase.EXECUTION, msg, *args, **kwargs)
    
    def result(self, msg: str, *args, **kwargs):
        self.phase(Phase.RESULT, msg, *args, **kwargs)
    
    def skill_load(self, msg: str, *args, **kwargs):
        self.phase(Phase.SKILL_LOAD, msg, *args, **kwargs)
    
    def complete(self, msg: str, *args, **kwargs):
        self.phase(Phase.COMPLETE, msg, *args, **kwargs)
    
    def error_phase(self, msg: str, *args, **kwargs):
        self.phase(Phase.ERROR, msg, *args, **kwargs)


# ============ 全局配置函数 ============

_configured = False


def configure_logging(
    level: int = logging.INFO,
    use_colors: bool = True,
    log_file: Optional[str] = None
):
    """配置全局日志系统
    
    Args:
        level: 日志级别
        use_colors: 是否使用彩色输出
        log_file: 可选的日志文件路径
    """
    global _configured
    
    if _configured:
        return
    
    # 获取根日志器
    root_logger = logging.getLogger()
    root_logger.setLevel(level)
    
    # 移除现有处理器
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
    
    # 控制台处理器 - 强制每条日志后刷新
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    console_handler.setFormatter(ColoredFormatter(use_colors=use_colors))
    
    # 确保立即刷新输出
    class FlushingHandler(logging.StreamHandler):
        def emit(self, record):
            super().emit(record)
            self.flush()
    
    console_handler = FlushingHandler(sys.stdout)
    console_handler.setLevel(level)
    console_handler.setFormatter(ColoredFormatter(use_colors=use_colors))
    root_logger.addHandler(console_handler)
    
    # 文件处理器（可选）
    if log_file:
        file_handler = logging.FileHandler(log_file, encoding='utf-8')
        file_handler.setLevel(level)
        file_handler.setFormatter(ColoredFormatter(use_colors=False))
        root_logger.addHandler(file_handler)
    
    _configured = True


def get_phase_logger(name: str) -> PhaseLoggerAdapter:
    """获取支持阶段标记的日志器
    
    Args:
        name: 日志器名称（通常使用 __name__）
    
    Returns:
        PhaseLoggerAdapter 实例
    """
    # 确保日志已配置
    if not _configured:
        configure_logging()
    
    logger = logging.getLogger(name)
    return PhaseLoggerAdapter(logger, {})


# ============ 便捷函数 ============

def print_separator(char: str = "═", length: int = 60):
    """打印分隔线"""
    print(char * length)


def format_skill_info(skill_name: str, skill_path: str, size_bytes: int = 0) -> str:
    """格式化技能信息字符串"""
    size_str = f" ({size_bytes} bytes)" if size_bytes else ""
    return f"{skill_name}: {skill_path}{size_str}"


def format_command_result(command: str, success: bool, changed: bool, result: str = "") -> str:
    """格式化命令执行结果"""
    status = "✅ 成功" if success else "❌ 失败"
    change_str = ", 有变更" if changed else ", 无变更"
    result_str = f"\n   结果: {result[:100]}..." if result and len(result) > 100 else f"\n   结果: {result}" if result else ""
    return f"{status}{change_str}{result_str}"
