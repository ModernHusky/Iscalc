"""工具函数模块

提供表达式验证、转换等工具函数。
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__))))

from integral import parser


def validate_expression(expression: str) -> tuple[bool, str]:
    """验证表达式是否有效
    
    Args:
        expression: 表达式字符串
    
    Returns:
        (是否有效, 错误信息或解析结果)
    """
    try:
        parsed = parser.parse_expr(expression)
        return True, str(parsed)
    except Exception as e:
        return False, str(e)


def expression_to_latex(expression: str) -> tuple[bool, str]:
    """将表达式转换为LaTeX
    
    Args:
        expression: 表达式字符串
    
    Returns:
        (是否成功, LaTeX字符串或错误信息)
    """
    try:
        parsed = parser.parse_expr(expression)
        if hasattr(parsed, 'export_latex'):
            return True, parsed.export_latex()
        elif hasattr(parsed, 'latex'):
            return True, parsed.latex()
        else:
            return True, str(parsed)
    except Exception as e:
        return False, str(e)


def validate_command(command: str) -> tuple[bool, str]:
    """验证命令是否有效
    
    Args:
        command: 命令字符串
    
    Returns:
        (是否有效, 错误信息或命令类型)
    """
    try:
        action = parser.parse_action(command)
        return True, type(action).__name__
    except Exception as e:
        return False, str(e)


def format_step_latex(step_num: int, command: str, result_latex: str) -> str:
    """格式化单步为LaTeX
    
    Args:
        step_num: 步骤编号
        command: 执行的命令
        result_latex: 结果的LaTeX表示
    
    Returns:
        格式化的LaTeX字符串
    """
    cmd_escaped = command.replace("_", "\\_").replace("&", "\\&")
    return f"\\text{{Step {step_num}: {cmd_escaped}}} \\\\\n&= {result_latex}"


def format_solution_markdown(steps: list, initial_expr: str, final_expr: str) -> str:
    """格式化完整求解过程为Markdown
    
    Args:
        steps: 步骤列表
        initial_expr: 初始表达式
        final_expr: 最终表达式
    
    Returns:
        Markdown格式的求解过程
    """
    lines = [
        "## 求解过程",
        "",
        f"**初始表达式:** `{initial_expr}`",
        "",
        "### 步骤",
        ""
    ]
    
    for i, step in enumerate(steps, 1):
        if step.get("success"):
            lines.append(f"{i}. **命令:** `{step['command']}`")
            if step.get("explanation"):
                lines.append(f"   - 说明: {step['explanation']}")
            lines.append(f"   - 结果: `{step.get('expr_after', '')}`")
            if step.get("latex_after"):
                lines.append(f"   - LaTeX: ${step['latex_after']}$")
        else:
            lines.append(f"{i}. **命令:** `{step['command']}` ❌")
            lines.append(f"   - 错误: {step.get('error', '未知错误')}")
        lines.append("")
    
    lines.extend([
        "### 结果",
        "",
        f"**最终表达式:** `{final_expr}`"
    ])
    
    return "\n".join(lines)
