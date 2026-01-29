"""LLM 提示词定义

采用Agent Skills渐进式披露理念，模块化提示词构建。
"""

from typing import List, Dict, Any, Optional

from .skills import (
    get_all_skill_metadata,
    get_relevant_skills,
    get_skill_details,
    get_state_skills,
    COMMAND_SKILLS,
)


# ============ 基础提示词（始终加载） ============

BASE_PROMPT = """
你是一个数学表达式化简专家，专门使用Iscalc工具来求解和化简数学表达式。

## 你的任务
分析用户给出的数学表达式，思考如何一步步化简它，然后生成Iscalc命令来执行化简操作。

## Iscalc表达式语法

### 基本元素
- 变量: x, y, z, u, v, t, n, m 等
- 常数: 1, 2, 3.14, pi, i (虚数单位), G (欧拉常数)
- 无穷: inf, oo (正无穷), -inf, -oo (负无穷)

### 运算符
- 加减乘除: +, -, *, /
- 幂运算: ^ (例如 x^2)
- 绝对值: |x|

### 函数
- 三角函数: sin(x), cos(x), tan(x), cot(x), sec(x), csc(x)
- 反三角函数: arcsin(x), arccos(x), arctan(x), arccot(x)
- 指数对数: exp(x), log(x), sqrt(x)
- 特殊函数: gamma(x), factorial(n), binom(n,k)

### 积分
- 定积分: INT x:[a,b]. f(x)  例如: INT x:[0,1]. x^2
- 不定积分: INT x. f(x)  例如: INT x. sin(x)
- 广义积分: INT x:[0,oo]. exp(-x)

### 极限
- 普通极限: LIM {x -> a}. f(x)
- 左极限: LIM {x -> a-}. f(x)
- 右极限: LIM {x -> a+}. f(x)

### 求和与求导
- 求和: SUM(n, a, b, f(n))  例如: SUM(n, 0, oo, 1/n^2)
- 求导: D x. f(x)  例如: D x. x^2

### 求值
- 代入求值: [f(x)]_x=a,b  表示 f(b) - f(a)

## 输出格式

你必须以JSON格式输出，包含以下字段:
{
    "thinking": "你的分析和推理过程。必须明确说明使用了哪个技能文件中的命令或策略（例如：'参考 skills/states/prove/SKILL.md' 或 '使用 rewrite-goal-proof 技能'）。这有助于用户进行调试。",
    "command": "要执行的Iscalc命令",
    "explanation": "这个命令会做什么",
    "is_final": false
}

当你认为表达式已经是最简形式时，设置 is_final 为 true，此时 command 可以为空字符串。
"""


# ============ 策略指南（已迁移至 Skills） ============

# ERROR_RECOVERY_GUIDE 保留在 Python 代码中，因为它不是基于表达式的技能，而是通用错误处理
ERROR_RECOVERY_GUIDE = """
## 错误恢复策略
如果命令执行失败:
1. 分析错误信息
2. 尝试不同的方法
3. 检查表达式语法是否正确
4. 考虑是否需要先进行其他变换
"""


# ============ 动态提示词构建 ============

def build_dynamic_system_prompt(
    expression: str,
    include_examples: bool = True,
    include_all_commands: bool = False,

    current_state: str = "CALCULATE",
    user_instruction: Optional[str] = None
) -> str:
    """构建动态系统提示词
    
    基于渐进式披露原则：
    - 第一层：基础提示词 + 所有命令摘要（始终加载）
    - 第二层：相关命令详情 + 策略（通过match_rules正则自动加载）
    - 第三层：状态相关技能（根据current_state加载）
    
    Args:
        expression: 当前表达式
        include_examples: 是否包含示例
        include_all_commands: 是否包含所有命令详情（fallback模式）
        current_state: 当前状态名称（如 'CALCULATE', 'PROVE', 'INDUCTION'）
    
    Returns:
        构建的系统提示词
    """
    parts = [BASE_PROMPT]
    
    # 1. 第一层：所有技能摘要
    parts.append(get_all_skill_metadata())
    
    # 2. 第二层：相关技能详情（含策略）
    if include_all_commands:
        # Fallback模式：加载所有命令
        # 注意：这里我们只加载 COMMAND_SKILLS 中的技能（兼容性）或者全部发现的技能
        from .skills import get_skill_loader
        all_skills = get_skill_loader().discover_skills()
        parts.append(get_skill_details(all_skills, include_examples=include_examples))
    else:
        # 渐进式披露：基于表达式匹配
        relevant_skills = get_relevant_skills(expression, user_instruction)
        if relevant_skills:
            parts.append(get_skill_details(relevant_skills, include_examples=include_examples))
    
    # 2.5. 状态相关技能（始终加载）
    state_skills = get_state_skills(current_state)
    if state_skills:
        parts.append("\n## 当前状态相关指南\n")
        parts.append(get_skill_details(state_skills, include_examples=include_examples))

    # 3. 示例 (按需加载 strategy-examples 技能)
    if include_examples:
        from .skills import get_skill_loader
        example_skill = get_skill_loader().load_skill_content("strategy-examples")
        if example_skill:
            parts.append(example_skill.full_content)
    
    # 错误恢复指南始终包含（或者也可以做成技能）
    parts.append(ERROR_RECOVERY_GUIDE)
    
    return "\n".join(parts)


# ============ 保留原有模板的兼容性 ============

# 为了向后兼容，保留原有的SYSTEM_PROMPT（完整版）
SYSTEM_PROMPT = build_dynamic_system_prompt("", include_all_commands=True)

USER_MESSAGE_TEMPLATE = """
当前状态: {current_state}
当前表达式: {expression}
附加条件: {conditions}
用户指导: {user_instruction}

{history_section}

请根据当前状态，分析表达式并生成下一个Iscalc命令来继续化简。
注意：不同状态下可用的命令不同，请参考技能指南中"当前状态相关指南"部分。
"""

HISTORY_TEMPLATE = """
历史步骤:
{steps}
"""

ERROR_FEEDBACK_TEMPLATE = """
上一个命令执行失败:
命令: {command}
错误: {error}

请分析错误原因，尝试其他方法。
"""
