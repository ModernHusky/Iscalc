"""LLM 提示词定义

采用Agent Skills渐进式披露理念，模块化提示词构建。
优化技能选择策略，引导 LLM 主动加载所需技能。
"""

from typing import List, Dict, Any, Optional

# ============ 技能加载特殊标记（Search-o1 风格） ============
# LLM 可以在 thinking 中使用这些标记来触发技能加载
SKILL_LOAD_BEGIN = "<|load_skill|>"
SKILL_LOAD_END = "<|end_load_skill|>"

from .skills import get_skills_categorized_xml


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
    "command": "要执行的Iscalc命令，必须严格遵守skills中的命令书写格式",
    "explanation": "这个命令会做什么",
    "is_final": false
}

当你认为表达式已经是最简形式时，设置 is_final 为 true，此时 command 可以为空字符串。
"""


# ============ 策略指南（已迁移至 Skills） ============

ERROR_RECOVERY_GUIDE = """
## 错误恢复策略
如果命令执行失败:
1. 分析错误信息
2. 尝试不同的方法
3. 检查表达式语法是否正确
4. 考虑是否需要先进行其他变换
"""


# ============ 技能选择策略指南 ============

SKILL_SELECTION_GUIDE = """
## 🧠 技能系统（单轮次边思考边加载模式）

**核心机制**：你可以在**同一个思考过程中**加载技能并继续推理，无需等待下一轮。

### 工作流程（单轮次完成）

1. **分析问题**：识别当前表达式的类型和需要的操作
2. **加载技能**：输出 `<|load_skill|>技能名<|end_load_skill|>` 
3. **立即参考**：技能加载后**立即**查阅技能内容继续思考
4. **输出命令**：根据技能文档给出 JSON 格式命令

### ⭐ 关键示例（边思考边加载）

**正确示范**：
```
分析当前表达式：这是一个广义积分 INT x:[0,oo]. 1/(1+exp(a*x))，包含无穷上限，
我需要了解如何处理这类积分。让我查阅积分策略技能：

<|load_skill|>strategy-integral<|end_load_skill|>
[✓strategy-integral技能已加载]

参考加载的技能文档，对于广义积分（improper integral），需要先将无穷替换为变量 t，
然后取极限。技能文档中的命令是：`improper integral to limit creating t`

{"thinking": "这是广义积分，需要先替换无穷为变量 t，再取极限", 
 "command": "improper integral to limit creating t", 
 "explanation": "将积分上限的无穷替换为变量 t，转化为极限形式", 
 "is_final": false}
```

### 必须遵守的规则

- ✅ **在同一思考过程中完成**：加载技能后立即参考并输出命令
- ✅ **技能已加载标记后继续思考**：看到 `[✓xxx技能已加载]` 后，立即参考技能内容
- ✅ **严格按技能文档格式**：命令格式必须与技能示例完全一致
- ❌ **禁止中途停止**：加载技能后必须继续推理直到输出 JSON 命令
- ❌ **禁止猜测命令**：不确定语法时必须先加载对应技能

### 常用技能对照

| 问题类型 | 推荐技能 | 关键命令示例 |
|---------|---------|-------------|
| 广义积分 | `strategy-integral` | `improper integral to limit creating t` |
| 分部积分 | `integrate-by-parts` | `integrate by parts, u = ..., v = ...` |
| 换元积分 | `substitute` | `substitute u=..., u_range=...` |
| 极限问题 | `strategy-limit` | `rewrite to limit at ...` |
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
    - 第二层：相关命令详情 + 策略等（通过match_rules正则自动加载）
    
    Args:
        expression: 当前表达式
        include_examples: 是否包含示例
        include_all_commands: 是否包含所有命令详情（fallback模式）
        current_state: 当前状态名称（如 'CALCULATE', 'PROVE', 'INDUCTION'）
    
    Returns:
        构建的系统提示词
    """
    parts = [BASE_PROMPT]
    
    # 0. 技能选择策略指南 (最重要，放在最前面)
    parts.append(SKILL_SELECTION_GUIDE)
    
    # 1. 第一层：按类别分组的技能列表 (XML格式)
    parts.append("## 可用技能库 (Available Skills)")
    parts.append("以下是所有可用技能的分类清单。技能内容**默认不加载**，你只能看到名称和描述。")
    parts.append("> **需要详细指令时**：输出 `<|load_skill|>技能名<|end_load_skill|>` 加载完整内容。")
    parts.append("")
    parts.append(get_skills_categorized_xml())
    
    # 2. 当前状态提示
    parts.append(f"""
## 当前求解状态

当前状态: **{current_state}**

> 建议: 如果你不熟悉 {current_state} 状态下的可用操作，请加载状态技能:
> `<|load_skill|>state-{current_state.lower()}<|end_load_skill|>`
""")
    
    # 错误恢复指南
    parts.append(ERROR_RECOVERY_GUIDE)
    
    return "\n".join(parts)


# ============ 保留原有模板的兼容性 ============

# 为了向后兼容，保留原有的 SYSTEM_PROMPT 访问方式，但避免在 import 阶段触发技能扫描。
_SYSTEM_PROMPT_CACHE: Optional[str] = None

def get_system_prompt() -> str:
    global _SYSTEM_PROMPT_CACHE
    if _SYSTEM_PROMPT_CACHE is None:
        _SYSTEM_PROMPT_CACHE = build_dynamic_system_prompt("", include_all_commands=True)
    return _SYSTEM_PROMPT_CACHE

# 兼容旧代码：不要在 import 时构建完整提示词
SYSTEM_PROMPT = ""  # Deprecated: use get_system_prompt()

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
