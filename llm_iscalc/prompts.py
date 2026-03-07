"""LLM 提示词定义

采用Agent Skills渐进式披露理念，模块化提示词构建。
优化技能选择策略，引导 LLM 主动加载所需技能。
"""

from typing import List, Dict, Any, Optional

# ============ 技能加载特殊标记（Search-o1 风格） ============
# LLM 可以在 thinking 中使用这些标记来触发技能加载
SKILL_LOAD_BEGIN = "<|load_skill|>"
SKILL_LOAD_END = "<|end_load_skill|>"

# ============ 动作执行特殊标记（ReAct 风格） ============
# LLM 可以在 thinking 中使用这些标记来立即执行命令并获取反馈
EXECUTE_BEGIN = "<|execute|>"
EXECUTE_END = "<|end_execute|>"

# ============ 计划与步骤标记（Plan-and-Execute 风格） ============
PLAN_BEGIN = "<|plan|>"
PLAN_END = "<|end_plan|>"
PLAN_STEP_BEGIN = "<|step|>"
PLAN_STEP_END = "<|end_step|>"

# ============ 层次化子计划标记 ============
SUB_PLAN_BEGIN = "<|sub_plan|>"
SUB_PLAN_END = "<|end_sub_plan|>"
SUB_STEP_BEGIN = "<|sub_step|>"
SUB_STEP_END = "<|end_sub_step|>"

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

### 输出格式

你必须以JSON格式输出，包含以下字段:
{
    "thinking": "你的分析和推理过程。必须明确说明使用了哪个技能文件中的命令或策略（例如：'参考 skills/states/prove/SKILL.md' 或 '使用 rewrite-goal-proof 技能'）。这有助于用户进行调试。",
    "command": "要执行的Iscalc命令，必须严格遵守skills中的命令书写格式",
    "explanation": "这个命令会做什么",
    "is_final": false
}

当你认为表达式已经是最简形式，并且**当前的计划的所有剩余步骤都已经通过 `<|execute|>` 标签执行完毕**时，才能设置 is_final 为 true，此时 command 可以为空字符串。绝不要在计划未全部执行完成时结束执行！
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
## 🧠 层次化双层 Plan-and-Execute 系统

**核心机制**：这允许你在**同一个思考连贯过程中**制定宏观外层计划，并在执行每一项外层任务时，进一步制定和执行具体的**内层子计划**！你的任务是一次性解决整个问题（规划宏观路径，再逐步拆解攻克）。

### 工作流程：层次化双层执行

你的工作分层次耦合，将在这**一次持续的回答中**连贯完成：

#### 阶段一：初遇问题，深入分析后制定外层宏观计划
当看到表达式时，你**必须先用自然语言深入分析问题**：
- 这个表达式的结构是什么？（类型、特征、难点）
- 有哪些可能的求解/证明方向？
- 哪个方向最合适？为什么？

分析完成后，再输出宏观步骤数组（使用 `<|plan|>` 标记）：
<|plan|>
[
  {"step": 1, "description": "将广义积分转换为极限形式"},
  {"step": 2, "description": "换元处理"},
  {"step": 3, "description": "部分分式分解，化简并求值"}
]
<|end_plan|>

#### 阶段二：攻克具体步骤，制定子计划并执行
使用 `<|step|>` 标记宣示你正在攻击哪个外层步骤。
紧接着，针对当前这一个步骤，制定其所需的具体操作 **内层子计划** （使用 `<|sub_plan|>` 标记）。
随后，使用 `<|sub_step|>` 逐一执行这个子计划中的项。

<|step|>1<|end_step|>
思考第一步如何实现...我们需要先加载积分策略技能，然后应用极限转换。
<|sub_plan|>
[
  {"step": 1, "description": "加载 strategy-integral 技能"},
  {"step": 2, "description": "执行 improper integral to limit 创建极限变量"}
]
<|end_sub_plan|>

<|sub_step|>1<|end_sub_step|>
<|load_skill|>strategy-integral<|end_load_skill|>
[✓ 技能已加载]

<|sub_step|>2<|end_sub_step|>
<|execute|>improper integral to limit creating t<|end_execute|>
[✓ 执行结果: LIM {t -> oo}. INT ...]

<|step|>2<|end_step|>
进入宏观第二步，思考子计划...
<|sub_plan|>
[
  {"step": 1, "description": "使用 substitute 换元"}
]
<|end_sub_plan|>

<|sub_step|>1<|end_sub_step|>
<|execute|>substitute u for exp(a*x)<|end_execute|>
[✓ 执行结果: ...]

#### 计划可变性 (动态更新)
- **修改外层计划**：如果你发现宏观思路不通，重新输出 `<|plan|>` ，但**只包含尚未执行的剩余宏观步骤**。然后在后续使用对应的 `<|step|>` 标记。
- **修改内层子计划**：如果你在某个外层步骤内部遇到了错误需要调整，在当前 `<|step|>` 内重新输出 `<|sub_plan|>`，必须**只包含剩余未执行的子项**。然后接着使用新的 `<|sub_step|>`。

### 规则
1. **严格分层**：无论是 `<|plan|>` 还是 `<|sub_plan|>`，制定或修改时绝不能包含 `<|execute|>`。先计划，再执行。
2. **无限次动作**：每次 `execute` 获取到结果后，只要不能达到最终目的，请【立刻开启下一项 `<|sub_step|>`】或【进入下一个宏观 `<|step|>`】。你拥有在这个回合持续行动到底的权利！
3. **不能提前结束**：除非你确定外层计划的所有步骤、以及最后一个子计划全部彻底运行成功，**否则绝不能直接输出 JSON 结束（不能输出 `"is_final": true` 或 false）**。必须不断用 `<|execute|>` 和 `sub_step` 推进直到完成证明或化简任务。
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
