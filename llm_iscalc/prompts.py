"""LLM 提示词定义

采用Agent Skills渐进式披露理念，模块化提示词构建。
优化技能选择策略，引导 LLM 主动加载所需技能。
"""

from typing import Optional

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

### 求值
- 代入求值: [f(x)]_x=a,b  表示 f(b) - f(a)

## 关键语法与状态约束（必须遵守）

这些约束直接来自 parser/state 实现，违反后命令会失败：

1. `subgoal` 语法
   - `subgoal <id>: <expr>` 或 `subgoal <id>: <expr> for <cond1>, <cond2>`
   - `<id>` 只能使用数字（如 `1`, `2`, `3`），便于后续 `apply 1 on ...`

2. `define` 语法
   - `define <lhs = rhs>` 或 `define <lhs = rhs> for <conds>`
   - `define` 的表达式必须是等式；左边应是变量或函数调用（如 `C(t,r)`）
   - `define` 只能在 `INITIAL/PROVE` 使用，不要在 `CALCULATE` 中使用

3. `from` / `apply` 语法
   - `from <id>:`（一次只能从一个 subgoal 开始）
   - `apply <id> on <expr>`（`<id>` 只能是数字）
   - `apply` 要求 `<expr>` 在当前目标中可匹配到子表达式

4. `done` 规则
   - 只有当前证明/计算真正闭合时才使用 `done`，否则会触发 `CheckFinishedException`

## 围道积分专用策略（高优先级）

当表达式出现 `CINT` / `com(...)` / `residue theorem` 时，优先采用：

1. 先 `define` 路径（例如 `C(t,r)`、`L(t,r)`）
2. 用多个 `subgoal` 拆分证明（留数主项、弧线极限为0、与实轴积分关系）
3. 对 `CINT z:com(path1,path2). ...`，先 `rewrite` 成单路径和
4. 再对单路径 `CINT` 使用 `apply cintegral identity` 或 `apply residue theorem`

注意：`apply cintegral identity` 对多路径 `com(...)` 往往无效，必须先拆分。

## 输出格式

1. **自由思考与技能加载（纯文本）**：
   在输出任何 JSON 之前，你可以先用纯文本输出你的分析过程。如果你需要加载技能，请在这个纯文本阶段输出 `<|load_skill|>技能名<|end_load_skill|>`。
   系统会自动拦截并为你注入技能内容。

2. **最终决定（JSON）**：
   当技能加载完成（看到系统返回 `[✓xxx技能已加载]`）并且你做出了最终决定后，**在输出的最后面**提供一个单独的 JSON 代码块。JSON 格式如下：

```json
{
    "thinking": "简短总结你的最终决定和参考了哪个技能。",
    "command": "要执行的Iscalc命令，必须严格遵守技能文档中的格式",
    "explanation": "这个命令会做什么",
    "is_final": false
}
```

当你认为表达式已经是最简形式时，设置 is_final 为 true，此时 command 可以为空字符串。
**注意：绝对不要把 `<|load_skill|>` 放在 JSON 的字段中！必须在 JSON 代码块之前的纯文本中进行技能加载。**
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


# ============ subgoal/define强化指南 ============

SUBGOAL_DEFINE_GUIDE = """
## subgoal / define 优先工作流（PROVE状态）

复杂证明（尤其是围道积分）优先使用以下顺序：
1. `define` 路径或辅助函数（如 `define C(t,r)=...`, `define I(a)=...`）
2. `subgoal 1`, `subgoal 2`, ... 拆分可独立证明的事实
3. 每个 subgoal 内用 `lhs:` / `rhs:` 做局部计算并 `done`
4. 回到主目标后，使用数字 id：`apply <id> on <expr>` 或 `from <id>:` 汇总
5. 最后 `simplify` 并 `done`

轮次策略：单轮只输出一个可执行命令，不要一次输出多行脚本。
"""


# ============ 技能选择策略指南 ============

SKILL_SELECTION_GUIDE = """
## 🧠 层次化双层 Plan-and-Execute 系统

**核心机制**：这允许你在**同一个思考连贯过程中**制定宏观外层计划，并在执行每一项外层任务时，进一步制定和执行具体的**内层子计划**！你的任务是一次性解决整个问题（规划宏观路径，再逐步拆解攻克）。

### 工作流程：层次化双层执行

1. **纯文本分析**：在纯文本区域自由思考，识别当前表达式类型。
2. **加载技能**：在纯文本中输出 `<|load_skill|>技能名<|end_load_skill|>`。
3. **立即参考**：技能加载后（看到技能已加载的系统提示），继续在纯文本中根据技能内容推理。
4. **输出命令（JSON代码块）**：在得出最终结论后，在回答的最后面单独输出 JSON 格式命令。

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

- ✅ **在同一思考过程中完成**：加载技能后立即参考并输出命令
- ✅ **技能已加载标记后继续思考**：看到 `[✓xxx技能已加载]` 后，立即参考技能内容
- ✅ **严格按技能文档格式**：命令格式必须与技能示例完全一致
- ❌ **禁止中途停止**：加载技能后必须继续推理直到输出 JSON 命令
- ❌ **禁止直接使用未经加载的命令**：即使你在其他文档（如 state-calculate）的提示中看到了某个命令（如 substitute, apply integral identity），你也**必须先使用 `<|load_skill|>` 加载该命令的专属技能文档**，仔细阅读其详细语法和约束后，才能生成该命令。**绝对禁止**仅凭简略提示或记忆就直接输出命令！

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
    parts.append(SUBGOAL_DEFINE_GUIDE)

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

> 若当前状态为 PROVE 且目标复杂，请优先加载:
> `<|load_skill|>subgoal-definition<|end_load_skill|>`
> `<|load_skill|>state-prove<|end_load_skill|>`
> 若包含围道积分，再加载:
> `<|load_skill|>strategy-contour-subgoal<|end_load_skill|>`
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

