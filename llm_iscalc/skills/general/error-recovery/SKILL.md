---
name: error-recovery
description: 常见错误诊断与修复指南。
match_rules:
- .*
---

# error-recovery
> 常见错误诊断与修复指南。
## 使用时机
当 iscalc 系统返回错误时，参考本技能诊断问题：
- **Parsing Error (解析错误)**：语法问题。
- **State Error (状态错误)**：状态转换问题（如在错误的状态下使用了某个命令）。
- **Rule Error (规则错误)**：规则应用失败（如参数不匹配）。
## 指令
### 常见解析错误 (Parsing Error)

1.  **条件格式错误**：
    - **错误**: `for x in [0,1]` 或 `for 0 <= x <= 1`
    - **错误**: `x:real` (不支持冒号类型声明，请使用 `x>=0` 或其他关系式隐含)
    - **正确**: `for x >= 0, x <= 1`

2.  **`from` 语法错误**：
    - **错误**: `from 2, 3:` (不支持多个 subgoal)
    - **正确**: `from 1:`（一次只能从一个 subgoal 开始）

3.  **`rewrite` 后加 `for`**：
    - **错误**: `rewrite a to b for x > 0`
    - **正确**: `rewrite a to b`（条件应在定义 subgoal 时声明）

4.  **`differentiate` 语法**：
    - **错误**: `differentiate both sides with respect to a`
    - **正确**: `differentiate both sides at a`

5.  **`subgoal` 格式**：
    - **错误**: `subgoal converges(SUM(...))`
    - **正确**: `subgoal 1: converges(SUM(...))` (需要 name/id 和冒号)

6.  **级数收敛**：
    - 系统不支持 `apply alternating series test` 等规则。收敛性由模式匹配自动判断。

7.  **系列求和变量约束**：
    - **错误**: `subgoal converges(SUM(n, 0, oo, ...)) for x >= 0, x <= 1, n >= 0, isInt(n)`
    - **正确**: `subgoal converges(SUM(n, 0, oo, ...)) for x >= 0, x <= 1`
    - 原因：求和变量 `n` 的类型和范围已由 SUM 的上下限隐含指定，无需额外约束。

8.  **无效的语句**：
    - `eval-at expression` 不是有效的 DSL 语句。

9.  **未实现的收敛规则**：
    - `apply Weierstrass M-test with M_n=...` 未实现。
    - 通过构造已知收敛级数的形式，让系统通过模式匹配判断收敛性。

### 常见状态错误 (State Error)

1.  **在 PROVE 状态使用 `calculate`**：
    - 不能在 ProveState 中使用 `calculate`。
    - 使用 `lhs:`, `rhs:`, 或 `arg:` 进入 CALCULATE 状态。

2.  **在 CALCULATE 状态使用 `lhs:` 或 `arg:`**：
    - 不能嵌套使用。先 `done` 结束当前计算块。

3.  **`arg` 嵌套**：
    - 不能在由 `arg:` 展开的证明中再次使用 `arg`。

### 常见规则错误 (Rule Error)

1.  **`u * dv does not equal body`**: 分部积分的 u, v 选择错误。检查参数。
2.  **`old expression not found`**: `rewrite` 的左侧表达式未精确匹配。检查括号和结合律。
3.  **`Applying the rule has no effect`**: 规则不适用于当前表达式。尝试先化简或重写。

### 关于 `done` 的使用

- 用于结束当前的计算或证明。
- **计算状态**: 当前表达式为闭合形式时使用。
- **证明状态**: 两边可通过化简证明相等时使用。
- **归纳/案例分析**: 所有分支完成后使用。

### CheckFinishedException (未完成异常)

当抛出 `CheckFinishedException` 时，表示证明或计算未完成：

1.  **计算证明未完成**：
    - 等式/不等式无法通过计算证明。
    - 检查是否需要更多步骤化简或重写。

2.  **缺少良构性条件 (Wellformedness Conditions)**：
    - 表达式依赖的条件未被证明。
    - 错误信息会显示缺少的条件，需要用 `subgoal` 证明这些条件。

3.  **计算未达到闭合形式**：
    - 结果中仍有未处理的积分、极限或级数。
    - 继续使用 `apply integral identity`, `simplify` 等规则直到闭合。

### 调试技巧

1.  **解析错误**：仔细检查语法，参考上述常见解析错误列表。
2.  **状态错误**：确认当前状态，避免在错误状态下使用命令。
3.  **规则错误**：检查参数是否精确匹配，尝试分步应用规则。
4.  **完成检查失败**：检查表达式是否真正化简完毕，是否有未证明的条件。
