---
name: proof-states
description: 证明状态指南 - INITIAL, PROVE, CALCULATE, INDUCTION, CASE 状态及其转换。
match_rules:
  - ".*"
applicable_types:
  - general
---

# 五种证明状态

1.  **INITIAL**: 初始状态，无活跃目标。
2.  **PROVE**: 当前有证明目标（如 `prove a = b`）。
3.  **CALCULATE**: 正在进行表达式变换。**这是唯一可以使用计算规则（如 simplify, rewrite）的状态。**
4.  **INDUCTION**: 归纳法开始，等待选择分支。
5.  **CASE**: 案例分析开始，等待选择分支。

# 状态转换命令

## INITIAL -> PROVE / CALCULATE
- `prove <expr>`: 开始证明。
- `calculate <expr>`: 开始计算。

## PROVE -> CALCULATE
- `lhs:`: 对等式目标的**左边**进行计算。
- `rhs:`: 对等式目标的**右边**进行计算。
- `arg:`: 对不等式或收敛性目标的**参数**进行计算。
- `from <name>:`: 从已证明的 subgoal 开始重写。

## PROVE -> INDUCTION
- `induction on <var>`: 开始对变量 var 进行归纳。
- `induction on <var> starting from <n>`: 指定归纳起点。

## PROVE -> CASE
- `case analysis on <bool_expr>`: 布尔表达式 (true/false 两个分支)。
- `case analysis on <num_expr>`: 数值表达式 (positive/zero/negative 三个分支)。

## INDUCTION -> PROVE
- `base:`: 进入基础情况。
- `induct:`: 进入归纳步骤。

## CASE -> PROVE
- 布尔: `case true:` / `case false:`
- 数值: `case positive:` / `case zero:` / `case negative:`

## CALCULATE -> (退出)
- `done`: 结束当前计算。
    - -> INITIAL (如果用 `calculate` 开始)
    - -> PROVE (如果用 `lhs:` 等开始)
    - -> INDUCTION / CASE (如果在某个分支内)
- `rhs:`: 从左边计算切换到右边（仅当目标是等式时）。

# 注意事项

- **不能在 PROVE 状态使用 `calculate`**。应使用 `lhs:` 或 `rhs:`。
- **不能在 CALCULATE 状态嵌套 `lhs:` 或 `arg:`**。
