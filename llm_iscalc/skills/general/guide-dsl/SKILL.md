---
name: guide-dsl
description: iscalc 系统 DSL 语言和交互协议指南
match_rules:
  - ".*"
---

# 核心概念

你正在使用一种特定领域语言 (DSL) 来表达和检查积分计算。系统的目标是通过多轮交互获得积分的计算结果。

## 证明状态

系统有五种状态：
1. **INITIAL**: 初始状态，无活跃目标。
2. **PROVE**: 当前有证明目标。
3. **CALCULATE**: 正在进行计算（变形）。**只有在此状态下才能应用 rewrite, simplify 等计算规则。**
4. **INDUCTION**: 归纳法开始。
5. **CASE**: 案例分析开始。

## 常用命令速查

### 开始计算/证明
- `calculate <expr>`: 开始计算表达式。
- `prove <expr>`: 开始证明表达式（通常用于证明相等 `a = b`）。

### 状态转换
- `lhs:`: 对当前证明目标的左边进行计算 (PROVE -> CALCULATE)。
- `rhs:`: 对当前证明目标的右边进行计算 (PROVE -> CALCULATE)。
- `arg:`: 对参数进行计算（用于证明不等式或收敛性）。
- `done`: 结束当前的计算或证明步骤 (CALCULATE -> INITIAL/PROVE)。

### 交互协议

你只需要输出**单个 Iscalc 命令字符串**。系统会自动解析并执行。

- **正确**: `rewrite log(x) to log(x) + 0`
- **正确**: `simplify`
- **错误**: `{ "lines": [...] }` (不要输出 JSON 结构的命令)

### 常用命令提示

- **化简**: `simplify`
- **重写**: `rewrite <old> to <new>`
- **积分**: `integrate by parts ...`, `substitute ...`
- **展开**: `expand polynomial`, `expand definition`

## 常见错误排查

1. **Syntax Error**:
   - 确保命令格式正确。
   - 不要输出 JSON 对象作为命令。
   - `calculate` 只能用于开始新的计算，**不能**在 `rewrite` 或 `simplify` 之后作为步骤使用。

2. **Rewriting Error**:
   - 如果提示找不到表达式，注意运算符优先级。尝试先重写外层结构添加括号。

3. **State Error**:
   - 不要在 `CALCULATE` 状态下使用 `lhs:` 或 `prove`。
   - 必须先 `done` 结束当前计算。

详细规则请参考各独立 Skill (rewrite, integrate-by-parts 等)。
