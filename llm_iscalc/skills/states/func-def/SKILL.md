---
name: state-func-def
description: define 命令指南。用于在 INITIAL/PROVE 中引入函数、常量或围道路径定义，并在后续 expand/fold 或 residue/cintegral 规则中复用。
match_rules:
- (?i)define
---

# state-func-def
> define 命令指南。

## 语法

- `define <lhs = rhs>`
- `define <lhs = rhs> for <cond1>, <cond2>`

## 左侧合法形式

- 常量定义：`define I = INT x:[0,1]. f(x)`
- 函数定义：`define f(x) = x^2 + 1`
- 多参数函数定义：`define C(t,r) = (r*exp(i*pi*(1-t)))_(t:[0,1])`

## 围道路径 define 的精确格式（按 parser.py）

围道路径 define 必须写成“函数 = 单个参数化路径”：

- `define C(t,r) = (r*exp(i*pi*(1-t)))_(t:[0,1])`
- `define L(t,r) = (r*(1-2*t))_(t:[0,1])`

参数语义（用于围道积分）：

- `t` 是方向参数（路径参数），决定从 `t=0` 到 `t=1` 的走向。
- `r` 是尺度参数，控制半径大小或线段长度。
- 改变 `t` 的参数化形式会改变方向，可能影响实轴项前的符号。

必须遵守：

- `define` 右边若用于围道规则，应是一个 `CINTPath`（即 `(...)_(...:[...,...])`）。
- 不要在 `define` 中写 `where`（DSL 不支持）。
- 不要把集合/几何描述写进 `define`（如 `{z: ...}`、`[-R,R]` 片段对象）。
- 不要把 `com(...)` 写在 `define` 右边（`com` 只在 `CINT` 中用于组合路径）。
- 路径引用要用函数调用形式（如 `C(t,r)`），不要只写 `C`。
- 建议让左侧第一个参数与右侧路径参数一致（如都用 `t`），避免替换和求导错位。

## 约束（来自运行时）

- `define` 必须是等式。
- 函数参数必须是变量，且不能重复。
- `for` 条件只能引用左侧参数。
- `define` 不要放在 `CALCULATE` 状态中执行。

## 与后续规则的配合

- `expand definition for <name>`：展开定义。
- `fold definition for <name>`：折叠回定义。
- 对围道路径定义，`apply residue theorem` / `apply cintegral identity` 可解析路径引用。
- 组合路径在积分处写：`CINT z:com(C(t,r),L(t,r)). ...`

## 围道路径定义示例

```text
define C(t,r) = (r*exp(i*pi*(1-t)))_(t:[0,1])
define L(t,r) = (r*(1-2*t))_(t:[0,1])

# 之后可写:
# CINT z:C(t,r). 1/(z^2+1)
# CINT z:com(C(t,r),L(t,r)). 1/(z^2+1)
```

