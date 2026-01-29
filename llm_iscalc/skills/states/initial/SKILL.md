---
name: state-initial
description: 初始状态 - 开始新计算或证明的入口点
match_rules:
  - ".*"
applicable_types:
  - general
---

# 初始状态 (Initial State)

在初始状态下，可以执行以下操作：

## 可用命令

- **开始计算**: `calculate <expr>` - 对表达式进行求值/化简
- **开始证明**: `prove <expr>` - 证明一个等式或不等式
- **定义函数**: `define <name> = <expr>` - 创建新的函数定义

## 注意事项

- 初始状态是系统的默认起点
- 不能在此状态下直接使用计算规则（如 `simplify`, `rewrite`）
- 必须先进入 `calculate` 或 `prove` 状态才能开始变换

# 完整示例 (来自 theories/)

## 示例1: 开始计算
```
calculate INT x:[0, 1]. x*exp(x)
```
进入 Calculate 状态后可以使用 `integrate by parts`, `simplify` 等规则。

## 示例2: 开始证明
```
prove (INT x:[1,oo]. 1 / ((x+a)*sqrt(x-1))) = pi / sqrt(a+1) for a > -1
```
进入 Prove 状态后可以使用 `lhs:`, `rhs:`, `subgoal` 等命令。

## 示例3: 定义函数
```
define I(a) = INT x:[0,oo]. exp(-a*x) for a > 0
```
定义后可以在后续计算中使用 `I(a)`。
