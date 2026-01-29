---
name: integrate-by-equation
description: 通过解方程求积分（用于循环积分等情况）
match_rules:
  - "(?i)int"
applicable_types:
  - definite_integral
  - indefinite_integral
---

# 何时使用

当当前的积分在早期的计算步骤中已经出现过，并且可以通过解线性方程来获得积分值时使用。常见于循环积分（如 $\int e^x \sin(x) dx$）。

# 使用说明

命令格式：`solve integral <expr>`

# 示例

假设计算过程如下（循环积分例子）：

```
  INT x:[0,pi / 2]. cos(x) * exp(2 * x)
= ...
= -(4 * (INT x:[0,pi / 2]. cos(x) * exp(2 * x))) + exp(pi) - 2
```

此时，表达式中包含原始积分（带有系数 -4）。要解出该积分，使用：
`solve integral INT x:[0, pi/2]. cos(x) * exp(2*x)`

这将自动解方程并给出最终答案：`(exp(pi) - 2) / 5`。

# 注意事项

*   **必须确保**被求解的积分完全出现在早期的步骤中。
*   当前表达式必须是该积分的线性函数（例如 `a * I + b`，其中 `a != 1`）。
