---
name: expand-polynomial
description: 展开多项式 - 将被积函数中的多项式展开为标准形式。
match_rules:
  - "(?i)int"
  - "\\^\\s*\\d"
applicable_types:
  - definite_integral
  - indefinite_integral
---

# 何时使用

当被积函数包含需要展开的多项式（如 $(x+1)^3$）时使用。

# 快速开始

```
expand polynomial
```

# 说明

- 此规则会展开**整个被积函数**中的多项式。
- 如果只需要展开某个子表达式，请使用 `rewrite` 手动指定。

# 示例

表达式: `INT x. (x + 1)^2`

```
expand polynomial
```

结果: `INT x. x^2 + 2*x + 1`

# 注意事项

- **仅适用于积分内的多项式**。
- 展开后通常可以直接使用 `apply integral identity`。
- 如果多项式过大，考虑使用换元法代替展开。
