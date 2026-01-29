---
name: simplify
description: 自动化简 - 自动简化当前表达式。在大多数步骤后使用，处理线性性、导数简化等。
keywords:
  - 简化
  - 化简
  - simplify
applicable_types:
  - definite_integral
  - indefinite_integral
  - improper_integral
  - limit
  - derivative
  - summation
  - general
match_rules:
  - ".*"                    # 匹配所有表达式（simplify 是通用命令）
---

# 何时使用

在几乎所有操作后都应使用：
- 换元后化简积分
- 积分后化简结果
- 极限计算后简化
- 任何需要整理表达式的时候

# 快速开始

```
simplify
```

# 功能

- 线性性处理（拆分加法）
- 常数提取
- 代数化简
- 导数计算
- 极限求值
- 表达式规范化

# 示例

应用于任何复杂表达式：
```
simplify
```

# 注意事项

# 注意事项

- **不做分式化简**：`simplify` 不会自动处理分式通分或拆分。请使用 `rewrite` (如 `x/(x+1) -> 1 - 1/(x+1)`) 或 `partial fraction decomposition`。
- **不做积分计算**：它只化简表达式，不计算积分。请用 `apply integral identity`。
- **条件依赖**：某些化简（如 `abs(x) -> x`）依赖于变量条件（如 `x >= 0`）。如果发现未化简，可能是条件不足。
- **多次使用**：复杂的表达式可能需要多次 `simplify` 或与其他规则交替使用。

# 完整示例 (来自 theories/)

## 示例: 通用工作流程
`simplify` 几乎在每个步骤后都会使用：
```
calculate INT x:[0, 1]. x*exp(x)
    integrate by parts with u = x, v = exp(x)
    apply integral identity
    simplify              <-- 化简积分结果
done

calculate INT x:[0,pi]. x * sin(x) / (1 + cos(x) ^ 2)
    substitute y for pi - x
    expand polynomial
    simplify              <-- 换元后化简
    ...
done
```
