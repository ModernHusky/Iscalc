---
name: state-calculate
description: 计算状态 - 对表达式应用计算规则进行变换
match_rules:
  - ".*"
applicable_types:
  - general
---

# 计算状态 (Calculate State)

在计算状态下，可以执行以下操作：

## 可用命令

- **使用计算规则**: `rewrite`, `simplify`, `substitute`, `integrate by parts` 等
- **完成计算**: `done` - 当表达式已达到闭合形式时使用
- **切换到右边**: `rhs:` - 在证明等式时，完成左边计算后切换到右边

## 工作流程

1. 使用 `calculate <expr>` 或 `lhs:` 进入此状态
2. 应用各种计算规则进行变换
3. 达到目标后使用 `done` 退出

## 注意事项

- 只有在此状态下才能使用计算规则
- 不能在此状态下使用 `calculate`（已经在计算中）
- 不能在此状态下使用 `lhs:` 或 `arg:`（需要先 `done`）

# 完整示例 (来自 theories/)

## 示例: 标准计算流程
```
calculate INT x:[0, 1]. x*exp(x)
    integrate by parts with u = x, v = exp(x)
    apply integral identity
    simplify
done
```

## 示例: 多步骤计算
```
calculate INT x:[0, 1]. 2*x*arctan(x)
    integrate by parts with u = arctan(x), v = x^2
    simplify
    rewrite x ^ 2 / (x ^ 2 + 1) to 1 - 1 / (x ^ 2 + 1)
    apply integral identity
    simplify
done
```
