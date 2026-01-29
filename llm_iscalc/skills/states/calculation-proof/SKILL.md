---
name: state-calculation-proof
description: 计算证明 - 使用 lhs/rhs/arg 进行的证明
match_rules:
  - "="
  - "converges"
applicable_types:
  - general
---

# 计算证明 (Proof by Calculation)

使用 `lhs:`, `rhs:`, 或 `arg:` 开始的证明方式。

## 适用目标类型

- **等式**: `a = b`
- **非等式**: `a != b`
- **不等式**: `a < b`, `a <= b`, `a > b`, `a >= b`
- **收敛性**: `converges(<expr>)`

## 工作流程

1. 在 Prove 状态使用 `lhs:` 或 `rhs:` 进入计算
2. 对表达式应用变换规则
3. 使用 `done` 完成一边的计算
4. 如有需要，使用 `rhs:` 切换到另一边
5. 最终使用 `done` 完成证明

## 注意事项

- 目标必须是等式、非等式、不等式或收敛性谓词
- 对于收敛性目标，使用 `arg:` 而非 `lhs:`

# 完整示例 (来自 theories/)

## 示例: lhs/rhs 证明等式
```
prove (INT x:[0,oo]. log(1 + a^2 / x^2)) = a * pi for a > 0
lhs:
    integrate by parts with u = log(1+a^2/x^2), v = x
    simplify
    rewrite x^2 * (a^2 / x^2 + 1) to a^2 + x^2
    apply integral identity
    simplify
done
```

## 示例: 两边都需要计算
```
prove (INT x:[0,pi / 2]. sin(x) ^ 2 / (sin(x) + cos(x))) = ...
lhs:
    substitute y for pi / 2 - x
done
rhs:
    simplify
    rewrite 1 to sin(x) ^ 2 + cos(x) ^ 2
    ...
done
```
