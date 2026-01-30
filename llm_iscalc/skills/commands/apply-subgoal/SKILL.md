---
name: apply-subgoal
description: 应用已证明的子目标到当前表达式。
match_rules:
- (?i)apply
- (?i)subgoal
---

# apply-subgoal
> 应用已证明的子目标到当前表达式。
## 使用时机
当需要使用之前证明的 subgoal 来替换或化简当前表达式中的某部分时。
## 指令
### 快速开始

```
apply <subgoal_name> on <expr>
```

### 示例

假设已证明：
```
subgoal 1: I(a) = pi / (2 * a) for a > 0
```

当前表达式: `... I(1) ...`

```
apply 1 on I(1)
```

结果: `... pi / 2 ...`

### 注意事项

- subgoal 必须已经被证明（出现在当前步骤之前）。
- `<expr>` 必须是可以应用该 subgoal 的表达式（即匹配 subgoal 左边的模式）。
- 这与 `rewrite` 不同：`apply` 用于已证明的恒等式，`rewrite` 用于代数变换。
