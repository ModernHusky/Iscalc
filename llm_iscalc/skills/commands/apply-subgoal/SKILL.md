---
name: apply-subgoal
description: 使用已证明的子目标替换当前表达式。典型用于多 subgoal 汇总、常数回代与围道积分最终拼接。
match_rules:
- (?i)apply\s+\w+\s+on
- (?i)subgoal
---

# apply-subgoal
> 使用已证明的子目标替换当前表达式。

## 语法

```text
apply <subgoal_id> on <expr>
```

- `<subgoal_id>`：子目标数字编号（只能是数字，如 `1`、`2`）。
- `<expr>`：当前式中要替换的子表达式。

## 实际匹配逻辑

1. 在上下文中查找 `<subgoal_id>` 对应的等式。
2. 尝试双向匹配：
   - `<expr>` 匹配子目标左侧，则替换为右侧；
   - `<expr>` 匹配子目标右侧，则替换为左侧。
3. 检查子目标 `for ...` 条件在当前上下文是否成立。

## 使用建议

- 优先写精确、最小的 `<expr>`，减少匹配歧义。
- 对复杂项先 `simplify` 再 `apply`，通常更稳定。
- 需要链式替换时按依赖顺序连续 `apply 1 on ...`、`apply 2 on ...`。

## 常见报错

- `lemma ... not found`：数字 id 不存在或对应子目标未证明完成（仅支持数字 id）。
- `source expression ... not found`：`on` 后表达式与当前式不一致。
- `Applying the rule has no effect`：子目标可用，但当前位置不匹配或条件不满足。

## 围道积分汇总示例

```text
lhs:
    apply 3 on (INT x:[-oo,oo]. 1/(x^2+1))
    apply 1 on (LIM {r -> oo}. CINT z:com(C(t,r),L(t,r)). 1/(z^2+1))
    apply 2 on (LIM {r -> oo}. CINT z:C(t,r). 1/(z^2+1))
    simplify
```


