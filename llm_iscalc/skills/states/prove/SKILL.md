---
name: state-prove
description: PROVE状态下的操作指南。用于组织子目标、进入lhs/rhs计算、使用from/apply汇总结果并完成证明。
match_rules:
- .*
---

# state-prove
> PROVE状态下的操作指南。

## 可用命令

- `lhs:` / `rhs:` / `arg:`
- `subgoal <id>: <expr> [for ...]`
- `from <id>:`
- `define <lhs = rhs> [for ...]`
- `induction on <var> [starting from <n>]`
- `case analysis on <cond_or_expr>`
- `done`

> **【⚠️严肃警告：必须先加载技能】**
> 虽然上面列出了可用操作（如 induction, case analysis 等），但在决定执行具体命令前，**绝对禁止**直接凭简略提示输出！
> 你**必须**先通过 `<|load_skill|>对应的技能名称<|end_load_skill|>` 加载专属技能文档，仔细阅读其详细语法和约束后，方可输出对应的 JSON 命令。

## 推荐决策顺序

1. 目标是简单等式：先 `lhs:` 或 `rhs:` 直接计算。
2. 目标含多个结构块（极限、积分、参数函数、常数项）：先拆 `subgoal`。
3. 需要复用结构（路径、记号、参数函数）：先 `define`。
4. 子目标完成后，用 `from <id>:` 或 `apply <id> on ...` 汇总。

## 关键约束

- 不要在 PROVE 里用 `calculate ...`，应使用 `lhs:` / `rhs:` / `arg:`。
- 子目标不能嵌套。
- `subgoal` 与 `from/apply` 的 id 必须是整数数字。
- `from <id>:` 一次只能引用一个 subgoal。
- `done` 前必须确保当前目标已闭合。

## 围道积分工作流（推荐）

当目标包含 `CINT` / `LIM` / `residue theorem`：

1. `define` 围道路径（如 `C(t,r)`, `L(t,r)`）。
2. 拆分 2~4 个 subgoal：
   - 留数定理主项；
   - 弧线项极限；
   - 与实轴积分关系；
   - （可选）常数或方向修正。
3. 子目标逐个 `done` 后，在主证明 `lhs:` 中连续 `apply` 汇总。
4. 最后 `simplify` 收束到目标。

## 例子（骨架）

```text
prove TARGET

define C(t,r) = ...
subgoal 1: A
lhs:
    ...
done
subgoal 2: B
lhs:
    ...
done
lhs:
    apply 1 on ...
    apply 2 on ...
    simplify
done
```


