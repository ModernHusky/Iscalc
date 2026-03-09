---
name: subgoal-definition
description: 在复杂证明中使用 define + subgoal + from/apply 进行分解与汇总。适用于参数积分、Skolem常数求解、围道积分与极限拼接等多阶段证明。
match_rules:
- (?i)subgoal
- (?i)from\s+\w+:
- (?i)define
- (?i)apply\s+\w+\s+on
---

# subgoal-definition
> 在复杂证明中使用 define + subgoal + from/apply 进行分解与汇总。

## 核心语法（与 parser 一致）

- `subgoal <int_id>: <expr>`
- `subgoal <int_id>: <expr> for <cond1>, <cond2>`
- `from <int_id>:`
- `define <lhs = rhs>`
- `define <lhs = rhs> for <cond1>, <cond2>`
- `apply <int_id> on <expr>`

规则：`<int_id>` 必须是数字（`1`, `2`, `3`），不要使用字符串名称。

## 状态约束（必须遵守）

- `define` 只能在 `INITIAL` 或 `PROVE` 使用。
- `subgoal` 用于进入一个新的 `PROVE` 子目标。
- 子目标不能嵌套（在 subgoal 内再开 subgoal 会失败）。
- `from <int_id>:` 只能从一个已存在的 subgoal 开始。
- `apply <int_id> on <expr>` 要求 `<expr>` 在当前式子中能精确定位。

## `define` 的合法形式

`define` 的表达式必须是等式，左边必须是：
- 一个变量：`define I = INT x:[0,1]. f(x)`
- 一个函数应用（参数是变量）：`define C(t,r) = ...`

条件 `for ...` 只能使用左边参数中出现的变量。

## 复杂证明推荐骨架

1. 先 `define` 辅助对象（路径、参数函数、记号）。
2. 建立多个 `subgoal`：
   - 子目标A：主公式（如留数定理结果）。
   - 子目标B：误差项/弧线项极限。
   - 子目标C：把主公式与原目标拼接。
3. 每个子目标内部用 `lhs:` 或 `rhs:` 局部证明并 `done`。
4. 回主目标后 `apply <id> on ...` 逐个替换，最后 `simplify` + `done`。

## 围道积分模板（针对 `INT 1/(x^2+1)` 型）

```text
prove (INT x:[-oo,oo]. 1/(x^2+1)) = pi

define C(t,r) = (r*exp(i*pi*(1-t)))_(t:[0,1])
define L(t,r) = (r*(1-2*t))_(t:[0,1])

subgoal 1: (LIM {r->oo}. CINT z:com(C(t,r),L(t,r)). 1/(z^2+1)) = -pi
lhs:
    apply residue theorem
done

subgoal 2: (LIM {r->oo}. CINT z:C(t,r). 1/(z^2+1)) = 0
lhs:
    apply cintegral identity
    simplify
done

subgoal 3: (LIM {r->oo}. ((CINT z:com(C(t,r),L(t,r)). 1/(z^2+1)) - (CINT z:C(t,r). 1/(z^2+1)))) = -(INT x:[-oo,oo]. 1/(x^2+1))
lhs:
    rewrite (CINT z:com(C(t,r),L(t,r)). 1/(z^2+1)) to (CINT z:C(t,r). 1/(z^2+1)) + (CINT z:L(t,r). 1/(z^2+1))
    rewrite (CINT z:L(t,r). 1/(z^2+1)) to (-(INT x:[-r,r]. 1/(x^2+1)))
    simplify
done

lhs:
    apply 3 on (INT x:[-oo,oo]. 1/(x^2+1))
    apply 1 on (LIM {r->oo}. CINT z:com(C(t,r),L(t,r)). 1/(z^2+1))
    apply 2 on (LIM {r->oo}. CINT z:C(t,r). 1/(z^2+1))
    simplify
done
```

## 常见失败点

- 缺冒号：`subgoal 1 ...`（错误）→ `subgoal 1: ...`
- `define` 非等式：`define C(t,r)`（错误）
- 围道定义误写为 `define C(R) = ... where ...` 或 `define C(R) = com(...)`（错误）；应先分别定义单路径（如 `define C(t,R) = (...)_(t:[0,1])`、`define L(t,R) = (...)_(t:[0,1])`），再在积分中写 `CINT z:com(C(t,R),L(t,R)). ...`
- `from 1,2:`（错误，一次只能一个）
- `apply` 目标表达式不精确，导致 `source expression ... not found`
- 对 `com(path1,path2)` 直接 `apply cintegral identity`（通常无效，先拆分）


