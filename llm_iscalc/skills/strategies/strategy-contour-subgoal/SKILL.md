---
name: strategy-contour-subgoal
description: 围道积分求解策略。用于含 CINT/com/留数定理 的证明，指导使用 define 路径、subgoal 拆分、apply 汇总完成实积分求值。
match_rules:
- (?i)cint
- (?i)contour
- (?i)residue
- (?i)cintegral
- (?i)com\(
- (?i)围道
- (?i)留数
---

# strategy-contour-subgoal
> 围道积分求解策略（define + subgoal 主导）。

## 触发场景

- 目标是 `INT x:[-oo,oo]` 或 `INT x:[0,oo]`，且出现复变方法线索。
- 推导中包含 `CINT`、`com(path1,path2)`、`apply residue theorem`。
- 需要把围道积分极限与实积分联系起来。

## 路径 define 语法检查（先过这一关）

- `define C(t,r) = (r*exp(i*pi*(1-t)))_(t:[0,1])`
- `define L(t,r) = (r*(1-2*t))_(t:[0,1])`

参数语义（建议固定此约定）：

- `t` 是方向参数，用来定义围道方向（从 `t=0` 到 `t=1`）。
- `r` 是尺度参数，用来定义半径大小或路径长短。
- 做实轴改写时方向会影响符号，先确认 `L(t,r)` 的起点和终点方向。

并且：

- 不要写 `define C(R) = ... where ...`（`where` 不在 DSL 中）。
- 不要写 `define C(R) = com(...)`（`com` 不是单路径对象定义）。
- 不要写集合/几何描述（如 `{z: ...}`、`[-R,R]` 段对象）。
- 组合路径只在积分里写：`CINT z:com(C(t,r),L(t,r)). f(z)`。

## 标准工作流

1. `define` 路径：
   - 弧线/半圆路径（如 `C(t,r)`）
   - 实轴线段路径（如 `L(t,r)`）
2. 拆 `subgoal`：
   - `subgoal 1`: 闭合围道积分极限（通常用 `apply residue theorem`）
   - `subgoal 2`: 弧线积分极限为 0（常用 `apply cintegral identity` + `simplify`）
   - `subgoal 3`: 围道与实积分关系（`rewrite` 拆分 `com(...)`，再转换线段积分）
3. 主目标汇总：
   - 在 `lhs:` 连续 `apply 1 on ...`, `apply 2 on ...`, `apply 3 on ...`
   - `simplify` 后 `done`

## 关键实现细节（必须遵守）

- `apply residue theorem` 需要闭合路径；路径引用必须能由 `define` 解析。
- `apply cintegral identity` 对单路径 `CINT` 更稳定。
- 对 `CINT z:com(path1,path2). ...` 先做：
  - `rewrite (CINT ... com(...). f) to (CINT ... path1. f) + (CINT ... path2. f)`
- 把 `CINT z:L(t,r). f(z)` 改写到实积分时注意方向：
  - 路径反向会引入负号。

## 推荐命令骨架

```text
prove TARGET

define C(t,r) = ...
define L(t,r) = ...

subgoal 1: ...
lhs:
    apply residue theorem
done

subgoal 2: ...
lhs:
    apply cintegral identity
    simplify
done

subgoal 3: ...
lhs:
    rewrite (CINT z:com(C(t,r),L(t,r)). F(z)) to (CINT z:C(t,r). F(z)) + (CINT z:L(t,r). F(z))
    rewrite (CINT z:L(t,r). F(z)) to (-(INT x:[-r,r]. G(x)))
    simplify
done

lhs:
    apply 1 on ...
    apply 2 on ...
    apply 3 on ...
    simplify
done
```

## 组合建议

- 命令细节：加载 `apply-subgoal`、`state-func-def`、`state-prove`。
- 复杂重写：加载 `rewrite`。
- 失败恢复：加载 `error-recovery`。

