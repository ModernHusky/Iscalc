---
name: state-prove
description: 证明状态操作指南。
match_rules:
- .*
---

# state-prove
> 证明状态操作指南。
## 使用时机
参考 description 描述。
## 指令
### 证明状态 (Prove State)

在证明状态下，可以执行以下操作：

### 可用命令



### 开始计算证明

- `lhs:` - 对等式目标的**左边**进行计算 (进入 Calculate 状态)
- `rhs:` - 对等式目标的**右边**进行计算
- `arg:` - 对收敛性目标的**参数**进行计算

### 从已证子目标开始

- `from <name>:` - 从已证明的 subgoal 开始重写

### 归纳法

- `induction on <var>` - 对变量进行归纳
- `induction on <var> starting from <n>` - 指定归纳起点

### 案例分析

- `case analysis on <cond>` - 布尔条件 (true/false 两个分支)
- `case analysis on <expr>` - 数值表达式 (positive/zero/negative 三个分支)

### 子目标和定义

- `subgoal <name>: <expr>` - 声明子目标
- `define <expr>` - 创建定义

### 完成证明

- `done` - 当目标已被证明时使用

### 注意事项

- **不能**在 Prove 状态使用 `calculate` 命令
- 条件可用 `for <cond1>, <cond2>` 指定
- **等式证明的绝对红线（极其严厉的警告！）**：iscalc **绝对不具备**判定两个表达式"逻辑等价"或"数学等价"的能力，它只认识字面上的完全相等！你**必须**选择等式的某一边（`lhs:` 或 `rhs:`），通过一系列操作，**完整地**将该边推导至与等式另一边**字面上完全一致**的形式，然后才能用 `done` 完成证明。
  - ✅ 正确做法：`lhs:` → 连续变换 → 获得的表达式与原始等式的 RHS **每一个字符、每一项**都完全等同 → 输入 `done`
  - ❌ **绝对禁止的错误做法**：将一边变换为另一边的**子表达式**！例如：
    - 错误情况 A（带有被积函数或导数）：试图证明积分等式时，对 RHS 求导得到了 LHS 的"被积函数"，就以为证明结束了。不！你必须让两边带着积分号完全一致，或者让两边都算出积分结果后完全一致。
    - 错误情况 B（带有不定积分常数）：LHS 算出结果没有 `SKOLEM_CONST(C)`，但 RHS 带有 `SKOLEM_CONST(C)`，这属于**不相等**！你必须想办法在某一边也引入常数（或利用 `rewrite` 补齐常数），使两边在形式上完全匹配。
  如果未能达到字面完全一致就输入 `done`，系统将直接报错拒绝！

### 完整示例 (来自 theories/)



### 示例1: 基本等式证明

```
prove (INT x:[1,oo]. 1 / ((x+a)*sqrt(x-1))) = pi / sqrt(a+1) for a > -1
lhs:
    substitute t for sqrt(x - 1)
    simplify
    substitute y for t / sqrt(a + 1)
    apply integral identity
    simplify
done
```

### 示例2: 使用子目标

```
prove (INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)) = pi / 8 * log(2)
subgoal 1: (INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)) = (INT x:[0,pi / 4]. log(tan(x) + 1))
lhs:
    substitute tan(u) for x
    ...
done
subgoal 2: ...
from 2:
    solve equation for INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)
done
```

### 示例3: 案例分析

```
prove (INT x:[0,oo]. 1/(x^4+2*x^2*cos(2*a)+1)) = pi/abs((4*cos(a))) for cos(a) != 0
...
case analysis on cos(a)
    case negative:
    lhs:
        apply 8 on ...
    rhs:
        simplify
    done
    case positive:
    ...
done
```
