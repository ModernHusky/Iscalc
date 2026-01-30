---
name: integrate-by-parts
description: 分部积分法（用于乘积形式函数积分）。
match_rules:
- (?i)int
---

# integrate-by-parts
> 分部积分法（用于乘积形式函数积分）。
## 使用时机
当被积函数是两个不同类型函数的乘积时：
- `x * exp(x)` - 多项式 × 指数
- `x * sin(x)` - 多项式 × 三角
- `log(x) * x^n` - 对数 × 多项式
- `exp(x) * sin(x)` - 指数 × 三角
## 指令
### 快速开始

```
integrate by parts with u = <expr>, v = <expr>
```

其中 u 是要求导的部分，v 是要积分后的结果。

### 分步指南

1. 将被积函数分解为 `u' * v` 的形式
2. 选择 u'（将被求导）和 v（将被积分）
3. 执行 `integrate by parts with u = u', v = v`
4. 公式转换: ∫u' * v = u' * V - ∫u'' * V（其中V是v的积分）
5. 继续化简新积分

### LIATE 选择规则

选择 u'（被求导）的优先级：
1. **L**ogarithmic - 对数函数 log(x)
2. **I**nverse trig - 反三角函数 arctan(x)
3. **A**lgebraic - 代数函数 x, x^2
4. **T**rigonometric - 三角函数 sin(x), cos(x)
5. **E**xponential - 指数函数 exp(x)

### 示例



### 示例1: 多项式 × 指数

表达式: `INT x:[0,1]. x * exp(x)`

```
integrate by parts with u = x, v = exp(x)
```

### 示例2: 对数积分

表达式: `INT x. log(x)`

```
integrate by parts with u = log(x), v = x
```

### 示例3: 需要两次分部积分

表达式: `INT x. exp(x) * sin(x)`

第一次分部后出现原积分，使用 `solve integral` 求解。

### 注意事项 (常见错误)

*   **最常见的错误**是给定的 `u` 和 `v` 不满足 `u * dv = 被积函数`。
    *   错误信息通常形式为：`u * dv does not equal body: t1 != t2`。
    *   请仔细检查 `u` 和 `v` 的选择。这里的 `v` 参数应该是积分后的结果（即 $dv$ 的原函数），而不是被积函数的一部分。
    *   如果由 `integrate by parts` 产生的 $\int v du$ 比原积分更复杂，说明 $u, v$ 选择可能不当，或者此方法不适用。
*   有时候两边实际上是相等的，但工具无法验证。这种情况下，请先尝试对被积函数使用 `rewrite` 或 `simplify`。
*   如果分部积分后出现原积分（循环积分），请在下一步使用 `solve integral` 命令。

### 完整示例 (来自 theories/)



### 示例: exp(x) * cos(x) 循环积分

```
calculate INT x:[0, 1]. exp(x) * cos(x)
    integrate by parts with u = exp(x), v = sin(x)
    simplify
    integrate by parts with u = exp(x), v = -cos(x)
    simplify
    solve integral INT x:[0, 1]. exp(x) * cos(x)
done
```

### 示例: log(x)^2

```
calculate INT x:[1, exp(1)]. log(x) ^ 2
    integrate by parts with u = log(x) ^ 2, v = x
    simplify
    integrate by parts with u = log(x), v = x
    apply integral identity
    simplify
done
```
