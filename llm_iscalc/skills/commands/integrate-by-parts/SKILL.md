---
name: integrate-by-parts
description: 分部积分法（用于乘积形式函数积分）。
match_rules:
- (?i)int
---

# integrate-by-parts
> 分部积分法（用于乘积形式函数积分）。

## 指令

### 语法格式

```
integrate by parts with u = <expr>, v = <expr>
```

### 参数说明

- `u`: 要求导的部分（将变成 `du`）
- `v`: 被积函数中的另一部分的**积分结果**（即 `∫dv`）

**重要**：这里的 `v` 是积分后的结果，不是被积函数本身！

## 使用时机

当被积函数是两个不同类型函数的乘积时：
- `x * exp(x)` — 多项式 × 指数
- `x * sin(x)` — 多项式 × 三角
- `log(x) * x^n` — 对数 × 多项式
- `exp(x) * sin(x)` — 指数 × 三角（循环积分）

## 适用目标类型

| 目标类型 | 适用性 |
|---------|--------|
| 定积分 `INT x:[a,b]. f(x)*g(x)` | ✅ |
| 不定积分 `INT x. f(x)*g(x)` | ✅ |
| 非乘积形式 | ⚠️ 需先变形 |

## 工作流程

1. **识别乘积结构**：将被积函数分解为 `u' * v'` 形式
2. **选择 u 和 dv**：根据 LIATE 规则选择
3. **验证匹配**：检查 `u * v' = 被积函数`
4. **公式应用**：`∫u'·v dx = u·v - ∫u·v' dx`
5. **继续求解**：对新积分继续化简

### LIATE 选择规则

选择 `u`（被求导）的优先级：
1. **L**ogarithmic — 对数函数 `log(x)`
2. **I**nverse trig — 反三角函数 `arctan(x)`
3. **A**lgebraic — 代数函数 `x`, `x^2`
4. **T**rigonometric — 三角函数 `sin(x)`, `cos(x)`
5. **E**xponential — 指数函数 `exp(x)`

## 注意事项

- **最常见的错误**是给定的 `u` 和 `v` 不满足 `u * dv = 被积函数`。
    - 错误信息通常形式为：`u * dv does not equal body: t1 != t2`。
    - 请仔细检查 `u` 和 `v` 的选择。**`v` 参数应该是积分后的结果**，而不是被积函数的一部分。
- **循环积分**：如果分部积分后出现原积分（如 `exp(x)*sin(x)`），请在下一步使用 `solve integral` 命令。
- **验证失败**：有时候两边实际上是相等的，但工具无法验证。这种情况下，请先尝试对被积函数使用 `rewrite` 或 `simplify`。
- **复杂度增加**：如果由 `integrate by parts` 产生的 `∫v du` 比原积分更复杂，说明 `u, v` 选择可能不当。

## 示例

### 示例1: 多项式 × 指数

表达式: `INT x:[0,1]. x * exp(x)`

```
integrate by parts with u = x, v = exp(x)
```

### 示例2: 对数积分

表达式: `INT x. log(x)`
（将 `log(x)` 视为 `log(x) * 1`）

```
integrate by parts with u = log(x), v = x
```

### 示例3: exp(x) * cos(x) 循环积分 (来自 theories/)

```
calculate INT x:[0, 1]. exp(x) * cos(x)
    integrate by parts with u = exp(x), v = sin(x)
    simplify
    integrate by parts with u = exp(x), v = -cos(x)
    simplify
    solve integral INT x:[0, 1]. exp(x) * cos(x)
done
```

### 示例4: log(x)^2 (来自 theories/)

```
calculate INT x:[1, exp(1)]. log(x) ^ 2
    integrate by parts with u = log(x) ^ 2, v = x
    simplify
    integrate by parts with u = log(x), v = x
    apply integral identity
    simplify
done
```

## 相关技能

> 💡 `<|load_skill|>apply-integral-identity<|end_load_skill|>` — 分部积分后应用恒等式

> 💡 `<|load_skill|>simplify<|end_load_skill|>` — 简化分部积分结果

> 💡 `<|load_skill|>integrate-by-equation<|end_load_skill|>` — 循环积分时使用 solve integral
