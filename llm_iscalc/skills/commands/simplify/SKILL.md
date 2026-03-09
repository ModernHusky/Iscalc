---
name: simplify
description: 通用表达式化简与计算。
match_rules:
- .*
---

# simplify
> 通用表达式化简与计算。

## 指令

### 语法格式

```
simplify
```

### 参数说明

无参数。该命令对当前整个表达式进行化简。

## 使用时机

在几乎所有操作后都应使用：
- 换元后化简积分
- 积分后化简结果
- 极限计算后简化
- 任何需要整理表达式的时候

## 适用目标类型

| 目标类型 | 适用性 |
|---------|--------|
| 任意数学表达式 | ✅ |
| 积分表达式 | ✅（化简被积函数）|
| 极限表达式 | ✅（尝试求值）|
| 导数表达式 | ✅（执行求导）|

## 工作流程

1. **线性性处理**：拆分加法，提取公因子
2. **常数提取**：将常数因子移到算符外
3. **代数化简**：合并同类项，约分等，以及**合并同底数的指数（如 `exp(A) * exp(B) -> exp(A+B)`）**和**聚拢离散的常数系数乘积（如 `-1 * (1/a) -> -1/a`）**
4. **导数计算**：对导数表达式执行求导
5. **极限求值**：对简单极限直接求值
6. **表达式规范化**：统一表达式的标准形式

## 注意事项

- **不做分式化简**：`simplify` 不会自动处理分式通分或拆分。请使用 `rewrite` (如 `x/(x+1)` 到 `1 - 1/(x+1)`) 或 `partial fraction decomposition`。
- **不做积分计算**：它只化简表达式，不计算积分。请用 `apply integral identity`。
- **条件依赖**：某些化简（如 `abs(x) -> x`）依赖于变量条件（如 `x >= 0`）。如果发现未化简，可能是条件不足。
- **多次交替使用**：复杂的表达式可能需要多次 `simplify` 与 `rewrite` 交替使用。
- **深层整理认知**：当结果看似正确，但系数散落（如含有负号和倒数），或者存在多个同底数项的乘积（如多个 `exp` 项的乘法）时，绝不要立刻认为“已是最简”。此时调用一次 `simplify` 能进一步聚合到底。例如 LLM 在推理出 `-exp(-b) * (1 / a * exp(a * x))` 时，应立刻想到还能用 `simplify` 一步化简成 `-exp(a * x - b) / a`。
- **酌情done**：当表达式看似最简，且系数已聚合到底，可以酌情使用 `done`。

> ⚠️ **重要：simplify 会重排表达式项的顺序**
>
> `simplify` 内部使用规范化算法，会按照表达式复杂度从大到小排序。例如：
> - 你 rewrite 得到：`log(t) - log(t^2+1)/2`
> - simplify 后变成：`-log(t^2+1)/2 + log(t)`
>
> **这是设计行为，无法避免。** 正确的工作流是：
> 1. **先完成所有 rewrite 操作**（包括对数合并等）
> 2. **最后再调用 simplify**
>
> 📝 **示例：极限中处理对数差**
> ```
> # ❌ 错误做法：rewrite 后立即 simplify
> rewrite -log(t^2+1)/2 + log(t) to log(t) - log(t^2+1)/2
> simplify  # 又会变回 -log(t^2+1)/2 + log(t)
>
> # ✅ 正确做法：先合并对数，再 simplify
> rewrite log(t) - log(t^2+1)/2 to log(t / sqrt(t^2+1))
> simplify  # 不会改变 log(t / sqrt(t^2+1)) 的形式
> ```

## 示例

### 示例1: 基本使用

应用于任何复杂表达式：
```
simplify
```

### 示例2: 换元后化简

```
calculate INT x:[0,pi]. x * sin(x) / (1 + cos(x) ^ 2)
    substitute y for pi - x
    expand polynomial
    simplify              # 换元后化简
    ...
done
```

### 示例3: 积分后化简 (来自 theories/)

```
calculate INT x:[0, 1]. x*exp(x)
    integrate by parts with u = x, v = exp(x)
    apply integral identity
    simplify              # 化简积分结果
done
```

### 示例4: 极限求值

```
calculate LIM {x -> 0}. sin(x)/x
    l'Hopital's rule
    simplify              # 求值得到 1
done
```

### 示例5: 代数式的深度整理聚合

应用 `simplify` 处理散落的常数项系数以及同底数的指数合并：

```
calculate ...
    ... # 假如你算到了 -exp(-b) * (1 / a * exp(a * x))，此时不应停滞
    simplify              # 通过系数聚合与指数合并，得到 -exp(a * x - b) / a
done
```
