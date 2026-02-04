---
name: integrate-by-equation
description: 通过解方程求积分（用于循环积分等情况）。
match_rules:
  - "(?i)int"
applicable_types:
  - definite_integral
  - indefinite_integral
---

# integrate-by-equation
> 通过解方程求积分（用于循环积分等情况）。

## 指令

### 语法格式

```
solve integral <expr>
```

### 参数说明

- `<expr>`: 需要求解的积分表达式，必须与计算过程中出现的积分完全匹配

## 使用时机

当当前的积分在早期的计算步骤中已经出现过，并且可以通过解线性方程来获得积分值时使用：
- **循环积分**：如 `∫ e^x sin(x) dx` 在两次分部积分后出现原积分
- **自我指向积分**：计算结果包含原积分的线性组合

## 适用目标类型

| 目标类型 | 适用性 |
|---------|--------|
| 定积分 `INT x:[a,b]. f(x)` | ✅ |
| 不定积分 `INT x. f(x)` | ✅ |
| 表达式包含原积分的线性组合 | ✅ |

## 工作流程

1. **识别循环**：检测当前表达式是否包含被求解的积分
2. **构建方程**：设原积分为 `I`，当前表达式形如 `a*I + b`（其中 `a ≠ 1`）
3. **求解方程**：解出 `I = b / (1 - a)`
4. **代入结果**：将解得的值作为最终结果

**数学原理**：
如果计算过程得到 `I = g(x) + c*I`（其中 `c ≠ 1`），则：
`I - c*I = g(x)` → `I(1-c) = g(x)` → `I = g(x)/(1-c)`

## 注意事项

- **必须确保**被求解的积分完全出现在早期的步骤中。
- 当前表达式必须是该积分的线性函数（例如 `a * I + b`，其中 `a != 1`）。
- **积分匹配**：`<expr>` 必须与原积分精确匹配，包括积分变量和边界。

## 示例

### 示例1: 循环积分 exp(x)*sin(x)

假设计算过程如下：

```
  INT x:[0,pi / 2]. cos(x) * exp(2 * x)
= ...
= -(4 * (INT x:[0,pi / 2]. cos(x) * exp(2 * x))) + exp(pi) - 2
```

此时，表达式中包含原始积分（带有系数 -4）。要解出该积分，使用：

```
solve integral INT x:[0, pi/2]. cos(x) * exp(2*x)
```

这将自动解方程并给出最终答案：`(exp(pi) - 2) / 5`

### 示例2: 完整工作流程 (来自 theories/)

```
calculate INT x:[0, 1]. exp(x) * cos(x)
    integrate by parts with u = exp(x), v = sin(x)
    simplify
    integrate by parts with u = exp(x), v = -cos(x)
    simplify
    # 此时表达式包含原积分
    solve integral INT x:[0, 1]. exp(x) * cos(x)
done
```

### 示例3: 不定积分循环

```
calculate INT x. exp(x) * sin(x)
    integrate by parts with u = exp(x), v = -cos(x)
    simplify
    integrate by parts with u = exp(x), v = sin(x)
    simplify
    solve integral INT x. exp(x) * sin(x)
    # 结果：exp(x) * (sin(x) - cos(x)) / 2
done
```
