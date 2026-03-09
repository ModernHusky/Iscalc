---
name: rewrite
description: 表达式代数恒等变换（三角、对数、分式等）。
match_rules:
- sin|cos|tan
- log
- /\(
- \^2
---

# rewrite
> 表达式代数恒等变换（三角、对数、分式等）。

## 指令

### 语法格式

```
rewrite <old_expr> to <new_expr>
```

将表达式中的 `old_expr` 替换为等价的 `new_expr`。

```
rewrite to <new_expr>
```

将整个当前表达式重写为 `new_expr`（需确保等价）。

> ⚠️ **重要：`rewrite to` 需要完整的表达式**
>
> 使用 `rewrite to` 时，`<new_expr>` 必须是**上一步的完整表达式**，包括外层结构（如 `LIM`、`INT` 等）。
>
> **示例**：
> - ✅ 正确：`rewrite to LIM {t -> oo}. log(2) / 2 + log(t / sqrt(t ^ 2 + 1))`
> - ❌ 错误：`rewrite to log(2) / 2 + log(t / sqrt(t ^ 2 + 1))`（缺少 `LIM {t -> oo}.`）
>
> 如果当前表达式是：
> - `LIM {x -> oo}. ...` → `rewrite to` 必须以 `LIM {x -> oo}.` 开头
> - `INT x:[a,b]. ...` → `rewrite to` 必须以 `INT x:[a,b].` 开头
> - 普通表达式 `a + b + c` → `rewrite to` 直接写变换后的表达式

### 参数说明

- `<old_expr>`: 需要被替换的子表达式（必须在原表达式中精确出现）
- `<new_expr>`: 替换后的等价表达式

## 使用时机

需要进行代数恒等变换时：
- 三角恒等式: `sin^2 + cos^2 = 1`
- 对数性质: `log(a) - log(b) = log(a/b)`
- 分数化简: `x/(x+1) = 1 - 1/(x+1)`
- 处理负号顺序

> 💡 **提示**：如果是表达式的多项式拆解展开（例如将 `(x+1)^2` 展开，或 `(x+1)(x-1)` 展开），**不需要**使用 `rewrite`，直接使用无参数的 `expand polynomial` 即可，它会自动扫描并展开函数中所有可展开的多项式结构。见`<|load_skill|>expand-polynomial<|end_load_skill|>` — 多项式乘法展开与整理。

## 适用目标类型

| 目标类型 | 适用性 |
|---------|--------|
| 任意数学表达式 | ✅ |
| 积分内部表达式 | ✅ |
| 极限表达式 | ✅ |
| 等式两边 | ✅ |

## 工作流程

1. **查找匹配**：在当前表达式中搜索 `old_expr` 的精确出现
2. **等价验证**：验证 `old_expr` 与 `new_expr` 数学上等价
3. **执行替换**：将找到的子表达式替换为 `new_expr`
4. **规范化**：对替换后的表达式进行规范化处理

## 注意事项

1. **括号与结合律**：
   - 请注意隐藏的括号，运算通常左结合。
   - 例如：在 `a * b * c` 中重写 `b * c` 会失败，因为系统将其视为 `(a * b) * c`。
   - **解决方法**：先重写 `a * b * c` to `a * (b * c)`，然后再重写内部的 `b * c`。

2. **重写范围**：
   - 旧表达式 (`old_expr`) 必须**精确**出现在原表达式中。
   - 两边必须严格相等。如果工具提示不相等，尝试分解为更小的步骤。
   - 例如：不要直接写 `sqrt(1-sin(x)^2)^(1/2)` to `cos(x)`。
   - 而是分步：先 `1-sin(x)^2` to `cos(x)^2`，再 `sqrt(cos(x)^2)` to `cos(x)`。

3. **避免滥用**：
   - **不要**用 `rewrite` 来做标准的积分恒等变换，请使用 `apply integral identity`。
   - **不要**用它来展开求和（SUM），请用级数展开规则。
   - **视情况**当进行多项式的乘法展开和恒等拆解时。推荐使用无参数的 `expand polynomial` 命令，它能自动一次性展开并在内部完成同类项合并。见`<|load_skill|>expand-polynomial<|end_load_skill|>` — 多项式乘法展开与整理。

4. **处理负号**：
   - 遇到 `-b` 找不到的情况（如 `a - b`），可能是表达式结构是 `a + (-b)`。

5. **极限技巧**：
   - 遇到 `a/(a+b)` 形式的极限，考虑重写为 `1 - b/(a+b)`。

6. **多项表达式的部分合并技巧**：

> ⚠️ **重要：表达式的左结合结构**
>
> 系统内部表达式采用**左结合**的二叉树结构：
> - `a + b + c` 实际是 `(a + b) + c`
> - `a - b + c` 实际是 `(a - b) + c`
> - `a + b - c` 实际是 `(a + b) - c`
>
> **问题**：当你想合并后两项 `b` 和 `c` 时，直接匹配 `b + c` 或 `b - c` 会失败，因为它们不是完整的子表达式。

**示例：合并对数表达式中的后两项(注意下述的是例子，可能还会有多个子表达式的情况，可能会有4个、5个或者更多，都要严格依照下述说明来给出命令)**

假设当前表达式为：`log(a) + log(b) + log(c)`（内部结构：`(log(a) + log(b)) + log(c)`）

```
# ❌ 错误做法：只匹配后两项（会失败！）
rewrite log(b) + log(c) to log(b * c)
# 失败原因：表达式树中不存在 "log(b) + log(c)" 这个完整子节点

# ✅ 正确做法1：使用 rewrite to 改写整个表达式
rewrite to log(a) + log(b * c)

# ✅ 正确做法2：匹配完整表达式进行转换
rewrite log(a) + log(b) + log(c) to log(a) + log(b * c)
```

**加减混合情况**

对于 `log(a) - log(b) + log(c)`（结构：`(log(a) - log(b)) + log(c)`）：

```
# ❌ 错误：尝试合并 -log(b) + log(c)（不存在这个子表达式）
rewrite -log(b) + log(c) to log(c/b)

# ✅ 正确：改写整个表达式
rewrite to log(a) + log(c/b)
# 或
rewrite log(a) - log(b) + log(c) to log(a) + log(c/b)
```


## 示例

### 常用变换模式

**三角恒等式**：
```
rewrite sin(x)^2 to 1 - cos(x)^2
rewrite sin(x)^2 to (1 - cos(2*x))/2
rewrite cos(x)^2 to (1 + cos(2*x))/2
rewrite tan(x) to sin(x)/cos(x)
```

**对数性质**：
```
rewrite log(a) - log(b) to log(a/b)
rewrite log(a) + log(b) to log(a*b)
rewrite n*log(x) to log(x^n)
```

**分数变换**：
```
rewrite x/(x+1) to 1 - 1/(x+1)
rewrite (x+1)/x to 1 + 1/x
```

### 示例1: 极限计算技巧

表达式: `LIM {x -> oo}. x/(x+1)`

```
rewrite x/(x+1) to 1 - 1/(x+1)
simplify
```

结果: `1`

### 示例2: 三角恒等变换链 (来自 theories/)

```
prove (INT x:[0,oo]. 1 / (x ^ 4 + 2 * x ^ 2 * cosh(2 * a) + 1)) = pi / (4 * cosh(a))
lhs:
    expand definition for cosh (all)
    rewrite x ^ 4 + 2 * x ^ 2 * ((exp(-(2 * a)) + exp(2 * a)) / 2) + 1 to (x ^ 2 + exp(2 * a)) * (x ^ 2 + exp(-(2 * a)))
    rewrite 1 / ((x ^ 2 + exp(2 * a)) * (x ^ 2 + exp(-(2 * a)))) to 1 / (exp(2 * a) - exp(-(2 * a))) * (1 / (x ^ 2 + exp(-(2 * a))) - 1 / (x ^ 2 + exp(2 * a)))
    simplify
    ...
done
```

### 示例3: 分数变换与部分分式准备 (来自 theories/)

```
calculate INT x:[0, 1]. 1 / (exp(x) + 1)
    rewrite 1 / (exp(x) + 1) to exp(x) / (exp(x) * (exp(x) + 1))
    substitute u for exp(x)
    partial fraction decomposition
    apply integral identity
    simplify
done
```
