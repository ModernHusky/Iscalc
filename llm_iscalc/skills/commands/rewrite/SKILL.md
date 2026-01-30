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
## 使用时机
需要进行代数恒等变换时：
- 三角恒等式: `sin^2 + cos^2 = 1`
- 对数性质: `log(a) - log(b) = log(a/b)`
- 分数化简: `x/(x+1) = 1 - 1/(x+1)`
- 处理负号顺序
## 指令
### 快速开始

```
rewrite <old_expr> to <new_expr>
```

将表达式中的 old_expr 替换为等价的 new_expr。

### 常用变换模式



### 三角恒等式

```
rewrite sin(x)^2 to 1 - cos(x)^2
rewrite sin(x)^2 to (1 - cos(2*x))/2
rewrite cos(x)^2 to (1 + cos(2*x))/2
rewrite tan(x) to sin(x)/cos(x)
```

### 对数性质

```
rewrite log(a) - log(b) to log(a/b)
rewrite log(a) + log(b) to log(a*b)
rewrite n*log(x) to log(x^n)
```

### 极限中的分数变换（重要！）

```
rewrite x/(x+1) to 1 - 1/(x+1)
rewrite x/(x-1) to 1 + 1/(x-1)
rewrite (x+1)/x to 1 + 1/x
rewrite (x-1)/x to 1 - 1/x
```

### 处理负号

```
rewrite -log(x+1) + log(x) to log(x) - log(x+1)
```

### 示例



### 极限计算关键技巧

表达式: `LIM {x -> oo}. x/(x+1)`

**先重写分数形式**：
```
rewrite x/(x+1) to 1 - 1/(x+1)
```

然后 `simplify` 得到结果 1。

### 对数合并

表达式: `log(x) - log(x+1)`

```
rewrite log(x) - log(x+1) to log(x/(x+1))
```

### 重要提醒 (注意事项)

1.  **括号与结合律**：
    - 请注意隐藏的括号，运算通常左结合。
    - 例如：在 `a * b * c` 中重写 `b * c` 会失败，因为系统将其视为 `(a * b) * c`。
    - **解决方法**：先重写 `a * b * c` to `a * (b * c)`，然后再重写内部的 `b * c`。
    - 同理，在 `c * a * b / d` 中找不到 `a * b`，应先重写为 `c * (a * b) / d`。

2.  **重写范围**：
    - 旧表达式 (`old_expr`) 必须**精确**出现在原表达式中。
    - 计算两边必须严格相等。如果工具提示不相等，尝试分解为更小的步骤。
      - 例如：不要直接写 `sqrt(1-sin(x)^2)^(1/2)` to `cos(x)`。
      - 而是分步：先 `1-sin(x)^2` to `cos(x)^2`，再 `sqrt(cos(x)^2)` to `cos(x)`。

3.  **避免滥用**：
    - **不要**用 `rewrite` 来做标准的积分恒等变换（如 `1/(x^2+1)` to `atan`），请使用 `apply integral identity`。
    - **不要**用它来展开求和（SUM），请用级数展开规则。
    - **不要**用它处理求和的线性性（提公因数），先 rewrite 内部项，再 `simplify`。

4.  **处理负号**：
    - 遇到 `-b` 找不到的情况（如 `a - b`），可能是因为表达式结构是 `a + (-b)` 或其他。尝试重写 `b` 或者先调整符号。
    - 示例：`rewrite -log(x+1) + log(x) to log(x) - log(x+1)` 再合并。

5.  **极限技巧**：
    - 遇到 `a/(a+b)` 形式的极限，考虑重写为 `1 - b/(a+b)`。

### 完整示例 (来自 theories/)



### 示例1: 三角恒等变换链

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

### 示例2: 对数合并

```
prove (INT x:[0,pi / 2]. log(a * sin(x))) = pi/2 * log(a/2) for a > 0
...
    rewrite log(a*cos(x)) to log(a) + log(cos(x))
    rewrite log(a * sin(x)) + log(a) to log(a * sin(x) * a)
    ...
done
```

### 示例3: 分数变换与部分分式准备

```
calculate INT x:[0, 1]. 1 / (exp(x) + 1)
    rewrite 1 / (exp(x) + 1) to exp(x) / (exp(x) * (exp(x) + 1))
    substitute u for exp(x)
    partial fraction decomposition
    apply integral identity
    simplify
done
```
