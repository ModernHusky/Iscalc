---
name: substitute
description: 积分换元法（令 u=expr 化简复合函数）。
match_rules:
- (?i)int
- improper_integral
---

# substitute
> 积分换元法（令 u=expr 化简复合函数）。
## 使用时机
当被积函数是复合函数形式时，用换元法将其化简为基本形式：
- 含有 `f(g(x))` 形式的复合函数
- 根式、三角函数、指数函数的复合
## 指令
### 快速开始

```
substitute u for <expr>
```

将令 u = expr 进行换元。

### 分步指南

1. 识别被积函数中的复合结构 `f(g(x))`
2. 令 `u = g(x)`，即内层函数
3. 执行 `substitute u for g(x)`
4. 化简后应用积分恒等式
5. 使用 `replace substitution` 代回原变量

### 示例



### 示例1: 多项式换元

表达式: `INT x:[0,1]. (x+1)^3`

```
substitute u for x + 1
```

结果: `INT u:[1,2]. u^3`

### 示例2: 三角函数换元

表达式: `INT x. sin(x)^2 * cos(x)`

```
substitute u for sin(x)
```

结果: `INT u. u^2`

### 示例3: 指数函数换元

表达式: `INT x. exp(x) / (1 + exp(x))`

```
substitute u for exp(x)
```

### 常见换元模式

| 被积函数形式 | 换元建议 |
|------------|---------|
| 根式 √(expr) | u = expr 或 u = √(expr) |
| sin(x), cos(x) | u = sin(x) 或 u = cos(x) |
| exp(x) | u = exp(x) |
| log(x) | u = log(x) |
| 线性复合 f(ax+b) | u = ax+b |

### 注意事项

1.  **修饰符语法（多个积分时）**：
    - 当表达式中有多个积分时，使用 `(at n)` 修饰符指定对第几个积分换元
    - **正确格式**: `substitute u for expr (at 2)` （注意括号）
    - **错误格式**: `substitute u for expr at 2` （缺少括号会导致解析错误）
    - **示例**: 
      ```
      表达式: (INT x:[0,1]. x^2) + (INT y:[0,1]. y^3)
      命令: substitute u for y (at 2)  # 对第二个积分换元
      ```

2.  **关于换元与反解**：
    - 当令 `u = g(x)` 时，被积函数最好能写成 `f(g(x)) * g'(x)` 的形式。
    - 如果方程 `u = g(x)` 难以反解出 `x`，请先尝试提取 `g'(x)`。

3.  **反向换元 (Inverse Substitution)**：
    - 使用 `substitute f(u) for x` 将积分变量 `x` 替换为含新变量 `u` 的表达式。
    - 此时积分变成 `INT u:[a',b']. g(f(u)) * f'(u)`。
    - **注意**：`x` 必须是积分变量。

4.  **变量恢复 (Recover Variable)**：
    - **仅限不定积分**：在计算完成后，**必须**使用 `replace substitution` 将中间变量 `u` 换回原变量 `x`。
    - **定积分**：不需要此步骤，因为积分限已经变换。

5.  **常见错误**：
    - `incorrect old variable y, should be x`：试图替换非积分变量。
    - 确保新变量 (如 `u`) 在之前的步骤中未被定义或占用。

### 完整示例 (来自 theories/)

### 示例1: 高斯积分缩放换元

```
# 将 exp(-(a*x^2)) 转换为标准形式 exp(-(u^2))
calculate INT x:[0,oo]. exp(-(a * x ^ 2))
    substitute u for sqrt(a) * x
    # 结果: INT u:[0,oo]. exp(-(u ^ 2)) / sqrt(a)
    simplify
    # 提取常数: (1/sqrt(a)) * INT u:[0,oo]. exp(-(u ^ 2))
done
```

### 示例2: 链式换元

```
calculate INT x:[0, 1]. x^3/(1+x^4)^(1/4)
    substitute u for x ^ 4
    simplify
    substitute v for 1 + u
    apply integral identity
    simplify
done
```

### 示例3: 三角函数与反向换元

```
prove (INT x:[1,oo]. 1 / ((x+a)*sqrt(x-1))) = pi / sqrt(a+1) for a > -1
lhs:
    substitute t for sqrt(x - 1)
    simplify
    substitute y for t / sqrt(a + 1)
    rewrite y ^ 2 * (a + 1) + a + 1 to (a + 1) * (y^2 + 1)
    apply integral identity
    simplify
done
```

### 示例4: 对数换元

```
calculate INT x:[exp(1), exp(2)]. 3/(x*log(x))
    substitute u for log(x)
    apply integral identity
    simplify
done
```
