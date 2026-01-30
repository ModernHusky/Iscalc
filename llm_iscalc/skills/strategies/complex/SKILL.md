---
name: complex
description: 复数运算与欧拉公式变换。
match_rules:
- (?i)complex
- (?i)imaginary
- (?i)复数
- (?i)虚数
- (?i)partial.*fraction
- (?i)decompose
- (?i)factor
---

# complex
> 复数运算与欧拉公式变换。
## 使用时机
- 当用户明确要求“使用复数方法”、“在复数域分解”时。
- 当有理函数在实数域无法进一步分解（例如分母为 `x^2+1`, `x^2+a^2`），但在复数域可以进一步处理时。
- 当涉及到虚数单位 `i` 的化简或重写时。
- **注意**：底层的 `partial fraction decomposition` 命令仅支持实数域分解。你需要使用 `rewrite` 命令手动进行复数分解。
## 指令
### 核心策略：使用 Rewrite 进行复数分解

对于形如 `1/(x^2 + a^2)` 的表达式，通过恒等式进行重写：

$$ \frac{1}{x^2+a^2} = \frac{1}{2ai} \left( \frac{1}{x-ai} - \frac{1}{x+ai} \right) $$

### 命令示例



### 示例 1: `1/(x^2 + 1)`

**策略**: 分解为 `(x-i)` 和 `(x+i)`。

```
rewrite 1/(x^2 + 1) to 1/(2*i) * (1/(x-i) - 1/(x+i))
```

### 示例 2: `1/(x^2 + 4)`

这里 `a=2`。

```
rewrite 1/(x^2 + 4) to 1/(4*i) * (1/(x-2*i) - 1/(x+2*i))
```

### 示例 3: `1/(u^2 + 1)`

```
rewrite 1/(u^2 + 1) to 1/(2*i) * (1/(u-i) - 1/(u+i))
```

### 其他策略



### 欧拉公式

处理复指数与三角函数转换：

- `rewrite exp(i*x) to cos(x) + i*sin(x)`
- `rewrite cos(x) to (exp(i*x) + exp(-i*x))/2`
- `rewrite sin(x) to (exp(i*x) - exp(-i*x))/(2*i)`

### 虚数化简

使用 `simplify` 自动处理：
- `i^2` -> `-1`
- `1/i` -> `-i`
