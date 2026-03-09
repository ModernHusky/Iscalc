---
name: apply-integral-identity
description: 应用基本积分公式（幂、三角、指数等），但仅限被积函数已经是标准公理型；若是复合函数、线性复合或差一个导数因子，先考虑配合 substitute 做换元。
match_rules:
- (?i)int
- improper_integral
---

# apply-integral-identity
> 应用基本积分公式（幂、三角、指数等）。

## 指令

### 语法格式

```
apply integral identity
```

自动识别并应用匹配的积分公式。

### 参数说明

无参数。该命令会自动扫描被积函数，匹配已知的积分恒等式。

## 使用时机

仅当被积函数已经是或已经被你化简成可直接套用积分公式的基本形式时使用：
- 多项式: `x^n`
- 三角函数: `sin(x)`, `cos(x)`, `tan(x)`
- 指数对数: `exp(x)`, `1/x`
- 反三角函数相关: `1/(1+x^2)`, `1/sqrt(1-x^2)`

如果被积函数是复合形式、链式结构或“差一个导数因子”的形式，不要直接无脑调用本命令；先判断是否应加载并使用 `substitute`，把积分改写成标准公理型后再回来执行 `apply integral identity`。

### 先做这个判断

在调用本命令前，先快速检查下面两件事：

1. **是否已经是标准公式表中的变量形状**
   - 例如 `x^n`、`sin(x)`、`exp(x)`、`1/(1+x^2)`
2. **是否本质上是某个内层表达式 `g(x)` 的函数**
   - 例如 `f(g(x)) * g'(x)`、`1/sqrt(a+bx)`、`cos(ax+b)`、`exp(x^2) * x`

若第 1 条不明显成立，而第 2 条成立或接近成立，优先考虑 `substitute`。

**相关技能**：
- `<|load_skill|>substitute<|end_load_skill|>` — 换元积分详细语法和示例

## 适用目标类型

| 目标类型 | 适用性 |
|---------|--------|
| 定积分 `INT x:[a,b]. f(x)` | ✅ |
| 不定积分 `INT x. f(x)` | ✅ |
| 广义积分 `INT x:[0,oo]. f(x)` | ✅ |
| 普通表达式 | ❌ |

## 工作流程

1. **先判定是否需要预处理**：检查被积函数是否已经是标准公理型；若不是，优先考虑 `simplify`、`rewrite`，以及尤其是 `substitute`
2. **模式匹配**：仅对已成为基本形式的被积函数匹配已知积分公式表
3. **线性性处理**：自动拆分加法 `INT (f+g) = INT f + INT g`
4. **常数提取**：自动提取常数因子 `INT c*f = c * INT f`
5. **公式应用**：将匹配的被积函数替换为其原函数
6. **边界代入**（定积分）：计算 `F(b) - F(a)`

### 推荐决策顺序

遇到积分时，优先按以下顺序思考：

1. 被积函数是否已经是公式表中的基本形状
2. 如果不是，能否通过简单代数化简变成基本形状
3. 如果仍不是，是否存在明显内层表达式 `g(x)`，使其适合先做 `substitute`
4. 换元后若得到 `u^n`、`sin(u)`、`1/(1+u^2)`、`1/u` 等标准形式，再执行 `apply integral identity`
5. 不定积分在换元完成后记得恢复变量

### 支持的积分公式

| 被积函数 | 积分结果 |
|---------|---------|
| `x^n` | `x^(n+1)/(n+1)` |
| `1/x` | `log(x)` |
| `exp(x)` | `exp(x)` |
| `sin(x)` | `-cos(x)` |
| `cos(x)` | `sin(x)` |
| `1/(1+x^2)` | `arctan(x)` |
| `1/sqrt(1-x^2)` | `arcsin(x)` |
| `sec(x)^2` | `tan(x)` |
| `csc(x)^2` | `-cot(x)` |

## 注意事项

- **线性性**：该规则会自动处理积分的线性性（即 `INT (f+g) = INT f + INT g`），无需手动拆分。
- **不要把“匹配失败”当成死路**：很多失败不是公式不存在，而是当前变量形状不对。若被积函数看起来像 `f(g(x))` 或只差一个内层导数因子，先尝试 `substitute`。
- **优先识别换元信号**：以下情况通常应先考虑 `substitute`，而不是立即套公式：
  - `sin(ax+b)`、`cos(ax+b)`、`exp(ax+b)`、`1/(1+(ax+b)^2)`
  - `(ax+b)^n`
  - `f(g(x)) * g'(x)` 或可通过常数整理得到这种形状
  - 含根式、分式、对数的复合表达式，如 `1/sqrt(1-x^2)` 的线性或多项式复合版本
- **匹配失败后的默认动作**：优先尝试 `simplify`、`rewrite`；若仍不是标准型，加载 `substitute` 并先换元。
- **不要用于复杂分式**：对于有理函数，先尝试 `partial fraction decomposition`。
- **不定积分的换元链**：如果先做了 `substitute`，在应用完积分恒等式后，通常还需要 `replace substitution` 恢复原变量。

## 与 `substitute` 的配合

把本命令视为“换元后的收尾器”而不是默认第一步。

推荐配合模式：

1. 先识别内层表达式 `g(x)`
2. 用 `substitute` 把积分改写为关于新变量的标准形式
3. 对新积分调用 `apply integral identity`
4. 必要时 `simplify`
5. 若是不定积分，最后恢复原变量

常见联动示例：

| 原积分 | 更合适的做法 |
|-------|-------------|
| `INT x. cos(3*x+1)` | 先 `substitute u for 3*x + 1`，再 `apply integral identity` |
| `INT x. x*exp(x^2)` | 先 `substitute u for x^2`，再 `apply integral identity` |
| `INT x. 1/(1+(2*x-5)^2)` | 先 `substitute u for 2*x - 5`，再 `apply integral identity` |
| `INT x. (x+1)^7` | 先视情况直接展开，或更自然地 `substitute u for x+1` 后再套公式 |

## 示例

### 示例1: 幂函数

表达式: `INT x. x^3`

```
apply integral identity
```

结果: `x^4/4`

### 示例2: 反三角函数

表达式: `INT x. 1/(1+x^2)`

```
apply integral identity
```

结果: `arctan(x)`

### 示例3: 先换元再套公式

```
calculate INT x. x * exp(x^2)
    substitute u for x ^ 2
    apply integral identity
    replace substitution
    simplify
done
```

结果: `exp(x^2) / 2`

### 示例4: 线性复合函数先换元

表达式: `INT x. 1/(1+(2*x-5)^2)`

```
substitute u for 2*x - 5
apply integral identity
replace substitution
simplify
```

### 示例5: 组合使用 (来自 theories/)

```
calculate INT x:[0, 1]. x*exp(x)
    integrate by parts with u = x, v = exp(x)
    apply integral identity    # 这里的被积函数已经是标准型 INT exp(x)
    simplify
done
```

## 相关技能

> `substitute` — 当被积函数不是标准公理型、而是复合函数或线性复合时，优先先做换元

> `rewrite` — 将表达式改写为更容易识别的标准形状

> `simplify` — 在套公式前后做必要的代数整理
