---
name: apply-integral-identity
description: 应用基本积分公式（幂、三角、指数等）。
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

当被积函数已化简为可直接积分的基本形式时使用：
- 多项式: `x^n`
- 三角函数: `sin(x)`, `cos(x)`, `tan(x)`
- 指数对数: `exp(x)`, `1/x`
- 反三角函数相关: `1/(1+x^2)`, `1/sqrt(1-x^2)`

## 适用目标类型

| 目标类型 | 适用性 |
|---------|--------|
| 定积分 `INT x:[a,b]. f(x)` | ✅ |
| 不定积分 `INT x. f(x)` | ✅ |
| 广义积分 `INT x:[0,oo]. f(x)` | ✅ |
| 普通表达式 | ❌ |

## 工作流程

1. **模式匹配**：扫描被积函数，尝试匹配已知积分公式表
2. **线性性处理**：自动拆分加法 `INT (f+g) = INT f + INT g`
3. **常数提取**：自动提取常数因子 `INT c*f = c * INT f`
4. **公式应用**：将匹配的被积函数替换为其原函数
5. **边界代入**（定积分）：计算 `F(b) - F(a)`

### 支持的积分公式

| 被积函数 | 积分结果 |
|---------|---------
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
- **匹配失败**：如果提示无效，通常需要通过 `simplify`、`expand polynomial` 或 `rewrite` 先将被积函数变形为基本形式。
- **不要用于复杂分式**：对于有理函数，先尝试 `partial fraction decomposition`。

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

### 示例3: 组合使用 (来自 theories/)

```
calculate INT x:[0, 1]. x*exp(x)
    integrate by parts with u = x, v = exp(x)
    apply integral identity    # 计算 INT exp(x) 得到 exp(x)
    simplify
done
```
