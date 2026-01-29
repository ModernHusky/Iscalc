---
name: partial-fraction
description: 部分分式分解 - 将有理函数分解为简单分式之和。当被积函数是分式且分母可分解时使用。
keywords:
  - 部分分式
  - 分式分解
  - 有理函数
applicable_types:
  - definite_integral
  - indefinite_integral
  - rational_function
match_rules:
  - "(?i)int"
---

# 何时使用

当被积函数是有理函数且分母可分解时：
- 分母是多项式的乘积
- 分母可因式分解
- 如 `1/(x^2-1)`, `1/(x*(x+1))`

# 快速开始

```
partial fraction decomposition
```

自动将有理函数分解为部分分式。

# 示例

## 示例1: 二次分母
表达式: `INT x. 1/(x^2-1)`

```
partial fraction decomposition
```

结果: `INT x. 1/(2*(x-1)) - 1/(2*(x+1))`

## 示例2: 线性因子乘积
表达式: `INT x. 1/(x*(x+1))`

```
partial fraction decomposition
```

结果: `INT x. 1/x - 1/(x+1)`

# 分步指南

1. 确认被积函数是真分式（分子次数 < 分母次数）
2. 如果是假分式，先做多项式除法
3. 应用 `partial fraction decomposition`
4. 对各项分别积分
5. `simplify` 合并结果

# 注意事项

# 注意事项

- **仅限有理函数**：被积函数必须是两个多项式的商。
- **应用范围**：此规则会自动判断，如果不是有理函数会报错。
- **后续处理**：分解通常会产生对数 (`log`) 或反正切 (`arctan`) 形式的积分，务必紧接着使用 `apply integral identity`。

# 完整示例 (来自 theories/)

## 示例1: 基本部分分式
```
calculate INT x:[3, 4]. 1 / (x ^ 2 - 4)
    partial fraction decomposition
    simplify
    substitute u for 4 * x + 8
    substitute u for 4 * x - 8 (at 2)
    apply integral identity
    simplify
done
```

## 示例2: 换元后的部分分式
```
calculate INT x:[1, 2]. exp(x) / ((exp(x) - 1) * (exp(x) + 3))
    substitute u for exp(x)
    partial fraction decomposition
    simplify
    substitute v for 4 * u + 12
    substitute v for 4 * u - 4 (at 2)
    apply integral identity
    simplify
done
```

## 示例3: 高次多项式分解
```
calculate INT x:[3, 4]. (x^4 + x^3 + x^2 + 1)/(x^2 + x - 2)
    partial fraction decomposition
    simplify
    apply integral identity
    simplify
done
```
