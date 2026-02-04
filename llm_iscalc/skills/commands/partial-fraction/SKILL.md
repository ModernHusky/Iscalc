---
name: partial-fraction
description: 部分分式分解（用于有理函数积分）。
match_rules:
- (?i)int
---

# partial-fraction
> 部分分式分解（用于有理函数积分）。

## 指令

### 语法格式

```
partial fraction decomposition
```

### 参数说明

无参数。系统自动识别有理函数并进行分解。

## 使用时机

当被积函数是有理函数且分母可分解时：
- 分母是多项式的乘积
- 分母可因式分解
- 如 `1/(x^2-1)`, `1/(x*(x+1))`

## 适用目标类型

| 目标类型 | 适用性 |
|---------|--------|
| 有理函数积分 `INT x. P(x)/Q(x)` | ✅ |
| 真分式（分子次数 < 分母次数）| ✅ |
| 假分式 | ⚠️ 需先做多项式除法 |
| 非有理函数 | ❌ |

## 工作流程

1. **识别有理函数**：确认被积函数是两个多项式的商 `P(x)/Q(x)`
2. **因式分解**：对分母 `Q(x)` 进行因式分解
3. **待定系数**：将分数拆分为简单分式之和
4. **求解系数**：通过待定系数法求出各项系数
5. **输出结果**：生成分解后的表达式

**分解规则**：
- 线性因子 `(x-a)` → `A/(x-a)`
- 重复线性因子 `(x-a)^n` → `A₁/(x-a) + A₂/(x-a)² + ... + Aₙ/(x-a)ⁿ`
- 不可约二次因子 `(x²+bx+c)` → `(Ax+B)/(x²+bx+c)`

## 注意事项

- **仅限有理函数**：被积函数必须是两个多项式的商。
- **应用范围**：此规则会自动判断，如果不是有理函数会报错。
- **后续处理**：分解通常会产生对数 (`log`) 或反正切 (`arctan`) 形式的积分，务必紧接着使用 `apply integral identity`。

## 示例

### 示例1: 二次分母

表达式: `INT x. 1/(x^2-1)`

```
partial fraction decomposition
```

结果: `INT x. 1/(2*(x-1)) - 1/(2*(x+1))`

### 示例2: 线性因子乘积

表达式: `INT x. 1/(x*(x+1))`

```
partial fraction decomposition
```

结果: `INT x. 1/x - 1/(x+1)`

### 示例3: 基本部分分式 (来自 theories/)

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

### 示例4: 换元后的部分分式 (来自 theories/)

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

### 示例5: 高次多项式分解 (来自 theories/)

```
calculate INT x:[3, 4]. (x^4 + x^3 + x^2 + 1)/(x^2 + x - 2)
    partial fraction decomposition
    simplify
    apply integral identity
    simplify
done
```
