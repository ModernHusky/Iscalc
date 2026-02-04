---
name: expand-polynomial
description: 多项式乘法展开与整理。
match_rules:
- (?i)int
- \^\s*\d
---

# expand-polynomial
> 多项式乘法展开与整理。

## 指令

### 语法格式

```
expand polynomial
```

### 参数说明

无参数。该命令会自动展开整个被积函数中的多项式。

## 使用时机

当被积函数包含需要展开的多项式时使用：
- 幂次表达式如 `(x+1)^3`
- 乘积表达式如 `(x+1)(x-1)`
- 展开后可直接用 `apply integral identity` 进行逐项积分

## 适用目标类型

| 目标类型 | 适用性 |
|---------|--------|
| 定积分 `INT x:[a,b]. poly(x)` | ✅ |
| 不定积分 `INT x. poly(x)` | ✅ |
| 普通多项式表达式 | ✅ |

## 工作流程

1. **识别多项式**：扫描被积函数，找到可展开的多项式结构
2. **分配律应用**：对 `(a+b)^n` 使用二项式定理展开
3. **乘法展开**：对 `(a+b)(c+d)` 应用分配律
4. **同类项合并**：将展开结果整理为标准多项式形式

## 注意事项

- **仅适用于积分内的多项式**：此规则主要针对被积函数。
- **后续处理**：展开后通常可以直接使用 `apply integral identity` 逐项积分。
- **复杂度考虑**：如果多项式次数很高，考虑使用换元法代替展开。
- **部分展开**：如果只需要展开某个子表达式，请使用 `rewrite` 手动指定。

## 示例

### 示例1: 二次多项式展开

表达式: `INT x. (x + 1)^2`

```
expand polynomial
```

结果: `INT x. x^2 + 2*x + 1`

### 示例2: 三次多项式展开

表达式: `INT x:[0,1]. (x+1)^3`

```
expand polynomial
```

结果: `INT x:[0,1]. x^3 + 3*x^2 + 3*x + 1`

### 示例3: 完整工作流程 (来自 theories/)

```
calculate INT x:[0,pi]. x * sin(x) / (1 + cos(x) ^ 2)
    substitute y for pi - x
    expand polynomial      # 展开 (pi - y) 等多项式
    simplify
    ...
done
```
