---
name: exchange-operators
description: 交换积分与求和/极限的顺序。
match_rules:
- (?i)int
- (?i)sum
- (?i)deriv
---

# exchange-operators
> 交换积分与求和/极限的顺序。

## 指令

### 语法格式

**交换导数与积分**：
```
exchange derivative and integral
```

**交换积分与求和**：
```
exchange integral and sum
```

**交换两个积分**：
```
exchange integral and integral
```

### 参数说明

无参数。系统自动识别表达式中的算符并执行交换。

## 使用时机

需要交换算符顺序以简化计算时：
1. **积分与求导**：将 `D x. INT y. f(x,y)` 变为 `INT y. D x. f(x,y)`（费曼技巧）
2. **积分与求和**：将 `INT x. SUM(n, ...)` 变为 `SUM(n, ... INT x ...)`
3. **两个积分**：交换二重积分的顺序（Fubini定理）

## 适用目标类型

| 目标类型 | 命令 | 适用性 |
|---------|------|--------|
| `D x. INT y. f(x,y)` | `exchange derivative and integral` | ✅ |
| `INT y. D x. f(x,y)` | `exchange derivative and integral` | ✅ |
| `INT x. SUM(n, l, u, f)` | `exchange integral and sum` | ✅ |
| `SUM(n, l, u, INT x. f)` | `exchange integral and sum` | ✅ |
| `INT x. INT y. f(x,y)` | `exchange integral and integral` | ✅ |

## 工作流程

### 交换导数与积分 (DerivIntExchange)

1. **识别结构**：找到形如 `D x. INT y. f(x,y)` 或 `INT y. D x. f(x,y)` 的表达式
2. **变量分离**：确保求导变量和积分变量不同
3. **执行交换**：交换导数和积分算符的位置

### 交换积分与求和 (IntSumExchange)

1. **识别结构**：找到形如 `INT x. SUM(n, l, u, f(n))` 的表达式
2. **收敛性检查**：验证级数一致收敛（可能需要先证明 subgoal）
3. **执行交换**：生成 `SUM(n, l, u, INT x. f(n))`

### 交换两个积分 (IntExchange)

1. **识别结构**：找到二重积分
2. **Fubini条件检查**：验证可积性条件
3. **执行交换**：交换积分顺序

## 注意事项

- **收敛性前提**：对于 `exchange integral and sum`，级数必须一致收敛。如果系统需要证明收敛性，请先建立 `subgoal` 证明 `converges(...)`。
- **无效操作**：如果提示 "Applying the rule has no effect"，请检查表达式形式是否严格匹配所需结构。
- **变量冲突**：确保交换后不会产生变量名冲突。

## 示例

### 示例1: 交换积分与求和

表达式: `INT x:[0,1]. SUM(n, 1, oo, x^n)`

```
exchange integral and sum
```

结果: `SUM(n, 1, oo, INT x:[0,1]. x^n)`

（之后可以对内部积分应用 `apply integral identity`）

### 示例2: 交换导数与积分（费曼技巧）

表达式: `D a. INT x:[0,oo]. exp(-a*x^2)`

```
exchange derivative and integral
```

结果: `INT x:[0,oo]. D a. exp(-a*x^2)`

### 示例3: 完整工作流程 (来自 theories/)

```
prove (INT x:[0,oo]. x^2 * exp(-a*x^2)) = sqrt(pi)/(4*a^(3/2)) for a > 0
# 使用费曼积分技巧
from 1:
    differentiate both sides at a
    expand definition for I
    exchange derivative and integral   # 关键步骤：交换导数和积分
    simplify
    ...
done
```
