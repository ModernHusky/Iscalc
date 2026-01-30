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
## 使用时机
需要交换算符顺序以简化计算时：
1. **积分与求导**：`D x. INT x. f(x)` 或 `INT x. D x. f(x)`
2. **积分与求和**：`INT x. SUM(n, ...)` 或 `SUM(n, ... INT x ...)`
## 指令
### 快速开始

- 交换导数与积分：
```
exchange derivative and integral
```

- 交换积分与求和：
```
exchange integral and sum
```

### 详细说明



### 交换导数与积分 (DerivIntExchange)

允许交换微分算子和积分算子的顺序。
- 适用：
    1. `D x. INT x. f(x)` -> `INT x. D x. f(x)`
    2. `INT x. D x. f(x)` -> `D x. INT x. f(x)`

### 交换积分与求和 (IntSumExchange)

允许交换积分算子和求和算子的顺序。
- 适用：`INT x. SUM(n, l, u, f(n))` -> `SUM(n, l, u, INT x. f(n))`
- **重要前提**：级数必须一致收敛。如果系统需要证明收敛性，请先建立 `subgoal` 证明 `converges(...)`。
- 如果提示 "Applying the rule has no effect"，请检查表达式形式是否严格匹配。

### 示例

表达式: `INT x:[0,1]. SUM(n, 1, oo, x^n)`

```
exchange integral and sum
```

结果: `SUM(n, 1, oo, INT x:[0,1]. x^n)`
（之后可以对内部积分应用 `apply integral identity`）
