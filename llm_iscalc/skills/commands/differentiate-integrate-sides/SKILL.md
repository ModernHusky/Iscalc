---
name: differentiate-integrate-sides
description: 对等式两边求导或积分。
match_rules:
- (?i)deriv
- (?i)from
---

# differentiate-integrate-sides
> 对等式两边求导或积分。
## 使用时机
在处理参数积分或已证明的等式时，有时需要对等式两边同时求导或积分。

**重要前提**：这些规则只能在 `from <name>:` 引入的 CALCULATE 状态下使用（即从一个已证明的 subgoal 开始）。
## 指令
### 命令



### 两边求导

```
differentiate both sides at <var>
```

**示例**：
假设 `subgoal 1: I(a) = pi / (2 * a) for a > 0` 已证明。
```
from 1:
    differentiate both sides at a
    ...
```
结果：`(D a. I(a)) = D a. (pi / (2 * a))`

### 两边积分

```
integrate both sides
```

**示例**：
假设 `subgoal 1: (D x. I(x)) = f(x)` 已证明。
```
from 1:
    integrate both sides
    ...
```
结果：`I(x) = INT x. f(x)`

### 注意事项

- **必须使用 `from <name>:`**：不能直接在 `lhs:` 或 `rhs:` 块中使用这些规则。
- **语法陷阱**：`differentiate both sides with respect to a` 是**错误的**，正确写法是 `differentiate both sides at a`。
- 如果提示 "Applying the rule has no effect"，检查是否正确使用了 `from:` 进入状态。
