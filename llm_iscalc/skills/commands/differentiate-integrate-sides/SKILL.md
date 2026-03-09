---
name: differentiate-integrate-sides
description: 对等式两边求导或积分。
match_rules:
- (?i)deriv
- (?i)from
---

# differentiate-integrate-sides
> 对等式两边求导或积分。

## 指令

### 语法格式

**两边求导**：
```
differentiate both sides at <var>
```

**两边积分**：
```
integrate both sides
```

### 参数说明

- `<var>`: 求导变量（仅对 `differentiate both sides` 命令）

## 使用时机

在处理参数积分或已证明的等式时，需要对等式两边同时求导或积分的场景：
- 费曼积分技巧（Feynman's trick）：对参数求导
- 从导数方程恢复原函数
- 微分方程的求解过程

**重要前提**：这些规则只能在 `from <id>:` 引入的 CALCULATE 状态下使用（即从一个已证明的 subgoal 开始）。

## 适用目标类型

| 目标类型 | 适用性 |
|---------|--------|
| `from <id>:` 状态 | ✅ |
| `lhs:` / `rhs:` 状态 | ❌ |
| `calculate` 状态 | ❌ |

## 工作流程

### 两边求导 `differentiate both sides at <var>`

1. **获取等式**：从当前状态获取等式 `L = R`
2. **对两边求导**：生成 `(D <var>. L) = (D <var>. R)`
3. **化简**（可选）：后续通常使用 `simplify` 进行求导计算

### 两边积分 `integrate both sides`

1. **获取等式**：从当前状态获取等式 `L = R`（通常 L 包含导数形式）
2. **对两边积分**：生成对应的积分等式
3. **化简**（可选）：后续通常使用 `simplify` 或 `apply integral identity`

## 注意事项

- **必须使用 `from <id>:`**：不能直接在 `lhs:` 或 `rhs:` 块中使用这些规则。
- **语法陷阱**：`differentiate both sides with respect to a` 是**错误的**，正确写法是 `differentiate both sides at a`。
- **无效操作**：如果提示 "Applying the rule has no effect"，检查是否正确使用了 `from:` 进入状态。

## 示例

### 示例1: 两边求导（费曼技巧）

假设 `subgoal 1: I(a) = pi / (2 * a) for a > 0` 已证明。

```
from 1:
    differentiate both sides at a
    # 结果：(D a. I(a)) = D a. (pi / (2 * a))
    simplify
    # 结果：I'(a) = -pi / (2 * a^2)
done
```

### 示例2: 两边积分

假设 `subgoal 1: (D x. I(x)) = f(x)` 已证明。

```
from 1:
    integrate both sides
    # 结果：I(x) = INT x. f(x)
    apply integral identity
    simplify
done
```

### 示例3: 完整工作流程 (来自 theories/)

```
prove (INT x:[0,oo]. x^2 * exp(-a*x^2)) = sqrt(pi)/(4*a^(3/2)) for a > 0
subgoal 1: (INT x:[0,oo]. exp(-a*x^2)) = sqrt(pi)/(2*sqrt(a)) for a > 0
lhs:
    ...
done

from 1:
    differentiate both sides at a
    expand definition for I
    exchange derivative and integral
    simplify
done
```

