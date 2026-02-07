---
name: strategy-limit
description: 极限求解通用策略。
match_rules:
- (?i)lim
---

# strategy-limit
> 极限求解通用策略。

## 使用时机

当遇到 `LIM {x -> a}. f(x)` 形式的极限表达式时使用。

## 极限求解基本策略

### 1. 直接代入法

首先尝试直接代入求值:
```
simplify
```
如果代入后得到确定值（非不定式），则求解完成。

### 2. 不定式类型识别与处理

#### 2.1 0/0 型或 ∞/∞ 型 → 洛必达法则
```
l'Hopital's rule
simplify
```

**适用场景**（来自 theories/lhopital.thy）:
- `LIM {x -> 1}. (x^2 - 1) / (x^2 + 3*x - 4)` → 分子分母同时趋向0
- `LIM {x -> 0}. sin(x) / x` → 经典 0/0 型
- `LIM {x -> 0}. (3^x - 2^x) / (x^2 - x)` → 指数函数的 0/0 型

#### 2.2 ∞ - ∞ 型 → 通分或对数合并

对于形如 `log(a) - log(b)` 的表达式，**先合并再求极限**：
```
rewrite log(a) - log(b) to log(a/b)
```

**为什么合并更优**：
- `log(x) - log(x+1)` → 合并为 `log(x/(x+1))`
- 当 `x → ∞` 时，`x/(x+1) → 1`，因此 `log(x/(x+1)) → log(1) = 0`
- 避免了 `∞ - ∞` 的不定式处理

> ⚠️ **注意：先完成所有 rewrite 再 simplify**
>
> `simplify` 会重排表达式项顺序（按复杂度排序），可能撤销你的 rewrite 结果。
> **正确做法**：在一个 rewrite 中完成对数合并，然后再 simplify。

**示例**（来自 theories/base.thy）:
```
axiom (LIM {x -> oo}. log(x / sqrt(x ^ 2 + 1))) = log(1)
```

#### 2.3 0 · ∞ 型 → 转换为分式

通过代数变换转为 0/0 或 ∞/∞ 型，然后用洛必达法则。

## 无穷极限处理技巧

### 分数重写模式

对于 `x → ∞` 的极限，常用以下重写技巧：

| 原形式 | 重写为 | 适用条件 |
|--------|--------|----------|
| `x/(x+1)` | `1 - 1/(x+1)` | 分子分母同阶 |
| `(x+1)/x` | `1 + 1/x` | 分子分母同阶 |
| `a/(a+b)` | `1 - b/(a+b)` | 通用形式 |
| `x^2/(x^2+1)` | `1 - 1/(x^2+1)` | 高阶情况 |

**示例命令**:
```
rewrite x/(x+1) to 1 - 1/(x+1)
simplify
```

### 对数表达式处理

#### 规则1: 负值加正值先调整顺序

当需要合并 `- log(a) + log(b)` 形式时，先调整为减法形式：
```
rewrite - log(t^2 + 1)/2 + log(t) to log(t) - log(t^2 + 1)/2
```

#### 规则2: 对数差合并为商的对数
```
rewrite log(a) - log(b) to log(a/b)
```

#### 规则3: 系数移入对数
```
rewrite n * log(x) to log(x^n)
rewrite log(x)/2 to log(sqrt(x))
```

**完整示例**:
```
# 原表达式: - log(t^2 + 1)/2 + log(t)
rewrite - log(t^2 + 1)/2 + log(t) to log(t) - log(t^2 + 1)/2
rewrite log(t) - log(t^2 + 1)/2 to log(t) - log(sqrt(t^2 + 1))
rewrite log(t) - log(sqrt(t^2 + 1)) to log(t / sqrt(t^2 + 1))
```

## 特殊极限技巧

### arctan 极限
```
# 当 x → ∞ 时
# arctan(x) → π/2
# arctan(-x) → -π/2
simplify
```

### 指数极限
```
# exp(-x) 当 x → ∞ 时趋向 0
# exp(x) 当 x → -∞ 时趋向 0
simplify
```

## 注意事项

1. **对数合并优先**：遇到 `log(a) - log(b)` 形式，优先合并为 `log(a/b)`，避免后续复杂计算。

2. **分数变换时机**：当极限变量趋向无穷时，分数重写（如 `x/(x+1) → 1 - 1/(x+1)`）通常能简化计算。

3. **洛必达法则限制**：仅适用于 0/0 或 ∞/∞ 型，使用前确认类型。

4. **simplify 收尾**：大多数极限操作后需要 `simplify` 来完成最终计算。

## 常见错误避免

- ❌ 直接对非不定式使用洛必达法则
- ❌ 忽略对数合并机会，导致步骤冗长
- ❌ 忘记处理负号顺序（`-a+b` 应先转为 `b-a`）

## 相关技能

在执行极限求解策略时，可能需要加载以下技能获取详细命令：

> 💡 `<|load_skill|>lhopital<|end_load_skill|>` — 洛必达法则详细用法和限制条件

> 💡 `<|load_skill|>rewrite<|end_load_skill|>` — 对数合并、分数变换等代数恒等变换

> 💡 `<|load_skill|>simplify<|end_load_skill|>` — 表达式简化规则
