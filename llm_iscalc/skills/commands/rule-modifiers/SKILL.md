---
name: rule-modifiers
description: 规则修饰符（at, on, for等）。
match_rules:
- .*
---

# rule-modifiers
> 规则修饰符（at, on, for等）。
## 使用时机
当表达式中有多个位置可以应用同一规则时，使用修饰符指定位置。
## 指令
### `(at n)` 修饰符

语法：`<rule> (at n)`，其中 `n` 从 1 开始。

**⚠️ 格式要求**：
- 必须使用**括号**：`(at n)` 而不是 `at n`
- 括号与规则之间有**空格**

**支持的规则**：
- `rewrite A to B (at n)`: 应用于第 n 次出现的 A。
- `substitute u for expr (at n)`: 应用于第 n 个积分。
- `substitute f(u) for x (at n)`: 逆换元，应用于第 n 个积分。
- `integrate by parts ... (at n)`: 应用于第 n 个积分。
- `split region at c (at n)`: 应用于第 n 个积分。
- `expand definition for f (at n)`: 应用于第 n 次出现的 f。

**不支持 `(at n)`**：
- `apply integral identity`: 自动尝试所有积分。
- `simplify`: 全局应用。

**错误示例对比**：

| ❌ 错误（导致解析失败） | ✅ 正确 |
|---------------------|--------|
| `substitute u for -x at 2` | `substitute u for -x (at 2)` |
| `rewrite A to B at 1` | `rewrite A to B (at 1)` |
| `substitute u for x(at 2)` | `substitute u for x (at 2)` |

### `(all)` 修饰符

语法：`<rule> (all)`，将规则应用于所有匹配位置。

**常用组合**：
- `rewrite A to B (all)`: 替换所有 A。
- `expand definition for f (all)`: 展开所有 f。
- `apply induction hypothesis (all)`: 应用归纳假设到所有匹配位置。

### 示例

表达式: `sqrt(2) + sqrt(2) * x`

- `rewrite sqrt(2) to 2^(1/2)`: 只替换第一个。
- `rewrite sqrt(2) to 2^(1/2) (at 2)`: 只替换第二个。
- `rewrite sqrt(2) to 2^(1/2) (all)`: 全部替换。
