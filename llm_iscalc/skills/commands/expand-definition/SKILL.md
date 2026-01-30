---
name: expand-definition
description: 展开函数或变量的定义。
match_rules:
- .*
---

# expand-definition
> 展开函数或变量的定义。
## 使用时机
当表达式中包含之前定义的函数（如 `let I(a) = ...`）或公理化定义的符号，且需要将其展开以便进行具体的代数运算或积分时。
## 指令
### 快速开始

```
expand definition for <name>
```

### 示例

假设之前定义了 `let I(t) = INT x. x^t`。
当前表达式: `D t. I(t)`

```
expand definition for I
```

结果: `D t. INT x. x^t`

### 注意事项

- **参数约束**：如果定义时有约束（如 `for t > 0`），而当前环境不满足该约束，展开可能会失败，提示 "Applying the rule has no effect"。
- **公理符号**：如果符号是公理化定义的（没有具体表达式），此规则无效。请尝试 `rewrite`。
- **全部展开**：如果要展开表达式中出现的所有该符号，可加 `(all)` 修饰符：`expand definition for f (all)`。

### 完整示例 (来自 theories/)



### 示例: 展开 cosh 定义

```
prove (INT x:[0,oo]. 1 / (x ^ 4 + 2 * x ^ 2 * cosh(2 * a) + 1)) = pi / (4 * cosh(a))
lhs:
    expand definition for cosh (all)
    rewrite x ^ 4 + 2 * x ^ 2 * ((exp(-(2 * a)) + exp(2 * a)) / 2) + 1 to ...
    ...
done
```

### 示例: 配合 fold definition

```
prove (INT x:[-oo,oo]. 1 / cosh(x)) = pi
lhs:
    expand definition for cosh (all)
    substitute t for exp(x)
    rewrite t * (1 / t + t) to 1 + t ^ 2
    apply integral identity
    simplify
done
```
