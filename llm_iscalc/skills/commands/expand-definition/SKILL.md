---
name: expand-definition
description: 展开函数或变量的定义。
match_rules:
- .*
---

# expand-definition
> 展开函数或变量的定义。

## 指令

### 语法格式

```
expand definition for <name>
```

展开单个出现的定义。

```
expand definition for <name> (all)
```

展开所有出现的该定义。

### 参数说明

- `<name>`: 需要展开的函数或符号名称

## 使用时机

当表达式中包含之前定义的函数（如 `let I(a) = ...`）或公理化定义的符号，且需要将其展开以便进行具体的代数运算或积分时：
- 展开用户定义的函数 `I(a)`
- 展开特殊函数定义如 `cosh`, `sinh`, `erf` 等
- 在费曼技巧中展开参数化积分

## 适用目标类型

| 目标类型 | 适用性 |
|---------|--------|
| 包含 `let` 定义函数的表达式 | ✅ |
| 包含特殊函数的表达式 | ✅ |
| 公理化定义（无具体表达式）| ❌ |

## 工作流程

1. **查找定义**：在上下文中查找 `<name>` 的定义
2. **条件检查**：验证当前参数是否满足定义条件（如 `for t > 0`）
3. **模式匹配**：将 `<name>(<args>)` 与定义模式匹配，提取参数绑定
4. **替换执行**：用定义体（代入参数后）替换函数调用

## 注意事项

- **参数约束**：如果定义时有约束（如 `for t > 0`），而当前环境不满足该约束，展开可能会失败，提示 "Applying the rule has no effect"。
- **公理符号**：如果符号是公理化定义的（没有具体表达式），此规则无效。请尝试 `rewrite`。
- **全部展开**：如果要展开表达式中出现的所有该符号，可加 `(all)` 修饰符。

## 示例

### 示例1: 展开用户定义函数

假设之前定义了 `let I(t) = INT x. x^t`。
当前表达式: `D t. I(t)`

```
expand definition for I
```

结果: `D t. INT x. x^t`

### 示例2: 展开 cosh 定义 (来自 theories/)

```
prove (INT x:[0,oo]. 1 / (x ^ 4 + 2 * x ^ 2 * cosh(2 * a) + 1)) = pi / (4 * cosh(a))
lhs:
    expand definition for cosh (all)
    # cosh(x) 展开为 (exp(x) + exp(-x)) / 2
    rewrite x ^ 4 + 2 * x ^ 2 * ((exp(-(2 * a)) + exp(2 * a)) / 2) + 1 to ...
    ...
done
```

### 示例3: 配合 fold definition (来自 theories/)

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
