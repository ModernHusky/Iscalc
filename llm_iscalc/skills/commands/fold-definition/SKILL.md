---
name: fold-definition
description: 将表达式折叠为函数定义形式。
match_rules:
- (?i)let
- (?i)fold
---

# fold-definition
> 将表达式折叠为函数定义形式。

## 指令

### 语法格式

```
fold definition for <name>
```

### 参数说明

- `<name>`: 之前用 `let` 定义的函数名称

## 使用时机

当需要将一个复杂表达式替换为之前使用 `let` 定义的函数名时：
- 这是 `expand definition` 的逆操作
- 用于简化最终结果的表示
- 在递归证明中引用定义

## 适用目标类型

| 目标类型 | 适用性 |
|---------|--------|
| 包含与 `let` 定义匹配的子表达式 | ✅ |
| 任意表达式 | ⚠️ 需精确匹配定义 |

## 工作流程

1. **查找定义**：在上下文中查找 `<name>` 的 `let` 定义
2. **模式匹配**：尝试将当前表达式的子部分与定义体匹配
3. **参数推断**：从匹配中推断函数的参数值
4. **条件检查**：验证推断的参数是否满足定义条件
5. **替换执行**：将匹配的表达式替换为函数调用 `<name>(<args>)`

## 注意事项

- **精确匹配**：表达式必须**精确匹配**定义中的结构。可能需要先使用 `rewrite` 或 `simplify` 调整形式。
- **条件约束**：如果定义有条件约束（如 `for t != -1`），确保当前参数满足条件。
- **与 expand 对称**：`fold` 和 `expand` 是互逆操作。

## 示例

### 示例1: 基本折叠

假设之前定义了：
```
let I(t) = INT x:[0, 1]. (x ^ t - 1) / log(x) for t != -1
```

当前表达式: `... INT x:[0, 1]. (x ^ 1 - 1) / log(x) ...`

```
fold definition for I
```

结果: `... I(1) ...`

### 示例2: 参数化折叠 (来自 theories/)

```
let I(a) = INT x:[0,oo]. exp(-a * x^2) for a > 0

subgoal 2: ...
lhs:
    ...
    # 当前表达式包含 INT x:[0,oo]. exp(-x^2)
    fold definition for I
    # 结果：I(1)
done
```

### 示例3: 递归引用

```
let I(n) = INT x:[0,1]. x^n * log(x)^2 for n : int, n >= 0

calculate I(3)
    expand definition for I
    integrate by parts with u = log(x)^2, v = x^4/4
    ...
    fold definition for I   # 将出现的积分识别为 I(4) 形式
    ...
done
```
