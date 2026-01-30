---
name: fold-definition
description: 将表达式折叠为函数定义形式。
match_rules:
- (?i)let
- (?i)fold
---

# fold-definition
> 将表达式折叠为函数定义形式。
## 使用时机
当需要将一个复杂表达式替换为之前使用 `let` 定义的函数名时。这是 `expand definition` 的逆操作。
## 指令
### 快速开始

```
fold definition for <name>
```

### 示例

假设之前定义了：
```
let I(t) = INT x:[0, 1]. (x ^ t - 1) / log(x) for t != -1
```

当前表达式: `... INT x:[0, 1]. (x ^ 1 - 1) / log(x) ...`

```
fold definition for I
```

结果: `... I(1) ...`

### 注意事项

- 表达式必须**精确匹配**定义中的结构（可能需要先 `rewrite` 或 `simplify`）。
- 如果定义有条件约束（如 `for t != -1`），确保当前参数满足条件。
