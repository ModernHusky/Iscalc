---
name: subgoal-definition
description: 子目标与定义 - 使用 subgoal 和 let 构建复杂证明。
match_rules:
  - "(?i)let"
  - "(?i)subgoal"
  - "(?i)define"
applicable_types:
  - proof
  - general
---

# 何时使用

当证明目标非常复杂，需要拆分为多个中间步骤，或者需要引入参数化函数 `I(t)` 来简化问题时。

# 定义参数化函数 (`let`)

语法：`let <expr> for <conditions>`

**示例**：
```
let I(t) = INT x:[0, 1]. (x ^ t - 1) / log(x) for t != -1
```

定义后，可以：
- `expand definition for I`: 将 `I(t)` 展开为其定义的积分表达式。
- `fold definition for I`: 将匹配的积分表达式折叠为 `I(t)`。

# 创建子目标 (`subgoal`)

语法：`subgoal <name>: <goal> for <conditions>`

**示例**：
```
subgoal 1: (D t. I(t)) = 1 / (t + 1) for t != -1
```

**从子目标开始证明** (`from`):
```
from 1:
    differentiate both sides at t
    ...
```

# 复杂证明模式

对于需要多步推导的积分（如参数微分法），典型流程为：

1.  定义 `I(t)` (参数化积分)。
2.  创建 `subgoal 1` 证明 `D t. I(t) = ...`。
3.  创建 `subgoal 2` 证明初值 `I(0) = ...`。
4.  创建 `subgoal 3`，使用 `from 1:` 对导数表达式积分，得到 `I(t)` 的显式形式。
5.  创建 `subgoal 4`，利用 `subgoal 2` 和 `subgoal 3` 确定积分常数。
6.  最终用 `fold definition for I` 和 `apply ? on I(?)` 得到答案。

**参见**: `dsl.md` 中的 Example 6 (Construct subgoal for proving a complex goal) 提供了完整示例。

# 注意事项

- **依赖关系**：子目标必须按依赖顺序证明（被依赖的先证）。
- **应用子目标**: 使用 `apply <name> on <expr>` 将已证明的子目标应用到当前表达式。
