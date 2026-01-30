---
name: merge-evalat
description: 求值项合并策略（EvalAt）。
match_rules:
- regex: .*\]_.*\s*[\+\-]\s*\[.*\].*
  flags:
  - re.DOTALL
---

# merge-evalat
> 求值项合并策略（EvalAt）。
## 使用时机
参考 description 描述。
## 指令
### 合并 EvalAt 表达式

当遇到多个共享相同变量和上下限的求值项 (EvalAt) 时，例如：

`[f(x)]_x=a,b - [g(x)]_x=a,b`

应使用 `rewrite` 命令将它们合并为单个求值项。

> [!WARNING]
> **绝对不要使用 `simplify` 命令来尝试此合并**。`simplify` 在处理此类合并时极其容易产生幻觉或失败。必须明确使用 `rewrite to ...`（对整个表达式进行rewrite）或 `rewrite <expr> to <target>`（对指定表达式进行rewrite）。

`rewrite to ([f(x) - g(x)]_x=a,b)` 或 `rewrite [f(x)]_x=a,b - [g(x)]_x=a,b to ([f(x) - g(x)]_x=a,b)`

这种简化可以揭示进一步的消去机会，或者使表达式更易于计算。

### 示例



### 示例 1: 有限限

**表达式**: `[log(u)]_u=1,t + [log(u + 1)]_u=1,t`
**动作**:
```json
{
    "thinking": "两项都是从 u=1 到 t 的求值。利用求值的线性性质，我可以将它们合并。参考 skills/strategies/merge-evalat/SKILL.md",
    "command": "rewrite to ([log(u) + log(u + 1)]_u=1,t)",
    "explanation": "将两个求值项合并为一个。",
    "is_final": false
}
```

### 示例 2: 无穷限

**表达式**: `[1/x]_x=1,oo + [1/(x+1)]_x=1,oo`
**动作**:
```json
{
    "thinking": "两个求值项共享相同的上下限，例如 [1, oo]。合并它们可能有帮助。参考 skills/strategies/merge-evalat/SKILL.md",
    "command": "rewrite to ([1/x + 1/(x+1)]_x=1,oo)",
    "explanation": "将各项合并到同一个求值符号下。",
    "is_final": false
}
```

### 示例 3： 子表达式求值

{
  "thinking": "两个子求值表达式项共享相同的上下限，例如 [1, t]。合并它们可能有帮助。参考 skills/strategies/merge-evalat/SKILL.md",
    "command": "rewrite ([log(x)]_x=1,t) - ([log(x ^ 2 + 1) / 2]_x=1,t) to ([log(x) - log(x ^ 2 + 1) / 2]_x=1,t)",
    "explanation": "将各项合并到同一个求值符号下。",
    "is_final": false
}

### 策略

1. **识别**: 寻找多个具有相同上下限的 `[...]` 项。
2. **构建**: 创建一个目标表达式，其中函数体合并在同一个 `[...]` 中。
3. **执行**: 使用 `rewrite to <target>` 命令。
