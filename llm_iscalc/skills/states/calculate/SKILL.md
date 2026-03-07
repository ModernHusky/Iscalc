---
name: state-calculate
description: 计算状态操作指南。
match_rules:
- .*
---

# state-calculate
> 计算状态操作指南。
## 使用时机
参考 description 描述。
## 指令
### 计算状态 (Calculate State)

在计算状态下，可以执行以下操作：

### 可用命令

- **使用计算规则**: `rewrite`, `simplify`, `substitute`, `integrate by parts` 等
- **完成计算**: `done` - 当表达式已达到闭合形式时使用
- **切换到右边**: `rhs:` - 在证明等式时，完成左边计算后切换到右边

### 工作流程

1. 使用 `calculate <expr>` 或 `lhs:` 进入此状态
2. 应用各种计算规则进行变换
3. 达到目标后使用 `done` 退出

### 注意事项

- 只有在此状态下才能使用计算规则
- 不能在此状态下使用 `calculate`（已经在计算中）
- 不能在此状态下使用 `lhs:` 或 `arg:`（需要先 `done`）
- **绝不允许的证明手法**：在等式证明中进入 `lhs:` 或`rhs:`  计算时，如果你得到的表达式等于另一边的**子表达式**（例如，对 RHS 积分结果求导得到了 LHS 的被积函数），这是**绝对错误且无效**的！
  - 计算状态中的每一步必须保持该边表达式的值完全不变（如代数变形、计算积分）。如果通过求导改变了它的值，然后输入 `done`，系统回到 Prove 状态比对两边时会直接判你的证明失败！
  - `done` 的要求：如果在 `lhs:` 或 `rhs:` 中，你的计算结果必须与另一边**字面上完全吻合**（包括形式、积分常数等）。

### 完整示例 (来自 theories/)



### 示例: 标准计算流程

```
calculate INT x:[0, 1]. x*exp(x)
    integrate by parts with u = x, v = exp(x)
    apply integral identity
    simplify
done
```

### 示例: 多步骤计算

```
calculate INT x:[0, 1]. 2*x*arctan(x)
    integrate by parts with u = arctan(x), v = x^2
    simplify
    rewrite x ^ 2 / (x ^ 2 + 1) to 1 - 1 / (x ^ 2 + 1)
    apply integral identity
    simplify
done
```
