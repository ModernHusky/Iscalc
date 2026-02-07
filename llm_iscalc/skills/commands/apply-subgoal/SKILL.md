---
name: apply-subgoal
description: 应用已证明的子目标到当前表达式。
match_rules:
- (?i)apply
- (?i)subgoal
---

# apply-subgoal
> 应用已证明的子目标到当前表达式。

## 指令

### 语法格式

```
apply <subgoal_name> on <expr>
```

### 参数说明

- `<subgoal_name>`: 已证明的子目标名称（数字或标识符）
- `<expr>`: 当前表达式中需要被替换的部分，必须匹配子目标左侧的模式

## 使用时机

当需要使用之前证明的 subgoal 来替换或化简当前表达式中的某部分时。典型场景：
- 多步证明中引用前置结论
- 利用辅助引理简化计算
- 参数化公式的实例化应用

## 适用目标类型

| 目标类型 | 适用性 |
|---------|--------|
| 任意表达式 | ✅ |
| 等式证明 `lhs:/rhs:` | ✅ |
| 计算流程 `calculate` | ✅ |

## 工作流程

1. **查找子目标**：在已证明的 subgoal 列表中查找指定名称
2. **模式匹配**：将 `<expr>` 与子目标左侧进行模式匹配，提取参数绑定
3. **条件检查**：验证当前环境是否满足子目标的前置条件（`for` 子句）
4. **替换执行**：用子目标右侧（代入参数后）替换 `<expr>`

## 注意事项

- **subgoal 必须已证明**：只能引用当前步骤之前已 `done` 的子目标。
- **精确匹配**：`<expr>` 必须能够匹配子目标的左侧模式。
- **与 rewrite 的区别**：`apply` 用于已证明的恒等式，`rewrite` 用于代数变换。

## 示例

### 示例1: 基本应用

假设已证明：
```
subgoal 1: I(a) = pi / (2 * a) for a > 0
```

当前表达式: `... I(1) ...`

```
apply 1 on I(1)
```

结果: `... pi / 2 ...`

### 示例2: 参数化应用 (来自 theories/)

```
subgoal 1: (INT x:[0,oo]. exp(-a*x^2)) = sqrt(pi)/(2*sqrt(a)) for a > 0
lhs:
    ...
done

subgoal 2: (INT x:[0,oo]. x^2 * exp(-x^2)) = sqrt(pi)/4
lhs:
    ...
    apply 1 on INT x:[0,oo]. exp(-x^2)   # 这里 a=1
    simplify
done
```
