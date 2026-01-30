---
name: strategy-integral
description: 积分求解通用策略。
match_rules:
- (?i)int
---

# strategy-integral
> 积分求解通用策略。
## 使用时机
参考 description 描述。
## 指令
### 积分求解策略

1. **观察被积函数的形式**

2. **如果是复合函数，考虑换元**
   - **命令**: `substitute <u_expr> for <u>`
   - **示例**: `substitute x + 1 for u`

3. **如果是乘积形式，考虑分部积分**
   - **命令**: `integration_by_parts <u_part> <dv_part>`
   - **示例**: `integration_by_parts x exp(x)`

4. **如果是有理函数，考虑部分分式分解**
   - **命令**: `rewrite <expr> to <partial_fraction_form>`
   - **提示**: 使用 `<|load_skill|>partial-fraction<|end_load_skill|>` 加载具体的分解技巧，然后手动构造 rewrite 目标。
   - **注意**: 系统**没有** `partial-fraction` 命令，必须使用 `rewrite`。

5. **尝试应用积分恒等式**
   - **命令**: `apply integral identity`

6. **简化结果**
   - **命令**: `simplify`

### 常见换元技巧

- 根式: 令 u = sqrt(expr) -> `substitute sqrt(...) for u`
- 三角函数: 令 u = sin(x) -> `substitute sin(x) for u`
- 指数函数: 令 u = exp(x) -> `substitute exp(x) for u`
- 对数函数: 令 u = log(x) -> `substitute log(x) for u`
