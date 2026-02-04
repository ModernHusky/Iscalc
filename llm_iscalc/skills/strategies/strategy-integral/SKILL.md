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

2. **广义积分：先转换为极限**
   - **如果积分限包含无穷** (如 `[-oo, oo]`, `[0, oo]`)
   - **命令**: `improper integral to limit creating t`
   - **效果**: `INT x:[a,oo]. f(x)` → `LIM {t -> oo}. INT x:[a,t]. f(x)`
   - **示例**: 
     ```
     INT x:[-oo,oo]. exp(-(a * x ^ 2))
     → LIM {t -> oo}. INT x:[-t,t]. exp(-(a * x ^ 2))
     ```

3. **偶函数对称性**（高斯积分关键技巧）
   - **如果被积函数是偶函数** (即 `f(-x) = f(x)`)
   - **策略**: 利用 `INT x:[-oo,oo]. f(x) = 2 * INT x:[0,oo]. f(x)`
   - **步骤**:
     1. 先用 `split region at 0` 将积分拆分
     2. 对左半部分用 `substitute u for -x` 转换为右半部分
     3. 用 `rewrite` 合并两个相同的积分
   - **示例**: `exp(-(a*x^2))` 是偶函数

4. **如果是复合函数，考虑换元**
   - **命令**: `substitute <u> for <u_expr>`
   - **示例**: `substitute u for x + 1`
   - **缩放换元**（处理 `exp(-(a*x^2))`）:
     - `substitute u for sqrt(a) * x`
     - 效果: 将 `exp(-(a*x^2))` 变为 `exp(-(u^2)) / sqrt(a)`

5. **如果是乘积形式，考虑分部积分**
   - **命令**: `integrate by parts with u = <u_part>, v = <v_part>`
   - **示例**: `integrate by parts with u = x, v = exp(x)`

6. **如果是有理函数，考虑部分分式分解**
   - **命令**: `rewrite <expr> to <partial_fraction_form>`
   - **提示**: 使用 `<|load_skill|>partial-fraction<|end_load_skill|>` 加载具体的分解技巧，然后手动构造 rewrite 目标。
   - **注意**: 系统**没有** `partial-fraction` 命令，必须使用 `rewrite`。

7. **尝试应用积分恒等式**
   - **命令**: `apply integral identity`
   - **适用**: 基本函数的积分，如 `INT u. exp(-u^2)`

8. **简化结果**
   - **命令**: `simplify`

### 常见换元技巧

- 根式: 令 u = sqrt(expr) -> `substitute sqrt(...) for u`
- 三角函数: 令 u = sin(x) -> `substitute sin(x) for u`
- 指数函数: 令 u = exp(x) -> `substitute exp(x) for u`
- 对数函数: 令 u = log(x) -> `substitute log(x) for u`

## 相关技能

> 💡 `<|load_skill|>substitute<|end_load_skill|>` — 换元积分详细语法和示例

> 💡 `<|load_skill|>integrate-by-parts<|end_load_skill|>` — 分部积分详细语法

> 💡 `<|load_skill|>improper-integral<|end_load_skill|>` — 广义积分转极限

> 💡 `<|load_skill|>rewrite<|end_load_skill|>` — 代数变换技巧

> 💡 `<|load_skill|>partial-fraction<|end_load_skill|>` — 部分分式分解策略
