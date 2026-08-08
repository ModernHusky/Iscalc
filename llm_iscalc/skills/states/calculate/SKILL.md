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

> **【⚠️严肃警告：必须先加载技能】**
> 虽然本页面列出了许多可用命令（如 substitute, apply integral identity 等），但**绝对禁止**你直接根据这里的简略提示输出这些命令！
> 在你决定使用某个具体的高级命令前，你**必须**先通过 `<|load_skill|>对应的技能名称<|end_load_skill|>` 加载该命令的专属技能文档，仔细阅读其详细语法和约束后，才能输出对应的 JSON 命令。
### 常见一般命令

- **使用计算规则**: `rewrite`, `simplify`, `substitute`, `integrate by parts` 等
- **完成计算**: `done` - 当表达式已达到闭合形式时使用
- **切换到右边**: `rhs:` - 在证明等式时，完成左边计算后切换到右边

#### 高级求解策略

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
   - **相关技能**：
     - `<|load_skill|>improper-integral<|end_load_skill|>` — 广义积分转极限

3. **偶函数对称性**（高斯积分关键技巧）
   - **如果被积函数是偶函数** (即 `f(-x) = f(x)`)
   - **策略**: 利用 `INT x:[-oo,oo]. f(x) = 2 * INT x:[0,oo]. f(x)`
   - **步骤**:
     1. 先用 `split region at 0` 将积分拆分
     2. 对左半部分用 `substitute u for -x` 转换为右半部分
     3. 用 `rewrite` 合并两个相同的积分
   - **示例**: `exp(-(a*x^2))` 是偶函数
   - **相关技能**：
     - `<|load_skill|>split-region-at<|end_load_skill|>` — 区间拆分技巧
     - `<|load_skill|>substitute<|end_load_skill|>` — 换元积分详细语法和示例

4. **如果是复合函数，考虑换元**
   - **命令**: `substitute <u> for <u_expr>`
   - **示例**: `substitute u for x + 1`
   - **缩放换元**（处理 `exp(-(a*x^2))`）:
     - `substitute u for sqrt(a) * x`
     - 效果: 将 `exp(-(a*x^2))` 变为 `exp(-(u^2)) / sqrt(a)`
   - **相关技能**：
     - `<|load_skill|>substitute<|end_load_skill|>` — 换元积分详细语法和示例

5. **如果是乘积形式，考虑分部积分**
   - **命令**: `integrate by parts with u = <u_part>, v = <v_part>`
   - **示例**: `integrate by parts with u = x, v = exp(x)`
   - **相关技能**：
     - `<|load_skill|>integrate-by-parts<|end_load_skill|>` — 分部积分详细语法

6. **如果是有理函数，考虑部分分式分解**
   - **命令**: `rewrite <expr> to <partial_fraction_form>`
   - **提示**: 使用 `<|load_skill|>partial-fraction<|end_load_skill|>` 加载具体的分解技巧，然后手动构造 rewrite 目标。
   - **注意**: 系统**没有** `partial-fraction` 命令，必须使用 `rewrite`。
   - **相关技能**：
     - `<|load_skill|>partial-fraction<|end_load_skill|>` — 部分分式分解策略

7. **尝试应用积分恒等式**
   - **⚠️重要提示**：必须判断当前的积分表达式是否是公理型常见一般积分式，否则必须加载skill文档查看具体使用方法
   - **命令**: `apply integral identity`
   - **适用**: 基本函数的积分，如 `INT u. exp(-u^2)`
   - **相关技能**：
     - `<|load_skill|>apply-integral-identity<|end_load_skill|>` — 积分恒等式应用策略

8. **围道积分求解**
    - **命令**：`apply cintegral identity`
    - **适用**：围道积分，如 `CINT z:C(t,r). 1/(z^2+1)` 注意，标志是**CINT**

9. **多项式乘法展开与整理**
    - **命令**：`expand polynomial`
    - **适用**：多项式乘法展开与整理，如 `(x+1)^3`、`(x+1)(x-1)`
    - **相关技能**：
      - `<|load_skill|>expand-polynomial<|end_load_skill|>` — 多项式乘法展开与整理


### 工作流程

1. 使用 `calculate <expr>` 或 `lhs:` 进入此状态
2. 应用各种计算规则进行变换
3. 达到目标后使用 `done` 退出

### 注意事项

- 只有在此状态下才能使用计算规则
- 不能在此状态下使用 `calculate`（已经在计算中）
- 不能在此状态下使用 `lhs:` 或 `arg:`（需要先 `done`）

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


# state-calculate
> 计算状态操作指南。
## 使用时机
参考 description 描述。
## 指令
### 计算状态 (Calculate State)

在计算状态下，可以执行以下操作：

> **【⚠️严肃警告：必须先加载技能】**
> 虽然本页面列出了许多可用命令（如 substitute, apply integral identity 等），但**绝对禁止**你直接根据这里的简略提示输出这些命令！
> 在你决定使用某个具体的高级命令前，你**必须**先通过 `<|load_skill|>对应的技能名称<|end_load_skill|>` 加载该命令的专属技能文档，仔细阅读其详细语法和约束后，才能输出对应的 JSON 命令。
### 常见一般命令

- **使用计算规则**: `rewrite`, `simplify`, `substitute`, `integrate by parts` 等
- **完成计算**: `done` - 当表达式已达到闭合形式时使用
- **切换到右边**: `rhs:` - 在证明等式时，完成左边计算后切换到右边

#### 高级求解策略

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
   - **相关技能**：
     - `<|load_skill|>improper-integral<|end_load_skill|>` — 广义积分转极限

3. **偶函数对称性**（高斯积分关键技巧）
   - **如果被积函数是偶函数** (即 `f(-x) = f(x)`)
   - **策略**: 利用 `INT x:[-oo,oo]. f(x) = 2 * INT x:[0,oo]. f(x)`
   - **步骤**:
     1. 先用 `split region at 0` 将积分拆分
     2. 对左半部分用 `substitute u for -x` 转换为右半部分
     3. 用 `rewrite` 合并两个相同的积分
   - **示例**: `exp(-(a*x^2))` 是偶函数
   - **相关技能**：
     - `<|load_skill|>split-region-at<|end_load_skill|>` — 区间拆分技巧
     - `<|load_skill|>substitute<|end_load_skill|>` — 换元积分详细语法和示例

4. **如果是复合函数，考虑换元**
   - **命令**: `substitute <u> for <u_expr>`
   - **示例**: `substitute u for x + 1`
   - **缩放换元**（处理 `exp(-(a*x^2))`）:
     - `substitute u for sqrt(a) * x`
     - 效果: 将 `exp(-(a*x^2))` 变为 `exp(-(u^2)) / sqrt(a)`
   - **相关技能**：
     - `<|load_skill|>substitute<|end_load_skill|>` — 换元积分详细语法和示例

5. **如果是乘积形式，考虑分部积分**
   - **命令**: `integrate by parts with u = <u_part>, v = <v_part>`
   - **示例**: `integrate by parts with u = x, v = exp(x)`
   - **相关技能**：
     - `<|load_skill|>integrate-by-parts<|end_load_skill|>` — 分部积分详细语法

6. **如果是有理函数，考虑部分分式分解**
   - **命令**: `rewrite <expr> to <partial_fraction_form>`
   - **提示**: 使用 `<|load_skill|>partial-fraction<|end_load_skill|>` 加载具体的分解技巧，然后手动构造 rewrite 目标。
   - **注意**: 系统**没有** `partial-fraction` 命令，必须使用 `rewrite`。
   - **相关技能**：
     - `<|load_skill|>partial-fraction<|end_load_skill|>` — 部分分式分解策略

7. **尝试应用积分恒等式**
   - **⚠️重要提示**：必须判断当前的积分表达式是否是公理型常见一般积分式，否则必须加载skill文档查看具体使用方法
   - **命令**: `apply integral identity`
   - **适用**: 基本函数的积分，如 `INT u. exp(-u^2)`
   - **相关技能**：
     - `<|load_skill|>apply-integral-identity<|end_load_skill|>` — 积分恒等式应用策略

8. **围道积分求解**
    - **命令**：`apply cintegral identity`
    - **适用**：围道积分，如 `CINT z:C(t,r). 1/(z^2+1)` 注意，标志是**CINT**

9. **多项式乘法展开与整理**
    - **命令**：`expand polynomial`
    - **适用**：多项式乘法展开与整理，如 `(x+1)^3`、`(x+1)(x-1)`
    - **相关技能**：
      - `<|load_skill|>expand-polynomial<|end_load_skill|>` — 多项式乘法展开与整理


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
