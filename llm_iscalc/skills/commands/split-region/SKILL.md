---
name: split-region
description: 拆分积分区间。
match_rules:
- (?i)int.*:\[
---

# split-region
> 拆分积分区间。

## 指令

### 语法格式

```
split region at <c>
```

### 参数说明

- `<c>`: 分割点，必须在积分区间 `(a, b)` 内

## 使用时机

当需要在某个特定点 `c` 将积分区间 `[a, b]` 分成 `[a, c] + [c, b]` 时使用：
- 被积函数在 `c` 处有分段定义
- 需要单独处理某个区间（如利用对称性）
- 无穷区间需要分割后分别处理

## 适用目标类型

| 目标类型 | 适用性 |
|---------|--------|
| 定积分 `INT x:[a,b]. f(x)` | ✅ |
| 广义积分 `INT x:[-oo,oo]. f(x)` | ✅ |
| 不定积分 | ❌ |

## 工作流程

1. **识别积分**：找到当前表达式中的定积分
2. **验证分割点**：确认 `c` 在区间 `(a, b)` 内
3. **拆分区间**：生成两个积分的和
   - `INT x:[a,b]. f(x)` → `INT x:[a,c]. f(x) + INT x:[c,b]. f(x)`
4. **保持被积函数**：两个新积分的被积函数与原积分相同

**数学原理**：
积分的区间可加性：$\int_a^b f(x)dx = \int_a^c f(x)dx + \int_c^b f(x)dx$

## 注意事项

- **分割点范围**：`c` 必须在积分区间 `(a, b)` 内。
- **后续处理**：分裂后需要分别处理两个积分。
- **对称性利用**：常用于分割后对某个区间进行换元，利用对称性简化。

## 示例

### 示例1: 分段函数

表达式: `INT x:[0, 2]. |x - 1|`

```
split region at 1
```

结果: `INT x:[0, 1]. |x - 1| + INT x:[1, 2]. |x - 1|`

### 示例2: 利用对称性 (来自 theories/)

```
prove (INT x:[0,pi / 2]. log(a * sin(2*x))) = (INT x:[0,pi / 2]. log(a * sin(x)))
lhs:
    substitute t for 2*x
    simplify
    split region at pi/2
    simplify
    substitute u for pi-t (at 2)
    rewrite sin(pi-u) to sin(u)
    ...
done
```

### 示例3: 无穷区间分割 (来自 theories/)

```
subgoal 3: (INT x:[0,oo]. 1/(x^4+...)) = 1/4*(INT x:[-oo,oo]. ...)
rhs:
    split region at 0
    # 结果：INT x:[-oo,0]. ... + INT x:[0,oo]. ...
    substitute u for -x
    # 利用偶函数性质
    ...
done
```

### 示例4: 完整工作流程

```
calculate INT x:[-1, 1]. |x|
    split region at 0
    # 结果：INT x:[-1,0]. |x| + INT x:[0,1]. |x|
    rewrite |x| to -x (at 1)    # 在 [-1,0] 上 |x| = -x
    rewrite |x| to x (at 2)      # 在 [0,1] 上 |x| = x
    apply integral identity
    simplify
done
```
