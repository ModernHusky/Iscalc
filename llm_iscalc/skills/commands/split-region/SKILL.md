---
name: split-region
description: 分裂积分区间 - 将定积分按某点分成两个积分的和。
match_rules:
  - "(?i)int.*:\\["
applicable_types:
  - definite_integral
---

# 何时使用

当需要在某个特定点 `c` 将积分区间 $[a, b]$ 分成 $[a, c] + [c, b]$ 时使用。常见场景：
- 被积函数在 `c` 处有分段定义。
- 需要单独处理某个区间（如利用对称性）。

# 快速开始

```
split region at <c>
```

# 示例

表达式: `INT x:[0, 2]. |x - 1|`

```
split region at 1
```

结果: `INT x:[0, 1]. |x - 1| + INT x:[1, 2]. |x - 1|`

# 注意事项

- `c` 必须在积分区间 $(a, b)$ 内。
- 分裂后需要分别处理两个积分。

# 完整示例 (来自 theories/)

## 示例: 利用对称性
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

## 示例: 无穷区间分割
```
subgoal 3: (INT x:[0,oo]. 1/(x^4+...)) = 1/4*(INT x:[-oo,oo]. ...)
rhs:
    split region at 0
    substitute u for -x
    ...
done
```
