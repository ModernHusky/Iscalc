---
name: improper-integral
description: 广义积分处理 - 将含无穷的积分转换为极限形式。当积分上下界含有无穷时使用。
keywords:
  - 广义积分
  - 无穷
  - 极限
  - improper
applicable_types:
  - improper_integral
match_rules:
  - "oo"
  - "inf"
---

# 何时使用

当积分上下界含有无穷：
- `INT x:[0,oo]. f(x)`
- `INT x:[-oo,0]. f(x)`
- `INT x:[-oo,oo]. f(x)`

# 快速开始

```
improper integral to limit creating t
```

将广义积分转换为 `LIM {t -> oo}. INT x:[a,t]. f(x)`

# 分步指南

1. 识别积分含有无穷上下界
2. 执行 `improper integral to limit creating t`
3. 执行 `apply integral identity` 计算积分
4. 执行 `simplify` 计算极限

# 示例

## 示例: 指数衰减积分
表达式: `INT x:[0,oo]. exp(-x)`

```
improper integral to limit creating t
```

转换为: `LIM {t -> oo}. INT x:[0,t]. exp(-x)`

继续:
```
apply integral identity
simplify
```

最终结果: 1

# 注意事项

- 变量名 `t` 可自定义，避免与被积变量冲突
- 转换后需要依次执行积分和极限计算
- 积分可能发散（结果为无穷）
