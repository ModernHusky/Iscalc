---
name: strategy-integral
description: 定积分求解策略指南
match_rules:
  - "(?i)int"
---

## 定积分求解策略
1. 观察被积函数的形式
2. 如果是复合函数，考虑换元
3. 如果是乘积形式，考虑分部积分
4. 如果是有理函数，考虑部分分式分解
5. 尝试应用积分恒等式
6. 简化结果

### 常见换元技巧
- 根式: 令 u = sqrt(expr)
- 三角函数: 令 u = sin(x) 或 u = tan(x/2)
- 指数函数: 令 u = exp(x)
- 对数函数: 令 u = log(x)
