---
name: strategy-limit
description: 极限求解策略指南
match_rules:
  - "(?i)lim"
---

## 极限求解策略
1. 直接代入尝试求值
2. 如果是0/0或∞/∞型，考虑洛必达法则
3. 对于a/(a+b)形式，重写为1-b/(a+b)
4. 使用simplify化简结果

### 常见分数重写模式
- x/(x+1) → 1 - 1/(x+1)  当 x → ∞ 时
- (x+1)/x → 1 + 1/x  当 x → ∞ 时
