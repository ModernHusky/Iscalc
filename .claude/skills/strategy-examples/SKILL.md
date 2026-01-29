---
name: strategy-examples
description: 常见数学问题的解题示例库
version: "1.0.0"
author: "iscalc-team"
match_rules: [] 
# 示例通常通过 prompts.py 的显式逻辑加载，或者可以添加规则
---

## 示例

### 示例1: 简单换元
表达式: INT x:[0,1]. (x+1)^3
步骤:
1. substitute u for x + 1
2. apply integral identity
3. simplify
4. replace substitution
5. simplify

### 示例2: 分部积分
表达式: INT x:[0,1]. x * exp(x)
步骤:
1. integrate by parts with u = x, v = exp(x)
2. simplify
3. apply integral identity
4. simplify

### 示例3: 部分分式
表达式: INT x:[0,1]. 1/(x^2 - 1)
步骤:
1. partial fraction decomposition
2. simplify
3. apply integral identity
4. simplify
