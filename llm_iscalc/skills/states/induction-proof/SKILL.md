---
name: state-induction-proof
description: 归纳证明 - 使用 induction on 进行的证明
match_rules:
  - "SUM"
  - "n"
  - "induction"
applicable_types:
  - general
---

# 归纳证明 (Proof by Induction)

使用 `induction on <var>` 开始的证明方式。

## 前提条件

- 目标必须是**等式**（当前唯一支持的情况）
- 归纳变量必须是**整数**
- 必须能证明变量至少达到起始值

## 语法

```
induction on n
induction on n starting from 0
induction on n starting from 1
```

## 分支处理

进入归纳状态后：
1. `base:` - 证明基础情况 (n = 起始值)
2. `induct:` - 证明归纳步骤 (n -> n+1)

## 归纳假设

在归纳步骤中，可以使用：
```
apply induction hypothesis
```
来应用归纳假设。

## 完成

两个分支都完成后，使用 `done` 退出归纳证明。

# 完整示例 (来自 theories/standard.thy)

## 示例: 归纳假设的应用
```
prove (INT x:[0,1]. x ^ m * log(x) ^ n) = (-1) ^ n * factorial(n) / (m + 1) ^ (n + 1) for m >= 0, n >= 0, isInt(n)
induction on n
base:
    lhs:
        apply integral identity
        simplify
    done
induct:
    lhs:
        integrate by parts with u = log(x) ^ (n + 1), v = x ^ (m + 1) / (m + 1)
        simplify
        apply induction hypothesis (all)    <-- 应用归纳假设
        simplify
        rewrite to (-1) ^ (n + 1) * (m + 1) ^ (-n - 2) * ((n + 1) * factorial(n))
        rewrite (n + 1) * factorial(n) to factorial(n + 1)
        simplify
    done
done
```
