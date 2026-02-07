---
name: state-induction
description: 归纳状态操作指南。
match_rules:
- induction
- SUM
- n
---

# state-induction
> 归纳状态操作指南。
## 使用时机
参考 description 描述。
## 指令
### 归纳状态 (Induction State)

在归纳状态下，可以执行以下操作：

### 可用命令

- `base:` - 进入基础情况 (Base Case)
- `induct:` - 进入归纳步骤 (Inductive Case)

### 工作流程

1. 使用 `induction on <var>` 开始归纳
2. 先用 `base:` 证明基础情况
3. 再用 `induct:` 证明归纳步骤
4. 两个分支都完成后使用 `done`

### 注意事项

- 目标必须是等式（当前唯一支持的情况）
- 归纳变量必须是整数
- 必须能证明变量至少达到起始值
- 在归纳步骤中，可以使用 `apply induction hypothesis` 应用归纳假设

### 完整示例 (来自 theories/standard.thy)



### 示例: 证明 ∫x^m*log(x)^n

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
        apply induction hypothesis (all)
        simplify
        rewrite to (-1) ^ (n + 1) * (m + 1) ^ (-n - 2) * ((n + 1) * factorial(n))
        rewrite (n + 1) * factorial(n) to factorial(n + 1)
        simplify
    done
done
```
