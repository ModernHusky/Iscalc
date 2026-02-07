---
name: induction-proof
description: 归纳证明模式操作指南。
match_rules:
- (?i)factorial
- (?i)sum.*n
- n\s*>?=\s*\d
---

# induction-proof
> 归纳证明模式操作指南。
## 使用时机
当需要证明一个关于整数变量 `n` 的等式（尤其是含有阶乘、求和等结构）时，使用数学归纳法。
## 指令
### 快速开始

```
induction on n
base:
    lhs:
        ... # 证明 n=0 (或起点) 时成立
    done
induct:
    lhs:
        ... # 假设 n=k 成立，证明 n=k+1 成立
        apply induction hypothesis (all)
        ...
    done
done
```

### 关键命令

1.  **开始归纳**: `induction on <var>` 或 `induction on <var> starting from <m>`
2.  **基础情况**: `base:` 后进入 PROVE 状态处理 n=m 的情况。
3.  **归纳步骤**: `induct:` 后进入 PROVE 状态处理 n=k+1 的情况，假设 n=k 成立。
4.  **应用归纳假设**: `apply induction hypothesis` 或 `apply induction hypothesis (all)`
    - 将当前表达式中匹配归纳假设左边的部分替换为右边。

### 示例：证明 $\int_0^\infty x^n e^{-x} dx = n!$

```
prove (INT x:[0, oo]. x^n * exp(-x)) = factorial(n) for n>=0, isInt(n)
induction on n
base:
    lhs:
        apply integral identity
        simplify
    rhs:
        simplify
    done
induct:
    lhs:
        integrate by parts with u=x^(n+1), v=-exp(-x)
        simplify
        apply induction hypothesis (all)
        rewrite (n+1)*factorial(n) to factorial(n+1)
    done
done
```

### 注意事项

- **归纳变量必须是整数**。
- **`apply induction hypothesis`** 只能在 `induct:` 分支的 CALCULATE 状态下使用。
- 如果归纳假设的左边在表达式中未出现，则无法应用。
