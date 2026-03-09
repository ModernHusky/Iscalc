---
name: state-case
description: 案例分析状态操作指南。
match_rules:
- case
---

# state-case
> 案例分析状态操作指南。
## 使用时机
参考 description 描述。
## 指令
### 案例分析状态 (Case State)

在案例分析状态下，可以执行以下操作：

### 布尔条件分支 (两路)

当使用 `case analysis on <condition>` 时：

- `case true:` - 进入条件为真的分支
- `case false:` - 进入条件为假的分支

### 数值表达式分支 (三路)

当使用 `case analysis on <expression>` 时：

- `case positive:` - 进入表达式为正的分支 (expr > 0)
- `case zero:` - 进入表达式为零的分支 (expr = 0)
- `case negative:` - 进入表达式为负的分支 (expr < 0)

### 注意事项

- 每个分支进入后会转为 Prove 状态
- 所有分支都必须完成证明
- 完成所有分支后使用 `done` 退出

### 完整示例 (来自 theories/interesting2.thy)



### 示例: 布尔条件分支 (x != 0)

```
subgoal 1: x^4 + 2*x^2*cos(2*a) + 1 != 0
case analysis on x != 0
    case true:
    lhs:
        rewrite to (x^2 - 1)^2 + 2*x^2*(1 + cos(2*a))
        rewrite cos(2*a) to 2*cos(a)^2 - 1
        simplify
    done
    case false:
    lhs:
        simplify                    
    done
done
```

### 示例: 数值表达式分支 (cos(a))

```
case analysis on cos(a)
    case negative:
    lhs:
        apply 8 on (INT x:[0,oo]. 1/(x^4+...))
    rhs:
        simplify
    done
    case positive:
    lhs:
        apply 7 on (INT x:[0,oo]. 1/(x^4+...))
    rhs:
        simplify
    done
done
```

