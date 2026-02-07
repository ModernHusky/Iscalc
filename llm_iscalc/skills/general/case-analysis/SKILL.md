---
name: case-analysis
description: 基于条件的案例分析证明。
match_rules:
- (?i)case
- (?i)if
- '>\s*0|<\s*0|=\s*0'
---

# case-analysis
> 基于条件的案例分析证明。
## 使用时机
当证明目标依赖于某个条件的真假（如 `x > 0` 或 `x = 0`）时，使用案例分析拆分证明。
## 指令
### 快速开始



### 布尔条件 (两个分支)

```
case analysis on <condition>
    case true:
        lhs:
            ... # 在 condition = true 下证明
        done
    case false:
        lhs:
            ... # 在 condition = false 下证明
        done
done
```

### 数值条件 (三个分支)

```
case analysis on <expression>
    case positive:
        ... # expr > 0
        done
    case zero:
        ... # expr = 0
        done
    case negative:
        ... # expr < 0
        done
done
```

### 示例

证明 `x^4 + 2*x^2*cos(2*a) + 1 != 0`:

```
prove x^4 + 2*x^2*cos(2*a) + 1 != 0
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

### 注意事项

- 每个分支必须以 `done` 结束。
- 外层的 `done` 用于结束整个案例分析。
