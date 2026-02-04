---
name: series-operations
description: 级数展开、求值与收敛性操作。
match_rules:
- (?i)sum
- (?i)series
---

# series-operations
> 级数展开、求值与收敛性操作。

## 指令

### 语法格式

**级数展开**：
```
apply series expansion on <expr> index <n>
```

**级数求值**：
```
apply series evaluation
```

### 参数说明

- `<expr>`: 需要展开为级数的表达式
- `<n>`: 级数的索引变量名

## 使用时机

处理涉及级数 (`SUM`) 的问题时使用：
- 将函数展开为幂级数
- 计算已知级数的值
- 交换积分与求和的顺序

## 适用目标类型

| 目标类型 | 命令 | 适用性 |
|---------|------|--------|
| 可展开为级数的函数 | `apply series expansion` | ✅ |
| 已知形式的级数 | `apply series evaluation` | ✅ |
| `SUM(n, l, u, f(n))` | 各种操作 | ✅ |

## 工作流程

### 级数展开 `apply series expansion on <expr> index <n>`

1. **识别函数**：分析 `<expr>` 的类型（指数、对数、三角等）
2. **查找展开式**：在已知级数展开表中查找匹配
3. **生成级数**：用索引变量 `<n>` 生成 `SUM(n, ...)` 形式
4. **替换**：将原表达式替换为级数形式

### 级数求值 `apply series evaluation`

1. **识别级数**：找到表达式中的 `SUM` 结构
2. **模式匹配**：与已知级数求和公式匹配
3. **计算结果**：用闭式表达式替换级数

## 注意事项

- **索引变量**：级数索引变量 `n` 默认为整数类型，无需额外声明 `isInt(n)`。
- **范围约束**：上下界已隐含了 `n` 的取值范围，无需再加 `n >= 0` 等条件。
- **收敛性**：系统通过模式匹配自动判断收敛性，不支持自定义判据（如 "apply alternating series test"）。
- **交换顺序**：交换积分与求和前，可能需要先证明收敛性（参见 `exchange-operators` 技能）。

## 示例

### 示例1: 级数展开

表达式: `exp(x)`

```
apply series expansion on exp(x) index n
```

结果: `SUM(n, 0, oo, x^n / n!)`

### 示例2: 级数求值

表达式: `SUM(n, 0, oo, x^n)` （其中 |x| < 1）

```
apply series evaluation
```

结果: `1 / (1 - x)`

### 示例3: 收敛性证明与交换 (来自 theories/)

```
prove (INT x:[0,1]. SUM(n, 0, oo, x^n * f(n))) = ...
lhs:
    # 先证明收敛性（如需要）
    subgoal 1: converges(SUM(n, 0, oo, ...))
    arg:
        simplify
        ...
    done
    
    # 然后交换积分与求和
    exchange integral and sum
    apply integral identity
    simplify
done
```

### 示例4: 完整级数操作流程

```
calculate SUM(n, 1, oo, 1/n^2)
    # 这是已知级数，直接求值
    apply series evaluation
    # 结果：pi^2 / 6
done
```
