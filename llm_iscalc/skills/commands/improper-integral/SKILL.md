---
name: improper-integral
description: 处理广义积分（无穷限或奇点）。
match_rules:
- oo
- inf
---

# improper-integral
> 处理广义积分（无穷限或奇点）。

## 指令

### 语法格式

```
improper integral to limit creating <var>
```

### 参数说明

- `<var>`: 新创建的极限变量名（如 `t`），用于替换无穷边界

## 使用时机

当积分上下界含有无穷时：
- `INT x:[0,oo]. f(x)` — 上界为正无穷
- `INT x:[-oo,0]. f(x)` — 下界为负无穷
- `INT x:[-oo,oo]. f(x)` — 双无穷区间

## 适用目标类型

| 目标类型 | 适用性 |
|---------|--------|
| `INT x:[a,oo]. f(x)` | ✅ |
| `INT x:[-oo,b]. f(x)` | ✅ |
| `INT x:[-oo,oo]. f(x)` | ✅ |
| 有限区间定积分 | ❌ |

## 工作流程

1. **识别无穷边界**：找到积分上下界中的 `oo` 或 `-oo`
2. **引入极限变量**：创建新变量 `<var>`（如 `t`）
3. **替换边界**：将无穷边界替换为新变量
4. **添加极限**：在积分外层包裹 `LIM {<var> -> oo}`

**转换规则**：
- `INT x:[a,oo]. f(x)` → `LIM {t -> oo}. INT x:[a,t]. f(x)`
- `INT x:[-oo,b]. f(x)` → `LIM {t -> -oo}. INT x:[t,b]. f(x)`

## 注意事项

- **变量命名**：变量名 `<var>` 应避免与被积变量或表达式中其他变量冲突。
- **后续步骤**：转换后需要依次执行积分（`apply integral identity`）和极限计算（`simplify`）。
- **收敛性**：积分可能发散（结果为无穷），系统会给出相应提示。

## 示例

### 示例1: 指数衰减积分

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

最终结果: `1`

### 示例2: 高斯积分 (来自 theories/)

```
calculate INT x:[0,oo]. exp(-(a * x ^ 2)) for a > 0
    improper integral to limit creating t
    # 结果：LIM {t -> oo}. INT x:[0,t]. exp(-(a * x ^ 2))
    substitute u for sqrt(a) * x
    simplify
    ...
done
```

### 示例3: 双无穷区间

```
prove (INT x:[-oo,oo]. 1 / cosh(x)) = pi
lhs:
    split region at 0
    # 分成两个广义积分
    improper integral to limit creating t
    improper integral to limit creating s (at 2)
    ...
done
```
