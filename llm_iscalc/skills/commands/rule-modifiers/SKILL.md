---
name: rule-modifiers
description: 规则修饰符（at, on, for等）。
match_rules:
- .*
---

# rule-modifiers
> 规则修饰符（at, on, for等）。

## 指令

### 语法格式

**指定第 n 个匹配**：
```
<rule> (at n)
```

**应用于所有匹配**：
```
<rule> (all)
```

### 参数说明

- `<rule>`: 任何支持修饰符的规则
- `n`: 从 1 开始的匹配位置编号

## 使用时机

当表达式中有多个位置可以应用同一规则时，使用修饰符指定位置：
- 多个积分需要分别换元
- 多个相同子表达式需要部分替换
- 批量展开所有定义

## 适用目标类型

### `(at n)` 修饰符支持的规则

| 规则 | 说明 |
|-----|------|
| `rewrite A to B (at n)` | 应用于第 n 次出现的 A |
| `substitute u for expr (at n)` | 应用于第 n 个积分 |
| `substitute f(u) for x (at n)` | 逆换元，应用于第 n 个积分 |
| `integrate by parts ... (at n)` | 应用于第 n 个积分 |
| `split region at c (at n)` | 应用于第 n 个积分 |
| `expand definition for f (at n)` | 应用于第 n 次出现的 f |

### 不支持 `(at n)` 的规则

| 规则 | 原因 |
|-----|------|
| `apply integral identity` | 自动尝试所有积分 |
| `simplify` | 全局应用 |

## 工作流程

### `(at n)` 修饰符

1. **扫描表达式**：按深度优先顺序遍历表达式树
2. **计数匹配**：记录规则可应用的位置
3. **定位第 n 个**：找到第 n 个匹配位置
4. **执行规则**：仅在该位置应用规则

### `(all)` 修饰符

1. **扫描表达式**：找到所有规则可应用的位置
2. **批量应用**：在所有匹配位置应用规则
3. **返回结果**：返回完全替换后的表达式

## 注意事项

**⚠️ 格式要求**：
- 必须使用**括号**：`(at n)` 而不是 `at n`
- 括号与规则之间有**空格**

**错误示例对比**：

| ❌ 错误（导致解析失败） | ✅ 正确 |
|---------------------|--------|
| `substitute u for -x at 2` | `substitute u for -x (at 2)` |
| `rewrite A to B at 1` | `rewrite A to B (at 1)` |
| `substitute u for x(at 2)` | `substitute u for x (at 2)` |

## 示例

### 示例1: 指定位置替换

表达式: `sqrt(2) + sqrt(2) * x`

- `rewrite sqrt(2) to 2^(1/2)`: 只替换第一个
- `rewrite sqrt(2) to 2^(1/2) (at 2)`: 只替换第二个
- `rewrite sqrt(2) to 2^(1/2) (all)`: 全部替换

### 示例2: 多积分换元 (来自 theories/)

```
calculate INT x:[3, 4]. 1 / (x ^ 2 - 4)
    partial fraction decomposition
    simplify
    substitute u for 4 * x + 8
    substitute u for 4 * x - 8 (at 2)   # 对第二个积分换元
    apply integral identity
    simplify
done
```

### 示例3: 全部展开定义 (来自 theories/)

```
prove (INT x:[0,oo]. 1 / (x ^ 4 + 2 * x ^ 2 * cosh(2 * a) + 1)) = ...
lhs:
    expand definition for cosh (all)    # 展开所有 cosh
    ...
done
```
