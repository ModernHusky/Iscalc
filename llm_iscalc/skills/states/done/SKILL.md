---
name: state-done
description: 完成条件 - 何时使用 done 命令
match_rules:
  - ".*"
applicable_types:
  - general
---

# 完成条件 (Using Done)

`done` 命令用于表示当前计算或证明已完成。

## 完成条件

### 计算 (Calculate)
- 当前表达式已达到**闭合形式** (closed form)
- 即：表达式不包含积分、极限等未求值的操作

### 等式证明 (Prove equality)
- 等式两边经过化简后可被证明相等

### 不等式/非等式证明
- 性质可通过标准自动化方法证明

### 归纳法/案例分析
- 所有分支都已完成证明

## 注意事项

- 过早使用 `done` 会导致错误
- 如果系统提示 "not finished"，需要继续变换
- 检查是否有未处理的积分、极限或级数

## 常见的"未完成"原因

### 1. 等式无法证明 (Equality Cannot Be Shown)
- 等式两边在计算后仍不相等
- **解决方案**：继续使用 `rewrite`, `simplify` 等规则化简

### 2. 缺少良构性条件 (Missing Wellformedness Conditions)
- 表达式依赖的某些条件未被证明
- **错误信息会显示缺少的条件**
- **解决方案**：使用 `subgoal` 命令证明这些条件

示例错误信息：
```
CheckFinishedException: proof not finished
Missing wellformedness conditions:
  - x > 0
  - a != 0
```

需要添加子目标：
```
subgoal c1: x > 0
...
done
subgoal c2: a != 0
...
done
```

### 3. 计算未达到闭合形式 (Not in Closed Form)
- 结果中仍包含未求值的积分 `INT`、极限 `LIM` 或级数 `SUM`
- **解决方案**：继续应用规则直到所有符号都被求值

检查清单：
- ✓ 所有积分是否已应用 `apply integral identity`？
- ✓ 所有极限是否已 `simplify` 求值？
- ✓ 所有级数是否已应用 `apply series evaluation`？

# 完整示例 (来自 theories/)

## 示例: 计算完成
```
calculate INT x:[0, 1]. x*exp(x)
    integrate by parts with u = x, v = exp(x)
    apply integral identity
    simplify
done     <-- 表达式已达到闭合形式
```

## 示例: 证明完成
```
prove (INT x:[1,oo]. 1 / ((x+a)*sqrt(x-1))) = pi / sqrt(a+1) for a > -1
lhs:
    substitute t for sqrt(x - 1)
    simplify
    substitute y for t / sqrt(a + 1)
    apply integral identity
    simplify
done     <-- 左右两边已证明相等
```
