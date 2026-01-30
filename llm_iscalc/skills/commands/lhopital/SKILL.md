---
name: lhopital
description: 洛必达法则求不定式极限。
match_rules:
- (?i)lim
---

# lhopital
> 洛必达法则求不定式极限。
## 使用时机
当极限是不定型时：
- 0/0 型
- ∞/∞ 型
## 指令
### 快速开始

```
l'Hopital's rule
```

对分子分母同时求导。

### 分步指南

1. 确认极限是 0/0 或 ∞/∞ 型
2. 应用 `l'Hopital's rule`
3. 如果仍是不定型，可再次应用
4. 最后 `simplify` 得结果

### 示例



### 示例1: 0/0 型

表达式: `LIM {x -> 0}. sin(x)/x`

```
l'Hopital's rule
```

变换为: `LIM {x -> 0}. cos(x)/1 = 1`

### 示例2: ∞/∞ 型

表达式: `LIM {x -> oo}. log(x)/x`

```
l'Hopital's rule
```

变换为: `LIM {x -> oo}. (1/x)/1 = 0`

### 注意事项

- 必须是分数形式的极限
- 直接代入得到不定型时才使用
- 某些情况下用 `rewrite` 更简单
