---
name: series-operations
description: 级数展开、求值与收敛性操作。
match_rules:
- (?i)sum
- (?i)series
---

# series-operations
> 级数展开、求值与收敛性操作。
## 使用时机
处理涉及级数 (`SUM`) 的问题时使用。
## 指令
### 常用命令



### 级数展开

```
apply series expansion on <expr> index <n>
```
将表达式 `expr` 展开为关于索引变量 `n` 的级数。

### 级数求值

```
apply series evaluation
```
尝试自动计算级数的值。

### 交换积分与求和

**参见技能**: `exchange-operators`

在交换之前，通常需要先证明级数收敛：
```
subgoal 1: converges(SUM(n, ...))
arg:
    simplify
    ...
done
```

### 注意事项

- 级数索引变量 `n` 默认为整数类型，无需额外声明 `isInt(n)`。
- 上下界已隐含了 `n` 的取值范围，无需再加 `n >= 0` 等条件。
- 系统通过模式匹配自动判断收敛性，不支持 "apply alternating series test" 等自定义规则。
