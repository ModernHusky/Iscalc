---
name: state-func-def
description: 函数定义 - 创建新的函数定义
match_rules:
  - "define"
applicable_types:
  - general
---

# 函数定义 (Function Definition)

在全局上下文或证明中引入新的函数定义。

## 语法

```
define f(x) = <expr>
define I = INT x:[0,1]. f(x)
```

## 规则

- 定义必须是**等式**形式
- 左边必须是：
  - 一个变量，或
  - 一个函数符号应用于不同变量的列表

## 示例

```
define f(x) = x^2 + 1
define I = INT x:[0,oo]. exp(-x^2)
define g(a, b) = a * b + 1
```

## 注意事项

- 定义后可以使用 `expand definition` 展开
- 也可以使用 `fold definition` 折叠表达式
- 条件可用 `for` 指定: `define f(x) = ... for x > 0`

# 完整示例 (来自 theories/)

## 示例: 定义积分函数
```
define I(t) = INT x:[0,1]. x^t for t > -1

prove D t. I(t) = (INT x:[0,1]. x^t * log(x)) for t > -1
lhs:
    expand definition for I
    ...
done
```

## 示例: 在证明中使用定义
```
prove (INT x:[0,oo]. 1 / (x ^ 4 + 2 * x ^ 2 * cosh(2 * a) + 1)) = pi / (4 * cosh(a))
lhs:
    expand definition for cosh (all)
    ...
    fold definition for cosh (all)
done
```
