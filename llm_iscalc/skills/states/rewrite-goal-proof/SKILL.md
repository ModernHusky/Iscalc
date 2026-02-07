---
name: state-rewrite-goal-proof
description: 重写目标证明模式操作指南。
match_rules:
- from
- subgoal
---

# state-rewrite-goal-proof
> 重写目标证明模式操作指南。
## 使用时机
参考 description 描述。
## 指令
### 重写目标证明 (Proof by Rewrite Goal)

使用 `from <name>:` 从已证明的子目标开始证明。

### 语法

```
from goal1:
```

### 前提条件

- 已证明的子目标必须是**等式**（当前唯一支持的情况）
- `from` 后的名称必须是已定义的子目标名称

### 工作流程

1. 先使用 `subgoal <name>: <expr>` 声明并证明子目标
2. 然后使用 `from <name>:` 从该子目标开始重写
3. 继续应用变换直到达到目标
4. 使用 `done` 完成证明

### 示例

```
subgoal lemma1: INT x:[0,1]. f(x) = 1/2
  lhs:
    simplify
  done
from lemma1:
  rewrite 1/2 to 2/4
  done
```

### 完整示例 (来自 theories/interesting2.thy)



### 示例: 从已证子目标开始证明

```
prove (INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)) = pi / 8 * log(2)
subgoal 1: (INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)) = (INT x:[0,pi / 4]. log(tan(x) + 1))
lhs:
    substitute tan(u) for x
    rewrite sec(u) ^ 2 to tan(u) ^ 2 + 1
    simplify
done
subgoal 2: (INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)) = pi / 4 * log(2) - (INT x:[0,1]. log(x + 1) / (x ^ 2 + 1))
lhs:
    apply 1 on INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)
    ...
done
from 2:                               <-- 从 subgoal 2 开始
    solve equation for INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)
done
```
