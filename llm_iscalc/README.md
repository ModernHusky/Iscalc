# LLM-Iscalc - 基于大语言模型的数学表达式求解器

LLM-Iscalc 是一个结合了大语言模型（LLM）和符号计算引擎（Iscalc）的智能数学求解系统。它能够自动分析数学表达式，生成求解策略，并逐步执行符号计算来化简和求解复杂的数学问题。

## ✨ 特性

- 🤖 **智能求解**: 使用 LLM 自动分析表达式并生成求解策略
- 🔄 **迭代优化**: 根据执行结果动态调整求解方案
- 📊 **实时反馈**: 流式显示思考过程和计算步骤
- 🎯 **多种积分技巧**: 支持换元、分部积分、部分分式分解等
- 🌐 **Web 界面**: 简洁友好的交互界面
- 🔍 **循环检测**: 自动检测并避免求解循环
- 📝 **LaTeX 输出**: 支持数学公式的 LaTeX 格式显示

## 🚀 快速开始

### 环境要求

- Python 3.10 或更高版本
- 已安装 iscalc 项目的依赖包

### 安装

1. 确保已安装主项目依赖：
```bash
pip install -r requirements.txt
```

2. 配置 API 密钥：

编辑 `llm_iscalc/config.py` 文件，设置你的 LLM API 配置：

```python
@dataclass
class LLMConfig:
    api_key: str = "your-api-key-here"
    api_base: str = "https://api.deepseek.com/v1"
    model: str = "deepseek-chat"
```

### 运行

启动 Flask 服务器：

```bash
python run_llm_iscalc_flask.py
```

或指定端口和主机：

```bash
python run_llm_iscalc_flask.py --port 8080 --host 0.0.0.0
```

然后在浏览器中访问：`http://127.0.0.1:7860`

## 📖 使用方法

### 基本用法

1. 在输入框中输入数学表达式，例如：
   - `INT x:[0,1]. x^2` - 定积分
   - `INT x. sin(x)` - 不定积分
   - `LIM {x -> 0}. sin(x)/x` - 极限
   - `SUM(n, 1, oo, 1/n^2)` - 级数求和

2. （可选）添加条件约束，例如：`x > 0, n > 1`

3. 点击"开始求解"或按 `Ctrl/Cmd + Enter`

4. 观察 LLM 的思考过程和执行步骤

5. 按 `Esc` 键可随时终止求解

### 表达式语法

#### 基本元素
- **变量**: `x`, `y`, `z`, `t`, `n`, `m` 等
- **常数**: `pi` (π), `e`, `i` (虚数单位), `G` (欧拉常数)
- **无穷**: `oo`, `inf` (正无穷), `-oo`, `-inf` (负无穷)

#### 运算符
- 四则运算: `+`, `-`, `*`, `/`
- 幂运算: `^` (例如 `x^2`)
- 绝对值: `|x|`

#### 函数
- 三角函数: `sin(x)`, `cos(x)`, `tan(x)`, `cot(x)`
- 反三角函数: `arcsin(x)`, `arccos(x)`, `arctan(x)`
- 指数对数: `exp(x)`, `log(x)`, `sqrt(x)`
- 特殊函数: `gamma(x)`, `factorial(n)`, `binom(n,k)`

#### 积分与极限
- 定积分: `INT x:[a,b]. f(x)`
- 不定积分: `INT x. f(x)`
- 极限: `LIM {x -> a}. f(x)`
- 左极限: `LIM {x -> a-}. f(x)`
- 右极限: `LIM {x -> a+}. f(x)`

### 示例

#### 示例 1: 简单定积分
```
表达式: INT x:[0,1]. x^2
结果: 1/3
```

#### 示例 2: 换元积分
```
表达式: INT x:[0,1]. 2*x*(x^2+1)^3
步骤:
1. substitute u for x^2 + 1
2. apply integral identity
3. simplify
4. replace substitution
```

#### 示例 3: 分部积分
```
表达式: INT x:[0,1]. x*exp(x)
步骤:
1. integrate by parts with u = x, v = exp(x)
2. simplify
3. apply integral identity
```

#### 示例 4: 极限计算
```
表达式: LIM {x -> oo}. x/(x+1)
步骤:
1. rewrite x/(x+1) to 1 - 1/(x+1)
2. simplify
```

## 🏗️ 架构设计

### 核心模块

```
llm_iscalc/
├── app_flask.py          # Flask Web 服务器
├── llm_engine.py         # LLM API 交互引擎
├── solver_loop.py        # 求解循环控制器
├── command_executor.py   # Iscalc 命令执行器
├── error_handler.py      # 错误处理器
├── session_manager.py    # 会话管理器
├── prompts.py           # LLM 提示词模板
├── config.py            # 配置管理
├── utils.py             # 工具函数
└── static/              # 前端资源
    ├── index.html       # Web 界面
    └── app.js           # 前端逻辑
```

### 工作流程

1. **用户输入** → 表达式和条件
2. **初始化** → CommandExecutor 解析并初始化计算状态
3. **LLM 分析** → LLMEngine 生成求解策略和命令
4. **命令执行** → CommandExecutor 执行 Iscalc 命令
5. **结果反馈** → 将结果返回给 LLM 进行下一步分析
6. **迭代求解** → 重复步骤 3-5 直到完成或达到最大迭代次数
7. **输出结果** → 显示最终化简结果

### 关键组件

#### LLMEngine
负责与 LLM API 交互，支持流式响应：
- 构建包含历史记录的提示词
- 处理错误反馈
- 解析 LLM 返回的 JSON 格式命令

#### SolverLoop
协调整个求解过程：
- 管理迭代循环
- 检测表达式循环
- 处理错误重试
- 触发事件回调

#### CommandExecutor
执行 Iscalc 命令：
- 维护计算状态
- 解析和执行命令
- 转换 LaTeX 格式
- 记录执行历史

## ⚙️ 配置选项

### LLM 配置 (LLMConfig)

```python
api_key: str           # API 密钥
api_base: str          # API 基础 URL
model: str             # 模型名称
temperature: float     # 温度参数 (0.0-1.0)
max_tokens: int        # 最大生成 token 数
timeout: float         # 请求超时时间（秒）
max_retries: int       # 最大重试次数
retry_delay: float     # 重试延迟（秒）
```

### 求解器配置 (SolverConfig)

```python
max_iterations: int              # 最大迭代次数
max_consecutive_errors: int      # 最大连续错误次数
timeout_seconds: float           # 总超时时间（秒）
loop_detection_window: int       # 循环检测窗口大小
```

### Iscalc 配置 (IscalcConfig)

```python
base_theory: str                 # 基础理论文件
default_conditions: list         # 默认条件列表
```

## 🔧 API 接口

### GET /api/config
获取当前配置信息

**响应**:
```json
{
  "model": "deepseek-chat",
  "api_base": "https://api.deepseek.com/v1"
}
```

### GET /api/solve
流式求解表达式（SSE）

**参数**:
- `expression`: 数学表达式（必需）
- `conditions`: 条件约束，逗号分隔（可选）

**响应**: Server-Sent Events 流
```javascript
data: {
  "commands": "执行的命令历史",
  "results": "计算结果",
  "errors": "错误信息",
  "thinking": {
    "step": 1,
    "thinking": "思考过程",
    "command": "生成的命令",
    "explanation": "命令说明",
    "is_final": false
  }
}
```

### POST /api/stop
终止当前求解

**响应**:
```json
{
  "status": "ok",
  "message": "已发送终止信号"
}
```

## 🎯 支持的 Iscalc 命令

### 变量替换
- `substitute u for <expr>` - 换元
- `substitute <expr> for u` - 逆替换
- `replace substitution` - 替换回代

### 积分技巧
- `apply integral identity` - 应用积分公式
- `integrate by parts with u = <expr>, v = <expr>` - 分部积分
- `split region at <point>` - 区域分割
- `improper integral to limit creating <var>` - 广义积分转极限

### 表达式变换
- `rewrite <old> to <new>` - 表达式重写
- `simplify` - 自动简化
- `expand polynomial` - 多项式展开
- `partial fraction decomposition` - 部分分式分解

### 高级技巧
- `l'Hopital's rule` - 洛必达法则
- `exchange integral and sum` - 交换积分与求和
- `exchange derivative and integral` - 交换导数与积分
- `apply series expansion on <expr> index <var>` - 级数展开
- `solve integral <expr>` - 求解积分方程

## 🐛 错误处理

系统具有完善的错误处理机制：

1. **语法错误**: 自动识别并提示正确的语法格式
2. **规则错误**: 分析规则应用失败的原因并尝试其他方法
3. **循环检测**: 检测到表达式循环时自动重试不同策略
4. **连续错误**: 超过阈值时建议检查输入表达式
5. **超时保护**: 防止无限循环，达到最大迭代次数自动停止

## 📊 性能优化

- **流式响应**: 实时显示 LLM 思考过程，提升用户体验
- **异步处理**: 使用 asyncio 实现高效的异步 I/O
- **状态缓存**: 缓存计算状态，避免重复计算
- **智能重试**: 指数退避策略处理 API 调用失败

## 🤝 贡献指南

欢迎贡献代码、报告问题或提出建议！

### 开发环境设置

1. Fork 并克隆仓库
2. 创建虚拟环境并安装依赖
3. 运行测试：`pytest tests/`
4. 提交 Pull Request

### 代码规范

- 遵循 PEP 8 代码风格
- 添加必要的类型注解
- 编写清晰的文档字符串
- 为新功能添加测试用例

## 📝 许可证

本项目遵循主项目 iscalc 的许可证。

## 🙏 致谢

- **Iscalc**: 强大的符号计算引擎
- **DeepSeek**: 提供高质量的 LLM API 服务
- **Flask**: 轻量级 Web 框架

## 📮 联系方式

如有问题或建议，请通过以下方式联系：

- 提交 Issue
- 发起 Discussion
- 提交 Pull Request

---

**注意**: 本项目仍在积极开发中，API 和功能可能会发生变化。
