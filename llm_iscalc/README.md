# LLM-Iscalc 使用说明

`llm_iscalc` 是 `iscalc` 项目的 LLM (大语言模型) 集成组件，允许通过自然语言交互和 LLM 辅助来执行数学推导和积分计算。

## 功能特点
- **交互式求解**: 使用 LLM 生成推导步骤。
- **自动化循环**: 自动执行 LLM 生成的 `iscalc` 命令。
- **Gradio 界面**: 友好的 Web 交互界面。
- **错误自愈**: 当命令执行失败时，LLM 会根据错误反馈调整策略。

## 快速开始

### 1. 安装依赖
确保你已经安装了 `iscalc` 的核心依赖，并添加 LLM 集成所需的库：
```bash
pip install gradio httpx aiohttp
```

### 2. 配置 API
编辑 `llm_iscalc/config.py` 文件，填写你的 DeepSeek API Key (或其他兼容 OpenAI 接口的服务)：
```python
@dataclass
class LLMConfig:
    api_key: str = "您的_API_KEY"
    api_base: str = "https://api.deepseek.com/v1"
    model: str = "deepseek-chat"
```

### 3. 运行服务
在项目根目录下运行启动脚本：
```bash
python run_llm_iscalc.py
```
默认会在 `http://127.0.0.1:7860` 启动 Gradio 界面。

## 参数说明
- `--port`: 指定运行端口（默认 7860）。
- `--share`: 创建一个可公开访问的 Gradio 链接。

示例：
```bash
python run_llm_iscalc.py --port 8080 --share
```

## 注意事项
- 确保网络能够访问指定的 `api_base`。
- LLM 生成的步骤会直接在本地 `iscalc` 环境中运行，请确保环境配置正确。
