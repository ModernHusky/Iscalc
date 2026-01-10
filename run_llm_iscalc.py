#!/usr/bin/env python
"""LLM-Iscalc 启动脚本

使用方法:
    python run_llm_iscalc.py [--port PORT] [--share]

示例:
    python run_llm_iscalc.py
    python run_llm_iscalc.py --port 8080
    python run_llm_iscalc.py --share
"""

import argparse
import os
import sys

# 确保可以导入llm_iscalc模块
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def main():
    parser = argparse.ArgumentParser(description="LLM-Iscalc 数学求解器")
    parser.add_argument("--port", type=int, default=7860, help="服务端口")
    parser.add_argument("--share", action="store_true", help="创建公共链接")
    
    args = parser.parse_args()
    
    from llm_iscalc.app_gradio import create_interface
    from llm_iscalc.config import default_config
    
    print(f"LLM-Iscalc 数学求解器")
    print(f"使用模型: {default_config.llm.model}")
    print(f"API地址: {default_config.llm.api_base}")
    print("-" * 40)
    
    # 启动界面
    app = create_interface()
    app.queue()
    app.launch(
        share=args.share,
        server_name="127.0.0.1",
        server_port=args.port
    )


if __name__ == "__main__":
    main()
