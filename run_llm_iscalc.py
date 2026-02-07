#!/usr/bin/env python
"""LLM-Iscalc Flask 版本启动脚本

使用方法:
    python run_llm_iscalc_flask.py [--port PORT] [--host HOST] [--debug]

示例:
    python run_llm_iscalc_flask.py
    python run_llm_iscalc_flask.py --port 8080
    python run_llm_iscalc_flask.py --host 0.0.0.0 --debug
"""

import argparse
import os
import sys

# 确保可以导入llm_iscalc模块
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def main():
    parser = argparse.ArgumentParser(description="LLM-Iscalc 数学求解器 (Flask 版本)")
    parser.add_argument("--port", type=int, default=7860, help="服务端口 (默认: 7860)")
    parser.add_argument("--host", type=str, default="127.0.0.1", help="服务地址 (默认: 127.0.0.1)")
    parser.add_argument("--debug", action="store_true", help="启用调试模式")
    
    args = parser.parse_args()
    
    import sys
    import os

    from llm_iscalc.app_flask import app
    from llm_iscalc.app_flask import main as flask_main
    from llm_iscalc.config import default_config
    
    print("=" * 60)
    print("🧮 LLM-Iscalc 数学求解器 (Flask 版本)")
    print("=" * 60)
    print(f"🤖 使用模型: {default_config.llm.model}")
    print(f"🌐 API地址: {default_config.llm.api_base}")
    print(f"📍 访问地址: http://{args.host}:{args.port}")
    print("=" * 60)
    print("\n提示:")
    print("  - 按 Ctrl+C 停止服务器")
    print("  - 使用 Ctrl/Cmd+Enter 快速提交表达式")
    print("  - 使用 Esc 键终止正在进行的求解")
    print("\n")
    
    # 启动 Flask 服务器
    flask_main(host=args.host, port=args.port, debug=args.debug)


if __name__ == "__main__":
    main()
