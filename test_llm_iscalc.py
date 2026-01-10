#!/usr/bin/env python
"""测试 LLM-Iscalc 集成系统"""

import asyncio
import sys

# 测试1: CommandExecutor
def test_command_executor():
    print("=" * 50)
    print("测试 CommandExecutor")
    print("=" * 50)
    
    from llm_iscalc.command_executor import CommandExecutor
    
    executor = CommandExecutor()
    
    # 测试初始化
    result = executor.initialize('INT x:[0,1]. x^2')
    print(f"初始化: {result.success}")
    print(f"表达式: {result.result}")
    
    # 测试执行命令
    result = executor.execute('apply integral identity')
    print(f"apply integral identity: {result.result}")
    
    result = executor.execute('simplify')
    print(f"simplify: {result.result}")
    
    print(f"完成: {executor.is_finished()}")
    print()


# 测试2: ErrorHandler
def test_error_handler():
    print("=" * 50)
    print("测试 ErrorHandler")
    print("=" * 50)
    
    from llm_iscalc.error_handler import ErrorHandler, ErrorCategory
    
    handler = ErrorHandler()
    
    # 测试错误分类
    try:
        raise ValueError("ParseException: invalid syntax")
    except Exception as e:
        error_info = handler.handle(e)
        print(f"错误类型: {error_info.category}")
        print(f"可恢复: {error_info.recoverable}")
        print(f"建议: {error_info.suggestion}")
    print()


# 测试3: SessionManager
def test_session_manager():
    print("=" * 50)
    print("测试 SessionManager")
    print("=" * 50)
    
    from llm_iscalc.session_manager import SessionManager
    
    manager = SessionManager()
    
    # 创建会话
    session = manager.create_session("INT x:[0,1]. x^2")
    print(f"会话ID: {session.id}")
    print(f"表达式: {session.expression}")
    print(f"状态: {session.status}")
    
    # 添加步骤
    manager.add_step(session.id, {"command": "test", "result": "ok"})
    print(f"历史步数: {len(session.history)}")
    print()


# 测试4: 完整流程（不调用LLM）
def test_full_flow_without_llm():
    print("=" * 50)
    print("测试完整流程（不调用LLM）")
    print("=" * 50)
    
    from llm_iscalc.command_executor import CommandExecutor
    
    executor = CommandExecutor()
    
    # 测试案例: INT x:[0,pi]. sin(x)
    print("表达式: INT x:[0,pi]. sin(x)")
    result = executor.initialize('INT x:[0,pi]. sin(x)')
    print(f"初始: {result.result}")
    
    commands = [
        'apply integral identity',
        'simplify'
    ]
    
    for cmd in commands:
        result = executor.execute(cmd)
        if result.success:
            print(f"{cmd}: {result.result}")
        else:
            print(f"{cmd}: 错误 - {result.error}")
    
    print(f"完成: {executor.is_finished()}")
    print()


# 测试5: 验证表达式
def test_validate_expression():
    print("=" * 50)
    print("测试表达式验证")
    print("=" * 50)
    
    from llm_iscalc.utils import validate_expression, expression_to_latex
    
    expressions = [
        "INT x:[0,1]. x^2",
        "sin(x) + cos(x)",
        "LIM {x -> 0}. sin(x)/x",
        "invalid expression @@#$"
    ]
    
    for expr in expressions:
        valid, result = validate_expression(expr)
        print(f"{expr}: {'有效' if valid else '无效'} - {result[:50]}...")
    print()


if __name__ == "__main__":
    test_command_executor()
    test_error_handler()
    test_session_manager()
    test_full_flow_without_llm()
    test_validate_expression()
    
    print("=" * 50)
    print("所有测试完成!")
    print("=" * 50)
