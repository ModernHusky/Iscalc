"""Gradio 前端界面

提供流式输出的Web界面。
"""

import asyncio
import gradio as gr
from typing import Optional

from .config import create_config, LLMConfig, Config
from .llm_engine import LLMEngine
from .command_executor import CommandExecutor
from .solver_loop import SolverLoop, SolveEvent, EventType


# 全局状态 - 使用默认配置自动初始化
_config: Config = create_config()
_llm_engine = LLMEngine(_config.llm)
_executor = CommandExecutor(_config.iscalc.base_theory)
_solver: SolverLoop = SolverLoop(_llm_engine, _executor, _config.solver)

# 终止标志
_stop_flag = False


def init_solver(api_key: str = None, api_base: str = None, model: str = None):
    """初始化求解器（可选，已有默认配置）"""
    global _solver, _config, _llm_engine, _executor
    
    if api_key or api_base or model:
        _config = create_config(
            api_key=api_key or _config.llm.api_key,
            api_base=api_base or _config.llm.api_base,
            model=model or _config.llm.model
        )
        _llm_engine = LLMEngine(_config.llm)
        _executor = CommandExecutor(_config.iscalc.base_theory)
        _solver = SolverLoop(_llm_engine, _executor, _config.solver)
    
    return f"求解器已初始化！使用模型: {_config.llm.model}"


def stop_solving():
    """终止求解"""
    global _stop_flag
    _stop_flag = True
    return "已发送终止信号"


async def solve_expression_async(expression: str, conditions: str = ""):
    """异步求解表达式"""
    global _solver, _config, _executor, _stop_flag
    
    _stop_flag = False  # 重置终止标志
    
    if not expression.strip():
        yield "", "", "", ""
        return
    
    # 解析条件
    cond_list = None
    if conditions.strip():
        cond_list = [c.strip() for c in conditions.split(",") if c.strip()]
    
    # 状态变量
    commands_text = ""
    results_text = ""
    errors_text = ""
    thinking_text = ""
    last_command = ""  # 记录上一个命令
    results_lines = []  # 用于支持覆盖操作
    current_step_thinking = ""  # 当前步骤的思考内容
    
    async def callback(event: SolveEvent):
        nonlocal commands_text, results_text, errors_text, thinking_text
        nonlocal last_command, results_lines, current_step_thinking
        
        if event.type == EventType.THINKING:
            # 保留每一步的思考过程
            if event.metadata.get("streaming"):
                # 流式输出时，更新当前步骤的思考内容
                current_step_thinking = event.content
                # 重建完整的思考文本
                if thinking_text:
                    # 找到最后一个分隔线，替换之后的内容
                    sep = "=" * 50
                    if sep in thinking_text:
                        parts = thinking_text.rsplit(sep, 1)
                        thinking_text = parts[0] + sep + "\n" + current_step_thinking
                    else:
                        thinking_text = current_step_thinking
                else:
                    thinking_text = current_step_thinking
            else:
                current_step_thinking = event.content
                thinking_text = event.content
        elif event.type == EventType.COMMAND:
            last_command = event.content
            commands_text += f"[步骤 {event.step}] {event.content}\n"
            if event.metadata.get("explanation"):
                commands_text += f"  说明: {event.metadata['explanation']}\n"
            # 命令生成后，保存当前思考并添加分隔符
            thinking_text += f"\n{'='*50}\n"
        elif event.type == EventType.RESULT:
            # 格式: = 表达式 (命令)
            if last_command:
                line = f"= {event.content} ({last_command})"
            else:
                line = f"  {event.content}"
            results_lines.append(line)
            results_text = "\n".join(results_lines) + "\n"
        elif event.type == EventType.ERROR:
            errors_text += f"[步骤 {event.step}] {event.content}\n"
        elif event.type == EventType.COMPLETE:
            results_text = "\n".join(results_lines) + f"\n\n✓ 求解完成\n"
        elif event.type == EventType.LOOP_RETRY:
            # 循环重试：移除最后一行结果，准备覆盖
            if results_lines:
                results_lines.pop()
                results_text = "\n".join(results_lines) + "\n"
            errors_text += f"⚠ {event.content}\n"
        elif event.type == EventType.LOOP_DETECTED:
            errors_text += f"⚠ {event.content}\n"
        elif event.type == EventType.MAX_ITERATIONS:
            errors_text += f"⚠ {event.content}\n"
    
    # 重新创建executor以重置状态
    _executor = CommandExecutor(_config.iscalc.base_theory)
    _solver.executor = _executor
    _solver.error_handler.clear_error_log()
    
    # 初始化时显示初始表达式
    init_result = _executor.initialize(expression, cond_list)
    if init_result.success:
        results_lines = [f"(calculate)", f"  {init_result.result}"]
        results_text = "\n".join(results_lines) + "\n"
    else:
        errors_text = f"初始化失败: {init_result.error}\n"
        yield commands_text, results_text, errors_text, thinking_text
        return
    
    yield commands_text, results_text, errors_text, thinking_text
    
    # 重新初始化executor用于solver
    _executor = CommandExecutor(_config.iscalc.base_theory)
    _solver.executor = _executor
    
    # 创建任务
    solve_task = asyncio.create_task(
        _solver.solve(expression, callback, cond_list)
    )
    
    # 定期更新UI，检查终止标志
    while not solve_task.done():
        if _stop_flag:
            solve_task.cancel()
            results_text += f"\n⚠ 已终止求解\n"
            break
        yield commands_text, results_text, errors_text, thinking_text
        await asyncio.sleep(0.1)
    
    # 获取最终结果
    if not _stop_flag:
        try:
            result = await solve_task
            if not result.success and result.error:
                errors_text += f"\n求解失败: {result.error}\n"
        except asyncio.CancelledError:
            pass
    
    yield commands_text, results_text, errors_text, thinking_text



def solve_expression(expression: str, conditions: str = ""):
    """同步包装器"""
    async def run():
        async for result in solve_expression_async(expression, conditions):
            yield result
    
    loop = asyncio.new_event_loop()
    asyncio.set_event_loop(loop)
    
    gen = run()
    try:
        while True:
            try:
                result = loop.run_until_complete(gen.__anext__())
                yield result
            except StopAsyncIteration:
                break
    finally:
        loop.close()


def create_interface() -> gr.Blocks:
    """创建Gradio界面"""
    
    with gr.Blocks(
        title="LLM-Iscalc 数学求解器",
        theme=gr.themes.Soft(),
        css="""
        .output-box { font-family: monospace; white-space: pre-wrap; }
        """
    ) as app:
        gr.Markdown("# 🧮 LLM-Iscalc 数学求解器")
        gr.Markdown("输入数学表达式，AI将自动分析并使用Iscalc工具进行化简求解。")
        gr.Markdown(f"*当前使用模型: {_config.llm.model}*")
        
        # 输入区域
        with gr.Row():
            with gr.Column(scale=3):
                expression_input = gr.Textbox(
                    label="📝 数学表达式",
                    placeholder="例如: INT x:[0,1]. x^2",
                    lines=2
                )
            with gr.Column(scale=2):
                conditions_input = gr.Textbox(
                    label="📋 附加条件 (可选，逗号分隔)",
                    placeholder="例如: x > 0, n >= 0",
                    lines=2
                )
        
        with gr.Row():
            solve_btn = gr.Button("🚀 开始求解", variant="primary", scale=2)
            stop_btn = gr.Button("⏹️ 终止", variant="stop", scale=1)
            clear_btn = gr.Button("🗑️ 清空", scale=1)
        
        # 思考过程
        with gr.Accordion("🤔 AI思考过程", open=True):
            thinking_output = gr.Textbox(
                label="",
                lines=8,
                interactive=False,
                elem_classes=["output-box"]
            )
        
        # 结果区域 - 单独一行
        gr.Markdown("### ✅ 结果")
        result_output = gr.Textbox(
            label="",
            lines=12,
            interactive=False,
            elem_classes=["output-box"]
        )
        
        # 命令和错误区域 - 两列布局
        with gr.Row():
            with gr.Column():
                gr.Markdown("### 📋 命令")
                command_output = gr.Textbox(
                    label="",
                    lines=12,
                    interactive=False,
                    elem_classes=["output-box"]
                )
            
            with gr.Column():
                gr.Markdown("### ❌ 错误")
                error_output = gr.Textbox(
                    label="",
                    lines=12,
                    interactive=False,
                    elem_classes=["output-box"]
                )
        
        # 示例
        gr.Markdown("### 📚 示例表达式")
        gr.Examples(
            examples=[
                ["INT x:[0,1]. x^2", ""],
                ["INT x:[0,pi]. sin(x)", ""],
                ["INT x:[0,1]. x * exp(x)", ""],
                ["INT x:[0,1]. 1/(1+x^2)", ""],
                ["LIM {x -> 0}. sin(x)/x", ""],
                ["D x. x^3 + 2*x", ""],
            ],
            inputs=[expression_input, conditions_input],
            label=""
        )
        
        # 事件绑定
        solve_btn.click(
            fn=solve_expression,
            inputs=[expression_input, conditions_input],
            outputs=[command_output, result_output, error_output, thinking_output]
        )
        
        stop_btn.click(
            fn=stop_solving,
            outputs=[]
        )
        
        def clear_outputs():
            return "", "", "", ""
        
        clear_btn.click(
            fn=clear_outputs,
            outputs=[command_output, result_output, error_output, thinking_output]
        )
    
    return app


def main():
    """主函数"""
    app = create_interface()
    app.queue()
    app.launch(share=False, server_name="127.0.0.1", server_port=7860)


if __name__ == "__main__":
    main()
