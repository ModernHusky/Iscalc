"""Flask 前端界面 - 简化版本"""

import asyncio
import json
from pathlib import Path
from flask import Flask, request, Response, jsonify, send_from_directory

from .config import create_config
from .llm_engine import LLMEngine
from .command_executor import CommandExecutor
from .solver_loop import SolverLoop, SolveEvent, EventType

# 获取当前文件所在目录
_current_dir = Path(__file__).parent

# 创建 Flask 应用
app = Flask(__name__, 
            static_folder=str(_current_dir / 'static'),
            static_url_path='/static')

# 全局状态
_config = create_config()
_llm_engine = LLMEngine(_config.llm)
_stop_flag = False


@app.route('/')
def index():
    """主页"""
    return send_from_directory(app.static_folder, 'index.html')


@app.route('/test')
def test():
    """测试页面"""
    return send_from_directory(app.static_folder, 'test.html')


@app.route('/api/config')
def get_config():
    """获取配置信息"""
    return jsonify({
        'model': _config.llm.model,
        'api_base': _config.llm.api_base
    })


@app.route('/api/stop', methods=['POST'])
def stop_solving():
    """终止求解"""
    global _stop_flag
    _stop_flag = True
    return jsonify({'status': 'ok', 'message': '已发送终止信号'})


@app.route('/api/solve')
def solve_expression():
    """求解表达式（SSE 流式输出）"""
    expression = request.args.get('expression', '').strip()
    conditions = request.args.get('conditions', '').strip()
    instruction = request.args.get('instruction', '').strip()
    
    if not expression:
        return jsonify({'error': '表达式不能为空'}), 400
    
    def generate():
        """生成 SSE 事件流"""
        global _stop_flag
        
        _stop_flag = False
        
        # 解析条件
        cond_list = None
        if conditions:
            cond_list = [c.strip() for c in conditions.split(",") if c.strip()]
        
        # 状态变量
        state = {
            'commands': "",
            'results': "",
            'errors': "",
            'thinking': {
                'thinking': '',
                'command': '',
                'explanation': '',
                'is_final': False,
                'step': -1  # 初始值为 -1，与前端一致
            }
        }
        
        def send_update():
            """发送更新"""
            data = {
                'commands': state['commands'],
                'results': state['results'],
                'errors': state['errors'],
                'thinking': state['thinking']
            }
            return f"data: {json.dumps(data, ensure_ascii=False)}\n\n"
        
        # 创建事件循环
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)
        
        # 创建消息队列
        queue = asyncio.Queue()
        
        try:
            # 创建 executor 和 solver
            executor = CommandExecutor(_config.iscalc.base_theory)
            solver = SolverLoop(_llm_engine, executor, _config.solver)
            
            # 初始化
            init_result = executor.initialize(expression, cond_list)
            if init_result.success:
                # 获取实际状态名称，而不是硬编码
                current_state = executor.get_current_state_name().lower()
                state['results'] = f"({current_state})\n  {init_result.result}\n"
            else:
                state['errors'] = f"初始化失败: {init_result.error}\n"
                loop.run_until_complete(queue.put(send_update()))
                yield send_update()
                yield "event: complete\ndata: {}\n\n"
                return
            
            # 初始更新
            loop.run_until_complete(queue.put(send_update()))
            
            # 回调函数
            # 获取实际状态名称
            current_state = executor.get_current_state_name().lower()
            results_lines = [f"({current_state})", f"  {init_result.result}"]
            last_command = ""
            
            # 状态切换命令列表（这些命令应该单独成行显示）
            state_switch_commands = ['lhs:', 'rhs:', 'arg:', 'base:', 'induct:', 'done',
                                     'case true:', 'case false:', 'case positive:', 
                                     'case negative:', 'case zero:']
            
            async def callback(event: SolveEvent):
                nonlocal last_command
                
                if event.type == EventType.THINKING:
                    # 如果步骤号变化，重置思考状态
                    if state['thinking']['step'] != event.step:
                         state['thinking'] = {
                            'step': event.step,
                            'thinking': '',
                            'command': '',
                            'explanation': '',
                            'is_final': False,
                        }
                    else:
                        state['thinking']['step'] = event.step
                    
                    # 尝试解析JSON并提取字段
                    try:
                        import re
                        import json
                        content = event.content
                        
                        # 预处理：去除可能存在的 Markdown 代码块标记
                        clean_content = content
                        if "```" in content:
                           clean_content = re.sub(r'```(?:json)?', '', content).strip()
                        
                        # 优先尝试正则解析（适合流式传输，能容忍 JSON 未闭合）
                        thinking_match = re.search(r'"thinking"\s*:\s*"((?:[^"\\]|\\.)*)(?:"|$)', clean_content, re.DOTALL)
                        if thinking_match:
                            val = thinking_match.group(1).replace('\\n', '\n').replace('\\"', '"')
                            state['thinking']['thinking'] = val
                        
                        command_match = re.search(r'"command"\s*:\s*"([^"]*)(?:"|$)', clean_content)
                        if command_match:
                            state['thinking']['command'] = command_match.group(1)
                            
                        explanation_match = re.search(r'"explanation"\s*:\s*"((?:[^"\\]|\\.)*)(?:"|$)', clean_content, re.DOTALL)
                        if explanation_match:
                            val = explanation_match.group(1).replace('\\n', '\n').replace('\\"', '"')
                            state['thinking']['explanation'] = val
                            
                        is_final_match = re.search(r'"is_final"\s*:\s*(true|false)', clean_content)
                        if is_final_match:
                            state['thinking']['is_final'] = is_final_match.group(1) == 'true'

                        if not state['thinking']['thinking'] and '{' not in clean_content:
                             if not clean_content.startswith('[DEBUG'):
                                state['thinking']['thinking'] = clean_content

                    except Exception:
                        pass
                    
                elif event.type == EventType.COMMAND:
                    is_intermediate = event.metadata.get("intermediate", False)
                    
                    if is_intermediate:
                        state['thinking'] = {
                            'step': event.step,
                            'thinking': '',
                            'command': event.content,
                            'explanation': '',
                            'is_final': False,
                        }
                    else:
                        last_command = event.content
                    
                    state['commands'] += f"[步骤 {event.step}] {event.content}\n"
                    if event.metadata.get("explanation"):
                        state['commands'] += f"  说明: {event.metadata['explanation']}\n"
                    
                elif event.type == EventType.RESULT:
                    # 处理中间步骤：为每个自动应用的规则创建thinking对象
                    intermediate_steps = event.metadata.get("intermediate_steps", [])
                    if intermediate_steps:
                        for step_info in intermediate_steps:
                            rule = step_info.get('rule', '系统自动操作')
                            # 为每个中间步骤发送thinking对象
                            state['thinking'] = {
                                'step': event.step,
                                'thinking': '',
                                'command': rule,
                                'explanation': '自动优化步骤',
                                'is_final': False,
                            }
                            # 同时记录到命令日志
                            state['commands'] += f"[自动] {rule}\n"
                            
                            await queue.put(send_update())
                    
                    # 处理主结果
                    # 检查是否是状态切换命令
                    is_state_switch = False
                    cmd_lower = last_command.lower().strip() if last_command else ''
                    for sc in state_switch_commands:
                        if cmd_lower.startswith(sc.rstrip(':')) or cmd_lower == sc.rstrip(':'):
                            is_state_switch = True
                            break
                    
                    if is_state_switch:
                        # 状态切换命令：单独成块显示
                        results_lines.append(f"({last_command})")
                        results_lines.append(f"= {event.content}")
                    elif last_command:
                        # 普通命令：显示在同一行
                        line = f"= {event.content} ({last_command})"
                        results_lines.append(line)
                    else:
                        line = f"  {event.content}"
                        results_lines.append(line)
                    state['results'] = "\n".join(results_lines) + "\n"
                    
                elif event.type == EventType.ERROR:
                    state['errors'] += f"[步骤 {event.step}] {event.content}\n"
                    
                elif event.type == EventType.COMPLETE:
                    # 只有当真正完成时才显示"求解完成"
                    if executor.is_finished():
                        state['results'] += "\n✓ 求解完成\n"
                    else:
                        # LLM认为完成但实际未完成
                        state['results'] += "\n⚠ LLM判断已完成，但证明尚未结束\n"

                # 每次状态更新都推送到队列
                await queue.put(send_update())
            
            # 创建求解任务
            async def solve_task():
                return await solver.solve(expression, callback, cond_list, user_instruction=instruction)
            
            task = asyncio.ensure_future(solve_task())
            
            # 队列消费循环
            while not task.done():
                if _stop_flag:
                    task.cancel()
                    state['results'] += "\n⚠ 已终止求解\n"
                    yield send_update()
                    break

                try:
                    # 等待队列消息，超时后发送心跳
                    # 注意：在 sync generator 中我们要用 run_until_complete 来等待 async 函数
                    data = loop.run_until_complete(asyncio.wait_for(queue.get(), timeout=0.3))
                    yield data
                    
                    # 尝试快速消费积压的消息
                    while not queue.empty():
                         data = loop.run_until_complete(queue.get())
                         yield data

                except asyncio.TimeoutError:
                    yield send_update()
                except Exception:
                    break
            
            # 获取最终结果
            if not _stop_flag:
                try:
                    result = loop.run_until_complete(task)
                    if not result.success and result.error:
                        state['errors'] += f"\n求解失败: {result.error}\n"
                except asyncio.CancelledError:
                    pass
            
            yield send_update()
            yield "event: complete\ndata: {}\n\n"
            
        except Exception as e:
            import traceback
            state['errors'] = f"错误: {str(e)}\n{traceback.format_exc()}"
            yield send_update()
            yield "event: complete\ndata: {}\n\n"
        finally:
            loop.close()
    
    return Response(generate(), mimetype='text/event-stream')


def main(host='127.0.0.1', port=7860, debug=False):
    """启动服务器"""
    print(f"🚀 LLM-IsCalc 服务器启动")
    print(f"📍 访问: http://{host}:{port}")
    print(f"🤖 模型: {_config.llm.model}")
    app.run(host=host, port=port, debug=debug, threaded=True)


if __name__ == "__main__":
    main()
