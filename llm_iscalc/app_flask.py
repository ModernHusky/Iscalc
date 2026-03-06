"""Flask 前端界面 - 简化版本"""

import asyncio
import json
import re
import datetime
from pathlib import Path
from flask import Flask, request, Response, jsonify, send_from_directory

from .config import create_config
from .llm_engine import LLMEngine
from .command_executor import CommandExecutor
from .solver_loop import SolverLoop, SolveEvent, EventType
from .logger_config import configure_logging, get_phase_logger
from .batch_tester import BatchTester

# 配置日志
configure_logging()
logger = get_phase_logger(__name__)

# 获取当前文件所在目录
_current_dir = Path(__file__).parent

# 创建 Flask 应用
app = Flask(__name__, 
            static_folder=str(_current_dir / 'static'),
            static_url_path='/static')

# 全局状态
_config = create_config()

_stop_flag = False
_batch_tester = None  # 当前运行的批量测试器


def _write_solve_log(expression: str, thinking_history: list, commands: str, results: str, errors: str, is_success: bool):
    """将本次求解过程写入日志文件到 test_logs/ 目录。

    文件名格式: 表达式(空格→_)_YYYYMMDDHHmm_pass/error.log
    """
    try:
        log_dir = _current_dir / 'test_logs'
        log_dir.mkdir(parents=True, exist_ok=True)

        # 生成文件名：表达式中空格→_，去掉非法字符
        expr_part = expression.replace(' ', '_')
        expr_part = re.sub(r'[<>:"/\\|?*]', '', expr_part)  # 去掉 Windows 非法字符
        expr_part = expr_part[:80]  # 限制长度避免路径过长

        now = datetime.datetime.now()
        timestamp = now.strftime('%Y%m%d%H%M')
        status = 'pass' if is_success else 'error'
        filename = f"{expr_part}_{timestamp}_{status}.log"

        # 大节分隔符（虚线）
        def _sec(title: str, width: int = 55) -> str:
            """生成形如  "--- 标题 ---"  的分隔行"""
            inner = f' {title} '
            pad = max(width - len(inner), 6)
            left = pad // 2
            right = pad - left
            return '-' * left + inner + '-' * right

        # 步骤内小分隔
        step_hr = '-' * 60 + '\n'

        # ── AI 思考过程 ──
        thinking_lines = [_sec('AI 思考过程'), '']
        if thinking_history:
            for i, t in enumerate(thinking_history, 1):
                thinking_lines.append(f"[步骤 {t.get('step', i)}]")
                thinking_lines.append(f"  thinking   : {t.get('thinking', '').strip()}")
                thinking_lines.append(f"  command    : {t.get('command', '').strip()}")
                thinking_lines.append(f"  explanation: {t.get('explanation', '').strip()}")
                thinking_lines.append(f"  is_final   : {t.get('is_final', False)}")
                thinking_lines.append(step_hr.rstrip('\n'))
        else:
            thinking_lines.append('（无思考过程记录）')

        # ── 执行命令完整过程 ──
        commands_lines = [
            _sec('执行命令完整过程'), '',
            commands.strip() if commands.strip() else '（无命令记录）'
        ]

        # ── 结果 ──
        results_lines_log = [
            _sec('结果'), '',
            results.strip() if results.strip() else '（无结果）'
        ]

        # ── 报错 ──
        errors_lines = [
            _sec('报错'), '',
            errors.strip() if errors.strip() else '（无报错）'
        ]

        header = (
            f"表达式 : {expression}\n"
            f"求解时间: {now.strftime('%Y-%m-%d %H:%M')}\n"
            f"状态   : {'PASS ✓' if is_success else 'ERROR ✗'}\n"
        )

        content = (
            header
            + '\n'.join(thinking_lines) + '\n'
            + '\n'.join(commands_lines) + '\n'
            + '\n'.join(results_lines_log) + '\n'
            + '\n'.join(errors_lines) + '\n'
        )

        log_path = log_dir / filename
        log_path.write_text(content, encoding='utf-8')
        logger.info(f"[日志] 已写入求解日志: {log_path}")

    except Exception as e:
        logger.warning(f"[日志] 写入日志失败: {e}")


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


@app.route('/api/test/theories')
def get_theories():
    """获取所有理论文件列表"""
    import os
    theories_dir = _current_dir.parent / 'theories'
    
    if not theories_dir.exists():
        return jsonify({'error': '理论文件目录不存在'}), 404
    
    # 获取所有 .thy 文件
    theory_files = []
    for file in sorted(theories_dir.glob('*.thy')):
        theory_name = file.stem  # 不带扩展名的文件名
        theory_files.append(theory_name)
    
    return jsonify(theory_files)


@app.route('/api/test/run')
def run_batch_test():
    """运行批量测试（SSE 流式输出）"""
    import threading
    
    theories = request.args.get('theories', '').strip()
    max_workers = int(request.args.get('max_workers', '5'))
    max_step = 25
    
    if not theories:
        return jsonify({'error': '理论列表不能为空'}), 400
    
    theory_list = [t.strip() for t in theories.split(',') if t.strip()]
    
    def generate():
        """生成 SSE 事件流"""
        global _batch_tester
        
        # 创建批量测试器
        _batch_tester = BatchTester(theory_list, max_workers=max_workers, max_step=max_step)
        
        # 在独立线程中运行异步事件循环
        def run_in_thread():
            loop = asyncio.new_event_loop()
            asyncio.set_event_loop(loop)
            try:
                loop.run_until_complete(_batch_tester.run_tests_async())
            finally:
                loop.close()
        
        test_thread = threading.Thread(target=run_in_thread)
        test_thread.start()
        
        # 从事件队列中读取事件并发送
        try:
            while True:
                try:
                    # 等待事件（超时 1 秒）
                    event = _batch_tester.event_queue.get(timeout=1)
                    
                    # 发送 SSE 事件
                    event_type = event['type']
                    event_data = json.dumps(event['data'], ensure_ascii=False)
                    yield f"event: {event_type}\ndata: {event_data}\n\n"
                    
                    # 如果是测试完成或错误，退出循环
                    if event_type in ['test_complete', 'error']:
                        break
                
                except Exception as e:
                    # 队列超时，检查线程是否还在运行
                    if not test_thread.is_alive():
                        break
        
        finally:
            # 等待测试线程结束
            test_thread.join(timeout=5)
            _batch_tester = None
    
    return Response(generate(), mimetype='text/event-stream')


@app.route('/api/test/stop', methods=['POST'])
def stop_batch_test():
    """停止批量测试"""
    global _batch_tester
    
    if _batch_tester:
        _batch_tester.stop()
        return jsonify({'status': 'ok', 'message': '已发送停止信号'})
    else:
        return jsonify({'status': 'ok', 'message': '没有正在运行的测试'})


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
        # 收集每轮 thinking 历史，用于日志（按步骤号索引，独立于 UI state）
        step_thinking_log = {}  # {step_num: {thinking, command, explanation, is_final}}
        
        def send_update():
            """发送更新"""
            data = {
                'commands': state['commands'],
                'results': state['results'],
                'errors': state['errors'],
                'thinking': state['thinking']
            }
            if 'command_success' in state:
                data['command_success'] = state['command_success']
                del state['command_success']  # 发送后清除，避免重复发送
            if 'command_failure' in state:
                data['command_failure'] = state['command_failure']
                del state['command_failure']  # 发送后清除，避免重复发送
            return f"data: {json.dumps(data, ensure_ascii=False)}\n\n"
        
        # 创建事件循环
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)
        
        # 创建消息队列
        queue = asyncio.Queue()
        
        try:
            # 创建 executor 和 solver
            llm_engine = LLMEngine(_config.llm)
            executor = CommandExecutor(_config.iscalc.base_theory)
            solver = SolverLoop(llm_engine, executor, _config.solver)
            
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
                    # 如果步骤号变化，重置 UI 思考状态
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
                        
                        # ── 续写事件：从 metadata 提取旧思考内容 ──
                        # 当 solver_loop 续写时，content = continuation_text（单 JSON），
                        # 旧思考通过 metadata["pre_skill_display"] 传递
                        _pre_skill_thinking = ""
                        pre_skill_raw = event.metadata.get("pre_skill_display", "")
                        if pre_skill_raw:
                            # 从旧 display_text 中提取纯思考文本（去掉 JSON 语法包装）
                            psd_match = re.search(r'"thinking"\s*:\s*"', pre_skill_raw)
                            if psd_match:
                                psd_start = psd_match.end()
                                psd_end = psd_start
                                while psd_end < len(pre_skill_raw):
                                    if pre_skill_raw[psd_end] == '\\' and psd_end + 1 < len(pre_skill_raw):
                                        psd_end += 2
                                    elif pre_skill_raw[psd_end] == '"':
                                        break
                                    else:
                                        psd_end += 1
                                _pre_skill_thinking = pre_skill_raw[psd_start:psd_end]
                                _pre_skill_thinking = _pre_skill_thinking.replace('\\n', '\n').replace('\\"', '"').replace('\\\\', '\\').strip()
                            else:
                                # 旧 display_text 不是 JSON 格式，直接用作前缀
                                _pre_skill_thinking = pre_skill_raw.strip()
                        
                        # 预处理：去除可能存在的 Markdown 代码块标记
                        clean_content = content
                        if "```" in content:
                           clean_content = re.sub(r'```(?:json)?', '', content).strip()
                        
                        # ============ 流式 thinking 解析优化 ============
                        # 策略1：用 rfind 尝试从最后一个 "thinking": " 之后提取内容
                        # 关键改进：Search-o1 技能加载时会有多个 JSON 拼接，rfind 能始终提取最新生成的思考
                        thinking_start = clean_content.rfind('"thinking"')
                        prefix_text = ""  # 最新 JSON 之前的文本
                        last_json_start = 0  # 记录最新 JSON 的起始位置
                        
                        if thinking_start != -1:
                            # 提取最新 JSON 之前的文本作为前缀
                            # 这包括旧的思考过程和技能加载确认标记
                            json_start = clean_content.rfind('{', 0, thinking_start)
                            if json_start >= 0:
                                last_json_start = json_start
                                
                            if json_start > 0:
                                raw_prefix = clean_content[:json_start].strip()
                                # 清洁前缀：如果前缀中包含旧 JSON 的 "thinking" 值，提取纯文本，去掉 JSON 语法残留
                                pf_match = re.search(r'"thinking"\s*:\s*"', raw_prefix)
                                if pf_match:
                                    pf_start = pf_match.end()
                                    pf_end = pf_start
                                    while pf_end < len(raw_prefix):
                                        if raw_prefix[pf_end] == '\\' and pf_end + 1 < len(raw_prefix):
                                            pf_end += 2
                                        elif raw_prefix[pf_end] == '"':
                                            break
                                        else:
                                            pf_end += 1
                                    extracted = raw_prefix[pf_start:pf_end]
                                    prefix_text = extracted.replace('\\n', '\n').replace('\\"', '"').replace('\\\\', '\\').strip()
                                else:
                                    prefix_text = raw_prefix
                            
                            # 找到最新 thinking 字段的值开始位置
                            colon_pos = clean_content.find(':', thinking_start)
                            if colon_pos != -1:
                                # 跳过冒号后的空格和引号
                                value_start = colon_pos + 1
                                while value_start < len(clean_content) and clean_content[value_start] in ' \t\n':
                                    value_start += 1
                                
                                if value_start < len(clean_content) and clean_content[value_start] == '"':
                                    value_start += 1  # 跳过开头的引号
                                    
                                    # 查找结束引号（需要处理转义）
                                    value_end = value_start
                                    while value_end < len(clean_content):
                                        char = clean_content[value_end]
                                        if char == '\\' and value_end + 1 < len(clean_content):
                                            value_end += 2  # 跳过转义字符
                                        elif char == '"':
                                            break  # 找到结束引号
                                        else:
                                            value_end += 1
                                    
                                    # 提取 thinking 内容（即使不完整）
                                    thinking_raw = clean_content[value_start:value_end]
                                    # 处理转义序列
                                    thinking_value = thinking_raw.replace('\\n', '\n').replace('\\"', '"').replace('\\\\', '\\')
                                    
                                    # 🔧 关键修复：合并前缀文本和最新 JSON thinking 字段
                                    # 优先使用 _pre_skill_thinking（从 metadata 提取的干净旧思考）
                                    effective_prefix = _pre_skill_thinking or prefix_text
                                    if effective_prefix:
                                        state['thinking']['thinking'] = effective_prefix + "\n\n" + thinking_value
                                    else:
                                        state['thinking']['thinking'] = thinking_value
                        
                        # 限定后续字段仅在最新的 JSON 内查找，防止意外匹配到旧的 JSON 属性
                        search_content = clean_content[last_json_start:]
                        
                        # 策略2：正则解析其他字段（在最新 JSON 范围内）
                        command_match = re.search(r'"command"\s*:\s*"([^"]*)(?:"|$)', search_content)
                        if command_match:
                            cmd_value = command_match.group(1)
                            # 过滤技能加载标记和 read_skill 命令（这些不应显示在"命令"字段）
                            if not cmd_value.startswith('<|load_skill|>') and not cmd_value.startswith('read_skill'):
                                state['thinking']['command'] = cmd_value
                            
                        explanation_match = re.search(r'"explanation"\s*:\s*"((?:[^"\\]|\\.)*)(?:"|$)', search_content, re.DOTALL)
                        if explanation_match:
                            val = explanation_match.group(1).replace('\\n', '\n').replace('\\"', '"')
                            state['thinking']['explanation'] = val
                            
                        is_final_match = re.search(r'"is_final"\s*:\s*(true|false)', search_content)
                        if is_final_match:
                            state['thinking']['is_final'] = is_final_match.group(1) == 'true'

                        # 策略3：如果没有找到 "thinking" 字段（如 Search-o1 的纯文本前缀），直接显示原始内容
                        if thinking_start == -1:
                            display_text_raw = clean_content.rstrip().rstrip('{').strip()
                            if _pre_skill_thinking:
                                if display_text_raw:
                                    state['thinking']['thinking'] = _pre_skill_thinking + "\n\n" + display_text_raw
                                else:
                                    state['thinking']['thinking'] = _pre_skill_thinking
                            else:
                                stripped = clean_content.strip()
                                is_incomplete_json_fragment = len(stripped) < 5 and (stripped.endswith('{') or stripped == '{')
                                if not clean_content.startswith('[DEBUG') and not is_incomplete_json_fragment:
                                    state['thinking']['thinking'] = clean_content

                    except Exception as e:
                        import traceback
                        print(f"DEBUG EXCEPTION in parse_response: {e}")
                        traceback.print_exc()
                        print(f"CONTENT WAS: {repr(content)}")
                        pass

                    # ── 日志独立积累（不受 UI state 覆盖影响）──
                    s = event.step
                    if s not in step_thinking_log:
                        step_thinking_log[s] = {
                            'step': s, 'thinking': '', 'command': '',
                            'explanation': '', 'is_final': False
                        }
                    log_entry = step_thinking_log[s]
                    # 取当前解析出的最新值（非空则覆盖）
                    cur = state['thinking']
                    if cur.get('thinking'):
                        log_entry['thinking'] = cur['thinking']
                    if cur.get('command') and not cur['command'].startswith('<|load_skill|'):
                        log_entry['command'] = cur['command']
                    if cur.get('explanation'):
                        log_entry['explanation'] = cur['explanation']
                    log_entry['is_final'] = cur.get('is_final', False)

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
                                'explanation': f'系统自动执行: {rule}',
                                'is_final': False,
                            }
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
                    is_really_finished = executor.is_finished()
                    current_state_name = executor.get_current_state_name()
                    
                    if is_really_finished:
                        state['results'] += "\n✓ 求解完成\n"
                    elif current_state_name == "CALCULATE":
                         # 计算模式下，用户没输入done，但LLM认为完成了，通常是计算出了结果
                         # 这种情况下显示完成，而不是告警
                         state['results'] += "\n✓ 计算完成 (LLM判定)\n"
                    else:
                        # LLM认为完成但实际未完成 (主要针对证明模式)
                        state['results'] += "\n⚠ LLM判断已完成，但证明尚未结束\n"

                elif event.type == EventType.SKILL_LOAD:
                    # 技能加载信息不再写入 commands 字段，仅通过前端右上角通知展示
                    # 保留事件处理以便日后扩展（如记录日志或发送专门的 skill_load 字段）
                    pass

                elif event.type == EventType.COMMAND_SUCCESS:
                    # 2026-02-01 架构重构：命令执行成功事件
                    # 将成功信息发送给前端，前端负责添加到 sidebar
                    state['command_success'] = {
                        'content': event.content,
                        'step': event.step,
                        'changed': event.metadata.get("changed", False),
                        'explanation': event.metadata.get("explanation", "")
                    }
                
                elif event.type == EventType.COMMAND_FAILURE:
                    # 2026-02-01 架构重构：命令执行失败事件
                    # 将失败信息发送给前端，前端负责不添加/删除
                    state['command_failure'] = {
                        'content': event.content,
                        'step': event.step,
                        'error': event.metadata.get("error", ""),
                        'explanation': event.metadata.get("explanation", "")
                    }

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

            # ── 写日志 ──
            # 从独立积累的 step_thinking_log 中按步骤顺序构建 thinking_history
            # 过滤：只保留有实质内容的步骤（thinking 或 command 非空）
            thinking_history = [
                step_thinking_log[k]
                for k in sorted(step_thinking_log.keys())
                if step_thinking_log[k].get('thinking') or step_thinking_log[k].get('command')
            ]

            is_success = '✓' in state['results']
            _write_solve_log(
                expression=expression,
                thinking_history=thinking_history,
                commands=state['commands'],
                results=state['results'],
                errors=state['errors'],
                is_success=is_success,
            )

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
