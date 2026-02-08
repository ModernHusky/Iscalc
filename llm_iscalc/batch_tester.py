"""批量测试模块

使用 SolverLoop 进行批量测试，支持并发和实时事件推送。
"""

import asyncio
import time
import queue
import re
from pathlib import Path
from typing import List, Dict, Any, Optional
from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor
from threading import Lock
from datetime import datetime
from collections import defaultdict

from .solver_loop import SolverLoop, SolveEvent, EventType
from .llm_engine import LLMEngine
from .command_executor import CommandExecutor
from .config import default_config
from .logger_config import get_phase_logger
from .models import TokenUsage, PromptComponent, RoundLog

# 添加 integral 模块导入以支持完整的上下文解析
import sys
from pathlib import Path as PathLib
sys.path.insert(0, str(PathLib(__file__).parent.parent))
from integral.context import Context
from integral import parser, action, expr


@dataclass
class ProblemInfo:
    """问题信息"""
    filename: str
    index: int
    problem: str
    context: Any = None  # Context 对象
    correct_answer: List[str] = None  # 正确答案步骤
    pre_problems_str: str = ""  # 前置问题
    pre_definitions_str: str = ""  # 前置定义


class BatchTester:
    """批量测试器
    
    使用 SolverLoop 进行批量测试，支持并发和实时事件推送。
    """
    
    def __init__(self, theories: List[str], max_workers: int = 5, max_step: int = 25):
        """初始化批量测试器
        
        Args:
            theories: 理论文件名列表
            max_workers: 最大并发数
            max_step: 每个题目的最大步数
        """
        self.theories = theories
        self.max_workers = max_workers
        self.max_step = max_step
        self.event_queue = queue.Queue()
        self.stop_flag = False
        self.logger = get_phase_logger(__name__)
        self._lock = Lock()
        self._futures = []  # 保存所有 futures 以便取消
        
        # 日志目录和统计
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        self.log_dir = Path(f"llm_iscalc/test_logs/{timestamp}")
        self.prompt_analytics = defaultdict(lambda: {"count": 0, "total_tokens": 0})
        
    def _parse_theory_file(self, theory_name: str) -> List[ProblemInfo]:
        """解析理论文件，提取所有问题
        
        参考旧版本 process_file() 的完整解析逻辑
        
        Args:
            theory_name: 理论文件名（不含 .thy 扩展名）
            
        Returns:
            问题信息列表
        """
        theory_file = Path(__file__).parent.parent / "theories" / f"{theory_name}.thy"
        
        if not theory_file.exists():
            self.logger.warning(f"理论文件不存在: {theory_file}")
            return []
        
        result = []
        problem_pattern = r"(prove|calculate) (\[.*\] )?(.+)"
        
        with open(theory_file, 'r', encoding='utf-8') as f:
            content = f.read()
            lines = [s for s in content.split('\n') if s.strip()]
        
        cur_goal = None
        steps = []
        i = 0
        ctx = Context()
        ctx.load_book("base")
        pre_problems, pre_definitions = [], []
        
        # 使用 for 循环遍历，确保有明确的终止条件
        for line in lines:
            line = line.strip()
            if line.startswith('#') or line.startswith('//'):
                # title or comment
                continue
            
            a = parser.parse_action(line)
            
            # 处理 imports
            if isinstance(a, action.ImportsAction):
                for theory in a.theories:
                    if theory != 'base':
                        ctx.load_book(theory)
            
            # 处理 prove/calculate
            if isinstance(a, (action.ProveAction, action.CalculateAction)):
                if cur_goal:
                    # 创建问题信息（不包含当前目标）
                    result.append(ProblemInfo(
                        filename=theory_name,
                        index=i,
                        problem=problem,
                        context=Context(ctx),  # 复制上下文
                        correct_answer=steps[:],  # 复制步骤
                        pre_problems_str="\n".join(pre_problems),
                        pre_definitions_str="\n".join(pre_definitions)
                    ))
                    
                    # 更新 pre_problems
                    match = re.search(problem_pattern, problem)
                    if match:
                        pre_problems.append(match.group(3).strip())
                    
                    # 创建新的上下文（继承当前上下文）
                    ctx = Context(ctx)
                    
                    # 将当前目标添加到上下文
                    if isinstance(cur_goal, action.ProveAction):
                        if cur_goal.expr.is_equals() and expr.is_indefinite_integral(cur_goal.expr.lhs):
                            ctx.add_indefinite_integral(cur_goal.expr, cur_goal.conditions, cur_goal.attrs)
                        elif cur_goal.expr.is_equals() and expr.is_integral(cur_goal.expr.lhs):
                            ctx.add_definite_integral(cur_goal.expr, cur_goal.conditions, cur_goal.attrs)
                        else:
                            ctx.add_other_identities(cur_goal.expr, cur_goal.conditions, cur_goal.attrs)
                
                cur_goal = a
                problem = line
                steps = []
                i += 1
            
            elif line and line != 'done':
                steps.append(line)
            
            # 处理 define
            if isinstance(a, action.DefineAction):
                ctx.add_definition(a.expr, a.conditions)
                s = str(a.expr)
                if a.conditions:
                    s = s + " for " + str(a.conditions)
                pre_definitions.append(s)
        
        # 添加最后一个问题
        if cur_goal:
            result.append(ProblemInfo(
                filename=theory_name,
                index=i,
                problem=problem,
                context=ctx,
                correct_answer=steps,
                pre_problems_str="\n".join(pre_problems),
                pre_definitions_str="\n".join(pre_definitions)
            ))
        
        self.logger.info(f"从 {theory_name}.thy 中解析出 {len(result)} 个问题")
        return result
    
    def _get_all_problems(self) -> List[ProblemInfo]:
        """获取所有理论文件的问题列表"""
        all_problems = []
        for theory in self.theories:
            problems = self._parse_theory_file(theory)
            all_problems.extend(problems)
        return all_problems
    
    def _save_round_log(self, problem_id: str, round_log: RoundLog):
        """保存单轮交互日志"""
        log_dir = self.log_dir / "temp" / problem_id
        log_dir.mkdir(parents=True, exist_ok=True)
        
        log_file = log_dir / f"round_{round_log.round}.md"
        
        # 更新提示词分析统计
        for component in round_log.prompt_components:
            self.prompt_analytics[component.name]["count"] += 1
            self.prompt_analytics[component.name]["total_tokens"] += component.token_count
        
        # 生成 Markdown 内容
        content = f"""# {problem_id} - Round {round_log.round}

**Timestamp**: {round_log.timestamp}

## 📊 Token & Cost

| Metric | Value |
|---|---|
| Prompt Tokens (Hit) | {round_log.token_usage.prompt_cache_hit_tokens} |
| Prompt Tokens (Miss) | {round_log.token_usage.prompt_cache_miss_tokens} |
| Completion Tokens | {round_log.token_usage.completion_tokens} |
| **Total Tokens** | **{round_log.token_usage.total_tokens}** |
| **Cost** | **￥{round_log.token_usage.cost:.6f}** |

### Prompt Components

| Component | Tokens | Percentage |
|---|---|---|
"""
        total_comp_tokens = sum(c.token_count for c in round_log.prompt_components)
        for comp in round_log.prompt_components:
            pct = (comp.token_count / total_comp_tokens * 100) if total_comp_tokens > 0 else 0
            content += f"| {comp.name} | {comp.token_count} | {pct:.1f}% |\n"
        
        content += f"""
---

## 🤖 LLM Interaction

### Prompt


{round_log.prompt_text}


### Response


{round_log.response_text}


---

## 🔄 State Change

### Before

```latex
{round_log.state_before}
```

### After

```latex
{round_log.state_after}
```
"""
        
        log_file.write_text(content, encoding='utf-8')
    
    def _save_output(self, problem_id: str, problem_text: str, result):
        """保存成功的题目"""
        output_dir = self.log_dir / "outputs"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        output_file = output_dir / f"{problem_id}.md"
        
        content = f"""# {problem_id} - SUCCESS ✅

## Problem

```
{problem_text}
```

## Solution

**Iterations**: {result.iterations}
**Total Cost**: ￥{result.total_token_usage.cost:.6f}
**Total Tokens**: {result.total_token_usage.total_tokens}

### Steps

"""
        for i, step in enumerate(result.steps, 1):
            content += f"{i}. {step}\n"
        
        content += f"""
## Result

```latex
{result.latex_solution}
```

## Token Statistics

| Metric | Value |
|---|---|
| Prompt Tokens (Hit) | {result.total_token_usage.prompt_cache_hit_tokens} |
| Prompt Tokens (Miss) | {result.total_token_usage.prompt_cache_miss_tokens} |
| Completion Tokens | {result.total_token_usage.completion_tokens} |
| **Total Tokens** | **{result.total_token_usage.total_tokens}** |
| **Total Cost** | **￥{result.total_token_usage.cost:.6f}** |
"""
        
        output_file.write_text(content, encoding='utf-8')
    
    def _save_error(self, problem_id: str, problem_text: str, error_message: str, result):
        """保存失败的题目"""
        error_dir = self.log_dir / "errors"
        error_dir.mkdir(parents=True, exist_ok=True)
        
        error_file = error_dir / f"{problem_id}.md"
        
        content = f"""# {problem_id} - FAILED ❌

## Problem

```
{problem_text}
```

## Error

{error_message}

## Statistics

**Iterations**: {result.iterations if result else 0}
**Total Cost**: ￥{result.total_token_usage.cost if result else 0:.6f}
**Total Tokens**: {result.total_token_usage.total_tokens if result else 0}
"""
        
        error_file.write_text(content, encoding='utf-8')
    
    async def _test_single_problem(self, problem: ProblemInfo):
        """测试单个问题
        
        Args:
            problem: 问题信息
        """
        # 检查是否已停止
        if self.stop_flag:
            self.logger.info(f"检测到停止信号，跳过问题 {problem.filename}_{problem.index}")
            return
        
        problem_id = f"{problem.filename}_{problem.index}"
        start_time = time.time()
        
        # 发送题目开始事件
        self.event_queue.put({
            "type": "problem_start",
            "data": {
                "problem_id": problem_id,
                "problem_text": problem.problem
            }
        })
        
        # 初始化求解器组件（每个问题独立的实例）
        llm_engine = LLMEngine(default_config.llm)
        # 传入问题的完整上下文
        if problem.context is not None:
            executor = CommandExecutor(context=problem.context)
        else:
            executor = CommandExecutor(base_theory=default_config.iscalc.base_theory)
        solver = SolverLoop(llm_engine, executor, default_config.solver)
        
        # 覆盖最大迭代次数
        solver.config.max_iterations = self.max_step
        
        # 记录步骤
        steps = []
        current_round = 0
        
        # 定义回调函数
        async def callback(event: SolveEvent):
            nonlocal current_round
            
            # 检查停止标志
            if self.stop_flag:
                raise Exception("用户停止测试")
            
            # 将 SolveEvent 转换为我们的事件格式
            if event.type == EventType.THINKING:
                # LLM 思考流式输出
                self.event_queue.put({
                    "type": "thinking_stream",
                    "data": {
                        "problem_id": problem_id,
                        "chunk": event.content
                    }
                })
            
            elif event.type in [EventType.COMMAND, EventType.COMMAND_SUCCESS]:
                # 命令执行
                current_round = event.step
                # round_update 事件已移除以避免卡死
                
                # 记录步骤
                steps.append({
                    "round": current_round,
                    "command": event.content,
                    "result": event.latex or ""
                })
            
            elif event.type == EventType.ERROR:
                # 错误
                self.event_queue.put({
                    "type": "thinking_stream",
                    "data": {
                        "problem_id": problem_id,
                        "chunk": f"\n❌ 错误: {event.content}\n"
                    }
                })
        
        # 调用 SolverLoop.solve
        try:
            result = await solver.solve(
                expression=problem.problem,
                callback=callback,
                conditions=None,
                user_instruction=None
            )
            
            # 保存每轮交互日志
            for round_log in result.round_logs:
                self._save_round_log(problem_id, round_log)
            
            # 保存最终结果
            if result.success:
                status = "success"
                self._save_output(problem_id, problem.problem, result)
            else:
                status = "failed"
                error_message = result.error or "Unknown error"
                self._save_error(problem_id, problem.problem, error_message, result)
            
            num_rounds = result.iterations
            error_message = result.error if not result.success else None
            
        except Exception as e:
            self.logger.error(f"测试 {problem_id} 时发生异常: {e}")
            status = "error"
            num_rounds = current_round
            error_message = str(e)
            steps = []
        
        finally:
            # 关闭 LLM 引擎
            await llm_engine.close()
        
        # 发送完成事件
        end_time = time.time()
        self.event_queue.put({
            "type": "problem_complete",
            "data": {
                "problem_id": problem_id,
                "status": status,
                "num_rounds": num_rounds,
                "time_elapsed": end_time - start_time,
                "problem_text": problem.problem,
                "steps": steps,
                "error_message": error_message,
                # Token 信息
                "token_usage": {
                    "total_tokens": result.total_token_usage.total_tokens if 'result' in locals() else 0,
                    "cost": result.total_token_usage.cost if 'result' in locals() else 0,
                    "prompt_cache_hit": result.total_token_usage.prompt_cache_hit_tokens if 'result' in locals() else 0,
                    "prompt_cache_miss": result.total_token_usage.prompt_cache_miss_tokens if 'result' in locals() else 0,
                    "completion": result.total_token_usage.completion_tokens if 'result' in locals() else 0
                }
            }
        })
    
    def _run_async_test(self, problem: ProblemInfo):
        """在新的事件循环中运行异步测试（用于线程池）"""
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)
        try:
            loop.run_until_complete(self._test_single_problem(problem))
        finally:
            loop.close()
    
    def run_tests(self):
        """运行批量测试（在独立线程中执行）"""
        try:
            # 获取所有问题
            all_problems = self._get_all_problems()
            total = len(all_problems)
            
            if total == 0:
                self.event_queue.put({
                    "type": "error",
                    "data": {"message": "没有找到任何问题"}
                })
                return
            
            # 发送测试开始事件
            self.event_queue.put({
                "type": "test_start",
                "data": {"total_problems": total}
            })
            
            # 使用线程池并发测试
            with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
                self._futures = []  # 重置 futures 列表
                for problem in all_problems:
                    if self.stop_flag:
                        break
                    future = executor.submit(self._run_async_test, problem)
                    self._futures.append(future)
                
                # 等待所有任务完成
                for future in self._futures:
                    if self.stop_flag:
                        break
                    try:
                        future.result()
                    except Exception as e:
                        self.logger.error(f"任务执行异常: {e}")
            
            # 发送测试完成事件
            if not self.stop_flag:
                # 生成最终报告
                self._generate_final_report()
                
                self.event_queue.put({
                    "type": "test_complete",
                    "data": {}
                })
        
        except Exception as e:
            self.logger.error(f"批量测试异常: {e}")
            self.event_queue.put({
                "type": "error",
                "data": {"message": str(e)}
            })
    
    def _generate_final_report(self):
        """生成最终测试报告"""
        report_file = self.log_dir / "report.md"
        
        content = f"""# 批量测试报告

**生成时间**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

## 🧠 提示词分析

### 成分调用频率 (Top 10)

| 成分名称 | 调用次数 | 占比 |
|---|---|---|
"""
        
        sorted_by_count = sorted(self.prompt_analytics.items(), 
                                 key=lambda x: x[1]["count"], reverse=True)[:10]
        total_count = sum(v["count"] for v in self.prompt_analytics.values())
        
        for name, stats in sorted_by_count:
            ratio = (stats["count"] / total_count * 100) if total_count > 0 else 0
            content += f"| {name} | {stats['count']} | {ratio:.1f}% |\n"
        
        content += """
### 成分 Token 消耗 (Top 10)

| 成分名称 | 总 Tokens | 占比 |
|---|---|---|
"""
        
        sorted_by_tokens = sorted(self.prompt_analytics.items(), 
                                  key=lambda x: x[1]["total_tokens"], reverse=True)[:10]
        total_tokens = sum(v["total_tokens"] for v in self.prompt_analytics.values())
        
        for name, stats in sorted_by_tokens:
            ratio = (stats["total_tokens"] / total_tokens * 100) if total_tokens > 0 else 0
            content += f"| {name} | {stats['total_tokens']:,} | {ratio:.1f}% |\n"
        
        report_file.write_text(content, encoding='utf-8')
        self.logger.info(f"最终报告已保存到: {report_file}")
    
    def stop(self):
        """停止测试"""
        self.stop_flag = True
        self.logger.info("收到停止信号")
        
        # 尝试取消所有未完成的任务
        if hasattr(self, '_futures'):
            for future in self._futures:
                if not future.done():
                    future.cancel()
                    self.logger.info(f"已取消一个未完成的任务")
