"""Command Executor 模块

负责执行Iscalc命令并管理计算状态。
"""

import sys
import os
from dataclasses import dataclass
from typing import Optional, List, Any
from enum import Enum

# 添加integral路径到sys.path
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__))))

from integral import parser
from integral.compstate import CompFile, Calculation
from integral import state as state_module
from integral.action import CalculateAction


class ErrorType(Enum):
    """错误类型枚举"""
    SYNTAX = "syntax"
    RULE = "rule"
    SYSTEM = "system"
    NONE = "none"


@dataclass
class ExecutionResult:
    """命令执行结果"""
    success: bool
    result: Optional[str] = None
    latex_result: Optional[str] = None
    error: Optional[str] = None
    error_type: ErrorType = ErrorType.NONE
    changed: bool = False
    raw_expr: Optional[Any] = None
    intermediate_steps: List[dict] = None  # 中间步骤列表


class CommandExecutor:
    """命令执行器
    
    使用Iscalc的状态机来执行命令。
    工作流程:
    1. initialize() - 创建CompFile和InitialState，然后用CalculateAction进入CalculateState
    2. execute() - 在CalculateState中执行RuleAction
    """
    
    def __init__(self, base_theory: str = "base"):
        self.base_theory = base_theory
        self.comp_file: Optional[CompFile] = None
        self.state: Optional[Any] = None
        self.history: List[dict] = []
        self.initial_expression: Optional[str] = None
        self._previous_expr_str: Optional[str] = None
    
    def initialize(self, expression: str, conditions: Optional[List[str]] = None) -> ExecutionResult:
        """初始化计算会话
        
        根据输入的命令（prove/calculate）创建对应的状态。
        """
        try:
            # 创建CompFile
            self.comp_file = CompFile(self.base_theory, "llm_session")
            
            # 创建InitialState
            initial_state = state_module.InitialState(self.comp_file.ctx)
            
            # 解析条件 (带智能预处理)
            cond_list = None
            if conditions:
                raw_conditions = conditions
                cond_list = []
                for c in raw_conditions:
                    c = c.strip()
                    if not c: continue
                    
                    # 处理 "var:type" 语法
                    if ":" in c:
                        parts = c.split(":")
                        if len(parts) == 2:
                            var_name = parts[0].strip()
                            type_name = parts[1].strip().lower()
                            
                            # 忽略 redundant real declaration
                            if type_name == "real":
                                continue
                            # 转换 int declaration
                            elif type_name == "int" or type_name == "integer":
                                c = f"isInt({var_name})"
                    
                    # 尝试解析
                    try:
                        cond_list.append(parser.parse_expr(c))
                    except Exception as e:
                        # 记录错误但不崩溃? 或者抛出更友好的错误
                        raise type(e)(f"无法解析条件 '{c}': {str(e)}")

            # 判断命令类型并创建对应的Action
            expr_lower = expression.lstrip().lower()
            
            if expr_lower.startswith("prove "):
                # prove 命令：使用parser.parse_action解析整个命令
                # 这会创建ProveAction并进入PROVE状态
                action = parser.parse_action(expression)
                self.state = initial_state.process_action(action)
            elif expr_lower.startswith("calculate "):
                # calculate 命令：去掉前缀，使用CalculateAction
                clean_expr = expression.lstrip()[len("calculate "):].strip()
                parsed_expr = parser.parse_expr(clean_expr)
                calc_action = CalculateAction(parsed_expr, cond_list)
                self.state = initial_state.process_action(calc_action)
            else:
                # 没有前缀：默认为calculate
                parsed_expr = parser.parse_expr(expression)
                calc_action = CalculateAction(parsed_expr, cond_list)
                self.state = initial_state.process_action(calc_action)
            
            self.initial_expression = expression
            current_expr = self._get_current_expr()
            self._previous_expr_str = str(current_expr) if current_expr else expression
            self.history = []
            
            return ExecutionResult(
                success=True,
                result=self._previous_expr_str,
                latex_result=self._expr_to_latex(current_expr),
                changed=False,
                raw_expr=current_expr
            )
            
        except Exception as e:
            error_type = ErrorType.SYNTAX if "parse" in str(type(e).__name__).lower() else ErrorType.SYSTEM
            return ExecutionResult(
                success=False,
                error=f"初始化错误: {str(e)}",
                error_type=error_type
            )
    
    def execute(self, command: str) -> ExecutionResult:
        """执行单个命令"""
        if self.state is None:
            return ExecutionResult(
                success=False,
                error="状态未初始化，请先调用initialize()",
                error_type=ErrorType.SYSTEM
            )
        
        expr_before = self._previous_expr_str or ""
        
        # 记录执行前的步骤数
        steps_before_count = 0
        if hasattr(self.state, 'calc') and hasattr(self.state.calc, 'steps'):
            steps_before_count = len(self.state.calc.steps)
        
        try:
            # 解析命令为Action
            action = parser.parse_action(command)
            
            # 执行命令
            self.state = self.state.process_action(action)
            
            # 获取执行后的表达式
            current_expr = self._get_current_expr()
            expr_after = str(current_expr) if current_expr else expr_before
            changed = expr_before != expr_after
            
            # 获取新增的中间步骤
            intermediate_steps = []
            if hasattr(self.state, 'calc') and hasattr(self.state.calc, 'steps'):
                current_steps = self.state.calc.steps
                if len(current_steps) > steps_before_count:
                    new_steps = current_steps[steps_before_count:]
                    for step in new_steps:
                        # 提取规则名称
                        rule_str = str(step.rule)
                        # 尝试更友好的显示
                        if hasattr(step.rule, 'name'):
                             rule_str = step.rule.name
                        
                        intermediate_steps.append({
                            'rule': str(step.rule), # 保持完整字符串形式
                            'res': str(step.res),
                            'latex': self._expr_to_latex(step.res)
                        })

            self.history.append({
                "command": command,
                "expr_before": expr_before,
                "expr_after": expr_after,
                "success": True,
                "changed": changed,
                "intermediate_steps": intermediate_steps
            })
            
            self._previous_expr_str = expr_after
            
            return ExecutionResult(
                success=True,
                result=expr_after,
                latex_result=self._expr_to_latex(current_expr),
                changed=changed,
                raw_expr=current_expr,
                intermediate_steps=intermediate_steps
            )
            
        except Exception as e:
            error_name = type(e).__name__.lower()
            if "parse" in error_name:
                error_type = ErrorType.SYNTAX
                error_msg = f"命令语法错误: {str(e)}"
            elif "rule" in error_name or "assertion" in error_name:
                error_type = ErrorType.RULE
                error_msg = f"规则应用错误: {str(e)}"
            else:
                error_type = ErrorType.SYSTEM
                error_msg = f"执行错误: {str(e)}"
            
            self.history.append({
                "command": command,
                "expr_before": expr_before,
                "success": False,
                "error": error_msg
            })
            
            return ExecutionResult(
                success=False,
                error=error_msg,
                error_type=error_type
            )
    
    def _get_current_expr(self) -> Optional[Any]:
        """获取当前表达式对象
        
        适配不同状态：
        - CalculateState: 有calc属性
        - ProveState: 返回goal.goal（实际表达式），而不是Goal对象
        - 其他状态: 尝试访问相关属性
        """
        if self.state is None:
            return None
        
        # CalculateState有calc属性
        if hasattr(self.state, 'calc'):
            calc = self.state.calc
            if calc.steps:
                return calc.steps[-1].res
            return calc.start
        
        # ProveState有goal属性（Goal对象）
        # 返回goal.goal（实际的Expr）而不是Goal对象
        # 避免LLM看到"Goal (finished)"这种误导信息
        if hasattr(self.state, 'goal'):
            goal_obj = self.state.goal
            if hasattr(goal_obj, 'goal'):
                return goal_obj.goal  # 返回实际表达式
            return goal_obj
        
        return None
    
    def get_current_expression(self) -> Optional[str]:
        """获取当前表达式字符串"""
        expr = self._get_current_expr()
        return str(expr) if expr else None
    
    def get_current_latex(self) -> Optional[str]:
        """获取当前表达式的LaTeX"""
        expr = self._get_current_expr()
        return self._expr_to_latex(expr) if expr else None
    
    def is_finished(self) -> bool:
        """检查是否已完成化简
        
        系统性修复：原始库的StateItem.is_finished()默认返回True，
        但Goal类没有覆盖此方法，导致没有proof时也返回True。
        这里在llm_iscalc层面做正确的判断。
        """
        if self.state is None:
            return False
        
        # 获取状态类型
        state_type = type(self.state).__name__
        
        # ProveState: 检查goal是否真正完成
        if state_type == "ProveState":
            goal = getattr(self.state, 'goal', None)
            if goal is not None:
                # Goal只有在有proof且proof完成时才算完成
                proof = getattr(goal, 'proof', None)
                if proof is None:
                    return False
                # 有proof时，再检查proof是否完成
                return proof.is_finished() if hasattr(proof, 'is_finished') else False
            return False
        
        # CalculateState: 计算状态总是被认为未完成（正在进行中）
        # 用户必须显式使用done命令退出此状态
        if state_type == "CalculateState":
            return False
            
        # InitialState: 初始状态被认为是"完成"（即处于空闲/就绪状态）
        # 当从最顶层状态执行done返回到InitialState时，也表示任务完成
        if state_type == "InitialState":
            return True
        
        # InductionState: 检查所有分支
        if state_type == "InductionState":
            induct_proof = getattr(self.state, 'induct_proof', None)
            if induct_proof is not None:
                return induct_proof.is_finished()
            return False
        
        # CaseAnalysisState: 检查所有case
        if state_type == "CaseAnalysisState":
            case_proof = getattr(self.state, 'case_proof', None)
            if case_proof is not None:
                return case_proof.is_finished()
            return False
        
        # 其他状态：回退到原始逻辑
        return self.state.is_finished()
    
    def get_history(self) -> List[dict]:
        return self.history.copy()
    
    def _expr_to_latex(self, e: Any) -> str:
        """将表达式转换为LaTeX"""
        if e is None:
            return ""
        try:
            if hasattr(e, 'export_latex'):
                return e.export_latex()
            elif hasattr(e, 'latex'):
                return e.latex()
            else:
                return str(e)
        except Exception:
            return str(e)
    
    def get_current_state_name(self) -> str:
        """获取当前状态名称
        
        返回状态的类名简称，用于技能匹配和提示词构建。
        可能的返回值: 'INITIAL', 'CALCULATE', 'PROVE', 'INDUCTION', 'CASE' 等
        """
        if self.state is None:
            return "INITIAL"
        
        state_class = type(self.state).__name__
        # 标准化状态名称：去掉 'State' 后缀，转大写
        if state_class.endswith("State"):
            state_class = state_class[:-5]
        return state_class.upper()
    
    def reset(self) -> None:
        """重置执行器状态"""
        self.comp_file = None
        self.state = None
        self.history = []
        self.initial_expression = None
        self._previous_expr_str = None
