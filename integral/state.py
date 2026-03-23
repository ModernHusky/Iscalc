"""State machine for processing the actions."""

from typing import Optional
import time

from integral import expr
from integral.rules import IntegrateByEquation, RuleException
from integral import compstate
from integral.compstate import Calculation, Goal, StateException
from integral.context import Context
from integral import poly
from integral.action import Action, CalculateAction, ProveAction, LHSAction, \
    RHSAction, DefineAction, ArgAction, RewriteGoalAction, InductionAction, \
    CaseAnalysisAction, SubgoalAction, DoneAction, RuleAction, SorryAction, \
    BaseCaseAction, InductCaseAction, CaseAction, ImportsAction, LetAction


def extract_solving_skolem_funcs(goal_expr: expr.Expr) -> list[str]:
    """Extract Skolem functions that are being solved for in the goal.

    For example, if the goal is: SKOLEM_FUNC(C(a)) = pi/2
    This returns ["C(a)"] to mark it as being solved.
    """
    solving_funcs = []

    # Check if it's an equation with Skolem func on LHS
    if goal_expr.is_equals():
        lhs = goal_expr.lhs
        if expr.is_skolem_func(lhs):
            if len(lhs.dependent_vars) == 0:
                solving_funcs.append(lhs.name)
            else:
                args_str = ",".join(str(arg) for arg in lhs.dependent_vars)
                solving_funcs.append(f"{lhs.name}({args_str})")

    return solving_funcs


class State:
    """Base class for states."""

    """The previous state this state comes from."""
    past: Optional["State"]

    def process_action(self, action: Action) -> "State":
        """Apply the given action, return new state."""
        pass

    def is_finished(self) -> bool:
        """Determine whether the given state is in finished form."""
        pass


class InitialState(State):
    """Initial state."""
    def __init__(self, ctx: Context, proof_finished: bool = False):
        self.ctx = ctx
        self.past = None
        self._proof_finished = proof_finished

    def process_action(self, action: Action) -> State:
        # Start a calculation
        if isinstance(action, CalculateAction):
            calc = Calculation(None, self.ctx, action.expr, conds=action.conditions)
            return CalculateState(self, calc)
        
        # Start a proof
        elif isinstance(action, ProveAction):
            if expr.is_equals(action.expr) and expr.is_indefinite_integral(action.expr.lhs):
                action.conditions.add_condition(expr.isReal(expr.Var(action.expr.lhs.var)))
            goal = Goal(None, self.ctx, action.expr, conds=action.conditions)
            # Auto-mark Skolem functions that are being solved for
            for skolem_str in extract_solving_skolem_funcs(action.expr):
                goal.ctx.mark_skolem_solving(skolem_str)
            return ProveState(self, goal)
        
        # Add a definition
        elif isinstance(action, DefineAction):
            self.ctx.add_definition(action.expr, conds=action.conditions)
            return self

        # Importing a theory, ignored for now        
        elif isinstance(action, ImportsAction):
            return self
        
        # Other actions are invalid
        else:
            raise StateException(
                "Initial",
                f"Action type {type(action).__name__} cannot be performed in initial state")
        
    def is_finished(self) -> bool:
        return self._proof_finished

    def __str__(self):
        return "(initial)"

class ProveState(State):
    """State when performing a proof."""
    def __init__(self, past: State, goal: Goal):
        self.past = past
        assert goal is not None
        self.goal = goal

    def process_action(self, action: Action) -> State:
        # Prove by calculating on both sides
        if isinstance(action, LHSAction):
            if not self.goal.proof:
                self.goal.proof_by_calculation()
            if not isinstance(self.goal.proof, compstate.CalculationProof):
                raise StateException("Prove", "lhs: can not appear in a calculation proof")
            return CalculateState(self, self.goal.proof.lhs_calc)
        
        elif isinstance(action, RHSAction):
            if not self.goal.proof:
                self.goal.proof_by_calculation()
            if not isinstance(self.goal.proof, compstate.CalculationProof):
                raise StateException("Prove", "rhs: not in calculation proof")
            return CalculateState(self, self.goal.proof.rhs_calc)
        
        elif isinstance(action, ArgAction):
            if self.goal.goal.is_equals():
                raise StateException(
                    "Prove",
                    f"Action type {type(action).__name__} cannot be performed in prove state because {str(self.goal.goal)} is a equation. If you wish to use calculation proof, use lhs: or rhs: to enter the calculate state.")
            if not self.goal.proof:
                self.goal.proof_by_calculation()
            if not isinstance(self.goal.proof, compstate.CalculationProof):
                raise StateException("Prove", "arg: not in calculation proof")
            return CalculateState(self, self.goal.proof.arg_calc)

        # Prove by rewriting goal
        elif isinstance(action, RewriteGoalAction):
            proof = self.goal.proof_by_rewrite_goal(begin=action.name)
            return CalculateState(self, proof.begin)
        
        # Prove by induction
        elif isinstance(action, InductionAction):
            proof = self.goal.proof_by_induction(action.var_name, start=action.start)
            return InductionState(self, proof)
        
        # Prove by case analysis
        elif isinstance(action, CaseAnalysisAction):
            proof = self.goal.proof_by_case(action.split_cond)
            return CaseAnalysisState(self, proof)

        # Start a subgoal
        elif isinstance(action, SubgoalAction):
            if self.goal.parent is not None and isinstance(self.goal.parent, Goal):
                raise StateException(
                "Prove",
                f"Sub-goals cannot be nested.")
            subgoal = self.goal.add_subgoal(action.name, action.expr, action.conditions)
            # Auto-mark Skolem functions that are being solved for
            for skolem_str in extract_solving_skolem_funcs(action.expr):
                subgoal.ctx.mark_skolem_solving(skolem_str)
            return ProveState(self, subgoal)
        
        # Done with current subgoal
        elif isinstance(action, DoneAction):
            self.goal.check_finished(stack=tuple())
            # Return to past, preserving proof completion status
            if isinstance(self.past, InitialState):
                return InitialState(self.past.ctx, proof_finished=True)
            return self.past

        # Make local definition
        elif isinstance(action, LetAction):
            self.goal.add_definition(action.expr, action.conditions)
            return self

        # Make definition
        elif isinstance(action, DefineAction):
            self.goal.add_definition(action.expr, action.conditions)
            return self

        # Abandon the current calculation or proof
        elif isinstance(action, SorryAction):
            if isinstance(self.past, InitialState):
                return self.past
            else:
                return self.past.process_action(action)

        # Other cases are invalid
        else:
            raise StateException(
                "Prove",
                f"Action type {type(action).__name__} cannot be performed in prove state")

    def is_finished(self) -> bool:
        return self.goal.is_finished()

    def __str__(self):
        return "(proof)\n%s" % self.goal


class CalculateState(State):
    """State when performing a calculation."""
    def __init__(self, past: State, calc: Calculation):
        self.past = past
        self.calc = calc

    def process_action(self, action: Action) -> State:
        # Perform a rule
        if isinstance(action, RuleAction):
            # Special check for IntegrateByEquation: lhs must appear exactly
            # as one of the steps.
            if isinstance(action.rule, IntegrateByEquation):
                # Check if lhs appears as a subexpression in any step
                def check_step(expr):
                    # Use the existing find_subexpr_pred method to check if lhs is a subexpression
                    subexprs = expr.find_subexpr_pred(lambda e: e == action.rule.lhs, is_nested=False)
                    return len(subexprs) > 0
                
                if not (check_step(self.calc.start) or
                        any(check_step(step.res) for step in self.calc.steps)):
                    # print("Current calculation is:")
                    # print(self.calc)
                    raise RuleException(
                        "IntegrateByEquation",
                        f"lhs {action.rule.lhs} must appear as one of the steps")
            self.calc.perform_rule(action.rule)
            return self
        
        # Handle rewrite goal action
        elif isinstance(action, RewriteGoalAction):
            if isinstance(self.past, ProveState):
                proof = self.past.goal.proof_by_rewrite_goal(begin=action.name)
                return CalculateState(self, proof.begin)
            else:
                raise StateException(
                    "Calculate",
                    "RewriteGoalAction can only be performed when past state is ProveState")
        
        # Done with current calculation or proof
        elif isinstance(action, DoneAction):
            if not self.is_finished():
                msg = "Use done when calculation is not finished\n"
                msg += f"Final expression {self.calc.steps[-1].res} is not closed"
                raise StateException("Done", msg)
            if isinstance(self.past, InitialState):
                return InitialState(self.past.ctx, proof_finished=True)
            else:
                return self.past.process_action(action)
            
        # Abandon the current calculation or proof
        elif isinstance(action, SorryAction):
            if isinstance(self.past, InitialState):
                return self.past
            else:
                return self.past.process_action(action)
        
        # Go to the other branch
        elif isinstance(action, RHSAction):
            return self.past.process_action(action)

        elif isinstance(action, SubgoalAction):
            subgoal = self.calc.add_subgoal(action.name, action.expr, action.conditions)
            return ProveState(self, subgoal)

        # Other cases are invalid
        else:
            raise StateException(
                "Calculate",
                f"Action type {type(action).__name__} cannot be performed in calculate state")

    def is_finished(self) -> bool:
        if isinstance(self.past, InitialState):
            if not self.calc.steps:
                return False

            res = self.calc.steps[-1].res
            ctx = self.calc.ctx
            cur_e = self.calc.start
            for step in self.calc.steps:
                ctx = step.rule.update_context(cur_e, ctx)
                cur_e = step.res
            substs = ctx.get_substs()
            for var, _, _ in substs:
                if res.contains_var(var) and not self.calc.start.contains_var(var):
                    return False

            self.calc.check_wellformed()

            return res.is_closed_form() and poly.normalize(res, self.calc.ctx) == res
        else:
            return self.past.is_finished()

    def __str__(self):
        return "(calculate)\n%s" % self.calc


class InductionState(State):
    """State when performing an induction."""
    def __init__(self, past: State, induct_proof: compstate.InductionProof):
        self.past = past
        self.induct_proof = induct_proof
    
    def process_action(self, action: Action) -> State:
        if isinstance(action, BaseCaseAction):
            return ProveState(self, self.induct_proof.base_case)
        elif isinstance(action, InductCaseAction):
            return ProveState(self, self.induct_proof.induct_case)
        elif isinstance(action, DoneAction):
            return self.past.process_action(action)
        else:
            raise StateException(
                "Induction",
                f"Action type {type(action).__name__} cannot be performed in induction state")
    
    def is_finished(self) -> bool:
        return self.induct_proof.is_finished()


class CaseAnalysisState(State):
    """State when performing case analysis."""
    def __init__(self, past: State, case_proof: compstate.CaseProof):
        self.past = past
        self.case_proof = case_proof

    def __str__(self):
        return str(self.case_proof)

    def process_action(self, action: Action) -> State:
        if isinstance(action, CaseAction):
            if action.mark == "true":
                if self.case_proof.split_type != "two-way":
                    raise StateException("Case", "case true when not in two-way analysis.")
                return ProveState(self, self.case_proof.cases[0])
            elif action.mark == "false":
                if self.case_proof.split_type != "two-way":
                    raise StateException("Case", "case false when not in two-way analysis.")
                return ProveState(self, self.case_proof.cases[1])
            elif action.mark == "negative":
                if self.case_proof.split_type != "three-way":
                    raise StateException("Case", "case negative when not in three-way analysis.")
                return ProveState(self, self.case_proof.cases[0])
            elif action.mark == "zero":
                if self.case_proof.split_type != "three-way":
                    raise StateException("Case", "case zero when not in three-way analysis.")
                return ProveState(self, self.case_proof.cases[1])
            elif action.mark == "positive":
                if self.case_proof.split_type != "three-way":
                    raise StateException("Case", "case positive when not in three-way analysis.")
                return ProveState(self, self.case_proof.cases[2])
            else:
                raise StateException("Case", f"Unknown case {action.mark}")
        elif isinstance(action, DoneAction):
            return self.past.process_action(action)
        else:
            raise StateException(
                "Case",
                f"Action type {type(action).__name__} cannot be performed in case state")
    
    def is_finished(self) -> bool:
        return self.case_proof.is_finished()



class ProblemInfo:
    """Information about a single problem.
    
    Attributes
    ----------
    filename : str
        name of the file, used for output
    index : int
        index of the problem, 1-based
    context : Context
        context for the problem
    problem : str
        statement of the problem
    steps : list[str]
        available answer for the problem

    """
    def __init__(self, filename: str, index: int, ctx: Context, problem: str, steps: list[str]):
        self.filename = filename
        self.index = index
        self.context = Context(ctx)
        self.problem = problem
        self.steps = steps

    def __str__(self):
        res = f"{self.filename}_{self.index}"
        if self.steps == ["sorry"]:
            res += "?"
        res += " " + self.problem
        return res


def process_file(filename: str) -> list[ProblemInfo]:
    """Process the content of a file containing calculations.
    
    Parameters
    ----------
    filename : str
        name of the file, in `iscalc/theories`.

    Returns
    -------
    list[ProblemInfo]
        list of problem infos contained in the file.

    """
    from integral import parser

    result = []
    with open(f'../iscalc/theories/{filename}.thy', 'r', encoding="utf-8") as problem_file:
        lines = [s for s in problem_file.read().split('\n') if s.strip()]
        cur_goal = None
        steps = []
        i = 0
        ctx = Context()
        ctx.load_book("base")
        for line in lines:
            line = line.strip()
            if line.startswith('#') or line.startswith('//'):
                # title or comment
                continue
            a = parser.parse_action(line)
            if isinstance(a, ImportsAction):
                for theory in a.theories:
                    if theory != 'base':
                        ctx.load_book(theory)
            if isinstance(a, (ProveAction, CalculateAction)):
                if cur_goal:
                    # First create problem using context *without* adding the current
                    # goal as theorem.
                    result.append(ProblemInfo(filename, i, ctx, problem, steps))
                    ctx = Context(ctx)

                    # Then add current theorem to context.
                    if isinstance(cur_goal, ProveAction):
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
        result.append(ProblemInfo(filename, i, ctx, problem, steps))
    return result


def check_actions(content: str, *, print_lines=False, print_state=False,
                  write_stats=False, filename=""):
    from integral import parser

    actions = content.split('\n')
    ctx = Context()
    st = InitialState(ctx)
    start_time = None
    cur_goal = None
    for i, act in enumerate(actions, 1):
        if print_lines:
            print(act)
        if not act.strip():
            # empty line
            continue
        if act.lstrip().startswith('#') or act.lstrip().startswith('//'):
            # title or comment
            continue
        a = parser.parse_action(act)
        if isinstance(a, ImportsAction):
            for thy_name in a.theories:
                ctx.load_book(thy_name)
        if isinstance(a, (ProveAction, CalculateAction)):
            cur_goal = a
            if write_stats:
                start_time = time.time()
                with open("stats.txt", "a", encoding='utf-8') as stats_file:
                    stats_file.write(f"{filename} {i} {cur_goal}\n")
        try:
            st = st.process_action(a)
            if isinstance(st, InitialState):
                if isinstance(cur_goal, ProveAction):
                    if cur_goal.expr.is_equals() and expr.is_indefinite_integral(cur_goal.expr.lhs):
                        ctx.add_indefinite_integral(cur_goal.expr, cur_goal.conditions, cur_goal.attrs)
                    elif cur_goal.expr.is_equals() and expr.is_integral(cur_goal.expr.lhs):
                        ctx.add_definite_integral(cur_goal.expr, cur_goal.conditions, cur_goal.attrs)
                    else:
                        ctx.add_other_identities(cur_goal.expr, cur_goal.conditions, cur_goal.attrs)
                if cur_goal and write_stats:
                    elapsed_time = time.time() - start_time
                    with open("stats.txt", "a", encoding='utf-8') as stats_file:
                        stats_file.write(f"{elapsed_time:.2f} seconds\n")
                cur_goal = None
        except Exception as e:
            print(cur_goal)
            print(st)
            raise e
    if print_state:
        print(st)
    if not print_state and not isinstance(st, InitialState):
        raise AssertionError("Does not end in initial state (add print_state=True to debug)")
