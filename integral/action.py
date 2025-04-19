"""Definition of internal language for actions."""

from typing import Optional

from integral import expr
from integral.expr import Expr
from integral.rules import Rule, IntegrateByEquation, RuleException
from integral import compstate
from integral.compstate import Calculation, Goal, CompFile, StateException
from integral.conditions import Conditions
from integral import poly

class Action:
    """Base class for actions."""
    def get_start_states(self) -> list[str]:
        raise NotImplementedError(f"get_start_states: {type(self)}")


class ProveAction(Action):
    """Start a proof."""
    def __init__(self, expr: Expr, conditions: Optional[Conditions] = None):
        self.expr = expr
        self.conditions = Conditions(conditions)

    def __str__(self):
        if self.conditions:
            return "prove %s for %s" % (self.expr, ', '.join(str(cond) for cond in self.conditions.data))
        else:
            return "prove %s" % self.expr

    def get_start_states(self) -> list[str]:
        return ["initial", "proof"]


class DefineAction(Action):
    """Make a definition."""
    def __init__(self, expr: Expr, conditions: Optional[Conditions] = None):
        self.expr = expr
        self.conditions = Conditions(conditions)

    def __str__(self):
        if self.conditions:
            return "define %s for %s" % (self.expr, ', '.join(str(cond) for cond in self.conditions.data))
        else:
            return "define %s" % self.expr

    def get_start_states(self) -> list[str]:
        return ["initial", "proof"]


class SubgoalAction(Action):
    """Start proof of a subgoal."""
    def __init__(self, name: str, expr: Expr, conditions: Optional[Conditions] = None):
        self.name = name
        self.expr = expr
        self.conditions = Conditions(conditions)

    def __str__(self):
        if self.conditions:
            return "subgoal %s: %s for %s" % (
                self.name, self.expr, ', '.join(str(cond) for cond in self.conditions))
        else:
            return "subgoal %s: %s" % (self.name, self.expr)

    def get_start_states(self) -> list[str]:
        return ["proof"]


class DoneAction(Action):
    """Done with current goal or subgoal."""
    def __init__(self):
        pass

    def __str__(self):
        return "done"

    def get_start_states(self) -> list[str]:
        return ["proof", "calculation", "case", "induction"]


class SorryAction(Action):
    """Abandon the current computation."""
    def __init__(self):
        pass

    def __str__(self):
        return "sorry"

    def get_start_states(self) -> list[str]:
        return ["proof"]


class RewriteGoalAction(Action):
    """Invoke rewriting goal."""
    def __init__(self, name: str):
        self.name = name

    def __str__(self):
        return "from %s:" % self.name

    def get_start_states(self) -> list[str]:
        return ["proof"]


class CalculateAction(Action):
    """Start a calculation."""
    def __init__(self, expr: Expr, conditions: Optional[Conditions] = None):
        self.expr = expr
        self.conditions = Conditions(conditions)

    def __str__(self):
        if self.conditions.data:
            return "calculate %s for %s" % (
                self.expr, ', '.join(str(cond) for cond in self.conditions.data))
        else:
            return "calculate %s" % self.expr

    def get_start_states(self) -> list[str]:
        return ["initial"]


class InductionAction(Action):
    """Start an induction."""
    def __init__(self, var_name: str, start: Expr):
        self.var_name = var_name
        self.start = start

    def __str__(self):
        if self.start == expr.Const(0):
            return f"induction on {self.var_name}"
        else:
            return f"induction on {self.var_name} starting from {self.start}"

    def get_start_states(self) -> list[str]:
        return ["proof"]


class CaseAnalysisAction(Action):
    """Start a case analysis."""
    def __init__(self, split_cond: Expr):
        self.split_cond = split_cond

    def __str__(self):
        return "case analysis on %s" % self.split_cond

    def get_start_states(self) -> list[str]:
        return ["proof"]


class LHSAction(Action):
    """Perform a proof by working on the left hand side."""
    def __init__(self):
        pass

    def __str__(self):
        return "lhs:"

    def get_start_states(self) -> list[str]:
        return ["proof"]


class RHSAction(Action):
    """Perform a proof by working on the right hand side."""
    def __init__(self):
        pass

    def __str__(self):
        return "rhs:"

    def get_start_states(self) -> list[str]:
        return ["proof"]


class ArgAction(Action):
    """Perform a proof by working on the argument."""
    def __init__(self):
        pass

    def __str__(self):
        return "arg:"

    def get_start_states(self) -> list[str]:
        return ["proof"]


class BaseCaseAction(Action):
    """Base case of an induction."""
    def __init__(self):
        pass

    def __str__(self):
        return "base:"
    
    def get_start_states(self) -> list[str]:
        return ["induction"]


class InductCaseAction(Action):
    """Induct case of an induction."""
    def __init__(self):
        pass

    def __str__(self):
        return "induct:"

    def get_start_states(self) -> list[str]:
        return ["induction"]


class CaseAction(Action):
    """Case in case analysis."""
    def __init__(self, mark: str):
        self.mark = mark

    def __str__(self):
        return "case %s" % self.mark

    def get_start_states(self) -> list[str]:
        return ["case"]


class RuleAction(Action):
    """Apply rule."""
    def __init__(self, rule: Rule):
        self.rule = rule

    def __str__(self):
        return str(self.rule)
    
    def get_start_states(self) -> list[str]:
        return ["calculation"]


"""State machine for processing the actions."""


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
    def __init__(self, comp_file: CompFile):
        self.comp_file = comp_file
        self.past = None

    def process_action(self, action: Action) -> State:
        # Start a calculation
        if isinstance(action, CalculateAction):
            calc = Calculation(self.comp_file, self.comp_file.ctx, action.expr, conds=action.conditions)
            return CalculateState(self, calc)
        
        # Start a proof
        elif isinstance(action, ProveAction):
            goal = compstate.Goal(self.comp_file, self.comp_file.ctx, action.expr, conds=action.conditions)
            return ProveState(self, goal)
        
        # Add a definition
        elif isinstance(action, DefineAction):
            self.comp_file.ctx.add_definition(action.expr, conds=action.conditions)
            return self
        
        # Other actions are invalid
        else:
            raise StateException(
                "Initial",
                f"Action type {type(action).__name__} cannot be performed in initial state")
        
    def is_finished(self) -> bool:
        return False

    def __str__(self):
        return "(initial)"

class ProveState(State):
    """State when performing a proof."""
    def __init__(self, past: State, goal: Goal):
        self.past = past
        self.goal = goal

    def process_action(self, action: Action) -> State:
        # Prove by calculating on both sides
        if isinstance(action, LHSAction):
            if not self.goal.proof:
                self.goal.proof_by_calculation()
            if not isinstance(self.goal.proof, compstate.CalculationProof):
                raise StateException("Prove", "lhs: not in calculation proof")
            return CalculateState(self, self.goal.proof.lhs_calc)
        
        elif isinstance(action, RHSAction):
            if not self.goal.proof:
                self.goal.proof_by_calculation()
            if not isinstance(self.goal.proof, compstate.CalculationProof):
                raise StateException("Prove", "rhs: not in calculation proof")
            return CalculateState(self, self.goal.proof.rhs_calc)
        
        elif isinstance(action, ArgAction):
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
            subgoal = self.goal.add_subgoal(action.name, action.expr, action.conditions)
            return ProveState(self, subgoal)
        
        # Done with current subgoal
        elif isinstance(action, DoneAction):
            self.goal.check_finished(stack=tuple())
            if isinstance(self.past, InitialState):
                if self.goal.goal.is_equals() and expr.is_integral(self.goal.goal.lhs):
                    self.past.comp_file.ctx.add_definite_integral(self.goal.goal, self.goal.conds)
                elif self.goal.goal.is_equals() and expr.is_indefinite_integral(self.goal.goal.lhs):
                    self.past.comp_file.ctx.add_indefinite_integral(self.goal.goal)
                else:
                    self.past.comp_file.ctx.add_lemma(self.goal.goal, self.goal.conds)
            return self.past

        # Make definition
        elif isinstance(action, DefineAction):
            self.goal.add_definition(action.expr, action.conditions)
            return self
        
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
                if not (self.calc.start == action.rule.lhs or
                        any(step.res == action.rule.lhs for step in self.calc.steps)):
                    # print("Current calculation is:")
                    # print(self.calc)
                    raise RuleException(
                        "IntegrateByEquation",
                        f"lhs {action.rule.lhs} must appear exactly as one of the steps")
            self.calc.perform_rule(action.rule)
            return self
        
        # Done with current calculation or proof
        elif isinstance(action, DoneAction):
            if isinstance(self.past, InitialState):
                if not self.is_finished():
                    msg = "Use done when calculation is not finished\n"
                    msg += f"Final expression {self.calc.steps[-1].res} is not closed"
                    raise StateException("Done", msg)
                return self.past
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
            for var, _ in substs:
                if res.contains_var(var) and not self.calc.start.contains_var(var):
                    return False
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
