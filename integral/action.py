"""Definition of internal language for actions."""

from typing import Optional, Tuple

from integral import expr
from integral.expr import Expr
from integral.rules import Rule
from integral.context import Context
from integral import compstate
from integral.compstate import Calculation, StateItem, Goal, CompFile
from integral.conditions import Conditions
from integral import poly

class Action:
    """Base class for actions."""
    pass


class ProveAction(Action):
    """Start a proof."""
    def __init__(self, expr: Expr, conditions: Optional[Conditions] = None):
        self.expr = expr
        self.conditions = Conditions(conditions)

    def __str__(self):
        if self.conditions:
            return "prove %s for %s" % (self.expr, ', '.join(str(cond) for cond in self.conditions))
        else:
            return "prove %s" % self.expr


class DefineAction(Action):
    """Make a definition."""
    def __init__(self, expr: Expr, conditions: Optional[Conditions] = None):
        self.expr = expr
        self.conditions = Conditions(conditions)

    def __str__(self):
        if self.conditions:
            return "define %s for %s" % (self.expr, ', '.join(str(cond) for cond in self.conditions))
        else:
            return "define %s" % self.expr


class SubgoalAction(Action):
    """Start a proof."""
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


class DoneAction(Action):
    """Done with current subgoal."""
    def __init__(self):
        pass

    def __str__(self):
        return "done"


class RewriteGoalAction(Action):
    """Invoke rewriting goal."""
    def __init__(self, name: str):
        self.name = name

    def __str__(self):
        return "from %s:" % self.name


class CalculateAction(Action):
    """Start a calculation."""
    def __init__(self, expr: Expr):
        self.expr = expr

    def __str__(self):
        return "calculate %s" % self.expr
    

class InductionAction(Action):
    """Start an induction."""
    def __init__(self, var_name: str, start: Expr):
        self.var_name = var_name
        self.start = start

    def __str__(self):
        if self.start == expr.Const(0):
            return "induction on %s" % self.var_name
        else:
            return "induction on %s starting from %s" % (self.var_name, self.expr)


class LHSAction(Action):
    """Perform a proof by working on the left hand side."""
    def __init__(self):
        pass

    def __str__(self):
        return "lhs:"

class RHSAction(Action):
    """Perform a proof by working on the right hand side."""
    def __init__(self):
        pass

    def __str__(self):
        return "rhs:"

class ArgAction(Action):
    """Perform a proof by working on the argument."""
    def __init__(self):
        pass

    def __str__(self):
        return "arg:"

class BaseCaseAction(Action):
    """Base case of an induction."""
    def __init__(self):
        pass

    def __str__(self):
        return "base:"
    
class InductCaseAction(Action):
    """Induct case of an induction."""
    def __init__(self):
        pass

    def __str__(self):
        return "induct:"

class RuleAction(Action):
    """Apply rule."""
    def __init__(self, rule: Rule):
        self.rule = rule

    def __str__(self):
        return str(self.rule)
    

"""State machine for processing the actions."""


class StateException(Exception):
    """Exception resulting from applying action to a state."""
    def __init__(self, msg: str):
        self.msg = msg

    def __str__(self):
        return self.msg


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
            calc = Calculation(self.comp_file, self.comp_file.ctx, action.expr)
            return CalculateState(self, calc)
        
        # Start a proof
        elif isinstance(action, ProveAction):
            goal = compstate.Goal(self.comp_file, self.comp_file.ctx, action.expr,
                                  conds=action.conditions)
            return ProveState(self, goal)
        
        # Add a definition
        elif isinstance(action, DefineAction):
            self.comp_file.add_definition(action.expr, conds=action.conditions)
            return self
        
        # Other actions are invalid
        elif isinstance(action, RuleAction):
            raise StateException("Cannot apply rule when at initial state.")
        else:
            raise StateException("Unknown action type %s" % type(action))
        
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
                raise StateException("lhs: not in calculation proof")
            return CalculateState(self, self.goal.proof.lhs_calc)
        
        elif isinstance(action, RHSAction):
            if not self.goal.proof:
                self.goal.proof_by_calculation()
            if not isinstance(self.goal.proof, compstate.CalculationProof):
                raise StateException("rhs: not in calculation proof")
            return CalculateState(self, self.goal.proof.rhs_calc)
        
        elif isinstance(action, ArgAction):
            if not self.goal.proof:
                self.goal.proof_by_calculation()
            if not isinstance(self.goal.proof, compstate.CalculationProof):
                raise StateException("arg: not in calculation proof")
            return CalculateState(self, self.goal.proof.arg_calc)

        # Prove by rewriting goal
        elif isinstance(action, RewriteGoalAction):
            proof = self.goal.proof_by_rewrite_goal(begin=action.name)
            return CalculateState(self, proof.begin)
        
        # Prove by induction
        elif isinstance(action, InductionAction):
            proof = self.goal.proof_by_induction(action.var_name, start=action.start)
            return InductionState(self, proof)

        # Start a subgoal
        elif isinstance(action, SubgoalAction):
            subgoal = self.goal.add_subgoal(action.name, action.expr, action.conditions)
            return ProveState(self, subgoal)
        
        # Done with current subgoal
        elif isinstance(action, DoneAction):
            if isinstance(self.past, InitialState):
                if self.goal.goal.is_equals() and expr.is_integral(self.goal.goal.lhs):
                    self.past.comp_file.ctx.add_definite_integral(self.goal.goal, self.goal.conds)
            return self.past

        # Make definition
        elif isinstance(action, DefineAction):
            self.goal.add_definition(action.expr, action.conditions)
            return self
        
        # Other cases are invalid
        else:
            raise StateException("Unknown action type %s" % type(action))

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
            self.calc.perform_rule(action.rule)
            return self
        
        # Done with current calculation or proof
        elif isinstance(action, DoneAction):
            return self.past.process_action(action)
        
        # Go to the other branch
        elif isinstance(action, RHSAction):
            return self.past.process_action(action)
        
        # Other cases are invalid
        elif isinstance(action, CalculateAction):
            raise StateException("Cannot start a new calculation within a calculation.")
        else:
            raise StateException("Unknown action type %s" % type(action))

    def is_finished(self) -> bool:
        if isinstance(self.past, InitialState):
            if not self.calc.steps:
                return False

            res = self.calc.steps[-1].res
            return res.is_evaluable() and poly.normalize(res, self.calc.ctx) == res
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
            raise StateException("Unknown action type %s" % type(action))
    
    def is_finished(self) -> bool:
        return self.induct_proof.is_finished()
