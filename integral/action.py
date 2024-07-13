"""Definition of internal language for actions."""

from typing import Optional

from integral.expr import Expr
from integral.rules import Rule
from integral.context import Context
from integral import compstate
from integral.compstate import Calculation, StateItem, Goal, CompFile
from integral import poly

class Action:
    """Base class for actions."""
    pass


class ProveAction(Action):
    """Start a proof."""
    def __init__(self, expr: Expr):
        self.expr = expr

    def __str__(self):
        return "prove %s" % self.expr


class CalculateAction(Action):
    """Start a calculation."""
    def __init__(self, expr: Expr):
        self.expr = expr

    def __str__(self):
        return "calculate %s" % self.expr
    

class LHSAction(Action):
    """Perform a proof by working on the left hand side."""
    def __init__(self):
        pass

    def __str__(self):
        return "lhs:"
    

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
        
        elif isinstance(action, ProveAction):
            goal = compstate.Goal(self.comp_file, self.comp_file.ctx, action.expr)
            return ProveState(self, goal)
        
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
        if isinstance(action, LHSAction):
            proof = self.goal.proof_by_calculation()
            return CalculateState(self, proof.lhs_calc)
        else:
            raise StateException("Unknown action type %s" % type(action))

    def is_finished(self) -> bool:
        return self.goal.proof is not None and self.goal.proof.is_finished()

    def __str__(self):
        return "(proof)\n%s" % self.goal


class CalculateState(State):
    """State when performing a calculation."""
    def __init__(self, past: State, calc: Calculation):
        self.past = past
        self.calc = calc

    def process_action(self, action: Action) -> State:
        if isinstance(action, CalculateAction):
            raise StateException("Cannot start a new calculation within a calculation.")
        elif isinstance(action, RuleAction):
            self.calc.perform_rule(action.rule)
            return self
        else:
            raise StateException("Unknown action type %s" % type(action))

    def is_finished(self) -> bool:
        if not self.calc.steps:
            return False

        res = self.calc.steps[-1].res
        return res.is_evaluable() and poly.normalize(res, self.calc.ctx) == res

    def __str__(self):
        return "(calculate)\n%s" % self.calc
