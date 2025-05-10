"""Definition of internal language for actions."""

from typing import Optional, Iterable

from integral import expr
from integral.expr import Expr
from integral.conditions import Conditions


class Action:
    """Base class for actions."""
    def get_start_states(self) -> list[str]:
        raise NotImplementedError(f"get_start_states: {type(self)}")


class ImportsAction(Action):
    """Import other theories."""
    def __init__(self, theories: Iterable[str]):
        self.theories = tuple(theories)

    def __str__(self):
        return "imports " + ", ".join(self.theories)
    
    def get_start_states(self):
        return ["initial"]


class AxiomAction(Action):
    """State an axiom."""
    def __init__(self, expr: Expr, conditions: Iterable[Expr], attrs: Iterable[str]):
        self.expr = expr
        self.conditions = Conditions(conditions)
        self.attrs = tuple(attrs)

    def __str__(self):
        res = "axiom "
        if self.attrs:
            res += "[" + ', '.join(self.attrs) + "] "
        res += str(self.expr)
        if self.conditions:
            res += " for " + ', '.join(str(cond) for cond in self.conditions.data)
        return res

    def get_start_states(self) -> list[str]:
        return ["initial"]


class ProveAction(Action):
    """Start a proof."""
    def __init__(self, expr: Expr, conditions: Iterable[Expr], attrs: Iterable[str]):
        self.expr = expr
        self.conditions = Conditions(conditions)
        self.attrs = tuple(attrs)

    def __str__(self):
        res = "prove "
        if self.attrs:
            res += "[" + ', '.join(self.attrs) + "] "
        res += str(self.expr)
        if self.conditions:
            res += " for " + ', '.join(str(cond) for cond in self.conditions.data)
        return res

    def get_start_states(self) -> list[str]:
        return ["initial"]


class LetAction(Action):
    """Make a local definition."""
    def __init__(self, expr: Expr, conditions: Optional[Conditions] = None):
        self.expr = expr
        self.conditions = Conditions(conditions)

    def __str__(self):
        if self.conditions:
            return "let %s for %s" % (self.expr, ', '.join(str(cond) for cond in self.conditions.data))
        else:
            return "let %s" % self.expr

    def get_start_states(self) -> list[str]:
        return ["proof"]


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
        return ["initial"]


class AxiomDefineAction(Action):
    """Make an axiomatic definition."""
    def __init__(self, expr: Expr, conditions: Optional[Conditions] = None):
        self.expr = expr
        self.conditions = Conditions(conditions)

    def __str__(self):
        if self.conditions:
            return "axiom_define %s for %s" % (self.expr, ', '.join(str(cond) for cond in self.conditions.data))
        else:
            return "axiom_define %s" % self.expr

    def get_start_states(self) -> list[str]:
        return ["initial"]


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
    def __init__(self, rule):
        self.rule = rule

    def __str__(self):
        return str(self.rule)
    
    def get_start_states(self) -> list[str]:
        return ["calculation"]
