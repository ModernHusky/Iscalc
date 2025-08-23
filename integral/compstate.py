"""State of computation"""

from typing import List, Optional, Union

from integral.expr import Expr, Var, Const
from integral import rules, expr
from integral.rules import Rule, check_wellformed, ProofObligation
from integral.conditions import Conditions
from integral import condprover
from integral.context import Context, Definition, Identity
from integral.poly import normalize
from integral import utils


class CheckFinishedException(expr.IscalcException):
    """Exception raised when check finished fails."""
    def __init__(self, stack: tuple[str], msg: str):
        self.stack = stack + (msg,)

    def __str__(self):
        return "Use done when goal is not finished, goal stack:\n" + '\n'.join(self.stack)

    def to_json(self) -> dict:
        return {
            "class": "CheckFinishedException",
            "stack": list(self.stack)
        }
    
    @staticmethod
    def from_json(data: dict):
        stack = data['stack']
        return CheckFinishedException(tuple(stack[:-1]), stack[-1])


class StateException(expr.IscalcException):
    """Exception resulting from applying action to a state."""
    def __init__(self, kind: str, msg: str):
        self.kind = kind
        self.msg = msg

    def __str__(self):
        return f"{self.kind}: {self.msg}"

    def to_json(self) -> dict:
        return {
            "class": "StateException",
            "kind": self.kind,
            "msg": self.msg
        }
    
    @staticmethod
    def from_json(data: dict):
        return StateException(data["kind"], data["msg"])


class Label:
    def __init__(self, data):
        self.data = []
        if isinstance(data, str):
            split = data.split(".")
            for n in split:
                if n == '':
                    continue
                assert int(n) >= 1, "Label: non-positive value"
                self.data.append(int(n) - 1)
        elif isinstance(data, list):
            assert all(n >= 0 for n in data), "Label: negative value"
            self.data = list(data)
        elif isinstance(data, Label):
            self.data = data.data
        else:
            raise AssertionError("Label: unexpected type")

    @property
    def head(self):
        return self.data[0]

    @property
    def tail(self):
        return Label(self.data[1:])

    def empty(self):
        return len(self.data) == 0

    def __str__(self):
        res = ""
        for n in self.data:
            res += str(n + 1) + "."
        return res

    def __eq__(self, other):
        return isinstance(other, Label) and self.data == other.data

    def append(self, i: int) -> "Label":
        return Label(self.data + [i])


class StateItem:
    """Items in a state of computation"""
    ctx: Context

    def get_by_label(self, label: Label) -> "StateItem":
        """Return the object at the given label."""
        raise NotImplementedError

    def get_facts(self):
        """Return the list of facts in this item."""
        return []

    def clear(self):
        """Clear itself."""
        pass

    def is_finished(self):
        """Whether the proof in the item is finished. Default to true."""
        return True


class Goal(StateItem):
    """Goal to be proved.
    
    Attributes
    ----------
    parent: Optional[StateItem]
        parent of the goal
    ctx: Context
        initial context of the goal. Note this may be different from context
        in parent, as earlier subgoals and definitions in the parent are added.
    goal: Expr
        statement to be proved
    conds: Optional[Conditions]
        additional conditions of the goal

    """
    def __init__(self, parent: Optional[StateItem], ctx: Context, goal: Expr, *,
                 conds: Optional[Conditions] = None):
        self.parent = parent

        # Statement to be proved
        self.goal = goal

        # List of assumptions for the goal
        if conds is None:
            conds = Conditions()
        self.conds = conds

        # Initialize context
        self.ctx = Context(ctx)
        self.ctx.extend_vars(goal.get_vars())
        self.ctx.extend_condition(self.conds)

        # Initialize proof
        self.proof = None

        # List of local definitions
        self.definitions: list[Definition] = list()

        # List of subgoals, as (name, goal) pairs
        self.subgoals: list[tuple[str, Goal]] = list()

    def __str__(self):
        if self.is_finished():
            res = "Goal (finished)\n"
        else:
            res = "Goal\n"
        res += "  %s" % self.goal
        if self.conds.data:
            res += " for %s" % (", ".join(str(cond) for cond in self.conds.data))
        res += "\n"
        for definition in self.definitions:
            res += str(definition) + "\n"
        for n, subgoal in self.subgoals:
            res += "subgoal %s\n" % n
            res += str(subgoal)
        if self.proof is not None:
            res += str(self.proof)
        return res

    def print_entry(self, is_toplevel=True):
        if is_toplevel:
            if self.conds and self.conds.data:
                print("prove %s for %s" % (self.goal, ', '.join(str(cond) for cond in self.conds.data)))
            else:
                print("prove %s" % self.goal)
        for func_def in self.definitions:
            print(func_def)
        for n, subgoal in self.subgoals:
            if subgoal.conds and subgoal.conds.data:
                print("subgoal %s: %s for %s" % (n, subgoal.goal, ', '.join(str(cond) for cond in subgoal.conds.data)))
            else:
                print("subgoal %s: %s" % (n, subgoal.goal))
            subgoal.print_entry(is_toplevel=False)
        if isinstance(self.proof, CalculationProof):
            self.proof.print_entry()
        elif isinstance(self.proof, RewriteGoalProof):
            self.proof.print_entry()
        elif isinstance(self.proof, InductionProof):
            self.proof.print_entry()
        elif isinstance(self.proof, CaseProof):
            self.proof.print_entry()

    def __eq__(self, other):
        if not isinstance(other, Goal):
            return False
        return self.proof == other.proof

    def check_finished(self, stack: tuple[str]):
        goal_str = str(self.goal)
        if self.conds:
            goal_str += " for " + ', '.join(str(cond) for cond in self.conds.data)
        if self.proof is None:
            raise CheckFinishedException(stack, f"goal {goal_str} has no proof")
        
        proof_obligs: list[ProofObligation] = check_wellformed(self.goal, self.ctx)
        if proof_obligs:
            msg = f"goal {self.goal} is not wellformed."
            for i, obligation in enumerate(proof_obligs, 1):
                msg += f"\nObligation {i}\n"
                msg += utils.indent(str(obligation))
            raise CheckFinishedException(stack, msg)
        for n, subgoal in self.subgoals:
            subgoal.check_finished(stack + (f"subgoal {n}: {subgoal.goal}",))
        self.proof.check_finished(stack + (f"proof of {goal_str}",))

    def clear(self):
        self.proof = None

    def add_subgoal(self, name: str, expr: Expr, conds: Optional[list[Expr]] = None) -> "Goal":
        """Add subgoal with given name and expression."""

        # Form context of the subgoal by adding existing subgoal and definitions
        # in the current goal.
        conds = Conditions(conds)
        goal = Goal(self, self.ctx, expr, conds=conds)
        self.subgoals.append((name, goal))
        self.ctx = Context(self.ctx)
        self.ctx.add_subgoal(name, Identity(expr, conds=conds))
        return goal

    def add_definition(self, eq: Expr, conds: Optional[list[Expr]] = None):
        if not eq.is_equals():
            raise AssertionError(f"define: {eq}")

        self.ctx = Context(self.ctx)
        self.ctx.add_definition(eq, conds)

    def proof_by_rewrite_goal(self, *, begin: str) -> "RewriteGoalProof":
        if not isinstance(begin, str):
            raise StateException("RewriteGoalProof", "begin should be a string")
        self.proof = RewriteGoalProof(self, self.ctx, self.goal, start=begin)
        return self.proof

    def proof_by_calculation(self) -> "CalculationProof":
        self.proof = CalculationProof(self, self.ctx, self.goal)
        return self.proof

    def proof_by_induction(self, induct_var: str, start: int = 0) -> "InductionProof":
        self.proof = InductionProof(self, self.ctx, self.goal, induct_var, start=start)
        return self.proof

    def proof_by_case(self, split_cond: Expr) -> "CaseProof":
        self.proof = CaseProof(self, self.ctx, self.goal, split_cond=split_cond)
        return self.proof

    def get_by_label(self, label: Label):
        if label.empty():
            return self
        else:
            if self.proof is None:
                raise AssertionError("get_by_label: goal %s has no proof" % str(self.goal))
            return self.proof.get_by_label(label)

    def get_facts(self):
        return [self.goal]


class CalculationStep(StateItem):
    """A step in the calculation.
    
    Attributes
    ----------
    parent: Calculation
        the calculation this step is contained in.
    rule: Rule
        rule to be applied in this calculation.
    res: Expr
        result of this calculation step.
    id: int
        index of this step within the calculation.
    ctx: Context
        context of the calculation step.

    """
    def __init__(self, parent: "Calculation", rule: Rule, res: Expr, id: int):
        self.parent = parent
        self.rule = rule
        self.res = res
        self.id = id
        self.ctx = parent.ctx

    def __str__(self):
        return "%s (%s)" % (self.res, self.rule)

    def __eq__(self, other):
        if not isinstance(other, CalculationStep):
            return False
        return self.rule.name == other.rule.name and str(self.res) == str(other.res)

    def clear(self):
        self.parent.clear(id=self.id)

    def perform_rule(self, rule: Rule):
        self.parent.perform_rule(rule, self.id)

    def perform_rules(self, calc_rules: tuple[Rule]):
        self.parent.perform_rules(calc_rules, self.id)


class Calculation(StateItem):
    """Calculation starting from an expression.

    Attributes
    ----------
    parent: Optional[StateItem]
        parent of the calculation
    ctx: Context
        context of the calculation
    start: Expr
        starting expression.
    connection_symbol: str
        one of '=' (for equality chaining) and '==>' (for rewriting goal)
    steps: list[CalculationStep]:
        list of steps in the calculation.
    conds: Conditions
        conditions under which the calculation is carried out.

    """
    def __init__(self, parent: Optional[StateItem], ctx: Context, start: Expr, *,
                 connection_symbol: str = '=', conds: Optional[Conditions] = None):
        self.parent = parent
        self.start = start

        self.steps: list[CalculationStep] = []
        if conds is None:
            conds = Conditions()
        self.conds = conds
        self.connection_symbol = connection_symbol

        self.ctx = Context(ctx)
        self.ctx.extend_vars(start.get_vars())
        if conds is not None:
            self.ctx.extend_condition(self.conds)

        self.subgoals = list()

    def check_wellformed(self):
        proof_obligs: list[ProofObligation] = check_wellformed(self.start, self.ctx)
        if proof_obligs:
            msg = f"start {self.start} is not wellformed."
            for i, obligation in enumerate(proof_obligs, 1):
                msg += f"\nObligation {i}\n"
                msg += utils.indent(str(obligation))
            raise CheckFinishedException(tuple(), msg)

    def add_subgoal(self, name: str, expr: Expr, conds: Optional[list[Expr]] = None) -> Goal:
        """Add subgoal with given name and expression."""

        # Form context of the subgoal by adding existing subgoal and definitions
        # in the current goal.
        conds = Conditions(conds)
        goal = Goal(self, self.ctx, expr, conds=Conditions(conds))
        self.subgoals.append((name, goal))
        self.ctx = Context(self.ctx)
        self.ctx.add_subgoal(name, Identity(expr, conds=conds))
        return goal

    def __eq__(self, other):
        if not isinstance(other, Calculation):
            return False
        return self.steps == other.steps

    def __str__(self):
        res = "  " + str(self.start) + "\n"
        for step in self.steps:
            res += self.connection_symbol + " %s\n" % step
        return res

    def print_entry(self):
        print("calculate %s" % self.start)
        for step in self.steps:
            print(str(step.rule))

    def clear(self, id: int = 0):
        self.steps = self.steps[:id]

    def add_step(self, step: CalculationStep):
        """Add the given step to the computation."""
        self.steps.append(step)

    @property
    def last_expr(self) -> Expr:
        """Last expression of the calculation."""
        if self.steps:
            return self.steps[-1].res
        else:
            return self.start

    def perform_rule(self, rule: Rule, id: Optional[int] = None):
        """Perform the given rule on the current expression."""
        if id is not None:
            # Cut off later steps
            self.steps = self.steps[:id + 1]
        else:
            id = len(self.steps) - 1

        e = self.last_expr
        ctx = Context(self.ctx)
        for n, subgoal in self.subgoals:
            ctx.subgoals[n] = Identity(subgoal.goal, conds=subgoal.conds)
        cur_e = self.start
        for step in self.steps:
            ctx = step.rule.update_context(cur_e, ctx)
            cur_e = step.res
        new_e = rule.eval(e, ctx)
        if str(new_e) == str(e):  # check equality ignoring alpha equivalence
            raise rules.RuleException(rules.get_rule_name(rule), f"Applying the rule has no effect: {str(rule)}")
        step = CalculationStep(self, rule, new_e, id + 1)
        self.add_step(step)

    def perform_rules(self, calc_rules: tuple[Rule], id: Optional[int] = None):
        for rule in calc_rules:
            self.perform_rule(rule)

    def get_by_label(self, label: Label) -> "StateItem":
        if label.empty():
            return self
        elif label.tail.empty():
            return self.steps[label.head]
        else:
            raise AssertionError("get_by_label: invalid label")


class CalculationProof(StateItem):
    """Proof for an equation by calculation.

    The proof consists of calculation of left and right sides.

    """
    def __init__(self, parent, ctx: Context, goal: Expr):
        self.parent = parent
        self.goal = goal
        self.ctx = Context(ctx)
        self.calcs: list[Calculation] = []
        if expr.is_compare(goal):
            self.predicate = goal.op
            if isinstance(parent, Goal):
                self.calcs.append(Calculation(self, self.ctx, goal.args[0], conds=parent.conds))
                self.calcs.append(Calculation(self, self.ctx, goal.args[1], conds=parent.conds))
            else:
                raise NotImplementedError
        elif expr.is_fun(goal) and goal.func_name == "converges":
            self.predicate = goal.func_name
            assert isinstance(parent, Goal)
            self.calcs.append(Calculation(self, self.ctx, goal.args[0], conds=parent.conds))
        else:
            # INT u:[0,oo]. 1 / (u ^ 2 + 1) ^ 2 = pi / 4
            if expr.is_integral(goal) and expr.is_equals(goal.body):
                new_expr = expr.Op("=", \
                                   expr.Integral(goal.var, goal.lower, goal.upper, goal.body.lhs), \
                                    goal.body.rhs)
                raise StateException("CalculationProof", f"The equality operator (=) has higher precedence than the integral operator (INT), so the integral must be enclosed in parentheses. The goal {goal} should be modified to ({new_expr.lhs})={new_expr.rhs}.")
            # D x. f(x) = g(x) ===> (D x. f(x)) = g(x)
            elif expr.is_deriv(goal) and expr.is_equals(goal.body):
                new_expr = expr.Op("=", \
                                   expr.Deriv(goal.var, goal.body.lhs), \
                                   goal.body.rhs)
                raise StateException("CalculationProof",
                                     f"The equality operator (=) has higher precedence than the derivative operator (D), so the integral must be enclosed in parentheses. The goal {goal} should be modified to ({new_expr.lhs})={new_expr.rhs}.")
            else:
                raise StateException("CalculationProof", "unknown form of goal.")

    def __eq__(self, other):
        return isinstance(other, CalculationProof) and \
            self.calcs == other.calcs and self.goal == other.goal

    def __str__(self):
        if self.is_finished():
            res = "Proof by calculation (finished)\n"
        else:
            res = "Proof by calculation\n"
        for calc in self.calcs:
            if calc.steps:
                res += str(calc)
        return res
    
    def print_entry(self):
        if expr.is_fun(self.goal) and self.goal.func_name == "converges":
            print("arg:")
            for step in self.arg_calc.steps:
                print("    " + str(step.rule))
        else:
            if self.lhs_calc and self.lhs_calc.steps:
                print("lhs:")
                for step in self.lhs_calc.steps:
                    print("    " + str(step.rule))
            if self.rhs_calc and self.rhs_calc.steps:
                print("rhs:")
                for step in self.rhs_calc.steps:
                    print("    " + str(step.rule))
        print("done")

    @property
    def lhs_calc(self) -> Calculation:
        if not self.goal.is_compare():
            if expr.is_fun(self.goal) and self.goal.func_name == "converges":
                raise StateException(
                    "Calculate",
                    f"Action type lhs: cannot be performed in calculate state when proving the convergence goal({str(self.goal)}).")
            else:
                raise StateException("CalculationProof", "currently only support equation goals.")
        return self.calcs[0]

    @property
    def rhs_calc(self) -> Calculation:
        assert self.goal.is_compare()
        return self.calcs[1]

    @property
    def arg_calc(self) -> Calculation:
        if not expr.is_fun(self.goal):
            raise StateException("CalculationProof", "Cannot perform the arg: operation in the lhs: context.")
        return self.calcs[0]

    def is_finished(self):
        if self.predicate == '=':
            return normalize(self.lhs_calc.last_expr, self.ctx) == \
                   normalize(self.rhs_calc.last_expr, self.ctx)
        elif self.predicate == '>':
            return self.ctx.is_greater(self.lhs_calc.last_expr, self.rhs_calc.last_expr)
        elif self.predicate == '<':
            return self.ctx.is_less(self.lhs_calc.last_expr, self.rhs_calc.last_expr)
        elif self.predicate == '<=':
            return self.ctx.is_less_eq(self.lhs_calc.last_expr, self.rhs_calc.last_expr)
        elif self.predicate == '>=':
            return self.ctx.is_greater_eq(self.lhs_calc.last_expr, self.rhs_calc.last_expr)
        elif self.predicate == '!=':
            return self.ctx.is_not_equal(self.lhs_calc.last_expr, self.rhs_calc.last_expr)
        elif self.predicate == 'converges':
            return rules.check_converge(self.arg_calc.last_expr, self.ctx)
        raise NotImplementedError(f"predicate: {self.predicate}")

    def check_finished(self, stack: tuple[str]):
        if self.predicate in ('=', '>', '<', '<=', '>=', '!='):
            lhs = self.lhs_calc.last_expr
            rhs = self.rhs_calc.last_expr
            if self.predicate == '=' and normalize(lhs, self.ctx) != normalize(rhs, self.ctx):
                raise CheckFinishedException(stack, f"calculation: {lhs} != {rhs}")
            if self.predicate == '>' and not self.ctx.is_greater(lhs, rhs):
                raise CheckFinishedException(stack, f"calculation: {lhs} > {rhs}")
            if self.predicate == '<' and not self.ctx.is_less(lhs, rhs):
                raise CheckFinishedException(stack, f"calculation: {lhs} < {rhs}")
            if self.predicate == '<=' and not self.ctx.is_less_eq(lhs, rhs):
                raise CheckFinishedException(stack, f"calculation: {lhs} <= {rhs}")
            if self.predicate == '>=' and not self.ctx.is_greater_eq(lhs, rhs):
                raise CheckFinishedException(stack, f"calculation: {lhs} >= {rhs}")
            if self.predicate == '!=' and not self.ctx.is_not_equal(lhs, rhs):
                raise CheckFinishedException(stack, f"calculation: {lhs} != {rhs}")
        elif self.predicate == 'converges':
            e1 = normalize(self.arg_calc.last_expr, self.ctx)
            e2 = normalize(-e1, self.ctx)
            if not rules.check_converge(e1, self.ctx) and not rules.check_converge(e2, self.ctx):
                raise CheckFinishedException(stack, f"calculation: {e1} does not converge")
        else:
            raise NotImplementedError(f"predicate: {self.predicate}")

    def clear(self):
        for calc in self.calcs:
            calc.clear()

    def get_by_label(self, label: Label):
        if label.empty():
            return self
        elif label.head < len(self.calcs):
            return self.calcs[label.head].get_by_label(label.tail)
        else:
            raise AssertionError("get_by_label: invalid label")


class InductionProof(StateItem):
    """Proof for an equation by induction on natural numbers.

    This breaks the equation goal into two goals, corresponding to the
    base case and inductive case.

    """

    def __init__(self, parent: Goal, ctx: Context, goal: Expr, induct_var: str,
                 *, start: Union[int, Expr] = 0):
        if not goal.is_equals():
            raise StateException("InductionProof", "currently only support equation goals.")

        self.parent = parent
        self.goal = goal
        self.induct_var = induct_var
        self.ctx = Context(ctx)

        if isinstance(start, int):
            self.start = Const(start)
        elif isinstance(start, Expr):
            self.start = start
        else:
            raise NotImplementedError

        if not ctx.check_condition(expr.Fun("isInt", expr.Var(self.induct_var))):
            raise StateException("InductionProof", f"induction variable {self.induct_var} is not integer")

        if not ctx.check_condition(expr.Op(">=", expr.Var(self.induct_var), self.start)):
            raise StateException("InductionProof", f"condition {self.induct_var} >= {self.start} does not hold")

        # Base case: n = start
        base_goal_ctx = Context(self.ctx)
        eq0 = normalize(goal.subst(induct_var, self.start), base_goal_ctx)
        self.base_case = Goal(self, base_goal_ctx, eq0)

        # Induction case
        n1 = Var(induct_var) + Const(1)
        induct_goal_ctx = Context(self.ctx)
        eqI = normalize(goal.subst(induct_var, n1), induct_goal_ctx)
        induct_goal_ctx.add_induct_hyp(self.goal)
        self.induct_case = Goal(self, induct_goal_ctx, eqI)

    def __eq__(self, other):
        if not isinstance(other, InductionProof):
            return False
        return self.start == other.start and self.induct_case == other.induct_case and \
               self.base_case == other.base_case

    def __str__(self):
        if self.is_finished():
            res = "Proof by induction on %s (finished)\n" % self.induct_var
        else:
            res = "Proof by induction on %s\n" % self.induct_var
        res += "Base case: %s\n" % self.base_case.goal
        res += str(self.base_case)
        res += "Induct case: %s\n" % self.induct_case.goal
        res += str(self.induct_case)
        return res
    
    def print_entry(self):
        if self.start == Const(0):
            print("induction on %s" % self.induct_var)
        else:
            print("induction on %s starting from %s" % (self.induct_var, self.start))
        if self.base_case.proof:
            print("base:")
            self.base_case.print_entry(is_toplevel=False)
        if self.induct_case.proof:
            print("induct:")
            self.induct_case.print_entry(is_toplevel=False)
        print("done")

    def is_finished(self):
        return self.base_case.is_finished() and self.induct_case.is_finished()

    def check_finished(self, stack: tuple[str]):
        self.base_case.check_finished(stack + ("base case",))
        self.induct_case.check_finished(stack + ("induct case",))

    def clear(self):
        self.base_case.clear()
        self.induct_case.clear()

    def get_by_label(self, label: Label):
        if label.empty():
            return self
        elif label.head == 0:
            return self.base_case.get_by_label(label.tail)
        elif label.head == 1:
            return self.induct_case.get_by_label(label.tail)
        else:
            raise AssertionError("get_by_label: invalid label")


class CaseProof(StateItem):
    """Prove an equation by cases.

    If split_cond is a condition, the two cases correspond to split_cond
    being true and false, respectively.

    If split_cond is an expression a, the three cases correspond to
    a > 0, a = 0, and a < 0.

    """
    def __init__(self, parent: Goal, ctx: Context, goal: Expr, *, split_cond: Expr):
        self.parent = parent
        self.goal = goal
        self.ctx = Context(ctx)
        self.split_cond = split_cond
        self.split_type = ""
        self.cases: List[Goal] = []
        assert isinstance(parent, Goal)

        if split_cond.is_compare():
            self.split_type = "two-way"
            # Case 1:
            conds1 = Conditions()
            conds1.add_condition(split_cond)
            conds1.update(parent.conds)
            self.cases.append(Goal(self, self.ctx, goal, conds=conds1))

            # Case 2:
            conds2 = Conditions()
            conds2.add_condition(expr.neg_expr(split_cond))
            conds2.update(parent.conds)
            self.cases.append(Goal(self, self.ctx, goal, conds=conds2))

        else:
            self.split_type = "three-way"
            # Case 1:
            conds1 = Conditions()
            conds1.add_condition(expr.Op("<", split_cond, Const(0)))
            conds1.update(parent.conds)
            self.cases.append(Goal(self, self.ctx, goal, conds=conds1))

            # Case 2:
            conds2 = Conditions()
            conds2.add_condition(expr.Op("=", split_cond, Const(0)))
            conds2.update(parent.conds)
            self.cases.append(Goal(self, self.ctx, goal, conds=conds2))

            # Case 3:
            conds3 = Conditions()
            conds3.add_condition(expr.Op(">", split_cond, Const(0)))
            conds3.update(parent.conds)
            self.cases.append(Goal(self, self.ctx, goal, conds=conds3))

    def __eq__(self, other):
        if not isinstance(other, CaseProof):
            return False
        return self.goal == other.goal and self.split_type == other.split_type and \
               self.split_cond == other.split_cond and self.cases == other.cases

    def print_entry(self):
        if self.split_type == "two-way":
            print("case analysis on %s" % self.split_cond)
            if self.cases[0].proof:
                print("case true:")
                self.cases[0].print_entry(is_toplevel=False)
            if self.cases[1].proof:
                print("case false:")
                self.cases[1].print_entry(is_toplevel=False)
            print("done")
        elif self.split_type == "three-way":
            print("case analysis on %s" % self.split_cond)
            if self.cases[0].proof:
                print("case negative:")
                self.cases[0].print_entry(is_toplevel=False)
            if self.cases[1].proof:
                print("case zero:")
                self.cases[1].print_entry(is_toplevel=False)
            if self.cases[2].proof:
                print("case positive:")
                self.cases[2].print_entry(is_toplevel=False)
            print("done")
        else:
            raise NotImplementedError(f"split_type = {self.split_type}")

    def __str__(self):
        if self.is_finished():
            res = "Proof by cases (finished)\n"
        else:
            res = "Proof by cases\n"
        for i, case in enumerate(self.cases):
            res += "case %d: %s for %s\n" % (i + 1, case.goal, case.conds)
            res += str(case)
        return res

    def is_finished(self):
        if self.split_type == "two-way":
            return self.cases[0].is_finished() and self.cases[1].is_finished()
        elif self.split_type == "three-way":
            all_conds = condprover.init_all_conds(self.parent.conds)
            if not condprover.check_cond(expr.Op(">=", self.split_cond, Const(0)), all_conds, dict()) and \
                not self.cases[0].is_finished():
                return False
            if not condprover.check_cond(expr.Op("!=", self.split_cond, Const(0)), all_conds, dict()) and \
                not self.cases[1].is_finished():
                return False
            if not condprover.check_cond(expr.Op("<=", self.split_cond, Const(0)), all_conds, dict()) and \
                not self.cases[2].is_finished():
                return False
            return True
        else:
            raise NotImplementedError(f"split_type = {self.split_type}")

    def check_finished(self, stack: tuple[str]):
        if self.split_type == "two-way":
            self.cases[0].check_finished(stack + (f"{self.split_cond} true branch",))
            self.cases[1].check_finished(stack + (f"{self.split_cond} false branch",))

        elif self.split_type == "three-way":
            all_conds = condprover.init_all_conds(self.parent.conds)
            if not condprover.check_cond(expr.Op(">=", self.split_cond, Const(0)), all_conds, dict()):
                self.cases[0].check_finished(stack + (f"{self.split_cond} < 0 case",))
            if not condprover.check_cond(expr.Op("!=", self.split_cond, Const(0)), all_conds, dict()):
                self.cases[1].check_finished(stack + (f"{self.split_cond} = 0 case",))
            if not condprover.check_cond(expr.Op("<=", self.split_cond, Const(0)), all_conds, dict()):
                self.cases[2].check_finished(stack + (f"{self.split_cond} > 0 case",))

        else:
            raise NotImplementedError(f"split_type = {self.split_type}")

    def clear(self):
        for case in self.cases:
            case.clear()

    def get_by_label(self, label: Label):
        if label.empty():
            return self
        elif label.head < len(self.cases):
            return self.cases[label.head].get_by_label(label.tail)
        else:
            raise AssertionError("get_by_label: invalid label")


class RewriteGoalProof(StateItem):
    """Prove an equation by transforming an initial equation.

    Attributes
    ----------
    parent: StateItem
        parent of the proof
    ctx: Context
        context of the proof
    goal: Expr
        goal to be proved
    start: str
        start of rewriting, should be the name of a previous subgoal
        
    """
    def __init__(self, parent: StateItem, ctx: Context, goal: Expr, *, start: str):
        if not goal.is_equals():
            raise StateException("RewriteGoalProof", f"goal {goal} is not an equality.")
        self.parent = parent
        self.goal = goal
        self.ctx = Context(ctx)
        self.start = start
        start_goal = ctx.get_subgoal(start)
        if not start_goal:
            raise StateException("RewriteGoalProof", f"start {start} not found")
        self.begin = Calculation(self, ctx, start_goal.expr, connection_symbol='==>',
                                 conds=start_goal.conds)

    def __eq__(self, other):
        if not isinstance(other, RewriteGoalProof):
            return False
        return self.goal == other.goal and self.begin == other.begin

    def print_entry(self):
        print("from %s:" % self.start)
        for step in self.begin.steps:
            print("    " + str(step.rule))
        print("done")

    def is_finished(self):
        f1 = normalize(self.begin.last_expr.lhs, self.ctx) == normalize(self.goal.lhs, self.ctx)
        f2 = normalize(self.begin.last_expr.rhs, self.ctx) == normalize(self.goal.rhs, self.ctx)
        return f1 and f2
    
    def check_finished(self, stack: tuple[str]):
        calc_lhs = normalize(self.begin.last_expr.lhs, self.ctx)
        goal_lhs = normalize(self.goal.lhs, self.ctx)
        if calc_lhs != goal_lhs:
            raise CheckFinishedException(stack, f"rewrite goal lhs not equal: {calc_lhs} != {goal_lhs}")

        calc_rhs = normalize(self.begin.last_expr.rhs, self.ctx)
        goal_rhs = normalize(self.goal.rhs, self.ctx)
        if calc_rhs != goal_rhs:
            raise CheckFinishedException(stack, f"rewrite goal rhs not equal: {calc_rhs} != {goal_rhs}")

    def clear(self):
        self.begin.clear()

    def __str__(self):
        if self.is_finished():
            res = "Proof by rewriting equation (finished)\n"
        else:
            res = "Proof by rewriting equation\n"

        res += str(self.begin)
        return res

    def get_by_label(self, label: Label):
        if label.empty() or len(label.data) == 1:
            return self
        elif not label.tail.empty():
            return self.begin.steps[label.tail.head]
        else:
            raise AssertionError("get_by_label: invalid label")


class CompFile:
    """Represent a file containing multiple StateItem objects.

    ctx - initial context of the file.
        either a Context or a string, specifying the base context or
        file name.

    """
    def __init__(self, ctx: Union[Context, str]):
        if isinstance(ctx, str):
            self.ctx = Context()
            self.ctx.load_book(ctx)
        else:
            self.ctx = ctx
        self.content: list[StateItem] = []

    def __eq__(self, other):
        return isinstance(other, CompFile) and self.content == other.content

    def __str__(self):
        res = ""
        for st in self.content:
            res += str(st)
        return res

    def add_definition(self, funcdef: Expr, *, conds: list[Expr] = None) -> Definition:
        """Add a function definition."""
        self.content.append(Definition(funcdef, Conditions(conds)))
        return self.content[-1]

    def add_calculation(self, calc: Expr, *, conds: list[Expr] = None) -> Calculation:
        """Add a calculation."""
        self.content.append(Calculation(self, self.ctx, calc, conds=Conditions(conds)))
        return self.content[-1]

    def add_goal(self, goal: Union[str, Expr, Goal], *,
                 conds: Optional[List[Union[str, Expr]]] = None) -> Goal:
        """Add a goal."""
        self.content.append(Goal(self, self.ctx, goal, conds))
        return self.content[-1]

    def add_item(self, item: StateItem):
        """Add item of arbitrary type"""
        self.content.append(item)

    def get_item_label(self, item: StateItem):
        res = None

        def rec(root: Union[CompFile, StateItem], loc: Label):
            nonlocal res, item
            if res != None:
                return
            if root == item:
                res = Label(loc)
            elif isinstance(root, CompFile):
                for idx, st in enumerate(self.content):
                    rec(st, loc.append(idx))
            elif isinstance(root, Goal):
                for n, subgoal in root.subgoals:
                    rec(subgoal, loc.append(int(n)))
                rec(root.proof, loc.append(0))
            elif isinstance(root, RewriteGoalProof):
                rec(root.begin, loc.append(0))
            elif isinstance(root, CalculationProof):
                rec(root.calcs[0], loc.append(0))
                rec(root.calcs[1], loc.append(1))
            elif isinstance(root, InductionProof):
                rec(root.base_case, loc.append(0))
                rec(root.induct_case, loc.append(1))
            elif isinstance(root, CaseProof):
                for i, c in enumerate(root.cases):
                    rec(c, loc.append(i))
            elif isinstance(root, Definition) or isinstance(root, CalculationStep):
                pass
            elif isinstance(root, Calculation):
                for i, step in enumerate(root.steps):
                    rec(root.steps[i], loc.append(i))
            else:
                print(type(root))
                raise NotImplementedError

        rec(self, Label(""))
        return res

    def get_by_label(self, label: Label):
        def rec(root: Union[CompFile, StateItem], loc: Label):
            if loc == Label(""):
                return root
            if isinstance(root, CompFile):
                return rec(root.content[loc.head], loc.tail)
            elif isinstance(root, Goal):
                return rec(root.proof, loc.tail)
            else:
                raise NotImplementedError

        return rec(self, label)
