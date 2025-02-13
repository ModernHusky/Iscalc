"""State of computation"""
from typing import Dict, List, Optional, Tuple, Union

from integral.expr import Expr, Var, Const
from integral import rules, expr
from integral.rules import Rule, check_wellformed
from integral.conditions import Conditions
from integral.context import Context, Identity
from integral import latex
from integral import parser
from integral.poly import normalize


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

    def append(self, i: int) -> "Location":
        return Label(self.data + [i, ])


class StateItem:
    """Items in a state of computation"""
    ctx: Context

    def export(self):
        """Obtain the JSON representation of the item."""
        raise NotImplementedError

    def export_book(self):
        """Obtain the JSON representation of the item in the book file."""
        raise NotImplementedError

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


class FuncDef(StateItem):
    """Introduce a new function definition."""

    def __init__(self, parent: "CompFile", ctx: Context, eq: Expr, conds: Optional[Conditions] = None):
        if not eq.is_equals():
            raise AssertionError("FuncDef: input should be an equation")

        self.parent = parent
        self.ctx = ctx
        self.eq = eq
        if expr.is_fun(self.eq.lhs):
            self.symb = self.eq.lhs.func_name
            self.args = self.eq.lhs.args
        elif expr.is_var(self.eq.lhs):
            self.symb = self.eq.lhs.name
            self.args = []
        else:
            raise AssertionError("FuncDef: left side of equation must be variable or function")
        self.body = self.eq.rhs

        if any(not expr.is_var(arg) for arg in self.args) or len(self.args) != len(set(self.args)):
            raise AssertionError("FuncDef: arguments should be distinct variables")

        if conds is None:
            conds = Conditions()
        self.conds = conds

    def __str__(self):
        res = "Definition\n"
        res += "  %s\n" % self.eq
        return res

    def __eq__(self, other):
        return isinstance(other, FuncDef) and self.eq == other.eq and self.conds == other.conds

    def export(self):
        res = {
            "type": "FuncDef",
            "eq": str(self.eq),
            "latex_lhs": latex.convert_expr(self.eq.lhs),
            "latex_eq": latex.convert_expr(self.eq)
        }
        if self.conds.data:
            res["conds"] = self.conds.export()
        return res

    def export_book(self):
        p = self.parent
        while (not isinstance(p, CompFile)):
            p = p.parent
        res = {
            "type": "definition",
            "expr": str(self.eq),
            "path": p.name
        }
        if self.conds.data:
            res["conds"] = [str(cond) for cond in self.conds.data]
        return res

    def get_by_label(self, label: Label):
        if not label.empty():
            raise AssertionError("get_by_label: invalid label")
        return self

    def get_facts(self):
        return [self.eq]


class Goal(StateItem):
    """Goal to be proved."""

    def __init__(self, parent: Union['CompFile', StateItem], ctx: Context, goal: Expr, *,
                 conds: Optional[Conditions] = None):
        self.parent = parent

        # Statement to be proved
        self.goal = goal
        # List of assumptions for the goal
        if conds is None:
            conds = Conditions()
        self.conds = conds

        self.proof = None
        self.ctx = ctx

        self.ctx.extend_vars(goal.get_vars())
        self.ctx.extend_condition(self.conds)

        # Check well-formedness of the goal
        proof_obligations_raw = check_wellformed(goal, self.ctx)
        self.proof_obligations = []
        for oblig in proof_obligations_raw:
            found = False
            for _, subgoal in self.ctx.subgoals.items():
                if subgoal.covers_obligation(oblig):
                    found = True
                    break
            if not found:
                self.proof_obligations.append(oblig)
        self.wellformed = (len(self.proof_obligations) == 0)

        # List of subgoals
        self.subgoals: List[Tuple[str, Goal]] = list()

        # List of temporary definitions
        self.definitions: List[FuncDef] = list()

    def __str__(self):
        if self.is_finished():
            res = "Goal (finished)\n"
        else:
            res = "Goal\n"
        res += "  %s\n" % self.goal
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
            if func_def.conds and func_def.conds.data:
                print("define %s for %s" % (func_def.eq, ', '.join(str(cond) for cond in func_def.conds.data)))
            else:
                print("define %s" % func_def.eq)
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

    def is_finished(self):
        # all conds are satisfied under context of proof
        if self.proof == None:
            return False
        if not self.wellformed:
            print("%s, %s" % (self.goal, self.wellformed))
            return False
        for n, subgoal in self.subgoals:
            if not subgoal.is_finished():
                return False
        return self.proof.is_finished()

    def clear(self):
        self.proof = None

    def export(self):
        res = {
            "type": "Goal",
            "goal": str(self.goal),
            "latex_goal": latex.convert_expr(self.goal),
            "finished": self.is_finished(),
        }
        if self.proof:
            res['proof'] = self.proof.export()
        if self.conds.data:
            res['conds'] = self.conds.export()
        if self.subgoals:
            res['subgoals'] = [{'name': name, 'goal': goal.export()}
                               for name, goal in self.subgoals]
        if not self.wellformed:
            res['wellformed'] = False
            res['obligations'] = [p.export() for p in self.proof_obligations]
        return res

    def export_book(self):
        p = self.parent
        while (not isinstance(p, CompFile)):
            p = p.parent
        res = {
            "type": "problem",
            "expr": str(self.goal),
            "path": p.name
        }
        if self.conds.data:
            res["conds"] = [str(cond) for cond in self.conds.data]
        return res
    
    def covers_obligation(self, oblig: rules.ProofObligation) -> bool:
        # List of conditions is a subset of conditions on obligation
        for cond in self.conds.data:
            if cond not in oblig.conds.data:
                return False

        # Satisfies the goal in one branch
        for branch in oblig.branches:
            if len(branch.exprs) == 1 and self.goal == branch.exprs[0]:
                return True
            
        return False

    def add_subgoal(self, name: str, expr: Union[str, Expr],
                    conds: Optional[List[Union[str, Expr]]] = None) -> "Goal":
        ctx = Context(self.ctx)
        for n, subgoal in self.subgoals:
            ctx.subgoals[n] = Identity(subgoal.goal, conds=subgoal.conds)
        for funcdef in self.definitions:
            ctx.add_definition(funcdef.eq, funcdef.conds)
        if isinstance(expr, str):
            expr = parser.parse_expr(expr)
        goal = Goal(self, ctx, expr, conds=Conditions(conds))
        self.subgoals.append((name, goal))

        # Recheck wellformedness conditions
        ctx = Context(ctx)
        ctx.subgoals[name] = Identity(expr, conds=Conditions(conds))
        proof_obligations_raw = check_wellformed(self.goal, ctx)
        self.proof_obligations = []
        for oblig in proof_obligations_raw:
            found = False
            for _, subgoal in self.subgoals:
                if subgoal.covers_obligation(oblig):
                    found = True
                    break
            if not found:
                self.proof_obligations.append(oblig)
        self.wellformed = (len(self.proof_obligations) == 0)

        return self.subgoals[-1][1]

    def add_definition(self, expr: Union[str, Expr],
                       conds: Optional[List[Union[str, Expr]]] = None) -> FuncDef:
        if isinstance(expr, str):
            expr = parser.parse_expr(expr)
        self.definitions.append(FuncDef(self, self.ctx, expr, conds=Conditions(conds)))
        return self.definitions[-1]

    def proof_by_rewrite_goal(self, *, begin):
        if not isinstance(begin, str):
            raise AssertionError("RewriteGoalProof: begin should be a string")
        ctx = Context(self.ctx)
        for n, subgoal in self.subgoals:
            ctx.subgoals[n] = Identity(subgoal.goal, conds=subgoal.conds)
        for funcdef in self.definitions:
            ctx.add_definition(funcdef.eq, funcdef.conds)
        self.proof = RewriteGoalProof(self, ctx, self.goal, start=begin)
        return self.proof

    def proof_by_calculation(self):
        ctx = Context(self.ctx)
        for n, subgoal in self.subgoals:
            ctx.subgoals[n] = Identity(subgoal.goal, conds=subgoal.conds)
        for funcdef in self.definitions:
            ctx.add_definition(funcdef.eq, funcdef.conds)
        self.proof = CalculationProof(self, ctx, self.goal)
        return self.proof

    def proof_by_induction(self, induct_var: str, start: int = 0):
        ctx = Context(self.ctx)
        for n, subgoal in self.subgoals:
            ctx.subgoals[n] = Identity(subgoal.goal, conds=subgoal.conds)
        for funcdef in self.definitions:
            ctx.add_definition(funcdef.eq, funcdef.conds)
        self.proof = InductionProof(self, ctx, self.goal, induct_var, start=start)
        return self.proof

    def proof_by_case(self, split_cond: Expr):
        ctx = Context(self.ctx)
        for n, subgoal in self.subgoals:
            ctx.subgoals[n] = Identity(subgoal.goal, conds=subgoal.conds)
        for funcdef in self.definitions:
            ctx.add_definition(funcdef.eq, funcdef.conds)
        self.proof = CaseProof(self, ctx, self.goal, split_cond=split_cond)
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
    parent (Calculation): the calculation this step is contained in.
    rule (Rule): rule to be applied in this calculation.
    res (Expr): result of this calculation step.
    id (int): index of this step within the calculation.
    ctx (Context): context of the calculation step.    

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

    def export(self):
        res = {
            "type": "CalculationStep",
            "rule": self.rule.export(),
            "res": str(self.res),
            "latex_res": latex.convert_expr(self.res)
        }
        return res

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
    parent: parent of the calculation, either a StateItem or CompFile.
    start (Expr): starting expression.
    steps: list of steps in the calculation.
    conds: (optional) a list of conditions under which the calculation
        is carried out.
    connection_symbol: one of '=' and '==>'
    ctx: current context (including existing identities, conditions,
        fixed variables, etc).

    """
    def __init__(self, parent, ctx: Context, start: Expr, *,
                 connection_symbol='=', conds: Optional[Conditions] = None):
        self.parent = parent
        self.start = start
        self.steps: List[CalculationStep] = []
        if conds is None:
            conds = Conditions()
        self.conds = conds
        self.connection_symbol = connection_symbol

        self.ctx = ctx
        self.ctx.extend_vars(start.get_vars())
        if conds is not None:
            self.ctx.extend_condition(self.conds)

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

    def export(self):
        res = {
            "type": "Calculation",
            "start": str(self.start),
            "latex_start": latex.convert_expr(self.start),
            "steps": [step.export() for step in self.steps],
        }
        if self.conds.data:
            res["conds"] = self.conds.export()
        return res

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
        cur_e = self.start
        for step in self.steps:
            ctx = step.rule.update_context(cur_e, ctx)
            cur_e = step.res
        new_e = rule.eval(e, ctx)
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

    def parse_expr(self, s: str) -> Expr:
        return parser.parse_expr(s)


class CalculationProof(StateItem):
    """Proof for an equation by calculation.

    The proof consists of calculation of left and right sides.

    """

    def __init__(self, parent, ctx: Context, goal: Expr):
        self.parent = parent
        self.goal = goal
        self.ctx = ctx
        self.calcs = []
        if goal.is_compare():
            self.predicate = goal.op
            if isinstance(parent, Goal):
                self.calcs.append(Calculation(self, self.ctx, self.goal.args[0]))
                self.calcs.append(Calculation(self, self.ctx, self.goal.args[1]))
            else:
                raise NotImplementedError
        elif expr.is_fun(goal) and goal.func_name == "converges":
            self.predicate = goal.func_name
            assert isinstance(parent, Goal)
            goal = parent
            self.calcs.append(Calculation(self, self.ctx, self.goal.args[0], conds=goal.conds))
        else:
            raise AssertionError("CalculationProof: unknown form of goal.")

    def __eq__(self, other):
        if not isinstance(other, CalculationProof):
            return False
        return self.calcs == other.calcs and self.goal == other.goal

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
        assert self.goal.is_compare()
        return self.calcs[0]

    @property
    def rhs_calc(self) -> Calculation:
        assert self.goal.is_compare()
        return self.calcs[1]

    @property
    def arg_calc(self) -> Calculation:
        assert expr.is_fun(self.goal)
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
        raise NotImplementedError

    def export(self):
        return {
            "type": "CalculationProof",
            "goal": str(self.goal),
            "latex_goal": latex.convert_expr(self.goal),
            "finished": self.is_finished(),
            "calcs": [calc.export() for calc in self.calcs]
        }

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


def get_comp_file(p):
    while not isinstance(p, CompFile):
        p = p.parent
    return p


def conds_subst(conds: Conditions, var: str, e: Expr):
    res = Conditions()
    for cond in conds.data:
        res.add_condition(parser.parse_expr(str(cond.subst(var, e))))
    return res


class InductionProof(StateItem):
    """Proof for an equation by induction on natural numbers.

    This breaks the equation goal into two goals, corresponding to the
    base case and inductive case.

    """

    def __init__(self, parent: Goal, ctx: Context, goal: Expr, induct_var: str,
                 *, start: Union[int, Expr] = 0):
        if not goal.is_equals():
            raise AssertionError("InductionProof: currently only support equation goals.")

        self.parent = parent
        self.goal = goal
        self.induct_var = induct_var
        self.ctx = ctx

        if isinstance(start, int):
            self.start = Const(start)
        elif isinstance(start, Expr):
            self.start = start
        else:
            raise NotImplementedError
        # Base case: n = start
        base_goal_ctx = Context(self.ctx)
        eq0 = normalize(goal.subst(induct_var, self.start), base_goal_ctx)
        self.base_case = Goal(self, base_goal_ctx, eq0)

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

    def export(self):
        return {
            "type": "InductionProof",
            "goal": str(self.goal),
            "latex_goal": latex.convert_expr(self.goal),
            "induct_var": self.induct_var,
            "base_case": self.base_case.export(),
            "induct_case": self.induct_case.export(),
            'start': str(self.start),
            "finished": self.is_finished()
        }

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

    def __init__(self, parent, ctx, goal: Expr, *, split_cond: Expr):
        self.parent = parent
        self.goal = goal
        self.ctx = ctx
        self.split_cond = split_cond
        self.split_type = ""
        self.cases: List[Goal] = []
        assert isinstance(parent, Goal)

        if split_cond.is_compare():
            self.split_type = "two-way"
            # Case 1:
            conds1 = Conditions()
            case1_ctx = self.ctx
            conds1.add_condition(split_cond)
            conds1.update(parent.conds)
            self.cases.append(Goal(self, case1_ctx, goal, conds=conds1))

            # Case 2:
            conds2 = Conditions()
            case2_ctx = self.ctx
            conds2.add_condition(expr.neg_expr(split_cond))
            conds2.update(parent.conds)
            self.cases.append(Goal(self, case2_ctx, goal, conds=conds2))

        else:
            self.split_type = "three-way"
            # Case 1:
            conds1 = Conditions()
            conds1.add_condition(expr.Op("<", split_cond, Const(0)))
            case1_ctx = self.ctx
            conds1.update(parent.conds)
            self.cases.append(Goal(self, case1_ctx, goal, conds=conds1))

            # Case 2:
            conds2 = Conditions()
            case2_ctx = self.ctx
            conds2.add_condition(expr.Op("=", split_cond, Const(0)))
            conds2.update(parent.conds)
            self.cases.append(Goal(self, case2_ctx, goal, conds=conds2))

            # Case 3:
            conds3 = Conditions()
            case3_ctx = self.ctx
            conds3.add_condition(expr.Op(">", split_cond, Const(0)))
            conds3.update(parent.conds)
            self.cases.append(Goal(self, case3_ctx, goal, conds=conds3))

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
            raise AssertionError

    def __str__(self):
        if self.is_finished():
            res = "Proof by cases (finished)\n"
        else:
            res = "Proof by cases\n"
        for i, case in enumerate(self.cases):
            res += "case%d: %s for %s\n" % (i + 1, case.goal, case.conds)
            res += str(case)
        return res

    def is_finished(self):
        for case in self.cases:
            if not case.is_finished():
                return False
        return True

    def export(self):
        return {
            "type": "CaseProof",
            "goal": str(self.goal),
            "latex_goal": latex.convert_expr(self.goal),
            "cases": [case.export() for case in self.cases],
            "split_cond": str(self.split_cond),
            "latex_split_cond": latex.convert_expr(self.split_cond),
            "finished": self.is_finished()
        }

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
    """

    def __init__(self, parent: StateItem, ctx: Context, goal: Expr, *, start: str):
        if not goal.is_equals():
            raise AssertionError("RewriteGoalProof: goal is not an equality.")
        self.parent = parent
        self.goal = goal
        self.ctx = ctx
        self.start = start
        start_goal = ctx.get_subgoal(start)
        if not start_goal:
            raise AssertionError("RewriteGoalProof: start %s not found" % start)
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

    def export(self):
        res = {
            "type": "RewriteGoalProof",
            "goal": str(self.goal),
            "latex_goal": latex.convert_expr(self.goal),
            "start": self.start,
            "finished": self.is_finished()
        }
        return res

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
    name - name of the file.

    """
    def __init__(self, ctx: Union[Context, str], name: str):
        if isinstance(ctx, str):
            self.ctx = Context()
            self.ctx.load_book(ctx, upto=name)
        else:
            self.ctx = ctx
        self.name: str = name
        self.content: List[StateItem] = []

    def __eq__(self, other):
        return isinstance(other, CompFile) and \
               self.name == other.name and self.content == other.content

    def __str__(self):
        res = "File %s\n" % self.name
        for st in self.content:
            res += str(st)
        return res

    def get_context(self, index: int = -1) -> Context:
        """Obtain the context up to the particular index (exclusive).

        If index = -1, return the context after processing all the content.

        """
        ctx = Context(self.ctx)
        for item in (self.content if index == -1 else self.content[:index]):
            if isinstance(item, FuncDef):
                ctx.add_definition(item.eq, item.conds)
                ctx.add_lemma(item.eq, item.conds)
            elif isinstance(item, Goal):
                ctx.add_lemma(item.goal, item.conds)
                ctx.extend_by_item(item.export_book())
        return ctx

    def add_definition(self, funcdef: Union[str, Expr], *, conds: List[Union[str, Expr]] = None) -> FuncDef:
        """Add a function definition.

        funcdef: statement of the definition.
        conds: list of conditions for the definition. This is ignored if input
               is already of type FuncDef.

        """
        if conds is not None:
            if isinstance(conds, Conditions):
                pass
            else:
                for i in range(len(conds)):
                    if isinstance(conds[i], str):
                        conds[i] = parser.parse_expr(conds[i])
        else:
            conds = []
        if isinstance(funcdef, str):
            funcdef = parser.parse_expr(funcdef)
        if isinstance(funcdef, Expr):
            if funcdef.is_equals():
                self.content.append(FuncDef(self, self.ctx, funcdef, Conditions(conds)))
            else:
                raise NotImplementedError
        else:
            raise NotImplementedError

        return self.content[-1]

    def add_calculation(self, calc: Union[str, Expr], *, conds: List[Union[str, Expr]] = None) -> Calculation:
        """Add a calculation."""
        ctx = self.get_context()
        if conds is not None:
            for i in range(len(conds)):
                if isinstance(conds[i], str):
                    conds[i] = parser.parse_expr(conds[i])
        else:
            conds = []
        conds = Conditions(conds)
        if isinstance(calc, str):
            self.content.append(Calculation(self, ctx, parser.parse_expr(calc), conds=conds))
        elif isinstance(calc, Expr):
            self.content.append(Calculation(self, ctx, calc, conds=conds))
        else:
            raise NotImplementedError
        return self.content[-1]

    def make_goal(self, goal: Union[str, Expr, Goal], *,
                  conds: Optional[List[Union[str, Expr]]] = None) -> Goal:
        if isinstance(goal, Goal):
            self.content.append(goal)
            return self.content[-1]
        # Parse goal statement
        if isinstance(goal, str):
            goal = parser.parse_expr(goal)
        assert isinstance(goal, Expr)

        # Parse conditions
        if conds is not None:
            for i in range(len(conds)):
                if isinstance(conds[i], str):
                    conds[i] = parser.parse_expr(conds[i])
        else:
            conds = []

        conds = Conditions(conds)
        ctx = self.get_context()
        return Goal(self, ctx, goal, conds=conds)

    def add_goal(self, goal: Union[str, Expr, Goal], *,
                 conds: Optional[List[Union[str, Expr]]] = None) -> Goal:
        """Add a goal.

        goal: statement of the goal.
        conds: list of conditions for the goal. This is ignored if input goal
               is already of type Goal.

        """
        self.content.append(self.make_goal(goal, conds=conds))
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
            elif isinstance(root, FuncDef) or isinstance(root, CalculationStep):
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

    def export(self):
        self.name = self.name
        return {
            "name": self.name,
            "content": [item.export() for item in self.content]
        }


def parse_rule(item, parent) -> Rule:
    if 'loc' in item:
        if item['loc'] == 'subterms':
            del item['loc']
            return rules.OnSubterm(parse_rule(item, parent))
        else:
            loc = item['loc']
            del item['loc']
            if loc == '' or loc == '.':
                return parse_rule(item, parent)
            else:
                return rules.OnLocation(parse_rule(item, parent), loc)
    elif item['name'] == 'ExpandDefinition':
        func_name = item['func_name']
        return rules.ExpandDefinition(func_name=func_name)
    elif item['name'] == 'FoldDefinition':
        func_name = item['func_name']
        return rules.FoldDefinition(func_name=func_name)
    elif item['name'] == 'DerivIntExchange':
        return rules.DerivIntExchange()
    elif item['name'] == 'Simplify':
        return rules.Simplify()
    elif item['name'] == 'ElimInfInterval':
        a = Const(0)
        if 'a' in item:
            a = parser.parse_expr(item['a'])
        return rules.ElimInfInterval(a)
    elif item['name'] == 'Substitution':
        var_name = item['var_name']
        var_subst = parser.parse_expr(item['var_subst'])
        return rules.Substitution(var_name, var_subst)
    elif item['name'] == 'SubstitutionInverse':
        var_name = item['var_name']
        old_var = item['old_var']
        var_subst = parser.parse_expr(item['var_subst'])
        return rules.SubstitutionInverse(var_name, old_var, var_subst)
    elif item['name'] == 'IntegrationByParts':
        u = parser.parse_expr(item['u'])
        v = parser.parse_expr(item['v'])
        return rules.IntegrationByParts(u, v)
    elif item['name'] == 'Equation':
        new_expr = parser.parse_expr(item['new_expr'])
        old_expr = parser.parse_expr(item['old_expr']) if ('old_expr' in item) else None
        return rules.Equation(old_expr, new_expr)
    elif item['name'] == 'ApplyEquation':
        eq = parser.parse_expr(item['eq'])
        if 'source' in item:
            source = parser.parse_expr(item['source'])
        else:
            source = None
        return rules.ApplyEquation(eq, source)
    elif item['name'] == 'ExpandPolynomial':
        return rules.ExpandPolynomial()
    elif item['name'] == 'SplitRegion':
        c = parser.parse_expr(item['c'])
        return rules.SplitRegion(c)
    elif item['name'] == 'IntegrateByEquation':
        lhs = parser.parse_expr(item['lhs'])
        return rules.IntegrateByEquation(lhs)
    elif item['name'] == 'LHopital':
        return rules.LHopital()
    elif item['name'] == 'ApplyInductHyp':
        return rules.ApplyInductHyp()
    elif item['name'] == 'DerivativeSimplify':
        return rules.DerivativeSimplify()
    elif item['name'] == 'IntegrateBothSide':
        return rules.IntegralEquation()
    elif item['name'] == 'LimitEquation':
        var = item['var']
        lim = parser.parse_expr(item['lim'])
        return rules.LimitEquation(var, lim)
    elif item['name'] == 'IntSumExchange':
        return rules.IntSumExchange()
    elif item['name'] == 'DerivEquation':
        var = item['var']
        return rules.DerivEquation(var)
    elif item['name'] == 'SolveEquation':
        solve_for = parser.parse_expr(item['solve_for'])
        return rules.SolveEquation(solve_for)
    elif item['name'] == 'VarSubsOfEquation':
        subst = item['subst']
        return rules.VarSubsOfEquation(subst)
    elif item['name'] == 'ApplyIdentity':
        source = parser.parse_expr(item['source'])
        target = parser.parse_expr(item['target'])
        return rules.ApplyIdentity(source, target)
    elif item['name'] == 'IndefiniteIntegralIdentity':
        return rules.IndefiniteIntegralIdentity()
    elif item['name'] == 'DefiniteIntegralIdentity':
        return rules.DefiniteIntegralIdentity()
    elif item['name'] == 'SeriesExpansionIdentity':
        index_var = item['index_var']
        old_expr = None
        if 'old_expr' in item:
            old_expr = parser.parse_expr(item['old_expr'])
        return rules.SeriesExpansionIdentity(old_expr=old_expr, index_var=index_var)
    elif item['name'] == 'SeriesEvaluationIdentity':
        return rules.SeriesEvaluationIdentity()
    elif item['name'] == 'ReplaceSubstitution':
        return rules.ReplaceSubstitution()
    elif item['name'] == 'ChangeSummationIndex':
        e = parser.parse_expr(item['new_lower'])
        return rules.ChangeSummationIndex(e)
    elif item['name'] == 'SummationEquation':
        idx_v = item['index_var']
        lower = parser.parse_expr(item['lower'])
        upper = parser.parse_expr(item['upper'])
        return rules.SummationEquation(idx_v, lower, upper)
    elif item['name'] == 'FunEquation':
        func_name = item['func_name']
        return rules.FunEquation(func_name)
    elif item['name'] == 'PartialFractionDecomposition':
        return rules.PartialFractionDecomposition()
    else:
        print(item['name'], flush=True)
        raise NotImplementedError


def parse_step(calc: Calculation, item, id: int) -> CalculationStep:
    assert item['type'] == 'CalculationStep'
    assert isinstance(calc, Calculation), "it should belong to a Calculation"
    rule = parse_rule(item['rule'], calc)
    ctx = Context(calc.ctx)
    cur_e = calc.start
    for step in calc.steps:
        ctx = step.rule.update_context(cur_e, ctx)
        cur_e = step.res
    new_e = rule.eval(calc.last_expr, ctx)
    step = CalculationStep(calc, rule, new_e, id)
    return step


def parse_conds(item) -> Conditions:
    res = Conditions()
    if 'conds' in item:
        for subitem in item['conds']:
            res.add_condition(parser.parse_expr(subitem['cond']))
    return res


def parse_calculatioin(parent, item) -> Calculation:
    assert item['type'] == 'Calculation'
    file = get_comp_file(parent)
    cur_id = len(file.content)
    if isinstance(parent, CompFile):
        ctx = parent.get_context()
        start = parser.parse_expr(item['start'])
        conds = parse_conds(item)
        res = Calculation(parent, ctx, start, conds=conds)
    elif isinstance(parent, CalculationProof):
        ctx = file.get_context(cur_id)
        goal = parent.parent
        assert isinstance(goal, Goal), "this calculation should belong to a Goal"
        if len(goal.ctx.induct_hyps) > 0:
            ctx.add_induct_hyp(goal.ctx.induct_hyps[0].expr)
        start = parser.parse_expr(item['start'])
        conds = goal.conds
        res = Calculation(parent, ctx, start, conds=conds)
    elif isinstance(parent, RewriteGoalProof):
        ctx = file.get_context(cur_id)
        start = item['start']
        conds = parse_conds(item)
        res = Calculation(parent, ctx, start, conds=conds)
    else:
        raise NotImplementedError
    for i, step in enumerate(item['steps']):
        e = res.last_expr
        res.add_step(parse_step(res, step, i))
    return res


def parse_goal(parent, item, ih=None) -> Goal:
    assert item['type'] == 'Goal'
    file = get_comp_file(parent)
    ctx = file.get_context(len(file.content) - 1)
    if ih is not None:
        ctx.add_induct_hyp(ih)
    goal = parser.parse_expr(item['goal'])
    conds = parse_conds(item)
    res = Goal(parent, ctx, goal, conds=conds)
    if 'subgoals' in item:
        res.subgoals = []
        for subgoal in item['subgoals']:
            res.subgoals.append((subgoal['name'], parse_goal(res, subgoal['goal'])))
    if 'proof' in item:
        res.proof = parse_item(res, item['proof'])
    if 'wellformed' in item:
        res.wellformed = item['wellformed']
        if not res.wellformed and 'obligations' in item:
            res.proof_obligations = list()
            for obligation in item['obligations']:
                branches = list()
                for b in obligation['branches']:
                    tmp = list()
                    for e in b['exprs']:
                        tmp.append(parser.parse_expr(e))
                    branches.append(rules.ProofObligationBranch(tmp))
                c = parse_conds(obligation)
                res.proof_obligations.append(rules.ProofObligation(branches, c))
    return res


def parse_item(parent, item) -> StateItem:
    file = get_comp_file(parent)
    if item['type'] == 'FuncDef':
        conds = parse_conds(item)
        eq = parser.parse_expr(item['eq'])
        return FuncDef(parent, parent.ctx, eq, conds=conds)
    elif item['type'] == 'CalculationProof':
        goal = parser.parse_expr(item['goal'])
        res = CalculationProof(parent, parent.ctx, goal)
        for i, calc_item in enumerate(item['calcs']):
            res.calcs[i] = parse_calculatioin(res, calc_item)
        return res
    elif item['type'] == 'Goal':
        return parse_goal(parent, item)
    elif item['type'] == 'Calculation':
        return parse_calculatioin(parent, item)
    elif item['type'] == 'InductionProof':
        file = get_comp_file(parent)
        assert isinstance(parent, Goal)
        goal = parser.parse_expr(item['goal'])
        induct_var = item['induct_var']
        res = InductionProof(parent, goal, induct_var)
        res.start = parser.parse_expr(item['start'])
        res.base_case = parse_goal(res, item['base_case'])
        res.induct_case = parse_goal(res, item['induct_case'], ih=goal)
        return res
    elif item['type'] == 'CaseProof':
        ctx = parent.ctx
        goal = parser.parse_expr(item['goal'])
        split_cond = parser.parse_expr(item['split_cond'])
        res = CaseProof(parent, goal, split_cond=split_cond)
        assert len(res.cases) == len(item['cases'])
        for i, case in enumerate(item['cases']):
            res.cases[i] = parse_goal(res, case)
        return res
    elif item['type'] == 'RewriteGoalProof':
        goal = parser.parse_expr(item['goal'])
        begin_goal=parser.parse_expr(item['start'])
        res = RewriteGoalProof(parent, goal, begin=begin_goal)
        res.begin = parse_calculatioin(res, item['start'])
        return res
    else:
        print(item['type'])
        raise NotImplementedError


def get_next_step_label(step: Union[Calculation, CalculationStep], label: Label) -> Label:
    if isinstance(step, Calculation):
        return Label(label.data + [0])
    elif isinstance(step, CalculationStep):
        return Label(label.data[:-1] + [label.data[-1] + 1])
    elif isinstance(step, RewriteGoalProof):
        return Label(label.data + [0])
    else:
        raise NotImplementedError
