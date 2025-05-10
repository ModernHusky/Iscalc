"""Context of integral calculations"""

from typing import Iterable, Optional, Callable
import os

from integral import expr
from integral.expr import Expr, Eq, Op, Const, expr_to_pattern
from integral.conditions import Conditions
from integral import action


class Definition:
    """Introduce a new function definition.
    
    Attributes
    ----------
    symbol: str
        symbol to be defined
    args: tuple[str]
        list of names of arguments
    conds: Conditions
        conditions on the symbol
    define_eq: Optional[Expr]
        equality of definition, does not contain patterns

    """
    def __init__(self, symbol: str, args: Iterable[str], conds: Conditions,
                 define_eq: Optional[Expr] = None):
        self.symbol = symbol
        self.args = tuple(args)
        if len(self.args) != len(set(self.args)):
            raise AssertionError("Definition: arguments must be distinct")
        self.conds = conds
        if define_eq is None:
            self.define_eq = None
        else:
            self.define_eq = define_eq

    def __str__(self):
        res = "define "
        if not self.define_eq:
            res += f"{self.symbol}({', '.join(self.args)})"
        else:
            res += str(self.define_eq)
        if self.conds:
            res += " for " + str(self.conds)
        return res
    
    def __repr__(self):
        return f"Definition({repr(self.symbol)}, {repr(self.args)}, [{str(self.conds)}], {str(self.define_eq)})"

    def __eq__(self, other):
        return isinstance(other, Definition) and self.symbol == other.symbol and \
            self.conds == other.conds and self.define_eq == other.define_eq


class Identity:
    """Introduces a theorem.
    
    Attributes
    ----------
    expr: Expr
        statement of the theorem
    conds: Conditions
        conditions on the theorem
    attributes: tuple[str]
        list of attributes of the theorem

            """
    def __init__(self, expr: Expr, *, conds: Conditions = None,
                 attrs: Iterable[str] = tuple()):
        self.expr = expr
        if conds is None:
            conds = Conditions()
        self.conds = conds
        self.attrs = attrs

    def __eq__(self, other: "Identity"):
        return isinstance(other, Identity) and self.expr == other.expr and \
            self.conds == other.conds and self.attrs == other.attrs

    @property
    def lhs(self):
        return self.expr.lhs

    @property
    def rhs(self):
        return self.expr.rhs

    def __str__(self):
        if self.conds:
            return f"{self.expr} for {self.conds}"
        else:
            return str(self.expr)

    def __repr__(self):
        return str(self)


def symb_identity(e: Expr, conds: Conditions, attrs: Iterable[str] = tuple()) -> Identity:
    """Obtain identity with symbolic variables."""
    if expr.is_compare(e):
        symb_e = expr.Op(e.op, expr_to_pattern(e.args[0]), expr_to_pattern(e.args[1]))
    else:
        symb_e = expr_to_pattern(e)
    symb_conds = [expr_to_pattern(cond) for cond in conds.data]
    return Identity(symb_e, conds=Conditions(symb_conds), attrs=attrs)

class Context:
    """Maintains the current context of calculation.

    Information kept in context include the following:

    - List of function definitions

    - List of existing identities (for indefinite integrals, definite integrals,
      trigonometric identities, etc).

    - Assumptions for the current computation.

    - List of variable substitutions.

    """
    def __init__(self, parent: Optional["Context"] = None):

        # Parent context
        if parent is not None:
            assert isinstance(parent, Context)
        self.parent = parent

        # List of definitions
        self.definitions: dict[str, Definition] = dict()

        # List of indefinite integral identities
        self.indefinite_integrals: list[Identity] = list()

        # List of definite integral identities
        self.definite_integrals: list[Identity] = list()

        # List of series expansions
        self.series_expansions: list[Identity] = list()

        # List of series evaluations
        self.series_evaluations: list[Identity] = list()

        # List of other identities (trigonometric, etc)
        self.other_identities: list[Identity] = list()

        # List of simplification rules
        self.simp_identities: list[Identity] = list()

        # List of inequalities
        self.inequalities: list[Identity] = list()

        # Inductive hypothesis
        self.induct_hyps: list[Identity] = list()

        # List of variables
        self.vars: list[str] = list()

        # List of assumptions
        self.conds: Conditions = Conditions()

        # List of substitutions
        self.substs: list[tuple[str, Expr]] = list()

        # List of subgoals
        self.subgoals: dict[str, Identity] = dict()

        # List of identities of summation/product split
        self.split_identities: list[Identity] = list()

    def __str__(self):
        res = ""
        res += "Definitions\n"
        for _, definition in self.get_definitions().items():
            res += str(definition) + "\n"
        res += "Indefinite integrals\n"
        for identity in self.get_indefinite_integrals():
            res += str(identity) + "\n"
        res += "Definite integrals\n"
        for identity in self.get_definite_integrals():
            res += str(identity) + "\n"
        res += "Series expansions\n"
        for identity in self.get_series_expansions():
            res += str(identity) + "\n"
        res += "Series evaluations\n"
        for identity in self.get_series_evaluations():
            res += str(identity) + "\n"
        res += "Other identities\n"
        for identity in self.get_other_identities():
            res += str(identity) + "\n"
        res += "Simplification rules\n"
        for identity in self.get_simp_identities():
            res += str(identity) + "\n"
        res += "Inequalities\n"
        for identity in self.get_inequalities():
            res += str(identity) + "\n"
        res += "Inductive hypothesis\n"
        for identity in self.get_induct_hyps():
            res += str(identity) + "\n"
        res += "Conditions\n"
        for cond in self.get_conds().data:
            res += str(cond) + "\n"
        res += "Substitutions\n"
        for var_name, var_subst in self.get_substs():
            res += "%s: %s" % (var_name, var_subst) + "\n"

        return res

    def get_definitions(self) -> dict[str, Definition]:
        """Obtain all definitions in the context."""        
        res = self.parent.get_definitions() if self.parent is not None else dict()
        for name, definition in self.definitions.items():
            res[name] = definition
        return res

    def get_split_identities(self) -> list[Identity]:
        res = self.parent.get_split_identities() if self.parent is not None else []
        res.extend(self.split_identities)
        return res

    def get_indefinite_integrals(self) -> list[Identity]:
        res = self.parent.get_indefinite_integrals() if self.parent is not None else []
        res.extend(self.indefinite_integrals)
        return res

    def get_definite_integrals(self) -> list[Identity]:
        res = self.parent.get_definite_integrals() if self.parent is not None else []
        res.extend(self.definite_integrals)
        return res

    def get_series_expansions(self) -> list[Identity]:
        res = self.parent.get_series_expansions() if self.parent is not None else []
        res.extend(self.series_expansions)
        return res

    def get_series_evaluations(self) -> list[Identity]:
        res = self.parent.get_series_evaluations() if self.parent is not None else []
        res.extend(self.series_evaluations)
        return res

    def get_other_identities(self) -> list[Identity]:
        res = self.parent.get_other_identities() if self.parent is not None else []
        res.extend(self.other_identities)
        return res

    def get_simp_identities(self) -> list[Identity]:
        res = self.parent.get_simp_identities() if self.parent is not None else []
        res.extend(self.simp_identities)
        return res

    def get_inequalities(self) -> list[Identity]:
        res = self.parent.get_inequalities() if self.parent is not None else []
        res.extend(self.inequalities)
        return res

    def get_induct_hyps(self) -> list[Identity]:
        res = self.parent.get_induct_hyps() if self.parent is not None else []
        res.extend(self.induct_hyps)
        return res

    def get_vars(self) -> list[str]:
        res = self.parent.get_vars() if self.parent is not None else []
        res.extend(self.vars)
        return res

    def get_conds(self) -> Conditions:
        res = self.parent.get_conds() if self.parent is not None else Conditions()
        for cond in self.conds.data:
            res.add_condition(cond)
        return res

    def get_subgoal(self, name: str) -> Optional[Identity]:
        """Obtain subgoal with the given name."""
        res = self.parent.get_subgoal(name) if self.parent is not None else None
        if res:
            return res
        if name in self.subgoals:
            return self.subgoals[name]
        else:
            return None
        
    def get_all_subgoals(self) -> dict[str, Identity]:
        """Obtain all subgoals in the context."""
        res = self.parent.get_all_subgoals() if self.parent is not None else dict()
        for name, identity in self.subgoals.items():
            res[name] = identity
        return res

    def get_eq_conds(self) -> Conditions:
        parent_conds = self.parent.get_conds() if self.parent is not None else Conditions()
        res = Conditions()
        for cond in parent_conds.data:
            if cond.is_equals():
                res.add_condition(cond)
        for cond in self.conds.data:
            if cond.is_equals():
                res.add_condition(cond)
        return res

    def get_substs(self) -> list[str, Expr]:
        res = self.parent.get_substs() if self.parent is not None else list()
        res.extend(self.substs)
        return res

    def add_axiom_definition(self, e: Expr, conds: Conditions):
        """Add axiomatic definition."""
        if isinstance(e, expr.Var):
            if conds:
                raise AssertionError(f"axiom_define: no condition allowed for {e}")
            self.definitions[e.name] = Definition(e.name, tuple(), Conditions())
        elif isinstance(e, expr.Fun):
            if not all(expr.is_var(arg) for arg in e.args):
                raise AssertionError(f"axiom_define: all arguments of {e} must be variables")
            args = list(arg.name for arg in e.args)
            for cond in conds.data:
                if not cond.get_vars().issubset(args):
                    raise AssertionError(f"axiom_define: condition {cond} contains extra variable")
            self.definitions[e.func_name] = Definition(e.func_name, args, conds)
        else:
            raise AssertionError(f"axiom define: unrecognized expression {e}")

    def add_definition(self, eq: Expr, conds: Conditions):
        """Add definition."""
        if not eq.is_equals():
            raise AssertionError(f"define: {eq}")

        e = eq.lhs
        if isinstance(e, expr.Var):
            if conds:
                raise AssertionError(f"define: no condition allowed for {e}")
            self.definitions[e.name] = Definition(e.name, tuple(), Conditions(), eq)
        elif isinstance(e, expr.Fun):
            if not all(expr.is_var(arg) for arg in e.args):
                raise AssertionError(f"axiom_define: all arguments of {e} must be variables")
            args = list(arg.name for arg in e.args)
            for cond in conds.data:
                if not cond.get_vars().issubset(args):
                    raise AssertionError(f"define: condition {cond} for {e} contains extra variable")
            self.definitions[e.func_name] = Definition(e.func_name, args, conds, eq)
        else:
            raise AssertionError(f"axiom define: unrecognized expression {e}")

    def add_indefinite_integral(self, eq: Expr, conds: Conditions, attrs: Iterable[str]):
        assert isinstance(conds, Conditions)
        if not (eq.is_equals() and expr.is_indefinite_integral(eq.lhs)):
            raise AssertionError(f"add_indefinite_integral: {eq}")

        self.indefinite_integrals.append(symb_identity(eq, conds, attrs))

    def add_definite_integral(self, eq: Expr, conds: Conditions, attrs: Iterable[str]):
        assert isinstance(conds, Conditions)
        if not (eq.is_equals() and expr.is_integral(eq.lhs)):
            raise AssertionError(f"add_definite_integral: {eq}")

        self.definite_integrals.append(symb_identity(eq, conds, attrs))

    def add_series_expansion(self, eq: Expr, conds: Conditions):
        assert isinstance(conds, Conditions)
        if not (eq.is_equals() and not expr.is_summation(eq.lhs) and expr.is_summation(eq.rhs)):
            raise AssertionError(f"add_series_expansion: {eq}")

        self.series_expansions.append(symb_identity(eq, conds))

    def add_series_evaluation(self, eq: Expr, conds: Conditions):
        assert isinstance(conds, Conditions)
        if not (eq.is_equals() and expr.is_summation(eq.lhs) and not expr.is_summation(eq.rhs)):
            raise AssertionError(f"add_series_expansion: {eq}")

        self.series_evaluations.append(symb_identity(eq, conds))

    def add_other_identities(self, eq: Expr, conds: Conditions, attrs: Iterable[str]):
        assert isinstance(conds, Conditions)
        if not eq.is_equals():
            raise AssertionError(f"add_other_identities: {eq}")

        self.other_identities.append(symb_identity(eq, conds, attrs))
        attrs = tuple(attrs)
        if attrs and 'bidirectional' in attrs:
            rev_eq = Eq(eq.rhs, eq.lhs)
            self.other_identities.append(symb_identity(rev_eq, conds, attrs))

    def add_simp_identity(self, eq: Expr, conds: Conditions):
        assert isinstance(conds, Conditions)
        if not eq.is_equals():
            raise AssertionError(f"add_simp_identity: {eq}")

        self.simp_identities.append(symb_identity(eq, conds))

    def add_inequality(self, e: Expr, conds: Conditions):
        assert isinstance(conds, Conditions)
        if not expr.is_compare(e):
            raise AssertionError(f"add_inequality: {e}")

        self.inequalities.append(symb_identity(e, conds))

    def add_split_identities(self, e: Expr, conds: Conditions):            
        assert isinstance(conds, Conditions)
        self.split_identities.append(symb_identity(e, conds))

    def add_induct_hyp(self, e: Expr):
        # Note: no conversion to symbols for inductive hypothesis
        self.induct_hyps.append(Identity(e))

    def add_condition(self, cond: Expr):
        if cond not in self.conds.data:
            self.conds.add_condition(cond)

    def extend_vars(self, vars: Iterable[str]):
        cur_vars = self.get_vars()
        for var in vars:
            if var not in cur_vars:
                self.vars.append(var)

    def extend_condition(self, conds: Conditions):
        for cond in conds.data:
            self.add_condition(cond)

    def add_subst(self, var: str, expr: Expr):
        self.substs.append((var, expr))

    def add_subgoal(self, name: str, identity: Identity):
        self.subgoals[name] = identity

    def load_book(self, book_name: str):
        """Load the book with the given name.
        
        This function recursively loads imported books.
        
        """
        from integral import parser

        assert isinstance(book_name, str)
        root_dir = os.path.dirname(os.path.dirname(__file__))

        thy_filename = os.path.join(root_dir, 'theories', book_name + '.thy')
        with open(thy_filename, 'r', encoding='utf-8') as f:
            content = f.read()
        actions = [s for s in content.split('\n') if s.strip()]
        for act in actions:
            if act.lstrip().startswith("#") or act.lstrip().startswith("//"):
                # title of comment
                continue
            a = parser.parse_action(act)
            if isinstance(a, action.ImportsAction):
                for book_name in a.theories:
                    self.load_book(book_name)
            elif isinstance(a, action.DefineAction):
                self.add_definition(a.expr, conds=a.conditions)
            elif isinstance(a, action.AxiomDefineAction):
                self.add_axiom_definition(a.expr, conds=a.conditions)
            elif isinstance(a, (action.AxiomAction, action.ProveAction)):
                if a.expr.is_equals() and expr.is_indefinite_integral(a.expr.lhs):
                    self.add_indefinite_integral(a.expr, a.conditions, a.attrs)
                elif a.expr.is_equals() and expr.is_integral(a.expr.lhs):
                    self.add_definite_integral(a.expr, a.conditions, a.attrs)
                elif a.expr.is_equals() and not expr.is_summation(a.expr.lhs) and expr.is_summation(a.expr.rhs):
                    self.add_series_expansion(a.expr, a.conditions)
                elif a.expr.is_equals() and expr.is_summation(a.expr.lhs) and not expr.is_summation(a.expr.rhs):
                    self.add_series_evaluation(a.expr, a.conditions)
                elif 'simp' in a.attrs:
                    self.add_simp_identity(a.expr, a.conditions)
                    self.add_other_identities(a.expr, a.conditions, a.attrs)
                else:
                    self.add_other_identities(a.expr, a.conditions, a.attrs)

    def check_condition(self, e: Expr) -> bool:
        """Check the given condition under the extra conditions"""
        from integral import condprover
        return condprover.check_condition(e, self)

    def is_positive(self, e: Expr) -> bool:
        return self.check_condition(Op(">", e, Const(0)))
    
    def is_negative(self, e: Expr) -> bool:
        return self.check_condition(Op("<", e, Const(0)))

    def is_nonzero(self, e: Expr) -> bool:
        return self.check_condition(Op("!=", e, Const(0)))

    def is_not_negative(self, e: Expr) -> bool:
        return self.check_condition(Op(">=", e, Const(0)))

    def is_not_positive(self, e: Expr) -> bool:
        return self.check_condition(Op("<=", e, Const(0)))

    def is_greater(self, e1: Expr, e2: Expr) -> bool:
        return self.check_condition(Op(">", e1, e2))

    def is_less(self, e1: Expr, e2: Expr) -> bool:
        return self.check_condition(Op("<", e1, e2))

    def is_greater_eq(self, e1: Expr, e2: Expr) -> bool:
        return self.check_condition(Op(">=", e1, e2))

    def is_less_eq(self, e1: Expr, e2: Expr) -> bool:
        return self.check_condition(Op("<=", e1, e2))

    def is_not_equal(self, e1: Expr, e2: Expr) -> bool:
        return self.check_condition(Op("!=", e1, e2))


def body_conds(e: Expr, ctx: Context) -> Context:
    """Return the conditions in the body."""
    ctx2 = Context(ctx)
    if expr.is_integral(e):
        ctx2.add_condition(expr.isReal(expr.Var(e.var)))
        if e.lower != expr.NEG_INF:
            ctx2.add_condition(Op(">", expr.Var(e.var), e.lower))
        if e.upper != expr.POS_INF:
            ctx2.add_condition(Op("<", expr.Var(e.var), e.upper))
    elif expr.is_indefinite_integral(e):
        ctx2.add_condition(expr.isReal(expr.Var(e.var)))
    elif expr.is_limit(e):
        ctx2.add_condition(expr.isReal(expr.Var(e.var)))
        if e.lim == expr.POS_INF:
            ctx2.add_condition(expr.Op(">", expr.Var(e.var), Const(0)))
    elif expr.is_summation(e) or expr.is_product(e):
        ctx2.add_condition(expr.Op(">=", expr.Var(e.index_var), e.lower))
        if e.upper != expr.POS_INF:
            ctx2.add_condition(expr.Op("<=", expr.Var(e.index_var), e.upper))
            ctx2.add_condition(expr.Op(">=", e.upper - expr.Var(e.index_var), Const(0)))
        ctx2.add_condition(expr.Fun("isInt", expr.Var(e.index_var)))
    else:
        raise TypeError
    return ctx2

def apply_subterm(e: Expr, f: Callable[[Expr, Context], Expr], ctx: Context) -> Expr:
    def rec(e: Expr, ctx: Context):
        if expr.is_var(e) or expr.is_const(e) or expr.is_inf(e) or expr.is_skolem_func(e):
            return f(e, ctx)
        elif expr.is_op(e):
            args = [rec(arg, ctx) for arg in e.args]
            return f(expr.Op(e.op, *args), ctx)
        elif expr.is_fun(e):
            args = [rec(arg, ctx) for arg in e.args]
            return f(expr.Fun(e.func_name, *args), ctx)
        elif expr.is_deriv(e):
            return f(expr.Deriv(e.var, rec(e.body, ctx)), ctx)
        elif expr.is_integral(e):
            lower = rec(e.lower, ctx)
            upper = rec(e.upper, ctx)
            body = rec(e.body, body_conds(e, ctx))
            return f(expr.Integral(e.var, lower, upper, body), ctx)
        elif expr.is_evalat(e):
            lower = rec(e.lower, ctx)
            upper = rec(e.upper, ctx)
            body = rec(e.body, ctx)
            return f(expr.EvalAt(e.var, lower, upper, body), ctx)
        elif expr.is_limit(e):
            return f(expr.Limit(e.var, rec(e.lim, ctx), rec(e.body, body_conds(e, ctx))), ctx)
        elif expr.is_indefinite_integral(e):
            body = rec(e.body, body_conds(e, ctx))
            return f(expr.IndefiniteIntegral(e.var, body, e.skolem_args), ctx)
        elif expr.is_summation(e):
            lower = rec(e.lower, ctx)
            upper = rec(e.upper, ctx)
            body = rec(e.body, body_conds(e, ctx))
            return f(expr.Summation(e.index_var, lower, upper, body), ctx)
        elif expr.is_product(e):
            lower = rec(e.lower, ctx)
            upper = rec(e.upper, ctx)
            body = rec(e.body, body_conds(e, ctx))
            return f(expr.Product(e.index_var, lower, upper, body), ctx)
        elif expr.is_symbol(e):
            return e
        else:
            raise NotImplementedError
    return rec(e, ctx)
