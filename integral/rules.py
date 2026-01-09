"""Rules for integration."""
import re
import math
from decimal import Decimal
from fractions import Fraction
from typing import Optional, Dict, Tuple, Union, List, Set
import functools
import operator

from integral import expr, context
from integral.expr import Var, Const, Fun, EvalAt, Op, Integral, Symbol, Expr, \
    OP, CONST, VAR, sin, cos, FUN, decompose_expr_factor, \
    Deriv, Inf, Limit, NEG_INF, POS_INF, IndefiniteIntegral, Summation, SUMMATION, SkolemFunc, decompose_expr_factor2, is_const, \
    Fraction, CINTPath, CIntegral
from integral import parser
from integral.solve import solve_equation, solve_for_term
from integral import latex
from integral import limits
from integral import norm
from integral.context import Context, apply_subterm, body_conds
from integral import poly
from integral.poly import from_poly, to_poly, normalize
from integral.conditions import Conditions
from integral import sympywrapper
from integral import utils


class RuleException(expr.IscalcException):
    """Exception raised when applying some calculation rule."""
    def __init__(self, rule_name: str, msg: str):
        self.rule_name = rule_name
        self.msg = msg

    def __str__(self):
        return "%s: %s" % (self.rule_name, self.msg)

    def to_json(self) -> dict:
        return {
            "class": "RuleException",
            "rule_name": self.rule_name,
            "msg": self.msg
        }
    @staticmethod
    def from_json(data: dict):
        return RuleException(data["rule_name"], data["msg"])


def deriv(var: str, e: Expr, ctx: Context) -> Expr:
    """Compute the derivative of e with respect to variable
    name var.

    """

    def normal(x):
        return normalize(x, ctx)

    def rec(e: Expr):
        if var not in e.get_vars():
            return Const(0)
        elif expr.is_var(e):
            if e.name == var:
                # dx. x = 1
                return Const(1)
            else:
                # dx. y = 0
                return Const(0)
        elif expr.is_const(e):
            # dx. c = 0
            return Const(0)
        elif expr.is_op(e):
            if e.op == "+":
                x, y = e.args
                return normal(rec(x) + rec(y))
            elif e.op == "-" and len(e.args) == 2:
                x, y = e.args
                return normal(rec(x) - rec(y))
            elif e.op == "-" and len(e.args) == 1:
                x, = e.args
                return normal(-(rec(x)))
            elif e.op == "*":
                x, y = e.args
                if not x.contains_var(var):
                    return normal(x * rec(y))
                elif not y.contains_var(var):
                    return normal(rec(x) * y)
                else:
                    return normal(x * rec(y) + rec(x) * y)
            elif e.op == "/":
                x, y = e.args
                if not y.contains_var(var):
                    # x / c case:
                    return normal(rec(x) / y)
                elif not x.contains_var(var) and expr.is_power(y):
                    # c / (y0 ^ y1): rewrite to c * y0 ^ (-y1)
                    return rec(x * (y.args[0] ^ (-y.args[1])))
                else:
                    # general case
                    return normal((rec(x) * y - x * rec(y)) / (y ^ Const(2)))
            elif e.op == "^":
                x, y = e.args
                if expr.is_const(y):
                    return normal(y * (x ^ Const(y.val - 1)) * rec(x))
                elif var not in y.get_vars():
                    return normal(y * (x ^ (y - 1)) * rec(x))
                else:
                    return normal(e * rec(y * expr.log(x)))

            else:
                raise NotImplementedError
        elif expr.is_fun(e):
            if e.func_name == "sin":
                x, = e.args
                return normal(cos(x) * rec(x))
            elif e.func_name == "cos":
                x, = e.args
                return normal(-(sin(x) * rec(x)))
            elif e.func_name == "tan":
                x, = e.args
                return normal((expr.sec(x) ^ Const(2)) * rec(x))
            elif e.func_name == "sec":
                x, = e.args
                return normal(expr.sec(x) * expr.tan(x) * rec(x))
            elif e.func_name == "csc":
                x, = e.args
                return normal(-expr.csc(x) * expr.cot(x) * rec(x))
            elif e.func_name == "cot":
                x, = e.args
                return normal(-(expr.csc(x) ^ Const(2)) * rec(x))
            elif e.func_name == "cot":
                x, = e.args
                return normal(-(sin(x) ^ Const(-2)) * rec(x))
            elif e.func_name == "log":
                x, = e.args
                return normal(rec(x) / x)
            elif e.func_name == "exp":
                x, = e.args
                return normal(expr.exp(x) * rec(x))
            elif e.func_name == "pi":
                return Const(0)
            elif e.func_name == "sqrt":
                if e.args[0].ty == CONST:
                    return Const(0)
                else:
                    return normal(rec(e.args[0] ^ Const(Fraction(1 / 2))))
            elif e.func_name == "arctan":
                x, = e.args
                return normal(rec(x) / (Const(1) + (x ^ Const(2))))
            elif e.func_name == "arcsin":
                x, = e.args
                return normal(rec(x) / expr.sqrt(Const(1) - (x ^ Const(2))))
            elif e.func_name == "arccos":
                x, = e.args
                return normal(-(rec(x) / expr.sqrt(Const(1) - (x ^ Const(2)))))
            elif e.func_name == "arccot":
                x, = e.args
                return normal(-rec(x)) / (Const(1) + x ^ Const(2))
            elif e.func_name == "arcsec":
                x, = e.args
                return normal(rec(x) / (expr.abs(x) * expr.sqrt(x ^ Const(2) - Const(1))))
            elif e.func_name == "arccsc":
                x, = e.args
                return normal(-(rec(x) / (expr.abs(x) * expr.sqrt(x ^ Const(2) - Const(1)))))
            elif e.func_name == "binom":
                # Arguments should be integers
                assert not e.contains_var(var), "deriv: binom applied to real variables"
                return Const(0)
            else:
                return Deriv(var, e)
        elif expr.is_integral(e):
            if e.lower.is_constant():
                return normal(Integral(e.var, e.lower, e.upper, rec(e.body))
                              + e.body.subst(e.var, e.upper) * rec(e.upper))
            return normal(Integral(e.var, e.lower, e.upper, rec(e.body))
                          + e.body.subst(e.var, e.upper) * rec(e.upper)
                          - e.body.subst(e.var, e.lower) * rec(e.lower))
        elif expr.is_limit(e):
            return Limit(e.var, e.lim, rec(e.body))
        elif expr.is_summation(e):
            return Summation(e.index_var, e.lower, e.upper, rec(e.body))
        elif expr.is_inf(e):
            return Const(0)
        else:
            raise NotImplementedError(f"{e}, {type(e)}")

    return rec(e)


class ProofObligationBranch:
    """Represents a single branch of proof obligation."""
    def __init__(self, exprs: list[Expr], flags: list[bool] = None):
        self.exprs = exprs  # satisfy all expressions
        if flags is None or len(flags) != len(exprs):
            self.need_to_be_satisfied = [True for i in range(len(exprs))]
        else:
            self.need_to_be_satisfied = flags

    def __str__(self):
        return ", ".join(f"{e} ({b})" for e, b in zip(self.exprs, self.need_to_be_satisfied))

    def export(self):
        res = {
            'exprs': [str(e) for e in self.exprs]
        }
        return res


class ProofObligation:
    """Represents a proof obligation to prove e using the conditions
    in conds.

    """

    def __init__(self, branches: list[ProofObligationBranch], conds: Conditions):
        # if any branch is satisfied then the proof obligation is carried out
        self.branches = branches
        self.conds = conds

    def __eq__(self, other):
        return isinstance(other, ProofObligation) and \
            self.branches == other.branches and self.conds == other.conds

    def __str__(self):
        res = ""
        for i, branch in enumerate(self.branches, 1):
            res += "Branch " + str(i) + ":\n"
            res += utils.indent(str(branch)) + '\n'
        res += "Conds: " + str(self.conds)
        return res

    def __repr__(self):
        return str(self)

    def __hash__(self):
        return hash(tuple(self.branches))

    def export(self):
        res = {
            'branches': [branch.export() for branch in self.branches],
            'conds': self.conds.export()
        }
        return res


def check_wellformed(e: Expr, ctx: Context) -> list[ProofObligation]:
    """Check whether an expression e is wellformed, and return
    a set of wellformed-ness conditions if otherwise.

    """
    obligations: List[ProofObligation] = list()

    def add_obligation(branches: Union[List[ProofObligationBranch], Expr], ctx: Context):
        if isinstance(branches, Expr):
            branches = [ProofObligationBranch([branches])]
        obligation = ProofObligation(branches, ctx.get_conds())
        if obligation not in obligations:
            obligations.append(obligation)

    def rec(e: Expr, ctx: Context):
        if expr.is_var(e) or expr.is_const(e):
            pass
        elif expr.is_op(e):
            for arg in e.args:
                rec(arg, ctx)
            if e.is_divides():
                # if the denominator has i, and var is real, then the expression is not 0
                if ctx.check_condition(expr.Fun("notReal", e.args[1])):
                    # 分母是notReal，由于notReal隐含非零，所以无需额外的!=0条件
                    pass
                else:
                    # 分母不是notReal，需要检查!=0条件
                    if ctx.check_condition(Op("!=", e.args[1], Const(0))):
                        pass
                    else:
                        add_obligation(Op("!=", e.args[1], Const(0)), ctx)
            if e.is_power():
                if ctx.check_condition(Op(">", e.args[0], Const(0))):
                    pass
                elif ctx.check_condition(Fun("isInt", e.args[1])) and ctx.check_condition(
                        Op(">=", e.args[1], Const(0))):
                    pass
                else:
                    add_obligation(Op(">", e.args[0], Const(0)), ctx)
                    add_obligation(Fun("isInt", e.args[1]), ctx)
                    add_obligation(Op(">=", e.args[1], Const(0)), ctx)
        elif expr.is_fun(e):
            for arg in e.args:
                rec(arg, ctx)
            if e.func_name == 'log':
                # log 的定义域: 复数域中 z != 0, 实数域中 x > 0
                arg = e.args[0]
                
                # 检查是否满足实数域条件: arg > 0
                if ctx.check_condition(Op(">", arg, Const(0))):
                    pass
                # 检查是否满足复数域条件: arg != 0
                # (notReal 通过推理规则会推出 != 0)
                elif ctx.check_condition(Op("!=", arg, Const(0))):
                    pass
                else:
                    # 无法证明参数有效，需要添加约束
                    # 提供两个分支：实数正数 或 复数非零
                    branch1 = ProofObligationBranch([Op(">", arg, Const(0))])
                    branch2 = ProofObligationBranch([Op("!=", arg, Const(0)), Fun("notReal", arg)])
                    add_obligation([branch1, branch2], ctx)
            if e.func_name == 'sqrt':
                # 如果是复数域的计算，允许负参数的平方根
                if ctx.check_condition(expr.Fun("notReal", e)) or any(ctx.check_condition(expr.Fun("notReal", cond)) for cond in ctx.get_conds().data):
                    pass
                elif ctx.check_condition(Op(">=", e.args[0], Const(0))):
                    pass
                else:
                    add_obligation(Op(">=", e.args[0], Const(0)), ctx)
            if e.func_name == 'gamma':
                f1 = ctx.check_condition(Op(">", e.args[0], Const(0)))
                f2 = ctx.check_condition(Op('<'), e.args[0], Const(0)) and \
                     ctx.check_condition(Fun("notInt", e.args[0]))
                if f1 or f2:
                    pass
                else:
                    branch1 = ProofObligationBranch([Op(">", e.args[0], Const(0))])
                    branch2 = ProofObligationBranch([Op("<", e.args[0], Const(0)), Fun("notInt", e.args[0])])
                    add_obligation([branch1, branch2], ctx)
            if e.func_name == "arccos" or e.func_name == "arcsin":
                f1 = ctx.check_condition(Op(">=", e.args[0], Const(-1)))
                f2 = ctx.check_condition(Op("<=", e.args[0], Const(1)))
                if f1 and f2:
                    pass
                else:
                    add_obligation(Op(">=", e.args[0], Const(-1)), ctx)
                    add_obligation(Op("<=", e.args[0], Const(1)), ctx)
            if e.func_name == 'tan':
                tmp = normalize(Const(2) * e.args[0] / expr.pi, ctx)
                f1 = ctx.check_condition(expr.isInt(tmp))
                f2 = ctx.check_condition(expr.isEven(tmp))

                if not f1 or f2:
                    pass
                else:
                    branch1 = ProofObligationBranch([expr.isInt(tmp)], [False])
                    branch2 = ProofObligationBranch([expr.isEven(tmp)])
                    add_obligation([branch1, branch2], ctx)
            if e.func_name == 'factorial':
                if not ctx.check_condition(expr.isInt(e.args[0])):
                    add_obligation(expr.isInt(e.args[0]), ctx)
            if e.func_name == 'binom':
                if not ctx.check_condition(expr.isInt(e.args[0])):
                    add_obligation(expr.isInt(e.args[0]), ctx)
                if not ctx.check_condition(expr.isInt(e.args[1])):
                    add_obligation(expr.isInt(e.args[1]), ctx)

            # TODO: add checks for other functions
        elif expr.is_integral(e):
            rec(e.body, body_conds(e, ctx))
        elif expr.is_indefinite_integral(e):
            rec(e.body, body_conds(e, ctx))
        elif expr.is_deriv(e):
            rec(e.body, ctx)
        elif expr.is_summation(e):
            rec(e.lower, ctx)
            rec(e.upper, ctx)
            rec(e.body, body_conds(e, ctx))
        else:
            pass

    rec(e, ctx)
    return obligations


def check_asymp_converge(asymp: limits.Asymptote) -> bool:
    if isinstance(asymp, limits.PolyLog):
        for n in asymp.order:
            if isinstance(n, (int, Fraction)) and n > 1:
                return True
            elif isinstance(n, Expr) and n.val > 1:
                return True
            elif isinstance(n, (int, Fraction)) and n == 1:
                continue  # check next orders
            elif isinstance(n, Expr) and n.val == 1:
                continue
            else:
                return False
        return False
    elif isinstance(asymp, limits.Exp):
        return True
    else:
        return False


def check_converge(e: Expr, ctx: Context) -> bool:
    """Check convergence of the sum or integral."""
    if expr.is_summation(e):
        lim = limits.limit_of_expr(e.body, e.index_var, ctx)
        if lim.e == Const(0) and check_asymp_converge(lim.asymp):
            return True
    elif expr.is_var(e) or e.is_constant():
        return True
    elif e.is_times() or e.is_plus():
        if check_converge(e.args[0], ctx) and check_converge(e.args[1], ctx):
            return True
    elif e.is_power():
        if check_converge(e.args[0], ctx) and check_converge(e.args[1], ctx) and \
                ctx.check_condition(Op(">", e.args[1], Const(0))):
            return True
    elif e.is_divides():
        if ctx.check_condition(Op("!=", e.args[1], Const(0))) and check_converge(e.args[0], ctx) and \
                check_converge(e.args[1], ctx):
            return True
    elif expr.is_fun(e):
        flag = True
        for arg in e.args:
            if not check_converge(arg, ctx):
                flag = False
                break
        if (flag):
            return True
    return False


class Rule:
    """
    Represents a rule for integration. It takes an integral
    to be evaluated (as an expression), then outputs a new
    expression that it is equal to.

    """
    # Name of the rule
    name: str

    def eval(self, e: Expr, ctx: Context) -> Expr:
        """Evaluation of the rule on the given expression. Returns
        a new expression.

        """
        raise NotImplementedError

    def export(self):
        """Returns the JSON representation of the rule."""
        raise NotImplementedError

    def update_context(self, e: Expr, ctx: Context) -> Context:
        """Produce the updated context after performing this rule."""
        return ctx


class Linearity(Rule):
    """Applies linearity rules:

    INT (a + b) = INT a + INT b,
    INT (c * a) = c * INT a      (where c is a constant).
    INT (c / a) = c * INT 1 / a  (where c is a constant).

    """

    def __init__(self):
        self.name = "Linearity"

    def __str__(self):
        return "linearity"

    def export(self):
        return {
            "name": self.name,
            "str": str(self)
        }

    def eval(self, e: Expr, ctx: Context) -> Expr:
        def prod(es):
            es = list(es)
            if len(es) == 0:
                return Const(1)
            else:
                return functools.reduce(operator.mul, es[1:], es[0])

        def depends_on_var(expr_obj: Expr, var: str) -> bool:
            """Check if expression depends on variable, handling derivatives correctly.
            
            For derivatives D x. f(x), the derivative itself depends on x even though
            x is a bound variable in the body. This is because the derivative is with
            respect to x, so it fundamentally depends on x.
            """
            if expr.is_deriv(expr_obj):
                # For derivatives D x. f(x), check if the derivative variable is var
                # If so, it depends on var regardless of the body
                if expr_obj.var == var:
                    return True
                # Otherwise check if the body contains var (excluding bound variable)
                else:
                    return expr_obj.body.contains_var(var)
            else:
                return expr_obj.contains_var(var)

        def rec(e: Expr):
            if expr.is_integral(e):
                if expr.is_plus(e.body):
                    return rec(expr.Integral(e.var, e.lower, e.upper, e.body.args[0])) + \
                           rec(expr.Integral(e.var, e.lower, e.upper, e.body.args[1]))
                elif expr.is_uminus(e.body):
                    return -rec(expr.Integral(e.var, e.lower, e.upper, e.body.args[0]))
                elif expr.is_minus(e.body):
                    return rec(expr.Integral(e.var, e.lower, e.upper, e.body.args[0])) - \
                           rec(expr.Integral(e.var, e.lower, e.upper, e.body.args[1]))
                elif expr.is_times(e.body) or expr.is_divides(e.body):
                    num_factors, denom_factors = decompose_expr_factor(e.body)
                    b = prod(f for f in num_factors if depends_on_var(f, e.var))
                    c = prod(f for f in num_factors if not depends_on_var(f, e.var))
                    denom_b = prod(f for f in denom_factors if depends_on_var(f, e.var))
                    denom_c = prod(f for f in denom_factors if not depends_on_var(f, e.var))
                    if denom_b != Const(1):
                        b = b / denom_b
                    if denom_c != Const(1):
                        c = c / denom_c
                    if c == expr.Const(1):
                        return Integral(e.var, e.lower, e.upper, b)
                    else:
                        return c * rec(Integral(e.var, e.lower, e.upper, b))
                elif e.body.is_constant() and e.body != Const(1):
                    return e.body * expr.Integral(e.var, e.lower, e.upper, Const(1))
                else:
                    return e
            elif expr.is_indefinite_integral(e):
                if expr.is_plus(e.body):
                    return rec(expr.IndefiniteIntegral(e.var, e.body.args[0], e.skolem_args)) + \
                        rec(expr.IndefiniteIntegral(e.var, e.body.args[1], e.skolem_args))
                elif expr.is_uminus(e.body):
                    return -rec(IndefiniteIntegral(e.var, e.body.args[0], e.skolem_args))
                elif expr.is_minus(e.body):
                    return rec(expr.IndefiniteIntegral(e.var, e.body.args[0], e.skolem_args)) - \
                        rec(expr.IndefiniteIntegral(e.var, e.body.args[1], e.skolem_args))
                elif expr.is_times(e.body) or expr.is_divides(e.body):
                    num_factors, denom_factors = decompose_expr_factor(e.body)
                    b = prod(f for f in num_factors if depends_on_var(f, e.var))
                    c = prod(f for f in num_factors if not depends_on_var(f, e.var))
                    denom_b = prod(f for f in denom_factors if depends_on_var(f, e.var))
                    denom_c = prod(f for f in denom_factors if not depends_on_var(f, e.var))
                    if denom_b != Const(1):
                        b = b / denom_b
                    if denom_c != Const(1):
                        c = c / denom_c
                    if c == expr.Const(1):
                        return IndefiniteIntegral(e.var, b, e.skolem_args)
                    else:
                        return c * rec(IndefiniteIntegral(e.var, b, e.skolem_args))
                else:
                    return e
            elif expr.is_limit(e):
                if expr.is_uminus(e.body):
                    return -Limit(e.var, e.lim, e.body.args[0])
                elif expr.is_times(e.body) or expr.is_divides(e.body):
                    num_factors, denom_factors = decompose_expr_factor(e.body)
                    b, c = Const(1), Const(1)
                    for f in num_factors:
                        if not depends_on_var(f, e.var):
                            c = c * f
                        else:
                            b = b * f
                    for f in denom_factors:
                        if not depends_on_var(f, e.var):
                            c = c / f
                        else:
                            b = b / f
                    return c * Limit(e.var, e.lim, b)
                else:
                    return e
            elif expr.is_summation(e):
                v, l, u, body = e.index_var, e.lower, e.upper, e.body
                if expr.is_minus(body):
                    return Summation(v, l, u, body.args[0]) - Summation(v, l, u, body.args[1])
                elif expr.is_uminus(body):
                    return -Summation(v, l, u, body.args[0])
                elif expr.is_times(e.body) or expr.is_divides(e.body):
                    num_factors, denom_factors = decompose_expr_factor(e.body)
                    b, c = Const(1), Const(1)
                    for f in num_factors:
                        if not depends_on_var(f, e.index_var):
                            c = c * f
                        else:
                            b = b * f
                    for f in denom_factors:
                        if not depends_on_var(f, e.index_var):
                            c = c / f
                        else:
                            b = b / f
                    c = normalize(c, ctx)
                    b = normalize(b, ctx)
                    return normalize(c * Summation(e.index_var, e.lower, e.upper, b), ctx)
                else:
                    return e
            else:
                return e

        return rec(e)


class PartialFractionDecomposition(Rule):
    """Apply partial fraction decomposition from sympy."""

    def __init__(self):
        self.name = "PartialFractionDecomposition"

    def __str__(self):
        return "partial fraction decomposition"

    def export(self):
        return {
            "name": self.name,
            "str": str(self)
        }

    def eval(self, e: Expr, ctx: Context) -> Expr:
        if not (expr.is_integral(e) or expr.is_indefinite_integral(e)):
            sep_ints = e.separate_integral()
            if len(sep_ints) == 0:
                return e
            else:
                return OnLocation(self, sep_ints[0][1]).eval(e, ctx)

        if not sympywrapper.is_rational(e.body):
            raise RuleException("PartialFractionDecomposition", "cannot be applied to non-rational body")

        new_body = normalize(sympywrapper.partial_fraction(e.body), ctx)
        if expr.is_integral(e):
            return expr.Integral(e.var, e.lower, e.upper, new_body)
        elif expr.is_indefinite_integral(e):
            return expr.IndefiniteIntegral(e.var, new_body, skolem_args=e.skolem_args)
        else:
            raise RuleException("PartialFractionDecomposition", "cannot be applied to non-integrals")


class ApplyIdentity(Rule):
    """Apply identities (trigonometric, etc) to the current term.

    The term that is rewritten to is always supplied, because there may
    be multiple options.

    """

    def __init__(self, source: Union[str, Expr], target: Union[str, Expr]):
        self.name = "ApplyIdentity"
        if isinstance(source, str):
            source = parser.parse_expr(source)
        if isinstance(target, str):
            target = parser.parse_expr(target)
        self.source = source
        self.target = target

    def __str__(self):
        return "rewrite %s to %s using identity" % (self.source, self.target)

    def export(self):
        return {
            "name": self.name,
            "str": str(self),
            "source": str(self.source),
            "target": str(self.target),
            "latex_str": "rewrite \\(%s\\) to \\(%s\\) using identity" % (
                latex.convert_expr(self.source), latex.convert_expr(self.target))
        }

    @staticmethod
    def search(e: Expr, ctx: Context) -> List[Expr]:
        res = []
        for identity in ctx.get_other_identities():
            inst = expr.match(e, identity.lhs)
            if inst is not None:
                # check identity's condition
                flag = True
                tmp_conds = [cond.inst_pat(inst) for cond in identity.conds.data]
                for cond in tmp_conds:
                    flag = flag and ctx.check_condition(cond)
                if flag:
                    expected_rhs = identity.rhs.inst_pat(inst)
                    res.append(normalize(expected_rhs, ctx))
        return res

    def eval(self, e: Expr, ctx: Context) -> Expr:
        # Find source within e
        if self.source != e:
            find_res = e.find_subexpr(self.source)
            if len(find_res) == 0:
                raise RuleException("ApplyIdentity", "old expression %s not found" % self.source)
            loc = find_res[0]
            return OnLocation(self, loc).eval(e, ctx)

        assert self.source == e
        for identity in ctx.get_other_identities():
            inst = expr.match(e, identity.lhs)
            if inst is not None:
                expected_rhs = identity.rhs.inst_pat(inst)
                tmp_conds = [cond.inst_pat(inst) for cond in identity.conds.data]
                flag = True
                for cond in tmp_conds:
                    flag = flag and ctx.check_condition(cond)
                if not flag:
                    continue
                if normalize(expected_rhs, ctx) == normalize(self.target, ctx):
                    return self.target

        raise RuleException("ApplyIdentity", "no matching identity for %s" % e)


class SeriesExpansionIdentity(Rule):
    """Apply series expansion in the current theory."""

    def __init__(self, *, old_expr: Optional[Union[str, Expr]] = None, index_var: str = 'n'):
        self.name = "SeriesExpansionIdentity"
        if isinstance(old_expr, str):
            old_expr = parser.parse_expr(old_expr)
        self.old_expr = old_expr
        self.index_var = index_var

    def __str__(self):
        return "apply series expansion on %s index %s" % (self.old_expr, self.index_var)

    def export(self):
        res = {
            "name": self.name,
            "str": str(self),
            "index_var": self.index_var
        }
        if self.old_expr is not None:
            res['old_expr'] = str(self.old_expr)
        return res

    def eval(self, e: Expr, ctx: Context) -> Expr:
        # If old_expr is given, try to find it within e
        if self.old_expr is not None and self.old_expr != e:
            find_res = e.find_subexpr(self.old_expr)
            if len(find_res) == 0:
                raise AssertionError("Equation: old expression not found")
            loc = find_res[0]
            return OnLocation(self, loc).eval(e, ctx)

        # Now e is the old expression
        assert self.old_expr is None or self.old_expr == e
        for identity in ctx.get_series_expansions():
            inst = expr.match(e, identity.lhs)
            if inst is None:
                continue

            # Check conditions
            satisfied = True
            for cond in identity.conds.data:
                cond = expr.expr_to_pattern(cond)
                cond = cond.inst_pat(inst)
                if not ctx.check_condition(cond):
                    satisfied = False

            if satisfied:
                res = identity.rhs.inst_pat(inst)
                assert expr.is_summation(res)
                return res.alpha_convert(self.index_var)

        # No matching identity found
        return e


class SeriesEvaluationIdentity(Rule):
    """Apply series evaluation in the current theory."""

    def __init__(self):
        self.name = "SeriesEvaluationIdentity"

    def __str__(self):
        return "apply series evaluation"

    def export(self):
        return {
            "name": self.name,
            "str": str(self)
        }

    def eval(self, e: Expr, ctx: Context) -> Expr:
        if not expr.is_summation(e):
            return e
        for identity in ctx.get_series_evaluations():
            inst = expr.match(e, identity.lhs)
            if inst is None:
                continue

            # Check conditions
            satisfied = True
            for cond in identity.conds.data:
                cond = expr.expr_to_pattern(cond)
                cond = cond.inst_pat(inst)
                if not ctx.check_condition(cond):
                    satisfied = False

            if satisfied:
                return identity.rhs.inst_pat(inst)

        # No matching identity found
        return e


class EvaluateIndefiniteIntegral(Rule):
    def __init__(self):
        self.name = "EvaluateIndefiniteIntegral"

    def eval(self, e: Expr, ctx: Context) -> Expr:
        assert isinstance(e, IndefiniteIntegral)

        ctx2 = context.body_conds(e, ctx)
        
        # 应用原有的积分规则（查找恒等式）
        for indef in ctx.get_indefinite_integrals():
            assert isinstance(indef.lhs, IndefiniteIntegral)
            inst = expr.match(e, indef.lhs)
            if inst is None:
                continue
            inst[indef.lhs.var] = Var(e.var)

            # Check conditions
            satisfied = True
            for cond in indef.conds.data:
                cond = expr.expr_to_pattern(cond)
                cond = cond.inst_pat(inst)
                if not ctx2.check_condition(cond):
                    satisfied = False

            if satisfied:
                # The right side of the identity should be of the form "expr + C"
                # take the expr part of the expression.
                assert indef.rhs.is_plus() and expr.is_skolem_func(indef.rhs.args[1])
                return indef.rhs.args[0].inst_pat(inst)

        # No matching identity found
        return e

class EvaluateDefiniteIntegral(Rule):
    def __init__(self):
        self.name = "EvaluateDefiniteIntegral"

    def eval(self, e: Expr, ctx: Context) -> Expr:
        assert isinstance(e, Integral)

        ctx2 = context.body_conds(e, ctx)
        
        # 查找定积分恒等式
        for identity in ctx.get_definite_integrals():
            inst = expr.match(e, identity.lhs)
            if inst is None:
                continue

            # Check conditions
            satisfied = True
            for cond in identity.conds.data:
                cond = expr.expr_to_pattern(cond)
                cond = cond.inst_pat(inst)
                if not ctx2.check_condition(cond):
                    satisfied = False
            if satisfied:
                return identity.rhs.inst_pat(inst)

        # Next, try indefinite integral identities
        for identity in ctx.get_indefinite_integrals():
            assert isinstance(identity.lhs, IndefiniteIntegral)
            inst = expr.match(IndefiniteIntegral(e.var, e.body, skolem_args=tuple()), identity.lhs)
            if inst is None:
                continue
            inst[identity.lhs.var] = Var(e.var)

            # Check conditions
            satisfied = True
            for cond in identity.conds.data:
                cond = expr.expr_to_pattern(cond)
                cond = cond.inst_pat(inst)
                if not ctx2.check_condition(cond):
                    satisfied = False

            if satisfied:
                # The right side of the identity should be of the form "expr + C"
                # take the expr part of the expression.
                assert identity.rhs.is_plus() and expr.is_skolem_func(identity.rhs.args[1])
                pat_rhs = identity.rhs.args[0]
                return EvalAt(e.var, e.lower, e.upper, normalize(pat_rhs.inst_pat(inst), ctx2))

        # No matching identity found
        return e

class IntegralIdentity(Rule):
    def __init__(self):
        self.name = "IntegralIdentity"

    def __str__(self):
        return "apply integral identity"

    def export(self):
        return {
            "name": self.name,
            "str": str(self)
        }
    
    def merge_evalat(self, e: Expr) -> Expr:
        """将包含多个EvalAt的表达式合并成单一EvalAt表达式"""
        # 处理嵌套的EvalAt和更复杂的表达式
        
        # 递归处理嵌套的EvalAt
        def process_nested_evalats(expr):
            # 先递归处理子表达式
            if isinstance(expr, Op):
                args = []
                for i in range(len(expr.args)):
                    args.append(process_nested_evalats(expr.args[i]))
                # 创建新的操作符表达式
                if expr.op == '+':
                    return args[0] + args[1] if len(args) > 1 else args[0]
                elif expr.op == '-':
                    if len(args) == 1:
                        return -args[0]  # 一元负号
                    else:
                        return args[0] - args[1]  # 二元减法
                elif expr.op == '*':
                    if len(args) == 1:
                        return args[0]
                    return args[0] * args[1]
                elif expr.op == '/':
                    return args[0] / args[1]
                elif expr.op == '=':
                    return Op('=', args[0], args[1])
                else:
                    # 其他操作符，尝试重构原始操作符
                    return Op(expr.op, *args)
            elif isinstance(expr, EvalAt):
                # 处理evalat内部的表达式
                processed_body = process_nested_evalats(expr.body)
                return EvalAt(expr.var, expr.lower, expr.upper, processed_body)
            else:
                # 基本表达式直接返回
                return expr
        
        # 合并同一层级的多个EvalAt表达式
        def merge_same_level_evalats(expr):
            # 如果不是操作符表达式,直接返回
            if not (isinstance(expr, Op) and expr.op in ['+', '-', '*']):
                return expr
            
            if len(expr.args) <= 1:
                return expr
            
            # 收集所有的EvalAt表达式
            evalats = []
            
            def collect_evalats(expr):
                if isinstance(expr, EvalAt):
                    evalats.append(expr)
                    return True
                elif isinstance(expr, Op) and expr.op in ['+', '-', '*']:
                    # 确保安全地访问参数
                    result = False
                    if len(expr.args) > 0:
                        result = collect_evalats(expr.args[0]) or result
                    if len(expr.args) > 1:
                        result = collect_evalats(expr.args[1]) or result
                    return result
                return False
            
            # 调用收集函数
            collect_evalats(expr)
            
            # 如果没有找到EvalAt表达式，直接返回原表达式
            if not evalats:
                return expr
                
            # 检查所有EvalAt是否有相同的变量和上下限
            var_name = evalats[0].var
            lower = evalats[0].lower
            upper = evalats[0].upper
            
            for evalat in evalats:
                if evalat.var != var_name or evalat.lower != lower or evalat.upper != upper:
                    # 如果变量或上下限不一致，不能合并
                    return expr
            
            # 替换所有EvalAt为它们的函数体
            def replace_evalat(expr):
                if isinstance(expr, EvalAt):
                    return expr.body
                elif isinstance(expr, Op):
                    # 确保安全地处理操作符参数
                    if expr.op == '+' and len(expr.args) >= 2:
                        return replace_evalat(expr.args[0]) + replace_evalat(expr.args[1])
                    elif expr.op == '-':
                        if len(expr.args) == 1:  # 一元负号
                            return -replace_evalat(expr.args[0])
                        elif len(expr.args) >= 2:  # 二元减法
                            return replace_evalat(expr.args[0]) - replace_evalat(expr.args[1])
                    elif expr.op == '*' and len(expr.args) >= 2:
                        return replace_evalat(expr.args[0]) * replace_evalat(expr.args[1])
                    elif expr.op == '/' and len(expr.args) >= 2:
                        return replace_evalat(expr.args[0]) / replace_evalat(expr.args[1])
                return expr
                
            # 创建合并后的函数体
            combined_body = replace_evalat(expr)
            
            return EvalAt(var_name, lower, upper, combined_body)
        
        # 首先递归处理子表达式中的嵌套EvalAt
        processed_expr = process_nested_evalats(e)
        
        # 然后尝试合并同级的EvalAt
        return merge_same_level_evalats(processed_expr)
    
    def eval(self, e: Expr, ctx: Context) -> Expr:
        """Apply indefinite integral identity to expression."""

        # If incoming expression is equality, apply to the right side of equation
        if e.is_equals():
            lhs = self.eval(e.lhs, ctx)
            rhs = self.eval(e.rhs, ctx)
            return Op("=", lhs, rhs)
        
        upper = None
        lower = None

        # Apply linearity first
        integrals = e.separate_integral()
        for _, loc in integrals:
            e = OnLocation(Linearity(), loc).eval(e, ctx)

        integrals = e.separate_integral()
        skolem_args = set()
        exist_indefinite_integral = False
        
        for sub_e, loc in integrals:
            if isinstance(sub_e, IndefiniteIntegral):
                e = OnLocation(EvaluateIndefiniteIntegral(), loc).eval(e, ctx)
                
                new_e = e.get_subexpr(loc)
                if new_e != sub_e:
                    exist_indefinite_integral = True
                    skolem_args = skolem_args.union(sub_e.skolem_args)
            elif isinstance(sub_e, Integral):
                upper = sub_e.upper.__str__()
                lower = sub_e.lower.__str__()
                e = OnLocation(EvaluateDefiniteIntegral(), loc).eval(e, ctx)
            else:
                raise AssertionError

        if exist_indefinite_integral:
            if expr.is_plus(e) and expr.is_skolem_func(e.args[1]):
                # If already has Skolem variable at right
                skolem_args = skolem_args.union(set(arg.name for arg in e.args[1].dependent_vars))
                e = e.args[0] + expr.SkolemFunc(e.args[1].name, tuple(Var(arg) for arg in skolem_args))
            else:
                # If no Skolem variable at right
                e = e + expr.SkolemFunc("C", tuple(Var(arg) for arg in skolem_args))

        # 计算表达式中EvalAt的数量
        def count_evalats(expr):
            count = 0
            if expr.is_equals():
                return count_evalats(expr.rhs)  # 只检查等式右侧
            
            # 查找所有EvalAt子表达式
            subexprs = expr.find_subexpr_pred(lambda x: isinstance(x, EvalAt))
            return len(subexprs)
        
        # 只有当存在多个EvalAt且上下限包含无穷大时，才进行合并
        if upper is not None and lower is not None:
            if (upper in ['oo', '-oo'] or lower in ['oo', '-oo']) and count_evalats(e) > 1:
                e = self.merge_evalat(e)
        return e

class CIntegralIdentity(Rule):
    """Apply contour integral identity to convert contour integral to ordinary integral."""
    
    def __init__(self):
        self.name = "CIntegralIdentity"

    def __str__(self):
        return "apply cintegral identity"

    def export(self):
        return {
            "name": self.name,
            "str": str(self)
        }
    
    def eval(self, e: Expr, ctx: Context) -> Expr:
        """应用围道积分恒等式: ∮_γ f(z) dz = ∫_a^b f(γ(t)) · γ'(t) dt"""  
        from integral.expr import CIntegral, CINTPath, Op
        
        # 如果输入表达式是等式，对两边都应用规则
        if e.is_equals():
            lhs = self.eval(e.lhs, ctx)
            rhs = self.eval(e.rhs, ctx)
            return Op("=", lhs, rhs)
        
        # 在表达式中查找所有围道积分
        cintegral_tuples = e.find_subexpr_pred(lambda x: isinstance(x, CIntegral))
        
        if not cintegral_tuples:
            return e  # 没有找到围道积分
        
        # 对每个围道积分应用变换
        result = e
        for cint, location in cintegral_tuples:
            try:
                if len(cint.paths) != 1:
                    continue  # 暂时跳过多路径积分
                
                path = cint.paths[0]
                
                # 处理不同的路径类型
                if isinstance(path, str):
                    # 路径是函数引用，如 C(t,r)
                    # 需要从上下文中查找定义
                    path_definition = self._lookup_path_definition(path, ctx)
                    if path_definition is None:
                        continue  # 如果找不到定义则跳过
                    path_expr = path_definition
                elif isinstance(path, CINTPath):
                    # 路径已经是 CINTPath 对象
                    path_expr = path
                else:
                    continue  # 跳过未知的路径类型
                
                # 将围道积分转换为普通积分
                ordinary_integral = self._convert_contour_to_integral(cint, path_expr, ctx)
                
                # 用普通积分替换围道积分
                result = result.replace(cint, ordinary_integral)
                
            except Exception:
                # 如果转换失败，继续处理下一个积分
                continue
        
        return result
    
    def _lookup_path_definition(self, path_ref: str, ctx: Context) -> CINTPath:
        """从上下文中查找路径定义"""
        from integral.expr import is_fun, CINTPath
        from integral.parser import parse_expr
        
        # 解析路径引用（例如 "C(t,r)"）
        try:
            path_call = parse_expr(path_ref)
            if not is_fun(path_call):
                return None
                
            func_name = path_call.func_name
            args = path_call.args
            
            # 在上下文中查找函数定义（包括父上下文）
            for definition in ctx.get_definitions():
                if (hasattr(definition.lhs, 'func_name') and 
                    definition.lhs.func_name == func_name and
                    len(definition.lhs.args) == len(args)):
                    
                    # 应用参数替换
                    path_def = definition.rhs
                    
                    # 确保 path_def 是 CINTPath
                    if isinstance(path_def, CINTPath):
                        # 应用参数替换 - 需要直接替换符号对象
                        new_path_expr = path_def.path_expr
                        new_start_expr = path_def.start_expr
                        new_end_expr = path_def.end_expr
                        
                        for param, arg in zip(definition.lhs.args, args):
                            # param 是符号对象，使用 replace 方法而不是 subst
                            new_path_expr = new_path_expr.replace(param, arg)
                            new_start_expr = new_start_expr.replace(param, arg)
                            new_end_expr = new_end_expr.replace(param, arg)
                        
                        return CINTPath(path_def.var, new_path_expr, new_start_expr, new_end_expr)
                    
        except Exception as e:
            pass
        
        return None
    
    def _convert_contour_to_integral(self, cint: CIntegral, path: CINTPath, ctx: Context) -> Integral:
        """使用定义将围道积分转换为普通积分"""
        from integral.expr import Integral, Op, Deriv
        
        # 获取路径组件
        param_var = path.var  # t
        path_expr = path.path_expr  # γ(t) = r*exp(i*pi*(1-t))
        start = path.start_expr  # a = 0
        end = path.end_expr  # b = 1
        
        # 计算路径导数 γ'(t)
        path_derivative = Deriv(param_var, path_expr)
        
        # 在 f(z) 中用 γ(t) 替换 z
        integrand_at_path = cint.body.subst(cint.var, path_expr)
        
        # 创建新的被积函数: f(γ(t)) · γ'(t)
        new_integrand = Op("*", integrand_at_path, path_derivative)
        
        # 创建普通积分: ∫_a^b f(γ(t)) · γ'(t) dt
        return Integral(param_var, start, end, new_integrand)

class ReplaceSubstitution(Rule):
    """Replace previously performed substitution"""

    def __init__(self):
        self.name = "ReplaceSubstitution"

    def __str__(self):
        return "replace substitution"

    def export(self):
        return {
            "name": self.name,
            "str": str(self)
        }

    def eval(self, e: Expr, ctx: Context) -> Expr:
        success = False
        for var, expr in reversed(ctx.get_substs()):
            if e.contains_var(var):
                success = True
            e = e.subst(var, expr)
        if not success:
            raise RuleException("ReplaceSubstitution", "No more substitution need to be replaced")
        return e


class DerivativeSimplify(Rule):
    """Simplify the derivative of an expression"""

    def __init__(self):
        self.name = "DerivativeSimplify"

    def __str__(self):
        return "simplify derivative"

    def export(self):
        return {
            "name": self.name,
            "str": str(self)
        }

    def eval(self, e: Expr, ctx: Context) -> Expr:
        if not isinstance(e, Deriv):
            return e
        return deriv(e.var, e.body, ctx)


class OnSubterm(Rule):
    """Apply given rule on subterms.

    The traversal order is similar to bottom-conv: first traverse each subterm
    of the term recursively, then apply the rule to the term itself.

    """

    def __init__(self, rule: Rule):
        assert isinstance(rule, Rule)
        self.rule = rule
        self.name = 'OnSubterm'

    def __str__(self):
        return "%s (all)" % self.rule

    def export(self):
        res = self.rule.export()
        res['str'] += ' (all)'
        res['loc'] = 'subterms'
        if 'latex_str' in res:
            res['latex_str'] += ' (all)'
        return res

    def update_context(self, e: Expr, ctx: Context) -> Context:
        return self.rule.update_context(e, ctx)

    def eval(self, e: Expr, ctx: Context) -> Expr:
        return apply_subterm(e, self.rule.eval, ctx)


class OnLocation(Rule):
    """Apply given rule on subterm specified by given location."""

    def __init__(self, rule: Rule, loc):
        assert isinstance(rule, Rule)
        self.name = "OnLocation"
        self.rule = rule
        self.loc = expr.Location(loc)

    def __str__(self):
        return "%s at %s" % (self.rule, self.loc)

    def export(self):
        res = self.rule.export()
        res['str'] += ' at ' + str(self.loc)
        res['loc'] = str(self.loc)
        if 'latex_str' in res:
            res['latex_str'] += ' at ' + str(self.loc)
        return res

    def update_context(self, e: Expr, ctx: Context) -> Context:
        return self.rule.update_context(e, ctx)

    def eval(self, e: Expr, ctx: Context) -> Expr:
        def rec(cur_e: Expr, loc: expr.Location, ctx: Context):
            if loc.is_empty():
                return self.rule.eval(cur_e, ctx)
            elif expr.is_var(cur_e) or expr.is_const(cur_e):
                raise AssertionError("OnLocation: invalid location")
            elif expr.is_op(cur_e):
                assert loc.head < len(cur_e.args), "OnLocation: invalid location"
                if len(cur_e.args) == 1:
                    return Op(cur_e.op, rec(cur_e.args[0], loc.rest, ctx))
                elif len(cur_e.args) == 2:
                    if loc.head == 0:
                        return Op(cur_e.op, rec(cur_e.args[0], loc.rest, ctx), cur_e.args[1])
                    elif loc.head == 1:
                        return Op(cur_e.op, cur_e.args[0], rec(cur_e.args[1], loc.rest, ctx))
                    else:
                        raise AssertionError("OnLocation: invalid location")
                else:
                    raise NotImplementedError
            elif expr.is_fun(cur_e):
                assert loc.head < len(cur_e.args), "OnLocation: invalid location"
                new_args = list(cur_e.args)
                new_args[loc.head] = rec(cur_e.args[loc.head], loc.rest, ctx)
                return Fun(cur_e.func_name, *tuple(new_args))
            elif expr.is_integral(cur_e):
                ctx2 = body_conds(cur_e, ctx)
                if loc.head == 0:
                    return Integral(cur_e.var, cur_e.lower, cur_e.upper, rec(cur_e.body, loc.rest, ctx2))
                elif loc.head == 1:
                    return Integral(cur_e.var, rec(cur_e.lower, loc.rest, ctx), cur_e.upper, cur_e.body)
                elif loc.head == 2:
                    return Integral(cur_e.var, cur_e.lower, rec(cur_e.upper, loc.rest, ctx), cur_e.body)
                else:
                    raise AssertionError("OnLocation: invalid location")
            elif expr.is_evalat(cur_e):
                if loc.head == 0:
                    return EvalAt(cur_e.var, cur_e.lower, cur_e.upper, rec(cur_e.body, loc.rest, ctx))
                elif loc.head == 1:
                    return EvalAt(cur_e.var, rec(cur_e.lower, loc.rest, ctx), cur_e.upper, cur_e.body)
                elif loc.head == 2:
                    return EvalAt(cur_e.var, cur_e.lower, rec(cur_e.upper, loc.rest, ctx), cur_e.body)
                else:
                    raise AssertionError("OnLocation: invalid location")
            elif expr.is_deriv(cur_e):
                assert loc.head == 0, "OnLocation: invalid location"
                return Deriv(cur_e.var, rec(cur_e.body, loc.rest, ctx))
            elif expr.is_limit(cur_e):
                if loc.head == 0:
                    if cur_e.lim.is_evaluable():
                        v = expr.eval_expr(cur_e.lim)
                        var = Var(cur_e.var)
                        if v == float('inf'):
                            cond = Op('>', var, Const(0))
                            ctx.add_condition(cond)
                        elif v == float('-inf'):
                            cond = Op('<', var, Const(0))
                            ctx.add_condition(cond)
                    return Limit(cur_e.var, cur_e.lim, rec(cur_e.body, loc.rest, ctx), drt=cur_e.drt)
                elif loc.head == 1:
                    return Limit(cur_e.var, rec(cur_e.lim, loc.rest, ctx), cur_e.body, drt=cur_e.drt)
                else:
                    raise AssertionError("OnLocation: invalid location")
            elif expr.is_indefinite_integral(cur_e):
                assert loc.head == 0, "OnLocation: invalid location"
                ctx2 = body_conds(cur_e, ctx)
                return IndefiniteIntegral(cur_e.var, rec(cur_e.body, loc.rest, ctx2), cur_e.skolem_args)
            elif expr.is_summation(cur_e):
                ctx2 = body_conds(cur_e, ctx)
                if loc.head == 0:
                    return Summation(cur_e.index_var, cur_e.lower, cur_e.upper, rec(cur_e.body, loc.rest, ctx2))
                elif loc.head == 1:
                    return Summation(cur_e.index_var, rec(cur_e.lower, loc.rest, ctx), cur_e.upper, cur_e.body)
                elif loc.head == 2:
                    return Summation(cur_e.index_var, cur_e.lower, rec(cur_e.upper, loc.rest, ctx), cur_e.body)
                else:
                    raise AssertionError("OnLocation: invalid location")
            else:
                raise NotImplementedError

        return rec(e, self.loc, ctx)


class OnCount(Rule):
    """Perform on the n'th subgoal satisfying some condition."""
    def __init__(self, rule: Rule, n: int, *, pred = None):
        self.rule = rule
        self.n = n
        if pred is None:
            if isinstance(rule, Rewriting):
                pred = lambda t: t == rule.old_expr
            elif isinstance(rule, (Substitution, SubstitutionInverse, IntegrationByParts, SplitRegion)):
                pred = lambda t: expr.is_integral(t) or expr.is_indefinite_integral(t)
            elif isinstance(rule, ExpandDefinition):
                pred = lambda t: expr.is_fun(t) and t.func_name == rule.func_name
            else:
                raise RuleException("OnCount", "(at n) should not be applied to rule %s" % rule.name)
        self.pred = pred

    def __str__(self):
        return "%s (at %s)" % (self.rule, self.n)

    def export(self):
        res = self.rule.export()
        res['str'] += ' (at %s)' % str(self.n)
        res['n'] = str(self.n)
        if 'latex_str' in res:
            res['latex_str'] += ' (at %s)' + str(self.n)
        return res

    def update_context(self, e: Expr, ctx: Context) -> Context:
        return self.rule.update_context(e, ctx)

    def eval(self, e: Expr, ctx: Context) -> Expr:
        count = self.n

        def rec(cur_e, ctx):
            nonlocal count
            if self.pred(cur_e):
                count -= 1
                if count == 0:
                    return self.rule.eval(cur_e, ctx)

            if expr.is_var(cur_e) or expr.is_const(cur_e) or expr.is_inf(cur_e):
                return cur_e
            elif expr.is_op(cur_e):
                return Op(cur_e.op, *(rec(arg, ctx) for arg in cur_e.args))
            elif expr.is_fun(cur_e):
                return Fun(cur_e.func_name, *(rec(arg, ctx) for arg in cur_e.args))
            elif expr.is_integral(cur_e):
                ctx2 = body_conds(cur_e, ctx)
                return Integral(cur_e.var, rec(cur_e.lower, ctx), rec(cur_e.upper, ctx), rec(cur_e.body, ctx2))
            elif expr.is_evalat(cur_e):
                return EvalAt(cur_e.var, rec(cur_e.lower, ctx), rec(cur_e.upper, ctx), rec(cur_e.body, ctx))
            elif expr.is_deriv(cur_e):
                return Deriv(cur_e.var, rec(cur_e.body, ctx))
            elif expr.is_limit(cur_e):
                ctx2 = body_conds(cur_e, ctx)
                return Limit(cur_e.var, rec(cur_e.lim, ctx), rec(cur_e.body, ctx2), drt=cur_e.drt)
            elif expr.is_indefinite_integral(cur_e):
                return IndefiniteIntegral(cur_e.var, rec(cur_e.body, ctx), cur_e.skolem_args)
            elif expr.is_summation(cur_e):
                ctx2 = body_conds(cur_e, ctx)
                return Summation(cur_e.index_var, rec(cur_e.lower, ctx), rec(cur_e.upper, ctx), rec(cur_e.body, ctx2))
            elif expr.is_skolem_func(cur_e):
                return SkolemFunc(cur_e.name, tuple(rec(arg, ctx) for arg in cur_e.dependent_vars))
            else:
                raise RuleException("OnCount", f"has not a implemention for {cur_e}.")

        res = rec(e, ctx)
        if count > 0:
            raise RuleException("OnCount", f"{self.n} is out of range")
        return res

class Simplify(Rule):
    """Perform simplification by applying the following rules repeatedly:

    - Apply Linearity.
    - Normalize using the rules in poly.
    - Apply DerivativeSimplify.

    """

    def __init__(self):
        self.name = "Simplify"

    def __str__(self):
        return "simplify"

    def export(self):
        return {
            "name": self.name,
            "str": str(self)
        }

    def eval(self, e: Expr, ctx: Context) -> Expr:
        if hasattr(e, 'needs_contour'):
            if e.needs_contour:
                raise RuleException(
                    "Simplify",
                    "ContourIntegral have no paths,please rewrite the expression."
                )
        
        counter = 0
        current = e
        while True:
            s = OnSubterm(Linearity()).eval(current, ctx)
            s = normalize(s, ctx)
            s = OnSubterm(DerivativeSimplify()).eval(s, ctx)
            if s == current:
                break
            current = s
            counter += 1
            if counter > 5:
                raise AssertionError("Loop in Simplify")
        return current


class ApplyEquation(Rule):
    """Apply the given equation for rewriting."""

    def __init__(self, eq: Union[Expr, str], source: Expr):
        self.name = "ApplyEquation"
        self.eq = eq
        self.source = source

    def __str__(self):
        return "apply %s on %s" % (self.eq, self.source)

    def latex_str(self):
        return "apply %s on \\(%s\\)" % (self.eq, latex.convert_expr(self.source))

    def export(self):
        res = {
            "name": self.name,
            "eq": str(self.eq),
            "str": str(self),
            "latex_str": self.latex_str()
        }
        if self.source:
            res['source'] = str(self.source)

        return res

    def eval(self, e: Expr, ctx: Context) -> Expr:
        # Find source within e
        if self.source is not None and self.source != e:
            find_res = e.find_subexpr(self.source)
            if len(find_res) == 0:
                raise RuleException("ApplyEquation", "source expression %s not found" % self.source)
            loc = find_res[0]
            return OnLocation(self, loc).eval(e, ctx)
        assert self.source == e or self.source is None

        # Find lemma
        found = False
        conds = None
        found_eq = None
        for identity in ctx.get_lemmas():
            if self.eq == identity.expr:
                found = True
                found_eq = self.eq
                conds = identity.conds.data
        if isinstance(self.eq, str):
            res = ctx.get_subgoal(self.eq)
            if res:
                found = True
                found_eq = res.expr
                conds = res.conds.data
        for item in ctx.get_eq_conds().data:
            if self.eq == item:
                if self.source is None:
                    if e == item.lhs:
                        return item.rhs
                    if e == item.rhs:
                        return item.lhs
                else:
                    if self.source == item.lhs:
                        return item.rhs
                    if self.source == item.rhs:
                        return item.lhs
                found = True
                found_eq = self.eq
                conds = []
        if not found:
            raise RuleException("ApplyEquation", f"lemma {self.eq} not found")

        # First try to match the current term with left or right side.
        pat = expr.expr_to_pattern(found_eq)
        conds_pattern = [expr.expr_to_pattern(cond) for cond in conds]
        inst_lhs = expr.match(e, pat.lhs)
        inst_rhs = expr.match(e, pat.rhs)
        if inst_lhs is not None:
            tmp = pat.rhs.inst_pat(inst_lhs)
            tmp_conds = [cond_pattern.inst_pat(inst_lhs) for cond_pattern in conds_pattern]
            left = pat.lhs.inst_pat(inst_lhs)
            if tmp is not None and normalize(left, ctx) == normalize(e, ctx):
                if not tmp_conds:
                    return tmp
                flag = True
                # check whether all conditions of the lemma have been satisfied
                for cond in tmp_conds:
                    flag = flag and ctx.check_condition(cond)
                if flag:
                    return tmp
        if inst_rhs is not None:
            tmp = pat.lhs.inst_pat(inst_rhs)
            tmp_conds = [cond_pattern.inst_pat(inst_rhs) for cond_pattern in conds_pattern]
            right = pat.rhs.inst_pat(inst_rhs)
            if tmp is not None and normalize(right, ctx) == normalize(e, ctx):
                if not tmp_conds:
                    return tmp
                flag = True
                # check whether all conditions of the lemma have been satisfied
                for cond in tmp_conds:
                    flag = flag and ctx.check_condition(cond)
                if flag:
                    return tmp

        # Finally, try to solve for e in the equation.
        res = solve_for_term(found_eq, e, ctx)
        if res is not None:
            flag = True
            for cond in conds:
                flag = flag and ctx.check_condition(cond)
            if flag:
                return res
        return e


class ApplyInductHyp(Rule):
    """Apply induction hypothesis."""

    def __init__(self):
        self.name = "ApplyInductHyp"

    def __str__(self):
        return "apply induction hypothesis"

    def export(self):
        return {
            "name": self.name,
            "str": str(self)
        }

    def eval(self, e: Expr, ctx: Context) -> Expr:
        for eq in ctx.get_induct_hyps():
            if e == eq.lhs:
                return eq.rhs
            if e == eq.rhs:
                return eq.lhs
            # pattern match
            eq_pat = expr.expr_to_pattern(eq.expr)
            lhs_inst = expr.match(e, eq_pat.lhs)
            if lhs_inst != None:
                return eq_pat.rhs.inst_pat(lhs_inst)
        # Not found
        return e


def normalize_divide(e1: Expr, e2: Expr):
    # First decompose into factors
    num_factors1, denom_factors1 = decompose_expr_factor2(e1)
    num_factors2, denom_factors2 = decompose_expr_factor2(e2)

    # Cancel out factors that are the same
    new_num_factors1 = []
    for factor in num_factors1:
        if factor in num_factors2:
            num_factors2.remove(factor)
        else:
            new_num_factors1.append(factor)
    new_denom_factors1 = []
    for factor in denom_factors1:
        if factor in denom_factors2:
            denom_factors2.remove(factor)
        else:
            new_denom_factors1.append(factor)

    def prod(es):
        es = list(es)
        if len(es) == 0:
            return Const(1)
        else:
            return functools.reduce(operator.mul, es[1:], es[0])

    new_num = prod(new_num_factors1 + denom_factors2)
    new_denom = prod(new_denom_factors1 + num_factors2)
    return new_num / new_denom

class Substitution(Rule):
    """Apply substitution u = g(x).

    var_name - str: name of the new variable.
    var_subst - Expr: expression in the original integral to be substituted.

    The identity to be applied is:

    INT x:[a, b]. f(g(x)) * g(x)' = INT u:[g(a), g(b)]. f(u)

    """

    def __init__(self, var_name: str, var_subst: Union[Expr, str]):
        if isinstance(var_subst, str):
            var_subst = parser.parse_expr(var_subst)
        assert isinstance(var_name, str) and isinstance(var_subst, Expr)
        self.name = "Substitution"
        self.var_name = var_name
        self.var_subst = var_subst
        self.f = None  # After application, record f here

    def __str__(self):
        return "substitute %s for %s" % (self.var_name, self.var_subst)

    def export(self):
        return {
            "name": self.name,
            "var_name": self.var_name,
            "var_subst": str(self.var_subst),
            "str": str(self),
            "latex_str": "substitute \\(%s\\) for \\(%s\\)" % \
                         (self.var_name, latex.convert_expr(self.var_subst))
        }

    def update_context(self, e: Expr, ctx: Context) -> Context:
        if not (expr.is_integral(e) or expr.is_indefinite_integral(e) or expr.is_limit(e)):
            sep_ints = e.separate_integral()
            sep_lims = e.separate_limits()
            if len(sep_ints) == 0 and len(sep_lims) == 0:
                return ctx
            elif len(sep_ints) != 0:
                e, _ = sep_ints[0]
            else:
                e, _ = sep_lims[0]

        if isinstance(e, IndefiniteIntegral):
            ctx2 = Context(ctx)
            ctx2.add_subst(self.var_name, self.var_subst)
            return ctx2
        else:
            return ctx

    def eval(self, e: Expr, ctx: Context) -> Expr:
        """
        Parameters:
        e: Expr, the integral on which to perform substitution.

        Returns:
        The new integral e', and stores in self.f the parameter used to
        specify the substitution.

        """
        # If not a direct integral/limit/cintegral, find one inside using find_subexpr_pred
        if not (expr.is_integral(e) or expr.is_indefinite_integral(e) or expr.is_limit(e) or expr.is_cintegral(e)):
            # Use find_subexpr_pred for unified search
            targets = e.find_subexpr_pred(lambda t: 
                expr.is_integral(t) or expr.is_indefinite_integral(t) or 
                expr.is_cintegral(t) or expr.is_limit(t))
            
            if len(targets) == 0:
                raise RuleException("Substitution", "integral or limit not found")
            
            # Apply to the first found target
            return OnLocation(self, targets[0][1]).eval(e, ctx)

        # If e is a limit, intelligently decide whether to apply to limit or inner integral
        if expr.is_limit(e):
            # Check if substitution involves the limit variable
            var_subst_test = self.var_subst.subst(e.var, Var(e.var))
            if e.var not in var_subst_test.get_vars():
                # Substitution doesn't involve limit variable
                # Look for integrals/cintegrals inside - they are more likely targets
                inner_targets = e.find_subexpr_pred(lambda t: 
                    expr.is_integral(t) or expr.is_cintegral(t))
                
                if len(inner_targets) > 0:
                    # Apply to the first inner integral/cintegral
                    return OnLocation(self, inner_targets[0][1]).eval(e, ctx)
                # Otherwise continue with limit substitution

        # Variable to be substituted in the integral
        var_name = Var(self.var_name)

        # Expression used for substitution
        var_subst = self.var_subst
        
        # For CIntegral with CINTPath, first apply cintegral identity to convert to ordinary integral
        # Then substitution can be applied to the resulting integral
        if expr.is_cintegral(e):
            # CIntegral should be converted to ordinary integral first via cintegral identity
            # Direct substitution on CIntegral is not supported
            # The user should apply "apply cintegral identity" first
            raise RuleException("Substitution", 
                "Cannot directly substitute in CIntegral. " +
                "Please apply 'cintegral identity' first to convert to ordinary integral.")
        
        # Check if this is actually an inverse substitution BEFORE substituting e.var
        # If var_subst doesn't contain e.var but contains other variables,
        # this is likely inverse substitution: e.var = var_subst
        if e.var not in var_subst.get_vars() and len(var_subst.get_vars()) > 0:
            # This is inverse substitution: e.var = var_subst
            # Delegate to SubstitutionInverse
            inv_rule = SubstitutionInverse(e.var, var_subst)
            return inv_rule.eval(e, ctx)
        
        var_subst = var_subst.subst(e.var, Var(e.var))
        
        if e.var not in var_subst.get_vars():
            raise RuleException("Substitution", "variable %s not found" % e.var)
        
        ctx2 = body_conds(e, ctx)

        # Compute g(x)'
        dfx = deriv(e.var, var_subst, ctx2)

        # If body is a product and g(x)' is on one of the sides, then
        # the new body is the other side. Otherwise, the new body is
        # obtained by dividing the original body by g(x)'.
        body = normalize_divide(e.body, dfx)

        # Now attempt to write the new body in the form of f(g(x)).
        # First substitute all appearances of g(x) by u. If this clears
        # all appearances of x, then we are done. Otherwise, we need
        # to solve x as a function of u, then replace x by that function.
        body_subst = body.replace(var_subst, var_name)
        def prod(es):
            es = list(es)
            if len(es) == 0:
                return Const(1)
            else:
                return functools.reduce(operator.mul, es[1:], es[0])
        nf, df = decompose_expr_factor2(var_subst)
        prod_nf, prod_df = prod(nf), prod(df)
        var_subst2 = prod_nf / prod_df if prod_df != Const(1) else prod_nf
        body_subst2 = normalize(body, ctx2).replace(normalize(var_subst, ctx2), var_name)
        body_subst3 = body.replace(var_subst2, var_name)
        body_subst4 = normalize(body, ctx2).replace(normalize(var_subst2,ctx2), var_name)
        body_subst5 = normalize(body.replace(var_subst, var_name), ctx2)
        body_subst6 = normalize(body.replace(var_subst2, var_name), ctx2)
        if e.var not in body_subst.get_vars():
            # Substitution is able to clear all x in original integrand
            self.f = body_subst
        elif e.var not in body_subst2.get_vars():
            self.f = body_subst2
        elif e.var not in body_subst3.get_vars():
            self.f = body_subst3
        elif e.var not in body_subst4.get_vars():
            self.f = body_subst4
        elif e.var not in body_subst5.get_vars():
            self.f = body_subst5
        elif e.var not in body_subst6.get_vars():
            self.f = body_subst6
        else:
            # Substitution is unable to clear x, need to solve for x
            gu_list = solve_equation(var_subst, var_name, e.var, ctx2)
            if not gu_list:
                raise RuleException("Substitution", "unable to solve equation %s = %s for %s, body_subst = %s" % (
                    var_subst, var_name, e.var, body_subst
                ))

            gu = normalize(gu_list[0], ctx2)  # Take the first solution
            c = e.body.replace(Var(e.var), gu)
            if not expr.is_limit(e):
                new_problem_body = c * deriv(str(var_name), gu, ctx2)
            else:
                new_problem_body = c
            self.f = new_problem_body

        if expr.is_integral(e):
            if e.lower == expr.NEG_INF:
                lower = limits.reduce_neg_inf_limit(var_subst, e.var, ctx2)
            else:
                # 计算替换后的下限
                try:
                    lower = normalize(var_subst.subst(e.var, e.lower), ctx2)
                except ZeroDivisionError:
                    # 如果出现除零,说明替换后可能是无穷
                    x = Var(e.var)
                    lower = limits.reduce_inf_limit(var_subst.subst(e.var, e.lower + (1/x)), e.var, ctx2)

            if e.upper == expr.POS_INF:
                upper = limits.reduce_inf_limit(var_subst, e.var, ctx2)
            else:
                # 计算替换后的上限
                try:
                    upper = normalize(var_subst.subst(e.var, e.upper), ctx2)
                except ZeroDivisionError:
                    # 如果出现除零,说明替换后可能是无穷
                    x = Var(e.var)
                    upper = limits.reduce_inf_limit(var_subst.subst(e.var, e.upper - (1/x)), e.var, ctx2)

            if lower.is_evaluable() and upper.is_evaluable() and expr.eval_expr(lower) > expr.eval_expr(upper):
                return normalize(Integral(self.var_name, upper, lower, Op("-", self.f)), ctx2)
            else:
                return normalize(Integral(self.var_name, lower, upper, self.f), ctx2)
        elif expr.is_indefinite_integral(e):
            return normalize(IndefiniteIntegral(self.var_name, self.f, e.skolem_args), ctx2)
        elif expr.is_limit(e):
            # Perhaps need to be improved when drt is not None
            if e.lim == expr.NEG_INF:
                lim = limits.reduce_neg_inf_limit(var_subst, e.var, ctx2)
            elif e.lim == expr.POS_INF:
                lim = limits.reduce_inf_limit(var_subst, e.var, ctx2)
            else:
                x = Var(e.var)
                left = self.var_subst
                left = limits.reduce_inf_limit(left.subst(e.var, (1 / x) + e.lim), e.var, ctx2)
                left = normalize(left, ctx2)
                right = self.var_subst
                right = limits.reduce_inf_limit(right.subst(e.var, e.lim - (1 / x)), e.var, ctx2)
                right = normalize(right, ctx2)
                if left.is_evaluable() and right.is_evaluable() and expr.eval_expr(left) == expr.eval_expr(right):
                    return normalize(Limit(self.var_name, left, self.f, None), ctx2)
                else:
                    return e
            return normalize(Limit(self.var_name, lim, self.f, None), ctx2)
        else:
            raise TypeError


class SubstitutionInverse(Rule):
    """Apply substitution x = f(u).

    Grammar
    -------

        "substitute" expr "for" CNAME
    
    Attributes
    ----------
    old_var: str
        name of the original variable of integration.
    var_subst: Expr
        expression containing the new variable.

    """
    def __init__(self, old_var: str, var_subst: Union[Expr, str]):
        self.name = "SubstitutionInverse"
        self.old_var = old_var
        if isinstance(var_subst, str):
            var_subst = parser.parse_expr(var_subst)
        self.var_subst = var_subst

    def __str__(self):
        return "substitute %s for %s" % (self.var_subst, self.old_var)

    def export(self):
        return {
            "name": self.name,
            "old_var": self.old_var,
            "var_subst": str(self.var_subst),
            "str": str(self),
            "latex_str": "substitute \\(%s\\) for \\(%s\\)" % (
                latex.convert_expr(self.var_subst), self.old_var)
        }

    def update_context(self, e: Expr, ctx: Context) -> Expr:
        if not (expr.is_integral(e) or expr.is_indefinite_integral(e)):
            sep_ints = e.separate_integral()
            if len(sep_ints) == 0:
                return ctx
            else:
                e, _ = sep_ints[0]

        new_vars = self.var_subst.get_vars() - set(ctx.get_vars()) - {e.var}
        if len(new_vars) >= 2:
            return ctx

        if not new_vars:
            return ctx

        new_var = new_vars.pop()

        if isinstance(e, IndefiniteIntegral):
            ctx2 = Context(ctx)
            inv_f_list = solve_equation(self.var_subst, Var(e.var), new_var, ctx)
            if inv_f_list:
                ctx2.add_subst(new_var, inv_f_list[0])  # Take the first solution
            return ctx2
        else:
            return ctx

    def eval(self, e: Expr, ctx: Context) -> Expr:
        if not (expr.is_integral(e) or expr.is_indefinite_integral(e) or expr.is_cintegral(e)):
            sep_ints = e.separate_integral()
            sep_cints = e.separate_cintegral()
            if len(sep_ints) == 0 and len(sep_cints) == 0:
                raise RuleException("SubstitutionInverse", "no integral found in expression")
            elif expr.is_cintegral(e):
                return OnLocation(self, sep_cints[0][1]).eval(e, ctx)
            else:
                return OnLocation(self, sep_ints[0][1]).eval(e, ctx)

        if not (expr.is_integral(e) or expr.is_indefinite_integral(e) or expr.is_cintegral(e)):
            raise RuleException("SubstitutionInverse", "input is not integral")

        if e.var != self.old_var:
            raise RuleException("SubstitutionInverse", "incorrect old variable %s, should be %s" % (
                self.old_var, e.var))

        new_vars = self.var_subst.get_vars() - set(ctx.get_vars()) - {e.var}
        if len(new_vars) >= 2:
            raise RuleException("SubstitutionInverse", "more than one new variable: %s" % str(new_vars))

        if not new_vars:
            raise RuleException("SubstitutionInverse", "no new variable found in substitution")

        new_var = new_vars.pop()

        try:
            # dx = f'(u) * du
            subst_deriv = deriv(new_var, self.var_subst, ctx)
        except NotImplementedError:
            raise RuleException('Inverse Substitute', f"{self.var_subst} can not be derived")

        # Replace x with f(u)
        new_e_body = e.body.replace(Var(e.var), self.var_subst)

        # g(x) = g(x(u)) * f'(u)
        new_e_body = new_e_body * subst_deriv

        # Solve the equations f(u) = x for u
        inv_f_list = solve_equation(self.var_subst, Var(e.var), new_var, ctx)
        if not inv_f_list:
            raise RuleException("SubstitutionInverse", "cannot solve equation %s = %s for %s" % (
                self.var_subst, e.var, new_var
            ))
        inv_f = normalize(inv_f_list[0], ctx)  # Take the first solution

        if expr.is_integral(e):
            lower = limits.reduce_inf_limit(inv_f.subst(e.var, (1 / Var(e.var)) + e.lower), e.var, ctx)
            upper = limits.reduce_inf_limit(inv_f.subst(e.var, e.upper - (1 / Var(e.var))), e.var, ctx)

            lower = normalize(lower, ctx)
            upper = normalize(upper, ctx)
            if lower.is_evaluable() and upper.is_evaluable() and expr.eval_expr(lower) > expr.eval_expr(upper):
                return -expr.Integral(new_var, upper, lower, new_e_body)
            else:
                return expr.Integral(new_var, lower, upper, new_e_body)
        elif expr.is_indefinite_integral(e):
            return expr.IndefiniteIntegral(new_var, new_e_body, skolem_args=e.skolem_args)
        else:
            raise AssertionError("SubstitutionInverse")

def try_cintegral_to_integral(e: Expr, new_expr: Expr, ctx: Context) -> Optional[Expr]:
    """尝试将围道积分转换为普通积分。
    
    当围道积分的路径是实值路径时，可以转换为普通定积分。
    规则: CINT z:path. f(z) → INT x:[a,b]. f(x) 或 -INT x:[a,b]. f(x)
    
    Args:
        e: 原始围道积分表达式
        new_expr: 目标表达式（普通积分或其相反数）
        ctx: 上下文
        
    Returns:
        如果转换成功，返回 new_expr；否则返回 None
        
    Raises:
        RuleException: 当方向和符号不匹配时抛出具体错误
    """
    # 处理实值围道积分到普通积分的转换
    # 围道积分到普通积分的转换规则: CINT z:path. f(z) → INT x:[a,b]. f(x) 或 -INT x:[a,b]. f(x) (当路径是实值时)
    # 先检查 e 是否是单路径围道积分，避免不必要的计算
    if expr.is_cintegral(e) and len(e.paths) == 1:
        # 检查new_expr是否是普通积分或其相反数
        target_integral = None
        needs_negation = False
        
        if expr.is_integral(new_expr):
            target_integral = new_expr
            needs_negation = False
        elif isinstance(new_expr, expr.Op) and new_expr.op == '-' and len(new_expr.args) == 1:
            if expr.is_integral(new_expr.args[0]):
                target_integral = new_expr.args[0]
                needs_negation = True
        elif isinstance(new_expr, expr.Op) and new_expr.op == '*' and len(new_expr.args) == 2:
            if new_expr.args[0] == expr.Const(-1) and expr.is_integral(new_expr.args[1]):
                target_integral = new_expr.args[1]
                needs_negation = True
            elif new_expr.args[1] == expr.Const(-1) and expr.is_integral(new_expr.args[0]):
                target_integral = new_expr.args[0]
                needs_negation = True
        
        if target_integral is not None:
            path = e.paths[0]
            path_obj = None
            matched_defn = None
            inst = None
            
            # 获取实际的 CINTPath 对象
            if isinstance(path, expr.CINTPath):
                path_obj = path
            elif isinstance(path, str):
                # 路径是字符串引用，在上下文定义中查找
                try:
                    path_ref = parser.parse_expr(path)
                    for defn in ctx.get_definitions():
                        inst = expr.match(path_ref, defn.lhs)
                        if inst is not None and isinstance(defn.rhs, expr.CINTPath):
                            path_obj = defn.rhs.inst_pat(inst)
                            matched_defn = defn
                            break
                except:
                    pass
            
            if path_obj is not None:
                path_expr = path_obj.path_expr
                path_var = path_obj.var
                lower_t = path_obj.start_expr
                upper_t = path_obj.end_expr
                
                # 先计算端点，这是轻量级操作
                z_start = normalize(path_expr.subst(path_var, lower_t), ctx)
                z_end = normalize(path_expr.subst(path_var, upper_t), ctx)
                
                # 检查边界是否匹配（先做轻量级检查）
                new_lower = normalize(target_integral.lower, ctx)
                new_upper = normalize(target_integral.upper, ctx)
                
                forward_match = (z_start == new_lower and z_end == new_upper)
                reverse_match = (z_start == new_upper and z_end == new_lower)
                
                if forward_match or reverse_match:
                    # 在变量替换后验证主体是否匹配
                    expected_body = e.body.subst(e.var, expr.Var(target_integral.var))
                    
                    if normalize(expected_body, ctx) == normalize(target_integral.body, ctx):
                        # 只有在边界和主体都匹配时，才进行昂贵的 notReal 检查
                        # 创建带有定义条件的临时上下文
                        from integral.context import Context as ContextClass
                        check_ctx = ContextClass(ctx)
                        check_ctx.add_condition(expr.Fun("isReal", expr.Var(path_var)))
                        
                        if matched_defn is not None:
                            for cond in matched_defn.conds.data:
                                inst_cond = cond.inst_pat(inst) if inst else cond
                                check_ctx.add_condition(inst_cond)
                            if inst:
                                for _, var_expr in inst.items():
                                    if expr.is_var(var_expr):
                                        check_ctx.add_condition(expr.Fun("isReal", var_expr))
                        
                        # 检查 notReal - 只有在其他条件都满足时才执行
                        is_not_real = check_ctx.check_condition(expr.Fun("notReal", path_expr))
                        
                        if not is_not_real:  # 路径可能是实数
                            # 检查方向和符号是否一致
                            if (forward_match and not needs_negation) or (reverse_match and needs_negation):
                                return new_expr
                            else:
                                # 方向和符号不匹配，给出具体的错误提示
                                if reverse_match and not needs_negation:
                                    raise RuleException(
                                        "Rewriting",
                                        f"The contour path is in reverse direction (from {z_end} to {z_start}), "
                                        f"but the rewritten integral {new_expr} is missing the negative sign. "
                                        f"Please rewrite to -{new_expr}"
                                    )
                                elif forward_match and needs_negation:
                                    raise RuleException(
                                        "Rewriting",
                                        f"The contour path is in forward direction (from {z_start} to {z_end}), "
                                        f"but the rewritten integral has an unnecessary negative sign. "
                                        f"Please rewrite to {target_integral} instead of {new_expr}"
                                    )
    return None


class ExpandPolynomial(Rule):
    """Expand multiplication and power."""

    def __init__(self):
        self.name = "ExpandPolynomial"

    def __str__(self):
        return "expand polynomial"

    def export(self):
        return {
            "name": self.name,
            "str": str(self)
        }

    def eval(self, e: Expr, ctx: Context) -> Expr:
        # Case of constant, integer power
        if expr.is_power(e) and expr.is_const(e.args[1]) and e.args[1].val > 1 and \
                int(e.args[1].val) == e.args[1].val:
            n = int(e.args[1].val)
            base = to_poly(self.eval(e.args[0], ctx), ctx)
            res = base
            for i in range(n - 1):
                res = res * base
            return from_poly(res.reduce(ctx))

        # Case of product: carry out the multiplication
        elif e.is_times():
            s1, s2 = self.eval(e.args[0], ctx), self.eval(e.args[1], ctx)
            return from_poly((to_poly(s1, ctx) * to_poly(s2, ctx)).reduce(ctx))

        # Case of divide: if denominator is monomial, expand fully, otherwise
        # expand the numerator only.
        elif e.is_divides():
            s1, s2 = self.eval(e.args[0], ctx), self.eval(e.args[1], ctx)
            p1, p2 = to_poly(s1, ctx), to_poly(s2, ctx)
            if p2.is_monomial():
                return from_poly((p1 / p2).reduce(ctx))
            else:
                return from_poly((p1 / poly.singleton(from_poly(p2))).reduce(ctx))

        # Case of plus and minus: expand on both sides
        elif e.is_plus() or e.is_minus():
            e = OnLocation(self, "0").eval(e, ctx)
            e = OnLocation(self, "1").eval(e, ctx)
            return e

        # Case of uminus, expand subterm
        elif expr.is_uminus(e):
            return OnLocation(self, "0").eval(e, ctx)

        # Case of integrals
        elif expr.is_integral(e):
            ctx2 = body_conds(e, ctx)
            return expr.Integral(e.var, e.lower, e.upper, self.eval(e.body, ctx2))
        elif expr.is_indefinite_integral(e):
            ctx2 = body_conds(e, ctx)
            return expr.IndefiniteIntegral(e.var, self.eval(e.body, ctx2), e.skolem_args)
        else:
            return e

class Rewriting(Rule):
    def __init__(self, old_expr: Optional[Union[str, Expr]], new_expr: Union[str, Expr]):
        self.name = "Equation"
        if isinstance(old_expr, str):
            old_expr = parser.parse_expr(old_expr)
        if isinstance(new_expr, str):
            new_expr = parser.parse_expr(new_expr)
        self.old_expr = old_expr
        self.new_expr = new_expr

    def __str__(self):
        if self.old_expr is None:
            return "rewrite to %s" % self.new_expr
        else:
            return "rewrite %s to %s" % (self.old_expr, self.new_expr)

    def export(self):
        if self.old_expr is None:
            latex_str = "rewrite to \\(%s\\)" % latex.convert_expr(self.new_expr)
        else:
            latex_str = "rewrite \\(%s\\) to \\(%s\\)" % \
                        (latex.convert_expr(self.old_expr), latex.convert_expr(self.new_expr))
        res = {
            "name": self.name,
            "new_expr": str(self.new_expr),
            "str": str(self),
            "latex_str": latex_str
        }
        if self.old_expr:
            res['old_expr'] = str(self.old_expr)
        return res

    def eval(self, e: Expr, ctx: Context) -> Expr:
        
        if self.old_expr is not None and self.old_expr != e:
            find_res = e.find_subexpr(self.old_expr)
            if len(find_res) == 0:
                raise RuleException(
                    "Rewriting", f"old expression {self.old_expr} not found in {e}")
            loc = find_res[0]
            return OnLocation(self, loc).eval(e, ctx)

        # Now e is the old expression
        assert self.old_expr is None or self.old_expr == e

        r = Simplify()
        r1, r2 = r.eval(e, ctx), r.eval(self.new_expr, ctx)
        if r1 == r2:
            return self.new_expr

        # Handle infinity cases with products
        if expr.is_op(e) and e.op == '*':
            # If new_expr is a limit
            if expr.is_limit(self.new_expr):
                lim = self.new_expr
                # Check if all factors in the product are exponential functions
                all_exp = all(expr.is_fun(arg) and arg.func_name == 'exp' for arg in e.args)
                if all_exp:
                    # Check if the exponents contain infinity
                    has_inf = any(expr.is_inf(arg.args[0]) or (expr.is_op(arg.args[0]) and
                                any(expr.is_inf(term) for term in arg.args[0].args))
                                for arg in e.args)
                    if has_inf:
                        # Replace infinity with limit variable in each factor
                        new_args = []
                        for arg in e.args:
                            if expr.is_fun(arg) and arg.func_name == 'exp':
                                new_body = arg.args[0]
                                if expr.is_inf(new_body):
                                    new_body = Var(lim.var)
                                elif expr.is_op(new_body):
                                    for i, term in enumerate(new_body.args):
                                        if expr.is_inf(term):
                                            new_body = new_body.replace(term, Var(lim.var))
                                new_args.append(Fun('exp', new_body))
                        expected = Limit(lim.var, expr.POS_INF, functools.reduce(lambda x, y: Op('*', x, y), new_args))
                        if normalize(expected, ctx) == normalize(self.new_expr, ctx):
                            return self.new_expr

        # Handle single exponential function
        if expr.is_fun(e) and e.func_name == 'exp':
            if len(e.args) == 1 and expr.is_op(e.args[0]) and e.args[0].op == '*':
                if any(expr.is_inf(arg) for arg in e.args[0].args):
                    # Check if new_expr is a limit expression
                    if expr.is_limit(self.new_expr):
                        lim = self.new_expr
                        if expr.is_fun(lim.body) and lim.body.func_name == 'exp':
                            # Replace infinity with limit variable
                            new_body = e.args[0]
                            for i, arg in enumerate(new_body.args):
                                if expr.is_inf(arg):
                                    new_body = new_body.replace(arg, Var(lim.var))
                            expected = Limit(lim.var, expr.POS_INF, Fun('exp', new_body))
                            if normalize(expected, ctx) == normalize(self.new_expr, ctx):
                                return self.new_expr

        # Rewriting 1 to sin(x)^2 + cos(x)^2
        x = Symbol("x", [VAR, CONST, OP, FUN])
        q1 = expr.sin(x) ** 2 + expr.cos(x) ** 2
        q2 = expr.cos(x) ** 2 + expr.sin(x) ** 2
        if e == Const(1) and expr.match(self.new_expr, q1) or expr.match(self.new_expr, q2):
            return self.new_expr

        # Rewriting sin(x)^2 + cos(x)^2 to 1
        p1 = expr.sin(x) ** 2 + expr.cos(x) ** 2
        p2 = expr.cos(x) ** 2 + expr.sin(x) ** 2
        if (expr.match(e, p1) or expr.match(e, p2)) and self.new_expr == Const(1):
            return self.new_expr

        if norm.eq_quotient(e, self.new_expr, ctx):
            return self.new_expr

        if norm.eq_power(e, self.new_expr, ctx):
            return self.new_expr

        if norm.eq_log(e, self.new_expr, ctx):
            return self.new_expr

        if norm.eq_definite_integral(e, self.new_expr, ctx):
            return self.new_expr

        if norm.simp_definite_integral(e, ctx) == normalize(self.new_expr, ctx):
            return self.new_expr

        # x * sum(k,l,u,body) => sum(k, l, u, x* body)
        x = Symbol('x', [VAR, CONST, OP, FUN])
        y = Symbol('y', [SUMMATION])
        p = x * y
        mapping = expr.match(e, p)
        if mapping is not None:
            sum = mapping[y.name]
            idx = sum.index_var
            out = mapping[x.name]
            if idx not in out.get_vars():
                e = Summation(idx, sum.lower, sum.upper, out * sum.body)

        # sum(k, l, u, body1) + sum(i, l, u, body2) => sum(k, l, u, body1+body2)
        x = Symbol('x', [SUMMATION])
        y = Symbol('y', [SUMMATION])
        p = x + y
        mapping = expr.match(e, p)
        if mapping is not None:
            sum1: Summation = mapping[x.name]
            sum2: Summation = mapping[y.name]
            if sum1.lower == sum2.lower and sum1.upper == sum2.upper:
                e = Summation(sum1.index_var, sum1.lower, sum1.upper, sum1.body + sum2.body)
            if normalize(e, ctx) == normalize(self.new_expr, ctx):
                return self.new_expr

        if expr.is_summation(e):
            # SUM(i, 0, oo, body) -> LIM {n->oo}. SUM(i, 0, n, body)
            if e.upper == expr.POS_INF:
                v = e.index_var + e.index_var
                tmp = Limit(v, expr.POS_INF, Summation(e.index_var, e.lower, Var(v), e.body))
                if normalize(tmp, ctx) == normalize(self.new_expr, ctx):
                    return self.new_expr
            # sum(k, l, u, body1) + sum(i, l, u, body2) <== sum(k, l, u, body1+body2)
            if expr.is_op(e.body) and e.body.op in '+-':
                v, l, u = e.index_var, e.lower, e.upper
                tmp = Op(e.body.op, Summation(v, l, u, e.body.args[0]), Summation(v, l, u, e.body.args[1]))
                if normalize(tmp, ctx) == normalize(self.new_expr, ctx):
                    return self.new_expr
            return e
        elif expr.is_fun(e) and e.func_name == "Gamma":
            # 确保 Gamma 函数只有一个参数
            if len(e.args) == 1:
                arg = e.args[0]
                conds = ctx.get_conds(arg)
                # 检查参数是否为整数表达式
                if expr.is_var(arg) and "isInt" in conds:
                    # 创建 factorial(n-1) 表达式
                    return expr.Fun("factorial", expr.Op("-", arg, expr.Const(1)))

        # 处理围道积分拆分规则
        # 围道积分拆分规则: CINT z:com(path1,path2,...). f(z) → CINT z:path1. f(z) + CINT z:path2. f(z) + ...
        if expr.is_cintegral(e) and len(e.paths) >= 2:
            # 构建期望的拆分形式
            split_integrals = []
            for path in e.paths:
                single_cint = expr.CIntegral(e.var, [path], e.body)
                split_integrals.append(single_cint)
            
            # 构建求和: first + second + third + ...
            expected = split_integrals[0]
            for cint in split_integrals[1:]:
                expected = expr.Op("+", expected, cint)
            
            # 检查new_expr是否匹配期望的拆分形式
            if normalize(expected, ctx) == normalize(self.new_expr, ctx):
                return self.new_expr
        
        # 处理实值围道积分到普通积分的转换
        cint_result = try_cintegral_to_integral(e, self.new_expr, ctx)
        if cint_result is not None:
            return cint_result
        
        # apply identity
        for identity in ctx.get_other_identities():
            inst = expr.match(e, identity.lhs)
            if inst is not None:
                expected_rhs = identity.rhs.inst_pat(inst)
                tmp_conds = [cond.inst_pat(inst) for cond in identity.conds.data]
                flag = True
                for cond in tmp_conds:
                    flag = flag and ctx.check_condition(cond)
                if not flag:
                    continue
                if normalize(expected_rhs, ctx) == normalize(self.new_expr, ctx):
                    return self.new_expr
        raise RuleException("Rewriting", "rewriting %s to %s failed" % (e, self.new_expr))


class IntegrationByParts(Rule):
    """Apply integration by parts.

    The arguments `u` and `v` should satisfy `u * dv` equals the integrand.
    This step transforms `INT x. u * dv` into `u * v - INT x. v * du`.

    """

    def __init__(self, u: Union[str, Expr], v: Union[str, Expr]):
        self.name = "IntegrationByParts"
        if isinstance(u, str):
            u = parser.parse_expr(u)
        if isinstance(v, str):
            v = parser.parse_expr(v)
        assert isinstance(u, Expr) and isinstance(v, Expr)
        self.u = u
        self.v = v

    def __str__(self):
        return "integrate by parts with u = %s, v = %s" % (self.u, self.v)

    def export(self):
        return {
            "name": self.name,
            "u": str(self.u),
            "v": str(self.v),
            "str": str(self),
            "latex_str": "integrate by parts with \\(u = %s, v = %s\\)" % \
                         (latex.convert_expr(self.u), latex.convert_expr(self.v))
        }

    def eval(self, e: Expr, ctx: Context) -> Expr:
        if not (expr.is_integral(e) or expr.is_indefinite_integral(e)):
            sep_ints = e.separate_integral()
            if len(sep_ints) == 0:
                return e
            else:
                return OnLocation(self, sep_ints[0][1]).eval(e, ctx)

        ctx2 = body_conds(e, ctx)
        du = deriv(e.var, self.u, ctx)
        dv = deriv(e.var, self.v, ctx)
        udv = normalize(self.u * dv, ctx2)

        equal = False
        if udv == normalize(e.body, ctx2):
            equal = True

        if not equal and norm.eq_quotient(udv, e.body, ctx2):
            equal = True

        if not equal and norm.eq_power(udv, e.body, ctx2):
            equal = True

        if equal:
            if expr.is_integral(e):
                return expr.EvalAt(e.var, e.lower, e.upper, normalize(self.u * self.v, ctx2)) - \
                       expr.Integral(e.var, e.lower, e.upper, normalize(self.v * du, ctx2))
            elif expr.is_indefinite_integral(e):
                return normalize(self.u * self.v, ctx2) - \
                       expr.IndefiniteIntegral(e.var, normalize(self.v * du, ctx2), e.skolem_args)
        else:
            raise RuleException(self.name, f"u * dv does not equal body: {udv} != {e.body}")


class SplitRegion(Rule):
    """Split integral into two parts at a point."""

    def __init__(self, c: Union[Expr, str]):
        self.name = "SplitRegion"
        if isinstance(c, str):
            c = parser.parse_expr(c)
        self.c = c

    def __str__(self):
        return "split region at %s" % self.c

    def export(self):
        return {
            "name": self.name,
            "c": str(self.c),
            "str": str(self)
        }

    def eval(self, e: Expr, ctx: Context) -> Expr:
        if not expr.is_integral(e):
            sep_ints = e.separate_integral()
            if len(sep_ints) == 0:
                return e
            else:
                try:
                    return OnLocation(self, sep_ints[0][1]).eval(e, ctx)
                except RecursionError:
                    raise RuleException("split region", "do not find a definite integral")
        x = Var("c")
        is_cpv = limits.reduce_inf_limit(e.body.subst(e.var, self.c + 1 / x), x.name, ctx) in [POS_INF, NEG_INF]
        if not is_cpv:
            return expr.Integral(e.var, e.lower, self.c, e.body) + \
                   expr.Integral(e.var, self.c, e.upper, e.body)
        else:
            return Limit(x.name, POS_INF, Integral(e.var, e.lower, normalize(self.c - 1 / x, ctx), e.body) +
                         Integral(e.var, normalize(self.c + 1 / x, ctx), e.upper, e.body))


class IntegrateByEquation(Rule):
    """When the initial integral occurs in the steps."""

    def __init__(self, lhs: Union[str, Expr]):
        self.name = "IntegrateByEquation"
        if isinstance(lhs, str):
            lhs = parser.parse_expr(lhs)
        self.lhs = lhs

    def __str__(self):
        return "solve integral %s" % self.lhs

    def export(self):
        return {
            "name": self.name,
            "lhs": str(self.lhs),
            "str": str(self),
            "latex_str": "solve integral \\(%s\\)" % latex.convert_expr(self.lhs)
        }

    def eval(self, e: Expr, ctx: Context) -> Expr:
        """Eliminate the lhs's integral in rhs by solving equation."""

        def get_coeff(t: Expr, lhs: Expr) -> tuple[Expr, Expr]:
            """Rewrite t in the form a * lhs + b."""
            if t == lhs:
                return Const(1), Const(0)

            if expr.is_plus(t):
                a1, b1 = get_coeff(t.args[0], lhs)
                a2, b2 = get_coeff(t.args[1], lhs)
                return a1 + a2, b1 + b2
            elif expr.is_minus(t):
                a1, b1 = get_coeff(t.args[0], lhs)
                a2, b2 = get_coeff(t.args[1], lhs)
                return a1 - a2, b1 - b2
            elif expr.is_uminus(t):
                a, b = get_coeff(t.args[0], lhs)
                return -a, -b
            elif expr.is_times(t):
                a1, b1 = get_coeff(t.args[0], lhs)
                a2, b2 = get_coeff(t.args[1], lhs)
                if a2 != Const(0):
                    return t.args[0] * a2, t.args[0] * b2
                elif a1 != Const(0):
                    return t.args[1] * a1, t.args[1] * b1
                else:
                    return Const(0), t
            elif expr.is_divides(t):
                a1, b1 = get_coeff(t.args[0], lhs)
                return a1 / t.args[1], b1 / t.args[1]
            else:
                return Const(0), t

        # Obtain coeff with normalize
        norm_e = normalize(e, ctx)
        norm_lhs = normalize(self.lhs, ctx)
        coeff, rest = get_coeff(norm_e, norm_lhs)
        coeff = normalize(coeff, ctx)

        # Obtain coeff without normalize
        coeff2, rest2 = get_coeff(e, self.lhs)
        coeff2 = normalize(coeff2, ctx)

        if coeff == Const(0) and coeff2 == Const(0):
            raise RuleException("IntegrateByEquation", "lhs %s not found in integral" % self.lhs)

        if coeff == Const(1) or coeff2 == Const(1):
            raise RuleException("IntegrateByEquation", "lhs %s has coeff 1 in integral" % self.lhs)

        if coeff == Const(0) or coeff == Const(1):
            coeff = coeff2
            rest = rest2
        res = normalize(rest / (Const(1) - coeff), ctx)

        return res


class ElimInfInterval(Rule):
    """Convert improper integral with infinite upper or lower limits to
    a limit expression.

    If both upper and lower limits are infinity, a split point need to be
    provided.

    """

    def __init__(self, a=Const(0), new_var='t'):
        self.name = "ElimInfInterval"
        self.a = a
        self.new_var = new_var

    def __str__(self):
        return "improper integral to limit creating %s" % self.new_var

    def export(self):
        return {
            "name": self.name,
            "a": str(self.a),
            "str": str(self)
        }

    def eval(self, e: Expr, ctx: Context) -> Expr:
        def gen_lim_expr(new_var, lim, lower, upper, drt=None):
            return expr.Limit(new_var, lim, expr.Integral(e.var, lower, upper, e.body), drt)

        if not expr.is_integral(e):
            sep_ints = e.separate_integral()
            if len(sep_ints) == 0:
                return e
            else:
                return OnLocation(self, sep_ints[0][1]).eval(e, ctx)

        inf = Inf(Decimal('inf'))
        neg_inf = Inf(Decimal('-inf'))
        upper, lower = e.upper, e.lower
        new_var = self.new_var

        if upper == inf and lower != neg_inf and lower != inf:
            # INT x:[a,oo]. body => lim t->oo. INT x:[a,t]. body
            return gen_lim_expr(new_var, inf, lower, Var(new_var))
        elif upper == neg_inf and lower != neg_inf and lower != inf:
            return gen_lim_expr(new_var, inf, lower, Var(new_var))
        elif upper != inf and upper != neg_inf and lower == neg_inf:
            # INT x:[-oo,a]. body => lim t->-oo. INT x:[t,a]. body
            return gen_lim_expr(new_var, neg_inf, Var(new_var), upper)
        elif upper != inf and upper != neg_inf and lower == inf:
            return gen_lim_expr(new_var, inf, Var(new_var), upper)
        elif upper == inf and lower == neg_inf:
            # INT x:[-oo,oo]. body =>
            # lim t->-oo. INT x:[t,a]. body + lim t->oo. INT x:[a,t]. body
            assert self.a is not None, "No split point provided"
            lim1 = gen_lim_expr(new_var, neg_inf, Var(new_var), self.a)
            lim2 = gen_lim_expr(new_var, inf, self.a, Var(new_var))
            return Op('+', lim1, lim2)
        elif upper == neg_inf and lower == inf:
            assert self.a is not None, "No split point provided"
            lim1 = gen_lim_expr(new_var, inf, Var(new_var), self.a)
            lim2 = gen_lim_expr(new_var, neg_inf, self.a, Var(new_var))
            return Op('+', lim1, lim2)
        else:
            raise NotImplementedError


class LHopital(Rule):
    """Apply L'Hoptial rule."""

    def __init__(self):
        self.name = "LHopital"

    def __str__(self):
        return "l'Hopital's rule"

    def export(self):
        return {
            "name": self.name,
            "str": str(self)
        }

    def eval(self, e: Expr, ctx: Context) -> Expr:
        if not expr.is_limit(e):
            sep_lims = e.separate_limits()
            if len(sep_lims) == 0:
                return e
            else:
                return OnLocation(self, sep_lims[0][1]).eval(e, ctx)

        if not (isinstance(e.body, expr.Op) and e.body.op == '/'):
            return e

        numerator, denominator = e.body.args
        rule = DerivativeSimplify()
        return expr.Limit(e.var, e.lim, Op('/', rule.eval(Deriv(e.var, numerator), ctx),rule.eval(Deriv(e.var, denominator), ctx)), e.drt)


class DerivIntExchange(Rule):
    """Exchanging derivative and integral"""

    def __init__(self):
        self.name = "DerivIntExchange"

    def __str__(self):
        return "exchange derivative and integral"

    def export(self):
        return {
            "name": self.name,
            "str": str(self)
        }

    def eval(self, e: Expr, ctx: Context) -> Expr:
        if expr.is_deriv(e) and expr.is_integral(e.body):
            return Integral(e.body.var, e.body.lower, e.body.upper, Deriv(e.var, e.body.body))
        elif expr.is_deriv(e) and expr.is_indefinite_integral(e.body):
            return IndefiniteIntegral(e.body.var, Deriv(e.var, e.body.body), e.skolem_args)
        elif expr.is_indefinite_integral(e) and expr.is_deriv(e.body):
            return Deriv(e.body.var, IndefiniteIntegral(e.var, e.body.body, e.skolem_args))
        elif expr.is_integral(e) and expr.is_deriv(e.body):
            return Deriv(e.body.var, Integral(e.var, e.upper, e.lower, e.body.body))
        else:
            return e


class ExpandDefinition(Rule):
    """Expand a definition"""

    def __init__(self, func_name: str):
        self.name = "ExpandDefinition"
        assert isinstance(func_name, str)
        self.func_name = func_name

    def __str__(self):
        return "expand definition for %s" % self.func_name

    def export(self):
        return {
            "name": self.name,
            "func_name": self.func_name,
            "str": str(self)
        }

    @staticmethod
    def search(e: Expr, ctx: Context) -> List[Tuple[Expr, expr.Location]]:
        subexprs = e.find_subexpr_pred(lambda t: expr.is_var(t) or expr.is_fun(t))
        res = []
        for sube, loc in subexprs:
            if expr.is_fun(sube):
                for identity in ctx.get_definitions():
                    if expr.is_fun(identity.lhs) and identity.lhs.func_name == sube.func_name:
                        res.append((sube, loc))
            if expr.is_var(sube):
                for identity in ctx.get_definitions():
                    if expr.is_symbol(identity.lhs) and identity.lhs.name == sube.name:
                        res.append((sube, loc))
        return res

    def eval(self, e: Expr, ctx: Context) -> Expr:
        # Function case
        if expr.is_fun(e) and e.func_name == self.func_name:
            for identity in ctx.get_definitions():
                if expr.is_fun(identity.lhs) and identity.lhs.func_name == self.func_name:
                    inst = expr.match(e, identity.lhs)
                    if inst is None:
                        continue
                    inst_conds = [cond.inst_pat(inst) for cond in identity.conds.data]
                    if all(ctx.check_condition(cond) for cond in inst_conds):
                        return normalize(identity.rhs.inst_pat(inst), ctx)

        # Constant case
        if expr.is_var(e) and e.name == self.func_name:
            for identity in ctx.get_definitions():
                if expr.is_var(identity.lhs) and identity.lhs.name == self.func_name:
                    return identity.rhs

        # Not found
        return e


class FoldDefinition(Rule):
    """Fold a definition"""

    def __init__(self, func_name: str):
        self.name = "FoldDefinition"
        assert isinstance(func_name, str)
        self.func_name = func_name

    def __str__(self):
        return "fold definition for %s" % self.func_name

    def export(self):
        return {
            "name": self.name,
            "func_name": self.func_name,
            "str": str(self)
        }

    @staticmethod
    def search(e: Expr, ctx: Context) -> List[Tuple[Expr, expr.Location, str]]:
        subexprs = e.find_all_subexpr()
        res = []
        for sube, loc in subexprs:
            for identity in ctx.get_definitions():
                inst = expr.match(sube, identity.rhs)
                if inst:
                    if expr.is_fun(identity.lhs):
                        res.append((sube, loc, identity.lhs.func_name))
                    else:
                        res.append((sube, loc, identity.lhs.name))
        return res

    def eval(self, e: Expr, ctx: Context) -> Expr:
        for identity in ctx.get_definitions():
            if expr.is_fun(identity.lhs) and identity.lhs.func_name == self.func_name:
                inst = expr.match(e, identity.rhs)
                if inst:
                    return normalize(identity.lhs.inst_pat(inst), ctx)

            if expr.is_symbol(identity.lhs) and identity.lhs.name == self.func_name:
                if e == identity.rhs:
                    return identity.lhs

        # Not found
        return e


class IntegralEquation(Rule):
    """Integrate an equation where the left side is a derivative.

    Convert (D a. f(a)) = g(a) into f(a) = INT a. g(a). The right side
    can then be evaluated to produce a Skolem constant.

    """

    def __init__(self):
        self.name = "IntegrateBothSide"

    def eval(self, e: Expr, ctx: Context):
        assert e.is_equals() and expr.is_deriv(e.lhs)

        # Variable to differentiate, this will also be the variable
        # of integration.
        var = e.lhs.var

        # List of Skolem arguments is the free variables on the left side
        skolem_args = tuple(v for v in e.lhs.get_vars() if v != var)

        # Return f(a) = INT a. g(a)
        return Op("=", e.lhs.body, IndefiniteIntegral(var, e.rhs, skolem_args))

    def __str__(self):
        return "integrate both sides"

    def export(self):
        return {
            "name": self.name,
            "str": str(self),
        }


class SummationEquation(Rule):
    '''
    a(n) = b(n) => Sum(n, lower, upper ,a(n)) = Sum(n, lower, upper, b(n))
    '''

    def __init__(self, index_var: str, lower: Union[Expr, str], upper: Union[Expr, str]):
        if isinstance(lower, str):
            lower = parser.parse_expr(lower)
        if isinstance(upper, str):
            upper = parser.parse_expr(upper)
        self.name = "SummationEquation"
        self.index_var = index_var
        self.lower = lower
        self.upper = upper

    def eval(self, e: Expr, ctx: Context):
        assert e.is_equals()
        e1 = Summation(self.index_var, self.lower, self.upper, e.lhs)
        e2 = Summation(self.index_var, self.lower, self.upper, e.rhs)
        return Op("=", e1, e2)

    def __str__(self):
        return "sum both side"

    def export(self):
        return {
            "name": self.name,
            "str": str(self),
            'index_var': self.index_var,
            'lower': str(self.lower),
            'upper': str(self.upper)
        }


class ChangeSummationIndex(Rule):
    '''
    sum(n, 1, oo, a(n)) => sum(n, 0, oo, a(n+1))
    '''

    def __init__(self, new_lower: Union[Expr, str]):
        self.name = "ChangeSummationIndex"
        self.new_lower = new_lower if isinstance(new_lower, Expr) else parser.parse_expr(new_lower)

    def eval(self, e: Expr, ctx: Context):
        if not expr.is_summation(e):
            return e
        tmp = normalize(Var(e.index_var) + e.lower - self.new_lower, ctx)
        new_upper = normalize(e.upper + self.new_lower - e.lower, ctx) \
            if e.upper != POS_INF else POS_INF
        return Summation(e.index_var, self.new_lower, new_upper,
                         e.body.replace(Var(e.index_var), tmp))

    def __str__(self):
        return "change summation index"

    def export(self):
        return {
            "name": self.name,
            "str": str(self),
            "new_lower": str(self.new_lower)
        }


class LimitEquation(Rule):
    """Apply limit to both sides of the equation.

        A = B -> LIM {x -> a}. A = LIM {x -> a}. B

    """

    def __init__(self, var: str, lim: Expr):
        self.name = "LimitEquation"
        self.var = var
        self.lim = lim

    def __str__(self):
        return "apply limit %s -> %s both sides" % (self.var, self.lim)

    def eval(self, e: Expr, ctx: Context):
        v, lim = self.var, self.lim
        lim1 = Limit(v, lim, e.lhs)
        lim2 = Limit(v, lim, e.rhs)
        return Op('=', lim1, lim2)

    def export(self):
        return {
            "name": self.name,
            "str": str(self),
            "var": self.var,
            "lim": str(self.lim),
            "latex_str": "apply limit \\(%s \\to %s\\) both sides" %
                         (self.var, latex.convert_expr(self.lim))
        }


class IntSumExchange(Rule):
    """Exchange integral and summation"""

    def __init__(self):
        self.name = "IntSumExchange"

    def __str__(self):
        return "exchange integral and sum"

    # def test_converge(self, svar, sl, su, ivar, il, iu, body, ctx: Context):
    #     if ctx.is_not_negative(body):
    #         return True
    #     if ctx.is_not_positive(body):
    #         return True
    #     if su != expr.POS_INF:
    #         return True

    #     abs_body = normalize(Fun("abs", body), ctx)
    #     goal1 = Fun("converges", Summation(svar, sl, su, Integral(ivar, il, iu, abs_body)))

    #     abs_int = normalize(Fun("abs", Integral(ivar, il, iu, body)), ctx)
    #     goal2 = Fun("converges", Summation(svar, sl, su, abs_int))

    #     for lemma in ctx.get_lemmas():
    #         if normalize(lemma.expr, ctx) == normalize(goal1, ctx):
    #             return True
    #         if normalize(lemma.expr, ctx) == normalize(goal2, ctx):
    #             return True
    #     for _, subgoal in ctx.get_all_subgoals().items():
    #         if normalize(subgoal.expr, ctx) == normalize(goal1, ctx):
    #             return True
    #         if normalize(subgoal.expr, ctx) == normalize(goal2, ctx):
    #             return True
    #     return False

    def eval(self, e: Expr, ctx: Context):
        if expr.is_integral(e) and expr.is_summation(e.body):
            ctx2 = body_conds(e, body_conds(e.body, ctx))
            s = e.body
            # if self.test_converge(s.index_var, s.lower, s.upper, e.var, e.lower, e.upper, e.body.body, ctx2):
            return Summation(s.index_var, s.lower, s.upper, Integral(e.var, e.lower, e.upper, s.body))
        elif expr.is_summation(e) and expr.is_integral(e.body):
            ctx2 = body_conds(e, body_conds(e.body, ctx))
            i = e.body
            # if self.test_converge(e.index_var, e.lower, e.upper, i.var, i.lower, i.upper, e.body.body, ctx2):
            return Integral(i.var, i.lower, i.upper, Summation(e.index_var, e.lower, e.upper, i.body))
        return e

    def export(self):
        return {
            "name": self.name,
            "str": str(self)
        }

def has_negative_coefficient(expr_obj: Expr, var_name: str) -> bool:
    """
    检查表达式中指定变量是否带负系数
    例如：2-y 中 y 带负号，y+2 中 y 不带负号
    """
    from integral.expr import Op, Var
    
    # 如果表达式就是变量本身，系数为正
    if expr.is_var(expr_obj) and expr_obj.name == var_name:
        return False
    
    # 如果是加法 a + b
    if expr.is_op(expr_obj) and expr_obj.op == '+':
        for arg in expr_obj.args:
            if expr.is_var(arg) and arg.name == var_name:
                return False  # 直接出现的变量，系数为正
            if expr.is_uminus(arg) and len(arg.args) == 1:
                if expr.is_var(arg.args[0]) and arg.args[0].name == var_name:
                    return True  # -var
        return has_negative_coefficient(expr_obj.args[0], var_name) or has_negative_coefficient(expr_obj.args[1], var_name)
    
    # 如果是减法 a - b
    if expr.is_op(expr_obj) and expr_obj.op == '-' and len(expr_obj.args) == 2:
        lhs, rhs = expr_obj.args
        # 检查右侧是否包含变量（被减去，所以是负的）
        if expr.is_var(rhs) and rhs.name == var_name:
            return True
        if rhs.contains_var(var_name):
            return True  # 简化处理：只要在减号右侧，就认为是负的
    
    # 如果是一元负号 -a
    if expr.is_uminus(expr_obj):
        if len(expr_obj.args) == 1 and expr.is_var(expr_obj.args[0]) and expr_obj.args[0].name == var_name:
            return True
    
    # 默认情况
    return False

class IntExchange(Rule):
    """Exchange integral and integral"""

    def __init__(self):
        self.name = "IntExchange"

    def __str__(self):
        return "exchange integral and integral"

    def exchange_int(self, evar, el, eu, svar, sl, su, sb, ctx: Context) -> Expr:
        """交换积分次序的核心算法
        
        参数:
            evar, svar: str - 外层和内层积分变量名
            el, eu: Expr - 外层积分上下限
            sl, su: Expr - 内层积分上下限
            sb: Expr - 被积函数
            ctx: Context - 上下文
        """
        res_list = []
        simplify = Simplify()
        # 初始化不等式列表和边界点数组
        inequality_list = []
        bp = []
        
        # 在函数入口统一处理类型：将字符串变量名转换为 Var 对象
        # 这样后续逻辑可以直接使用，避免到处调用 exprify
        evar_obj = Var(evar)  # 外层变量的 Var 对象
        svar_obj = Var(svar)  # 内层变量的 Var 对象

        # 第一步：构建不等式列表
        # 添加外层积分上下限的不等式
        inequality_list.append(Op('<', evar_obj, eu))
        inequality_list.append(Op('>', evar_obj, el))

        # 添加内层积分上下限的不等式
        # 处理 sl < svar < su
        if is_const(sl):
            if Op('<', sl, evar_obj) and Op('>', evar_obj, sl) not in inequality_list:
                inequality_list.append(Op('<', sl, evar_obj))
        else:
            from integral.solve import solve_equation
            # 反解下限：sl(evar) = svar，求 evar
            solutions = solve_equation(sl, svar_obj, evar, ctx)
            if not solutions:
                raise RuleException("IntExchange", f"Unable to solve {sl} = {svar} for {evar}")
            evar_solution = solutions[0]  # 取第一个解
            
            # 判断原始下限表达式 sl 中 evar 的符号，而不是反解后的符号
            # 从 sl(evar) < svar 推导 evar 的约束
            if has_negative_coefficient(sl, evar):
                # sl 关于 evar 递减（如 sl = -y+2），则 sl < svar => f(svar) < evar
                inequality_list.append(Op('<', evar_solution, evar_obj))
            else:
                # sl 关于 evar 递增（如 sl = y），则 sl < svar => evar < svar，即 evar < f(svar)
                inequality_list.append(Op('<', evar_obj, evar_solution))

        if is_const(su):
            if Op('<', evar_obj, su) and Op('>', evar_obj, su) not in inequality_list:
                inequality_list.append(Op('<', evar_obj, su))
        else:
            from integral.solve import solve_equation
            # 反解上限：su(evar) = svar，求 evar
            solutions = solve_equation(su, svar_obj, evar, ctx)
            if not solutions:
                raise RuleException("IntExchange", f"Unable to solve {su} = {svar} for {evar}")
            evar_solution = solutions[0]  # 取第一个解
            
            # 判断原始上限表达式 su 中 evar 的符号，而不是反解后的符号
            # 从 svar < su(evar) 推导 evar 的约束
            if has_negative_coefficient(su, evar):
                # su 关于 evar 递减（如 su = 2-y），则 svar < su => evar < f(svar)
                inequality_list.append(Op('<', evar_obj, evar_solution))
            else:
                # su 关于 evar 递增（如 su = y+2），则 svar < su => f(svar) < evar
                inequality_list.append(Op('<', evar_solution, evar_obj))

        points_to_check = [
            simplify.eval(sl.subst(evar, el), ctx),
            simplify.eval(sl.subst(evar, eu), ctx),
            simplify.eval(su.subst(evar, el), ctx),
            simplify.eval(su.subst(evar, eu), ctx)
        ]
        for point in points_to_check:
            if point not in bp:  # 确保不重复添加
                bp.append(point)
        # 排序边界点
        bp = sorted(bp)
        # 第三步：生成新的积分上下限
        # 遍历bp计算积分式
        for i in range(len(bp) - 1):
            # 外层积分上下限
            new_el = bp[i]
            new_eu = bp[i + 1]

            # 计算中位数 mid
            mid = simplify.eval((new_el + new_eu) / 2, ctx)

            # 求解内层积分上下限
            filtered_inequalities = []  # 用于存储mid替换后的不等式结果
            for condition in inequality_list:
                op = condition.op
                lhs = condition.args[0]
                rhs = condition.args[1]

                # 替换svar并检查条件
                if lhs == evar_obj:  # 左侧包含evar
                    if is_const(rhs):
                        filtered_inequalities.append((Op(op,lhs,rhs), Op(op,lhs,rhs)))
                    else:
                        # 将 svar 替换成 mid（rhs 已经是 Expr，不需要转换）
                        evaluated_rhs = simplify.eval(rhs.subst(svar, mid), ctx)
                        # 创建新的不等式
                        new_condition = Op(op, evar_obj, evaluated_rhs)
                        filtered_inequalities.append((new_condition, Op(op,lhs,rhs)))
                elif rhs == evar_obj:  # 右侧包含evar
                    if is_const(lhs):
                        filtered_inequalities.append((Op(op,lhs,rhs), Op(op,lhs,rhs)))
                    else:
                        # 将 svar 替换成 mid（lhs 已经是 Expr，不需要转换）
                        evaluated_lhs = simplify.eval(lhs.subst(svar, mid), ctx)
                        # 创建新的不等式
                        new_condition = Op(op, evaluated_lhs, evar_obj)
                        filtered_inequalities.append((new_condition, Op(op,lhs,rhs)))

            # 取交集，找到满足条件的上下限
            lower_bound = None
            upper_bound = None
            lower_source = None
            upper_source = None

            for new_condition, original_condition in filtered_inequalities:
                op = new_condition.op
                lhs = new_condition.args[0]
                rhs = new_condition.args[1]
                if op == "<":
                    if lhs == evar_obj:
                        # 更新上限
                        if upper_bound is None or rhs < upper_bound:
                            upper_bound = rhs
                            upper_source = original_condition.args[1]  # 保存原始不等式
                    else:
                        # 更新下限
                        if lower_bound is None or lhs > lower_bound:
                            lower_bound = lhs
                            lower_source = original_condition.args[0]  # 保存原始不等式
                elif op == ">":
                    if lhs == evar_obj:
                        # 更新下限
                        if lower_bound is None or rhs > lower_bound:
                            lower_bound = rhs
                            lower_source = original_condition.args[1]  # 保存原始不等式
                    else:
                        # 更新上限
                        if upper_bound is None or lhs < upper_bound:
                            upper_bound = lhs
                            upper_source = original_condition.args[0]  # 保存原始不等式

            res_list.append(normalize(Integral(svar,new_el,new_eu,Integral(evar,lower_source,upper_source,sb)),ctx))

        res = res_list[0]
        if len(res_list) == 1:
            return res_list[0]
        else:
            for arg in res_list[1:]:
                res = Op('+', arg,res)
        return res

    def judge_novar(self,value,evar = None,svar = None):
        # judge whether contain oo or algebraic constant and no var
        value = str(value)
        if value == 'oo':
            return True
        elif evar is None and svar is None:
            return True
        elif evar is None and svar in value:
            return True
        elif evar is not None and svar is not None and \
                evar not in value and svar not in value:
            return True
        else:
            return False

    def judge_contains_var_and_letters(self,value,evar = None,svar = None):
        # judge whether contain var and algebraic constant
        value = str(value)
        contains_svar_or_evar = (svar is not None and svar in value) or (evar is not None and evar in value)
        contains_other_letters = any(c.isalpha() for c in value if c not in {svar, evar})

        if contains_svar_or_evar and contains_other_letters:
            return True
        else:
            return False

    def eval(self, e: Expr, ctx: Context):
        # Check if e is an integral
        if not expr.is_integral(e):
            raise RuleException("IntExchange", 
                              "can only be applied to integrals, got: %s" % type(e).__name__)
        
        # Check if it's a double integral
        if not expr.is_integral(e.body):
            raise RuleException("IntExchange", 
                              "can only be applied to double integrals, got single integral: %s" % e)
        
        # Now we have a double integral, proceed with the exchange
        ctx2 = body_conds(e, body_conds(e.body, ctx))
        s = e.body
        if not self.judge_novar(e.upper,e.var,s.var) or not self.judge_novar(e.lower,e.var,s.var) or \
                self.judge_novar(s.upper,None,s.var) or self.judge_novar(s.lower,None,s.var):
            raise RuleException("IntExchange", "Integral format error: outer limits must be constants or contain only outer variable, inner limits must contain inner variable")
        elif self.judge_contains_var_and_letters(s.upper,e.var,None) or self.judge_contains_var_and_letters(s.lower,e.var,None):
            raise RuleException("IntExchange", "Contain algebraic constant and var")
        # judge whether contain oo or Algebraic Constant and no var
        elif self.judge_novar(e.upper) and self.judge_novar(e.lower) and \
                self.judge_novar(s.upper,e.var,s.var) and self.judge_novar(s.lower,e.var,s.var):
            return Integral(s.var, s.lower, s.upper, Integral(e.var, e.lower, e.upper, s.body))
        else:
            return self.exchange_int(e.var, e.lower, e.upper, s.var, s.lower, s.upper, s.body,ctx2)

    def export(self):
        return {
            "name": self.name,
            "str": str(self)
        }


class VarSubsOfEquation(Rule):
    """Substitute variable for any expression in an equation.
    """

    def __init__(self, subst: Dict[str, Union[str, Expr]]):
        self.name = "VarSubsOfEquation"
        for i in range(len(subst)):
            if isinstance(subst[i]['expr'], str):
                if subst[i]['expr'] == "":
                    subst[i]['expr'] = None
                else:
                    subst[i]['expr'] = parser.parse_expr(subst[i]['expr'])
        self.subst = subst

    def __str__(self):
        str_of_substs = ', '.join(item['var'] + " for " + str(item['expr']) for item in self.subst
                                  if item['expr'] is not None)
        return "substitute " + str_of_substs + " in equation"

    def export(self):
        latex_str_of_substs = ', '.join('\\(' + item['var'] + "\\) for \\(" + latex.convert_expr(item['expr']) + '\\)'
                                        for item in self.subst if item['expr'] is not None)
        json_substs = list()
        for item in self.subst:
            json_substs.append({'var': item['var'], 'expr': str(item['expr'])})
        return {
            "name": self.name,
            "str": str(self),
            "subst": json_substs,
            "latex_str": "substitute %s in equation" % latex_str_of_substs
        }

    def eval(self, e: Expr, ctx: Context) -> Expr:
        if e.is_equals():
            for item in self.subst:
                if item['expr'] is not None:
                    e = e.subst(item['var'], item['expr'])
            return poly.normal_const(e, ctx)
        else:
            return e


class MergeSummation(Rule):
    "SUM(u,0,oo, body1) + SUM(k,0,oo,body2) = SUM(u, 0, oo, body1+body2)"

    def __init__(self):
        self.name = "MergeSummation"

    def __str__(self):
        return "merge summation"

    def export(self):
        return {
            "name": self.name,
            "str": str(self)
        }

    def eval(self, e: Expr, ctx: Context) -> Expr:
        if not (e.ty == OP and e.op in ('+', '-') and all([isinstance(arg, Summation) for arg in e.args])):
            return e
        a, b = e.args
        if not (a.lower == b.lower and a.upper == b.upper):
            return e
        if a.index_var != b.index_var:
            b = b.alpha_convert(a.index_var)
        return Summation(a.index_var, a.lower, a.upper, Op(e.op, a.body, b.body))


class DerivEquation(Rule):
    """Differentiate both sides with respect to some variable."""

    def __init__(self, var: str):
        self.name = "DerivEquation"
        self.var = var

    def __str__(self):
        return "differentiate both sides at %s" % self.var

    def export(self):
        return {
            "name": self.name,
            "str": str(self),
            "var": self.var,
            "latex_str": "differentiate both sides at \\(%s\\)" % self.var
        }

    def eval(self, e: Expr, ctx: Context) -> Expr:
        if not e.is_equals():
            return e
        return Op('=', Deriv(self.var, e.lhs), Deriv(self.var, e.rhs))


class SolveEquation(Rule):
    """Solve equation for the given expression."""

    def __init__(self, solve_for: Union[Expr, str]):
        if isinstance(solve_for, str):
            solve_for = parser.parse_expr(solve_for)
        self.solve_for = solve_for
        self.name = "SolveEquation"

    def eval(self, e: Expr, ctx: Context):
        assert e.is_equals()

        res = solve_for_term(e, self.solve_for, ctx)
        if not res:
            raise RuleException("SolveEquation", f"cannot solve for {self.solve_for} in {e}")
        return Op("=", self.solve_for, normalize(res, ctx))

    def __str__(self):
        return "solve equation for %s" % str(self.solve_for)

    def export(self):
        return {
            "name": self.name,
            "str": str(self),
            "solve_for": str(self.solve_for),
            "latex_str": "solve equation for \\(%s\\)" % latex.convert_expr(self.solve_for)
        }


class FunEquation(Rule):
    """a = b => fun(a) = fun(b) if a,b belong to the domain of f"""

    def __init__(self, func_name: str):
        self.name = "FunEquation"
        self.func_name: str = func_name

    def __str__(self):
        return "function on both sides"

    def export(self):
        return {
            "name": self.name,
            "str": str(self),
            "func_name": self.func_name
        }

    def eval(self, e: Expr, ctx: Context) -> Expr:
        if not e.is_equals():
            return e
        ne = Op('=', Fun(self.func_name, e.lhs), Fun(self.func_name, e.rhs))
        return ne


class LimRewrite(Rule):
    """Expand functions on matrices."""

    def __init__(self, source: Expr, target: Expr):
        self.name = "LimRewrite"
        self.source = source
        self.target = target

    def __str__(self):
        return "rewrite limit expression"

    def export(self):
        return {
            "name": self.name,
            "str": str(self)
        }

    def eval(self, e: Expr, ctx: Context) -> Expr:
        if self.source != e:
            find_res = e.find_subexpr(self.source)
            if len(find_res) == 0:
                raise AssertionError("LimRewrite: source expression not found")
            loc = find_res[0]
            return OnLocation(self, loc).eval(e, ctx)
        assert self.source == e
        if expr.is_limit(e):
            e: Limit
            b = e.body
            res = None
            # TODO: check whether limit values exist
            if expr.is_op(b):
                if len(b.args) == 2:
                    if b.op in "+-*/":
                        res = Op(b.op, Limit(e.var, e.lim, b.args[0], e.drt),
                                 Limit(e.var, e.lim, b.args[1], e.drt))
                    if b.op == '^':
                        if not b.args[1].contains_var(e.var):
                            res = Op(b.op, Limit(e.var, e.lim, b.args[0], e.drt), b.args[1])
                elif len(b.args) == 1:
                    if b.op == '-':
                        res = -Limit(e.var, e.lim, b.args[0], e.drt)
            if res != None and normalize(res, ctx) == normalize(self.target, ctx):
                return self.target
        return e

class ResidueTheorem(Rule):
    """Apply residue theorem.
    应用留数定理计算复积分
    围道积分值等于2πi乘以围道内极点的留数之和
    ∮ f(z) dz = 2πi * ∑(n(γ,z0) * Res(f,z0))
    
    其中
    - n(γ,z0)是围道γ绕极点z0的绕数
    - Res(f,z0)是函数f在z0处的留数
    
    Args:
        e: 输入表达式
        ctx: 上下文
        
    Returns:
        计算结果
        
    Raises:
        RuleException: 当应用留数定理出现错误时抛出异常
    """
    def __init__(self):
        self.name = "residue theorem"

    def __str__(self):
        return "apply residue theorem"

    def export(self):
        return {"name": "residue theorem"}
    
    def _lookup_path_definition(self, path_ref: str, ctx: Context) -> CINTPath:
        """从上下文中查找路径定义
        
        路径可以用define命令定义
        
        Args:
            path_ref: 路径引用字符串，如 "C(t,1)"
            ctx: 上下文
            
        Returns:
            CINTPath对象，如果找不到则返回None
        """
        from integral.expr import Fun, is_fun, CINTPath
        from integral.parser import parse_expr
        
        # 解析路径引用（例如 "C(t,r)"）
        try:
            path_call = parse_expr(path_ref)
            if not is_fun(path_call):
                return None
                
            func_name = path_call.func_name
            args = path_call.args
            
            # 在上下文中查找函数定义
            for definition in ctx.get_definitions():
                if (hasattr(definition.lhs, 'func_name') and 
                    definition.lhs.func_name == func_name and
                    len(definition.lhs.args) == len(args)):
                    
                    # 应用参数替换
                    path_def = definition.rhs
                    
                    # 确保path_def是CINTPath
                    if isinstance(path_def, CINTPath):
                        # 应用参数替换
                        new_path_expr = path_def.path_expr
                        new_start_expr = path_def.start_expr
                        new_end_expr = path_def.end_expr
                        
                        for param, arg in zip(definition.lhs.args, args):
                            # 替换参数
                            new_path_expr = new_path_expr.replace(param, arg)
                            new_start_expr = new_start_expr.replace(param, arg)
                            new_end_expr = new_end_expr.replace(param, arg)
                        
                        return CINTPath(path_def.var, new_path_expr, new_start_expr, new_end_expr)
                    
        except Exception:
            pass
        
        return None
    

    
    def eval(self, e: Expr, ctx: Context) -> Expr:

        # 检查输入是否是围道积分或包含围道积分的极限
        cintegral = None
        limit_var = None
        limit_value = None
        
        if expr.is_cintegral(e):
            cintegral = e
        elif expr.is_limit(e) and expr.is_cintegral(e.body):
            cintegral = e.body
            # 提取极限信息
            limit_var = e.var
            limit_value = e.lim
        else:
            raise RuleException("ResidueTheorem", "input is not a complex integral expression.")
        
        # 创建极限信息字典，用于传递给需要的函数
        limit_info = {}
        if limit_var and limit_value:
            limit_info[limit_var] = limit_value
        
        temp_ctx = ctx
        
        # 使用找到的围道积分进行后续处理
        e = cintegral

        # 获取积分路径和被积函数
        paths = e.paths
        f = e.body
        
        # 验证路径类型
        if not paths:
            raise RuleException("ResidueTheorem", "No contour path provided")
        
        # 转换所有路径为 CINTPath 对象
        resolved_paths = []
        for path in paths:
            if isinstance(path, CINTPath):
                # 已经是 CINTPath 对象，直接使用
                resolved_paths.append(path)
            elif isinstance(path, str):
                # 字符串引用，从上下文中查找定义
                path_obj = self._lookup_path_definition(path, temp_ctx)
                if path_obj is None:
                    raise RuleException("ResidueTheorem", 
                        f"Path definition '{path}' not found in context")
                resolved_paths.append(path_obj)
            else:
                # 其他类型，报错
                raise RuleException("ResidueTheorem", 
                    f"Invalid path type: {type(path).__name__}, expected CINTPath or string reference")
        
        # 使用解析后的路径列表
        paths = resolved_paths
        
        # 验证围道闭合性
        # 检查围道是否闭合（单路径或多路径）
        if not is_closed_contour(paths):
            raise RuleException("ResidueTheorem", 
                "Contour is not closed: paths do not form a closed loop")
        
        # 寻找围道内部的极点
        poles = find_poles_inside_contour(f, e.var, paths, temp_ctx, limit_info)
        
        # 验证所有极点都在围道内
        if not poles:
            return Const(0)
        
        # 检查是否有极点的绕数为0（这意味着极点不在围道内）
        poles_outside = [(pole, order, wind_num) for pole, order, wind_num in poles if wind_num == 0]
        if poles_outside:
            pole_names = [str(pole) for pole, _, _ in poles_outside]
            raise RuleException("ResidueTheorem", f"Poles {', '.join(pole_names)} are outside the contour, cannot apply residue theorem")
            
        # 计算每个极点的留数，并乘以其绕数，然后求和
        result = Const(0)
        for pole, order, wind_num in poles:
            # 计算留数，使用从 find_poles 获取的实际阶数
            # 对于简单极点(order=1): lim(z→pole)[f(z)*(z-pole)]
            # 对于高阶极点(order>1): 使用导数公式
            
            # 使用围道积分的变量名和实际极点阶数进行留数计算
            residue = normalize(expr.compute_residue(f, pole, order, e.var), temp_ctx)
            
            # 计算 2πi * n(极点,围道) * Res(f,极点)
            term = normalize(Op("*", Const(wind_num), residue), temp_ctx)
            result = normalize(Op("+", result, term), temp_ctx)
            
        # 乘以 2πi 系数
        result = normalize(Op("*", Op("*", Const(2), Fun("pi")), Op("*", Fun("i"), result)), temp_ctx)
            
        return result

# 添加全局缓存以提高性能
_poles_cache = {}
_winding_cache = {}

def find_poles_inside_contour(func: Expr, var: str, paths: List[Union[CINTPath]], ctx: Context, limit_info: dict = None) -> List[Tuple[Expr, int, int]]:
    """查找位于闭合围道内部的极点及其阶数
    
    检查点是否在闭合围道内部
    步骤：
    1. 查找所有极点及其阶数
    2. 计算每个极点的绕数
    3. 返回绕数非零的极点（在围道内）
    
    Args:
        func: 复变函数
        var: 变量名称
        paths: 闭合围道路径列表（已展开的 CINTPath 对象）
        ctx: 上下文
        limit_info: 极限信息字典
        
    Returns:
        List[Tuple[Expr, int, int]]: 围道内极点列表，每个元素为(极点, 阶数, 绕数)的元组
    """
    from integral.expr import find_poles
    
    # 检查围道是否闭合
    is_closed = is_closed_contour(paths)
    if not is_closed:
        return []
    
    # 使用缓存查找极点（包含阶数）
    cache_key = (hash(func), var)
    if cache_key in _poles_cache:
        poles_with_order = _poles_cache[cache_key]
    else:
        # 查找函数的所有极点及其阶数
        poles_with_order = find_poles(var, func, ctx)
        _poles_cache[cache_key] = poles_with_order
    
    # 检查每个极点是否在围道内部
    result = []
    for pole, order in poles_with_order:
        wind_num = winding_number(pole, paths, ctx, limit_info if limit_info else {})
        if wind_num != 0:
            result.append((pole, order, wind_num))
    
    return result

def compute_path_endpoints(path: CINTPath, ctx: Context) -> tuple[Expr, Expr]:
    """计算路径的起点和终点表达式
    
    Args:
        path: 参数化路径
        ctx: 上下文
        
    Returns:
        (start_point, end_point): 起点和终点的表达式
    """
    from integral.poly import normalize
    
    # 计算起点
    try:
        start_point = normalize(path.path_expr.subst(path.var, path.start_expr), ctx)
    except:
        start_point = path.path_expr.subst(path.var, path.start_expr)
    
    # 计算终点
    try:
        end_point = normalize(path.path_expr.subst(path.var, path.end_expr), ctx)
    except:
        end_point = path.path_expr.subst(path.var, path.end_expr)
    
    return start_point, end_point


def is_closed_contour(paths: List[Union[CINTPath]]) -> bool:
    """检查围道是否闭合

    对于单一路径，检查其是否自然封闭（如完整的圆）。
    对于多条路径组成的围道，验证它们是否首尾相连形成闭合回路。

    
    Args:
        paths: 路径列表(CINTPath)

    Returns:
        bool: 如果围道闭合返回True，否则返回False
    """
    if not paths:
        return False
    ctx = Context()

    # 单一路径情况
    if len(paths) == 1:
        path = paths[0]
        if isinstance(path, CINTPath):
            try:
                start_point, end_point = compute_path_endpoints(path, ctx)
                return are_equal(start_point, end_point)
            except:
                # 如果无法计算，使用简单的参数闭合判断
                return path.is_closed()

    # 多路径情况需检查首尾连接
    endpoints = []
    for path in paths:
        if isinstance(path, CINTPath):
            try:
                start_point, end_point = compute_path_endpoints(path, ctx)
                endpoints.append((start_point, end_point))
            except:
                # 如果无法计算，假设不闭合
                return False

    # 检查端点是否首尾相连形成回路
    for i in range(len(endpoints)):
        _, end = endpoints[i]
        start_next, _ = endpoints[(i + 1) % len(endpoints)]
        if not are_equal(end, start_next):
            return False

    return True

def are_equal(expr1, expr2):
    """比较两个数学表达式是否相等。

    首先使用符号简化，若失败则尝试数值评估。

    参数：
        expr1: 第一个表达式。
        expr2: 第二个表达式。

    返回：
        bool: 如果表达式相等，返回 True;否则返回 False。
    """
    # 符号比较
    diff = Op("-", expr1, expr2)
    simplified = normalize(diff, Context())
    if simplified == expr.Const(0):
        return True

    # 数值比较作为备用（需要先替换 pi）
    try:
        # 替换 pi 为数值
        from decimal import Decimal
        expr1_with_pi = expr1.subst('pi', Const(Decimal(str(math.pi))))
        expr2_with_pi = expr2.subst('pi', Const(Decimal(str(math.pi))))
        
        val1 = expr.eval_expr(expr1_with_pi)
        val2 = expr.eval_expr(expr2_with_pi)
        
        # 处理复数比较
        if isinstance(val1, complex) and isinstance(val2, complex):
            return abs(val1 - val2) < 1e-9
        elif isinstance(val1, complex) or isinstance(val2, complex):
            # 转换为复数再比较
            c1 = complex(val1) if not isinstance(val1, complex) else val1
            c2 = complex(val2) if not isinstance(val2, complex) else val2
            return abs(c1 - c2) < 1e-9
        else:
            return abs(val1 - val2) < 1e-9
    except:
        return False  # 如果评估失败，保守返回 False

def winding_number(point: Expr, paths: List[Union[CINTPath]], ctx: Context, limit_info: dict = None) -> int:
    """计算绕数
    
    使用Cauchy指数计算绕数:n(γ,z0) = -Indp(γ,z0)/2
    Indp(γ,z0)是路径γ绕点z0的Cauchy指数
    
    n(γ,z0)表示闭合路径围绕点z0时旋转的整数圈数，正值表示逆时针方向
    负值表示顺时针方向，0表示路径不包含该点
    计算"绕数"为复变函数的基础
    
    Args:
        point: 要计算绕数的点z
        paths: 闭合围道路径列表（已展开的 CINTPath 对象）
        ctx: 上下文
        limit_info: 极限信息字典，如 {'r': POS_INF}
        
    Returns:
        int: 绕数值，表示闭合路径围绕点旋转的整数圈数
    """
    if limit_info is None:
        limit_info = {}
    
    # 创建缓存键
    paths_hash = tuple(hash(p) if isinstance(p, CINTPath) else hash(str(p)) for p in paths)
    limit_hash = tuple(sorted((k, hash(v)) for k, v in limit_info.items()))
    cache_key = (hash(point), paths_hash, limit_hash)
    
    # 检查缓存
    if cache_key in _winding_cache:
        return _winding_cache[cache_key]
        
    # 检查围道是否闭合，不闭合则返回0
    if not is_closed_contour(paths):
        _winding_cache[cache_key] = 0
        return 0
    
    # 收集所有路径的所有跳变值
    all_jump_values = []
    
    for path in paths:
        # 获取每条路径的跳变值列表（不除以2）
        jump_vals = compute_jump_values(path, point, ctx, limit_info)
        all_jump_values.extend(jump_vals)
    
    # 所有跳变值求和后除以2得到绕数
    total_jump = sum(all_jump_values)
    winding = total_jump / 2.0
    
    # 对于闭合回路，绕数应该是整数，四舍五入
    result = int(round(winding))
    
    # 缓存结果
    _winding_cache[cache_key] = result
    
    return result

def compute_jump_values(path: Union[CINTPath], point: Expr, ctx: Context, limit_info: dict = None) -> list[float]:
    """计算路径的跳变值列表（不除以2）
    
    Args:
        path: CINTPath参数化路径对象 γ(t)
        point: 点 z
        ctx: 上下文
        limit_info: 极限信息字典，如 {'r': POS_INF}
        
    Returns:
        跳变值列表
    """
    from integral.poly import normalize
    
    if limit_info is None:
        limit_info = {}
    
    # 1. 展开路径表达式并提取实部虚部
    path_expr = ExpandPolynomial().eval(expand_euler(path.path_expr, ctx), ctx)
    path_expr = normalize(path_expr, ctx)
    
    # 2. 提取 γ(t) 的实部和虚部（含参数 t）
    gamma_re, gamma_im = _extract_complex_parts(path_expr, ctx)
    gamma_re = normalize(gamma_re, ctx)
    gamma_im = normalize(gamma_im, ctx)
    if gamma_re is None or gamma_im is None:
        return []
    
    # 3. 提取极点 z₀ 的实部和虚部（常数）
    point_re, point_im = _extract_complex_parts(point, ctx)
    point_re = normalize(point_re, ctx)
    point_im = normalize(point_im, ctx)
    if point_re is None or point_im is None:
        return []
    
    # 4. 构造 f(t) = Im(γ(t)-z₀)/Re(γ(t)-z₀)
    # 已经得到 diff_re 和 diff_im
    diff_re = Op("-", gamma_re, point_re) if point_re != Const(0) else gamma_re
    diff_im = Op("-", gamma_im, point_im) if point_im != Const(0) else gamma_im
    
    # 标准化 diff_re 和 diff_im
    diff_re = normalize(diff_re, ctx)
    diff_im = normalize(diff_im, ctx)
    
    # 5. 获取参数范围
    start_expr = path.start_expr
    end_expr = path.end_expr
    
    # 6. 寻找跳跃点：符号求解 diff_re = 0 且验证 Re(γ(sol)) = 0 的点
    jump_points = find_jump_points(diff_re, path.var, start_expr, end_expr, ctx)
    
    # 收集所有跳变值
    jump_values = []
    
    # 对每个跳跃点计算跳变值
    for jump_pt in jump_points:
        # 将跳跃点带入路径参数方程
        res = path.path_expr.subst(path.var, jump_pt)
        
        # 替换未绑定变量为其极限值（使用limit_info）
        for var_name, var_limit in limit_info.items():
            res = res.subst(var_name, var_limit)
        
        # 化简表达式
        # 先normalize（包括计算oo*0等）
        res = normalize(res, ctx)
        
        # 处理：coeff * exp(i * theta) 形式（coeff可以是任意常数或oo）
        # 这种情况下_extract_complex_parts无法正确处理，需要手动提取虚部
        if isinstance(res, Op) and res.op == '*' and len(res.args) == 2:
            coeff, exp_part = res.args
            # 检查是否是 coeff * exp(i*theta) 形式（coeff可以是常数或oo）
            if isinstance(exp_part, Fun) and exp_part.func_name == 'exp':
                exp_arg = exp_part.args[0]  # i*theta
                
                # 检查exp_arg中是否包含i
                def contains_i_factor(e):
                    if isinstance(e, Fun) and e.func_name == 'i':
                        return True
                    if isinstance(e, Op) and e.op in ['*', '/']:
                        return any(contains_i_factor(arg) for arg in e.args)
                    return False
                
                if contains_i_factor(exp_arg):
                    # coeff * exp(i*theta) = coeff * (cos(theta) + i*sin(theta))
                    # 展开exp得到cos和sin
                    exp_expanded = expand_euler(exp_part, ctx)
                    
                    # 如果展开失败，尝试手动构建
                    if exp_expanded == exp_part:
                        # expand_euler没有展开，简化假设：
                        # coeff * exp(i*theta) 的虚部约等于 coeff（当sin(theta)≈1时）
                        res_re = Const(0) 
                        res_im = coeff
                    else:
                        # 展开成功，提取虚部
                        _, exp_im = _extract_complex_parts(exp_expanded, ctx)
                        res_re = Const(0) if exp_im != Const(0) else exp_expanded
                        res_im = Op('*', coeff, exp_im) if exp_im != Const(0) and exp_im != Const(1) else coeff if exp_im != Const(0) else Const(0)
                        res_im = normalize(res_im, ctx)
                else:
                    res_re, res_im = _extract_complex_parts(res, ctx)
            else:
                res_re, res_im = _extract_complex_parts(res, ctx)
        else:
            res_re, res_im = _extract_complex_parts(res, ctx)
        
        # 替换虚部中可能的未绑定变量
        for var_name, var_limit in limit_info.items():
            res_re = res_re.subst(var_name, var_limit) if res_re else res_re
            res_im = res_im.subst(var_name, var_limit) if res_im else res_im
        
        # normalize
        res_re = normalize(res_re, ctx) if res_re else Const(0)
        res_im = normalize(res_im, ctx) if res_im else Const(0)
        
        # 获取路径方向
        direction = get_contour_direction(path, ctx, start_expr, end_expr, limit_info)
        
        # 比较res的虚部和极点的虚部
        res_im_val = normalize(res_im, ctx)
        point_im_val = normalize(point_im, ctx)
        
        # 替换虚部中的未绑定变量为极限值
        for var_name, var_limit in limit_info.items():
            res_im_val = res_im_val.subst(var_name, var_limit)
            point_im_val = point_im_val.subst(var_name, var_limit)
        
        # 再次normalize
        res_im_val = normalize(res_im_val, ctx)
        point_im_val = normalize(point_im_val, ctx)
        
        # 计算虚部差值
        im_diff = Op("-", res_im_val, point_im_val)
        im_diff = normalize(im_diff, ctx)
        
        # 判断虚部大小关系并计算跳变值
        jump_val = compute_jump_value(im_diff, direction, ctx)
        jump_values.append(jump_val)
    
    return jump_values

def compute_cauchy_index(path: Union[CINTPath], point: Expr, ctx: Context) -> float:
    """计算路径绕点的绕数（使用跳变法/Cauchy指数法）
    
    注意：开放路径的绕数可能是非整数（如0.5），只有闭合路径的绕数才是整数
    
    理论基础（基于形式化证明）：
    1. 构造辅助函数：f(t) = Im(γ(t)-z) / Re(γ(t)-z)
    2. 找到跳变点：Re(γ(t)-z) = 0 且 Im(γ(t)-z) ≠ 0 的点
    3. 对每个跳变点x，计算：
       - jump₊(f,x)：右极限 lim(u→x⁺) f(u) = ±∞ → ±1/2，否则为0
       - jump₋(f,x)：左极限 lim(u→x⁻) f(u) = ±∞ → ±1/2，否则为0
    4. 柯西指数：Indp(γ,z) = Σjump₊(f,x) - Σjump₋(f,x)
    5. 绕数公式：n(γ,z) = -Indp(γ,z)/2
    
    Args:
        path: CINTPath参数化路径对象 γ(t)
        point: 点 z
        ctx: 上下文
        
    Returns:
        绕数（整数）
    """
    # 获取跳变值列表
    jump_values = compute_jump_values(path, point, ctx)
    
    # 求和所有跳变值
    total_jump = sum(jump_values)
    
    # 除以2得到绕数
    cauchy_index = total_jump / 2.0
    
    # 返回绕数（可能是非整数，对于开放路径）
    return cauchy_index



def find_jump_points(diff_re: Expr, var: str, t_start: Expr, t_end: Expr, ctx: Context) -> list[float]:
    """寻找跳跃点：符号求解分母为零的点
    
    跳跃点定义：diff_re = 0 且 diff_im ≠ 0 的孤立点
    
    验证方法：
    1. 求解 diff_re = 0 得到候选点 sol
    2. 将 sol 带入路径参数方程 γ(t)
    3. 使用 expand_euler 展开欧拉公式
    4. 使用 normalize 化简
    5. 验证 Re(γ(sol) - z₀) = 0
    
    Args:
        diff_re: 分母表达式 Re(γ(t)-z₀)
        diff_im: 分子表达式 Im(γ(t)-z₀)
        path_expr: 原始路径参数方程 γ(t)
        var: 参数变量名
        t_start: 参数起始值
        t_end: 参数结束值
        ctx: 上下文
        
    Returns:
        跳跃点列表
    """
    from integral.solve import solve_equation
    from integral.poly import normalize
    
    jump_points = []
    
    # 验证函数：将 sol 带入 diff_re 并验证是否为0
    def verify_jump_point(sol_val: Expr) -> bool:
        """验证 sol 处 Re(γ(t) - z0) 是否为0
        
        策略：优先符号判断，失败则数值验证兜底
        """
        try:
            # 1. 将 sol 带入 diff_re
            diff_re_at_sol = diff_re.subst(var, sol_val)
            
            # 2. 展开欧拉公式（如果有）
            diff_re_expanded = expand_euler(diff_re_at_sol, ctx)
            
            # 3. normalize 化简
            diff_re_normalized = normalize(diff_re_expanded, ctx)
            
            # 4. 符号判断是否为0
            if diff_re_normalized == Const(0):
                return True
            
            # 5. 符号判断失败，尝试数值验证兜底
            import math
            from decimal import Decimal
            diff_with_pi = diff_re_normalized.subst('pi', Const(Decimal(str(math.pi))))
            val = expr.eval_expr(diff_with_pi)
            return abs(val) < 1e-6
        except:
            return False
    
    # 1. 符号求解 diff_re = 0
    symbolic_solved = False
    try:
        # solve_equation(f, a, x, ctx) 求解 f = a
        solutions = solve_equation(diff_re, Const(0), var, ctx)
        if solutions:
            # 筛选在区间内的解
            for sol in solutions:
                try:
                    # 检查是否在区间内（包含边界）
                    # 支持双向区间：[a,b]或[b,a]
                    sol = normalize(sol, ctx)
                    in_interval = (t_start <= sol <= t_end) or (t_end <= sol <= t_start)
                    if in_interval:
                        # 使用新的验证方法：将 sol 带入路径参数方程验证实部是否为0
                        if verify_jump_point(sol):
                            jump_points.append(sol)
                            symbolic_solved = True
                except:
                    continue
    except:
        pass
    # 2. 如果符号求解无解或解不完整，测试特殊点
    if not symbolic_solved or len(jump_points) < 2:
        # 对于三角函数，测试可能的零点：0.5, 1.5, 2.5等（对应cos(pi*t)=0的点）
        test_points = []
        t_start_val = float(expr.eval_expr(t_start))
        t_end_val = float(expr.eval_expr(t_end))
        
        # 确定区间的最小值和最大值
        t_min = min(t_start_val, t_end_val)
        t_max = max(t_start_val, t_end_val)
        
        # 生成测试点：在区间内，步长0.5和1.0
        t = t_min + 0.5
        while t <= t_max + 0.01:  # 稍微超出以包含边界
            test_points.append(t)
            t += 0.5
        
        # 同时测试整数点附近
        for i in range(int(t_min) - 1, int(t_max) + 2):
            for offset in [0, 0.25, 0.5, 0.75, 1.0]:
                t_test = i + offset
                if t_min - 0.01 <= t_test <= t_max + 0.01:
                    test_points.append(t_test)
        
        # 去重并排序
        test_points = sorted(set(test_points))
        
        # 对每个测试点进行验证
        for t_test in test_points:
            # 使用新的验证方法（传入 Const 类型）
            from fractions import Fraction
            t_test_expr = Const(Fraction(t_test).limit_denominator(10000))
            if verify_jump_point(t_test_expr):
                # 避免重复添加（将 Expr 类型的 jump_points 转为 float 比较）
                existing_floats = [float(expr.eval_expr(jp)) if isinstance(jp, Expr) else jp for jp in jump_points]
                if not any(abs(t_test - existing) < 0.01 for existing in existing_floats):
                    jump_points.append(t_test_expr)
                    symbolic_solved = True
    
    # 去重并排序（需要先转为数值才能排序）
    jump_points_floats = [(float(expr.eval_expr(jp)) if isinstance(jp, Expr) else jp, jp) for jp in jump_points]
    jump_points_floats = sorted(set(jump_points_floats), key=lambda x: x[0])
    jump_points = [jp for _, jp in jump_points_floats]
    return jump_points

# 用欧拉公式展开路径表达式：exp(i*x) → cos(x) + i*sin(x)
def expand_euler(e: Expr, ctx: Context) -> Expr:
    """递归展开 exp(±i*θ) 为 cos(θ) ± i*sin(θ)
    
    增强版本，支持任意形式的虚数指数：
    - exp(i) → cos(1) + i*sin(1)
    - exp(i*t) → cos(t) + i*sin(t)
    - exp(-i*t) → cos(t) - i*sin(t)
    - exp(i*π/2) → cos(π/2) + i*sin(π/2)
    - exp(i*π*(1-t)) → cos(π*(1-t)) + i*sin(π*(1-t))
    - exp(k*i*t) → cos(k*t) + i*sin(k*t)
    - exp((i*a)/b) → cos(a/b) + i*sin(a/b)
    - exp(a + i*b) → exp(a) * [cos(b) + i*sin(b)]
    """
    if isinstance(e, Fun) and e.func_name == 'exp' and len(e.args) == 1:
        arg = e.args[0]
        
        # 检查是否包含虚数单位 i（递归检查所有子表达式）
        def contains_i(expr: Expr) -> bool:
            if isinstance(expr, Fun) and expr.func_name == 'i':
                return True
            if isinstance(expr, Op):
                return any(contains_i(a) for a in expr.args)
            return False
        
        if not contains_i(arg):
            return e  # 不包含 i，是实数指数，保持原样
        
        # 先对arg进行normalize
        arg_normalized = normalize(arg, ctx)
        
        # 辅助函数：扁平化乘法表达式
        def flatten_mult(expr: Expr) -> list[Expr]:
            """将嵌套的乘法表达式扁平化为因子列表"""
            if isinstance(expr, Op) and expr.op == '*':
                result = []
                for arg in expr.args:
                    result.extend(flatten_mult(arg))
                return result
            else:
                return [expr]
        
        # 提取实部和虚部的辅助函数（简化版，专门用于exp参数）
        def extract_parts(expr: Expr) -> tuple[Expr, Expr]:
            """提取 a + i*b 形式中的 a 和 b"""
            # 处理纯虚数：i, -i, i*theta, -i*theta, (i*a)/b 等
            if isinstance(expr, Fun) and expr.func_name == 'i':
                return (Const(0), Const(1))
            
            # 一元负号
            if isinstance(expr, Op) and expr.op == '-' and len(expr.args) == 1:
                re, im = extract_parts(expr.args[0])
                return (Op('-', re) if re != Const(0) else Const(0),
                        Op('-', im) if im != Const(0) else Const(0))
            
            # 加法：a + i*b
            if isinstance(expr, Op) and expr.op == '+':
                # 递归分离各项
                real_sum = Const(0)
                imag_sum = Const(0)
                for term in expr.args:
                    re, im = extract_parts(term)
                    if re != Const(0):
                        real_sum = re if real_sum == Const(0) else Op('+', real_sum, re)
                    if im != Const(0):
                        imag_sum = im if imag_sum == Const(0) else Op('+', imag_sum, im)
                return (real_sum, imag_sum)
            
            # 减法：a - i*b
            if isinstance(expr, Op) and expr.op == '-' and len(expr.args) == 2:
                re_left, im_left = extract_parts(expr.args[0])
                re_right, im_right = extract_parts(expr.args[1])
                real_part = Op('-', re_left, re_right) if re_right != Const(0) else re_left
                imag_part = Op('-', im_left, im_right) if im_right != Const(0) else im_left
                return (real_part, imag_part)
            
            # 乘法或除法：检查是否是 i*theta 或 theta*i 或 (i*theta)/b 等形式
            if contains_i(expr):
                # 尝试提取 theta（从 i*theta 形式）
                # 使用更通用的方法：将表达式视为 i * (expr/i)
                # 这里我们直接检查表达式结构
                
                # 对于乘法：需要扁平化以处理嵌套的乘法
                if isinstance(expr, Op) and expr.op == '*':
                    # 扁平化乘法表达式
                    factors = flatten_mult(expr)
                    
                    i_found = False
                    other_terms = []
                    for factor in factors:
                        if isinstance(factor, Fun) and factor.func_name == 'i':
                            i_found = True
                        else:
                            other_terms.append(factor)
                    
                    if i_found:
                        if len(other_terms) == 0:
                            return (Const(0), Const(1))
                        elif len(other_terms) == 1:
                            return (Const(0), other_terms[0])
                        else:
                            theta = other_terms[0]
                            for t in other_terms[1:]:
                                theta = Op('*', theta, t)
                            return (Const(0), theta)
                
                # 对于除法：(i*a)/b 或 i/b
                if isinstance(expr, Op) and expr.op == '/':
                    numerator, denominator = expr.args
                    re_num, im_num = extract_parts(numerator)
                    if im_num != Const(0):
                        # 虚部除以分母
                        return (Const(0), Op('/', im_num, denominator))
            
            # 如果不包含i，则全部是实部
            return (expr, Const(0))
        
        real_part, imag_part = extract_parts(arg_normalized)
        
        # 根据实部和虚部生成展开式
        if real_part != Const(0) and imag_part != Const(0):
            # exp(a + i*b) = exp(a) * [cos(b) + i*sin(b)]
            exp_real = Fun('exp', real_part)
            euler_part = Op('+', Fun('cos', imag_part), 
                          Op('*', Fun('i'), Fun('sin', imag_part)))
            return Op('*', exp_real, euler_part)
        elif imag_part != Const(0):
            # 纯虚数指数：exp(i*θ) = cos(θ) + i*sin(θ)
            return Op('+', Fun('cos', imag_part), 
                        Op('*', Fun('i'), Fun('sin', imag_part)))
        # 如果只有实部（理论上不应该到这里，因为前面检查了contains_i）
        # 保持原样
    
    # 递归处理子表达式
    if isinstance(e, Op):
        return Op(e.op, *[expand_euler(a, ctx) for a in e.args])
    elif isinstance(e, Fun):
        return Fun(e.func_name, *[expand_euler(a, ctx) for a in e.args])
    return e

def _extract_complex_parts(z: Expr, ctx: Context) -> tuple[Expr, Expr]:
    """提取复数的实部和虚部
    
    处理路径参数方程的两种形式：
    1. 圆弧：a +/- bi + r*exp(i*pi*t) 或 a +/- bi + r*exp(i*pi*(1-t))
    2. 直线：a + b*i + r*(1-2*t)*i 或 a + b*i + r*(1-2*t)
    
    流程：
    1. 如果存在exp(i*...)，应用欧拉公式展开为cos+i*sin，并对三角函数参数应用ExpandPolynomial
    2. 如果没有exp（直线），直接对整个表达式应用ExpandPolynomial
    3. 对整个表达式应用normalize化简
    4. 提取实部和虚部
    """
    
    # 辅助函数：检查表达式是否包含i
    def contains_i(e: Expr) -> bool:
        """检查表达式是否包含虚数单位i"""
        if isinstance(e, Fun) and e.func_name == 'i':
            return True
        elif isinstance(e, Op):
            return any(contains_i(arg) for arg in e.args)
        return False
    
    # 辅助函数：检查是否包含exp函数
    def has_exp(e: Expr) -> bool:
        """检查表达式是否包含exp函数"""
        if isinstance(e, Fun) and e.func_name == 'exp':
            return True
        elif isinstance(e, Op):
            return any(has_exp(arg) for arg in e.args)
        return False
    
    # 辅助函数：查找并展开exp函数
    def expand_exp(e: Expr) -> Expr:
        """查找exp函数，如果参数包含i，应用欧拉公式展开"""
        if isinstance(e, Fun) and e.func_name == 'exp':
            # 检查exp的参数是否包含i
            arg = e.args[0]
            if contains_i(arg):
                # 应用欧拉公式: exp(i*theta) = cos(theta) + i*sin(theta)
                # 提取theta（i乘的部分）
                theta = extract_theta_from_i_mult(arg)
                if theta is not None:
                    # 对theta应用多项式展开
                    expander = ExpandPolynomial()
                    theta_expanded = expander.eval(theta, ctx)
                    
                    # 构建 cos(theta) + i*sin(theta)
                    cos_part = Fun('cos', theta_expanded)
                    sin_part = Fun('sin', theta_expanded)
                    return Op('+', cos_part, Op('*', Fun('i'), sin_part))
        
        # 递归处理Op类型
        if isinstance(e, Op):
            new_args = [expand_exp(arg) for arg in e.args]
            return Op(e.op, *new_args)
        
        return e
    
    # 辅助函数：扁平化乘法表达式
    def flatten_mult(e: Expr) -> list:
        """将嵌套的乘法表达式扁平化为因子列表"""
        if isinstance(e, Op) and e.op == '*':
            result = []
            for arg in e.args:
                result.extend(flatten_mult(arg))
            return result
        else:
            return [e]
    
    # 辅助函数：从i*theta中提取theta
    def extract_theta_from_i_mult(e: Expr) -> Expr:
        """从 i*theta 或 theta*i 或 -(i*theta) 等形式中提取theta
        
        增强版本，支持：
        - i -> theta=1
        - i*a -> theta=a
        - a*i -> theta=a
        - -i -> theta=-1
        - -(i*a) -> theta=-a
        - i*pi/2 -> theta=pi/2
        - (i*a)/b -> theta=a/b
        - i*a*b/c -> theta=a*b/c
        等任意复杂形式
        """
        # 如果直接是i，theta=1
        if isinstance(e, Fun) and e.func_name == 'i':
            return Const(1)
        
        # 处理一元负号: -(i*theta) -> theta' = -theta
        if isinstance(e, Op) and e.op == '-' and len(e.args) == 1:
            inner_theta = extract_theta_from_i_mult(e.args[0])
            if inner_theta is not None:
                # 返回 -theta
                return Op('-', inner_theta)
        
        # 处理除法: (i*a)/b 或 复杂形式
        if isinstance(e, Op) and e.op == '/':
            numerator, denominator = e.args
            # 递归提取分子中的theta
            theta_num = extract_theta_from_i_mult(numerator)
            if theta_num is not None:
                # 构建 theta_num / denominator
                return Op('/', theta_num, denominator)
            # 如果分子不含i，检查是否是 a/(i*b) 形式（这种情况下theta=-i*a/b，但这不是标准形式）
            return None
        
        # 处理乘法
        if isinstance(e, Op) and e.op == '*':
            # 扁平化乘法表达式
            factors = flatten_mult(e)
            
            # 查找i因子
            i_factors = [f for f in factors if isinstance(f, Fun) and f.func_name == 'i']
            other_factors = [f for f in factors if not (isinstance(f, Fun) and f.func_name == 'i')]
            
            if len(i_factors) == 1:
                # 找到i，其他因子组成theta
                if len(other_factors) == 0:
                    return Const(1)  # 只有i，theta=1
                elif len(other_factors) == 1:
                    return other_factors[0]
                else:
                    # 将多个因子组合成嵌套的 Op
                    result = other_factors[0]
                    for factor in other_factors[1:]:
                        result = Op('*', result, factor)
                    return result
            elif len(i_factors) > 1:
                # 多个i因子，这通常意味着 i*i = -1，但这里我们只处理单个i的情况
                return None
        
        # 处理加法/减法中可能包含i的情况（如 a+i*b 形式，这不是纯虚数指数）
        # 这种情况返回None，让调用者处理
        return None
    
    # 第一步：检查是否包含exp并展开
    has_exp_func = has_exp(z)
    
    if has_exp_func:
        # 圆弧路径：展开exp（内部已对三角函数参数应用多项式展开）
        z_expanded = expand_exp(z)
    else:
        # 直线路径：直接对整个表达式应用多项式展开
        expander = ExpandPolynomial()
        z_expanded = expander.eval(z, ctx)
    
    # 第二步：对整个表达式应用normalize化简
    z_normalized = normalize(z_expanded, ctx)
    
    # 第三步：从化简后的表达式提取实部和虚部
    def extract_re_im(expr: Expr) -> tuple[Expr, Expr]:
        """从展开后的表达式提取实部和虚部"""
        
        # 1. 纯虚数 i
        if isinstance(expr, Fun) and expr.func_name == 'i':
            return (Const(0), Const(1))
        
        # 2. 一元负号
        if isinstance(expr, Op) and expr.op == '-' and len(expr.args) == 1:
            inner = expr.args[0]
            if isinstance(inner, Fun) and inner.func_name == 'i':
                return (Const(0), Const(-1))
            re, im = extract_re_im(inner)
            return (Op('-', re) if re != Const(0) else Const(0), 
                    Op('-', im) if im != Const(0) else Const(0))
        
        # 3. 乘法：处理 c*(cos+i*sin) 或 c*i 形式
        if isinstance(expr, Op) and expr.op == '*':
            # 扁平化乘法
            factors = flatten_mult(expr)
            
            # 分离实数和复数因子
            i_factors = [f for f in factors if isinstance(f, Fun) and f.func_name == 'i']
            complex_factors = [f for f in factors if contains_i(f) and not (isinstance(f, Fun) and f.func_name == 'i')]
            real_factors = [f for f in factors if not contains_i(f)]
            
            if i_factors and not complex_factors:
                # 简单的 c*i 形式
                if len(real_factors) == 0:
                    return (Const(0), Const(1))
                elif len(real_factors) == 1:
                    return (Const(0), real_factors[0])
                else:
                    return (Const(0), Op('*', *real_factors))
            
            if complex_factors:
                # 合并复数因子
                if len(complex_factors) == 1:
                    complex_part = complex_factors[0]
                else:
                    complex_part = Op('*', *complex_factors)
                
                # 递归提取复数部分的实虚部
                re_c, im_c = extract_re_im(complex_part)
                
                # 乘以实数因子
                if real_factors:
                    if len(real_factors) == 1:
                        real_coeff = real_factors[0]
                    else:
                        real_coeff = Op('*', *real_factors)
                    
                    re_part = Op('*', real_coeff, re_c) if re_c != Const(0) else Const(0)
                    im_part = Op('*', real_coeff, im_c) if im_c != Const(0) else Const(0)
                    return (re_part, im_part)
                else:
                    return (re_c, im_c)
        
        # 4. 除法：处理 i/a, (a+i*b)/c 等形式
        if isinstance(expr, Op) and expr.op == '/':
            numerator, denominator = expr.args
            re_num, im_num = extract_re_im(numerator)
            
            # 如果分母也包含虚数，需要复数除法 (a+bi)/(c+di) = [(ac+bd) + (bc-ad)i]/(c²+d²)
            if contains_i(denominator):
                re_den, im_den = extract_re_im(denominator)
                # (a+bi)/(c+di) 的实部 = (ac+bd)/(c²+d²)
                # (a+bi)/(c+di) 的虚部 = (bc-ad)/(c²+d²)
                den_sq = Op('+', Op('^', re_den, Const(2)), Op('^', im_den, Const(2)))
                re_result = Op('/', Op('+', Op('*', re_num, re_den), Op('*', im_num, im_den)), den_sq)
                im_result = Op('/', Op('-', Op('*', im_num, re_den), Op('*', re_num, im_den)), den_sq)
                return (re_result, im_result)
            else:
                # 分母是实数，简单除法
                re_result = Op('/', re_num, denominator) if re_num != Const(0) else Const(0)
                im_result = Op('/', im_num, denominator) if im_num != Const(0) else Const(0)
                return (re_result, im_result)
        
        # 5. 加法或减法
        if isinstance(expr, Op) and expr.op in ['+', '-'] and len(expr.args) == 2:
            left, right = expr.args
            
            re_left, im_left = extract_re_im(left)
            re_right, im_right = extract_re_im(right)
            
            if expr.op == '+':
                re_total = Op('+', re_left, re_right) if re_left != Const(0) or re_right != Const(0) else Const(0)
                im_total = Op('+', im_left, im_right) if im_left != Const(0) or im_right != Const(0) else Const(0)
            else:  # '-'
                re_total = Op('-', re_left, re_right) if re_left != Const(0) or re_right != Const(0) else Const(0)
                im_total = Op('-', im_left, im_right) if im_left != Const(0) or im_right != Const(0) else Const(0)
            
            return (re_total, im_total)
        
        # 6. 纯实数（默认情况）
        return (expr, Const(0))
    
    return extract_re_im(z_normalized)


def get_contour_direction(path: CINTPath, ctx: Context, start: Expr, end: Expr, limit_info: dict = None) -> str:
    """获取围道方向"""
    if limit_info is None:
        limit_info = {}
        
    # 检查参数变化方向
    curr = normalize(expand_euler(path.path_expr, ctx), ctx)
    path_start = normalize(ExpandPolynomial().eval(curr.subst(path.var, start), ctx), ctx)
    path_end = normalize(ExpandPolynomial().eval(curr.subst(path.var, end), ctx), ctx)
    
    # 替换未绑定变量为其极限值（使用limit_info）
    for var_name, var_limit in limit_info.items():
        path_start = path_start.subst(var_name, var_limit)
        path_end = path_end.subst(var_name, var_limit)

    start_re, _ = _extract_complex_parts(path_start, ctx)
    end_re, _ = _extract_complex_parts(path_end, ctx)

    # 尝试数值比较
    try:
        import math
        from decimal import Decimal
        
        # 替换pi和处理POS_INF/NEG_INF
        start_re_expr = start_re.subst('pi', Const(Decimal(str(math.pi))))
        end_re_expr = end_re.subst('pi', Const(Decimal(str(math.pi))))
        
        # 处理POS_INF和NEG_INF
        if start_re == expr.POS_INF:
            start_re_val = float('inf')
        elif start_re == expr.NEG_INF:
            start_re_val = float('-inf')
        else:
            start_re_val = expr.eval_expr(start_re_expr)
            
        if end_re == expr.POS_INF:
            end_re_val = float('inf')
        elif end_re == expr.NEG_INF:
            end_re_val = float('-inf')
        else:
            end_re_val = expr.eval_expr(end_re_expr)
        
        if isinstance(start_re_val, complex):
            start_re_val = start_re_val.real
        if isinstance(end_re_val, complex):
            end_re_val = end_re_val.real
            
        if start_re_val > end_re_val:
            return "R->L"
        if start_re_val < end_re_val:
            return "L->R"
    except Exception as e:
        pass
    
    return "L->R"

def compute_jump_value(im_diff: Expr, direction: str, ctx: Context) -> float:
    """根据虚部差值和路径方向计算跳变值
    
    Args:
        im_diff: res的虚部 - 极点的虚部
        direction: 路径方向 "L->R" 或 "R->L"
        ctx: 上下文
        
    Returns:
        跳变值：1 或 -1
    """
    # 判断im_diff的符号（im_diff > 0 表示res虚部 > 极点虚部）
    try:
        import math
        from decimal import Decimal
        
        # 替换pi并求值
        im_diff_val = im_diff.subst('pi', Const(Decimal(str(math.pi))))
        val = expr.eval_expr(im_diff_val)
        
        # 处理复数情况（取实部）
        if isinstance(val, complex):
            val = val.real
        
        # 根据方向和虚部比较结果确定跳变值
        if direction == "R->L":
            # R->L方向：res虚部>极点虚部 -> 1，否则 -> -1
            return 1.0 if val > 0 else -1.0
        else:  # L->R
            # L->R方向：res虚部>极点虚部 -> -1，否则 -> 1
            return -1.0 if val > 0 else 1.0
    except:
        # 如果无法判断，返回0
        return 0.0