"""Rules for integration."""
import re
from decimal import Decimal
from fractions import Fraction
from typing import Optional, Dict, Tuple, Union, List, Set
import functools
import operator

from sympy import false

from integral import expr, context
from integral.expr import Var, Const, Fun, EvalAt, Op, Integral, Symbol, Expr, \
    OP, CONST, VAR, sin, cos, FUN, decompose_expr_factor, \
    Deriv, Inf, Limit, NEG_INF, POS_INF, IndefiniteIntegral, Summation, SUMMATION, INTEGRAL, INF, \
    Product, SYMBOL, SkolemFunc, decompose_expr_factor2, is_const, exprify, Complex
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


class RuleException(Exception):
    def __init__(self, rule_name: str, msg: str):
        self.rule_name = rule_name
        self.msg = msg

    def __str__(self):
        return "%s: %s" % (self.rule_name, self.msg)


def deriv(var: str, e: Expr, ctx: Context) -> Expr:
    """Compute the derivative of e with respect to variable
    name var.

    """

    def normal(x):
        return normalize(x, ctx)

    def rec(e):
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
                elif not x.contains_var(var) and y.ty == OP and y.op == "^":
                    # c / (y0 ^ y1): rewrite to c * y0 ^ (-y1)
                    return rec(x * (y.args[0] ^ (-y.args[1])))
                else:
                    # general case
                    return normal((rec(x) * y - x * rec(y)) / (y ^ Const(2)))
            elif e.op == "^":
                x, y = e.args
                if y.ty == CONST:
                    return normal(y * (x ^ Const(y.val - 1)) * rec(x))
                elif var not in y.get_vars():
                    return normal(y * (x ^ (y - 1)) * rec(x))
                else:
                    return normal(rec(expr.exp(y * expr.log(x))))

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
        elif isinstance(e, Complex):
            # For a complex expression z = x + i * y
            # dz/dt = dx/dt + i * dy/dt
            real_deriv = rec(e.real)
            imag_deriv = rec(e.imag)
            return Complex(real_deriv, imag_deriv)
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
    def __init__(self, exprs: List[Expr], flags: List[bool] = None):
        self.exprs = exprs  # satisfy all expressions
        if flags is None or len(flags) != len(exprs):
            self.need_to_be_satisfied = [True for i in range(len(exprs))]
        else:
            self.need_to_be_satisfied = flags

    def __str__(self):
        res = ", ".join([str(e) for e in self.exprs])
        res += "\n"
        res += ", ".join([str(b) for b in self.need_to_be_satisfied])
        return res

    def export(self):
        res = {
            'exprs': [str(e) for e in self.exprs]
        }
        return res


class ProofObligation:
    """Represents a proof obligation to prove e using the conditions
    in conds.

    """

    def __init__(self, branches: List['ProofObligationBranch'], conds: Conditions):
        self.branches = branches  # if any branch is satisfied then the proof obligation is carried out
        self.conds = conds

    def __eq__(self, other: "ProofObligation"):
        return self.branches == other.branches

    def __le__(self, other: "ProofObligation"):
        if len(self.branches) < other.branches:
            return True
        elif len(self.branches > other.branches):
            return False
        return all(a < b for a, b in zip(self.branches, other.branches))

    def __str__(self):
        res = ""
        for i, b in enumerate(self.branches):
            res += "branch " + str(i) + ":\n"
            res += str(b) + '\n'
        res += str(self.conds) + "\n"
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


def check_wellformed(e: Expr, ctx: Context) -> List[ProofObligation]:
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
                if ctx.check_condition(Op(">", e.args[0], Const(0))):
                    pass
                else:
                    add_obligation(Op(">", e.args[0], Const(0)), ctx)
            if e.func_name == 'sqrt':
                if ctx.check_condition(Op(">=", e.args[0], Const(0))):
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
                f1 = ctx.check_condition(Fun("isInt", tmp))
                f2 = ctx.check_condition(Fun("isEven", tmp))

                if not f1 or f2:
                    pass
                else:
                    branch1 = ProofObligationBranch([Fun("isInt", tmp)], [False])
                    branch2 = ProofObligationBranch([Fun("isEven", tmp)])
                    add_obligation([branch1, branch2], ctx)
            # TODO: add checks for other functions
        elif expr.is_integral(e):
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
        def rec(e: Expr):
            if expr.is_integral(e):
                if isinstance(e.body, Complex):
                    # 分离实部和虚部进行线性积分
                    real_part = expr.Integral(e.var, e.lower, e.upper, e.body.real)
                    imag_part = expr.Integral(e.var, e.lower, e.upper, e.body.imag)
                    return Complex(rec(real_part), rec(imag_part))
                elif e.body.is_plus():
                    return rec(expr.Integral(e.var, e.lower, e.upper, e.body.args[0])) + \
                           rec(expr.Integral(e.var, e.lower, e.upper, e.body.args[1]))
                elif expr.is_uminus(e.body):
                    return -rec(expr.Integral(e.var, e.lower, e.upper, e.body.args[0]))
                elif e.body.is_minus():
                    return rec(expr.Integral(e.var, e.lower, e.upper, e.body.args[0])) - \
                           rec(expr.Integral(e.var, e.lower, e.upper, e.body.args[1]))
                elif e.body.is_times() or e.body.is_divides():
                    num_factors, denom_factors = decompose_expr_factor(e.body)
                    b = prod(f for f in num_factors if f.contains_var(e.var))
                    c = prod(f for f in num_factors if not f.contains_var(e.var))
                    denom_b = prod(f for f in denom_factors if f.contains_var(e.var))
                    denom_c = prod(f for f in denom_factors if not f.contains_var(e.var))
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
                if isinstance(e.body, Complex):
                    # 分离实部和虚部进行线性积分
                    real_part = expr.IndefiniteIntegral(e.var, e.body.real, e.skolem_args)
                    imag_part = expr.IndefiniteIntegral(e.var, e.body.imag, e.skolem_args)
                    return Complex(rec(real_part), rec(imag_part))
                if e.body.is_plus():
                    return rec(expr.IndefiniteIntegral(e.var, e.body.args[0], e.skolem_args)) + \
                        rec(expr.IndefiniteIntegral(e.var, e.body.args[1], e.skolem_args))
                elif expr.is_uminus(e.body):
                    return -rec(IndefiniteIntegral(e.var, e.body.args[0], e.skolem_args))
                elif e.body.is_minus():
                    return rec(expr.IndefiniteIntegral(e.var, e.body.args[0], e.skolem_args)) - \
                        rec(expr.IndefiniteIntegral(e.var, e.body.args[1], e.skolem_args))
                elif e.body.is_times() or e.body.is_divides():
                    num_factors, denom_factors = decompose_expr_factor(e.body)
                    b = prod(f for f in num_factors if f.contains_var(e.var))
                    c = prod(f for f in num_factors if not f.contains_var(e.var))
                    denom_b = prod(f for f in denom_factors if f.contains_var(e.var))
                    denom_c = prod(f for f in denom_factors if not f.contains_var(e.var))
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
                elif e.body.is_times() or e.body.is_divides():
                    num_factors, denom_factors = decompose_expr_factor(e.body)
                    b, c = Const(1), Const(1)
                    for f in num_factors:
                        if not f.contains_var(e.var):
                            c = c * f
                        else:
                            b = b * f
                    for f in denom_factors:
                        if not f.contains_var(e.var):
                            c = c / f
                        else:
                            b = b / f
                    return c * Limit(e.var, e.lim, b)
                else:
                    return e
            elif expr.is_summation(e):
                v, l, u, body = e.index_var, e.lower, e.upper, e.body
                if e.body.is_minus():
                    return Summation(v, l, u, body.args[0]) - Summation(v, l, u, body.args[1])
                elif expr.is_uminus(e.body):
                    return -Summation(v, l, u, body.args[0])
                elif e.body.is_times() or e.body.is_divides():
                    num_factors, denom_factors = decompose_expr_factor(e.body)
                    b, c = Const(1), Const(1)
                    for f in num_factors:
                        if not f.contains_var(e.index_var):
                            c = c * f
                        else:
                            b = b * f
                    for f in denom_factors:
                        if not f.contains_var(e.index_var):
                            c = c / f
                        else:
                            b = b / f
                    c = normalize(c, ctx)
                    b = normalize(b, ctx)
                    return normalize(c * Summation(e.index_var, e.lower, e.upper, b), ctx)
                else:
                    return e
            elif expr.is_complex(e):
                real_part = e.real
                imag_part = e.imag
                # 对实部进行化简
                if isinstance(real_part, Op) and real_part.op == "+":
                    # 如果实部是加法，递归化简
                    simplified_real = sum([rec(arg) for arg in real_part.args])
                elif isinstance(real_part, Op) and real_part.op == "-":
                    # 如果实部是减法，递归化简
                    simplified_real = rec(real_part.args[0]) - rec(real_part.args[1])
                else:
                    # 如果实部不是加法或减法，直接使用原实部
                    simplified_real = rec(real_part)
                # 对虚部进行化简
                if isinstance(imag_part, Op) and imag_part.op == "+":
                    simplified_imag = sum([rec(arg) for arg in imag_part.args])
                elif isinstance(imag_part, Op) and imag_part.op == "-":
                    simplified_imag = rec(imag_part.args[0]) - rec(imag_part.args[1])
                else:
                    simplified_imag = rec(imag_part)

                return Complex(simplified_real, simplified_imag)
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


class DefiniteIntegralIdentity(Rule):
    """Apply definite integral identity in current theory."""

    def __init__(self):
        self.name = "DefiniteIntegralIdentity"

    def __str__(self):
        return "apply integral identity"

    def export(self):
        return {
            "name": self.name,
            "str": str(self)
        }

    def eval(self, e: Expr, ctx: Context) -> Expr:
        # Apply linearity
        if expr.is_integral(e) or expr.is_indefinite_integral(e):
            e = Linearity().eval(e, ctx)

        if not (expr.is_integral(e) or expr.is_indefinite_integral(e)):
            sep_ints = e.separate_integral()
            for _, loc in sep_ints:
                e = OnLocation(self, loc).eval(e, ctx)
            return e

        # First, look for indefinite integrals identities
        for identity in ctx.get_indefinite_integrals():
            inst = expr.match(IndefiniteIntegral(e.var, e.body, skolem_args=tuple()), identity.lhs)
            if inst is None:
                continue

            inst[identity.lhs.var] = Var(e.var)
            assert identity.rhs.is_plus() and expr.is_skolem_func(identity.rhs.args[1])
            pat_rhs = identity.rhs.args[0]  # remove Skolem constant C
            return EvalAt(e.var, e.lower, e.upper, normalize(pat_rhs.inst_pat(inst), ctx))

        # Look for definite integral identities
        for identity in ctx.get_definite_integrals():
            inst = expr.match(e, identity.lhs)
            if inst is not None:
                # Check conditions
                satisfied = True
                for cond in identity.conds.data:
                    cond = expr.expr_to_pattern(cond)
                    cond = cond.inst_pat(inst)
                    if not ctx.check_condition(cond):
                        satisfied = False
                if satisfied:
                    return normalize(identity.rhs.inst_pat(inst), ctx)

        # No matching identity found
        return e


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
            res = identity.rhs.inst_pat(inst)
            assert expr.is_summation(res)
            res = res.alpha_convert(self.index_var)
            return res
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

            return identity.rhs.inst_pat(inst)
        # No matching identity found
        return e


class IndefiniteIntegralIdentity(Rule):
    """Apply indefinite integral identity in current theory."""

    def __init__(self):
        self.name = "IndefiniteIntegralIdentity"

    def __str__(self):
        return "apply indefinite integral"

    def export(self):
        return {
            "name": self.name,
            "str": str(self)
        }

    def eval(self, e: Expr, ctx: Context) -> Expr:
        """Apply indefinite integral identity to expression."""
        def apply(e: Expr):
            for indef in ctx.get_indefinite_integrals():
                inst = expr.match(e, indef.lhs)
                if inst is None:
                    continue

                inst['x'] = Var(e.var)
                assert indef.rhs.is_plus() and expr.is_skolem_func(indef.rhs.args[1])
                return indef.rhs.args[0].inst_pat(inst)

            # No matching identity found
            return e

        # Apply linearity
        if expr.is_integral(e) or expr.is_indefinite_integral(e):
            e = Linearity().eval(e, ctx)

        # Apply to right side of equation
        if e.is_equals():
            return OnLocation(self, "1").eval(e, ctx)

        integrals = e.separate_integral()
        skolem_args = set()
        for sub_e, loc in integrals:
            if expr.is_integral(sub_e):
                raise RuleException("apply indefinite integral", "Attempting to apply indefinite integral methods to a definite integral expression will not work")
            new_e = apply(sub_e)
            if new_e != sub_e:
                e = e.replace_expr(loc, new_e)
                skolem_args = skolem_args.union(set(sub_e.skolem_args))

        if e.is_plus() and expr.is_skolem_func(e.args[1]):
            # If already has Skolem variable at right
            skolem_args = skolem_args.union(set(arg.name for arg in e.args[1].dependent_vars))
            e = e.args[0] + expr.SkolemFunc(e.args[1].name, tuple(Var(arg) for arg in skolem_args))
        else:
            # If no Skolem variable at right
            e = e + expr.SkolemFunc("C", tuple(Var(arg) for arg in skolem_args))

        return e


class EvaluateIndefiniteIntegral(Rule):
    def __init__(self):
        self.name = "EvaluateIndefiniteIntegral"

    def eval(self, e: Expr, ctx: Context) -> Expr:
        assert isinstance(e, IndefiniteIntegral)
        for indef in ctx.get_indefinite_integrals():
            assert isinstance(indef.lhs, IndefiniteIntegral)
            inst = expr.match(e, indef.lhs)
            if inst is None:
                continue

            inst[indef.lhs.var] = Var(e.var)

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

        # First, try indefinite integral identities
        for identity in ctx.get_indefinite_integrals():
            assert isinstance(identity.lhs, IndefiniteIntegral)
            inst = expr.match(IndefiniteIntegral(e.var, e.body, skolem_args=tuple()), identity.lhs)
            if inst is None:
                continue

            inst[identity.lhs.var] = Var(e.var)

            # The right side of the identity should be of the form "expr + C"
            # take the expr part of the expression.
            assert identity.rhs.is_plus() and expr.is_skolem_func(identity.rhs.args[1])
            pat_rhs = identity.rhs.args[0]
            return EvalAt(e.var, e.lower, e.upper, normalize(pat_rhs.inst_pat(inst), ctx))

        # Next, look for definite integral identities
        for identity in ctx.get_definite_integrals():
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
                return normalize(identity.rhs.inst_pat(inst), ctx)

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

    def eval(self, e: Expr, ctx: Context) -> Expr:
        """Apply indefinite integral identity to expression."""

        # If incoming expression is equality, apply to the right side of equation
        if e.is_equals():
            lhs = self.eval(e.lhs, ctx)
            rhs = self.eval(e.rhs, ctx)
            return Op("=", lhs, rhs)

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
                e = OnLocation(EvaluateDefiniteIntegral(), loc).eval(e, ctx)
            else:
                raise AssertionError

        if exist_indefinite_integral:
            if e.is_plus() and expr.is_skolem_func(e.args[1]):
                # If already has Skolem variable at right
                skolem_args = skolem_args.union(set(arg.name for arg in e.args[1].dependent_vars))
                e = e.args[0] + expr.SkolemFunc(e.args[1].name, tuple(Var(arg) for arg in skolem_args))
            else:
                # If no Skolem variable at right
                e = e + expr.SkolemFunc("C", tuple(Var(arg) for arg in skolem_args))
        return e

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
        def rec(cur_e, loc, ctx):
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
                return IndefiniteIntegral(cur_e.var, rec(cur_e.body, loc.rest, ctx), cur_e.skolem_args)
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
                else:
                    return cur_e

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
            raise RuleException("OnCount",f"{self.n} is out of range")
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
        assert found, "ApplyEquation: lemma %s not found" % self.eq

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
        # print("solve equation %s for %s"%(self.eq, e))
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
        if not (expr.is_integral(e) or expr.is_indefinite_integral(e) or expr.is_limit(e)):
            sep_ints = e.separate_integral()
            sep_lims = e.separate_limits()
            if len(sep_ints) == 0 and len(sep_lims) == 0:
                raise RuleException("Substitution", "integral or limit not found")
            elif len(sep_ints) != 0:
                return OnLocation(self, sep_ints[0][1]).eval(e, ctx)
            else:
                return OnLocation(self, sep_lims[0][1]).eval(e, ctx)

        # Variable to be substituted in the integral
        var_name = Var(self.var_name)

        # Expression used for substitution
        var_subst = self.var_subst

        if e.var not in var_subst.get_vars():
            raise RuleException("Substitution", "variable %s not found" % e.var)

        # Compute g(x)'
        dfx = deriv(e.var, var_subst, ctx)

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
        body_subst2 = normalize(body, ctx).replace(normalize(var_subst, ctx), var_name)
        body_subst3 = body.replace(var_subst2, var_name)
        body_subst4 = normalize(body, ctx).replace(normalize(var_subst2,ctx), var_name)
        body_subst5 = normalize(body.replace(var_subst, var_name), ctx)
        body_subst6 = normalize(body.replace(var_subst2, var_name), ctx)
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
            gu = solve_equation(var_subst, var_name, e.var, ctx)
            if gu is None:
                raise RuleException("Substitution", "unable to solve equation %s = %s for %s, body_subst = %s" % (
                    var_subst, var_name, e.var, body_subst
                ))

            gu = normalize(gu, ctx)
            c = e.body.replace(Var(e.var), gu)
            if not expr.is_limit(e):
                new_problem_body = c * deriv(str(var_name), gu, ctx)
            else:
                new_problem_body = c
            self.f = new_problem_body

        if expr.is_integral(e):
            if e.lower == expr.NEG_INF:
                lower = limits.reduce_neg_inf_limit(var_subst, e.var, ctx)
            else:
                x = Var(e.var)
                lower = self.var_subst
                lower = limits.reduce_inf_limit(lower.subst(e.var, (1 / x) + e.lower), e.var, ctx)
                lower = normalize(lower, ctx)
            if e.upper == expr.POS_INF:
                upper = limits.reduce_inf_limit(var_subst, e.var, ctx)
            else:
                x = Var(e.var)
                upper = self.var_subst
                upper = limits.reduce_inf_limit(upper.subst(e.var, e.upper - (1 / x)), e.var, ctx)
                upper = normalize(upper, ctx)
            if lower.is_evaluable() and upper.is_evaluable() and expr.eval_expr(lower) > expr.eval_expr(upper):
                return normalize(Integral(self.var_name, upper, lower, Op("-", self.f)), ctx)
            else:
                return normalize(Integral(self.var_name, lower, upper, self.f), ctx)
        elif expr.is_indefinite_integral(e):
            return normalize(IndefiniteIntegral(self.var_name, self.f, e.skolem_args), ctx)
        elif expr.is_limit(e):
            # Perhaps need to be improved when drt is not None
            if e.lim == expr.NEG_INF:
                lim = limits.reduce_neg_inf_limit(var_subst, e.var, ctx)
            elif e.lim == expr.POS_INF:
                lim = limits.reduce_inf_limit(var_subst, e.var, ctx)
            else:
                x = Var(e.var)
                left = self.var_subst
                left = limits.reduce_inf_limit(left.subst(e.var, (1 / x) + e.lim), e.var, ctx)
                left = normalize(left, ctx)
                right = self.var_subst
                right = limits.reduce_inf_limit(right.subst(e.var, e.lim - (1 / x)), e.var, ctx)
                right = normalize(right, ctx)
                if left.is_evaluable() and right.is_evaluable() and expr.eval_expr(left) == expr.eval_expr(right):
                    return normalize(Limit(self.var_name, left, self.f, None), ctx)
                else:
                    return e
            return normalize(Limit(self.var_name, lim, self.f, None), ctx)
        else:
            raise TypeError


class SubstitutionInverse(Rule):
    """Apply substitution x = f(u).

    var_name - str: name of the new variable u.
    var_subst - Expr: expression containing the new variable.

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
            inv_f = solve_equation(self.var_subst, Var(e.var), new_var, ctx)
            ctx2.add_subst(new_var, inv_f)
            return ctx2
        else:
            return ctx

    def eval(self, e: Expr, ctx: Context) -> Expr:
        if not (expr.is_integral(e) or expr.is_indefinite_integral(e)):
            sep_ints = e.separate_integral()
            if len(sep_ints) == 0:
                raise RuleException("SubstitutionInverse", "no integral found in expression")
            else:
                return OnLocation(self, sep_ints[0][1]).eval(e, ctx)

        if not (expr.is_integral(e) or expr.is_indefinite_integral(e)):
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
        inv_f = solve_equation(self.var_subst, Var(e.var), new_var, ctx)
        if inv_f is None:
            raise RuleException("SubstitutionInverse", "cannot solve equation %s = %s for %s" % (
                self.var_subst, e.var, new_var
            ))

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
        if e.is_power() and expr.is_const(e.args[1]) and e.args[1].val > 1 and \
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
        # If old_expr is given, try to find it within e
        if self.old_expr is not None and self.old_expr != e:
            find_res = e.find_subexpr(self.old_expr)
            if len(find_res) == 0:
                print(e)
                raise RuleException("Rewriting", "old expression %s not found" % self.old_expr)
            loc = find_res[0]
            return OnLocation(self, loc).eval(e, ctx)

        # Now e is the old expression
        assert self.old_expr is None or self.old_expr == e

        r = Simplify()
        r1, r2 = r.eval(e, ctx), r.eval(self.new_expr, ctx)
        if r1 == r2:
            return self.new_expr

        # Rewriting 1 to sin(x)^2 + cos(x)^2
        x = Symbol("x", [VAR, CONST, OP, FUN])
        p = expr.sin(x) ** 2 + expr.cos(x) ** 2
        if e == Const(1) and expr.match(self.new_expr, p):
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

class Equation(Rule):
    """Apply substitution for equal expressions"""

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
        # If old_expr is given, try to find it within e
        if self.old_expr is not None and self.old_expr != e:
            find_res = e.find_subexpr(self.old_expr)
            if len(find_res) == 0:
                print(e)
                raise RuleException("Equation", "old expression %s not found" % self.old_expr)
            loc = find_res[0]
            return OnLocation(self, loc).eval(e, ctx)

        # Now e is the old expression
        assert self.old_expr is None or self.old_expr == e

        r = Simplify()
        r1, r2 = r.eval(e, ctx), r.eval(self.new_expr, ctx)
        if r1 == r2:
            return self.new_expr

        # Rewriting 1 to sin(x)^2 + cos(x)^2
        x = Symbol("x", [VAR, CONST, OP, FUN])
        p = expr.sin(x) ** 2 + expr.cos(x) ** 2
        if e == Const(1) and expr.match(self.new_expr, p):
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
                vars = e.get_vars(with_bd=True)
                all_index_vars = "ijklmn"
                flag = False
                for v in all_index_vars:
                    if v not in vars:
                        tmp = Limit(v, expr.POS_INF, Summation(e.index_var, e.lower, Var(v), e.body))
                        flag = True
                        break
                if not flag:
                    raise AssertionError("all variables are run out")
                if normalize(tmp, ctx) == normalize(self.new_expr, ctx):
                    return self.new_expr

            if expr.is_op(e.body) and e.body.op in '+-':
                v, l, u = e.index_var, e.lower, e.upper
                tmp = Op(e.body.op, Summation(v, l, u, e.body.args[0]), Summation(v, l, u, e.body.args[1]))
                if normalize(tmp, ctx) == normalize(self.new_expr, ctx):
                    return self.new_expr
        raise RuleException("Equation", "rewriting %s to %s failed" % (e, self.new_expr))


class IntegrationByParts(Rule):
    """Apply integration by parts.

    The arguments u and v should satisfy u * dv equals the integrand.

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
            raise RuleException("Integration by parts", "u * dv does not equal body: %s != %s" % (
                str(udv), str(e.body)))


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
        lhs = normalize(self.lhs, ctx)

        def get_coeff(t: Expr, lhs:Expr):
            """Obtain the coefficient of lhs within t."""
            if t == lhs:
                return Const(1)

            if t.is_plus():
                return get_coeff(t.args[0], lhs) + get_coeff(t.args[1], lhs)
            elif t.is_minus():
                return get_coeff(t.args[0], lhs) - get_coeff(t.args[1], lhs)
            elif expr.is_uminus(t):
                return -get_coeff(t.args[0], lhs)
            elif t.is_times():
                return t.args[0] * get_coeff(t.args[1], lhs)
            elif t.is_divides():
                return get_coeff(t.args[0], lhs) / t.args[1]
            else:
                return Const(0)

        norm_e = normalize(e, ctx)
        coeff = normalize(get_coeff(norm_e, lhs), ctx)
        lhs = self.lhs
        coeff2 = normalize(get_coeff(e, lhs), ctx)

        if coeff == Const(0) and coeff2 == Const(0):
            raise RuleException("IntegrateByEquation", "lhs %s not found in integral" % self.lhs)

        if coeff == Const(1) or coeff2 == Const(1):
            raise RuleException("IntegrateByEquation", "lhs %s has coeff 1 in integral" % self.lhs)

        if coeff == Const(0) or coeff == Const(1):
            coeff = coeff2
        res = normalize((e - (coeff * lhs)) / ((Const(1) - coeff)), ctx)

        if lhs.contains_indefinite_integral() and not lhs.contains_skolem_func():
            res = res + SkolemFunc("C", tuple())

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


def check_item(item, target=None, *, debug=False):
    """Check application of rules in the item."""
    problem = parser.parse_expr(item['problem'])

    if debug:
        print("\n%s: %s" % (item['name'], problem))

    current = problem
    prev_steps = []
    ctx = Context()

    for step in item['calc']:
        reason = step['reason']
        expected = parser.parse_expr(step['text'])

        if reason == 'Initial':
            result = current

        elif reason == 'Simplification':
            if "location" in step:
                result = OnLocation(Simplify(), step["location"]).eval(current, ctx)
            else:
                result = Simplify().eval(current, ctx)

        elif reason == 'Substitution':
            var_name = step['params']['var_name']
            f = parser.parse_expr(step['params']['f'])
            g = parser.parse_expr(step['params']['g'])
            rule = Substitution(var_name, g)
            if 'location' in step:
                result = OnLocation(rule, step['location']).eval(current, ctx)
            else:
                result = rule.eval(current, ctx)
            rule.f = parser.parse_expr(str(rule.f))  # trick to eliminate difference in printing
            if rule.f != f:
                print("Expected f: %s" % f)
                print("Actual f: %s" % rule.f)
                raise AssertionError("Unexpected value of f in substitution")

        elif reason == 'Integrate by parts':
            u = parser.parse_expr(step['params']['parts_u'])
            v = parser.parse_expr(step['params']['parts_v'])
            rule = IntegrationByParts(u, v)
            if 'location' in step:
                result = OnLocation(rule, step['location']).eval(current, ctx)
            else:
                result = rule.eval(current, ctx)

        elif reason == 'Rewrite':
            rhs = parser.parse_expr(step['params']['rhs'])
            if 'denom' in step['params']:
                rule = Equation(rhs, parser.parse_expr(step['params']['denom']))
            else:
                rule = Equation(rhs)
            if 'location' in step:
                result = OnLocation(rule, step['location']).eval(current, ctx)
            else:
                result = rule.eval(current, ctx)

        elif reason == 'Substitution inverse':
            var_name = step['params']['var_name']
            old_var = step['params']['old_var']
            g = parser.parse_expr(step['params']['g'])
            rule = SubstitutionInverse(var_name, old_var, g)
            if 'location' in step:
                result = OnLocation(rule, step['location']).eval(current, ctx)
            else:
                result = rule.eval(current, ctx)

        elif reason == 'Split region':
            c = parser.parse_expr(step['params']['c'])
            rule = SplitRegion(c)
            if 'location' in step:
                result = OnLocation(rule, step['location']).eval(current, ctx)
            else:
                result = rule.eval(current, ctx)

        elif reason == 'Solve equation':
            prev_id = int(step['params']['prev_id'])
            rule = IntegrateByEquation(prev_steps[prev_id])
            result = rule.eval(current, ctx)

        else:
            print("Reason: %s" % reason)
            raise NotImplementedError

        if result != expected:
            print("Expected: %s" % expected)
            print("Result: %s" % result)
            raise AssertionError("Error on intermediate step (%s)" % reason)

        current = result
        prev_steps.append(current)

    if target is not None:
        target = parser.parse_expr(target)
        if current != target:
            print("Target: %s" % target)
            print("Result: %s" % current)
            raise AssertionError("Error on final answer")


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
        if expr.is_fun(e) and e.func_name == self.func_name:
            for identity in ctx.get_definitions():
                if expr.is_fun(identity.lhs) and identity.lhs.func_name == self.func_name:
                    inst = expr.match(e, identity.lhs)
                    if inst == None:
                        continue
                    tmp_conds = [cond.inst_pat(inst) for cond in identity.conds.data]
                    flag = True
                    for cond in tmp_conds:
                        flag = flag and ctx.check_condition(cond)
                    if flag:
                        return normalize(identity.rhs.inst_pat(inst), ctx)
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

def is_negative_var(s: str, var: str) -> bool:
    # 将空格替换为空
    s = s.replace(" ", "")

    # 判断负号是否直接与变量相邻，或者负号在括号外作用于括号内的变量
    # 1. 负号与变量直接相邻：'-' + var
    # 2. 括号内的负号作用于变量，忽略括号
    pattern = r'(?<!\w)-' + re.escape(var) + r'(?!\w)|-\((.*' + re.escape(var) + r'.*)\)'

    # 使用正则表达式搜索
    match = re.search(pattern, s)

    # 如果匹配到符合条件的负号与变量相邻的情况，返回True
    return bool(match)

class IntExchange(Rule):
    """Exchange integral and integral"""

    def __init__(self):
        self.name = "IntExchange"

    def __str__(self):
        return "exchange integral and integral"

    def solve_var(self, f: Expr, a: Expr, x: str, ctx: Context) -> Optional[Expr]:
        from integral.expr import exprify   #  exprify 函数将非 Expr对象转换为 Expr

        if expr.is_var(f):
            if f.name == x:
                return a

        if f.is_plus():
            u, v = f.args
            if not u.contains_var(x):
                return self.solve_var(v, a - exprify(u), x, ctx)
            if not v.contains_var(x):
                return self.solve_var(u, a - exprify(v), x, ctx)

        if expr.is_uminus(f):
            u, = f.args
            return self.solve_var(u, -exprify(a), x, ctx)

        if f.is_minus():
            u, v = f.args
            if not u.contains_var(x):
                return self.solve_var(v, exprify(u) - a, x, ctx)
            if not v.contains_var(x):
                return self.solve_var(u, exprify(v) + a, x, ctx)

        if f.is_times():
            u, v = f.args
            if not u.contains_var(x) and ctx.is_nonzero(u):
                return self.solve_var(v, a / exprify(u), x, ctx)
            if not v.contains_var(x) and ctx.is_nonzero(v):
                return self.solve_var(u, a / exprify(v), x, ctx)

        if f.is_divides():
            u, v = f.args
            if not u.contains_var(x):
                rhs = exprify(u) / a
                if u.is_constant() and a in (POS_INF, NEG_INF):
                    rhs = Const(0)
                return self.solve_var(v, rhs, x, ctx)
            if not v.contains_var(x):
                return self.solve_var(u, exprify(v) * a, x, ctx)

        if f.is_power():
            u, v = f.args
            if not v.contains_var(x):
                return self.solve_var(u, a ^ (1 / exprify(v)), x, ctx)

    def exchange_int(self, evar, el, eu, svar, sl, su, sb, ctx: Context) -> Expr:
        res_list = []
        simplify = Simplify()
        # 初始化不等式列表和边界点数组
        inequality_list = []
        bp = []

        # 第一步：构建不等式列表
        # 添加外层积分上下限的不等式
        inequality_list.append(Op('<', Var(evar), eu))
        inequality_list.append(Op('>', Var(evar), el))

        # 添加内层积分上下限的不等式
        # 处理 sl < svar < su
        if is_const(sl):
            if Op('<', sl, Var(evar)) and Op('>', Var(evar), sl) not in inequality_list:
                inequality_list.append(Op('<', sl, Var(evar)))
        else:
            evar_expr = self.solve_var(sl, exprify(svar), evar, ctx)  # sl = y, svar = x, evar = y ==> y = x ==> y < x
            # 判断反解中变量是否有负号
            evar_expr_str = str(evar_expr)
            if is_negative_var(evar_expr_str, svar):
                inequality_list.append(Op('>', Var(evar), exprify(evar_expr)))  # 如果有负号，反转不等式方向
            else:
                inequality_list.append(Op('<', Var(evar), exprify(evar_expr)))

        if is_const(su):
            if Op('<', Var(evar), su) and Op('>', Var(evar), su) not in inequality_list:
                inequality_list.append(Op('<', Var(evar), su))
        else:
            evar_expr = self.solve_var(su, exprify(svar), evar, ctx)  # su = -y + 4, svar = x, evar = y ==> y = -x+4 ==> y < -x+4
            # 判断反解中变量是否有负号
            evar_expr_str = str(evar_expr)
            if is_negative_var(evar_expr_str, svar):
                inequality_list.append(Op('>', exprify(evar_expr), Var(evar)))  # 如果有负号，反转不等式方向
            else:
                inequality_list.append(Op('<', exprify(evar_expr), Var(evar)))

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
                if lhs == Var(evar):  # 左侧包含evar
                    if is_const(rhs):
                        filtered_inequalities.append((Op(op,lhs,rhs), Op(op,lhs,rhs)))
                    else:
                        # 将 evar 替换成 mid
                        evaluated_rhs = simplify.eval(exprify(rhs).subst(svar, mid), ctx)
                        # 创建新的不等式
                        new_condition = Op(op, Var(evar), evaluated_rhs)
                        filtered_inequalities.append((new_condition, Op(op,lhs,rhs)))
                elif rhs == Var(evar):  # 右侧包含evar
                    if is_const(lhs):
                        filtered_inequalities.append((Op(op,lhs,rhs), Op(op,lhs,rhs)))
                    else:
                        # 将 evar 替换成 mid
                        evaluated_lhs = simplify.eval(exprify(lhs).subst(svar, mid),ctx)
                        # 创建新的不等式
                        new_condition = Op(op, evaluated_lhs, Var(evar))
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
                    if lhs == Var(evar):
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
                    if lhs == Var(evar):
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
        if expr.is_integral(e) and expr.is_integral(e.body):
            ctx2 = body_conds(e, body_conds(e.body, ctx))
            s = e.body
            if not self.judge_novar(e.upper,e.var,s.var) or not self.judge_novar(e.lower,e.var,s.var) or \
                    self.judge_novar(s.upper,None,s.var) or self.judge_novar(s.lower,None,s.var):
                raise TypeError("Integral format error")
            elif self.judge_contains_var_and_letters(s.upper,e.var,None) or self.judge_contains_var_and_letters(s.lower,e.var,None):
                raise NotImplementedError("Contain algebraic constant and var")
            # judge whether contain oo or Algebraic Constant and no var
            elif self.judge_novar(e.upper) and self.judge_novar(e.lower) and \
                    self.judge_novar(s.upper,e.var,s.var) and self.judge_novar(s.lower,e.var,s.var):
                return Integral(s.var, s.lower, s.upper, Integral(e.var, e.lower, e.upper, s.body))
            else:
                return self.exchange_int(e.var, e.lower, e.upper, s.var, s.lower, s.upper, s.body,ctx2)
        return e

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
            raise AssertionError("SolveEquation: cannot solve")
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
        # if len(check_wellformed(ne, ctx)) == 0:
        #     return ne
        # return e



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
