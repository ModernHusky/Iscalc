"""Normalization of expressions."""

import math
from typing import Union
from fractions import Fraction

from integral.expr import Expr, Const, Integral, Var
from integral import poly, expr
from integral.poly import Polynomial, from_poly, to_poly
from integral.context import Context


def unfold_power(p: Polynomial, n: int, ctx: Context) -> Polynomial:
    """Unfold power of a polynomial."""
    assert n >= 0
    if n == 0:
        return poly.constant(1)

    res = p
    for i in range(n-1):
        res = (p * res).reduce(ctx)
    return res

class NormalQuotient:
    def __init__(self, num: Polynomial, denom: Polynomial, ctx: Context):
        self.num = to_poly(from_poly(num), ctx)
        self.denom = to_poly(from_poly(denom), ctx)
    
    def __str__(self):
        return "(%s, %s)" % (from_poly(self.num), from_poly(self.denom))

    def to_expr(self) -> Expr:
        denom = from_poly(self.denom)
        if denom == Const(1):
            return from_poly(self.num)
        else:
            return from_poly(self.num) / denom

def add_normal_quotient(n1: NormalQuotient, n2: NormalQuotient, ctx: NormalQuotient) -> NormalQuotient:
    num = (n1.num * n2.denom + n1.denom * n2.num).reduce(ctx)
    denom = (n1.denom * n2.denom).reduce(ctx)
    return NormalQuotient(num, denom, ctx)

def uminus_normal_quotient(n: NormalQuotient, ctx: Context) -> NormalQuotient:
    return NormalQuotient((-n.num).reduce(ctx), n.denom, ctx)

def minus_normal_quotient(n1: NormalQuotient, n2: NormalQuotient, ctx: Context) -> NormalQuotient:
    return add_normal_quotient(n1, uminus_normal_quotient(n2, ctx), ctx)

def mult_normal_quotient(n1: NormalQuotient, n2: NormalQuotient, ctx: Context) -> NormalQuotient:
    return NormalQuotient((n1.num * n2.num).reduce(ctx), (n1.denom * n2.denom).reduce(ctx), ctx)

def inverse_normal_quotient(n: NormalQuotient, ctx: Context) -> NormalQuotient:
    return NormalQuotient(n.denom, n.num, ctx)

def divide_normal_quotient(n1: NormalQuotient, n2: NormalQuotient, ctx: Context) -> NormalQuotient:
    return mult_normal_quotient(n1, inverse_normal_quotient(n2, ctx), ctx)

def exp_normal_quotient(base: NormalQuotient, val: int, ctx: Context):
    if val == 1:
        return base
    elif val == -1:
        return inverse_normal_quotient(base, ctx)
    elif isinstance(val, int):
        if val >= 0:
            return NormalQuotient(unfold_power(base.num, val, ctx),
                                  unfold_power(base.denom, val, ctx), ctx)
        else:
            return NormalQuotient(unfold_power(base.denom, -val, ctx),
                                  unfold_power(base.num, -val, ctx), ctx)
    elif isinstance(val, Fraction):
        if val >= 0:
            return NormalQuotient(poly.singleton(from_poly(base.num) ** val),
                                  poly.singleton(from_poly(base.denom) ** val), ctx)
        else:
            return NormalQuotient(poly.singleton(from_poly(base.denom) ** -val),
                                  poly.singleton(from_poly(base.num) ** -val), ctx)
    else:
        raise TypeError

def equal_normal_quotient(n1: NormalQuotient, n2: NormalQuotient, ctx: Context) -> bool:
    e1 = from_poly((n1.num * n2.denom).reduce(ctx))
    e2 = from_poly((n1.denom * n2.num).reduce(ctx))
    return e1 == e2

def normalize_quotient(e: Expr, ctx: Context) -> NormalQuotient:
    def rec(e: Expr) -> NormalQuotient:
        if e.is_plus():
            return add_normal_quotient(rec(e.args[0]), rec(e.args[1]), ctx)
        elif expr.is_uminus(e):
            return uminus_normal_quotient(rec(e.args[0]), ctx)
        elif e.is_minus():
            return minus_normal_quotient(rec(e.args[0]), rec(e.args[1]), ctx)
        elif e.is_times():
            return mult_normal_quotient(rec(e.args[0]), rec(e.args[1]), ctx)
        elif e.is_divides():
            if all(isinstance(arg, expr.Fun) and arg.func_name == "factorial" for arg in e.args):
                f1,f2 = e.args
                n1,n2 = f1.args[0], f2.args[0]
                if n1 == n2:
                    return NormalQuotient(poly.constant(1), poly.constant(1), ctx)
                try:
                    v = expr.eval_expr(poly.normalize(n1 - n2, ctx))
                    if v > 0:
                        res = Const(1)
                        for i in range(1, v+1):
                            res *= Const(i)
                        res_e = divide_normal_quotient(rec(res), rec(Const(1)), ctx)
                    else:
                        v = abs(v)
                        res = Const(1)
                        for i in range(1, v+1):
                            res *= (n1 + Const(i))
                        res_e = divide_normal_quotient(rec(Const(1)), rec(res), ctx)
                except NotImplementedError as e:
                    res = expr.Fun("factorial", n1 - n2)
                    res_e = divide_normal_quotient(rec(res), rec(Const(1)), ctx)
                finally:
                    return res_e
            return divide_normal_quotient(rec(e.args[0]), rec(e.args[1]), ctx)
        elif e.is_power():
            if expr.is_const(e.args[1]):
                return exp_normal_quotient(rec(e.args[0]), e.args[1].val, ctx)
        elif expr.is_fun(e):
            if e.func_name == 'sqrt':
                return rec(e.args[0] ** Const(Fraction(1,2)))
            elif e.func_name == "tan":
                return NormalQuotient(poly.singleton(expr.sin(e.args[0])), poly.singleton(expr.cos(e.args[0])), ctx)
            elif e.func_name == "cot":
                return NormalQuotient(poly.singleton(expr.cos(e.args[0])), poly.singleton(expr.sin(e.args[0])), ctx)
            elif e.func_name == "sec":
                return NormalQuotient(poly.constant(1), poly.singleton(expr.cos(e.args[0])), ctx)
            elif e.func_name == "csc":
                return NormalQuotient(poly.constant(1), poly.singleton(expr.sin(e.args[0])), ctx)
            else:
                e = expr.Fun(e.func_name, *(quotient_normalize(arg, ctx) for arg in e.args))

        # Un-handled cases
        return NormalQuotient(poly.singleton(e), poly.constant(1), ctx)

    return rec(e)

def quotient_normalize(t: Expr, ctx: Context) -> Expr:
    return normalize_quotient(t, ctx).to_expr()

def eq_quotient(t1: Expr, t2: Expr, ctx: Context) -> bool:
    n1 = normalize_quotient(t1, ctx)
    n2 = normalize_quotient(t2, ctx)
    return equal_normal_quotient(n1, n2, ctx)


class NormalPower:
    """A more general normal form for polynomials.
    
    NormalPower(p, q, m) represents the expression (p/q)^(1/m), where
    p and q are polynomials, and m is an integer >= 1.
    
    """
    def __init__(self, num: Polynomial, denom: Polynomial, root: int, ctx: Context):
        self.num = num
        self.denom = denom
        self.root = root
        self.ctx = ctx

    def __str__(self):
        return "(%s, %s, %s)" % (self.num, self.denom, self.root)

    def to_expr(self) -> Expr:
        denom = from_poly(self.denom)
        if denom == Const(1):
            inner = from_poly(self.num)
        else:
            inner = from_poly(self.num) / denom
        if self.root == Const(1):
            return inner
        else:
            return inner ** Const(Fraction(1, self.root))

def add_normal_power(n1: NormalPower, n2: NormalPower) -> NormalPower:
    """Add two normal forms.
    
    If both sides do not have roots, take common denominators.
    Not much is done otherwise.

    """
    if n1.root == 1 and n2.root == 1:
        num = n1.num * n2.denom + n1.denom * n2.num
        denom = n1.denom * n2.denom
        return NormalPower(num, denom, 1, n1.ctx)
    elif n1.root == 1 and n2.root > 1:
        # p/q + y^(1/n) = (p + q * y^(1/n)) / q
        num = n1.num + n1.denom * poly.singleton(n2.to_expr())
        denom = n1.denom
        return NormalPower(num, denom, 1, n1.ctx)
    elif n1.root > 1 and n2.root == 1:
        return add_normal_power(n2, n1)
    else:
        return NormalPower(poly.singleton(n1.to_expr()) + poly.singleton(n2.to_expr()),
                           poly.constant(1), 1, n1.ctx)

def uminus_normal_power(n: NormalPower) -> NormalPower:
    """Negation of a normal form.
    
    If argument has roots, not much is done.

    """
    if n.root == 1:
        return NormalPower(-n.num, n.denom, 1, n.ctx)
    else:
        return NormalPower(-poly.singleton(n.to_expr()), n.denom, 1, n.ctx)

def minus_normal_power(n1: NormalPower, n2: NormalPower) -> NormalPower:
    return add_normal_power(n1, uminus_normal_power(n2))

def mult_normal_power(n1: NormalPower, n2: NormalPower) -> NormalPower:
    """Multiply two normal forms.
    
    If the two sides are (p/q)^(1/m) and (r/s)^(1/n), take the
    lcm of m and n to be k, then the product is
    
      (p^(k/m) * r^(k/n) / q^(k/m) * s^(k/n)) ^ (1/k)

    """
    root = math.lcm(n1.root, n2.root)
    p1 = root // n1.root
    p2 = root // n2.root
    num = unfold_power(n1.num, p1, n1.ctx) * unfold_power(n2.num, p2, n1.ctx)
    denom = unfold_power(n1.denom, p1, n1.ctx) * unfold_power(n2.denom, p2, n1.ctx)
    return NormalPower(num, denom, root, n1.ctx)

def inverse_normal_power(n: NormalPower) -> NormalPower:
    """Inverse of normal form.
    
    The inverse of (p/q)^(1/m) is (q/p)^(1/m).

    """
    return NormalPower(n.denom, n.num, n.root, n.ctx)

def divide_normal_power(n1: NormalPower, n2: NormalPower) -> NormalPower:
    return mult_normal_power(n1, inverse_normal_power(n2))

def exp_normal_power(base: NormalPower, val: Union[int, Fraction]) -> NormalPower:
    if val == 1:
        return base
    elif val == -1:
        return inverse_normal_power(base)
    elif isinstance(val, int):
        if val >= 0:
            return NormalPower(unfold_power(base.num, val, base.ctx),
                               unfold_power(base.denom, val, base.ctx),
                               base.root, base.ctx)
        else:
            return NormalPower(unfold_power(base.denom, -val, base.ctx),
                               unfold_power(base.num, -val, base.ctx),
                               base.root, base.ctx)
    elif isinstance(val, Fraction):
        if val >= 0:
            return NormalPower(unfold_power(base.num, val.numerator, base.ctx),
                               unfold_power(base.denom, val.numerator, base.ctx),
                               base.root * val.denominator, base.ctx)
        else:
            return NormalPower(unfold_power(base.denom, -val.numerator, base.ctx),
                               unfold_power(base.num, -val.numerator, base.ctx),
                               base.root * val.denominator, base.ctx)
    else:
        raise TypeError

def equal_normal_power(n1: NormalPower, n2: NormalPower) -> bool:
    e1 = from_poly((n1.num * n2.denom).reduce(n1.ctx))
    e2 = from_poly((n1.denom * n2.num).reduce(n1.ctx))
    return n1.root == n2.root and e1 == e2

def normalize_power(e: Expr, ctx: Context) -> NormalPower:
    def rec(e: Expr) -> NormalPower:
        if e.is_plus():
            return add_normal_power(rec(e.args[0]), rec(e.args[1]))
        elif expr.is_uminus(e):
            return uminus_normal_power(rec(e.args[0]))
        elif e.is_minus():
            return minus_normal_power(rec(e.args[0]), rec(e.args[1]))
        elif e.is_times():
            return mult_normal_power(rec(e.args[0]), rec(e.args[1]))
        elif e.is_divides():
            return divide_normal_power(rec(e.args[0]), rec(e.args[1]))
        elif e.is_power():
            if expr.is_const(e.args[1]):
                return exp_normal_power(rec(e.args[0]), e.args[1].val)
        elif expr.is_fun(e):
            if e.func_name == 'sqrt':
                return rec(e.args[0] ** Const(Fraction(1,2)))

        # Un-handled cases
        return NormalPower(poly.singleton(e), poly.constant(1), 1, ctx)

    return rec(e)

def power_normalize(t: Expr, ctx: Context) -> Expr:
    return normalize_power(t, ctx).to_expr()

def eq_power(t1: Expr, t2: Expr, ctx: Context) -> bool:
    n1 = normalize_power(t1, ctx)
    n2 = normalize_power(t2, ctx)
    return equal_normal_power(n1, n2)


class NormalLog:
    """NormalLog(e) represents an expression exp(e), where e is a polynomial."""
    def __init__(self, e: Polynomial):
        self.e = e

    def __str__(self):
        return "(%s)" % from_poly(self.e)

    def to_expr(self) -> Expr:
        return from_poly(self.e)

def minus_normal_log(a: NormalLog, b: NormalLog) -> NormalLog:
    return NormalLog(a.e / b.e)

def add_normal_log(a: NormalLog, b: NormalLog) -> NormalLog:
    return NormalLog(a.e * b.e)

def normalize_log(e: Expr, ctx: Context) -> NormalLog:
    def rec(e: Expr) -> NormalLog:
        if e.is_minus():
            return minus_normal_log(rec(e.args[0]), rec(e.args[1]))
        elif e.is_plus():
            return add_normal_log(rec(e.args[0]), rec(e.args[1]))
        elif expr.is_fun(e) and e.func_name == 'log':
            return NormalLog(poly.singleton(e.args[0]))
        return NormalLog(poly.singleton(expr.Fun("exp", e)))

    return rec(e)

def equal_normal_log(t1: NormalLog, t2: NormalLog):
    e1 = from_poly(t1.e)
    e2 = from_poly(t2.e)
    return e1 == e2

def eq_log(t1: Expr, t2: Expr, ctx: Context) -> bool:
    n1 = normalize_log(t1, ctx)
    n2 = normalize_log(t2, ctx)
    return equal_normal_log(n1, n2)

def normalize_exp(t: Expr) -> Expr:
    def rec(t):
        if expr.is_fun(t) and t.func_name == 'exp':
            a = t.args[0]
            if a.is_plus():
                return rec(expr.exp(a.args[0])) * rec(expr.exp(a.args[1]))
            elif a.is_minus():
                return rec(expr.exp(a.args[0])) / rec(expr.exp(a.args[1]))
            elif expr.is_uminus(a):
                return rec(expr.exp(a.args[0])) ** (-1)
            elif a.is_divides() and expr.is_fun(a.args[0]) and a.args[0].func_name == 'log':
                return rec(a.args[0].args[0] ** (1 / a.args[1]))
            elif expr.is_fun(a) and a.func_name == 'log':
                return a.args[0]
            else:
                return t
        else:
            return t

    return rec(t)

class NormalDefiniteIntegral:
    def __init__(self, var: str, lower: Polynomial, upper: Polynomial, body: Polynomial, ctx: Context):
        self.var = var
        self.lower = to_poly(from_poly(lower), ctx)
        self.upper = to_poly(from_poly(upper), ctx)
        self.body = to_poly(from_poly(body), ctx)

    def __str__(self):
        return "INT(%s, %s, %s, %s)" % (self.var, from_poly(self.lower),from_poly(self.upper),from_poly(self.body))


def equal_normal_definite_integral(t1: NormalDefiniteIntegral, t2: NormalDefiniteIntegral):
    e1 = from_poly(t1.body), from_poly(t1.lower), from_poly(t1.upper)
    e2 = from_poly(t2.body), from_poly(t2.lower), from_poly(t2.upper)
    return e1 == e2

def add_normal_definite_integral(t1: NormalDefiniteIntegral, t2: NormalDefiniteIntegral, ctx: Context):
    tmp = from_poly(t2.body)
    tmp = to_poly(tmp.subst(t2.var, expr.Var(t1.var)), ctx)
    return NormalDefiniteIntegral(t1.var, t1.lower, t1.upper, t1.body + tmp, ctx)

def minus_normal_definite_integral(t1: NormalDefiniteIntegral, t2: NormalDefiniteIntegral, ctx: Context):
    tmp = from_poly(t2.body)
    tmp = to_poly(tmp.subst(t2.var, expr.Var(t1.var)), ctx)
    return NormalDefiniteIntegral(t1.var, t1.lower, t1.upper, t1.body - tmp, ctx)

def normalize_definite_integral(e: Expr, ctx: Context):
    def rec(e: Expr) -> NormalDefiniteIntegral:
        if e.is_plus():
            return add_normal_definite_integral(rec(e.args[0]), rec(e.args[1]), ctx)
        elif e.is_minus():
            return minus_normal_definite_integral(rec(e.args[0]), rec(e.args[1]), ctx)
        elif expr.is_integral(e):
            return NormalDefiniteIntegral(e.var, to_poly(e.lower, ctx), to_poly(e.upper, ctx), to_poly(e.body, ctx), ctx)
        else:
            return NormalDefiniteIntegral("_x", to_poly(Const(0), ctx), to_poly(Const(1), ctx), to_poly(e, ctx), ctx)
    return rec(e)

def eq_definite_integral(t1: Expr, t2: Expr, ctx: Context) -> bool:
    n1 = normalize_definite_integral(t1, ctx)
    n2 = normalize_definite_integral(t2, ctx)
    return equal_normal_definite_integral(n1, n2)

def is_odd(e, var, conds) -> bool:
    from integral import poly
    tmp1 = e
    tmp2 = e.subst(var, -Var(var))
    if poly.normalize(tmp1 + tmp2, conds) == Const(0):
        return True
    return False

def simp_definite_integral(e: Integral, ctx: Context) -> Expr:
    if not expr.is_integral(e):
        return e
    if is_odd(e.body, e.var, ctx) and \
        from_poly(to_poly(e.lower, ctx)) == from_poly(to_poly(expr.Op("-", e.upper), ctx)):
        return Const(0)
    return e
