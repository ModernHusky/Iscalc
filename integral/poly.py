"""Polynomials."""

from fractions import Fraction
from typing import Union, Dict, Tuple, Optional
import functools
import operator
import sympy
import math

from integral import expr, context
from integral.context import Context, apply_subterm
from integral.conditions import Conditions
from functools import lru_cache



def collect_pairs(ps, ctx:Context):
    """
    Reduce a list of pairs by collecting into groups according to
    first components, and adding the second component for each group.

    It is assumed that the first components are hashable. The second
    components are either Expr, ConstantPolynomial, or numbers. Pairs
    whose second component equals zero are removed.

    e.g.    

    - [("x", 1), ("y", 2), ("x", 3)] => [("x", 4), ("y", 2)]
    - [("x", 1), ("x", -1), ("y", 1)] => [("y", 1)]

    """
    res = {}
    for v, c in ps:
        if v in res:
            res[v] += c
        else:
            res[v] = c
    
    def zero_for(v):
        if isinstance(v, expr.Expr):
            return expr.Const(0)
        elif isinstance(v, Polynomial):
            return Polynomial(tuple())
        elif isinstance(v, (int, Fraction)):
            return 0
        else:
            raise NotImplementedError

    res_list = []
    for k, v in res.items():
        if v != zero_for(v):
            res_list.append((k, v))
    try:
        res = tuple(sorted(res_list))
    except:
        print(res_list)
        raise NotImplementedError
    return res


def collect_pairs_power(ps: Dict[expr.Expr, "Polynomial"], ctx: Context):
    res = {}
    res_list = []
    def is_non_negative(c: Polynomial):
        if c.is_fraction() and c.get_fraction() >= 0:
            return True
        e = from_poly(c)
        if expr.is_var(e) and ctx.is_not_negative(e):
            return True
        return False

    mat_list = []
    for v, c in ps:
        if v in res:
            if is_non_negative(c) and is_non_negative(res[v]):
                res[v] += c
                res[v].reduce(ctx)
            elif ctx.is_nonzero(v):
                res[v] += c
                res[v].reduce(ctx)
            else:
                res_list.append((v, c))
        else:  # v not in res
            res[v] = c

    for v, c in res.items():
        if c != Polynomial(tuple()):
            res_list.append((v, c))
    return tuple(sorted(res_list) + mat_list)

def reduce_power(n: expr.Expr, e: "Polynomial") -> Tuple[Tuple[expr.Expr, "Polynomial"]]:
    """Reduce n ^ e to normal form.
    
    Returns a list of (n_i, e_i), so that n ^ e equals (n_1 ^ e^1) * ... (n_k ^ e_k).

    If n has type Expr, n ^ e is left unchanged. If n is an integer,
    it is factored to simplify the representation.

    """
    if expr.is_const(n) and isinstance(n.val, int) and e.is_fraction():
        if n.val >= 0:
            # Compute factors of n. Let n = (n_1 ^ e_1) * ... * (n_k ^ e_k), then
            # n ^ e = (n_1 ^ (e * e_1)) * ... * (n_k ^ (e * e_k)).
            return tuple((expr.Const(ni), e * ei) for ni, ei in sympy.factorint(n.val).items())
        else:
            # If n is negative, the denominator of e must be odd.
            # If the numerator of e is also odd, add an extra -1 factor.

            # Handle negative integers
            abs_n = -n.val
            factors_abs = sympy.factorint(abs_n)
            factors_abs_tuples = tuple((expr.Const(ni), e * ei) for ni, ei in factors_abs.items())
            
            # Create the (-1)^e factor
            sqrt_neg_one = expr.Fun('sqrt', expr.Const(-1))
            e_twice = e * 2
            factor_neg = (sqrt_neg_one, e_twice)
            
            return (factor_neg,) + factors_abs_tuples

            # assert Fraction(e.get_fraction()).denominator % 2 == 1, \
            #     'reduce_power: exponent has even denominator'
            #     # TODO
            # if Fraction(e.get_fraction()).numerator % 2 == 0:
            #     return tuple((expr.Const(ni), e * ei) for ni, ei in sympy.factorint(-n.val).items())
            # else:
            #     return ((expr.Const(-1), constant(1)),) + tuple((expr.Const(ni), e * ei) for ni, ei in sympy.factorint(-n.val).items())
    else:
        return ((n, e),)

def extract_frac(ps: Tuple[Tuple[expr.Expr, "Polynomial"]]) -> Tuple[Tuple[Tuple[expr.Expr, "Polynomial"]], Fraction]:
    """Reduce (n_1 ^ e_1) * ... * (n_k ^ e_k) by extracting fractions.
    
    Collects the integer part of e_i into a separate coefficient. E.g.
    2 ^ (3/2) => 2 * 2^(1/2),
    2 ^ (-1/2) => 1/2 * 2^(1/2)
    
    """
    res = []
    coeff = 1

    for n, e in ps:
        if expr.is_const(n) and e.is_fraction():
            bval = n.val
            val = e.get_fraction()
            if val >= 1:
                coeff *= (bval ** math.floor(val))
            if val < 0:
                coeff *= Fraction(1, bval ** (-math.floor(val)))
            if val - math.floor(val) != 0:
                res.append((n, constant(val - math.floor(val))))
        else:
            res.append((n, e))

    return tuple(res), coeff


class Monomial:
    """Represents a monomial."""
    def __init__(self, coeff: Union[int, Fraction],
                 factors: Tuple[Tuple[expr.Expr, Union[int, Fraction, "Polynomial"]]]):
        """Construct a monomial from coefficient and tuple of factors,
        where each factor is associated its power. For example,

        (2, ()) -> 2
        (2, ((x, 1))) -> 2 * x
        (2, ((x, 2), (y, 1))) -> 2 * x^2 * y

        coeff: Union[int, Fraction] - coefficient of the monomial.

        """
        assert isinstance(coeff, (int, Fraction)), "Unexpected coeff: %s" % str(coeff)

        self.coeff = coeff

        # First, coerce type of exponent to Polynomial
        self.factors = []
        for base, power in factors:
            assert isinstance(base, expr.Expr)
            if isinstance(power, (int, Fraction)):
                power = constant(power)
            assert isinstance(power, Polynomial), "Unexpected power: %s" % str(power)
            self.factors.append((base, power))
        self.factors = tuple(self.factors)

    def reduce(self, ctx: Context) -> "Monomial":
        # Reduce power in factors
        reduced_factors = []
        for n, e in self.factors:
            e.reduce(ctx)
            reduced_factors.extend(reduce_power(n, e))

        # Here using collect power version
        self.factors = tuple((i, j) for i, j in collect_pairs_power(reduced_factors, ctx))

        # Extract fractions
        self.factors, coeff2 = extract_frac(self.factors)
        self.coeff = self.coeff * coeff2
        return self

    def __hash__(self):
        return hash(("MONO", self.coeff, self.factors))

    def __eq__(self, other):
        return isinstance(other, Monomial) and self.coeff == other.coeff and \
            self.factors == other.factors

    def __str__(self):
        res = ""
        if self.coeff != 1:
            res += str(self.coeff)
        for var, p in self.factors:
            s = str(var)
            if len(s) != 1:
                s = "(" + s + ")"
            if str(p) != "1":
                if isinstance(p, Fraction):
                    assert p.denominator != 0
                    if p.denominator == 1:
                        s = s + "^" + str(p.numerator)
                    else:
                        s = s + " ^ " + str(p)
                else:
                    s = s + "^" + str(p)
            if res:
                res += " * " + s
            else:
                res = s

        if res == "":
            res = "1"
        return res

    def __repr__(self):
        return "Monomial(%s)" % str(self)

    def __le__(self, other):
        if not isinstance(other, Monomial):
            raise TypeError
        return (self.factors, self.coeff) <= (other.factors, other.coeff)

    def __lt__(self, other):
        return self <= other and self != other

    def __mul__(self, other):
        if isinstance(other, (int, Fraction)):
            return Monomial(self.coeff * other, self.factors)
        elif isinstance(other, Monomial):
            return Monomial(self.coeff * other.coeff, self.factors + other.factors)
        else:
            raise NotImplementedError

    def __neg__(self):
        return Monomial(-self.coeff, self.factors)

    def __truediv__(self, other):
        if isinstance(other, Monomial):
            inv_factors = tuple((n, -e) for n, e in other.factors)
            return Monomial(Fraction(self.coeff) / other.coeff, self.factors + inv_factors)
        else:
            raise NotImplementedError

    def __pow__(self, exp):
        if isinstance(exp, int):
            return Monomial(Fraction(self.coeff) ** exp, [(n, e * exp) for n, e in self.factors])
        elif isinstance(exp, Fraction) and exp.denominator % 2 == 0:
            sqrt_factors = []
            for n, e in self.factors:
                if isinstance(n, expr.Expr) and e.is_fraction() and e.get_fraction() % 2 == 0:
                    sqrt_factors.append((expr.Fun('abs', n), e * exp))
                else:
                    sqrt_factors.append((n, e * exp))
            if self.coeff == 1:
                return Monomial(1, sqrt_factors)
            elif self.coeff < 0 and len(sqrt_factors) > 0:
                # Avoid taking square root of negative coefficient, move the negation
                # to the first sqrt_factor
                sqrt_factors = [(-sqrt_factors[0][0], sqrt_factors[0][1])] + sqrt_factors[1:]
                return Monomial(1, [(expr.Const(-self.coeff), exp)] + sqrt_factors)
            else:
                return Monomial(1, [(expr.Const(self.coeff), exp)] + sqrt_factors)
        elif isinstance(exp, Fraction) and exp.denominator % 2 == 1:
            factors = []
            for n, e in self.factors:
                factors.append((n, e * exp))
            if self.coeff == 1:
                return Monomial(1, factors)
            else:
                return Monomial(1, [(expr.Const(self.coeff), exp)] + factors)
        else:
            raise ValueError

    def is_constant(self) -> bool:
        return len(self.factors) == 0

    def get_constant(self) -> Union[int, Fraction]:
        if len(self.factors) == 0:
            return self.coeff
        else:
            raise AssertionError

    def is_fraction(self) -> bool:
        return len(self.factors) == 0

    def get_fraction(self) -> Union[int, Fraction]:
        return self.coeff


class Polynomial:
    """Represents a polynomial."""
    def __init__(self, monomials: Tuple[Monomial]):
        self.monomials = tuple(monomials)
        assert all(isinstance(mono, Monomial) for mono in self.monomials)

    def reduce(self, ctx: Context) -> "Polynomial":
        for mono in self.monomials:
            mono.reduce(ctx)
        assert all(isinstance(mono.coeff, (int, Fraction)) for mono in self.monomials)
        ts = collect_pairs([(mono.factors, mono.coeff) for mono in self.monomials], ctx)
        self.monomials = tuple(Monomial(coeff, factor) for factor, coeff in ts if coeff != 0)
        return self

    def __eq__(self, other):
        if isinstance(other, (int, Fraction)):
            return self.is_fraction() and self.get_fraction() == other

        return isinstance(other, Polynomial) and self.monomials == other.monomials

    def __le__(self, other):
        if not isinstance(other, Polynomial):
            raise TypeError
        return self.monomials <= other.monomials

    def __lt__(self, other):
        return self <= other and self != other

    def __hash__(self):
        return hash(("POLY", self.monomials))

    def __str__(self):
        if len(self.monomials) == 0:
            return "0"
        else:
            return " + ".join(str(mono) for mono in self.monomials)

    def __repr__(self):
        return "Polynomial(%s)" % str(self)

    def __len__(self):
        return len(self.monomials)

    def __add__(self, other):
        if not isinstance(other, Polynomial):
            raise TypeError
        return Polynomial(self.monomials + other.monomials)

    def __neg__(self):
        return Polynomial([-m for m in self.monomials])

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        if isinstance(other, (int, Fraction)):
            return Polynomial([m * other for m in self.monomials])
        elif isinstance(other, Polynomial):
            # Applies distributivity - could expand the number of terms exponentially
            return Polynomial([m1 * m2 for m1 in self.monomials for m2 in other.monomials])
        else:
            raise NotImplementedError

    def __truediv__(self, other):
        # Assume the denominator is a monomial
        if isinstance(other, Polynomial):
            if len(other.monomials) == 0:
                raise ZeroDivisionError
            elif len(other.monomials) == 1:
                return Polynomial([m / other.monomials[0] for m in self.monomials])
            else:
                raise ValueError
        else:
            raise NotImplementedError

    def __pow__(self, exp: "Polynomial") -> "Polynomial":
        # Assume self is a monomial and exp is a fraction
        if len(self.monomials) == 1 and isinstance(exp, (int, Fraction)):
            return Polynomial([self.monomials[0] ** exp])
        else:
            raise ValueError('%s, %s' % (self, exp))

    def is_monomial(self) -> bool:
        return len(self.monomials) == 1

    def get_monomial(self) -> Monomial:
        if self.is_monomial():
            return self.monomials[0]
        else:
            raise AssertionError("get_monomial")

    def is_fraction(self) -> bool:
        if len(self.monomials) == 0:
            return True
        return self.is_monomial() and self.get_monomial().is_fraction()

    def get_fraction(self) -> Union[int, Fraction]:
        if len(self.monomials) == 0:
            return 0
        else:
            return self.get_monomial().get_fraction()
        
    def get_constant(self) -> Union[int, Fraction]:
        """If self is a constant, return the constant. Otherwise raise an exception."""
        if len(self.monomials) == 0:
            return 0
        elif len(self.monomials) == 1 and self.monomials[0].is_constant():
            return self.monomials[0].get_constant()
        else:
            raise AssertionError

    def is_one(self) -> bool:
        return len(self.monomials) == 1 and self.monomials[0].is_one()


def constant(c: Union[int, Fraction]) -> Polynomial:
    """Polynomial for c (numerical constant)."""
    return Polynomial([Monomial(c, tuple())])

def singleton(s: expr.Expr) -> Polynomial:
    """Polynomial for 1*s^1."""
    if expr.is_const(s):
        return constant(s.val)
    else:
        return Polynomial([Monomial(1, [(s, 1)])])


"""
Conversion from expressions to polynomials.
"""

def to_poly_r(e: expr.Expr, ctx: Context) -> Polynomial:
    """Convert expression to polynomial."""
    if expr.is_var(e):
        return singleton(e)

    elif e.is_plus():
        return to_poly(e.args[0], ctx) + to_poly(e.args[1], ctx)

    elif expr.is_uminus(e):
        return -to_poly(e.args[0], ctx)

    elif e.is_minus():
        return to_poly(e.args[0], ctx) - to_poly(e.args[1], ctx)

    elif e.is_times():
        a, b = to_poly(e.args[0], ctx), to_poly(e.args[1], ctx)
        if a.is_monomial() and b.is_monomial() or a.is_fraction() or b.is_fraction():
            return a * b
        elif a.is_monomial():
            return a * singleton(from_poly(b))
        elif b.is_monomial():
            return b * singleton(from_poly(a))
        else:
            return singleton(from_poly(a)) * singleton(from_poly(b))

    elif e.is_divides():
        a, b = to_poly(e.args[0], ctx), to_poly(e.args[1], ctx)
        if a.is_fraction() and a.get_fraction() == 0:
            return constant(0)
        elif b.is_fraction() and b.get_fraction() == 1:
            return a
        elif a.is_monomial() and b.is_monomial():
            return a / b
        elif a.is_monomial():
            return a / singleton(from_poly(b))
        elif b.is_monomial():
            return singleton(from_poly(a)) / b
        else:
            return singleton(from_poly(a)) / singleton(from_poly(b))

    elif e.is_power():
        a, b = to_poly(e.args[0], ctx), to_poly(e.args[1], ctx)
        if a.is_fraction() and a.get_fraction() == 0:
            return singleton(expr.Const(0))
        elif a.is_fraction() and a.get_fraction() == 1:
            return singleton(expr.Const(1))
        elif a.is_monomial() and b.is_fraction():
            return a ** b.get_fraction()
        elif b.is_fraction():
            return Polynomial([Monomial(1, [(from_poly(a), b.get_fraction())])])
        else:
            return Polynomial([Monomial(1, [(from_poly(a), b)])])

    elif expr.is_fun(e) and e.func_name == "exp":
        a = e.args[0]
        if expr.is_fun(a) and a.func_name == "log":
            return to_poly(a.args[0], ctx)
        else:
            return Polynomial([Monomial(1, [(expr.E, to_poly(a, ctx))])])

    elif expr.is_fun(e) and e.func_name in ("sin", "cos", "tan", "cot", "csc", "sec"):
        a = e.args[0]
        if expr.is_fun(a) and e.func_name == "arc" + a.func_name:
            # sin(arcsin(x)) = x
            return to_poly(a.args[0], ctx)
        else:
            tmp = normalize(a, ctx)
            if e.func_name == "cos" and expr.is_uminus(tmp):
                return singleton(expr.Fun(e.func_name, tmp.args[0]), )
            else:
                return singleton(expr.Fun(e.func_name, tmp))

    elif expr.is_fun(e) and e.func_name in ("arcsin", "arccos", "arctan", "arccot", "arccsc", "arcsec"):
        a, = e.args
        if e.func_name in ("arcsin", "arctan", "arccot", "arccos") and expr.is_fun(a) and e.func_name == "arc" + a.func_name:
            return to_poly(a.args[0], ctx)
        else:
            return singleton(expr.Fun(e.func_name, normalize(a, ctx)))

    elif expr.is_fun(e) and e.func_name == "sqrt":
        return to_poly(expr.Op("^", e.args[0], expr.Const(Fraction(1, 2))), ctx)
    elif expr.is_fun(e) and e.func_name == "log":
        a = e.args[0]
        # log(a^b) ==> b * log(a)
        if expr.is_op(a) and a.op == '^':
            m1,m2=a.args
            return to_poly(m2, ctx) * to_poly(expr.Fun('log', m1), ctx)
        args_norm = [normalize(arg, ctx) for arg in e.args]
        return singleton(expr.Fun(e.func_name, *args_norm))
    elif expr.is_fun(e):
        args_norm = [normalize(arg, ctx) for arg in e.args]
        return singleton(expr.Fun(e.func_name, *args_norm))

    elif expr.is_evalat(e):
        if e.upper == expr.POS_INF:
            upper = expr.Limit(e.var, expr.POS_INF, e.body)
        elif e.upper == expr.NEG_INF:
            x = expr.Var(e.var)
            upper = expr.Limit(e.var, expr.POS_INF, e.body.subst(e.var, -x))
        else:
            try:
                upper = normalize(e.body.subst(e.var, e.upper), ctx)
            except:
                x = expr.Var(e.var)
                a = e.upper
                upper = expr.Limit(e.var, expr.POS_INF, e.body.subst(e.var, a - 1 / x))

        if e.lower == expr.POS_INF:
            lower = expr.Limit(e.var, expr.POS_INF, e.body)
        elif e.lower == expr.NEG_INF:
            x = expr.Var(e.var)
            lower = expr.Limit(e.var, expr.POS_INF, e.body.subst(e.var, -x))
        else:
            try:
                lower = normalize(e.body.subst(e.var, e.lower), ctx)
            except:
                x = expr.Var(e.var)
                a = e.lower
                lower = expr.Limit(e.var, expr.POS_INF, e.body.subst(e.var, a + 1 / x))
        return to_poly(normalize(upper, ctx) - normalize(lower, ctx), ctx)

    elif expr.is_integral(e):
        ctx2 = context.body_conds(e, ctx)
        body = normalize(e.body, ctx2)
        l, h = normalize(e.lower, ctx), normalize(e.upper, ctx)

        if l.is_evaluable() and h.is_evaluable() :
            ll, hh = expr.eval_expr(l), expr.eval_expr(h)
            if ll > hh:
                return singleton(-expr.Integral(e.var, h, l, body))
        return singleton(expr.Integral(e.var, normalize(e.lower, ctx), normalize(e.upper, ctx), body))
    elif expr.is_limit(e):
        ctx2 = context.body_conds(e,ctx)
        return singleton(expr.Limit(e.var, normalize(e.lim, ctx), normalize(e.body, ctx2)))
    elif expr.is_inf(e):
        if e == expr.POS_INF:
            return singleton(e)
        else:
            return -singleton(expr.POS_INF)
    elif expr.is_indefinite_integral(e):
        return singleton(expr.IndefiniteIntegral(e.var, normalize(e.body, ctx), e.skolem_args))

    elif expr.is_summation(e):
        l, u = normalize(e.lower, ctx), normalize(e.upper, ctx)
        if l == u:
            return to_poly(e.body.subst(e.index_var, l), ctx)
        ctx2 = context.body_conds(e, ctx)
        return singleton(expr.Summation(e.index_var, normalize(e.lower, ctx),
                                        normalize(e.upper, ctx), normalize(e.body, ctx2)))
    elif expr.is_product(e):
        l, u = normalize(e.lower, ctx), normalize(e.upper, ctx)
        if l == u:
            return to_poly(e.body.subst(e.index_var, l), ctx)
        ctx2 = context.body_conds(e, ctx)
        return singleton(expr.Product(e.index_var, l, u, normalize(e.body, ctx2)))
    elif expr.is_deriv(e):
        return singleton(expr.Deriv(e.var, normalize(e.body, ctx)))
    else:
        return singleton(e)

def to_poly(e: expr.Expr, ctx: Context) -> Polynomial:
    return to_poly_r(e, ctx).reduce(ctx)

def function_eval(e: expr.Expr, ctx: Context) -> expr.Expr:
    if expr.is_fun(e) and e.func_name == "binom":
        if expr.is_const(e.args[0]) and expr.is_const(e.args[1]):
            return expr.Const(math.comb(e.args[0].val, e.args[1].val))

    if expr.is_fun(e) and e.func_name == 'factorial':
        if expr.is_const(e.args[0]) and abs(round(e.args[0].val) - e.args[0].val) < 1e-15 and round(e.args[0].val) >= 0:
            return expr.Const(math.factorial(round(e.args[0].val)))

    if expr.is_fun(e) and e.func_name == 'floor':
        if expr.is_inf(e.args[0]):
            return e.args[0]
    if expr.is_fun(e) and e.func_name == 'sin':
        a = e.args[0]
        if expr.is_uminus(a):
            return -expr.Fun('sin', a.args[0])
    return e

def function_table(e: expr.Expr, ctx: Context) -> expr.Expr:
    if not expr.is_fun(e) or len(e.args) != 1:
        return e

    func_table = ctx.get_function_tables()
    if not e.func_name in func_table:
        return e
    if not e.args[0].is_constant():
        return e
    if e.args[0] in func_table[e.func_name]:
        return func_table[e.func_name][e.args[0]]
    else:
        return e

def simplify_identity(e: expr.Expr, ctx: Context) -> expr.Expr:
    for identity in ctx.get_simp_identities():
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
                return identity.rhs.inst_pat(inst)
    return e

def simplify_eq(e: expr.Expr, ctx: Context) -> expr.Expr:
    if not expr.is_var(e):
        return e

    for eq in ctx.get_conds().data:
        if eq.is_equals() and e == eq.lhs and eq.rhs.is_constant():
            return eq.rhs
    return e

def simplify_limit(e: expr.Expr, ctx: Context) -> expr.Expr:
    from integral import limits
    if not expr.is_limit(e):
        return e
    if e.var not in e.body.get_vars():
        return e.body

    if e.lim == expr.POS_INF:
        return limits.reduce_inf_limit(e.body, e.var, ctx)
    elif e.lim == expr.NEG_INF:
        raise limits.reduce_neg_inf_limit(e.body, e.var, ctx)
    else:
        return limits.reduce_finite_limit(e, ctx)

def simplify_integral(e: expr.Expr, ctx: Context) -> expr.Expr:
    if not expr.is_integral(e):
        return e

    if expr.is_deriv(e.body) and e.body.var == e.var:
        return expr.EvalAt(e.var, e.lower, e.upper, e.body.body)
    elif e.lower == e.upper:
        return expr.Const(0)
    else:
        return e

def simplify_power(e: expr.Expr, ctx: Context) -> expr.Expr:
    if not e.is_power():
        return e
    if e.args[1].is_plus() and expr.is_const(e.args[0]) and expr.is_const(e.args[1].args[1]):
        # c1 ^ (a + c2) => c1 ^ c2 * c1 ^ a
        return (e.args[0] ^ e.args[1].args[1]) * (e.args[0] ^ e.args[1].args[0])
    elif e.args[1].is_minus() and expr.is_const(e.args[0]) and expr.is_const(e.args[1].args[1]):
        # c1 ^ (a - c2) => c1 ^ -c2 * c1 ^ a
        return (e.args[0] ^ e.args[1].args[0]) * (e.args[0] ^ (-(e.args[1].args[1])))
    elif expr.is_uminus(e.args[0]):
        if expr.is_const(e.args[1]):
            # (-a) ^ n = (-1) ^ n * a ^ n
            return (expr.Const(-1) ^ e.args[1]) * (e.args[0].args[0] ^ e.args[1])
        else:
            return e
    elif e.args[0].is_minus() and expr.is_uminus(e.args[0].args[0]) and expr.is_const(e.args[1]):
        # (-a - b) ^ n = (-1) ^ n * (a + b) ^ n
        nega, negb = e.args[0].args
        return (expr.Const(-1) ^ e.args[1]) * ((nega.args[0] + negb) ^ e.args[1])
    else:
        return e

def simplify_trig(e: expr.Expr, ctx: Context) -> expr.Expr:
    if not (expr.is_fun(e) and e.func_name in ('sin', 'cos', 'tan', 'cot', 'csc', 'sec')):
        return e

    a = e.args[0]
    c = expr.Symbol('c', [expr.CONST])
    d = expr.Symbol('d', [expr.CONST])

    # Find constant coefficient of a
    coeff = None
    if expr.match(a, c * expr.pi):
        coeff = a.args[0].val
    elif expr.match(a, c * expr.pi / d):
        coeff = Fraction(a.args[0].args[0].val, a.args[1].val)
    elif expr.match(a, expr.pi / d):
        coeff = Fraction(1, a.args[1].val)
    elif expr.match(a, -(c * expr.pi)):
        coeff = -a.args[0].args[0].val
    elif expr.match(a, -(c * expr.pi / d)):
        coeff = Fraction(-a.args[0].args[0].args[0].val, a.args[0].args[1].val)
    elif expr.match(a, -(expr.pi / d)):
        coeff = Fraction(-1, a.args[0].args[1].val)
    elif a == -expr.pi:
        coeff = -1

    if coeff is None:
        return e

    # Normalize coefficient to the interval (-1, 1]
    n = int(coeff)
    if n % 2 == 0:
        coeff = coeff - n
    elif n > 0:
        coeff = coeff - (n + 1)
    else:
        coeff = coeff - (n - 1)

    def build(val):
        if isinstance(val, int):
            return expr.Fun(e.func_name, expr.Const(val) * expr.pi)
        elif isinstance(val, Fraction):
            return expr.Fun(e.func_name, expr.Const(val.numerator) * expr.pi / expr.Const(val.denominator))
        else:
            raise TypeError

    # Further normalize to [0, pi]
    if coeff < 0:
        if e.func_name in ('sin', 'tan', 'cot', 'csc'):
            return -build(-coeff)
        else:  # cos, sec: does not change sign
            return build(-coeff)
    else:
        return build(coeff)

def simplify_log(e: expr.Expr, ctx: Context) -> expr.Expr:
    if not (expr.is_fun(e) and e.func_name == 'log'):
        return e
    
    a = e.args[0]
    if not a.is_constant():
        return e

    if expr.is_const(a) and isinstance(a.val, int):
        int_factors = sympy.factorint(a.val)
        log_ints = []
        for b, e in int_factors.items():
            if e != 1:
                log_ints.append(e * expr.log(expr.Const(b)))
            else:
                log_ints.append(expr.log(expr.Const(b)))
        return sum(log_ints[1:], log_ints[0])
    elif expr.is_const(a) and isinstance(a.val, Fraction):
        return expr.log(expr.Const(a.val.numerator)) - expr.log(expr.Const(a.val.denominator))
    elif a.is_times():
        return expr.log(a.args[0]) + expr.log(a.args[1])
    elif a.is_divides():
        return expr.log(a.args[0]) - expr.log(a.args[1])
    elif expr.is_fun(a) and a.func_name == 'sqrt':
        return expr.log(a.args[0]) / 2
    else:
        return e

def simplify_sqrt(e: expr.Expr, ctx: Context) -> expr.Expr:
    if not (expr.is_fun(e) and e.func_name == 'sqrt'):
        return e

    if e.args[0] == expr.Const(0):
        return expr.Const(0)
    if e.args[0] == expr.Const(1):
        return expr.Const(1)
    if e.args[0] == expr.Const(4):
        return expr.Const(2)
    if e.args[0] == expr.Const(Fraction(1, 4)):
        return expr.Const(Fraction(1, 2))
    if e.args[0] == expr.Const(Fraction(1, 2)):
        return 1 / expr.sqrt(expr.Const(2))

    return e

def simplify_sum(e: expr.Expr, ctx:Context) -> expr.Expr:
    if expr.is_summation(e):
        if e.lower == e.upper and e.lower not in (expr.POS_INF, expr.NEG_INF):
            return e.body
        if e.body == expr.Const(0):
            return e.body
    elif e.is_plus():
        a = expr.Symbol('a',[expr.SUMMATION])
        b = expr.Symbol('b',[expr.SUMMATION])
        pat_list = [a+b, -(a)+b, a - b, -(a) - b]
        for pat in pat_list:
            inst = expr.match(e, pat)
            if inst != None:
                tmpa:expr.Summation = a.inst_pat(inst)
                tmpb:expr.Summation = b.inst_pat(inst)
                al, bl, au, bu, av, bv = tmpa.lower,tmpb.lower,tmpa.upper,tmpb.upper, tmpa.index_var, tmpb.index_var
                ab, bb = tmpa.body, tmpb.body
                if al == bl and au == bu and ab == bb:
                    new_inst = {'a':ab, 'b':bb}
                    new_body = pat.inst_pat(new_inst)
                    return expr.Summation(av, al, au, new_body)
    return e

def simplify_inf(e: expr.Expr, ctx: Context) -> expr.Expr:
    if e.is_plus():
        if e.args[0] == expr.POS_INF and e.args[1].is_constant():
            return e.args[0]
    elif e.is_minus():
        if e.args[0] == expr.POS_INF and e.args[1].is_constant():
            return e.args[0]
    elif e.is_divides():
        if e.args[0] == expr.POS_INF and e.args[1].is_constant():
            if expr.eval_expr(e) != 0:
                return e.args[0]
    return e

def simplify_skolem(e:expr.Expr, ctx:Context):
    if expr.is_skolem_func(e) or expr.is_const(e) or expr.is_var(e) or expr.is_inf(e):
        return e
    elif expr.is_op(e):
        op = e.op
        if op == '*':
            a, b = e.args
            if expr.is_const(a) and expr.is_skolem_func(b):
                if a.val != 0:
                    return b
                else:
                    return expr.Const(0)
            elif expr.is_const(b) and expr.is_skolem_func(a):
                if b.val != 0:
                    return a
                else:
                    return expr.Const(0)
        elif op == '/':
            a, b = e.args
            if expr.is_const(b) and expr.is_skolem_func(a):
                if b.val != 0:
                    return a
                else:
                    raise ValueError("The denominator cannot be 0")
        return expr.Op(op, *[simplify_skolem(arg, ctx) for arg in e.args])
    elif expr.is_fun(e):
        return expr.Fun(e.func_name, *[simplify_skolem(arg, ctx) for arg in e.args])
    elif expr.is_summation(e):
        return expr.Summation(e.index_var, simplify_skolem(e.lower,ctx), simplify_skolem(e.upper,ctx), simplify_skolem(e.body,ctx))
    elif expr.is_limit(e):
        return expr.Limit(e.var, simplify_skolem(e.lim, ctx), simplify_skolem(e.body, ctx))
    elif expr.is_integral(e):
        e:expr.Integral
        return expr.Integral(e.var, simplify_skolem(e.lower,ctx), simplify_skolem(e.upper,ctx), simplify_skolem(e.body,ctx))
    return e
def normal_const(e:expr.Expr, ctx:Context):
    if e.is_constant():
        return normalize(e, ctx)
    elif expr.is_var(e) or expr.is_inf(e) or expr.is_skolem_func(e):
        return e
    elif expr.is_op(e):
        return expr.Op(e.op, *[normal_const(arg, ctx) for arg in e.args])
    elif expr.is_fun(e):
        return expr.Fun(e.func_name, *[normal_const(arg, ctx) for arg in e.args])
    elif expr.is_summation(e):
        return expr.Summation(e.index_var, normal_const(e.lower,ctx), \
                              normal_const(e.upper,ctx), normal_const(e.body,ctx))
    elif expr.is_limit(e):
        return expr.Limit(e.var, normal_const(e.lim, ctx), normal_const(e.body, ctx))
    elif expr.is_integral(e):
        e:expr.Integral
        return expr.Integral(e.var, normal_const(e.lower,ctx), normal_const(e.upper,ctx),\
                             normal_const(e.body,ctx))
    raise NotImplementedError(str(e))

def normalize(e: expr.Expr, ctx: Context) -> expr.Expr:
    if e.is_equals():
        return expr.Eq(normalize(e.lhs, ctx), normalize(e.rhs, ctx))

    if expr.is_const(e):
        return e
    for i in range(5):
        old_e = e
        e = from_poly(to_poly(e, ctx))
        e = apply_subterm(e, function_table, ctx)
        e = apply_subterm(e, function_eval, ctx)
        e = apply_subterm(e, simplify_identity, ctx)
        e = apply_subterm(e, simplify_eq, ctx)
        e = apply_subterm(e, simplify_limit, ctx)
        e = apply_subterm(e, simplify_integral, ctx)
        e = apply_subterm(e, simplify_power, ctx)
        e = apply_subterm(e, simplify_trig, ctx)
        e = apply_subterm(e, simplify_log, ctx)
        e = apply_subterm(e, simplify_sqrt, ctx)
        e = apply_subterm(e, simplify_inf, ctx)
        e = apply_subterm(e, simplify_sum, ctx)
        e = apply_subterm(e, simplify_skolem, ctx)
        if is_complex_expr(e, ctx.get_conds()):
            e = simplify_complex(e, ctx.get_conds())
        if e == old_e:
            break

    return e

def is_complex_expr(e: expr.Expr, conds: Conditions) -> bool:
    """判断表达式是否为复数表达式"""
    if expr.is_fun(e):
        if e.func_name in ("Re", "Im", "conj", "abs"):
            return True
    elif expr.is_var(e):
        # 检查变量是否为i
        if str(e) == 'i':
            return True
        # 检查是否被声明为复数
        for cond in conds.data:
            if expr.is_fun(cond) and cond.func_name == "isComplex" and cond.args[0] == e:
                return True
    elif expr.is_op(e):
        # 检查是否包含虚数单位i
        if e.is_times() and any(str(arg) == 'i' for arg in e.args):
            return True
        # 递归检查子表达式
        return any(is_complex_expr(arg, conds) for arg in e.args)
    return False

def evaluate_const_expr(e: expr.Expr) -> expr.Expr:
    """尝试计算常量表达式的值"""
    if expr.is_const(e):
        return e
    elif expr.is_op(e):
        args = [evaluate_const_expr(arg) for arg in e.args]
        if e.is_plus():
            a, b = args[0], args[1]
            if expr.is_const(a) and expr.is_const(b):
                return expr.Const(a.val + b.val)
            return expr.Op("+", a, b)
        elif e.is_minus() and len(e.args) == 2:
            a, b = args[0], args[1]
            if expr.is_const(a) and expr.is_const(b):
                return expr.Const(a.val - b.val)
            return expr.Op("-", a, b)
        elif e.is_times():
            a, b = args[0], args[1]
            if expr.is_const(a) and expr.is_const(b):
                return expr.Const(a.val * b.val)
            # 处理平方根的乘法
            if expr.is_fun(a) and a.func_name == 'sqrt' and expr.is_fun(b) and b.func_name == 'sqrt':
                return expr.sqrt(a.args[0] * b.args[0])
            return expr.Op("*", a, b)
        elif e.is_divides():
            a, b = args[0], args[1]
            if expr.is_const(a) and expr.is_const(b):
                from fractions import Fraction
                new_val = Fraction(a.val) / Fraction(b.val)
                return expr.Const(new_val)
            return expr.Op("/", a, b)
        elif e.is_power():
            a, b = args[0], args[1]
            if expr.is_const(a) and expr.is_const(b):
                return expr.Const(a.val ** b.val)
            return expr.Op("^", a, b)
        else:
            return expr.Op(e.op, *args)
    elif expr.is_fun(e):
        args = [evaluate_const_expr(arg) for arg in e.args]
        if all(expr.is_const(arg) for arg in args):
            if e.func_name == "sqrt":
                from math import sqrt
                val = args[0].val
                result = sqrt(val)
                # 如果结果是整数
                if result.is_integer():
                    return expr.Const(int(result))
                else:
                    return expr.Fun("sqrt", expr.Const(val))
        return expr.Fun(e.func_name, *args)
    return e

def get_complex_parts(e: expr.Expr) -> tuple[Optional[expr.Expr], Optional[expr.Expr]]:
    """获取复数表达式的实部和虚部
    返回 (实部, 虚部),如果某部分为0则返回None
    """
    print("get_complex_parts",e.__repr__())
    print("get_complex_parts",e)
    if str(e) == 'i':
        return None, expr.Const(1)
    elif expr.is_const(e):
        return e if e.val != 0 else None, None
    elif expr.is_op(e):
        if e.is_times():
            print("e.args[0],e.args[1]",e.args[0],e.args[1])
            # 处理形如 a*i 的情况
            if str(e.args[1]) == 'i':
                print("e.args[0]",e.args[0])
                return None, e.args[0]
            elif str(e.args[0]) == 'i':
                print("e.args[1]",e.args[1])
                return None, e.args[1]
            elif len(str(e.args[0])) == 1 and expr.is_var(e.args[1]) and str(e.args[1]) != 'i':
                print("e.args[1]",e.args[1])
                return None, e
            elif len(str(e.args[1])) == 1 and expr.is_var(e.args[0]) and str(e.args[0]) != 'i':
                print("e.args[0]",e.args[0])
                return None, e
            else:
                
                
                # 处理一般乘法
                re1, im1 = get_complex_parts(e.args[0])
                re2, im2 = get_complex_parts(e.args[1])
                # 将None转换为Const(0)进行计算
                re1 = re1 if re1 is not None else expr.Const(0)
                re2 = re2 if re2 is not None else expr.Const(0)
                im1 = im1 if im1 is not None else expr.Const(0)
                im2 = im2 if im2 is not None else expr.Const(0)
                
                # (a + bi)(c + di) = (ac - bd) + (ad + bc)i
                re = evaluate_const_expr(re1 * re2 - im1 * im2)
                im = evaluate_const_expr(re1 * im2 + im1 * re2)
                
                # 如果结果为0则返回None
                return (re if not (expr.is_const(re) and re.val == 0) else None,
                       im if not (expr.is_const(im) and im.val == 0) else None)
                
        elif expr.is_uminus(e):
            # 处理负号
            re, im = get_complex_parts(e.args[0])
            if re is None and im is None:
                return None, None
            elif re is None:
                return None, evaluate_const_expr(-im)
            elif im is None:
                return evaluate_const_expr(-re), None
            else:
                return evaluate_const_expr(-re), evaluate_const_expr(-im)
            
        elif e.is_minus() and len(e.args) == 2:
            # 处理减法
            re1, im1 = get_complex_parts(e.args[0])
            re2, im2 = get_complex_parts(e.args[1])
            # 分别处理实部和虚部的减法
            re_diff = None
            im_diff = None

            # 计算实部
            if re1 is None:
                if re2 is None:
                    re_diff = expr.Const(0)
                else:
                    re_diff = -re2
            elif re2 is None:
                re_diff = re1
            elif expr.is_const(re1) and re1.val == 0:
                if expr.is_const(re2) and re2.val == 0:
                    re_diff = expr.Const(0)
                else:
                    re_diff = -re2
            elif expr.is_const(re2) and re2.val == 0:
                re_diff = re1
            else:
                re_diff = evaluate_const_expr(re1 - re2)

            # 计算虚部
            if im1 is None:
                if im2 is None:
                    im_diff = expr.Const(0)
                else:
                    im_diff = -im2
            elif im2 is None:
                im_diff = im1
            elif expr.is_const(im1) and im1.val == 0:
                if expr.is_const(im2) and im2.val == 0:
                    im_diff = expr.Const(0)
                else:
                    im_diff = -im2
            elif expr.is_const(im2) and im2.val == 0:
                im_diff = im1
            else:
                im_diff = evaluate_const_expr(im1 - im2)

            # 如果结果为0则返回None
            re_result = None if expr.is_const(re_diff) and re_diff.val == 0 else re_diff
            im_result = None if expr.is_const(im_diff) and im_diff.val == 0 else im_diff
            
            return (re_result if not (expr.is_const(re_result) and re_result.val == 0) else None,
                   im_result if not (expr.is_const(im_result) and im_result.val == 0) else None)
            
            
        elif e.is_plus():
            # 处理加法
            re1, im1 = get_complex_parts(e.args[0])
            re2, im2 = get_complex_parts(e.args[1])
            
            # 分别处理实部和虚部的减法
            re_diff = None
            im_diff = None

            # 计算实部
            if re1 is None:
                if re2 is None:
                    re_diff = expr.Const(0)
                else:
                    re_diff = re2
            elif re2 is None:
                re_diff = re1
            elif expr.is_const(re1) and re1.val == 0:
                if expr.is_const(re2) and re2.val == 0:
                    re_diff = expr.Const(0)
                else:
                    re_diff = re2
            elif expr.is_const(re2) and re2.val == 0:
                re_diff = re1
            else:
                re_diff = evaluate_const_expr(re1 + re2)

            # 计算虚部
            if im1 is None:
                if im2 is None:
                    im_diff = expr.Const(0)
                else:
                    im_diff = im2
            elif im2 is None:
                im_diff = im1
            elif expr.is_const(im1) and im1.val == 0:
                if expr.is_const(im2) and im2.val == 0:
                    im_diff = expr.Const(0)
                else:
                    im_diff = im2
            elif expr.is_const(im2) and im2.val == 0:
                im_diff = im1
            else:
                im_diff = evaluate_const_expr(im1 + im2)

            # 如果结果为0则返回None
            re_result = None if expr.is_const(re_diff) and re_diff.val == 0 else re_diff
            im_result = None if expr.is_const(im_diff) and im_diff.val == 0 else im_diff
            
            return (re_result if not (expr.is_const(re_result) and re_result.val == 0) else None,
                   im_result if not (expr.is_const(im_result) and im_result.val == 0) else None)
            
            
    # 其他情况，视为实数
    return e, None

def make_complex(re: Optional[expr.Expr], im: Optional[expr.Expr]) -> expr.Expr:
    print("!!! re, im  make_complex",re, im)
    """从实部和虚部构造复数表达式"""
    # 如果实部和虚部都为None，返回0
    if re is None and im is None:
        return expr.Const(0)
        
    # 如果虚部为None，只返回实部
    if im is None:
        print("make_complex:re",re)
        return re
        
    # 如果实部为None
    if re is None:
        if expr.is_const(im):
            if im.val == 1:
                return expr.Var('i')
            elif im.val == -1:
                return -expr.Var('i')
            return im * expr.Var('i')
        elif expr.is_var(im):  # 处理单个变量的情况
            if str(im) == 'i':
                return im
            return im * expr.Var('i')
        elif expr.is_op(im) and im.is_times():  # 处理乘法表达式
            if any(expr.is_var(arg) and str(arg) != 'i' for arg in im.args):
                print("!!1 im",im)
                return im
            return im * expr.Var('i')
        return im * expr.Var('i')
        
    # 处理负的虚部
    if expr.is_const(im) and im.val < 0:
        if im.val == -1:
            return re - expr.Var('i')
        neg_im = expr.Const(-im.val)
        return re - neg_im * expr.Var('i')
        
    # 一般情况
    if expr.is_var(im):  # 处理单个变量的情况
        if str(im) == 'i':
            return re + im
        return re + im * expr.Var('i')
    elif expr.is_op(im) and im.is_times():  # 处理乘法表达式
        if any(expr.is_var(arg) and str(arg) != 'i' for arg in im.args):
            print("!!2 im",im)
            return re + im
        return re + im * expr.Var('i')
    elif expr.is_const(im) and im.val == 1:
        return re + expr.Var('i')
    return re + im * expr.Var('i')

@lru_cache(maxsize=None)
def abs(arg: expr.Expr) -> expr.Expr:
    """计算复数的绝对值"""
    if str(arg) == 'i':
        return expr.Const(1)  # |i| = 1
    elif expr.is_const(arg):
        return expr.abs(arg)  # |c| = |c| for real constant c
    
    # 获取实部和虚部
    re_part, im_part = get_complex_parts(arg)
    re_part = evaluate_const_expr(re_part)
    im_part = evaluate_const_expr(im_part)
        
    # 计算平方和
    if expr.is_const(re_part) and expr.is_const(im_part):
        square_sum = re_part.val * re_part.val + im_part.val * im_part.val
        return expr.Fun("sqrt", expr.Const(square_sum))
    else:
        square_sum = evaluate_const_expr(re_part * re_part + im_part * im_part)
        return expr.Fun("sqrt", square_sum)

@lru_cache(maxsize=None)
def simplify_complex(e: expr.Expr, conds: Conditions) -> expr.Expr:
    print("simplify_complex",e.__repr__())
    """化简复数表达式"""
    # 转换为可哈希类型以便缓存
    def make_hashable(e):
        if isinstance(e, expr.Expr):
            return expr_to_hashable(e)
        elif isinstance(e, Conditions):
            return tuple(expr_to_hashable(cond) for cond in e.data)
        return e
    
    def helper(e_tuple, conds_tuple):
        # 将元组转回表达式
        e = expr_from_hashable(e_tuple)
        
        if expr.is_const(e):
            return e
        elif expr.is_var(e):
            if str(e) == 'i':
                return e
            # 检查变量是否有定义
            for cond_tuple in conds_tuple:
                cond = expr_from_hashable(cond_tuple)
                if expr.is_op(cond) and cond.op == "=" and cond.lhs == e:
                    return helper(make_hashable(cond.rhs), conds_tuple)
            return e
        elif expr.is_fun(e):
            if e.func_name == "abs":
                arg = helper(make_hashable(e.args[0]), conds_tuple)
                return abs(arg)
                
            elif e.func_name == "Re":
                arg = helper(make_hashable(e.args[0]), conds_tuple)
                re_part, _ = get_complex_parts(arg)
                return re_part if re_part is not None else expr.Const(0)
                
            elif e.func_name == "Im":
                arg = helper(make_hashable(e.args[0]), conds_tuple)
                _, im_part = get_complex_parts(arg)
                return im_part if im_part is not None else expr.Const(0)
                
            elif e.func_name == "conj":
                arg = helper(make_hashable(e.args[0]), conds_tuple)
                if expr.is_const(arg):
                    return arg
                elif str(arg) == 'i':
                    return -expr.Var('i')
            
                re_part, im_part = get_complex_parts(arg)
                print("!!! re_part, im_part",re_part, im_part)
                # 如果只有实部，直接返回
                if im_part is None:
                    return re_part
                # 如果只有虚部，返回其相反数
                if re_part is None:
                    return make_complex(None, evaluate_const_expr(-im_part))
                # 如果都有，返回实部减虚部
                return make_complex(re_part, evaluate_const_expr(-im_part))
                
            else:
                return expr.Fun(e.func_name, *(helper(make_hashable(arg), conds_tuple) for arg in e.args))
                
        elif expr.is_op(e):
            if e.op == "+":
                print("e.op == +",e)
                a = helper(make_hashable(e.args[0]), conds_tuple)
                b = helper(make_hashable(e.args[1]), conds_tuple)
                
                # 检查常量情况
                if expr.is_const(a) and expr.is_const(b):
                    return evaluate_const_expr(expr.Op("+", a, b))
                
                # 检查是否有一个操作数为0
                if expr.is_const(a) and a.val == 0:
                    return b
                if expr.is_const(b) and b.val == 0:
                    return a
                    
                # 检查是否都是纯虚数
                if str(a) == 'i' and str(b) == 'i':
                    return expr.Const(2) * expr.Var('i')
                
                print("+ a, b",a, b)
                # 获取实部和虚部
                re_a, im_a = get_complex_parts(a)
                re_b, im_b = get_complex_parts(b)
                print("+ re_a, im_a",re_a, "+++",im_a)
                print("+ re_b, im_b",re_b, "+++",im_b)
                # 分别处理实部和虚部的加法
                re_sum = None
                im_sum = None

                # 计算实部
                if re_a is None:
                    if re_b is None:
                        re_sum = expr.Const(0)
                    else:
                        re_sum = re_b
                elif re_b is None:
                    re_sum = re_a
                elif expr.is_const(re_a) and re_a.val == 0:
                    if expr.is_const(re_b) and re_b.val == 0:
                        re_sum = expr.Const(0)
                    else:
                        re_sum = re_b
                elif expr.is_const(re_b) and re_b.val == 0:
                    re_sum = re_a
                else:
                    re_sum = evaluate_const_expr(re_a + re_b)

                # 计算虚部
                if im_a is None:
                    if im_b is None:
                        im_sum = expr.Const(0)
                    else:
                        im_sum = im_b
                elif im_b is None:
                    im_sum = im_a
                elif expr.is_const(im_a) and im_a.val == 0:
                    if expr.is_const(im_b) and im_b.val == 0:
                        im_sum = expr.Const(0)
                    else:
                        im_sum = im_b
                elif expr.is_const(im_b) and im_b.val == 0:
                    im_sum = im_a
                else:
                    im_sum = evaluate_const_expr(im_a + im_b)
                print("+ re_sum, im_sum",re_sum, "+++",im_sum)

                # 如果结果为0则返回None
                re_result = None if expr.is_const(re_sum) and re_sum.val == 0 else re_sum
                im_result = None if expr.is_const(im_sum) and im_sum.val == 0 else im_sum
                
                print("+ re_result, im_sum",re_result, "+++",im_result)
                return make_complex(re_result, im_result)
                
            elif e.op == "-" and len(e.args) == 2:
                print("e.op == -",e)
                a = helper(make_hashable(e.args[0]), conds_tuple)
                b = helper(make_hashable(e.args[1]), conds_tuple)
                print("a, b",a, b)
                # 检查常量情况
                if expr.is_const(a) and expr.is_const(b):
                    return evaluate_const_expr(expr.Op("-", a, b))
                
                # 检查特殊情况
                if expr.is_const(b) and b.val == 0:
                    return a
                if a == b:
                    return expr.Const(0)
                    
                # 检查是否都是纯虚数
                if str(a) == 'i' and str(b) == 'i':
                    return expr.Const(0)
                
                re_a, im_a = get_complex_parts(a)
                re_b, im_b = get_complex_parts(b)
                
                # 分别处理实部和虚部的减法
                re_diff = None
                im_diff = None

                # 计算实部
                if re_a is None:
                    if re_b is None:
                        re_diff = expr.Const(0)
                    else:
                        re_diff = -re_b
                elif re_b is None:
                    re_diff = re_a
                elif expr.is_const(re_a) and re_a.val == 0:
                    if expr.is_const(re_b) and re_b.val == 0:
                        re_diff = expr.Const(0)
                    else:
                        re_diff = -re_b
                elif expr.is_const(re_b) and re_b.val == 0:
                    re_diff = re_a
                else:
                    re_diff = evaluate_const_expr(re_a - re_b)

                # 计算虚部
                if im_a is None:
                    if im_b is None:
                        im_diff = expr.Const(0)
                    else:
                        im_diff = -im_b
                elif im_b is None:
                    im_diff = im_a
                elif expr.is_const(im_a) and im_a.val == 0:
                    if expr.is_const(im_b) and im_b.val == 0:
                        im_diff = expr.Const(0)
                    else:
                        im_diff = -im_b
                elif expr.is_const(im_b) and im_b.val == 0:
                    im_diff = im_a
                else:
                    im_diff = evaluate_const_expr(im_a - im_b)

                # 如果结果为0则返回None
                re_result = None if expr.is_const(re_diff) and re_diff.val == 0 else re_diff
                im_result = None if expr.is_const(im_diff) and im_diff.val == 0 else im_diff
                print("- re_diff, im_diff",re_result, "---",im_result)
                return make_complex(re_result, im_result)
                
            elif e.op == "*":
                print("e.op == *",e)
                a = helper(make_hashable(e.args[0]), conds_tuple)
                b = helper(make_hashable(e.args[1]), conds_tuple)
                
                # 检查常量情况
                if expr.is_const(a) and expr.is_const(b):
                    return evaluate_const_expr(expr.Op("*", a, b))
                
                # 检查特殊情况
                if expr.is_const(a):
                    if a.val == 0:
                        return expr.Const(0)
                    if a.val == 1:
                        return b
                    if a.val == -1:
                        return -b
                if expr.is_const(b):
                    if b.val == 0:
                        return expr.Const(0)
                    if b.val == 1:
                        return a
                    if b.val == -1:
                        return -a
                        
                # 检查是否涉及i
                if str(a) == 'i':
                    if str(b) == 'i':
                        return expr.Const(-1)
                    re_b, im_b = get_complex_parts(b)
                    if im_b is None and re_b is not None:
                        return make_complex(im_b, re_b)
                    elif im_b is None and re_b is not None:
                        return make_complex(im_b, re_b)
                    return make_complex(-im_b, re_b)
                if str(b) == 'i':
                    re_a, im_a = get_complex_parts(a)
                    if im_a is None and re_a is None:
                        return make_complex(im_a, re_a)
                    elif im_a is None and re_a is not None:
                        return make_complex(im_a, re_a)
                    return make_complex(-im_a, re_a)
                
                # 获取实部和虚部
                re_a, im_a = get_complex_parts(a)
                re_b, im_b = get_complex_parts(b)

                re_a = re_a if re_a is not None else expr.Const(0)
                im_a = im_a if im_a is not None else expr.Const(0)
                re_b = re_b if re_b is not None else expr.Const(0)
                im_b = im_b if im_b is not None else expr.Const(0)
                
                # 计算实部和虚部的乘积
                # 优化乘法计算
                re_prod = None
                im_prod = None

                # 计算实部 (re_a * re_b - im_a * im_b)
                if expr.is_const(re_a) and re_a.val == 0:
                    if expr.is_const(im_a) and im_a.val == 0:
                        re_prod = expr.Const(0)
                    else:
                        re_prod = evaluate_const_expr(-im_a * im_b)
                elif expr.is_const(im_a) and im_a.val == 0:
                    re_prod = evaluate_const_expr(re_a * re_b)
                elif expr.is_const(re_b) and re_b.val == 0:
                    if expr.is_const(im_b) and im_b.val == 0:
                        re_prod = expr.Const(0)
                    else:
                        re_prod = evaluate_const_expr(-im_a * im_b)
                elif expr.is_const(im_b) and im_b.val == 0:
                    re_prod = evaluate_const_expr(re_a * re_b)
                else:
                    re_prod = evaluate_const_expr(re_a * re_b - im_a * im_b)

                # 计算虚部 (re_a * im_b + im_a * re_b)
                if expr.is_const(re_a) and re_a.val == 0:
                    if expr.is_const(im_a) and im_a.val == 0:
                        im_prod = expr.Const(0)
                    else:
                        im_prod = evaluate_const_expr(im_a * re_b)
                elif expr.is_const(im_a) and im_a.val == 0:
                    im_prod = evaluate_const_expr(re_a * im_b)
                elif expr.is_const(re_b) and re_b.val == 0:
                    if expr.is_const(im_b) and im_b.val == 0:
                        im_prod = expr.Const(0)
                    else:
                        im_prod = evaluate_const_expr(im_a * re_b)
                elif expr.is_const(im_b) and im_b.val == 0:
                    im_prod = evaluate_const_expr(re_a * im_b)
                else:
                    im_prod = evaluate_const_expr(re_a * im_b + im_a * re_b)
                print("* re_prod, im_prod",re_prod, "***",im_prod)
                
                # 如果结果为0则返回None
                re_result = None if expr.is_const(re_prod) and re_prod.val == 0 else re_prod
                im_result = None if expr.is_const(im_prod) and im_prod.val == 0 else im_prod
                
                print("* re_result, im_result",re_result, "***",im_result)
                return make_complex(re_result, im_result)
                
            elif e.op == "/":
                print("e.op == /",e)
                num = helper(make_hashable(e.args[0]), conds_tuple)
                denom = helper(make_hashable(e.args[1]), conds_tuple)
                print("1/:num, denom",num, denom)
                
                # 检查常量情况
                # 如果分母是实数
                if not is_complex_expr(denom, conds):
                    # 简单除法
                    print("2/:num, denom",num, denom)
                    return expr.Op("/", num, denom)
                
                # 检查特殊情况
                if expr.is_const(num):
                    if num.val == 0:
                        return expr.Const(0)  # 如果分子为0，直接返回0，不需要进行后续计算
                    if num.val == 1 and expr.is_const(denom) and denom.val == 1:
                        return expr.Const(1)
                    if num.val == -1 and expr.is_const(denom) and denom.val == 1:
                        return expr.Const(-1)
                if expr.is_const(denom):
                    if denom.val == 0:
                        raise ValueError("Division by zero in complex division")
                    if denom.val == 1:
                        return num
                    if denom.val == -1:
                        return -num
                if num == denom:
                    return expr.Const(1)
                    
                # 检查是否涉及i
                if str(denom) == 'i':
                    re_num, im_num = get_complex_parts(num)
                    re_num = re_num if re_num is not None else expr.Const(0)
                    im_num = im_num if im_num is not None else expr.Const(0)
                    return make_complex(im_num, -re_num)
                
                # 获取实部和虚部
                re_num, im_num = get_complex_parts(num)
                re_denom, im_denom = get_complex_parts(denom)
                re_num = re_num if re_num is not None else expr.Const(0)
                im_num = im_num if im_num is not None else expr.Const(0)
                re_denom = re_denom if re_denom is not None else expr.Const(0)
                im_denom = im_denom if im_denom is not None else expr.Const(0)

                # 计算分母的模的平方
                denom_squared = None
                if expr.is_const(re_denom) and re_denom.val == 0:
                    if expr.is_const(im_denom) and im_denom.val == 0:
                        raise ValueError("Division by zero in complex division")
                    denom_squared = evaluate_const_expr(im_denom * im_denom)
                elif expr.is_const(im_denom) and im_denom.val == 0:
                    denom_squared = evaluate_const_expr(re_denom * re_denom)
                else:
                    denom_squared = evaluate_const_expr(re_denom * re_denom + im_denom * im_denom)

                # 计算实部和虚部
                real_part = None
                imag_part = None

                # 计算实部 (re_num * re_denom + im_num * im_denom) / denom_squared
                if expr.is_const(re_num) and re_num.val == 0:
                    if expr.is_const(im_num) and im_num.val == 0:
                        real_part = expr.Const(0)
                    else:
                        real_part = evaluate_const_expr(im_num * im_denom / denom_squared)
                elif expr.is_const(im_num) and im_num.val == 0:
                    real_part = evaluate_const_expr(re_num * re_denom / denom_squared)
                elif expr.is_const(re_denom) and re_denom.val == 0:
                    if expr.is_const(im_denom) and im_denom.val == 0:
                        raise ValueError("Division by zero in complex division")
                    real_part = evaluate_const_expr(im_num * im_denom / denom_squared)
                elif expr.is_const(im_denom) and im_denom.val == 0:
                    real_part = evaluate_const_expr(re_num * re_denom / denom_squared)
                else:
                    real_part = evaluate_const_expr((re_num * re_denom + im_num * im_denom) / denom_squared)

                # 计算虚部 (im_num * re_denom - re_num * im_denom) / denom_squared
                if expr.is_const(re_num) and re_num.val == 0:
                    if expr.is_const(im_num) and im_num.val == 0:
                        imag_part = expr.Const(0)
                    else:
                        imag_part = evaluate_const_expr(im_num * re_denom / denom_squared)
                elif expr.is_const(im_num) and im_num.val == 0:
                    imag_part = evaluate_const_expr(-re_num * im_denom / denom_squared)
                elif expr.is_const(re_denom) and re_denom.val == 0:
                    if expr.is_const(im_denom) and im_denom.val == 0:
                        raise ValueError("Division by zero in complex division")
                    imag_part = evaluate_const_expr(-im_num * re_denom / denom_squared)
                elif expr.is_const(im_denom) and im_denom.val == 0:
                    imag_part = evaluate_const_expr(im_num * re_denom / denom_squared)
                else:
                    imag_part = evaluate_const_expr((im_num * re_denom - re_num * im_denom) / denom_squared)

                # 如果结果为0则返回None
                re_result = None if expr.is_const(real_part) and real_part.val == 0 else real_part
                im_result = None if expr.is_const(imag_part) and imag_part.val == 0 else imag_part
                print("re_result, im_result",re_result, "///",im_result)
                return make_complex(re_result, im_result)
                
            elif e.op == "^" and expr.is_const(e.args[1]):
                base = helper(make_hashable(e.args[0]), conds_tuple)
                power = e.args[1].val
                
                # 检查常量情况
                if expr.is_const(base):
                    if base.val == 0:
                        return expr.Const(0)
                    if base.val == 1:
                        return expr.Const(1)
                    if power == 0:
                        return expr.Const(1)
                    if power == 1:
                        return base
                    
                # 处理i的幂
                if str(base) == 'i':
                    if isinstance(power, int):
                        remainder = power % 4
                        if remainder == 0:
                            return expr.Const(1)
                        elif remainder == 1:
                            return expr.Var('i')
                        elif remainder == 2:
                            return expr.Const(-1)
                        else:  # remainder == 3
                            return -expr.Var('i')
                
                # 处理一般情况
                if isinstance(power, int) and power > 0:
                    result = base
                    for _ in range(power - 1):
                        result = helper(make_hashable(expr.Op("*", result, base)), conds_tuple)
                    return result
                elif isinstance(power, int) and power == 0:
                    return expr.Const(1)
                elif isinstance(power, int) and power < 0:
                    inv_base = helper(make_hashable(expr.Op("/", expr.Const(1), base)), conds_tuple)
                    result = inv_base
                    for _ in range(-power - 1):
                        result = helper(make_hashable(expr.Op("*", result, inv_base)), conds_tuple)
                    return result
                elif isinstance(power, Fraction):
                    if power.denominator == 2:
                        re_part, im_part = get_complex_parts(base)
                        if re_part is None and im_part is None:
                            return expr.Const(0)
                        
                        re_part = re_part if re_part is not None else expr.Const(0)
                        im_part = im_part if im_part is not None else expr.Const(0)
                        
                        r = evaluate_const_expr(expr.sqrt(re_part * re_part + im_part * im_part))
                        theta = expr.Fun("atan2", im_part, re_part)
                        
                        new_r = evaluate_const_expr(expr.sqrt(r))
                        new_theta = evaluate_const_expr(theta / expr.Const(2))
                        
                        new_re = evaluate_const_expr(new_r * expr.Fun("cos", new_theta))
                        new_im = evaluate_const_expr(new_r * expr.Fun("sin", new_theta))
                        
                        return make_complex(new_re, new_im)
                    else:
                        raise NotImplementedError(f"Complex power with fraction {power} is not supported")
                else:
                    raise NotImplementedError(f"Complex power with exponent type {type(power)} is not supported")
        return e
    
    # 转换输入为可哈希类型并调用helper函数
    e_tuple = make_hashable(e)
    conds_tuple = make_hashable(conds)
    return helper(e_tuple, conds_tuple)

"""
Conversion from polynomials to terms.
"""

def rsize(e: expr.Expr) -> int:
    """Find size of term without constants."""
    if expr.is_const(e):
        return 0
    elif expr.is_uminus(e):
        return rsize(e.args[0])
    elif e.is_times():
        return rsize(e.args[0]) + e.args[1].size() + 1
    else:
        return e.size()

def display_large(e: expr.Expr) -> bool:
    """Determine whether the expression requires large display."""
    def pred(e: expr.Expr) -> bool:
        return expr.is_integral(e) or e.is_divides() or (expr.is_fun(e) and e.func_name == "binom")
    return len(e.find_subexpr_pred(pred)) > 0

def from_mono(m: Monomial) -> expr.Expr:
    """Convert a monomial to an expression."""
    sign = 1
    num_factors = []
    denom_factors = []
    if isinstance(m.coeff, int) or (isinstance(m.coeff, Fraction) and m.coeff.denominator == 1):
        if m.coeff == 1:
            pass
        elif m.coeff == -1:
            sign = -1
        elif m.coeff >= 0:
            num_factors.append(expr.Const(m.coeff))
        else:
            sign = -1
            num_factors.append(expr.Const(-m.coeff))
    elif isinstance(m.coeff, Fraction):
        denom_factors.append(expr.Const(m.coeff.denominator))
        if m.coeff.numerator == 1:
            pass
        elif m.coeff.numerator == -1:
            sign = -1
        elif m.coeff.numerator >= 0:
            num_factors.append(expr.Const(m.coeff.numerator))
        else:
            sign = -1
            num_factors.append(expr.Const(-m.coeff.numerator))
    else:
        raise TypeError

    for base, power in m.factors:
        if isinstance(power, Polynomial) and power.is_fraction():
            power = power.get_fraction()
        if isinstance(base, expr.Expr) and base == expr.E:
            if isinstance(power, Polynomial):
                num_factors.append(expr.exp(from_poly(power)))
            else:
                num_factors.append(expr.exp(expr.Const(power)))
        else:
            if power == 1:
                num_factors.append(base)
            elif power == Fraction(1 / 2):
                num_factors.append(expr.sqrt(base))
            elif power == -1:
                denom_factors.append(base)
            elif power == Fraction(-1 / 2):
                denom_factors.append(expr.sqrt(base))
            elif isinstance(power, (int, Fraction)) and power > 0:
                num_factors.append(base ** expr.Const(power))
            elif isinstance(power, (int, Fraction)) and power < 0:
                denom_factors.append(base ** expr.Const(-power))
            elif isinstance(power, Polynomial):
                num_factors.append(base ** from_poly(power))
            else:
                raise TypeError("from_mono: unexpected type %s for power" % type(power))

    large_num_factors = [factor for factor in num_factors if display_large(factor)]
    num_factors = [factor for factor in num_factors if not display_large(factor)]

    def prod(factors):
        if len(factors) == 0:
            return expr.Const(1)
        else:
            return functools.reduce(operator.mul, factors[1:], factors[0])

    res = prod(num_factors)
    if len(denom_factors) != 0:
        res = res / prod(denom_factors)
    if len(large_num_factors) != 0:
        if res == expr.Const(1):
            res = prod(large_num_factors)
        else:
            res = res * prod(large_num_factors)
    if sign == -1:
        res = -res
    return res

def from_poly(p: Polynomial) -> expr.Expr:
    """Convert a polynomial to an expression."""
    if len(p.monomials) == 0:
        return expr.Const(0)
    else:
        monos = [from_mono(m) for m in p.monomials]
        monos = sorted(monos, key=lambda p: rsize(p), reverse=True)
        res = monos[0]
        for mono in monos[1:]:
            if expr.is_uminus(mono):
                res = res - mono.args[0]
            elif mono.is_times() and expr.is_uminus(mono.args[0]):
                res = res - mono.args[0].args[0] * mono.args[1]
            elif mono.is_times() and expr.is_const(mono.args[0]) and mono.args[0].val < 0:
                res = res - expr.Const(-mono.args[0].val) * mono.args[1]
            elif expr.is_const(mono) and mono.val < 0:
                res = res - expr.Const(-mono.val)
            else:
                res = res + mono
        return res

def expr_to_hashable(e: expr.Expr) -> tuple:
    """将表达式转换为可哈希的元组形式"""
    if expr.is_const(e):
        return ('CONST', e.val)
    elif expr.is_var(e):
        return ('VAR', str(e))
    elif expr.is_fun(e):
        return ('FUN', e.func_name, tuple(expr_to_hashable(arg) for arg in e.args))
    elif expr.is_op(e):
        return ('OP', e.op, tuple(expr_to_hashable(arg) for arg in e.args))
    elif expr.is_limit(e):
        return ('LIMIT', str(e.var), expr_to_hashable(e.lim), expr_to_hashable(e.body))
    elif expr.is_integral(e):
        return ('INTEGRAL', str(e.var), expr_to_hashable(e.lower), expr_to_hashable(e.upper), expr_to_hashable(e.body))
    elif expr.is_deriv(e):
        return ('DERIV', str(e.var), expr_to_hashable(e.body))
    else:
        return ('OTHER', str(e))

def expr_from_hashable(t: tuple) -> expr.Expr:
    """将可哈希的元组形式转换回表达式"""
    if t[0] == 'CONST':
        return expr.Const(t[1])
    elif t[0] == 'VAR':
        return expr.Var(t[1])
    elif t[0] == 'FUN':
        args = tuple(expr_from_hashable(arg) for arg in t[2])
        return expr.Fun(t[1], *args)
    elif t[0] == 'OP':
        args = tuple(expr_from_hashable(arg) for arg in t[2])
        return expr.Op(t[1], *args)
    elif t[0] == 'LIMIT':
        return expr.Limit(t[1], expr_from_hashable(t[2]), expr_from_hashable(t[3]))
    elif t[0] == 'INTEGRAL':
        return expr.Integral(t[1], expr_from_hashable(t[2]), expr_from_hashable(t[3]), expr_from_hashable(t[4]))
    elif t[0] == 'DERIV':
        return expr.Deriv(t[1], expr_from_hashable(t[2]))
    else:
        raise ValueError(f"Unknown expression type: {t[0]}")

