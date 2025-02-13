"""Functions for solving equations"""

from typing import Optional, Tuple

from integral import expr
from integral.expr import Expr, POS_INF, NEG_INF, Const, Var
from integral.poly import normalize
from integral.context import Context


def solve_equation(f: Expr, a: Expr, x: str, ctx: Context) -> Optional[Expr]:
    """Solve the equation f(x) = a for the variable x.
    
    First, try to isolate x on the left side by moving expressions
    independent of x to the right side.

    Next, several other heuristics are tried, such as using linearity.

    Note: for trigonometric or other special functions, may only produce
    one of the solutions.

    """
    if expr.is_var(f):
        if f.name == x:
            return a
    if f.is_plus():
        u, v = f.args
        if not u.contains_var(x):
            # u + v = a  ==>  v = a - u
            return solve_equation(v, a - u, x, ctx)
        if not v.contains_var(x):
            # u + v = a  ==>  u = a - v
            return solve_equation(u, a - v, x, ctx)
    if expr.is_uminus(f):
        # -u = a  ==>  u = -a
        u, = f.args
        return solve_equation(u, -a, x, ctx)
    if f.is_minus():
        u, v = f.args
        if not u.contains_var(x):
            # u - v = a  ==>  v = u - a
            return solve_equation(v, u - a, x, ctx)
        if not v.contains_var(x):
            # u - v = a  ==>  u = v + a
            return solve_equation(u, v + a, x, ctx)
    if f.is_times():
        u, v = f.args
        if not u.contains_var(x) and ctx.is_nonzero(u):
            # u * v = a  ==>  v = a / u
            return solve_equation(v, a / u, x, ctx)
        if not v.contains_var(x) and ctx.is_nonzero(v):
            # u * v = a  ==>  u = a / v
            return solve_equation(u, a / v, x, ctx)
    if f.is_divides():
        u, v = f.args
        if not u.contains_var(x):
            # u / v = a  ==>  v = a / u
            rhs = u / a
            if u.is_constant() and a in (POS_INF, NEG_INF):
                rhs = Const(0)
            return solve_equation(v, rhs, x, ctx)
        if not v.contains_var(x):
            # u / v = a  ==>  u = v * a
            return solve_equation(u, v * a, x, ctx)
    if f.is_power():
        u, v = f.args
        if not v.contains_var(x):
            # u ^ v = a  ==>  u = a ^ (1/v)
            return solve_equation(u, a ^ (1/v), x, ctx)
    if expr.is_fun(f):
        if f.func_name == "log":
            return solve_equation(f.args[0], expr.exp(a), x, ctx)
        elif f.func_name == "exp":
            return solve_equation(f.args[0], expr.log(a), x, ctx)
        elif f.func_name == "sin":
            return solve_equation(f.args[0], expr.arcsin(a), x, ctx)
        elif f.func_name == "cos":
            return solve_equation(f.args[0], expr.arccos(a), x, ctx)
        elif f.func_name == "tan":
            return solve_equation(f.args[0], expr.arctan(a), x, ctx)
        elif f.func_name == "cot":
            return solve_equation(f.args[0], expr.arccot(a), x, ctx)
        elif f.func_name == "sec":
            return solve_equation(f.args[0], expr.arcsec(a), x, ctx)
        elif f.func_name == "csc":
            return solve_equation(f.args[0], expr.arccsc(a), x, ctx)
        elif f.func_name == "arcsin":
            return solve_equation(f.args[0], expr.sin(a), x, ctx)
        elif f.func_name == "arccos":
            return solve_equation(f.args[0], expr.cos(a), x, ctx)
        elif f.func_name == "arctan":
            return solve_equation(f.args[0], expr.tan(a), x, ctx)
        elif f.func_name == "arccot":
            return solve_equation(f.args[0], expr.cot(a), x, ctx)
        elif f.func_name == "arcsec":
            return solve_equation(f.args[0], expr.sec(a), x, ctx)
        elif f.func_name == "arccsc":
            return solve_equation(f.args[0], expr.csc(a), x, ctx)
        elif f.func_name == "sqrt":
            return solve_equation(f.args[0], a ^ 2, x, ctx)

    # Try linearity
    extract_res = extract_linear(f, x)
    if extract_res:
        # b * x + c = a  ==>  x = (a - c) / b
        b, c = extract_res
        if ctx.is_nonzero(b):
            return normalize((a - c) / b, ctx)

def extract_linear(e: Expr, x: str) -> Optional[Tuple[Expr, Expr]]:
    """Attempt to write e in the form a * x + b.
    
    If this is possible, return the pair (a, b). Otherwise return None.
    The results should be normalized before use.
    
    """
    if not e.contains_var(x):
        return Const(0), e
    elif expr.is_var(e):
        assert e.name == x
        return Const(1), Const(0)
    elif e.is_plus():
        res1 = extract_linear(e.args[0], x)
        res2 = extract_linear(e.args[1], x)
        if res1 and res2:
            return res1[0] + res2[0], res1[1] + res2[1]
    elif expr.is_uminus(e):
        res = extract_linear(e.args[0], x)
        if res:
            return -res[0], -res[1]
    elif e.is_minus():
        res1 = extract_linear(e.args[0], x)
        res2 = extract_linear(e.args[1], x)
        if res1 and res2:
            return res1[0] - res2[0], res1[1] - res2[1]
    elif e.is_times():
        u, v = e.args
        if not u.contains_var(x):
            res = extract_linear(v, x)
            if res:
                return u * res[0], u * res[1]
        elif not v.contains_var(x):
            res = extract_linear(u, x)
            if res:
                return v * res[0], v * res[1]
    elif e.is_divides():
        u, v = e.args
        if not v.contains_var(x):
            res = extract_linear(u, x)
            if res:
                return res[0] / v, res[1] / v

def solve_for_term(eq: Expr, t: Expr, ctx: Context) -> Optional[Expr]:
    """A more general solving procedure for term t.
    
    Given equation of the form f = g, where both f and g may contain t.
    Try to derive an equation of the form t = t' from f = g.
    
    """
    if not eq.is_equals():
        raise AssertionError("solve_for_term: input should be an equation.")

    # Take variable name that have not appeared
    var_name = "_v"
    var = Var(var_name)

    # Replace all appearances of t in equation by var
    eq = eq.replace(t, var)

    # Now consider some simple cases
    if not eq.rhs.contains_var(var_name):
        return solve_equation(eq.lhs, eq.rhs, var_name, ctx)
    
    if not eq.lhs.contains_var(var_name):
        return solve_equation(eq.rhs, eq.lhs, var_name, ctx)
    
    # Finally, try transforming the equation to f = 0
    res = solve_equation(eq.lhs - eq.rhs, Const(0), var_name, ctx)
    if res:
        if res.contains_var(var_name):
            raise AssertionError("solve_equation returns %s" % res)
        else:
            return res
