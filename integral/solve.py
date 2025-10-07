"""Functions for solving equations"""

from typing import Optional, Tuple, List
from integral import expr
from integral.expr import Expr, POS_INF, NEG_INF, Const, Var, Op, Fun
from integral.poly import normalize
from integral.context import Context

def solve_equation(f: Expr, a: Expr, x: str, ctx: Context) -> List[Expr]:
    """Solve the equation f(x) = a for variable x, returning ALL solutions.
    
    Breaking change: Returns List[Expr] instead of Optional[Expr].
    
    Returns:
        [] - no solution found
        [sol1, sol2, ...] - all candidate solutions
    
    Examples:
        solve_equation(x, 5, 'x', ctx) => [5]
        solve_equation(x^2 + 1, 0, 'x', ctx) => [i, -i]
        solve_equation((x-1)*(x-2), 0, 'x', ctx) => [1, 2]
    """
    # Backward compatibility: handle Var object as x parameter
    if isinstance(x, Var):
        x = x.name
    elif isinstance(x, Expr) and expr.is_var(x):
        x = x.name
    
    # Note: We don't normalize f here as it might change structure (e.g., division to multiplication)
    # We only normalize a (the right-hand side) to simplify constants
    a = normalize(a, ctx)
    
    # ========== Base case: variable ==========
    if expr.is_var(f):
        if f.name == x:
            return [a]
        else:
            return []  # Unrelated variable
    
    # ========== Product equals zero: (u*v*w = 0) ==========
    if a == Const(0) and f.is_times():
        solutions = []
        for factor in f.args:
            if factor.contains_var(x):
                factor_sols = solve_equation(factor, Const(0), x, ctx)
                for sol in factor_sols:
                    # Avoid duplicates
                    if not any(_are_equal(sol, s, ctx) for s in solutions):
                        solutions.append(sol)
        return solutions
    
    # ========== Addition: u + v = a ==========
    if expr.is_plus(f):
        u, v = f.args
        if not u.contains_var(x):
            # u + v = a  =>  v = a - u
            return solve_equation(v, normalize(a - u, ctx), x, ctx)
        if not v.contains_var(x):
            # u + v = a  =>  u = a - v
            return solve_equation(u, normalize(a - v, ctx), x, ctx)
    
    # ========== Unary minus: -u = a ==========
    if expr.is_uminus(f):
        u, = f.args
        return solve_equation(u, normalize(-a, ctx), x, ctx)
    
    # ========== Subtraction: u - v = a ==========
    if f.is_minus():
        u, v = f.args
        if not u.contains_var(x):
            # u - v = a  =>  v = u - a
            return solve_equation(v, normalize(u - a, ctx), x, ctx)
        if not v.contains_var(x):
            # u - v = a  =>  u = v + a
            return solve_equation(u, normalize(v + a, ctx), x, ctx)
    
    # ========== Multiplication: u * v = a ==========
    if f.is_times():
        u, v = f.args
        if not u.contains_var(x) and ctx.is_nonzero(u):
            # u * v = a  =>  v = a / u
            return solve_equation(v, normalize(a / u, ctx), x, ctx)
        if not v.contains_var(x) and ctx.is_nonzero(v):
            # u * v = a  =>  u = a / v
            return solve_equation(u, normalize(a / v, ctx), x, ctx)
        # If neither is constant, can't solve directly
        return []
    
    # ========== Division: u / v = a ==========
    if f.is_divides():
        u, v = f.args
        if not u.contains_var(x):
            # u / v = a  =>  v = u / a
            if a == Const(0):
                return []  # Division by zero
            rhs = u / a
            if u.is_constant() and a in (POS_INF, NEG_INF):
                rhs = Const(0)
            return solve_equation(v, normalize(rhs, ctx), x, ctx)
        if not v.contains_var(x):
            # u / v = a  =>  u = v * a
            return solve_equation(u, normalize(v * a, ctx), x, ctx)
    
    # ========== Power: u^n = a ==========
    if f.is_power():
        u, v = f.args
        if not v.contains_var(x):
            # u^v = a  =>  u = a^(1/v)
            # Special handling for even powers: return both ±solutions
            if expr.is_const(v) and v.val % 2 == 0 and v.val > 0:
                # Even power: get one solution and add its negative
                one_sol_list = solve_equation(u, normalize(a ^ (Const(1) / v), ctx), x, ctx)
                if len(one_sol_list) == 1:
                    sol = one_sol_list[0]
                    neg_sol = _create_negative(sol, ctx)
                    # Check if they're different
                    if not _are_equal(sol, neg_sol, ctx):
                        return [sol, neg_sol]
                return one_sol_list
            else:
                # Odd power or non-integer: only one solution
                return solve_equation(u, normalize(a ^ (Const(1) / v), ctx), x, ctx)
    
    # ========== Functions ==========
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
            return solve_equation(f.args[0], normalize(a ^ Const(2), ctx), x, ctx)
    
    # ========== Linear extraction: b*x + c = a ==========
    extract_res = extract_linear(f, x)
    if extract_res:
        b, c = extract_res
        b = normalize(b, ctx)
        c = normalize(c, ctx)
        if ctx.is_nonzero(b):
            return [normalize((a - c) / b, ctx)]
        else:
            # b = 0: equation becomes c = a
            if _are_equal(c, a, ctx):
                # Identity: infinite solutions (can't represent, return empty)
                return []
            else:
                # Contradiction: no solution
                return []
    
    # ========== Default: unable to solve ==========
    return []


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
        solutions = solve_equation(eq.lhs, eq.rhs, var_name, ctx)
        return solutions[0] if solutions else None
    
    if not eq.lhs.contains_var(var_name):
        solutions = solve_equation(eq.rhs, eq.lhs, var_name, ctx)
        return solutions[0] if solutions else None
    
    # Finally, try transforming the equation to f = 0
    solutions = solve_equation(eq.lhs - eq.rhs, Const(0), var_name, ctx)
    if solutions:
        res = solutions[0]
        if res.contains_var(var_name):
            raise AssertionError("solve_equation returns %s" % res)
        else:
            return res


def _create_negative(e: Expr, ctx: Context) -> Expr:
    """Create negative of an expression with special handling for complex numbers.
    
    For b*i, creates (-b)*i instead of -(b*i) to maintain canonical form.
    """
    # Special case: b*i => (-b)*i
    if isinstance(e, Op) and e.op == "*" and len(e.args) == 2:
        if isinstance(e.args[1], Fun) and e.args[1].func_name == "i":
            b = e.args[0]
            if isinstance(b, Const):
                return Op("*", Const(-b.val), Fun("i"))
            else:
                return Op("*", Op("-", b), Fun("i"))
    
    # General case: -e
    return normalize(Op("-", e), ctx)


def _are_equal(e1: Expr, e2: Expr, ctx: Context) -> bool:
    """Check if two expressions are equal after normalization."""
    try:
        n1 = normalize(e1, ctx)
        n2 = normalize(e2, ctx)
        return n1 == n2
    except:
        return False
