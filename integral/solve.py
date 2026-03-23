"""Functions for solving equations"""

import sympy
from fractions import Fraction
from typing import Optional, Tuple, List
from integral import expr
from integral.expr import Expr, POS_INF, NEG_INF, Const, Var, Op, Fun
from integral import sympywrapper
from integral.poly import normalize
from integral.context import Context


def _is_pure_z_squared(e: Expr, x: str) -> bool:
    """Check if e is exactly z^2 (no coefficient)."""
    return (
        expr.is_power(e)
        and expr.is_var(e.args[0])
        and e.args[0].name == x
        and e.args[1] == Const(2)
    )


def _solve_quadratic(f: Expr, a: Expr, x: str, ctx: Context) -> List[Expr]:
    """Solve quadratic equations: f(x) = a, where f is ax^2+bx+c.

    Parsing invariants:
      x^2 + x + c   →  Op(+, (Op(+, Op(^,x,2), x)), c)
      x^2 - x + c   →  Op(-, (Op(-, Op(^,x,2), x)), c)
      x^2 + c       →  Op(+, Op(^,x,2), c)
      x^2 - c       →  Op(-, Op(^,x,2), c)
      a*x^2 + b*x + c → Op(+, (Op(+, Op(*,a,Op(^,x,2)), Op(*,b,x))), c)
      a*x^2 - b*x + c → Op(-, (Op(-, Op(*,a,Op(^,x,2)), Op(*,b,x))), c)
      -(x^2) + ...    → unary minus node
    """
    def _coeff_of(term, var_name):
        """Return (ax2, bx) where term = ax2*x^2 + bx*x.
        Coefficients as int/Fraction, or (None, None) if not parseable."""
        # Op(*, c, x^2) or Op(^, x, 2)
        if _is_pure_z_squared(term, var_name):
            return (1, 0)
        if term.is_times() and _is_pure_z_squared(term.args[1], var_name):
            coeff = term.args[0]
            if isinstance(coeff, Const):
                return (coeff.val, 0)
            return (None, None)
        # Op(*, c, x) or Var(x)
        if expr.is_var(term) and term.name == var_name:
            return (0, 1)
        if term.is_times() and len(term.args) == 2:
            coeff = term.args[0]
            var_part = term.args[1]
            if expr.is_var(var_part) and var_part.name == var_name:
                if isinstance(coeff, Const):
                    return (0, coeff.val)
        # Unary minus of any of the above
        if expr.is_uminus(term):
            ax2, bx = _coeff_of(term.args[0], var_name)
            if ax2 is not None:
                return (-ax2, -bx)
        return (None, None)

    A_num, B_num, C_num = None, None, None

    # Pure z^2 = a
    if _is_pure_z_squared(f, x):
        A_num, B_num, C_num = 1, 0, 0

    elif getattr(f, 'op', None) in ("+", "-") and len(f.args) == 2:
        inner, C_raw = f.args
        if isinstance(C_raw, Const):
            C_num = -C_raw.val if f.op == "-" else C_raw.val
        else:
            return []
        inner_op = getattr(inner, 'op', None)

        if inner_op in ("+", "-") and len(inner.args) == 2:
            # inner = left +/- right
            left, right = inner.args
            sign = -1 if inner_op == "-" else 1
            ax2_l, bx_l = _coeff_of(left, x)
            ax2_r, bx_r = _coeff_of(right, x)
            if ax2_l is None or ax2_r is None or bx_l is None or bx_r is None:
                return []
            A_num, B_num = ax2_l, sign * bx_r
        elif expr.is_uminus(inner):
            # -(...) + C, e.g. -(x^2) + 4*x
            ax2, bx = _coeff_of(inner, x)
            if ax2 is None:
                return []
            A_num, B_num = ax2, bx
        else:
            # inner = term with no x (e.g. Op(^, x, 2))
            ax2, bx = _coeff_of(inner, x)
            if ax2 is None:
                return []
            A_num, B_num = ax2, bx

    if A_num is None:
        return []

    # Degenerate case: A=0 → linear equation bx + C_adj = 0
    if A_num == 0:
        return []

    R_num = normalize(a, ctx)
    if not isinstance(R_num, Const):
        return []
    R_num = R_num.val

    C_adj = C_num - R_num
    disc = B_num * B_num - 4 * A_num * C_adj
    half = Op("/", Const(1), Const(2))
    neg_B = Op("-", Const(B_num))
    two_A = Op("*", Const(2), Const(A_num))
    denom = normalize(two_A, ctx)
    sqrt_disc = Op("^", Const(disc), half)

    if disc < 0:
        abs_disc = -disc
        sqrt_abs = Op("^", Const(abs_disc), half)
        i_sqrt_abs = Op("*", sqrt_abs, Fun("i"))
        num_plus = Op("+", neg_B, i_sqrt_abs)
        num_minus = Op("+", neg_B, Op("-", i_sqrt_abs))
    else:
        num_plus = Op("+", neg_B, sqrt_disc)
        num_minus = Op("+", neg_B, Op("-", sqrt_disc))

    sol1 = normalize(Op("/", num_plus, denom), ctx)
    sol2 = normalize(Op("/", num_minus, denom), ctx)

    if _are_equal(sol1, sol2, ctx):
        return [sol1]
    return [sol1, sol2]


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

    a = normalize(a, ctx)

    # ========== Quadratic equations: ax^2+bx+c = a ==========
    # NOTE: placed BEFORE addition handler to avoid decomposing quadratic into pieces
    # that can't be solved rec.
    quad_res = _solve_quadratic(f, a, x, ctx)
    if quad_res:
        return quad_res

    # ========== Base case: variable ==========
    if expr.is_var(f):
        if f.name == x:
            return [a]
        else:
            return []

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
            # u^v = a  =>  u = a^(1/v). 当 n 为正整数时返回全部 n 个复根，供留数/极点使用。
            n_val = None
            if expr.is_const(v):
                try:
                    n_val = int(v.val)
                except (TypeError, ValueError):
                    pass
            if n_val is not None and n_val > 0:
                if a == Const(0):
                    return solve_equation(u, Const(0), x, ctx)
                base = normalize(a ^ (Const(1) / v), ctx)
                one_sol_list = solve_equation(u, base, x, ctx)
                if len(one_sol_list) == 1:
                    sol_base = one_sol_list[0]
                    # 全部 n 个复根: sol_base * exp(2*pi*i*k/n), k=0..n-1
                    roots = []
                    for k in range(n_val):
                        if k == 0:
                            root_k = sol_base
                        else:
                            two_pi_k_n = Op("/", Op("*", Op("*", Const(2), Fun("pi")), Const(k)), Const(n_val))
                            angle = Op("*", two_pi_k_n, Fun("i"))
                            exp_factor = Fun("exp", angle)
                            root_k = Op("*", sol_base, exp_factor)
                            root_k = normalize(root_k, ctx)
                        if not any(_are_equal(root_k, r, ctx) for r in roots):
                            roots.append(root_k)
                    if roots:
                        return roots
                # 若 principal 解非单解或生成根失败，则回退到原逻辑（偶次时 ± 两个解）
                if n_val == 2:
                    if len(one_sol_list) == 1:
                        sol = one_sol_list[0]
                        neg_sol = _create_negative(sol, ctx)
                        if not _are_equal(sol, neg_sol, ctx):
                            return [sol, neg_sol]
                    return one_sol_list
                return one_sol_list if one_sol_list else solve_equation(u, normalize(a ^ (Const(1) / v), ctx), x, ctx)
            if expr.is_const(v) and v.val % 2 == 0 and v.val > 0:
                one_sol_list = solve_equation(u, normalize(a ^ (Const(1) / v), ctx), x, ctx)
                if len(one_sol_list) == 1:
                    sol = one_sol_list[0]
                    neg_sol = _create_negative(sol, ctx)
                    if not _are_equal(sol, neg_sol, ctx):
                        return [sol, neg_sol]
                return one_sol_list
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

    # ========== Polynomial factorization ==========
    if a == Const(0):
        try:
            factored = _factor_polynomial(f, x, ctx)
            if factored and len(factored) > 1:
                # 对每个因式递归求根，去重后返回
                solutions = []
                for factor in factored:
                    if factor.contains_var(x):
                        sols = solve_equation(factor, Const(0), x, ctx)
                        for sol in sols:
                            if not any(_are_equal(sol, s, ctx) for s in solutions):
                                solutions.append(sol)
                if solutions:
                    return solutions
        except Exception:
            pass

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
    except Exception:
        return False


def _factor_polynomial(f: Expr, x: str, ctx: Context) -> Optional[List[Expr]]:
    """尝试使用 sympy 对多项式 f(x) 进行因式分解。

    Returns:
        None  - 不是可分解的多项式
        []    - 空多项式
        [f1, f2, ...] - 分解后的因式列表（各因式均含变量 x）
    """
    try:
        sym_f = sympywrapper.convert_to_sympy(f)
        x_sym = sympy.symbols(x)
        factored_sym = sympy.factor(sym_f, x_sym)
        # 如果不可分解，factor 返回原式（与原式相等）
        if factored_sym == sym_f:
            return None
        # 展开并收集各因式
        factors = []
        if factored_sym.is_Mul:
            for factor_sym in factored_sym.args:
                factor_expr = sympywrapper.convert_from_sympy(factor_sym)
                if factor_expr.contains_var(x):
                    factors.append(factor_expr)
        elif factored_sym.is_Pow:
            base_sym, exp_sym = factored_sym.args
            base_expr = sympywrapper.convert_from_sympy(base_sym)
            if base_expr.contains_var(x):
                if isinstance(exp_sym, sympy.Integer) and int(exp_sym) > 0:
                    # 展开幂: x^4 -> x*x*x*x
                    for _ in range(int(exp_sym)):
                        factors.append(base_expr)
                else:
                    factors.append(sympywrapper.convert_from_sympy(factored_sym))
        else:
            if factored_sym.contains(x_sym):
                factors.append(sympywrapper.convert_from_sympy(factored_sym))
        return factors if factors else None
    except Exception:
        return None
