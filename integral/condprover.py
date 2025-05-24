"""Module for reasoning about conditions."""

from copy import copy

from integral import expr
from integral.expr import Expr, eval_expr, match, expr_to_pattern, Op, Const, Var, Fun
from integral.conditions import Conditions
from integral.context import Context, Identity
from integral.parser import parse_expr


def subject_of(cond: Expr) -> Expr:
    """Return the subject of a condition.

    The following rules are used to determine subject:

    - For equality and comparisons, left side of the operator
    - For predicates with a single argument (such as isEven), its only
      argument.
    
    """
    if expr.is_equals(cond) or expr.is_not_equals(cond):
        return cond.args[0]
    if expr.is_greater(cond) or expr.is_greater_eq(cond):
        return cond.args[0]
    if expr.is_less(cond) or expr.is_less_eq(cond):
        return cond.args[0]
    if expr.is_fun(cond):
        return cond.args[0]
    raise NotImplementedError(f"subject_of: {cond}")

# Tolerance for floating-point rounding errors
tol = 1e-15

# Comparison of floating-point numbers up to rounding error
def approx_equal(a: Expr, b: Expr) -> bool:
    try:
        a_val = float(eval_expr(a))
        b_val = float(eval_expr(b))
        return abs(a_val - b_val) < tol
    except:
        return False

def approx_not_equal(a: Expr, b: Expr) -> bool:
    try:
        a_val = float(eval_expr(a))
        b_val = float(eval_expr(b))
        return abs(a_val - b_val) > tol
    except:
        return True

def approx_greater(a: Expr, b: Expr) -> bool:
    try:
        a_val = float(eval_expr(a))
        b_val = float(eval_expr(b))
        return a_val - b_val > tol
    except:
        return False

def approx_greater_eq(a: Expr, b: Expr) -> bool:
    try:
        a_val = float(eval_expr(a))
        b_val = float(eval_expr(b))
        return a_val - b_val > -tol
    except:
        return False

def approx_less(a: Expr, b: Expr) -> bool:
    try:
        a_val = float(eval_expr(a))
        b_val = float(eval_expr(b))
        return b_val - a_val > tol
    except:
        return False

def approx_less_eq(a: Expr, b: Expr) -> bool:
    try:
        a_val = float(eval_expr(a))
        b_val = float(eval_expr(b))
        return b_val - a_val > -tol
    except:
        return False

def approx_integer(a: Expr) -> bool:
    try:
        a_val = float(eval_expr(a))
        return abs(round(a_val) - a_val) < tol
    except:
        return False

def approx_even(a: Expr) -> bool:
    try:
        if approx_integer(a):
            a_val = float(eval_expr(a))
            return round(a_val) % 2 == 0
        return False
    except:
        return False
    
def approx_real(a: Expr) -> bool:
    try:
        a_val = complex(eval_expr(a))
        return abs(a_val.imag) < tol
    except:
        return False

def approx_not_real(a: Expr) -> bool:
    try:
        a_val = complex(eval_expr(a))
        return abs(a_val.imag) > tol
    except:
        return False

def update_inst(k: str, v: Expr, inst: dict[str, Expr]) -> dict[str, Expr]:
    """Update instantiation without changing the original."""
    res = copy(inst)
    res[k] = v
    return res

def check_cond(cond: Expr, all_conds: dict[Expr, list[Expr]],
               inst: dict[str, Expr]) -> list[dict[str, Expr]]:
    """Determine whether cond is implied by the existing set of
    conditions. There may be uninstantiated symbols in cond. The
    list of valid instantiations are returned.
    
    The following checks are performed:

    - If subject of cond is a constant, and the right side is also constant,
      compare using numerical calculations.

    - If subject of cond appears in all_conds, try to use the conditions
      available to verify cond.

    - Perform pattern matching.

    Parameters
    ----------
    cond: Expr
        the condition to be checked
    all_conds: dict[Expr, list[Expr]]
        mapping from subject to list of conditions on the subject
    inst: dict[str, Expr]
        current instantiation

    Returns
    -------
    list[dict[str, Expr]]
        list of new instantiations 

    """
    # Get subject of the cond
    x = subject_of(cond)

    # Trivial case: cond already appears as a fact.
    if x in all_conds and cond in all_conds[x]:
        return [inst]

    # If subject of cond is a constant, evaluate using numerical
    # calculations (TODO: this is not guaranteed to be correct).
    if x.is_constant():
        if expr.is_equals(cond) and cond.args[1].is_constant():
            if approx_equal(x, cond.args[1]):
                return [inst]
        elif expr.is_not_equals(cond) and cond.args[1].is_constant():
            if approx_not_equal(x, cond.args[1]):
                return [inst]
        elif expr.is_greater(cond) and cond.args[1].is_constant():
            if approx_greater(x, cond.args[1]):
                return [inst]
        elif expr.is_greater_eq(cond) and cond.args[1].is_constant():
            if approx_greater_eq(x, cond.args[1]):
                return [inst]
        elif expr.is_less(cond) and cond.args[1].is_constant():
            if approx_less(x, cond.args[1]):
                return [inst]
        elif expr.is_less_eq(cond) and cond.args[1].is_constant():
            if approx_less_eq(x, cond.args[1]):
                return [inst]
        elif expr.is_fun(cond, "isInt"):
            if approx_integer(x):
                return [inst]
        elif expr.is_fun(cond, "isEven"):
            if approx_even(x):
                return [inst]
        elif expr.is_fun(cond, "isReal"):
            if approx_real(x):
                return [inst]
        elif expr.is_fun(cond, "isComplex"):
            return [inst]

    # If the goal is of form x ?= c, where c is a constant, try to
    # apply transitivity with facts in all_conds.
    if expr.is_compare(cond) and x in all_conds and cond.args[1].is_constant():
        for fact in all_conds[x]:
            if not (expr.is_compare(fact) and fact.args[1].is_constant()):
                continue
            if expr.is_greater_eq(cond):
                # x >= b --> b >= a --> x >= a
                if expr.is_greater(fact) or expr.is_greater_eq(fact):
                    if approx_greater_eq(fact.args[1], cond.args[1]):
                        return [inst]
            if expr.is_greater(cond):
                # x >= b --> b > a --> x > a
                if expr.is_greater_eq(fact) and approx_greater(fact.args[1], cond.args[1]):
                    return [inst]
                # x > b --> b >= a --> x > a
                if expr.is_greater(fact) and approx_greater_eq(fact.args[1], cond.args[1]):
                    return [inst]
            if expr.is_less_eq(cond):
                # x <= b --> b <= a --> x <= a
                if expr.is_less(fact) or expr.is_less_eq(fact):
                    if approx_less_eq(fact.args[1], cond.args[1]):
                        return [inst]
            if expr.is_less(cond):
                # x <= b --> b < a --> x < a
                if expr.is_less_eq(fact) and approx_less(fact.args[1], cond.args[1]):
                    return [inst]
                # x < b --> b <= a --> x < a
                if expr.is_less(fact) and approx_less_eq(fact.args[1], cond.args[1]):
                    return [inst]
            if expr.is_equals(cond):
                if expr.is_equals(fact) and approx_equal(fact.args[1], cond.args[1]):
                    return [inst]
            if expr.is_not_equals(cond):
                if expr.is_not_equals(fact) and approx_equal(fact.args[1], cond.args[1]):
                    return [inst]
                # x < a --> a <= b --> x != b
                if expr.is_less(fact) and approx_less_eq(fact.args[1], cond.args[1]):
                    return [inst]
                # x <= a --> a < b --> x != b
                if expr.is_less_eq(fact) and approx_less(fact.args[1], cond.args[1]):
                    return [inst]
                # x > a --> a >= b --> x != b
                if expr.is_greater(fact) and approx_greater_eq(fact.args[1], cond.args[1]):
                    return [inst]
                # x >= a --> a > b --> x != b
                if expr.is_greater_eq(fact) and approx_greater(fact.args[1], cond.args[1]):
                    return [inst]
                # x = a --> a != b --> x != b
                if expr.is_equals(fact) and approx_not_equal(fact.args[1], cond.args[1]):
                    return [inst]
        
    # If the other side of cond is a symbol, update with appropriate
    # instantiation.
    if expr.is_compare(cond) and x in all_conds and expr.is_symbol(cond.args[1]):
        symb = cond.args[1].name
        res = []
        for fact in all_conds[x]:
            if not expr.is_compare(fact):
                continue
            if expr.is_greater_eq(cond):
                if expr.is_greater_eq(fact) or expr.is_greater(fact):
                    res.append(update_inst(symb, fact.args[1], inst))
            if expr.is_greater(cond):
                if expr.is_greater(fact):
                    res.append(update_inst(symb, fact.args[1], inst))
            if expr.is_less_eq(cond):
                if expr.is_less_eq(fact) or expr.is_less(fact):
                    res.append(update_inst(symb, fact.args[1], inst))
            if expr.is_less(cond):
                if expr.is_less(fact):
                    res.append(update_inst(symb, fact.args[1], inst))
            if expr.is_equals(cond):
                if expr.is_equals(fact):
                    res.append(update_inst(symb, fact.args[1], inst))
            if expr.is_not_equals(cond):
                if expr.is_not_equals(fact):
                    res.append(update_inst(symb, fact.args[1], inst))
        return res

    # Not found
    return list()

def add_condition(all_conds: dict[Expr, list[Expr]], e: Expr, cond: Expr):
    """Helper function: add condition `cond` to `all_conds` only if
    it is not already implied by the existing conditions.
    
    """
    if e not in all_conds:
        all_conds[e] = list()
    if len(check_cond(cond, all_conds, dict())) == 0:
        all_conds[e].append(cond)

def init_all_conds(conds: Conditions) -> dict[Expr, list[Expr]]:
    """Initialize mapping from subject to list of facts from a given
    condition object.
    
    """
    all_conds: dict[Expr, list[Expr]] = dict()
    
    # Rewrite all absolute value conditions
    for cond in conds.data:
        if not (expr.is_compare(cond) or expr.is_fun(cond)):
            continue
        x = subject_of(cond)
        add_condition(all_conds, x, cond)
            
        # Handle absolute value conditions
        if expr.is_fun(x, 'abs') and expr.is_less(cond):
            # abs(x) < c  -->  -c < x < c
            add_condition(all_conds, x.args[0], Op("<", x.args[0], cond.args[1]))
            add_condition(all_conds, x.args[0], Op(">", x.args[0], -cond.args[1]))
        if expr.is_fun(x, 'abs') and expr.is_less_eq(cond):
            # abs(x) <= c  -->  -c <= x <= c
            add_condition(all_conds, x.args[0], Op("<=", x.args[0], cond.args[1]))
            add_condition(all_conds, x.args[0], Op(">=", x.args[0], -cond.args[1]))

    # add simple condition transition
    for k in all_conds:
        for x in all_conds[k]:
            if expr.is_less(x):
                if x.args[1] in all_conds:
                    for y in all_conds[x.args[1]]:
                        # x: k < b
                        # y: b < c or b <= c or b = c
                        # x and y ==> k < c
                        if expr.is_less(y) or expr.is_less_eq(y) or expr.is_equals(y):
                            add_condition(all_conds, k, Op('<', k, y.args[1]))
    return all_conds

def saturate_expr(e: Expr, ineq: Identity, all_conds: dict[Expr, list[Expr]]):
    """Use the rule `ineq` to saturate facts about `e`. New facts are
    added to `all_conds`.
    
    """
    pat = subject_of(ineq.expr)
    inst = match(e, pat)
    if inst is not None:
        # Check conditions of the inequality, instantiating schematic
        # variables in the inequality in the process.
        old_list = [inst]
        for cond in ineq.conds.data:
            new_list = []
            for inst in old_list:
                res = check_cond(cond.inst_pat(inst), all_conds, inst)
                new_list.extend(res)
            old_list = new_list
        # For each instantiation, apply the identity
        for mapping in old_list:
            res = ineq.expr.inst_pat(mapping)
            add_condition(all_conds, e, res)
    return

def saturate_once(e: Expr, ineqs: list[Identity], all_conds: dict[Expr, list[Expr]]):
    """Perform one round of saturation. New facts are added
    onto `all_conds`.
    
    """
    all_subs = e.find_all_subexpr()
    for sube, _ in all_subs:
        for ineq in ineqs:
            saturate_expr(sube, ineq, all_conds)

def all_conds_size(all_conds: dict[Expr, list[Expr]]) -> int:
    """Return number of facts in all_conds."""
    res = 0
    for _, conds in all_conds.items():
        res += len(conds)
    return res

def saturate(e: Expr, ineqs: list[Identity], all_conds: dict[Expr, list[Expr]], *,
             round_limit: int = 5, size_limit: int = 200):
    """Saturate up to given number of rounds and size limits. New facts
    are added onto `all_conds`.
    
    If number of rounds and size limits have been reached without
    saturation, warning is printed to alert possible problems.
    
    """
    i = 0
    while True:
        prev_size = all_conds_size(all_conds)
        saturate_once(e, ineqs, all_conds)
        i += 1
        next_size = all_conds_size(all_conds)
        if prev_size == next_size:
            # No new facts are added
            return
        if next_size > size_limit:
            print(f"Warning: size limit reached during saturation, size = {next_size}, limit = {size_limit}")
            print_all_conds(all_conds)
            raise AssertionError
        if i > round_limit:
            print(f"Warning: round limit reached during saturation")
            return

def print_all_conds(all_conds: dict[Expr, list[Expr]]):
    """Print all conditions (for debugging)."""
    for x, conds in all_conds.items():
        print("%s: %s\n" % (x, ', '.join(str(cond) for cond in conds)))

def get_standard_inequalities() -> list[Identity]:
    """List of standard inequalities."""
    data = [
        # Addition
        (["c > 0"], "a + c > a"),
        (["a > 0"], "a + c > c"),
        (["a > b"], "a + c > b + c"),
        (["a > b"], "c + a > c + b"),
        (["a < b"], "a + c < b + c"),
        (["a < b"], "c + a < c + b"),
        (["a != b"], "a + c != b + c"),
        (["a != b"], "c + a != c + b"),
        (["a >= b"], "a + c >= b + c"),
        (["a >= b"], "c + a >= c + b"),
        (["a <= b"], "a + c <= b + c"),
        (["a <= b"], "c + a <= c + b"),
        (["a != -b"], "a + b != 0"),
        (["a != b"], "a - b != 0"),
        (["a < 0", "isReal(b)"], "a + b * i != 0"),
        (["a < 0", "isReal(b)"], "a - b * i != 0"),
        (["a >= b", "c > d"], "a + c > b + d"),
        (["a > b", "c >= d"], "a + c > b + d"),
        (["a <= b", "c < d"], "a + c < b + d"),
        (["a < b", "c <= d"], "a + c < b + d"),
        (["a >= b", "c >= d"], "a + c >= b + d"),
        (["a <= b", "c <= d"], "a + c <= b + d"),

        # Unary minus
        (["x > a"], "-x < -a"),
        (["x < a"], "-x > -a"),
        (["x != a"], "-x != -a"),
        (["x >= a"], "-x <= -a"),
        (["x <= a"], "-x >= -a"),

        # Subtraction
        (["a > b"], "c - a < c - b"),
        (["a < b"], "c - a > c - b"),
        (["a > b"], "a - c > b - c"),
        (["a < b"], "a - c < b - c"),
        (["a >= b"], "c - a <= c - b"),
        (["a <= b"], "c - a >= c - b"),
        (["a >= b"], "a - c >= b - c"),
        (["a <= b"], "a - c <= b - c"),
        (["a > b", "c <= d"], "a - c > b - d"),
        (["a >= b", "c < d"], "a - c > b - d"),
        (["a < b", "c >= d"], "a - c < b - d"),
        (["a <= b", "c > d"], "a - c < b - d"),
        (["a >= b", "c <= d"], "a - c >= b - d"),
        (["a <= b", "c >= d"], "a - c <= b - d"),
        (["a < b"], "a - b < 0"),
        (["a < b"], "b - a > 0"),
        (["a > -b"], "a + b > 0"),
        (["a > -b"], "b + a > 0"),

        # Multiplication (simple)
        (["a != 0", "b != 0"], "a * b != 0"),
        (["a > 0", "b > 0"], "a * b > 0"),
        (["a < 0", "b > 0"], "a * b < 0"),
        (["a > 0", "b < 0"], "a * b < 0"),
        (["a < 0", "b < 0"], "a * b > 0"),

        # Multiplication (one side is constant)
        (["a > b", "c > 0"], "c * a > c * b"),
        (["a > b", "c > 0"], "a * c > b * c"),
        (["a < b", "c > 0"], "c * a < c * b"),
        (["a < b", "c > 0"], "a * c < b * c"),
        (["a >= b", "c >= 0"], "c * a >= c * b"),
        (["a >= b", "c >= 0"], "a * c >= b * c"),
        (["a <= b", "c >= 0"], "c * a <= c * b"),
        (["a <= b", "c >= 0"], "a * c <= b * c"),

        # Multiplication (left side > 0)
        (["a > b", "c >= d", "b > 0"], "a * c > b * d"),
        (["a >= b", "c > d", "b > 0"], "a * c > b * d"),
        (["a < b", "c < d", "a > 0"], "a * c < b * d"),
        (["a >= b", "c >= d", "b >= 0"], "a * c >= b * d"),
        (["a <= b", "c < d", "a > 0"], "a * c <= b * d"),
        (["a < b", "c <= d", "a > 0"], "a * c <= b * d"),

        # Multiplication (right side > 0)
        (["a > b", "c >= d", "b > 0"], "c * a > d * b"),
        (["a >= b", "c > d", "b > 0"], "c * a > d * b"),
        (["a <= b", "c < d", "a > 0"], "c * a < d * b"),
        (["a < b", "c <= d", "a > 0"], "c * a < d * b"),
        (["a >= b", "c >= d", "b >= 0"], "c * a >= d * b"),
        (["a <= b", "c <= d", "a > 0"], "c * a <= d * b"),

        # Division
        (["a > 0", "b > 0"], "a / b > 0"),
        (["a > 0", "b < 0"], "a / b < 0"),
        (["a < 0", "b > 0"], "a / b < 0"),
        (["a < 0", "b < 0"], "a / b > 0"),
        (["a > b", "c > 0"], "a / c > b / c"),
        (["a < b", "c > 0"], "a / c < b / c"),
        (["a > b", "c < 0"], "a / c < b / c"),
        (["a < b", "c < 0"], "a / c > b / c"),
        (["a >= b", "c > 0"], "a / c >= b / c"),
        (["a <= b", "c > 0"], "a / c <= b / c"),
        (["a >= b", "c < 0"], "a / c <= b / c"),
        (["a <= b", "c < 0"], "a / c >= b / c"),
        (["a > 0", "b != 0"], "a / b != 0"),
        (["a < 0", "b != 0"], "a / b != 0"),
        (["a != 0", "b != 0"], "a / b != 0"),
        (["x > 1"], "1 / x < 1"),
        (["x > 0"], "1 / x > 0"),
        (["x > -a", "x < a"], "x / a < 1"),
        (["x > -a", "x < a"], "x / a > -1"),

        # Square root
        ([], "sqrt(a) >= 0"),
        (["a > 0"], "sqrt(a) > 0"),
        (["a < 1", "a > 0"], "sqrt(a) < 1"),
        (["a > b", "b >= 0"], "sqrt(a) > sqrt(b)"),
        (["a >= b", "b >= 0"], "sqrt(a) >= sqrt(b)"),
        (["a >= 0"], "sqrt(a ^ 2) = a"),

        # Power
        (["a != 0"], "a ^ 2 > 0"),
        ([], "a ^ 2 >= 0"),
        (["x > 0"], "x ^ y > 0"),
        (["x != 0"], "x ^ n != 0"),
        (["x > y", "y >= 0", "z > 0"], "x ^ z > y ^ z"),
        (["x < a", "x > -a"], "x ^ 2 < a ^ 2"),
        (["x > a", "a >= 0"], "x ^ 2 > a ^ 2"),
        (["x < -1"], "x ^ 2 > 1"),
        (["x <= a", "x >= -a"], "x ^ 2 <= a ^ 2"),
        (["x >= a", "a >= 0"], "x ^ 2 >= a ^ 2"),
        (["x != y"], "x ^ 2 - y ^ 2 != 0"),
        (["y != x"], "x ^ 2 - y ^ 2 != 0"),
        (["x != y"], "x ^ 4 - y ^ 4 != 0"),
        (["y != x"], "x ^ 4 - y ^ 4 != 0"),
        (["x > 0", "x < 1", "y > 1"], "x ^ y < 1"),
        (["x > 0", "x <= 1", "isInt(n)", "n > 0"], "x - x^n >= 0"),
        (["x > 0", "x < 1", "y > 1"], "x ^ y > 0"),

        # Log
        (["x >= 1"], "log(x) >= 0"),
        (["x > 1"], "log(x) > 0"),
        (["x <= 1", "x > 0"], "log(x) <= 0"),
        (["x < 1", "x > 0"], "log(x) < 0"),
        (["x != 1"], "log(x) != 0"),
        (["a > b", "b > 0"], "log(a) > log(b)"),

        # Absolute value
        (["x != 0"], "abs(x) > 0"),

        # Exponential
        ([], "exp(x) > 0"),
        (["x > 0"], "exp(x) > 1"),
        (["x < 0"], "exp(x) < 1"),
        (["x >= 0"], "exp(x) >= 1"),
        (["x <= 0"], "exp(x) <= 1"),

        # Trigonometric
        (["x > -pi / 2", "x < pi / 2"], "cos(x) > 0"),
        (["x > -pi / 2", "x < pi / 2"], "cos(x) <= 1"),
        (["x > pi / 2", "x < 3 * pi / 2"], "cos(x) < 0"),
        (["x >= -pi / 2", "x <= pi / 2"], "cos(x) >= 0"),
        (["x >= pi / 2", "x <= 3 * pi / 2"], "cos(x) <= 0"),
        (["x > 0", "x < 2 * pi"], "cos(x) < 1"),
        (["x > -pi / 2", "x < pi / 2"], "sin(x) > -1"),
        (["x > -pi / 2", "x < pi / 2"], "sin(x) < 1"),
        (["x > 0", "x < pi"], "sin(x) > 0"),
        (["x > -pi", "x < 0"], "sin(x) < 0"),
        (["x >= 0", "x <= pi"], "sin(x) >= 0"),
        (["x >= -pi", "x <= 0"], "sin(x) <= 0"),
        (["x > -pi / 2", "x < pi / 4"], "tan(x) < 1"),
        (["x > 0", "x < pi / 2"], "tan(x) > 0"),
        (["cos(x) != 0"], "sin(x) > -1"),
        (["cos(x) != 0"], "sin(x) < 1"),
        (["sin(x) != 0"], "cos(x) > -1"),
        (["sin(x) != 0"], "cos(x) < 1"),
        (["x > -pi / 2", "x < pi / 2"], "sec(x) >= 1"),
        (["x > pi / 4", " x < pi / 2"], "sec(x) < sqrt(2)"),

        # Inverse trigonometric
        (["x >= -1", "x <= 1"], "arcsin(x) >= -pi / 2"),
        (["x >= -1", "x <= 1"], "arcsin(x) <= pi / 2"),
        (["x > -1", "x < 1"], "arcsin(x) > -pi / 2"),
        (["x > -1", "x < 1"], "arcsin(x) < pi / 2"),
        (["x > 0", "x <= 1"], "arcsin(x) > 0"),
        (["x >= 0", "x <= 1"], "arcsin(x) >= 0"),
        (["x < 0", "x >= -1"], "arcsin(x) < 0"),
        (["x <= 0", "x >= -1"], "arcsin(x) <= 0"),
        (["x > 0", "x < 1 / 2"], "arcsin(x) < pi / 6"),
        (["x != -1"], "arcsin(x) != -pi/2"),
        (["x != 1"], "arcsin(x) != pi/2"),
        (["x >= -1", "x <= 1"], "arccos(x) >= 0"),
        (["x >= -1", "x <= 1"], "arccos(x) <= pi"),
        (["x > 0", "x <= 1"], "arccos(x) < pi / 2"),
        (["x >= 0", "x <= 1"], "arccos(x) <= pi / 2"),
        (["x < 0", "x >= -1"], "arccos(x) > pi / 2"),
        (["x <= 0", "x >= -1"], "arccos(x) >= pi / 2"),
        (["x > -1", "x < 1"], "arccos(x) > 0"),
        (["x > -1", "x < 1"], "arccos(x) < pi"),
        ([], "arctan(x) > -pi / 2"),
        ([], "arctan(x) < pi / 2"),
        (["x >= 0"], "arctan(x) >= 0"),
        (["x > 0"], "arctan(x) > 0"),
        (["x <= 0"], "arctan(x) <= 0"),
        (["x < 0"], "arctan(x) < 0"),
        (["x >= 1"], "arcsec(x) >= 0"),
        (["x >= 1"], "arcsec(x) < pi / 2"),
        (["x > 1"], "arcsec(x) > 0"),
        (["x <= -1"], "arcsec(x) <= pi"),
        (["x <= -1"], "arcsec(x) > pi / 2"),
        (["x < -1"], "arcsec(x) < pi"),
        (["x >= 1"], "arccsc(x) <= pi / 2"),
        (["x >= 1"], "arcsec(x) > 0"),
        (["x > 1"], "arcsec(x) < pi / 2"),
        (["x <= -1"], "arcsec(x) >= -pi / 2"),
        (["x <= -1"], "arcsec(x) < 0"),
        (["x < -1"], "arcsec(x) > -pi / 2"),
        (["x > 1"], "arctan(x) > pi / 4"),
        (["x > 0"], "arctan(x) < pi / 2"),
        (["x > 0"], "arctan(x) > 0"),
        ([], "arctan(x) > -pi/2"),
        ([], "arctan(x) < pi/2"),

        # Value comparison of trig functions in an interval
        (["x > 0", "x < pi/4"], "cos(x) > sin(x)"),

        # Hyperbolic
        ([], "cosh(x) > 0"),

        # Factorial
        ([], "factorial(x) >= 1"),

        # isInt
        (["isInt(a)", "isInt(b)"], "isInt(a + b)"),
        (["isInt(a)", "isInt(b)"], "isInt(a - b)"),
        (["isInt(a)", "isInt(b)"], "isInt(a * b)"),

        # inequality of !=
        (["a > 0"], "a != 0"),
        (["a < 0"], "a != 0"),
        (["a > b"], "a - b != 0"),
        (["a > b"], "b - a != 0"),
        (["x != t", "b != 0"], "x * b != t * b"),
        (["x != t", "b != 0"], "b * x != b * t"),
        (["x != pi/2", "x != -pi/2"], "cos(x) != 0"),
        (["x != 0"], "arctan(x) != 0"),
        (["x != 0"], "arctan(x) > 0"),
        (["x != 0"], "arctan(x) > -pi/2"),
        (["x != 0"], "arctan(x) < pi /2 "),
        (["x != 0"], "arctan(x) < 0"),
        (["x > 0", "x < pi/2"], "sin(x) > 0"),
        (["x < 0", "x > -pi/2"], "sin(x) < 0"),
        (["x > 0", "x < pi / 6"], "sin(x) < 1 / 2"),

        (["a >= b", "a != b"], "a > b"),
        (["a <= b", "a != b"], "a < b"),
        (["a = b", "a > c"], "b > c"),
        (["a > b", "b > c"], "a > c"),

        # Complex number rules
        (["isReal(a)"], "isReal(cos(a))"),
        (["isReal(a)"], "isReal(sin(a))"),
        (["isReal(a)"], "isReal(tan(a))"),
        (["isReal(a)"], "isReal(cot(a))"),
        (["isReal(a)"], "isReal(csc(a))"),
        (["isReal(a)"], "isReal(sec(a))"),

        (["isReal(a)"], "isReal(arcsin(a))"),
        (["isReal(a)"], "isReal(arccos(a))"),
        (["isReal(a)"], "isReal(arctan(a))"),
        (["isReal(a)"], "isReal(arcsec(a))"),
        (["isReal(a)"], "isReal(arccsc(a))"),
        (["isReal(a)"], "isReal(arccot(a))"),

        (["isReal(a)", "a > 0"], "isReal(log(a))"),
        (["isReal(a)"], "isReal(exp(a))"),
        (["isReal(a)"], "exp(a) > 0"),
        (["isReal(a)"], "isReal(abs(a))"),
        (["isReal(a)"], "isReal(sqrt(a))"),

        (["isComplex(a)", "isComplex(i)"], "isComplex(a + i)"),
        (["isComplex(a)", "isComplex(i)"], "isComplex(a - i)"),
        (["isComplex(a)", "isComplex(i)"], "isComplex(a * i)"),
        (["isComplex(a)", "isComplex(i)"], "isComplex(a / i)"),
        (["isComplex(a)"], "isComplex(-a)"),

        (["isEven(a)"], "isInt(a)"),
        (["isInt(a)"], "isReal(a)"),
        (["isReal(a)"], "isComplex(a)"),

        (["isReal(a)", "isReal(b)"], "isReal(a + b)"),
        (["isReal(a)", "isReal(b)"], "isReal(a - b)"),
        (["isReal(a)", "isReal(b)"], "isReal(a * b)"),
        (["isReal(a)", "isReal(b)","b != 0"], "isReal(a / b)"),
        
        (["isReal(a)"], "isReal(-a)"),

        # Real number power rules
        (["isReal(x)"], "isReal(x ^ n)"),
        (["isReal(x)", "isReal(y)"], "isReal(x ^ y)"),

    ]

    ineqs = []
    for conds, e in data:
        symb_e = expr_to_pattern(parse_expr(e))
        symb_conds = [expr_to_pattern(parse_expr(cond)) for cond in conds]
        ineqs.append(Identity(symb_e, conds=Conditions(symb_conds)))
    return ineqs

standard_inequalities = get_standard_inequalities()

def check_condition(e: Expr, ctx: Context) -> bool:
    """Check whether e holds under the given context."""

    ### Some special checks ###

    if expr.is_conj(e):
        return all(check_condition(arg, ctx) for arg in e.args)
    if expr.is_disj(e):
        return any(check_condition(arg, ctx) for arg in e.args)

    # If integrand is non-negative, then the integral is non-negative
    if expr.is_greater_eq(e) and expr.is_integral(e.args[0]) and e.args[1] == Const(0):
        ctx2 = Context(ctx)
        ctx2.add_condition(Op(">", Var(e.args[0].var), e.args[0].lower))
        ctx2.add_condition(Op("<", Var(e.args[0].var), e.args[0].upper))
        return check_condition(Op(">=", e.args[0].body, Const(0)), ctx2)
    
    # abs(s) < t <-- -t < s < t &&
    if expr.is_less(e) and expr.is_fun(e.args[0], 'abs'):
        arg = e.args[0].args[0]
        e1 = Op("<", arg, e.args[1])
        e2 = Op(">", arg, -e.args[1])
        return check_condition(e1, ctx) and check_condition(e2,ctx)

    # real vs. non-real
    if expr.is_not_equals(e) and e.rhs.is_constant() and approx_not_real(e.rhs):
        if check_condition(expr.isReal(e.lhs), ctx):
            return True

    # Substitute for equations in the context
    if ctx.get_substs():
        new_e = e
        for var, subst_e in reversed(ctx.get_substs()):
            new_e = new_e.subst(var, subst_e)
        if new_e != e:
            if check_condition(new_e, ctx):
                return True

    # a <= inf or a < inf
    if (expr.is_less(e) or expr.is_less_eq(e)) and expr.is_pos_inf(e.args[1]):
        return True

    # INT Real Condition
    def contains_i(e: Expr):
        if expr.is_fun(e) and e.func_name == 'i':
            return True
        if e.ty in (expr.VAR, expr.CONST, expr.SYMBOL, expr.INF):
            return False
        if e.ty in (expr.OP, expr.FUN):
            return any(contains_i(arg) for arg in e.args)
        if expr.is_integral(e):
            return contains_i(e.body) or contains_i(e.lower) or contains_i(e.upper)
        return False

    def add_integral_real_cond(e: Expr, all_conds: dict[Expr, list[Expr]]):
        if expr.is_integral(e):
            if not contains_i(e.body) and not contains_i(e.lower) and not contains_i(e.upper):
                if e not in all_conds:
                    all_conds[e] = []
                all_conds[e].append(Fun('isReal', e))
                
    # Otherwise, perform saturation search
    conds = ctx.get_conds()
    for _, g in ctx.get_all_subgoals().items():
        if expr.is_compare(g.expr) and not expr.is_equals(g.expr):
            satisfied = True
            for cond in g.conds.data:
                if cond in ctx.get_conds().data:
                    satisfied = True
            if satisfied:
                conds.add_condition(g.expr)
    all_conds = init_all_conds(conds)
    
    # Check all subexpressions of e
    def check_subexpr(e: Expr):
        if expr.is_integral(e):
            add_integral_real_cond(e, all_conds)
        if expr.is_op(e) or expr.is_fun(e):
            for arg in e.args:
                check_subexpr(arg)
        if expr.is_integral(e):
            check_subexpr(e.body)
            check_subexpr(e.lower)
            check_subexpr(e.upper)
    
    check_subexpr(e)
    
    ineqs = copy(standard_inequalities)
    ineqs.extend(ctx.get_inequalities())

    saturate(subject_of(e), ineqs, all_conds)
    return len(check_cond(e, all_conds, dict())) == 1

