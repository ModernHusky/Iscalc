"""Module for reasoning about conditions."""

from copy import copy
from typing import Dict, List

from integral import expr
from integral.expr import Expr, eval_expr, match, expr_to_pattern, Op, Const, Var
from integral.conditions import Conditions
from integral.context import Context, Identity
from integral.parser import parse_expr
from integral.poly import normalize


def subject_of(cond: Expr) -> Expr:
    """Return the subject of a condition.
    
    This is usually the left side of an inequality, or the (only)
    argument of a predicate.
    
    """
    if cond.is_equals() or cond.is_not_equals():
        return cond.args[0]
    if cond.is_greater() or cond.is_greater_eq():
        return cond.args[0]
    if cond.is_less() or cond.is_less_eq():
        return cond.args[0]
    # if cond.is_fun() and cond.func_name == 'isInt':
    #     return cond.args[0]
    # if cond.is_fun() and cond.func_name == 'isEven':
    #     return cond.args[0]
    if expr.is_fun(cond):
        return cond.args[0]
    raise TypeError

# Tolerance for floating-point rounding errors
tol = 1e-15

# Comparison of floating-point numbers up to rounding error
def approx_equal(a: Expr, b: Expr) -> bool:
    a, b = eval_expr(a), eval_expr(b)
    return abs(a - b) < tol

def approx_not_equal(a: Expr, b: Expr) -> bool:
    a, b = eval_expr(a), eval_expr(b)
    return abs(a - b) > tol

def approx_greater(a: Expr, b: Expr) -> bool:
    a, b = eval_expr(a), eval_expr(b)
    return a - b > tol

def approx_greater_eq(a: Expr, b: Expr) -> bool:
    a, b = eval_expr(a), eval_expr(b)
    return a - b > -tol

def approx_less(a: Expr, b: Expr) -> bool:
    a, b = eval_expr(a), eval_expr(b)
    return b - a > tol

def approx_less_eq(a: Expr, b: Expr) -> bool:
    a, b = eval_expr(a), eval_expr(b)
    return b - a > -tol

def approx_integer(a: Expr) -> bool:
    a = eval_expr(a)
    return abs(round(a) - a) < tol

def approx_even(a: Expr) -> bool:
    if approx_integer(a):
        return eval_expr(a) % 2 == 0
    else:
        return False

def init_all_conds(conds: Conditions) -> Dict[Expr, List[Expr]]:
    """Initialize all_conds from a condition object."""
    all_conds: Dict[Expr, List[Expr]] = dict()
    for cond in conds.data:
        x = subject_of(cond)
        if x not in all_conds:
            all_conds[x] = list()
        all_conds[x].append(cond)
        if expr.is_fun(x) and x.func_name == 'abs' and cond.is_less():
            if x.args[0] not in all_conds:
                all_conds[x.args[0]] = list()
            all_conds[x.args[0]].append(Op("<", x.args[0], cond.args[1]))
            all_conds[x.args[0]].append(Op(">", x.args[0], -cond.args[1]))
        if expr.is_fun(x) and x.func_name == 'abs' and cond.is_less_eq():
            if x.args[0] not in all_conds:
                all_conds[x.args[0]] = list()
            all_conds[x.args[0]].append(Op("<=", x.args[0], cond.args[1]))
            all_conds[x.args[0]].append(Op(">=", x.args[0], -cond.args[1]))

    # add simple condition transition
    for k in all_conds:
        for x in all_conds[k]:
            if x.is_less():
                if x.args[1] in all_conds:
                    for y in all_conds[x.args[1]]:
                        # x: k < b
                        # y: b < c or b <= c or b = c
                        # x and y ==> k < c
                        if y.is_less() or y.is_less_eq() or y.is_equals():
                            all_conds[k].append(Op('<', k, y.args[1]))
    return all_conds

def update_inst(k: str, v: Expr, inst: Dict[str, Expr]) -> Dict[str, Expr]:
    """Update instantiation without changing the original."""
    res = copy(inst)
    res[k] = v
    return res

def check_cond(cond: Expr, all_conds: Dict[Expr, List[Expr]], inst: Dict[str, Expr]) -> List[Dict[str, Expr]]:
    """Determine whether cond is implied by the existing set of conditions.
    
    The following checks are performed:

    - If subject of cond is a constant, and the right side is also constant,
      compare using eval_expr.

    - If subject of cond appears in all_conds, try to use the conditions
      available to verify cond.

    - Perform pattern matching.

    """
    x = subject_of(cond)

    # Trivial case
    if x in all_conds and cond in all_conds[x]:
        return [inst]

    # If subject of cond is a constant
    if x.is_constant():
        if cond.is_equals() and cond.args[1].is_constant():
            if approx_equal(x, cond.args[1]):
                return [inst]
        elif cond.is_not_equals() and cond.args[1].is_constant():
            if approx_not_equal(x, cond.args[1]):
                return [inst]
        elif cond.is_greater() and cond.args[1].is_constant():
            if approx_greater(x, cond.args[1]):
                return [inst]
        elif cond.is_greater_eq() and cond.args[1].is_constant():
            if approx_greater_eq(x, cond.args[1]):
                return [inst]
        elif cond.is_less() and cond.args[1].is_constant():
            if approx_less(x, cond.args[1]):
                return [inst]
        elif cond.is_less_eq() and cond.args[1].is_constant():
            if approx_less_eq(x, cond.args[1]):
                return [inst]
        elif expr.is_fun(cond) and cond.func_name == 'isInt':
            if approx_integer(x):
                return [inst]
        elif expr.is_fun(cond) and cond.func_name == 'isEven':
            if approx_even(x):
                return [inst]

    # If subject of cond appears in all_conds
    if cond.is_compare() and x in all_conds and cond.args[1].is_constant():
        for fact in all_conds[x]:
            if not (fact.is_compare() and fact.args[1].is_constant()):
                continue
            if cond.is_greater_eq():
                # x >= b --> b >= a --> x >= a
                if fact.is_greater() or fact.is_greater_eq():
                    if approx_greater_eq(fact.args[1], cond.args[1]):
                        return [inst]
            if cond.is_greater():
                # x >= b --> b > a --> x > a
                if fact.is_greater_eq() and approx_greater(fact.args[1], cond.args[1]):
                    return [inst]
                # x > b --> b >= a --> x > a
                if fact.is_greater() and approx_greater_eq(fact.args[1], cond.args[1]):
                    return [inst]
            if cond.is_less_eq():
                # x <= b --> b <= a --> x <= a
                if fact.is_less() or fact.is_less_eq():
                    if approx_less_eq(fact.args[1], cond.args[1]):
                        return [inst]
            if cond.is_less():
                # x <= b --> b < a --> x < a
                if fact.is_less_eq() and approx_less(fact.args[1], cond.args[1]):
                    return [inst]
                # x < b --> b <= a --> x < a
                if fact.is_less() and approx_less_eq(fact.args[1], cond.args[1]):
                    return [inst]
            if cond.is_equals():
                if fact.is_equals() and approx_equal(fact.args[1], cond.args[1]):
                    return [inst]
            if cond.is_not_equals():
                if fact.is_not_equals() and approx_equal(fact.args[1], cond.args[1]):
                    return [inst]
                # x < a --> a <= b --> x != b
                if fact.is_less() and approx_less_eq(fact.args[1], cond.args[1]):
                    return [inst]
                # x <= a --> a < b --> x != b
                if fact.is_less_eq() and approx_less(fact.args[1], cond.args[1]):
                    return [inst]
                # x > a --> a >= b --> x != b
                if fact.is_greater() and approx_greater_eq(fact.args[1], cond.args[1]):
                    return [inst]
                # x >= a --> a > b --> x != b
                if fact.is_greater_eq() and approx_greater(fact.args[1], cond.args[1]):
                    return [inst]
                # x = a --> a != b --> x != b
                if fact.is_equals() and approx_not_equal(fact.args[1], cond.args[1]):
                    return [inst]
        
    # If the other side of cond is a pattern
    if cond.is_compare() and x in all_conds and expr.is_symbol(cond.args[1]):
        symb = cond.args[1].name
        res = []
        for fact in all_conds[x]:
            if not fact.is_compare():
                continue
            if cond.is_greater_eq():
                if fact.is_greater_eq() or fact.is_greater():
                    res.append(update_inst(symb, fact.args[1], inst))
            if cond.is_greater():
                if fact.is_greater():
                    res.append(update_inst(symb, fact.args[1], inst))
            if cond.is_less_eq():
                if fact.is_less_eq() or fact.is_less():
                    res.append(update_inst(symb, fact.args[1], inst))
            if cond.is_less():
                if fact.is_less():
                    res.append(update_inst(symb, fact.args[1], inst))
            if cond.is_equals():
                if fact.is_equals():
                    res.append(update_inst(symb, fact.args[1], inst))
            if cond.is_not_equals():
                if fact.is_not_equals():
                    res.append(update_inst(symb, fact.args[1], inst))
        return res

    # Not found
    return list()

def saturate_expr(e: Expr, ineq: Identity, all_conds: Dict[Expr, List[Expr]], ctx: Context):
    """Use the rule ineq to saturate facts about e, add to all_conds."""
    pat = subject_of(ineq.expr)
    inst = match(e, pat)
    if inst is not None:
        old_list = [inst]
        for cond in ineq.conds.data:
            new_list = []
            for inst in old_list:
                res = check_cond(cond.inst_pat(inst), all_conds, inst)
                new_list.extend(res)
            old_list = new_list
        for mapping in old_list:
            res = ineq.expr.inst_pat(mapping)
            res_norm = Op(res.op ,res.args[0], normalize(res.args[1], ctx)) if res.is_compare() else res
            # if check_cond(res, all_conds, mapping) == [mapping]:
            #     continue
            if e not in all_conds:
                all_conds[e] = list()
            if res not in all_conds[e]:
                all_conds[e].append(res)
            if res_norm not in all_conds[e]:
                all_conds[e].append(res_norm)
    return

def saturate_once(e: Expr, ineqs: List[Identity], all_conds: Dict[Expr, List[Expr]], ctx: Context):
    """Perform one round of saturation"""
    all_subs = e.find_all_subexpr()
    for sube, _ in all_subs:
        for ineq in ineqs:
            saturate_expr(sube, ineq, all_conds, ctx)

def all_conds_size(all_conds: Dict[Expr, List[Expr]]) -> int:
    """Return number of facts in all_conds."""
    res = 0
    for _, conds in all_conds.items():
        res += len(conds)
    return res        

def saturate(e: Expr, ineqs: List[Identity], all_conds: Dict[Expr, List[Expr]], ctx: Context, *,
             round_limit: int = 5, size_limit: int = 1000):
    """Saturate up to given number of rounds and size limits.
    
    If number of rounds and size limits have been reached without
    saturation, assertion is thrown to alert possible problems.
    
    """
    i = 0
    while True:
        prev_size = all_conds_size(all_conds)
        saturate_once(e, ineqs, all_conds, ctx)
        i += 1
        next_size = all_conds_size(all_conds)
        if prev_size == next_size:
            return
        if next_size > size_limit or i > round_limit:
            print_all_conds(all_conds)
            raise AssertionError("saturate: limit reached")

def print_all_conds(all_conds: Dict[Expr, List[Expr]]):
    for x, conds in all_conds.items():
        print("%s: %s" % (x, ', '.join(str(cond) for cond in conds)))

def get_standard_inequalities() -> List[Identity]:
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
        (["a >= b", "c > d"], "a + c > b + d"),
        (["a > b", "c >= d"], "a + c > b + d"),
        (["a <= b", "c < d"], "a + c < b + d"),
        (["a < b", "c <= d"], "a + c < b + d"),
        (["a >= b", "c >= d"], "a + c >= b + d"),
        (["a <= b", "c <= d"], "a + c <= b + d"),
        ([], "a + b = b + a"),

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
        (["x > 1"], "1 / x < 1"),
        (["x > 0"], "1 / x > 0"),

        # Square root
        (["a > 0"], "sqrt(a) > 0"),
        (["a < 1", "a > 0"], "sqrt(a) < 1"),
        (["a > b", "b >= 0"], "sqrt(a) > sqrt(b)"),
        (["a >= b", "b >= 0"], "sqrt(a) >= sqrt(b)"),

        # Power
        (["a != 0"], "a ^ 2 > 0"),
        ([], "a ^ 2 >= 0"),
        (["x > 0"], "x ^ y > 0"),
        (["x != 0"], "x ^ n != 0"),
        (["x > y", "y >= 0", "z > 0"], "x ^ z > y ^ z"),
        (["x < a", "x > -a"], "x ^ 2 < a ^ 2"),
        (["x > a", "a >= 0"], "x ^ 2 > a ^ 2"),
        (["x <= a", "x >= -a"], "x ^ 2 <= a ^ 2"),
        (["x >= a", "a >= 0"], "x ^ 2 >= a ^ 2"),
        (["x != y"], "x ^ 2 - y ^ 2 != 0"),
        (["y != x"], "x ^ 2 - y ^ 2 != 0"),
        (["x != y"], "x ^ 4 - y ^ 4 != 0"),
        (["y != x"], "x ^ 4 - y ^ 4 != 0"),

        # Log
        (["x >= 1"], "log(x) >= 0"),
        (["x > 1"], "log(x) > 0"),
        (["x <= 1", "x > 0"], "log(x) <= 0"),
        (["x < 1", "x > 0"], "log(x) < 0"),
        (["x != 1"], "log(x) != 0"),

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
        (["x > pi / 2", "x < 3 * pi / 2"], "cos(x) < 0"),
        (["x >= -pi / 2", "x <= pi / 2"], "cos(x) >= 0"),
        (["x >= pi / 2", "x <= 3 * pi / 2"], "cos(x) <= 0"),
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
        (["x > -pi / 2", "x < pi / 2"], "sec(x) > 0"),
        (["x > pi / 4", " x < pi / 2"], "sec(x) > 1"),
        (["x > pi / 4", " x < pi / 2"], "sec(x) < sqrt(2)"),
        # Inverse trigonometric
        (["x >= -1", "x <= 1"], "asin(x) >= -pi / 2"),
        (["x >= -1", "x <= 1"], "asin(x) <= pi / 2"),
        (["x > -1", "x < 1"], "asin(x) > -pi / 2"),
        (["x > -1", "x < 1"], "asin(x) < pi / 2"),
        (["x > 0", "x <= 1"], "asin(x) > 0"),
        (["x >= 0", "x <= 1"], "asin(x) >= 0"),
        (["x < 0", "x >= -1"], "asin(x) < 0"),
        (["x <= 0", "x >= -1"], "asin(x) <= 0"),
        (["x >= -1", "x <= 1"], "acos(x) >= 0"),
        (["x >= -1", "x <= 1"], "acos(x) <= pi"),
        (["x > 0", "x <= 1"], "acos(x) < pi / 2"),
        (["x >= 0", "x <= 1"], "acos(x) <= pi / 2"),
        (["x < 0", "x >= -1"], "acos(x) > pi / 2"),
        (["x <= 0", "x >= -1"], "acos(x) >= pi / 2"),
        (["x > -1", "x < 1"], "acos(x) > 0"),
        (["x > -1", "x < 1"], "acos(x) < pi"),
        ([], "atan(x) > -pi / 2"),
        ([], "atan(x) < pi / 2"),
        (["x >= 0"], "atan(x) >= 0"),
        (["x > 0"], "atan(x) > 0"),
        (["x <= 0"], "atan(x) <= 0"),
        (["x < 0"], "atan(x) < 0"),
        (["x >= 1"], "asec(x) >= 0"),
        (["x >= 1"], "asec(x) < pi / 2"),
        (["x > 1"], "asec(x) > 0"),
        (["x <= -1"], "asec(x) <= pi"),
        (["x <= -1"], "asec(x) > pi / 2"),
        (["x < -1"], "asec(x) < pi"),
        (["x >= 1"], "acsc(x) <= pi / 2"),
        (["x >= 1"], "asec(x) > 0"),
        (["x > 1"], "asec(x) < pi / 2"),
        (["x <= -1"], "asec(x) >= -pi / 2"),
        (["x <= -1"], "asec(x) < 0"),
        (["x < -1"], "asec(x) > -pi / 2"),
        (["x > 1"], "atan(x) > pi / 4"),
        (["x > 0"], "atan(x) < pi / 2"),
        (["x > 0"], "atan(x) > 0"),
        ([], "atan(x) > -pi/2"),
        ([], "atan(x) < pi/2"),
        # Hyperbolic
        ([], "cosh(x) > 0"),

        # Factorial
        ([], "factorial(x) >= 1"),

        # reflexivity
        ([], "x = x"),

        # isInt
        (["isInt(a)", "isInt(b)"], "isInt(a + b)"),
        (["isInt(a)", "isInt(b)"], "isInt(a - b)"),
        (["isInt(a)", "isInt(b)"], "isInt(a * b)"),

        (["a > 0"], "a != 0"),
        (["a < 0"], "a != 0"),
        (["a >= b", "a != b"], "a > b"),
        (["a <= b", "a != b"], "a < b"),
        (["a = b", "a > c"], "b > c"),
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

    # Some special checks
    if e.is_greater_eq() and expr.is_integral(e.args[0]) and e.args[1] == Const(0):
        ctx2 = Context(ctx)
        ctx2.add_condition(Op(">", Var(e.args[0].var), e.args[0].lower))
        ctx2.add_condition(Op("<", Var(e.args[0].var), e.args[0].upper))
        return check_condition(Op(">=", e.args[0].body, Const(0)), ctx2)
    if e.is_less() and expr.is_fun(e.args[0]) and e.args[0].func_name == 'abs':
        arg = e.args[0].args[0]
        e1 = Op("<", arg, e.args[1])
        e2 = Op(">", arg, -e.args[1])
        return check_condition(e1, ctx) and check_condition(e2,ctx)
    if ctx.get_substs():
        new_e = e
        for var, subst_e in ctx.get_substs().items():
            new_e = new_e.subst(var, subst_e)
        if new_e != e:
            if check_condition(new_e, ctx):
                return True
    # a <= inf or a < inf
    if (e.is_less() or e.is_less_eq()) and expr.is_pos_inf(e.args[1]):
        return True

    conds = ctx.get_conds()
    all_conds = init_all_conds(conds)
    ineqs = copy(standard_inequalities)
    ineqs.extend(ctx.get_inequalities())
    for lemma in ctx.get_lemmas():
        if lemma.expr.is_compare():
            ineqs.append(lemma)

    saturate(subject_of(e), ineqs, all_conds, ctx)
    return len(check_cond(e, all_conds, dict())) == 1
