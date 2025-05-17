"""Expressions."""
import math
import functools
import operator
from decimal import Decimal
from fractions import Fraction
from collections.abc import Iterable
from typing import Dict, List, Optional, Set, TypeGuard, Tuple, Union, Callable

import sympy

class IscalcException(Exception):
    """Parent class of all exceptions in Iscalc."""
    def to_json(self) -> dict:
        """Convert current object to json format."""
        raise NotImplementedError(f"to_json: {type(self)}")
    
    @staticmethod
    def from_json(data: dict):
        """Load object from json format."""
        raise NotImplementedError(f"from_json: {__class__.__name__}")


VAR, CONST, OP, FUN, DERIV, INTEGRAL,CINTEGRAL, EVAL_AT, SYMBOL, LIMIT, INF, INDEFINITEINTEGRAL, \
SKOLEMFUNC, SUMMATION, PRODUCT, MULTIPOLECONTOUR, COMPOUNDCONTOUR, CIRCLEPATH, POLEPATH, RECTANGLEPATH, LINEPATH = range(21)

op_priority = {
    "+": 65, "-": 65, "*": 70, "/": 70, "%": 70, "^": 75, "=": 50, "<": 50, ">": 50, "<=": 50, ">=": 50, "!=": 50
}

class Location:
    """Location within an expression."""

    def __init__(self, data):
        if isinstance(data, Iterable) and all(isinstance(n, int) for n in data):
            self.data = tuple(data)
        elif isinstance(data, str):
            if data in (".", ""):
                self.data = tuple([])
            else:
                self.data = tuple(int(n) for n in data.split('.'))
        elif isinstance(data, Location):
            self.data = data.data
        else:
            raise TypeError

    def __str__(self):
        if not self.data:
            return "."
        else:
            return ".".join(str(n) for n in self.data)

    def is_empty(self):
        return len(self.data) == 0

    @property
    def head(self):
        return self.data[0]

    @property
    def top2(self):
        return self.data[:2]
    @property
    def rest(self):
        return Location(self.data[1:])

    @property
    def rest2(self):
        return Location(self.data[2:])

    def append(self, i: int) -> "Location":
        return Location(self.data + (i,))


class Expr:
    """Expressions."""
    def __init__(self):
        pass

    def __add__(self, other):
        if isinstance(other, (int, Fraction)):
            other = Const(other)
        return Op("+", self, other)

    def __radd__(self, other):
        if isinstance(other, (int, Fraction)):
            other = Const(other)
        return Op("+", other, self)

    def __sub__(self, other):
        if isinstance(other, (int, Fraction)):
            other = Const(other)
        return Op("-", self, other)

    def __rsub__(self, other):
        if isinstance(other, (int, Fraction)):
            other = Const(other)
        return Op("-", other, self)

    def __mul__(self, other):
        if isinstance(other, (int, Fraction)):
            other = Const(other)
        return Op("*", self, other)

    def __rmul__(self, other):
        if isinstance(other, (int, Fraction)):
            other = Const(other)
        return Op("*", other, self)

    def __truediv__(self, other):
        if isinstance(other, (int, Fraction)):
            other = Const(other)
        if is_const(self) and is_const(other) and isinstance(self.val, int) and isinstance(other.val, int):
            return Const(Fraction(self.val, other.val))
        return Op("/", self, other)

    def __rtruediv__(self, other):
        if isinstance(other, (int, Fraction)):
            other = Const(other)
        return Op("/", other, self)

    def __xor__(self, other):
        if isinstance(other, (int, Fraction)):
            other = Const(other)
        return Op("^", self, other)

    def __pow__(self, other):
        if isinstance(other, (int, Fraction)):
            other = Const(other)
        return Op("^", self, other)

    def __mod__(self, other):
        if isinstance(other, (int, Fraction)):
            other = CONST(other)
        return Op("%", self, other)

    def __neg__(self):
        if self == POS_INF:
            return NEG_INF
        elif self == NEG_INF:
            return POS_INF
        elif is_const(self) and self.val > 0:
            return Const(-self.val)
        return Op("-", self)

    def size(self):
        if self.ty in (VAR, CONST, SYMBOL, INF):
            return 1
        elif is_op(self) or is_fun(self):
            return 1 + sum(arg.size() for arg in self.args)
        elif is_deriv(self):
            return 1 + self.body.size()
        elif is_integral(self) or is_evalat(self):
            return 1 + self.lower.size() + self.upper.size() + self.body.size()
        elif is_limit(self):
            return 1 + self.lim.size() + self.body.size()
        elif is_indefinite_integral(self):
            return 1 + self.body.size()
        elif is_skolem_func(self):
            return 1 + len(self.dependent_vars)
        elif is_summation(self):
            return 1 + self.lower.size() + self.upper.size() + self.body.size()
        elif is_product(self):
            return 1 + self.lower.size() + self.upper.size() + self.body.size()
        elif is_cintegral(self):
            return 1 + sum(path.size() for path in self.paths) + self.body.size()        
        else:
            raise NotImplementedError

    def is_zero(self) -> bool:
        return is_const(self) and self.val == 0

    def is_plus(self):
        return self.ty == OP and self.op == '+'

    def is_minus(self):
        return self.ty == OP and self.op == '-' and len(self.args) == 2

    def is_uminus(self):
        return self.ty == OP and self.op == '-' and len(self.args) == 1

    def is_times(self):
        return self.ty == OP and self.op == '*'

    def is_divides(self):
        return self.ty == OP and self.op == '/'

    def is_power(self):
        return self.ty == OP and self.op == '^'

    def is_mod(self):
        return self.ty == OP and self.op == '%'

    def is_equals(self):
        return self.ty == OP and self.op == '='

    def is_not_equals(self):
        return self.ty == OP and self.op == "!="

    def is_less(self):
        return self.ty == OP and self.op == "<"

    def is_less_eq(self):
        return self.ty == OP and self.op == "<="

    def is_greater(self):
        return self.ty == OP and self.op == ">"

    def is_greater_eq(self):
        return self.ty == OP and self.op == ">="

    def is_compare(self) -> bool:
        return self.ty == OP and self.op in ('=', '!=', '<', '<=', '>', '>=')

    def is_trig(self):
        return self.ty == FUN and self.func_name in ("sin", "cos", "tan", "cot", "csc", "sec")

    def is_inverse_trig(self):
        return self.ty == FUN and self.func_name in ("arcsin", "arccos", "arctan", "arccot", "arccsc", "arcsec")

    def is_skolem_term(self):
        if self.get_vars() != set():
            return False
        if not self.is_constant():
            return True
        else:
            return False

    def is_odd(self, var, conds) -> bool:
        from integral import poly
        tmp1 = self
        tmp2 = self.subst(var, -Var(var))
        if poly.normalize(tmp1 + tmp2, conds) == Const(0):
            return True
        return False

    @property
    def lhs(self) -> "Expr":
        if self.is_compare():
            return self.args[0]
        else:
            raise AssertionError(f"lhs: term {self} is not a comparison")

    @property
    def rhs(self) -> "Expr":
        if self.is_compare():
            return self.args[1]
        else:
            raise AssertionError(f"rhs: term {self} is not a comparison")

    def __le__(self, other):
        if isinstance(other, (int, Fraction)):
            return False

        if self.size() != other.size():
            return self.size() <= other.size()

        if self.ty != other.ty:
            return self.ty <= other.ty

        if is_var(self):
            return self.name <= other.name
        elif is_const(self):
            return self.val <= other.val
        elif is_op(self):
            return (self.op, self.args) <= (other.op, other.args)
        elif is_fun(self):
            return (self.func_name, self.args) <= (other.func_name, other.args)
        elif is_deriv(self) or is_indefinite_integral(self):
            return (self.body, self.var) <= (other.body, other.var)
        elif is_integral(self) or is_evalat(self):
            return (self.body, self.lower, self.upper, self.var) <= \
                   (other.body, other.lower, other.upper, other.var)
        elif is_symbol(self):
            return sum(self.ty) <= sum(other.ty)
        elif is_summation(self):
            return (self.body, self.lower, self.upper, self.index_var) <= \
                   (other.body, other.lower, other.upper, other.index_var)
        elif is_skolem_func(self):
            return (self.name, self.dependent_vars) <= (other.name, other.dependent_vars)
        elif is_limit(self):
            return (self.var, self.lim, self.body, self.drt) <= (other.var, other.lim, other.body, other.drt)
        else:
            print(type(self))
            raise NotImplementedError

    def __lt__(self, other):
        return self <= other and self != other

    def __gt__(self, other):
        return other <= self and self != other

    def __ge__(self, other):
        return not self < other

    def priority(self):
        if self.ty in (VAR, SYMBOL, INF, SKOLEMFUNC):
            return 100
        elif self.ty == CONST:
            if isinstance(self.val, Fraction) and self.val.denominator != 1:
                return op_priority['/']
            elif self.val < 0:
                # return 80  # priority of uminus
                return 74
            else:
                return 100
        elif self.ty == OP:
            if len(self.args) == 1:
                return 80  # priority of uminus
            elif self.op in op_priority:
                return op_priority[self.op]
            else:
                raise NotImplementedError
        elif self.ty in (FUN, SUMMATION, PRODUCT):
            return 95
        elif self.ty in (DERIV, INTEGRAL, EVAL_AT, INDEFINITEINTEGRAL, CINTEGRAL, COMPOUNDCONTOUR):
            return 10
        elif self.ty == LIMIT:
            return 5
        else:
            raise NotImplementedError

    def __lt__(self, other):
        return self <= other and self != other

    def get_subexpr(self, loc) -> "Expr":
        """Given an expression, return the subexpression at location."""
        if not isinstance(loc, Location):
            loc = Location(loc)
        if loc.is_empty():
            return self
        elif is_var(self) or is_const(self):
            raise AssertionError("get_subexpr: invalid location")
        elif is_op(self) or is_fun(self):
            assert loc.head < len(self.args), "get_subexpr: invalid location"
            return self.args[loc.head].get_subexpr(loc.rest)
        elif is_deriv(self):
            assert loc.head == 0, "get_subexpr: invalid location"
            return self.body.get_subexpr(loc.rest)
        elif is_integral(self) or is_evalat(self):
            if loc.head == 0:
                return self.body.get_subexpr(loc.rest)
            elif loc.head == 1:
                return self.lower.get_subexpr(loc.rest)
            elif loc.head == 2:
                return self.upper.get_subexpr(loc.rest)
            else:
                raise AssertionError("get_subexpr: invalid location")
        elif is_limit(self):
            assert loc.head == 0, "get_subexpr: invalid location"
            return self.body.get_subexpr(loc.rest)

        else:
            raise NotImplementedError

    def replace_expr(self, loc, new_expr: "Expr") -> "Expr":
        """Replace self's subexpr at location."""
        if not isinstance(loc, Location):
            loc = Location(loc)
        if loc.is_empty():
            return new_expr
        elif is_var(self) or is_const(self):
            raise AssertionError("replace_expr: invalid location")
        elif is_op(self):
            assert loc.head < len(self.args), "replace_expr: invalid location"
            if len(self.args) == 1:
                return Op(self.op, self.args[0].replace_expr(loc.rest, new_expr))
            elif len(self.args) == 2:
                if loc.head == 0:
                    return Op(self.op, self.args[0].replace_expr(loc.rest, new_expr), self.args[1])
                elif loc.head == 1:
                    return Op(self.op, self.args[0], self.args[1].replace_expr(loc.rest, new_expr))
                else:
                    raise AssertionError("replace_expr: invalid location")
            else:
                raise NotImplementedError
        elif is_fun(self):
            assert loc.head < len(self.args), "replace_expr: invalid location"
            arg = self.args[loc.head].replace_expr(loc.rest, new_expr)
            return Fun(self.func_name, arg)
        elif is_integral(self):
            if loc.head == 0:
                return Integral(self.var, self.lower, self.upper, self.body.replace_expr(loc.rest, new_expr))
            elif loc.head == 1:
                return Integral(self.var, self.lower.replace_expr(loc.rest, new_expr), self.upper, self.body)
            elif loc.head == 2:
                return Integral(self.var, self.lower, self.upper.replace_expr(loc.rest, new_expr), self.body)
            else:
                raise AssertionError("replace_expr: invalid location")
        elif is_evalat(self):
            if loc.head == 0:
                return EvalAt(self.var, self.lower, self.upper, self.body.replace_expr(loc.rest, new_expr))
            elif loc.head == 1:
                return EvalAt(self.var, self.lower.replace_expr(loc.rest, new_expr), self.upper, self.body)
            elif loc.head == 2:
                return EvalAt(self.var, self.lower, self.upper.replace_expr(loc.rest, new_expr), self.body)
            else:
                raise AssertionError("replace_expr: invalid location")
        elif is_deriv(self):
            assert loc.head == 0, "replace_expr: invalid location"
            return Deriv(self.var, self.body.replace_expr(loc.rest, new_expr))
        elif is_limit(self):
            assert loc.head == 0, "replace_expr: invalid location"
            return Limit(self.var, self.lim, self.body.replace_expr(loc.rest, new_expr), self.drt)
        elif is_summation(self):
            if loc.head == 0:
                return Summation(self.index_var, self.lower, self.upper, self.body.replace_expr(loc.rest, new_expr))
            elif loc.head == 1:
                return Summation(self.index_var, self.lower.replace_expr(loc.rest, new_expr), self.upper, self.body)
            elif loc.head == 2:
                return Summation(self.index_var, self.lower, self.upper.replace_expr(loc.rest, new_expr), self.body)
            else:
                raise AssertionError("replace_expr: invalid location")
        elif is_cintegral(self):
            return CIntegral(self.var, self.paths, self.body.replace(loc.rest, new_expr))
        else:
            raise NotImplementedError(self)

    def get_location(self) -> Location:
        """Returns the location at which the 'selected' field is True."""
        location = []

        def get(exp: Expr, loc=''):
            if hasattr(exp, 'selected') and exp.selected == True:
                location.append(loc[1:])
                exp.selected = False  # Once it is found, restore it.
            elif is_op(exp) or is_fun(exp):
                for i in range(len(exp.args)):
                    get(exp.args[i], loc + "." + str(i))
            elif is_integral(exp) or is_evalat(exp):
                get(exp.lower, loc + ".1")
                get(exp.upper, loc + ".2")
                get(exp.body, loc + ".0")
            elif is_deriv(exp) or is_summation(exp) or is_limit(exp):
                get(exp.body, loc + ".0")

        get(self)
        return location[0]
    def get_all_func_name(self) -> Set[str]:
        return set([pair[0].func_name for pair in self.find_all_subexpr() if is_fun(pair[0])])

    def get_all_symbols(self) -> Set["Symbol"]:
        return set([pair[0] for pair in self.find_all_subexpr() if is_symbol(pair[0])])

    def find_subexpr(self, subexpr: "Expr") -> List[Location]:
        """Returns the location of a subexpression."""
        locations = []

        def find(e: Expr, loc: Location):
            if e == subexpr:
                locations.append(Location(loc))
            elif is_op(e) or is_fun(e):
                for i, arg in enumerate(e.args):
                    find(arg, loc.append(i))
            elif is_integral(e) or is_evalat(e):
                find(e.lower, loc.append(1))
                find(e.upper, loc.append(2))
                find(e.body, loc.append(0))
            elif is_deriv(e) or is_limit(e) or is_indefinite_integral(e):
                find(e.body, loc.append(0))
            elif is_summation(e):
                find(e.body, loc.append(0))
                find(e.lower, loc.append(1))
                find(e.upper, loc.append(2))
        find(self, Location(""))
        return locations

    def find_subexpr_pred(self, pred: Callable[["Expr"], bool]) -> List[Tuple["Expr", Location]]:
        """Find list of subexpressions satisfying a given predicate.

        Larger expressions are placed later.

        """
        results = []

        def find(e: Expr, loc: Location):
            if is_op(e) or is_fun(e):
                for i, arg in enumerate(e.args):
                    find(arg, loc.append(i))
            elif is_integral(e) or is_evalat(e):
                find(e.lower, loc.append(1))
                find(e.upper, loc.append(2))
                find(e.body, loc.append(0))
            elif is_deriv(e) or is_limit(e) or is_indefinite_integral(e) or is_cintegral(e):
                find(e.body, loc.append(0))
            elif is_summation(e):
                find(e.body, loc.append(0))
                find(e.lower, loc.append(1))
                find(e.upper, loc.append(2))

            if pred(e):
                results.append((e, Location(loc)))

        find(self, Location(""))
        return results

    def find_all_subexpr(self) -> list[tuple["Expr", Location]]:
        return self.find_subexpr_pred(lambda t: True)

    def subst(self, var: str, e: "Expr") -> "Expr":
        """Substitute occurrence of var for e in self."""
        assert isinstance(var, str) and isinstance(e, Expr)
        if is_var(self):
            if self.name == var:
                return e
            else:
                return self
        elif is_const(self):
            return self
        elif is_skolem_func(self):
            return SkolemFunc(self.name, tuple(arg.subst(var, e) for arg in self.dependent_vars))
        elif is_symbol(self):
            return self
        elif is_op(self):
            return Op(self.op, *[arg.subst(var, e) for arg in self.args])
        elif is_fun(self):
            return Fun(self.func_name, *[arg.subst(var, e) for arg in self.args])
        elif is_deriv(self):
            return Deriv(self.var, self.body.subst(var, e))
        elif is_limit(self):
            return Limit(self.var, self.lim.subst(var, e), self.body.subst(var, e))
        elif is_inf(self):
            return self
        elif is_integral(self):
            return Integral(self.var, self.lower.subst(var, e), self.upper.subst(var, e), self.body.subst(var, e))
        elif is_indefinite_integral(self):
            return IndefiniteIntegral(self.var, self.body.subst(var, e), self.skolem_args)
        elif is_cintegral(self):
            return CIntegral(self.var, self.paths, self.body.subst(var,e))
        elif is_evalat(self):
            return EvalAt(self.var, self.lower.subst(var, e), self.upper.subst(var, e), self.body.subst(var, e))
        elif is_summation(self):
            return Summation(self.index_var, self.lower.subst(var, e), self.upper.subst(var, e),
                             self.body.subst(var, e))
        elif is_product(self):
            return Product(self.index_var, self.lower.subst(var, e), self.upper.subst(var, e),
                             self.body.subst(var, e))
        else:
            raise NotImplementedError(f"subst: {type(self)}")

    def is_constant(self):
        """Determine whether expr is a number.

        Note Inf is not considered to be constants.

        """
        if is_const(self):
            return True
        elif is_op(self):
            return all(arg.is_constant() for arg in self.args)
        elif is_fun(self):
            self: Fun
            if self.func_name in ('inv', 'unit_matrix', 'zero_matrix'):
                return False
            elif self.func_name == 'i':
                return True
            return all(arg.is_constant() for arg in self.args)
        else:
            return False

    def is_evaluable(self):
        return self.is_constant() or is_inf(self)
    
    def is_closed_form(self):
        """Determine whether expression is in closed form."""
        if is_const(self):
            return True
        elif is_var(self):
            return True
        elif is_op(self):
            return all(arg.is_closed_form() for arg in self.args)
        elif is_fun(self):
            return all(arg.is_closed_form() for arg in self.args)
        elif is_skolem_func(self):
            return True
        else:
            return False

    def get_vars(self) -> Set[str]:
        """Obtain the set of variables in self."""
        res = set()
        def rec(t, bd_vars):
            if is_var(t):
                if t.name not in bd_vars:
                    res.add(t.name)
            elif t.ty in (CONST, INF, SYMBOL):
                return
            elif is_op(t) or is_fun(t):
                for arg in t.args:
                    rec(arg, bd_vars)
            elif is_deriv(t):
                rec(t.body, bd_vars + [t.var])
            elif is_limit(t):
                rec(t.lim, bd_vars + [t.var])
                rec(t.body, bd_vars + [t.var])
            elif is_integral(t) or is_evalat(t):
                rec(t.lower, bd_vars + [t.var])
                rec(t.upper, bd_vars + [t.var])
                rec(t.body, bd_vars + [t.var])
            elif is_indefinite_integral(t):
                rec(t.body, bd_vars + [t.var])
            elif is_summation(t):
                rec(t.lower, bd_vars + [t.index_var])
                rec(t.upper, bd_vars + [t.index_var])
                rec(t.body, bd_vars + [t.index_var])
            elif is_product(t):
                rec(t.lower, bd_vars + [t.index_var])
                rec(t.upper, bd_vars + [t.index_var])
                rec(t.body, bd_vars + [t.index_var])
            elif t.is_equals():
                rec(t.lhs, bd_vars)
                rec(t.rhs, bd_vars)
            elif is_skolem_func(t):
                t:SkolemFunc
                for var in t.dependent_vars:
                    rec(var, bd_vars)
            elif is_cintegral(t):
                t:CompoundContourIntegral   #TODO: add other types of contour integrals
                rec(t.body, bd_vars + [t.var])
            else:
                print(t, type(t))
                raise NotImplementedError
        bd = []
        rec(self, bd)
        return res

    def contains_var(self, x: str) -> bool:
        """Whether self contains variable x."""
        assert isinstance(x, str)
        return x in self.get_vars()

    def contains_skolem_func(self):
        if is_skolem_func(self):
            return True
        elif is_integral(self) or is_indefinite_integral(self):
            return False
        elif is_op(self) or is_fun(self):
            return any(arg.contains_skolem_func() for arg in self.args)
        else:
            return False

    def replace(self, e: "Expr", repl_e: "Expr") -> "Expr":
        """Replace occurrences of e with repl_e."""
        assert isinstance(e, Expr) and isinstance(repl_e, Expr)
        if self == e:
            return repl_e
        elif self.ty in (VAR, CONST, INF, SYMBOL):
            return self
        elif is_op(self):
            return Op(self.op, *[arg.replace(e, repl_e) for arg in self.args])
        elif is_fun(self):
            return Fun(self.func_name, *[arg.replace(e, repl_e) for arg in self.args])
        elif is_deriv(self):
            return Deriv(self.var, self.body.replace(e, repl_e))
        elif is_integral(self):
            return Integral(self.var, self.lower.replace(e, repl_e), self.upper.replace(e, repl_e),
                            self.body.replace(e, repl_e))
        elif is_evalat(self):
            return EvalAt(self.var, self.lower.replace(e, repl_e), self.upper.replace(e, repl_e),
                          self.body.replace(e, repl_e))
        elif is_skolem_func(self):
            return SkolemFunc(self.name, tuple(var.replace(e, repl_e) for var in self.dependent_vars))
        elif is_summation(self):
            nl = self.lower.replace(e, repl_e)
            nu = self.upper.replace(e, repl_e)
            nbody = self.body.replace(e, repl_e)
            return Summation(self.index_var, nl, nu, nbody)
        elif is_product(self):
            nl = self.lower.replace(e, repl_e)
            nu = self.upper.replace(e, repl_e)
            nbody = self.body.replace(e, repl_e)
            return Product(self.index_var, nl, nu, nbody)
        elif is_limit(self):
            return Limit(self.var, self.lim.replace(e, repl_e), self.body.replace(e, repl_e), self.drt)
        elif is_cintegral(self):
            return CIntegral(self.var, self.paths, self.body.replace(e, repl_e))
        else:
            print(self, e, repl_e)
            raise NotImplementedError

    def separate_integral(self) -> list[tuple[Union["Integral", "IndefiniteIntegral"], Location]]:
        """Collect the list of all integrals appearing in self."""
        return self.find_subexpr_pred(lambda e: is_integral(e) or is_indefinite_integral(e))

    def separate_limits(self) -> list[tuple["Limit", Location]]:
        """Collect the list of all limits appearing in self."""
        return self.find_subexpr_pred(lambda e: is_limit(e))
    
    def separate_cintegral(self) -> List[Tuple["Expr", Location]]:
        """Collect the list of all integrals appearing in self."""
        return self.find_subexpr_pred(lambda e: is_cintegral(e))

    @property
    def depth(self):
        """Return the depth of expression as an estimate of problem difficulty."""
        def d(expr):
            if expr.ty in (VAR, CONST):
                return 0
            elif expr.ty in (OP, FUN):
                if len(expr.args) == 0:
                    return 1
                return 1 + max([d(expr.args[i]) for i in range(len(expr.args))])
            elif expr.ty in (EVAL_AT, INTEGRAL, DERIV):
                return d(expr.body)
            elif expr.ty == SYMBOL:
                raise TypeError

        return d(self)

    def is_spec_function(self, fun_name):
        """Return true iff e is formed by rational options of fun_name."""
        v = Symbol("v", [VAR, OP, FUN])
        if fun_name == "sin":
            pat1 = sin(v)
        elif fun_name == "cos":
            pat1 = cos(v)
        else:
            return False
        if len(find_pattern(self, pat1)) != 1:
            return False

        def rec(ex):
            if ex.ty == CONST:
                return True
            elif ex.ty == VAR:
                return False
            elif ex.ty == OP:
                return all(rec(arg) for arg in ex.args)
            elif ex.ty == FUN:
                return True if ex.func_name == fun_name else False
            else:
                return False

        return rec(self)

    def nonlinear_subexpr(self):
        """Return nonlinear & nonconstant subexpression."""
        subs = []
        a = Symbol('a', [CONST])
        b = Symbol('b', [CONST])
        x = Symbol('x', [VAR])
        patterns = [a * x, a * x + b, a * x - b, x, b + a * x, a + x, x + a]

        def traverse(exp):
            table = [match(exp, p) for p in patterns]
            is_linear = functools.reduce(lambda x, y: x or y, table)
            if not exp.is_constant() and not is_linear:
                if exp not in subs:
                    subs.append(exp)
            if exp.ty in (OP, FUN):
                for arg in exp.args:
                    traverse(arg)
            elif exp.ty in (INTEGRAL, EVAL_AT, DERIV):
                traverse(exp.body)

        traverse(self)
        if self in subs:
            subs.remove(self)
        return tuple(subs)

    def inst_pat(self, mapping: Dict) -> "Expr":
        """Instantiate by replacing symbols in term with mapping."""
        if is_var(self) or is_const(self) or is_inf(self):
            return self
        elif is_symbol(self):
            if self.name in mapping:
                res = mapping[self.name]
                return res
            else:
                return self
        elif is_op(self):
            return Op(self.op, *(arg.inst_pat(mapping) for arg in self.args))
        elif is_fun(self):
            return Fun(self.func_name, *(arg.inst_pat(mapping) for arg in self.args))
        elif is_skolem_func(self):
            return SkolemFunc(self.name, tuple(arg.inst_pat(mapping) for arg in self.dependent_vars))
        elif is_integral(self):
            return Integral(self.var, self.lower.inst_pat(mapping), self.upper.inst_pat(mapping),self.body.inst_pat(mapping))
        elif is_evalat(self):
            return EvalAt(self.var, self.lower.inst_pat(mapping), self.upper.inst_pat(mapping),
                          self.body.inst_pat(mapping))
        elif is_deriv(self):
            if self.var in mapping and is_var(mapping[self.var]):
                return Deriv(mapping[self.var].name, self.body.inst_pat(mapping))
            return Deriv(self.var, self.body.inst_pat(mapping))
        elif is_summation(self):
            return Summation(self.index_var, self.lower.inst_pat(mapping), self.upper.inst_pat(mapping), \
                             self.body.inst_pat(mapping))
        elif is_product(self):
            return Product(self.index_var, self.lower.inst_pat(mapping), self.upper.inst_pat(mapping), \
                             self.body.inst_pat(mapping))
        elif is_limit(self):
            return Limit(self.var, self.lim.inst_pat(mapping), self.body.inst_pat(mapping), self.drt)
        else:
            print(type(self))
            raise NotImplementedError

    def has_var(self, var):
        """Check if var occurs in self"""
        assert isinstance(var, Expr) and var.ty == VAR, \
            "%s is not a var" % var
        if self.ty in (VAR, CONST):
            return self == var
        elif self.ty == SKOLEMFUNC:
            return var in self.dependent_vars
        elif self.ty in (OP, FUN):
            return any(subexpr.has_var(var) for subexpr in self.args)
        elif self.ty == DERIV:
            return self.body.has_var(var)
        elif self.ty == INTEGRAL:
            return self.lower.has_var(var) or self.upper.has_var(var) or \
                   self.body.has_var(var)
        elif self.ty == EVAL_AT:
            return self.var != str(var) and (self.body.has_var(var) or \
                                             self.upper.has_var(var) or self.lower.has_var(var))
        else:
            raise NotImplementedError

def exprify(value):
    # judge whether the value is Expr

    if isinstance(value, Expr):
        return value

    if isinstance(value, (int, float)):
        return Const(value)

    if isinstance(value, str):
        return Var(value)
    # 对于其他类型的输入，抛出异常
    raise TypeError(f"无法将类型 {type(value).__name__} 的值 {value} 转换为 Expr")

def contains_i(e: Expr) -> bool:
    """Check if the expression contains the imaginary unit i."""
    if is_const(e):
        return False
    elif is_var(e):
        return False
    elif is_fun(e) and e.func_name == "i":
        return True
    elif is_inf(e):
        return False
    elif is_symbol(e):
        return False
    elif is_op(e):
        return any(contains_i(arg) for arg in e.args)
    elif is_integral(e) or is_deriv(e) or is_limit(e) or is_summation(e) or is_product(e):
        return contains_i(e.body)
    else:
        return False

def is_var(e: Expr) -> TypeGuard["Var"]:
    return e.ty == VAR

def is_const(e: Expr) -> TypeGuard["Const"]:
    return e.ty == CONST

def is_op(e: Expr) -> TypeGuard["Op"]:
    return e.ty == OP

def is_fun(e: Expr, name: str = "") -> TypeGuard["Fun"]:
    if name == "":
        return e.ty == FUN
    else:
        return e.ty == FUN and e.func_name == name

def is_deriv(e: Expr) -> TypeGuard["Deriv"]:
    return e.ty == DERIV

def is_skolem_func(e: Expr) -> TypeGuard["SkolemFunc"]:
    return e.ty == SKOLEMFUNC

def is_symbol(e: Expr) -> TypeGuard["Symbol"]:
    return e.ty == SYMBOL

def is_integral(e: Expr) -> TypeGuard["Integral"]:
    return e.ty == INTEGRAL

def is_cintegral(e: Expr) -> TypeGuard["CIntegral"]:
    return e.ty == COMPOUNDCONTOUR or e.ty == MULTIPOLECONTOUR or e.ty == CINTEGRAL

def is_circlepath(e: Expr) -> TypeGuard["CirclePath"]:
    return e.ty == CIRCLEPATH

def is_rectanglepath(e: Expr) -> TypeGuard["RectanglePath"]:
    return e.ty == RECTANGLEPATH

def is_polepath(e: Expr) -> TypeGuard["PolePath"]:
    return e.ty == POLEPATH

def is_indefinite_integral(e: Expr) -> TypeGuard["IndefiniteIntegral"]:
    return e.ty == INDEFINITEINTEGRAL

def is_evalat(e: Expr) -> TypeGuard["EvalAt"]:
    return e.ty == EVAL_AT

def is_limit(e: Expr) -> TypeGuard["Limit"]:
    return e.ty == LIMIT

def is_summation(e: Expr) -> TypeGuard["Summation"]:
    return e.ty == SUMMATION

def is_product(e: Expr) -> TypeGuard["Product"]:
    return e.ty == PRODUCT

def is_inf(e: Expr) -> TypeGuard["Inf"]:
    return e.ty == INF and (e.t == Decimal("inf") or e.t == Decimal("-inf"))

def is_pos_inf(e: Expr) -> TypeGuard["Inf"]:
    return e.ty == INF and e.t == Decimal("inf")

def is_neg_inf(e: Expr) -> TypeGuard["Inf"]:
    return e.ty == INF and e.t == Decimal("-inf")

def is_plus(e: Expr) -> TypeGuard["Op"]:
    return e.ty == OP and e.op == '+' and len(e.args) == 2

def is_minus(e: Expr) -> TypeGuard["Op"]:
    return e.ty == OP and e.op == '-' and len(e.args) == 2

def is_uminus(e: Expr) -> TypeGuard["Op"]:
    return e.ty == OP and e.op == '-' and len(e.args) == 1

def is_times(e: Expr) -> TypeGuard["Op"]:
    return e.ty == OP and e.op == '*' and len(e.args) == 2

def is_divides(e: Expr) -> TypeGuard["Op"]:
    return e.ty == OP and e.op == '/' and len(e.args) == 2

def is_power(e: Expr) -> TypeGuard["Op"]:
    return e.ty == OP and e.op == '^' and len(e.args) == 2

def is_mod(e: Expr) -> TypeGuard["Op"]:
    return e.ty == OP and e.op == '%' and len(e.args) == 2

def is_less(e: Expr) -> TypeGuard["Op"]:
    return is_op(e) and e.op == '<'

def is_greater(e: Expr) -> TypeGuard["Op"]:
    return is_op(e) and e.op == '>'

def is_less_eq(e: Expr) -> TypeGuard["Op"]:
    return is_op(e) and e.op == '<='

def is_greater_eq(e: Expr) -> TypeGuard["Op"]:
    return is_op(e) and e.op == '>='

def is_equals(e: Expr) -> TypeGuard["Op"]:
    return is_op(e) and e.op == '='

def is_not_equals(e: Expr) -> TypeGuard["Op"]:
    return is_op(e) and e.op == '!='

def is_compare(e: Expr) -> TypeGuard["Op"]:
    return is_op(e) and e.op in ('<', '>', '<=', '>=', '=', '!=')

def match(exp: Expr, pattern: Expr) -> Optional[Dict]:
    """Match expr with given pattern.

    If successful, return a dictionary mapping symbols to expressions.
    Otherwise returns None.

    """
    d = dict()

    def rec(exp: Expr, pattern: Expr, bd_vars: Dict[str, str]):
        if isinstance(pattern, Symbol):
            if pattern.name in d:
                return exp == d[pattern.name]
            # Check exp does not contain bound variables
            for var in exp.get_vars():
                if var in bd_vars.values():
                    return False
            if exp.ty in pattern.pat:
                d[pattern.name] = exp
                return True
            else:
                return False
        if exp.ty != pattern.ty:
            return False
        if is_var(exp):
            return pattern.name == exp.name or \
                (pattern.name in bd_vars and bd_vars[pattern.name] == exp.name)
        elif is_const(exp):
            return pattern.val == exp.val
        elif is_op(exp):
            if exp.op != pattern.op or len(exp.args) != len(pattern.args):
                return False
            for i in range(len(exp.args)):
                if not rec(exp.args[i], pattern.args[i], bd_vars):
                    return False
            return True
        elif is_fun(exp):
            if exp.func_name != pattern.func_name or len(exp.args) != len(pattern.args):
                return False
            for i in range(len(exp.args)):
                if not rec(exp.args[i], pattern.args[i], bd_vars):
                    return False
            return True
        elif is_skolem_func(exp):
            if exp.name != pattern.name or len(exp.dependent_vars) != len(pattern.dependent_vars):
                return False
            for i in range(len(exp.dependent_vars)):
                if not rec(exp.dependent_vars[i], pattern.dependent_vars[i], bd_vars):
                    return False
            return True
        elif is_indefinite_integral(exp):
            # Note this ignores set of skolem arguments
            bd_vars[pattern.var] = exp.var
            res = rec(exp.body, pattern.body, bd_vars)
            del bd_vars[pattern.var]
            return res
        elif is_integral(exp):
            bd_vars[pattern.var] = exp.var
            res1 = rec(exp.upper, pattern.upper, bd_vars)
            res2 = rec(exp.lower, pattern.lower, bd_vars)
            res3 = rec(exp.body, pattern.body, bd_vars)
            del bd_vars[pattern.var]
            return res1 and res2 and res3
        elif is_summation(exp):
            bd_vars[pattern.index_var] = exp.index_var
            res1 = rec(exp.upper, pattern.upper, bd_vars)
            res2 = rec(exp.lower, pattern.lower, bd_vars)
            res3 = rec(exp.body, pattern.body, bd_vars)
            del bd_vars[pattern.index_var]
            return res1 and res2 and res3
        elif is_product(exp):
            bd_vars[pattern.index_var] = exp.index_var
            res1 = rec(exp.upper, pattern.upper, bd_vars)
            res2 = rec(exp.lower, pattern.lower, bd_vars)
            res3 = rec(exp.body, pattern.body, bd_vars)
            del bd_vars[pattern.index_var]
            return res1 and res2 and res3
        elif is_inf(exp):
            return exp.t == pattern.t
        elif is_limit(exp):
            bd_vars[pattern.var] = exp.var
            res1 = rec(exp.body, pattern.body, bd_vars)
            res2 = rec(exp.lim, pattern.lim, bd_vars)
            del bd_vars[pattern.var]
            return res1 and res2
        elif is_deriv(exp):
            # TODO: think more about matching of derivatives
            res1 = pattern.var == exp.var
            res2 = rec(exp.body, pattern.body, bd_vars)
            return res1 and res2
        else:
            # Currently not implemented
            print("Match Failed for type:", type(exp))
            return False

    bd_vars = dict()
    if rec(exp, pattern, bd_vars):
        return d
    else:
        return None

def expr_to_pattern(e: Expr) -> Expr:
    """Convert an expression to pattern."""
    vars = e.get_vars()
    def rec(_e:Expr):
        if _e.ty in (CONST, SYMBOL, VAR, INF, SKOLEMFUNC):
            return _e
        elif _e.ty == OP:
            return Op(_e.op, *[rec(arg) for arg in _e.args])
        elif _e.ty == FUN:
            return Fun(_e.func_name, *[rec(arg) for arg in _e.args])
        elif _e.ty == SUMMATION:
            return Summation(_e.index_var, rec(_e.lower), rec(_e.upper), rec(_e.body))
        elif _e.ty == PRODUCT:
            return Product(_e.index_var, rec(_e.lower), rec(_e.upper), rec(_e.body))
        elif _e.ty == INTEGRAL:
            return Integral(_e.var, rec(_e.lower), rec(_e.upper), rec(_e.body))
        elif _e.ty == EVAL_AT:
            return EvalAt(_e.var, rec(_e.lower), rec(_e.upper), rec(_e.body))
        elif _e.ty == INDEFINITEINTEGRAL:
            return IndefiniteIntegral(_e.var, rec(_e.body), _e.skolem_args)
        elif _e.ty == LIMIT:
            return Limit(_e.var, rec(_e.lim), rec(_e.body), _e.drt)
        elif _e.ty == DERIV:
            return Deriv(_e.var, rec(_e.body))
        else:
            raise NotImplementedError(str(_e))
    e = rec(e)
    for var in vars:
        sym = Symbol(var[0], [VAR, CONST, OP, FUN, INTEGRAL, INF, SYMBOL])
        e = e.subst(var[0], sym)
    return e


def find_pattern(expr, pat, transform=None):
    """Find all subexpr can be matched with the given pattern.

    Return a list of: matched expression, location, mapping of symbols.
    If the transform function is provided, first apply it to the mapping
    of symbols.

    """
    c = []

    def rec(e, pat, cur_loc):
        mapping = match(e, pat)
        if mapping:
            if transform is None:
                c.append((e, cur_loc, mapping))
            else:
                c.append((e, cur_loc, transform(mapping)))
        if e.ty in (OP, FUN):
            for i in range(len(e.args)):
                rec(e.args[i], pat, cur_loc + (i,))
        elif e.ty in (INTEGRAL, DERIV, EVAL_AT):
            rec(e.body, pat, cur_loc + (0,))

    rec(expr, pat, tuple())
    return c


def collect_spec_expr(expr, symb):
    c = [p.args[0] for p, _, _ in find_pattern(expr, symb) if len(p.args) != 0]
    return c


def term_decomposition(e:Expr) -> list[Expr]:
    if e.is_plus():
        a, b = e.args
        return term_decomposition(a) + term_decomposition(b)
    elif e.is_minus():
        a, b = e.args
        return term_decomposition(a) + [(item, -1*sign) for item, sign in term_decomposition(b)]
    elif is_uminus(e):
        return [(item , -1*sign) for item, sign in term_decomposition(e.args[0])]
    else:
        return [(e,1)]

def common_factor_extraction(e:Expr) -> tuple[list[Expr], list[Expr], Expr]:
    term_list = term_decomposition(e)
    num_factors_list = []
    denom_factors_list = []
    sign_list = []
    is_all_negative = True
    for term, sign in term_list:
        a, b = decompose_expr_factor(term)
        num_factors_list.append(a)
        denom_factors_list.append(b)
        sign_list.append(sign)
        if sign == 1:
            is_all_negative = False
    def extraction(factors_list, flag_list):
        res = list()
        for i, factor in enumerate(factors_list[0]):
            flag = True
            for factors in factors_list[1:]:
                if factor not in factors:
                    flag = False
                    break
            if flag:
                flag_list[0][i] = False
                for j, factors in enumerate(factors_list[1:]):
                    for k, factor2 in enumerate(factors):
                        if factor2 == factor:
                            flag_list[j+1][k] = False
                res.append(factor)
        return res
    denom_flag_list = [[True for _ in factors] for factors in denom_factors_list]
    num_flag_list = [[True for _ in factors] for factors in num_factors_list]
    denom_comm_factors = extraction(denom_factors_list, denom_flag_list)
    num_common_factors = extraction(num_factors_list, num_flag_list)

    if is_all_negative:
        num_common_factors.append(Const(-1))
        sign_list = [1 for item in sign_list]

    num_factors_list = [[factor for j, factor in enumerate(num_factors) if num_flag_list[i][j] is True] for i, num_factors in enumerate(num_factors_list)]
    denom_factors_list = [[factor for j, factor in enumerate(denom_factors) if denom_flag_list[i][j] is True] for i, denom_factors in enumerate(denom_factors_list)]
    def prod(es):
        es = list(es)
        if len(es) == 0:
            return Const(1)
        else:
            return functools.reduce(operator.mul, es[1:], es[0])
    for i in range(len(num_factors_list)):
        num = prod(num_factors_list[i])
        denom = prod(denom_factors_list[i])
        if i == 0:
            if sign_list[0] == 1:
                s = num/denom if denom != Const(1) else num
            else:
                s = - (num / denom) if denom != Const(1) else -num
        else:
            if sign_list[i] == 1:
                s = s + num/denom if denom != Const(1) else s + num
            else:
                s = s - num/denom if denom != Const(1) else s - num
    if s != Const(1):
        return num_common_factors + [s], denom_comm_factors
    else:
        return num_common_factors, denom_comm_factors

def decompose_expr_factor2(e:Expr) -> tuple[list[Expr], list[Expr]]:
    if e.is_plus() or e.is_minus():
        return common_factor_extraction(e)
    elif e.is_times():
        a, b = e.args
        n1, d1 = decompose_expr_factor2(a)
        n2, d2 = decompose_expr_factor2(b)
        return n1+n2, d1+d2
    elif e.is_divides():
        a, b = e.args
        n1, d1 = decompose_expr_factor2(a)
        n2, d2 = decompose_expr_factor2(b)
        return n1 + d2, d1 + n2
    else:
        return decompose_expr_factor(e)

def decompose_expr_factor(e) -> tuple[list[Expr], list[Expr]]:
    """Get production factors from expr."""
    num_factors, denom_factors = [], []
    def rec(e: Expr, sign):
        nonlocal num_factors, denom_factors
        if e.is_times():
            rec(e.args[0], sign)
            rec(e.args[1], sign)
        elif is_uminus(e):
            num_factors.append(Const(-1))
            rec(e.args[0], sign)
        elif e.is_divides():
            rec(e.args[0], sign)
            rec(e.args[1], -1 * sign)
        elif is_const(e) and isinstance(e.val, int):
            res_list = []
            if e.val > 0:
                factor_dict = sympy.factorint(e.val)
                for k, v in factor_dict.items():
                    res_list = res_list + [Const(k)] * v
            else:
                factor_dict = sympy.factorint(e.val)
                for k, v in factor_dict.items():
                    res_list = res_list + [Const(k)] * v
            if sign == 1:
                num_factors = num_factors + res_list
            else:
                denom_factors = denom_factors + res_list
        elif sign == 1:
            num_factors.append(e)
        else:
            denom_factors.append(e)
    rec(e, 1)
    return num_factors, denom_factors



class Var(Expr):
    """Variable."""
    def __init__(self, name: str):
        assert isinstance(name, str)
        self.ty = VAR
        self.name = name

    def __hash__(self):
        return hash((VAR, self.name))

    def __eq__(self, other):
        return isinstance(other, Var) and self.name == other.name

    def __str__(self):
        return self.name

    def __repr__(self):
        return "Var(%s)" % self.name

class Const(Expr):
    """Constants."""

    def __init__(self, val: Union[bool, int, Fraction, Decimal]):
        assert isinstance(val, (bool, int, Fraction, Decimal))
        if isinstance(val, Decimal):
            val = Fraction(val)
        self.ty = CONST
        if isinstance(val, Fraction) and val.denominator == 1:
            self.val = val.numerator
        else:
            self.val = val

    def __hash__(self):
        return hash((CONST, self.val))

    def __eq__(self, other):
        return isinstance(other, Const) and self.val == other.val

    def __str__(self):
        return str(self.val)

    def __repr__(self):
        return "Const(%s)" % str(self.val)


class Op(Expr):
    """Operators."""
    def __init__(self, op: str, *args):
        assert isinstance(op, str)
        assert all(isinstance(arg, Expr) for arg in args), op +":"+ str(args)
        if len(args) == 1:
            assert op == "-"
        elif len(args) == 2:
            assert op in ["+", "-", "*", "/", "%", "^", "=", "!=", "<", "<=", ">", ">="]
        else:
            raise NotImplementedError
        self.ty = OP
        self.op = op
        self.args: tuple[Expr, ...] = tuple(args)

    def __hash__(self):
        return hash((OP, self.op, tuple(self.args)))

    def __eq__(self, other):
        return isinstance(other, Op) and self.op == other.op and self.args == other.args

    def __str__(self):
        if len(self.args) == 1:
            a, = self.args
            s = str(a)
            if a.priority() < self.priority():
                s = "(%s)" % s
            return "%s%s" % (self.op, s)
        elif len(self.args) == 2:
            a, b = self.args
            if self.op == '/' and is_const(a) and is_const(b) and isinstance(a.val, int) and isinstance(b.val, int):
                return "%s/%s" % (a.val, b.val)
            s1, s2 = str(a), str(b)
            if a.priority() < op_priority[self.op]:
                s1 = "(%s)" % s1
            if b.priority() <= op_priority[self.op]:
                s2 = "(%s)" % s2
            if a.priority() > op_priority[self.op]:
                if is_uminus(a) and self.op == '^':
                    s1 = "(%s)" % s1
            return "%s %s %s" % (s1, self.op, s2)
        else:
            raise NotImplementedError

    def __repr__(self):
        return "Op(%s,%s)" % (self.op, ",".join(repr(arg) for arg in self.args))

    def contains_var(self, x: str) -> bool:
        """Whether self contains variable x."""
        assert isinstance(x, str)
        return x in self.get_vars()

class Fun(Expr):
    """Functions."""

    def __init__(self, func_name: str, *args: Expr):
        if not isinstance(func_name, str):
            raise AssertionError("Fun:", func_name)
        if not all(isinstance(arg, Expr) for arg in args):
            raise AssertionError("Fun:", args)

        self.ty = FUN
        self.args: tuple[Expr, ...] = tuple(args)
        self.func_name = func_name

    def __hash__(self):
        return hash((FUN, self.func_name, self.args))

    def __eq__(self, other):
        return isinstance(other, Fun) and self.func_name == other.func_name and self.args == other.args

    def __str__(self):
        if len(self.args) > 0:
            return "%s(%s)" % (self.func_name, ",".join(str(arg) for arg in self.args))
        else:
            return self.func_name

    def __repr__(self):
        if len(self.args) > 0:
            return "Fun(%s,%s)" % (self.func_name, ",".join(repr(arg) for arg in self.args))
        else:
            return "Fun(%s)" % self.func_name


class Limit(Expr):
    """Limit expression.

    - var: variable which approaches the limit
    - lim: variable
    - body: expression
    - dir: limit side

    """

    def __init__(self, var: str, lim: Expr, body: Expr, drt=None):
        assert isinstance(var, str) and isinstance(lim, Expr) and isinstance(body, Expr), \
            "Illegal expression: %s %s %s" % (type(var), type(lim), type(body))
        self.ty = LIMIT
        self.var = var
        self.lim = lim
        self.body = body.subst(var, Var(var))
        self.drt = drt

    def alpha_convert(self, new_name: str):
        """Change the variable of limit expression to new_name."""
        assert isinstance(new_name, str), "alpha_convert"
        return Limit(new_name, self.lim, self.body.subst(self.var, Var(new_name)), self.drt)

    def __eq__(self, other):

        if not isinstance(other, Limit):
            return False
        if other.var == self.var:
            return other.drt == self.drt and \
                other.lim == self.lim and \
                other.body == self.body
        else:
            return other.alpha_convert(self.var) == self

    def __hash__(self):
        return hash((LIMIT, self.var, self.lim, self.body, self.drt))

    def __str__(self):
        if self.lim == POS_INF or self.lim == NEG_INF:
            return "LIM {%s -> %s}. %s" % (self.var, self.lim, self.body)
        else:
            return "LIM {%s -> %s %s}. %s" % (
                self.var, self.lim, self.drt if self.drt != None else "", self.body)

    def __repr__(self):
        if self.lim == POS_INF or self.lim == NEG_INF:
            return "Limit(%s, %s, %s)" % (self.var, self.lim, self.body)
        else:
            return "Limit(%s, %s%s, %s)" % (
                self.var, self.lim, "" if self.drt == None else self.drt, self.body)


class Inf(Expr):
    """The infinity."""

    def __init__(self, t):
        assert t in (Decimal("inf"), Decimal("-inf"))
        self.ty = INF
        self.t = t

    def __str__(self):
        if self.t == Decimal("inf"):
            return "oo"
        else:
            return "-oo"

    def __repr__(self):
        return "Inf(%s)" % self.t

    def __hash__(self):
        return hash((INF, self.t))

    def __eq__(self, other):
        return isinstance(other, Inf) and self.t == other.t

    def keys(self):
        return ('ty', 't')

    def __getitem__(self, item):
        return getattr(self, item)


class SkolemFunc(Expr):
    """Skolem variable or function"""
    def __init__(self, name: str, dep_vars: Iterable[Var]):
        self.ty = SKOLEMFUNC
        self.name = name
        self.dependent_vars: tuple[Var] = tuple(dep_vars)

    def __eq__(self, other):
        return isinstance(other, SkolemFunc) and \
            self.dependent_vars == other.dependent_vars and self.name == other.name

    def __str__(self):
        if not self.dependent_vars:
            return "SKOLEM_CONST(" + self.name + ")"
        else:
            return "SKOLEM_FUNC(" + self.name + "(" + ", ".join(str(var) for var in self.dependent_vars) + "))"

    def __hash__(self):
        return hash((self.name, tuple(self.dependent_vars), self.ty))


NEG_INF = Inf(Decimal('-inf'))
POS_INF = Inf(Decimal('inf'))
ZERO = Const(0)
TRUE = Const(True)
FALSE = Const(False)


def sin(e: Expr) -> Expr:
    return Fun("sin", e)

def cos(e: Expr) -> Expr:
    return Fun("cos", e)

def tan(e):
    return Fun("tan", e)

def cot(e):
    return Fun("cot", e)

def sec(e):
    return Fun("sec", e)

def csc(e):
    return Fun("csc", e)

def log(e):
    return Fun("log", e)

def exp(e):
    return Fun("exp", e)

def arcsin(e):
    return Fun("arcsin", e)

def arccos(e):
    return Fun("arccos", e)

def arctan(e):
    return Fun("arctan", e)

def arccot(e):
    return Fun("arccot", e)

def arcsec(e):
    return Fun("arcsec", e)

def arccsc(e):
    return Fun("arccsc", e)

def sqrt(e):
    return Fun("sqrt", e)

def abs(e):
    return Fun("abs", e)

def binom(e1: Expr, e2: Expr) -> Expr:
    """Binomial coefficients"""
    return Fun("binom", e1, e2)

def factorial(e: Expr) -> Expr:
    """Factorial of e"""
    return Fun('factorial', e)


pi = Fun("pi")
E = Fun("exp", Const(1))
G = Fun("G")
i = Fun("i")
int_type = Fun("int")
real_type = Fun("real")
complex_type = Fun("complex")

def Eq(s: Expr, t: Expr) -> Expr:
    return Op("=", s, t)

def isInt(t: Expr) -> Expr:
    return Fun("isInt", t)

def isReal(t: Expr) -> Expr:
    return Fun("isReal", t)

def isEven(t: Expr) -> Expr:
    return Fun("isEven", t)

class Deriv(Expr):
    """Derivative of an expression."""

    def __init__(self, var: str, body: Expr):
        assert isinstance(var, str) and isinstance(body, Expr)
        self.ty = DERIV
        self.var: str = var
        self.body: Expr = body.subst(var, Var(var))

    def __hash__(self):
        return hash((DERIV, self.var, self.body))

    def __eq__(self, other):
        return isinstance(other, Deriv) and self.var == other.var and self.body == other.body

    def __str__(self):
        return "D %s. %s" % (self.var, str(self.body))

    def __repr__(self):
        return "Deriv(%s,%s)" % (self.var, repr(self.body))


class IndefiniteIntegral(Expr):
    """Indefinite integral of an expression."""

    def __init__(self, var: str, body: Expr, skolem_args: Tuple[str]):
        assert isinstance(var, str) and isinstance(body, Expr)
        self.ty = INDEFINITEINTEGRAL
        self.var = var
        self.body = body.subst(var, Var(var))
        self.skolem_args = tuple(skolem_args)

    def __hash__(self):
        return hash((INDEFINITEINTEGRAL, self.var, self.body, self.skolem_args))

    def __eq__(self, other):
        return isinstance(other, IndefiniteIntegral) and self.body == other.alpha_convert(self.var).body and \
            self.skolem_args == other.skolem_args

    def __str__(self):
        if self.skolem_args:
            return "INT %s [%s]. %s" % (self.var, ', '.join(self.skolem_args), self.body)
        else:
            return "INT %s. %s" % (self.var, self.body)

    def __repr__(self):
        return "IndefiniteIntegral(%s,%s,%s)" % (self.var, repr(self.body), self.skolem_args)

    def alpha_convert(self, new_name: str):
        """Change the variable of integration to new_name."""
        assert isinstance(new_name, str), "alpha_convert"
        return IndefiniteIntegral(new_name, self.body.subst(self.var, Var(new_name)), self.skolem_args)


class Integral(Expr):
    """Integral of an expression.
    
    Note equality is with respect to alpha equivalence. The hash function
    is likewise.

    """
    def __init__(self, var: str, lower: Expr, upper: Expr, body: Expr):
        assert isinstance(var, str) and isinstance(lower, Expr) and \
               isinstance(upper, Expr) and isinstance(body, Expr)
        self.ty = INTEGRAL
        self.var = var
        self.lower = lower
        self.upper = upper
        self.body = body.subst(var, Var(var))

    def __hash__(self):
        # Convert to standard bound variable
        return hash((INTEGRAL, self.lower, self.upper, self.body.subst(self.var, Var("_u"))))

    def __eq__(self, other):
        return isinstance(other, Integral) and self.lower == other.lower and self.upper == other.upper and \
               self.body == other.alpha_convert(self.var).body

    def __str__(self):
        return "INT %s:[%s,%s]. %s" % (self.var, str(self.lower), str(self.upper), str(self.body))

    def __repr__(self):
        return "Integral(%s,%s,%s,%s)" % (self.var, repr(self.lower), repr(self.upper), repr(self.body))

    def alpha_convert(self, new_name):
        """Change the variable of integration to new_name."""
        assert isinstance(new_name, str), "alpha_convert"
        return Integral(new_name, self.lower, self.upper, self.body.subst(self.var, Var(new_name)))


class EvalAt(Expr):
    """Evaluation at upper and lower, then subtract."""

    def __init__(self, var: str, lower: Expr, upper: Expr, body: Expr):
        assert isinstance(var, str) and isinstance(lower, Expr) and \
               isinstance(upper, Expr) and isinstance(body, Expr)
        self.ty = EVAL_AT
        self.var = var
        self.lower = lower
        self.upper = upper
        self.body = body.subst(var, Var(var))

    def __hash__(self):
        return hash((EVAL_AT, self.var, self.lower, self.upper, self.body))

    def __eq__(self, other):
        return isinstance(other, EvalAt) and self.var == other.var and \
               self.lower == other.lower and self.upper == other.upper and self.body == other.body

    def __str__(self):
        return "[%s]_%s=%s,%s" % (str(self.body), self.var, str(self.lower), str(self.upper))

    def __repr__(self):
        return "EvalAt(%s,%s,%s,%s)" % (self.var, repr(self.lower), repr(self.upper), repr(self.body))


class Symbol(Expr):
    """Pattern expression.
    
    It can be used to find expression with the given specific structure.
    
    """
    def __init__(self, name: str, pat: List[str]):
        self.ty = SYMBOL
        self.name = name
        self.pat = tuple(pat)

    def __eq__(self, other):
        return isinstance(other, Symbol) and self.name == other.name and self.pat == other.pat

    def __hash__(self):
        return hash((SYMBOL, self.name, self.ty, sum(self.pat)))

    def __str__(self):
        return "?%s" % self.name

    def __repr__(self):
        return "Symbol(%s, %s)" % (self.name, self.pat)


class Summation(Expr):
    """Summation of integers over some range."""
    def __init__(self, index_var: str, lower: Expr, upper: Expr, body: Expr):
        self.ty = SUMMATION
        self.index_var: str = index_var
        self.lower: Expr = lower
        self.upper: Expr = upper
        self.body: Expr = body.subst(index_var, Var(index_var))

    def __str__(self):
        return "SUM(" + self.index_var + ", " + str(self.lower) + ", " + str(self.upper) + ", " + str(self.body) + ")"

    def __eq__(self, other):

        if isinstance(other, Summation):
            if self.index_var == other.index_var:
                return self.lower == other.lower and \
                self.upper == other.upper and \
                self.body == other.body
            else:
                return other.alpha_convert(self.index_var) == self
        return False

    def __hash__(self):
        return hash((SUMMATION, self.index_var, self.ty, self.lower, self.upper, self.body))

    def alpha_convert(self, new_var: str):
        """Rename the bound variable of a summation."""
        assert isinstance(new_var, str), "alpha_convert"
        return Summation(new_var, self.lower, self.upper, self.body.subst(self.index_var, Var(new_var)))


class Product(Expr):
    """Summation of integers over some range."""
    def __init__(self, index_var: str, lower: Expr, upper: Expr, body: Expr):
        self.ty = PRODUCT
        self.index_var: str = index_var
        self.lower: Expr = lower
        self.upper: Expr = upper
        self.body: Expr = body.subst(index_var, Var(index_var))

    def __str__(self):
        return "MUL(" + self.index_var + ", " + str(self.lower) + ", " + str(self.upper) + ", " + str(self.body) + ")"

    def __eq__(self, other):
        if isinstance(other, Product):
            if self.index_var == other.index_var:
                return self.lower == other.lower and \
                self.upper == other.upper and \
                self.body == other.body
            else:
                return other.alpha_convert(self.index_var) == self
        return False

    def __hash__(self):
        return hash((PRODUCT, self.index_var, self.ty, self.lower, self.upper, self.body))

    def alpha_convert(self, new_var: str):
        """Rename the bound variable of a product."""
        assert isinstance(new_var, str), "alpha_convert"
        return Product(new_var, self.lower, self.upper, self.body.subst(self.index_var, Var(new_var)))


def eval_expr(e: Expr):
    if is_inf(e):
        if e == POS_INF:
            return float('inf')
        else:
            return float('-inf')
    elif is_const(e):
        return e.val
    elif is_plus(e):
        return eval_expr(e.args[0]) + eval_expr(e.args[1])
    elif is_uminus(e):
        return -eval_expr(e.args[0])
    elif is_minus(e):
        return eval_expr(e.args[0]) - eval_expr(e.args[1])
    elif is_times(e):
        return eval_expr(e.args[0]) * eval_expr(e.args[1])
    elif is_divides(e):
        return eval_expr(e.args[0]) / eval_expr(e.args[1])
    elif is_mod(e):
        return eval_expr(e.args[0]) % eval_expr(e.args[1])
    elif is_power(e):
        return eval_expr(e.args[0]) ** eval_expr(e.args[1])
    elif is_fun(e):
        if e.func_name == 'sqrt':
            return math.sqrt(eval_expr(e.args[0]))
        elif e.func_name == 'exp':
            return math.exp(eval_expr(e.args[0]))
        elif e.func_name == 'i':
            return 1j  # return imaginary unit in Python
        elif e.func_name == 'abs':
            return abs(eval_expr(e.args[0]))
        elif e.func_name == 'pi':
            return math.pi
        elif e.func_name == 'sin':
            return math.sin(eval_expr(e.args[0]))
        elif e.func_name == 'cos':
            return math.cos(eval_expr(e.args[0]))
        elif e.func_name == 'tan':
            return math.tan(eval_expr(e.args[0]))
        elif e.func_name == 'cot':
            return 1.0 / math.tan(eval_expr(e.args[0]))
        elif e.func_name == 'sec':
            return 1.0 / math.cos(eval_expr(e.args[0]))
        elif e.func_name == 'csc':
            return 1.0 / math.sin(eval_expr(e.args[0]))
        elif e.func_name == 'arcsin':
            return math.asin(eval_expr(e.args[0]))
        elif e.func_name == 'arccos':
            return math.acos(eval_expr(e.args[0]))
        elif e.func_name == 'arctan':
            return math.atan(eval_expr(e.args[0]))
        elif e.func_name == 'log':
            a = eval_expr(e.args[0])
            if a <= 0.0:
                return -math.inf
            else:
                return math.log(a)
        elif e.func_name == 'factorial':
            arg = eval_expr(e.args[0])
            if int(arg) == arg:
                return math.factorial(arg)
            else:
                from scipy.special import gamma
                return gamma(float(arg) + 1)
        elif e.func_name == 'Gamma':
            arg = eval_expr(e.args[0])
            from scipy.special import gamma
            return gamma(float(arg))

    raise NotImplementedError(f"eval_expr on {e}")

def neg_expr(e: Expr):
    """Return the negation of the given formula."""
    if is_op(e):
        if e.op == "=":
            return Op("!=", e.args[0], e.args[1])
        elif e.op == "!=":
            return Op("=", e.args[0], e.args[1])
        elif e.op == ">":
            return Op("<=", e.args[0], e.args[1])
        elif e.op == "<":
            return Op(">=", e.args[0], e.args[1])
        elif e.op == ">=":
            return Op("<", e.args[0], e.args[1])
        elif e.op == "<=":
            return Op(">", e.args[0], e.args[1])
        else:
            raise NotImplementedError(f"neg_expr: {e}")
    else:
        raise NotImplementedError(f"neg_expr: {e}")

class MultiPoleContourIntegral(Expr):
    """Multi-pole contour integral of an expression.
    
    The contour consists of small circles around multiple poles.
    """
    def __init__(self, var: str, poles: List[Expr], radii: List[Expr], body: Expr):
        assert isinstance(var, str) and isinstance(poles, list) and \
               isinstance(radii, list) and isinstance(body, Expr)
        assert len(poles) == len(radii)
        self.ty = MULTIPOLECONTOUR
        self.var = var
        self.poles = poles  # List of pole points
        self.radii = radii  # List of radii
        self.body = body.subst(var, Var(var))

    def __hash__(self):
        return hash((MULTIPOLECONTOUR, self.var, tuple(self.poles), tuple(self.radii),
                    self.body.subst(self.var, Var("_u"))))

    def __eq__(self, other):
        return isinstance(other, MultiPoleContourIntegral) and \
               self.poles == other.poles and self.radii == other.radii and \
               self.body == other.alpha_convert(self.var).body

    def __str__(self):
        pole_str = ",".join("%s,%s" % (str(p), str(r)) for p, r in zip(self.poles, self.radii))
        return "CINT %s:poles(%s). %s" % (self.var, pole_str, str(self.body))

    def __repr__(self):
        return "MultiPoleContourIntegral(%s,%s,%s,%s)" % (
            self.var, repr(self.poles), repr(self.radii), repr(self.body))

    def alpha_convert(self, new_name):
        """Change the variable of integration to new_name."""
        assert isinstance(new_name, str), "alpha_convert"
        return MultiPoleContourIntegral(new_name, self.poles, self.radii,
                                      self.body.subst(self.var, Var(new_name)))

class CirclePath:
    """Circle path for compound contour integral."""
    def __init__(self, center: Expr, end_r: Expr, begin_a: Expr, end_a: Expr, direction: str):
        assert isinstance(center, Expr) and isinstance(end_r, Expr) and \
               isinstance(begin_a, Expr) and isinstance(end_a, Expr) and isinstance(direction, str)
        assert direction in ["cw", "ccw"]
        self.ty = CIRCLEPATH
        # 处理圆心的复数形式
        if center.is_plus() and len(center.args) == 2 and \
           center.args[1].is_times() and center.args[1].args[1] == i:
            self.center = center
        else:
            # 如果不是复数形式，尝试转换
            self.center = center
        self.end_r = end_r
        self.begin_a = begin_a
        self.end_a = end_a
        self.direction = direction

    def __hash__(self):
        return hash((self.center, self.end_r, self.begin_a, self.end_a, self.direction))

    def get_vars(self) -> Set[str]:
        """获取路径中的所有变量"""
        vars = set()
        vars.update(self.center.get_vars())
        vars.update(self.end_r.get_vars())
        vars.update(self.begin_a.get_vars())
        vars.update(self.end_a.get_vars())
        return vars

    def is_inside(self, point: Expr, ctx) -> bool:
        """判断一个点是否在圆形路径内部。
        
        Args:
            point: 要判断的点
            ctx: 上下文环境
        
        Returns:
            bool: 点在圆内返回True，否则返回False
        """
        from integral.poly import normalize
        
        # 如果点就是圆心，点在圆内
        if point == self.center:
            return True
            
        # 计算点到圆心的距离
        distance = Fun("abs", Op("-", point, self.center))
        distance = normalize(distance, ctx)
        radius = normalize(self.end_r, ctx)
        
        # 如果距离和半径都可以直接求值
        if distance.is_evaluable() and radius.is_evaluable():
            return distance < radius
        # 否则进行符号比较
        comparison = Op("<", distance, radius)
        return ctx.check_condition(comparison)

    def __eq__(self, other):
        return isinstance(other, CirclePath) and \
               self.center == other.center and self.end_r == other.end_r and \
               self.begin_a == other.begin_a and self.end_a == other.end_a and \
               self.direction == other.direction

    def __str__(self):
        # 正确处理圆心的复数形式
        if self.center.is_plus() and len(self.center.args) == 2 and \
           self.center.args[1].is_times() and self.center.args[1].args[1] == i:
            real_part = self.center.args[0]
            imag_part = self.center.args[1].args[0]
            return "circle((%s,%s),%s,%s,%s,%s)" % (
                real_part, imag_part, self.end_r, self.begin_a, self.end_a, self.direction)
        else:
            return "circle((%s,0),%s,%s,%s,%s)" % (
                self.center, self.end_r, self.begin_a, self.end_a, self.direction)

    def __repr__(self):
        return "CirclePath(%s,%s,%s,%s,%s)" % (
            repr(self.center), repr(self.end_r), repr(self.begin_a), repr(self.end_a), repr(self.direction))
    
    def size(self):
        return 1 + self.center.size() + self.end_r.size() + self.begin_a.size() + self.end_a.size()

class LinePath:
    """Line path for compound contour integral."""
    def __init__(self, start: Expr, end: Expr):
        assert isinstance(start, Expr) and isinstance(end, Expr)
        self.ty = LINEPATH
        self.start = start
        self.end = end

    def get_vars(self) -> Set[str]:
        """获取路径中的所有变量"""
        vars = set()
        vars.update(self.start.get_vars())
        vars.update(self.end.get_vars())
        return vars
    
    def __hash__(self):
        return hash((self.start, self.end))

    def __eq__(self, other):
        return isinstance(other, LinePath) and \
               self.start == other.start and self.end == other.end

    def __str__(self):
        return "line(%s,%s)" % (str(self.start), str(self.end))

    def __repr__(self):
        return "LinePath(%s,%s)" % (repr(self.start), repr(self.end))
    
    def size(self):
        return 1 + self.start.size() + self.end.size()

class PolePath:
    """Pole path for compound contour integral."""
    def __init__(self, Re: Expr, Im: Expr):
        assert isinstance(Re, Expr) and isinstance(Im, Expr)
        self.ty = POLEPATH
        self.Re = Re
        self.Im = Im
        self.point = Op("+", Re, Op("*", Im, Fun("i")))

    def get_vars(self) -> Set[str]:
        """获取路径中的所有变量"""
        vars = set()
        vars.update(self.Re.get_vars())
        vars.update(self.Im.get_vars())
        return vars

    def __eq__(self, other):
        return isinstance(other, PolePath) and \
               self.Re == other.Re and self.Im == other.Im
    
    def __hash__(self):
        return hash((self.Re, self.Im))

    def __str__(self):
        return "pole(%s,%s)" % (str(self.Re), str(self.Im))

    def __repr__(self):
        return "PolePath(%s,%s)" % (repr(self.Re), repr(self.Im))
    
    def size(self):
        return 1 + self.Re.size() + self.Im.size()

class RectanglePath:
    """Rectangle path for compound contour integral."""
    def __init__(self, Re1: Expr, Im1: Expr, Re2: Expr, Im2: Expr, Re3: Expr, Im3: Expr, Re4: Expr, Im4: Expr):
        assert isinstance(Re1, Expr) and isinstance(Im1, Expr) and \
               isinstance(Re2, Expr) and isinstance(Im2, Expr) and \
               isinstance(Re3, Expr) and isinstance(Im3, Expr) and \
               isinstance(Re4, Expr) and isinstance(Im4, Expr)
        self.ty = RECTANGLEPATH
        self.Re1 = Re1
        self.Im1 = Im1
        self.Re2 = Re2
        self.Im2 = Im2
        self.Re3 = Re3
        self.Im3 = Im3
        self.Re4 = Re4
        self.Im4 = Im4

    def get_vars(self) -> Set[str]:
        """获取路径中的所有变量"""
        vars = set()
        vars.update(self.Re1.get_vars())
        vars.update(self.Im1.get_vars())
        vars.update(self.Re2.get_vars())
        vars.update(self.Im2.get_vars())
        vars.update(self.Re3.get_vars())
        vars.update(self.Im3.get_vars())
        vars.update(self.Re4.get_vars())
        vars.update(self.Im4.get_vars())
        return vars

    def __eq__(self, other):
        return isinstance(other, RectanglePath) and \
               self.Re1 == other.Re1 and self.Im1 == other.Im1 and \
               self.Re2 == other.Re2 and self.Im2 == other.Im2 and \
               self.Re3 == other.Re3 and self.Im3 == other.Im3 and \
               self.Re4 == other.Re4 and self.Im4 == other.Im4
    
    def __hash__(self):
        return hash((self.Re1, self.Im1, self.Re2, self.Im2, self.Re3, self.Im3, self.Re4, self.Im4))

    def __str__(self):
        return "rectangle((%s,%s),(%s,%s),(%s,%s),(%s,%s))" % (
            str(self.Re1), str(self.Im1), str(self.Re2), str(self.Im2), str(self.Re3), str(self.Im3), str(self.Re4), str(self.Im4))

    def __repr__(self):
        return "RectanglePath((%s,%s),(%s,%s),(%s,%s),(%s,%s))" % (
            repr(self.Re1), repr(self.Im1), repr(self.Re2), repr(self.Im2), repr(self.Re3), repr(self.Im3), repr(self.Re4), repr(self.Im4))

    def size(self):
        return 1 + self.Re1.size() + self.Im1.size() + self.Re2.size() + self.Im2.size() + self.Re3.size() + self.Im3.size() + self.Re4.size() + self.Im4.size()

class CompoundContourIntegral(Expr):
    """Compound contour integral of an expression.
    
    The contour is composed of multiple paths.
    """
    def __init__(self, var: str, paths: List[Union[CirclePath, PolePath, RectanglePath]], body: Expr):
        assert isinstance(var, str) and isinstance(body, Expr)
        self.ty = COMPOUNDCONTOUR
        self.var = var
        # 确保paths是列表
        self.paths = [paths] if not isinstance(paths, list) else paths
        self.body = body.subst(var, Var(var))

    def __hash__(self):
        # 为每个路径创建一个唯一的哈希值
        path_hashes = []
        for p in self.paths:
            if isinstance(p, CirclePath):
                path_hashes.append(hash((CIRCLEPATH, str(p.center), str(p.end_r), str(p.begin_a), str(p.end_a), p.direction)))
            elif isinstance(p, LinePath):
                path_hashes.append(hash((LINEPATH, str(p.start), str(p.end))))
            elif isinstance(p, PolePath):
                path_hashes.append(hash((POLEPATH, str(p.point), str(p.radius))))
            elif isinstance(p, RectanglePath):
                path_hashes.append(hash((RECTANGLEPATH, str(p.z1), str(p.z2), str(p.z3), str(p.z4))))
        return hash((COMPOUNDCONTOUR, self.var, tuple(path_hashes),
                    self.body.subst(self.var, Var("_u"))))

    def __eq__(self, other):
        return isinstance(other, CompoundContourIntegral) and \
               self.paths == other.paths and \
               self.body == other.alpha_convert(self.var).body

    def __str__(self):
        paths_str = ",".join(str(p) for p in self.paths)
        return "CINT %s:com(%s). %s" % (self.var, paths_str, str(self.body))

    def __repr__(self):
        return "CompoundContourIntegral(%s,%s,%s)" % (
            self.var, repr(self.paths), repr(self.body))

    def alpha_convert(self, new_name):
        """Change the variable of integration to new_name."""
        assert isinstance(new_name, str), "alpha_convert"
        return CompoundContourIntegral(new_name, self.paths,
                                     self.body.subst(self.var, Var(new_name)))

class ComplexNumber:
    """Complex number representation."""
    def __init__(self, real, imag):
        self.real = real
        self.imag = imag

    def __str__(self):
        if self.imag == 0:
            return str(self.real)
        if self.real == 0:
            return str(self.imag) + "*"+"i"
        if self.imag < 0:
            return f"{self.real}{self.imag}*i"
        return f"{self.real}+{self.imag}*i"

    def __add__(self, other):
        return ComplexNumber(self.real + other.real, self.imag + other.imag)

    def __mul__(self, other):
        return ComplexNumber(
            self.real * other.real - self.imag * other.imag,
            self.real * other.imag + self.imag * other.real
        )

    def conjugate(self):
        return ComplexNumber(self.real, -self.imag)

def find_poles(var, expr, ctx=None):
    """
    查找复函数的极点。

    Args:
        var: 要查找极点的变量
        expr: 要查找极点的表达式
        ctx: 上下文环境(可选)
        
    Returns:
        极点列表
    """
    from integral import solve
    
    if ctx is None:
        from integral.context import Context
        ctx = Context()

    poles = []

    # 处理分式表达式
    if is_op(expr) and expr.op == '/':
        # 获取分母的零点
        denominator = expr.args[1]
        if denominator.contains_var(var):
            try:
                # 尝试解方程找到极点
                zeros = solve.solve_equation(denominator, Const(0), var, ctx)
                if isinstance(zeros, list):
                    poles.extend(zeros)
                else:
                    poles.append(zeros)
            except:
                pass
    
    # 处理特殊情况: 1/(z^2+1)
    if is_op(expr) and expr.op == '/' and len(expr.args) == 2:
        num, denom = expr.args
        if is_const(num) and num.val == 1:
            if denom.is_plus() and len(denom.args) == 2:
                # 检查是否是 z^2 + 1 的形式
                if (denom.args[0].is_power() and 
                    is_var(denom.args[0].args[0]) and
                    is_const(denom.args[0].args[1]) and
                    denom.args[0].args[1].val == 2 and
                    is_const(denom.args[1]) and
                    denom.args[1].val == 1):
                    
                    # 返回极点 i 和 -i
                    poles.extend([Fun("i"), Op("*", Const(-1), Fun("i"))])
    
    # 处理三角函数
    if is_fun(expr):
        if expr.func_name == "tan":
            # tan(z)在z = (n + 1/2)π处有极点
            poles.append(Op("*", Fun("pi"), Const(1)/Const(2)))
        elif expr.func_name == "cot":
            # cot(z)在z = nπ处有极点
            poles.append(Fun("pi"))
        elif expr.func_name == "csc":
            # csc(z)在z = nπ处有极点
            poles.append(Fun("pi"))
        elif expr.func_name == "sec":
            # sec(z)在z = (n + 1/2)π处有极点
            poles.append(Op("*", Fun("pi"), Const(1)/Const(2)))
    
    # 递归处理复合表达式
    if is_op(expr):
        if expr.op in ['+', '-', '*']:
            for arg in expr.args:
                poles.extend(find_poles(var, arg, ctx))
    
    # 去重
    unique_poles = []
    for pole in poles:
        if pole not in unique_poles:
            unique_poles.append(pole)
    
    return unique_poles

def compute_residue(expr, pole, order=1):
    """Compute residue at a pole.
    
    Args:
        expr: The expression to compute residue for.
        pole: The pole to compute residue at.
        order: Order of the pole (default is 1 for simple poles).
        
    Returns:
        The residue at the pole.
    """
    if order == 1:
        # 对于一阶极点，使用极限公式
        # Res(f,a) = lim(z->a) (z-a)f(z)
        z = Var('z')
        # 确保表达式中的变量被正确替换为z
        if isinstance(expr, Expr):
            vars = expr.get_vars()
            if len(vars) == 1:
                var = list(vars)[0]
                expr = expr.subst(var, z)
        residue = Limit('z', pole, (z - pole) * expr)
        return residue
    else:
        # 对于高阶极点，使用导数公式
        # Res(f,a) = 1/(n-1)! * lim(z->a) d^(n-1)/dz^(n-1) [(z-a)^n * f(z)]
        z = Var('z')
        # 确保表达式中的变量被正确替换为z
        if isinstance(expr, Expr):
            vars = expr.get_vars()
            if len(vars) == 1:
                var = list(vars)[0]
                expr = expr.subst(var, z)
        n = order
        expr_mult = ((z - pole) ** n) * expr
        
        # 计算n-1阶导数
        for _ in range(n-1):
            expr_mult = Deriv('z', expr_mult)
            
        # 计算极限
        residue = Limit('z', pole, expr_mult) / factorial(n-1)
        return residue

def factorial(n):
    """计算阶乘。"""
    if n <= 0:
        return Const(1)
    result = 1
    for i in range(1, n+1):
        result *= i
    return Const(result)

def is_inside_contour(point: Expr, contour: Union[MultiPoleContourIntegral, MultiPoleContourIntegral]) -> bool:
    """判断点是否在围道内。
    
    使用绕数(winding number)方法判断点是否在围道内。
    对于简单闭合曲线,如果绕数不为0,则点在曲线内部。
    
    Args:
        point: 要判断的点
        contour: 围道积分
        
    Returns:
        bool: 点是否在围道内
    """
    if isinstance(contour, MultiPoleContourIntegral):
        # 对于多极点围道，判断点是否在任一小圆内
        for pole, radius in zip(contour.poles, contour.radii):
            if abs(point - pole) < radius:
                return True
        return False
        
    elif isinstance(contour, InfinityContourIntegral):
        # 对于无穷远点围道，判断点是否在大圆内
        return abs(point) < contour.radius
        
    return False

def is_complex_analytic(e: Expr) -> bool:
    """判断表达式是否为解析函数。"""
    if is_const(e) or is_var(e):
        return True
        
    elif e.is_plus() or e.is_minus() or e.is_times():
        return all(is_complex_analytic(arg) for arg in e.args)
        
    elif e.is_divides():
        # 分母不为0时解析
        return is_complex_analytic(e.args[0]) and is_complex_analytic(e.args[1])
        
    elif e.is_power():
        base, exp = e.args
        if is_const(exp) and isinstance(exp.val, int) and exp.val >= 0:
            return is_complex_analytic(base)
        return False
        
    elif is_fun(e):
        if e.func_name in ["sin", "cos", "tan", "exp", "log"]:
            return all(is_complex_analytic(arg) for arg in e.args)
        return False
        
    return False

def get_branch_points(e: Expr) -> List[Expr]:
    """获取复变函数的分支点。"""
    from integral import solve
    points = []
    
    if e.is_divides():
        # 分母的零点是极点
        zeros = solve.solve_equation(e.args[1], Const(0))
        points.extend(zeros)
        
    elif e.is_power():
        base, exp = e.args
        if not (is_const(exp) and isinstance(exp.val, int) and exp.val >= 0):
            # 非整数幂的底数零点是分支点
            zeros = solve.solve_equation(base, Const(0))
            points.extend(zeros)
            
    elif is_fun(e):
        if e.func_name == "log":
            # 对数函数在0处有分支点
            points.append(Const(0))
            
        elif e.func_name == "sqrt":
            # 平方根在0处有分支点
            points.append(Const(0))
            
        elif e.func_name in ["arcsin", "arccos"]:
            # 反三角函数在±1处有分支点
            points.extend([Const(1), Const(-1)])
            
    return points

def get_singular_points(e: Expr) -> List[Tuple[Expr, str]]:
    """获取复变函数的奇点。
    
    返回值是一个列表，每个元素是(point, type)的元组，
    其中type可以是:
    - "pole": 极点
    - "essential": 本性奇点
    - "branch": 分支点
    """
    from integral import solve
    points = []
    
    if e.is_divides():
        # 分母的零点是极点
        zeros = solve.solve_equation(e.args[1], Const(0))
        for zero in zeros:
            points.append((zero, "pole"))
            
    elif is_fun(e):
        if e.func_name == "tan":
            # 正切函数在kπ处有极点
            points.append((Fun("pi"), "pole"))
            
        elif e.func_name == "log":
            # 对数在0处有分支点
            points.append((Const(0), "branch"))
            
        elif e.func_name == "exp":
            # 指数函数在无穷远处有本性奇点
            points.append((POS_INF, "essential"))
            
    # 递归处理子表达式
    if e.is_plus() or e.is_minus() or e.is_times() or e.is_divides():
        for arg in e.args:
            points.extend(get_singular_points(arg))
            
    return points

class CIntegral(Expr):
    """Contour integral of an expression."""
    def __init__(self, var: str, paths: List[Union[CirclePath, PolePath, RectanglePath]], body: Expr):
        assert isinstance(var, str) and isinstance(body, Expr)
        self.ty = CINTEGRAL
        self.var = var
        self.paths = paths
        self.body = body.subst(var, Var(var))
    
    def __hash__(self):
        # 确保paths是一个元组，并且每个path都是可哈希的
        path_tuple = tuple(hash(p) for p in self.paths)
        return hash((CINTEGRAL, self.var, path_tuple, self.body))
    
    def __eq__(self, other):
        return isinstance(other, CIntegral) and self.var == other.var and self.paths == other.paths and self.body == other.body

    def __str__(self):
        return "CINT %s:com(%s). %s" % (self.var, ",".join(str(p) for p in self.paths), str(self.body))
    
    def alpha_convert(self, new_name):
        """Change the variable of integration to new_name."""
        assert isinstance(new_name, str), "alpha_convert: new_name must be a string"
        return CIntegral(new_name, self.paths, self.body.subst(self.var, Var(new_name)))

def compute_residue(f: Expr, pole: Expr, order: int = 1) -> Expr:
    """计算函数f在极点pole处的留数
    
    对于简单极点（order=1），留数计算公式为：
    Res(f, pole) = lim_{z->pole} (z-pole) * f(z)
    
    对于高阶极点，留数计算公式为：
    Res(f, pole) = (1/(order-1)!) * lim_{z->pole} (d^(order-1)/dz^(order-1)) [(z-pole)^order * f(z)]
    
    Args:
        f: 复变函数表达式
        pole: 极点
        order: 极点的阶数，默认为1（简单极点）
    
    Returns:
        Expr: 留数表达式
    """
    var = "z"  # 假设变量是z
    z = Var(var)
    
    if order == 1:  # 简单极点
        # 留数 = lim_{z->pole} (z-pole) * f(z)
        expr_to_limit = Op("*", Op("-", z, pole), f)
        return Limit(var, pole, expr_to_limit)
    else:  # 高阶极点
        # 对于高阶极点，需要进行复杂的计算
        # (z-pole)^order * f(z)
        expr_to_diff = Op("*", Op("^", Op("-", z, pole), Const(order)), f)
        
        # 对上面的表达式求(order-1)阶导数
        expr_derivative = expr_to_diff
        for _ in range(order - 1):
            expr_derivative = Deriv(var, expr_derivative)
        
        # 计算常数系数 1/(order-1)!
        factorial = 1
        for i in range(1, order):
            factorial *= i
        coef = Op("/", Const(1), Const(factorial))
        
        # 计算留数表达式
        return Op("*", coef, Limit(var, pole, expr_derivative))

def is_polynomial(expr: Expr, var: str) -> bool:
    """判断表达式是否为关于var的多项式"""
    if is_const(expr):
        return True
    elif is_var(expr) and expr.name == var:
        return True
    elif expr.is_plus() or expr.is_minus():
        return all(is_polynomial(arg, var) for arg in expr.args)
    elif expr.is_times():
        return all(is_polynomial(arg, var) for arg in expr.args)
    elif expr.is_power():
        base, exp = expr.args
        return is_var(base) and base.name == var and is_const(exp) and isinstance(exp.val, int) and exp.val >= 0
    return False

def solve_polynomial_zeros(polynomial: Expr, var: str) -> List[Expr]:
    """求解多项式的零点"""
    # 这是一个简化版的求解函数，实际应用中可能需要更复杂的实现
    # 对于一阶多项式 ax + b = 0
    if is_polynomial(polynomial, var) and max_degree(polynomial, var) == 1:
        # 获取系数
        a = coefficient(polynomial, var, 1)
        b = coefficient(polynomial, var, 0)
        # 求解 x = -b/a
        return [Op("/", Op("-", b), a)]
    
    # 对于二阶多项式 ax^2 + bx + c = 0
    elif is_polynomial(polynomial, var) and max_degree(polynomial, var) == 2:
        # 获取系数
        a = coefficient(polynomial, var, 2)
        b = coefficient(polynomial, var, 1)
        c = coefficient(polynomial, var, 0)
        
        # 计算判别式
        delta = Op("-", Op("^", b, Const(2)), Op("*", Const(4), Op("*", a, c)))
        
        # 求解 x = (-b ± √Δ) / (2a)
        x1 = Op("/", Op("+", Op("-", b), Fun("sqrt", delta)), Op("*", Const(2), a))
        x2 = Op("/", Op("-", Op("-", b), Fun("sqrt", delta)), Op("*", Const(2), a))
        
        return [x1, x2]
    
    # 对于高阶多项式，返回空列表
    return []

def max_degree(polynomial: Expr, var: str) -> int:
    """计算多项式关于var的最高次数"""
    if is_const(polynomial):
        return 0
    elif is_var(polynomial) and polynomial.name == var:
        return 1
    elif polynomial.is_plus() or polynomial.is_minus():
        return max(max_degree(arg, var) for arg in polynomial.args)
    elif polynomial.is_times():
        return sum(max_degree(arg, var) for arg in polynomial.args)
    elif polynomial.is_power() and is_var(polynomial.args[0]) and polynomial.args[0].name == var:
        if is_const(polynomial.args[1]):
            return polynomial.args[1].val
    return 0

def coefficient(polynomial: Expr, var: str, degree: int) -> Expr:
    """计算多项式中var^degree项的系数"""
    # 这是一个简化版的系数提取函数，实际应用中可能需要更复杂的实现
    if is_const(polynomial):
        return polynomial if degree == 0 else Const(0)
    elif is_var(polynomial) and polynomial.name == var:
        return Const(1) if degree == 1 else Const(0)
    elif polynomial.is_plus():
        return Op("+", coefficient(polynomial.args[0], var, degree), 
                     coefficient(polynomial.args[1], var, degree))
    elif polynomial.is_minus():
        return Op("-", coefficient(polynomial.args[0], var, degree), 
                     coefficient(polynomial.args[1], var, degree))
    elif polynomial.is_times():
        # 检查是否有var^degree因子
        var_factor = None
        const_factors = []
        
        for arg in polynomial.args:
            if is_var(arg) and arg.name == var and degree == 1:
                var_factor = arg
            elif arg.is_power() and is_var(arg.args[0]) and arg.args[0].name == var:
                if is_const(arg.args[1]) and arg.args[1].val == degree:
                    var_factor = arg
            else:
                const_factors.append(arg)
        
        if var_factor is not None:
            # 返回其他因子的乘积作为系数
            if not const_factors:
                return Const(1)
            elif len(const_factors) == 1:
                return const_factors[0]
            else:
                result = const_factors[0]
                for factor in const_factors[1:]:
                    result = Op("*", result, factor)
                return result
    
    return Const(0)
