"""Conditions"""

from typing import Union, Iterable

from integral.expr import Expr
from integral import latex


class Conditions:
    """A condition is represented by a list of boolean expressions."""
    def __init__(self, conds: Union["Conditions", Iterable[Union[str, Expr]]] = None):
        from integral import parser
        self.data: list[Expr] = list()
        if conds is None:
            pass
        elif isinstance(conds, Conditions):
            self.data.extend(conds.data)
        else:
            parsed_conds = []
            for cond in list(conds):
                if isinstance(cond, str):
                    parsed_conds.append(parser.parse_condition(cond))
                else:
                    assert isinstance(cond, Expr)
                    parsed_conds.append(cond)
            self.data.extend(parsed_conds)

    def __bool__(self):
        return bool(self.data)

    def __hash__(self):
        return hash(tuple(self.data))

    def __str__(self):
        return ", ".join(str(cond) for cond in self.data)

    def add_condition(self, cond: Expr):
        assert isinstance(cond, Expr)
        if cond not in self.data:
            self.data.append(cond)

    def __eq__(self, other):
        return isinstance(other, Conditions) and self.data == other.data

    def export(self):
        res = list()
        for cond in self.data:
            res.append({
                "cond": str(cond),
                "latex_cond": latex.convert_expr(cond)
            })
        return res

    def update(self, other: 'Conditions'):
        for e in other.data:
            self.add_condition(e)
