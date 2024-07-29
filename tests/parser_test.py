"""Unit test for parsing."""

import unittest
from fractions import Fraction

from integral.expr import Var, Const, Op, Fun
from integral.parser import parse_expr, parse_action
from integral.action import Action


class ParserTest(unittest.TestCase):
    def testParseTerm(self):
        test_data = [
            "x", "1", "11/10", "-1", "-11/10",
            "x + y", "x - y", "-x", "x * y", "x / y", "x ^ y",
            "x + y * z", "(x + y) * z",
            "x * y + z", "x * (y + z)",
            "x * y ^ 2", "(x * y) ^ 2",
            "sin(x)", "cos(x)", "log(x)", "exp(x)",
            "D x. 3 * x",
            "INT x:[1,2]. 3 * x",
            "[3 * x]_x=1,2",
            "INT x:[0,pi / 4]. sin(x)",
            "x ^ (1/2)",
            '(-2) ^ n',
            '(-2) ^ (n + 1)',
            'a ^ b * c ^ d',
            '(-1) ^ n * factorial(n) / (m + 1) ^ (n + 1)',
            'x ^ m * log(x) ^ n',
            '(-c) ^ k / (k * a + 1) ^ (k + 1)',
            'x ^ (c * x ^ a)',
            '(c * x ^ a * log(x)) ^ k',
            '(c * x ^ a) ^ k * log(x) ^ k',
        ]

        for s in test_data:
            e = parse_expr(s)
            self.assertEqual(str(e), s)

    def testParseTerm2(self):
        test_data = [
            ("-x", -Var("x")),
            ("-2", Const(-2)),
            ("1/2", Const(Fraction(1) / 2)),
            ("-1/2", Const(Fraction(-1) / 2)),
            # ("0.5", Const(Decimal("0.5"))),
            ("pi", Fun("pi")),
            ("-x^2", Op("-", Op("^", Var("x"), Const(2)))),
            ('a ^ b * c ^ d', Op('*', Op('^', Var('a'), Var('b')), Op('^', Var('c'), Var('d')))),
            ('a * -b', Op('*', Var('a'), Op('-', Var('b')))),
            ('(-x) ^ k * log(x)', Op('*', Op('^', Op('-',Var('x')), Var('k')), Fun('log', Var('x')))),
            ("x^m * log(x) ^ n", Op('*', Op('^', Var('x'), Var('m')), Op('^', Fun('log',Var('x')), Var('n')))),
            ("(-1)^n * factorial(n) / (m+1) ^ (n+1)", Op('/', Op('*', Op('^', Const(-1), Var('n')), Fun("factorial", Var('n'))), Op('^', Op('+', Var('m'), Const(1)), Op('+', Var('n'), Const(1))))),
            ("(-c) ^ k / (k * a + 1) ^ (k + 1)", Op('/', Op('^', Op('-', Var('c')), Var('k')), Op('^', Op('+', Op('*', Var('k'), Var('a')), Const(1)), Op('+', Var('k'), Const(1))))),
        ]

        for s, e, in test_data:
            self.assertEqual(parse_expr(s), e)

    def testParseAction(self):
        test_data = [
            "calculate INT x:[0,1]. (3 * x + 1) ^ (-2)",
            "substitute u for 3 * x + 1",
            "substitute sin(t) for x",
            "apply integral identity",
            "simplify"
        ]

        for s in test_data:
            action = parse_action(s)
            self.assertIsInstance(action, Action)
            self.assertEqual(str(action), s)


if __name__ == "__main__":
    unittest.main()
