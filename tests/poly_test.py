"""Unit test for polynomial module."""

import unittest

from integral import parser
from integral import context
from integral import poly

class PolyTest(unittest.TestCase):
    def testNormalizeAlpha(self):
        t = parser.parse_expr("(INT x:[1,2]. x ^ 2) - (INT y:[1,2]. y ^ 2)")
        ctx = context.Context()
        p = poly.to_poly(t, ctx)
        self.assertEqual(p, poly.Polynomial(tuple()))  # equal to zero

    def testSimplifyLog(self):
        t = parser.parse_expr("log(10)")
        ctx = context.Context()
        simp_t = poly.simplify_log(t, ctx)
        self.assertEqual(simp_t, parser.parse_expr("log(2) + log(5)"))

    def testSimplifyAbs(self):
        t = parser.parse_expr("pi / abs(4 * cos(a))")
        ctx = context.Context()
        ctx.load_book("base")
        ctx.add_condition("cos(a) < 0")
        ctx.add_condition("cos(a) != 0")
        simp_t = poly.normalize(t, ctx)
        self.assertEqual(simp_t, parser.parse_expr("-(pi / (4 * cos(a)))"))


if __name__ == "__main__":
    unittest.main()
