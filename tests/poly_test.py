"""Unit test for polynomial module."""

import unittest

from sympy.matrices.expressions.kronecker import rules

from integral import parser
from integral import context
from integral import poly
from integral.conditions import Conditions

import os
os.chdir('E:\\=graduatelife======\\learn-git\\iscalc')

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

    def testComplexMultiplication(self):
        # Test complex number multiplication
        ctx = context.Context()

        # Test i * i = -1
        t1 = parser.parse_expr("i * i")
        self.assertEqual(poly.normalize(t1, ctx), parser.parse_expr("-1"))

if __name__ == "__main__":
    unittest.main()
