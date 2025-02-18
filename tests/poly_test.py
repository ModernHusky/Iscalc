"""Unit test for polynomial module."""

import unittest

from integral import parser
from integral import context
from integral import poly
from integral.conditions import Conditions

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

    def testComplexAbs(self):
        # Test abs of complex numbers
        ctx = context.Context()
        conds = Conditions([])
        
        # Test |i| = 1
        t1 = parser.parse_expr("abs(i)")
        self.assertEqual(poly.simplify_complex(t1, conds), parser.parse_expr("1"))
        
        # Test |3 + 4i| = 5
        t2 = parser.parse_expr("abs(3 + 4 * i)")
        self.assertEqual(poly.simplify_complex(t2, conds), parser.parse_expr("5"))
        
        # Test |2i| = 2
        t3 = parser.parse_expr("abs(2 * i)")
        self.assertEqual(poly.simplify_complex(t3, conds), parser.parse_expr("2"))

    def testComplexMultiplication(self):
        # Test complex number multiplication
        ctx = context.Context()
        conds = Conditions([])
        
        # Test (1 + i)(1 - i) = 2
        t1 = parser.parse_expr("(1 + i) * (1 - i)")
        self.assertEqual(poly.simplify_complex(t1, conds), parser.parse_expr("2"))
        
        # Test i * i = -1
        t2 = parser.parse_expr("i * i")
        self.assertEqual(poly.simplify_complex(t2, conds), parser.parse_expr("-1"))
        
        # Test (2 + 3i)(4 + 5i) = (8 - 15) + (12 + 20)i = -7 + 32i
        t3 = parser.parse_expr("(2 + 3 * i) * (4 + 5 * i)")
        self.assertEqual(poly.simplify_complex(t3, conds), parser.parse_expr("-7 + 32 * i"))
        
        # Test distributive property
        t4 = parser.parse_expr("2 * (3 + 4 * i)")
        self.assertEqual(poly.simplify_complex(t4, conds), parser.parse_expr("6 + 8 * i"))

if __name__ == "__main__":
    unittest.main()
