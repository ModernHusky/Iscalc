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

    def testSimplifyOp(self):
        t1 = parser.parse_expr("3/COMPLEX(1,1)")
        ctx1 = context.Context()
        simp_t1 = poly.simplify_op(t1, ctx1)
        self.assertEqual(simp_t1, parser.parse_expr("COMPLEX(3/2,-3/2)"))

        # 无差异但是执行失败
        # t2 = parser.parse_expr("COMPLEX(1,1)/3")
        # ctx2 = context.Context()
        # simp_t2 = poly.simplify_op(t2, ctx2)
        # self.assertEqual(simp_t2, parser.parse_expr("COMPLEX(1/3,1/3)"))

        t3 = parser.parse_expr("COMPLEX(1,1)*3")
        ctx3 = context.Context()
        simp_t3 = poly.simplify_op(t3, ctx3)
        self.assertEqual(simp_t3, parser.parse_expr("COMPLEX(3,3)"))

        t4 = parser.parse_expr("3*COMPLEX(1,1)")
        ctx4 = context.Context()
        simp_t4 = poly.simplify_op(t4, ctx4)
        self.assertEqual(simp_t4, parser.parse_expr("COMPLEX(3,3)"))

        t5 = parser.parse_expr("COMPLEX(1,1)^3")
        ctx5 = context.Context()
        simp_t5 = poly.simplify_op(t5, ctx5)
        self.assertEqual(simp_t5, parser.parse_expr("COMPLEX(-2,2)"))


if __name__ == "__main__":
    unittest.main()
