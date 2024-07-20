"""Unit test for normalization module."""

import unittest

from integral import norm
from integral.context import Context
from integral.parser import parse_expr


class NormTest(unittest.TestCase):
    def testQuotientTrig(self):
        ctx = Context()

        test_data = [
            ("tan(x)", "sin(x) / cos(x)"),
            ("cot(x)", "cos(x) / sin(x)"),
            ("sec(x)", "1 / cos(x)"),
            ("csc(x)", "1 / sin(x)"),
            ("tan(x) * sec(x)", "sin(x) / cos(x) ^ 2")
        ]

        for t1, t2 in test_data:
            t1 = parse_expr(t1)
            t2 = parse_expr(t2)
            self.assertEqual(norm.normalize_quotient(t1, ctx).to_expr(), t2)


if __name__ == "__main__":
    unittest.main()
