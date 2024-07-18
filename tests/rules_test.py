"""Unit test for rules."""

import unittest

from integral import parser
from integral import context
from integral import rules
from integral.rules import RuleException


class RulesTest(unittest.TestCase):
    def testSubstitutionInverse(self):
        # Basic correct case
        ctx = context.Context()
        t = parser.parse_expr("INT x:[0,1]. x ^ 2")
        rule = rules.SubstitutionInverse("u", "x", parser.parse_expr("u + 1"))
        res = "INT u:[-1,0]. (u + 1) ^ 2 * 1"
        self.assertEqual(rule.eval(t, ctx), parser.parse_expr(res))

    def testSubstitutionInverseWrong(self):
        # Substitution variable already used, case 1
        ctx = context.Context()
        t = parser.parse_expr("INT x:[0,1]. x + y")
        rule = rules.SubstitutionInverse("x", "x", parser.parse_expr("x + 1"))
        self.assertRaises(RuleException, rule.eval, t, ctx)

    def testSubstitutionInverseWrong2(self):
        # Substitution variable already used, case 2
        ctx = context.Context()
        t = parser.parse_expr("INT x:[0,1]. x + y")
        rule = rules.SubstitutionInverse("y", "x", parser.parse_expr("y + 1"))
        self.assertRaises(RuleException, rule.eval, t, ctx)

    def testSubstitutionInverseIndef(self):
        # Substitution for indefinite integrals
        ctx = context.Context()
        t = parser.parse_expr("INT x. 1 / (1 + sqrt(x))")
        rule = rules.SubstitutionInverse("u", "x", parser.parse_expr("u ^ 2"))
        res = "INT u. 1 / (1 + sqrt(u ^ 2)) * (2 * u)"
        self.assertEqual(rule.eval(t, ctx), parser.parse_expr(res))


if __name__ == "__main__":
    unittest.main()
