"""Unit test for rules."""

import unittest

from integral import parser
from integral.parser import parse_expr
from integral import context
from integral import rules
from integral import condprover
from integral import compstate
from integral.rules import RuleException


class RulesTest(unittest.TestCase):
    def testSubstitutionCondCheck(self):
        # After substituting u for sqrt(5 + sqrt(x)), should be able to derive
        # u ^ 2 - 5 >= 0, so simplify abs(u ^ 2 - 5) to u ^ 2 - 5.
        file = compstate.CompFile("base", "standard")
        ctx = file.ctx

        t = parse_expr("INT x. sqrt(5 + sqrt(x))")
        ctx.add_condition(parse_expr("x > 0"))

        rule = rules.Substitution("u", parse_expr("sqrt(5 + sqrt(x))"))
        t2 = rule.eval(t, ctx)
        self.assertEqual(t2, parse_expr("INT u. 4 * u * (u ^ 2 - 5) * sqrt(abs(u ^ 2 - 5) + 5)"))
        ctx2 = rule.update_context(ctx)

        t3 = rules.Simplify().eval(t2, ctx2)
        self.assertEqual(t3, parse_expr("4 * (INT u. u ^ 2 * (u ^ 2 - 5))"))

    def testSubstitutionCondCheck2(self):
        # After substituting u for (1 + sqrt(x - 3)) ^ (1/3), should be able to
        # derive u ^ 3 - 1 >= 0, so simplify abs(u ^ 3 - 1) to u ^ 3 - 1
        file = compstate.CompFile("base", "standard")
        ctx = file.ctx

        t = parse_expr("INT x. (1 + sqrt(x - 3)) ^ (1/3)")
        ctx.add_condition(parse_expr("x > 3"))

        rule = rules.Substitution("u", parse_expr("(1 + sqrt(x - 3)) ^ (1/3)"))
        t2 = rule.eval(t, ctx)
        self.assertEqual(t2, parse_expr("INT u. 6 * u ^ 2 * (u ^ 3 - 1) * (abs(u ^ 3 - 1) + 1) ^ (1/3)"))
        ctx2 = rule.update_context(ctx)

        t3 = rules.Simplify().eval(t2, ctx2)
        self.assertEqual(t3, parse_expr("6 * (INT u. u ^ 3 * (u ^ 3 - 1))"))

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
