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
    def testSubstitutionIndefinite(self):
        ctx = context.Context()

        t = parse_expr("INT x. sin(5 * x)")
        rule = rules.Substitution("u", parse_expr("5 * x"))
        self.assertEqual(rule.eval(t, ctx), parse_expr("INT u. sin(u) / 5"))

    def testSubstitutionDefinite(self):
        ctx = context.Context()

        t = parse_expr("INT x:[0, pi/5]. sin(5 * x)")
        rule = rules.Substitution("u", parse_expr("5 * x"))
        self.assertEqual(rule.eval(t, ctx), parse_expr("INT u:[0, pi]. sin(u) / 5"))

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

    def testSubstitutionCondCheck3(self):
        file = compstate.CompFile("base", "standard")
        ctx = file.ctx

        t = parse_expr("INT x. (x^2 - 1)^(3/2) / x")
        ctx.add_condition(parse_expr("x > 1"))

        rule = rules.Substitution("u", parse_expr("acos(1/x)"))
        t2 = rule.eval(t, ctx)
        self.assertEqual(t2, parse_expr("INT u. cos(u) * sin(u) / cos(u) ^ 2 * (1 / cos(u) ^ 2 - 1) ^ (3/2)"))
        ctx2 = rule.update_context(ctx)

        self.assertTrue(condprover.check_condition(parse_expr("sin(u) >= 0"), ctx2))
        self.assertTrue(condprover.check_condition(parse_expr("cos(u) >= 0"), ctx2))

    def testSubstitutionCondCheck3a(self):
        file = compstate.CompFile("base", "standard")
        ctx = file.ctx

        t = parse_expr("INT x. (x^2 - 1)^(3/2) / x")
        ctx.add_condition(parse_expr("x > 1"))

        rule = rules.Substitution("u", parse_expr("asec(x)"))
        t2 = rule.eval(t, ctx)
        self.assertEqual(t2, parse_expr("INT u. sec(u) * tan(u) * (sec(u) ^ 2 - 1) ^ (3/2) / sec(u)"))
        ctx2 = rule.update_context(ctx)

        self.assertTrue(condprover.check_condition(parse_expr("tan(u) >= 0"), ctx2))
        self.assertTrue(condprover.check_condition(parse_expr("sec(u) > 0"), ctx2))

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

    def testSubstitutionCondCheck4(self):
        file = compstate.CompFile("base", "simple_integral_01")
        ctx = file.ctx

        t = parse_expr("INT x. 1 / sqrt(x ^ 2 + 1)")
        ctx.add_condition(parse_expr("x > 1"))

        rule = rules.Substitution("u", parse_expr("atan(x)"))
        ctx2 = rule.update_context(ctx)

        self.assertTrue(condprover.check_condition(parse_expr("tan(u) > 0"), ctx2))
        self.assertTrue(condprover.check_condition(parse_expr("u > pi / 4"), ctx2))
        self.assertTrue(condprover.check_condition(parse_expr("u < pi / 2"), ctx2))
        self.assertTrue(condprover.check_condition(parse_expr("sec(u) < sqrt(2)"), ctx2))
        self.assertTrue(condprover.check_condition(parse_expr("sec(u) > 1"), ctx2))
    def testSubstitutionCondCheck5(self):
        file = compstate.CompFile("base", "simple_integral_01")
        ctx = file.ctx

        t = parse_expr("INT x. sqrt(x ^ 2 + 1)")

        rule = rules.Substitution("u", "atan(x)")
        ctx2 = rule.update_context(ctx)

        self.assertTrue(condprover.check_condition(parse_expr("u > -pi/2"), ctx2))
        self.assertTrue(condprover.check_condition(parse_expr("u < pi/2"), ctx2))
        self.assertTrue(condprover.check_condition(parse_expr("cos(u) > 0"), ctx2))
        self.assertTrue(condprover.check_condition(parse_expr("sec(u) > 0"), ctx2))

    def testRewriting1(self):
        file = compstate.CompFile("base", "test_rewriting1")
        ctx = file.ctx
        ctx.load_book("interesting")

        # "conds": [ "m + 1 > 0" ],
        # "expr": "(INT x:[0,1].x^m * log(x)^n) = (-1)^n * factorial(n) / (m+1)^(n+1)",
        # "path": "chapter4_practice02",
        # "type": "problem"

        e = parse_expr("SUM(k, 0, oo, c ^ k / factorial(k) * (INT x:[0,1]. x ^ (a * k) * log(x) ^ k))")
        ctx.add_condition("a > 0")
        ctx.add_condition("c != 0")

        e = rules.IntegralIdentity().eval(e, ctx)
        self.assertEqual(e, parser.parse_expr("SUM(k, 0, oo, c ^ k / factorial(k) * ((-1) ^ k * factorial(k) * (a * k + 1) ^ (-k - 1)))"))


    def testIntegralIdentity(self):
        ctx = context.Context()
        ctx.load_book("base")

        e = parser.parse_expr("INT x. exp(3*x)")
        e = rules.IntegralIdentity().eval(e, ctx)
        self.assertEqual(e, parser.parse_expr("exp(3 * x) / 3 + SKOLEM_CONST(C)"))

        e = parser.parse_expr("INT x. exp(3*x+2)")
        e = rules.IntegralIdentity().eval(e, ctx)
        self.assertEqual(e, parser.parse_expr("exp(3 * x + 2) / 3 + SKOLEM_CONST(C)"))

        e = parser.parse_expr("INT x. sin(3*x+2)")
        e = rules.IntegralIdentity().eval(e, ctx)
        self.assertEqual(e, parser.parse_expr("-cos(3 * x + 2) / 3 + SKOLEM_CONST(C)"))

        e = parser.parse_expr("INT x. cos(3*x+2)")
        e = rules.IntegralIdentity().eval(e, ctx)
        self.assertEqual(e, parser.parse_expr("sin(3 * x + 2) / 3 + SKOLEM_CONST(C)"))

        e = parser.parse_expr("INT x. cos(-x)")
        e = rules.IntegralIdentity().eval(e, ctx)
        self.assertEqual(e, parser.parse_expr("-sin(-x) + SKOLEM_CONST(C)"))

        e = parser.parse_expr("INT x. sin(-x)")
        e = rules.IntegralIdentity().eval(e, ctx)
        self.assertEqual(e, parser.parse_expr("cos(-x) + SKOLEM_CONST(C)"))

        e = parser.parse_expr("INT x. 1 / (2 * x - 3)")
        e = rules.IntegralIdentity().eval(e, ctx)
        self.assertEqual(e, parser.parse_expr("log(abs(2*x-3))/2 + SKOLEM_CONST(C)"))

if __name__ == "__main__":
    unittest.main()
