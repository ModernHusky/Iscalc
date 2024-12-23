"""Overall test for integrals."""

import unittest
import json
import sys,os
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(BASE_DIR)
from integral import expr
from integral import compstate
from integral import rules
from integral import parser


class IntegralTest(unittest.TestCase):
    def checkAndOutput(self, file: compstate.CompFile, omit_finish: bool = False):
        # Test parsing of json file
        json_file = file.export()
        for i, item in enumerate(json_file['content']):
            aa, bb = compstate.parse_item(file.content[i].parent, item), file.content[i]
            a = aa.export()
            b = bb.export()
            # if a != b:
            #     if isinstance(aa, compstate.Goal) and isinstance(bb, compstate.Goal):
            #         aa.is_finished()
            #     with open('examples/a.json', 'w', encoding='utf-8') as f:
            #         json.dump(a, f, indent=4, ensure_ascii=False, sort_keys=True)
            #     with open('examples/b.json', 'w', encoding='utf-8') as f:
            #         json.dump(b, f, indent=4, ensure_ascii=False, sort_keys=True)
            self.assertEqual(a, b)

        # Output to file
        with open('examples/' + file.name + '.json', 'w', encoding='utf-8') as f:
            json.dump(file.export(), f, indent=4, ensure_ascii=False, sort_keys=True)

        # Test goals are finished
        if not omit_finish:
            for content in file.content:
                self.assertTrue(content.is_finished())

    def testStandard(self):
        file = compstate.CompFile("base", "standard")

        goal1 = file.add_goal("(INT x. 1 / (x + a)) = log(abs(x + a)) + SKOLEM_CONST(C)", conds=["x + a != 0"])
        proof = goal1.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Substitution("u", parser.parse_expr("x + a")))
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.ReplaceSubstitution())
        assert goal1.is_finished()

        goal2 = file.add_goal("(INT x. exp(a * x)) = exp(a * x) / a + SKOLEM_CONST(C)", conds=["a != 0"])
        proof = goal2.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Substitution("u", parser.parse_expr("a * x")))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.ReplaceSubstitution())
        assert goal2.is_finished()

        goal3 = file.add_goal("(INT x. sin(a * x)) = - (cos(a * x) / a) + SKOLEM_CONST(C)", conds=["a != 0"])
        proof = goal3.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Substitution("u", parser.parse_expr("a * x")))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.ReplaceSubstitution())
        calc.perform_rule(rules.Simplify())
        assert goal3.is_finished()

        goal4 = file.add_goal("(INT x. cos(a * x)) = sin(a * x) / a + SKOLEM_CONST(C)", conds=["a != 0"])
        proof = goal4.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Substitution("u", parser.parse_expr("a * x")))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.ReplaceSubstitution())
        assert goal4.is_finished()

        goal5 = file.add_goal("(INT x. 1 / (a ^ 2 + x ^ 2)) = (1/a) * atan(x / a) + SKOLEM_CONST(C)", conds=["a != 0"])
        proof = goal5.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Substitution("u", parser.parse_expr("x / a")))
        calc.perform_rule(rules.Equation("a ^ 2 * u ^ 2 + a ^ 2", "a ^ 2 * (u ^ 2 + 1)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.ReplaceSubstitution())
        calc.perform_rule(rules.Simplify())
        assert goal5.is_finished()

        goal6 = file.add_goal(
            "(INT x. x ^ k * log(x)) = x ^ (k + 1) * log(x) / (k + 1) - x ^ (k + 1) / (k + 1) ^ 2 + SKOLEM_CONST(C)", conds=["x > 0", "k != -1"])
        proof = goal6.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.IntegrationByParts(
            u=parser.parse_expr("log(x)"),
            v=parser.parse_expr("x ^ (k+1) / (k+1)")))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        assert goal6.is_finished()

        goal7 = file.add_goal("(INT x:[0,1]. x ^ m * log(x) ^ n) = (-1)^n * factorial(n) / (m+1) ^ (n+1)", conds=["m >= 0", "n >= 0", "isInt(n)"])
        proof = goal7.proof_by_induction("n")
        proof_base = proof.base_case.proof_by_calculation()
        proof_induct = proof.induct_case.proof_by_calculation()

        # Base case
        calc = proof_base.lhs_calc
        # calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        calc = proof_base.rhs_calc
        calc.perform_rule(rules.Simplify())
        calc = proof_induct.rhs_calc
        calc.perform_rule(rules.Simplify())
        calc = proof_induct.lhs_calc
        # calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IntegrationByParts(
            u=parser.parse_expr("log(x) ^ (n + 1)"),
            v=parser.parse_expr("x ^ (m+1) / (m+1)")))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.OnSubterm(rules.ApplyInductHyp()))
        calc.perform_rule(rules.Simplify())
        s = parser.parse_expr("(-1) ^ (n + 1) * (m + 1) ^ (-n - 2) * ((n + 1) * factorial(n))")
        calc.perform_rule(rules.Equation(None, s))
        s1 = parser.parse_expr("(n + 1) * factorial(n)")
        s2 = parser.parse_expr("factorial(n + 1)")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        calc.perform_rule(rules.Simplify())
        assert goal7.is_finished()

        goal8 = file.add_goal("(INT x:[0,oo]. exp(-(x * y)) * sin(a * x)) = a / (a ^ 2 + y ^ 2)", conds=["y > 0"])
        proof = goal8.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(
            rules.IntegrationByParts(parser.parse_expr("exp(-(x * y))"), parser.parse_expr("-cos(a * x) / a")))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(
            rules.IntegrationByParts(parser.parse_expr("exp(-(x * y))"), parser.parse_expr("sin(a * x) / a")))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IntegrateByEquation(parser.parse_expr("INT x:[0,oo]. exp(-(x * y)) * sin(a * x)")))
        calc.perform_rule(rules.Equation(None, "a / (a ^ 2 + y ^ 2)"))
        assert goal8.is_finished()

        goal9 = file.add_goal("(INT x. a ^ x) = a ^ x / log(a) + SKOLEM_CONST(C)",
                              conds = ["a>0", "a!=1"])
        proof = goal9.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Equation("a ^ x", "exp(log(a) * x)"))
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        assert goal9.is_finished()

        goal10 = file.add_goal("(INT x. cos(x) ^ 2) = 1/2 * (sin(2 * x) / 2 + x) + SKOLEM_CONST(C)")
        proof = goal10.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.ApplyIdentity("cos(x)^2", "(1 + cos(2*x)) / 2"))
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        assert goal10.is_finished()
        # self.checkAndOutput(file)

    def testTongji(self):
        file = compstate.CompFile("tongji", "tongji7")

        calc = file.add_calculation("INT x:[2,3]. 2 * x + x ^ 2")
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "34/3")

        calc = file.add_calculation("INT x:[0,1]. (3 * x + 1) ^ (-2)")
        calc.perform_rule(rules.Substitution("u", parser.parse_expr("3 * x + 1")))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "1/4")

        calc = file.add_calculation("INT x:[0,1]. exp(6 * x)")
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "exp(6) / 6 - 1/6")

        calc = file.add_calculation("INT x:[-1,2]. x * exp(x)")
        calc.perform_rule(rules.IntegrationByParts(parser.parse_expr("x"), parser.parse_expr("exp(x)")))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "2 * exp(-1) + exp(2)")

        calc = file.add_calculation("INT x:[0,pi/4]. sin(x)")
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "-(sqrt(2) / 2) + 1")

        calc = file.add_calculation("INT x:[0,1]. 3*x^2 - x + 1")
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "3/2")

        calc = file.add_calculation("INT x:[1,2]. x^2 + 1/x^4")
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "21/8")

        calc = file.add_calculation("INT x:[pi/3, pi]. sin(2*x + pi/3)")
        calc.perform_rule(rules.Substitution("u", parser.parse_expr("2*x + pi/3")))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "-3/4")

        calc = file.add_calculation("INT x:[4, 9]. x ^ (1 / 3) * (x ^ (1 / 2) + 1)")
        calc.perform_rule(rules.ExpandPolynomial())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "-(81 * 2 ^ (2/3) / 11) + 945 * 3 ^ (2/3) / 44")

        calc = file.add_calculation("INT x:[-1, 0]. (3 * x ^ 4 + 3 * x ^ 2 + 1) / (x ^ 2 + 1)")
        calc.perform_rule(rules.PartialFractionDecomposition())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "pi / 4 + 1")

        calc = file.add_calculation("INT x:[4, exp(1) + 3]. (x ^ 3 - 12 * x ^ 2 - 42) / (x - 3)")
        calc.perform_rule(rules.PartialFractionDecomposition())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Substitution("u", parser.parse_expr("x - 3")))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.ExpandPolynomial())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "-(3 * exp(2) / 2) + exp(3) / 3 - 45 * exp(1) - 461/6")

        calc = file.add_calculation("INT x:[0, pi / 2]. sin(x) * cos(x) ^ 3")
        calc.perform_rule(rules.Substitution("u", parser.parse_expr("cos(x)")))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "1/4")

        calc = file.add_calculation("INT x:[0, pi]. (1 - sin(x)^3)")
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("sin(x)^3", "sin(x) * sin(x)^2"))
        calc.perform_rule(rules.ApplyIdentity("sin(x)^2", "1 - cos(x)^2"))
        calc.perform_rule(rules.Substitution("u", parser.parse_expr("cos(x)")))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "pi - 4/3")

        calc = file.add_calculation("INT x:[pi/6, pi/2]. cos(x) ^ 2")
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "-(sqrt(3) / 8) + pi / 6")

        calc = file.add_calculation("INT x:[0, 1]. (1 - x^2) ^ (1/2)")
        calc.perform_rule(rules.SubstitutionInverse("u", "x", parser.parse_expr("sin(u)")))
        calc.perform_rule(rules.ApplyIdentity("1-sin(u)^2", "cos(u)^2"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "pi / 4")

        calc = file.add_calculation("INT x:[0, sqrt(2)]. sqrt(2 - x^2)")
        calc.perform_rule(rules.SubstitutionInverse("u", "x", parser.parse_expr("sqrt(2) * sin(u)")))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.ApplyIdentity("sin(u)^2", "1 - cos(u)^2"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "pi / 2")

        calc = file.add_calculation("INT y:[-sqrt(2), sqrt(2)]. sqrt(8 - 2*y^2)")
        calc.perform_rule(rules.SubstitutionInverse("u", "y", parser.parse_expr("2 * sin(u)")))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.ApplyIdentity("sin(u)^2", "1 - cos(u)^2"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.ExpandPolynomial())
        self.assertEqual(str(calc.last_expr), "sqrt(2) * pi + 2 * sqrt(2)")

        calc = file.add_calculation("INT x:[1/sqrt(2), 1]. sqrt(1 - x^2) / x ^ 2")
        calc.perform_rule(rules.SubstitutionInverse("u", "x", parser.parse_expr("sin(u)")))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.ApplyIdentity("sin(u)^2", "1 - cos(u)^2"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.ApplyIdentity("cos(u)^2", "1 - sin(u)^2"))
        calc.perform_rule(rules.ExpandPolynomial())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.ApplyIdentity("1 / sin(u) ^ 2", "csc(u)^2"))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "-(pi / 4) + 1")

        calc = file.add_calculation("INT x:[-1, 1]. x / sqrt(5 - 4 * x)")
        calc.perform_rule(rules.Substitution("u", parser.parse_expr("5 - 4 * x")))
        calc.perform_rule(rules.ExpandPolynomial())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "1/6")

        calc = file.add_calculation("INT x:[1,4]. 1 / (1 + sqrt(x))")
        calc.perform_rule(rules.Substitution("u", parser.parse_expr("sqrt(x)")))
        calc.perform_rule(rules.Substitution("v", parser.parse_expr("u + 1")))
        calc.perform_rule(rules.ExpandPolynomial())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "2 * log(2) - 2 * log(3) + 2")

        calc = file.add_calculation("INT x:[3/4, 1]. 1 / (sqrt(1-x) - 1)")
        calc.perform_rule(rules.Substitution("u", parser.parse_expr("sqrt(1 - x)")))
        calc.perform_rule(rules.Substitution("v", parser.parse_expr("u - 1")))
        calc.perform_rule(rules.ExpandPolynomial())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "-(2 * log(2)) + 1")

        calc = file.add_calculation("INT t:[0, 1]. t * exp(-t ^ 2 / 2)")
        calc.perform_rule(rules.Substitution("u", parser.parse_expr("t ^ 2 / 2")))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "-exp(-1/2) + 1")

        calc = file.add_calculation("INT x:[1, exp(2)]. 1 / (x*sqrt(1+log(x)))")
        calc.perform_rule(rules.Substitution("u", parser.parse_expr("1 + log(x)")))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "2 * sqrt(3) - 2")

        calc = file.add_calculation("INT x:[-2, 0]. (x + 2)/(x^2 + 2*x + 2)")
        calc.perform_rule(rules.Equation("x^2 + 2*x + 2", "(x + 1) ^ 2 + 1"))
        calc.perform_rule(rules.Substitution("u", parser.parse_expr("x + 1")))
        calc.perform_rule(rules.ExpandPolynomial())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.SplitRegion(parser.parse_expr("0")))
        calc.perform_rule(rules.Substitution("v", parser.parse_expr("u ^ 2 + 1")))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Substitution("v", parser.parse_expr("u ^ 2 + 1")))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "pi / 2")

        calc = file.add_calculation("INT x:[-pi/2,pi/2]. cos(x)^4")
        calc.perform_rule(rules.Equation("cos(x)^4", "(cos(x)^2)^2"))
        calc.perform_rule(rules.ApplyIdentity("cos(x)^2", "(1 + cos(2*x)) / 2"))
        calc.perform_rule(rules.ExpandPolynomial())
        calc.perform_rule(rules.Substitution("u", parser.parse_expr("2 * x")))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "3 * pi / 8")

        # Need to check substitution is monotonic
        calc = file.add_calculation("INT x:[-pi/2, pi/2]. sqrt(cos(x) - cos(x)^3)")
        calc.perform_rule(rules.Equation("cos(x) - cos(x)^3", "cos(x) * (1 - cos(x)^2)"))
        calc.perform_rule(rules.ApplyIdentity("1 - cos(x)^2", "sin(x) ^ 2"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.SplitRegion(parser.parse_expr("0")))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Substitution("u", parser.parse_expr("cos(x)")))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Substitution("u", parser.parse_expr("cos(x)")))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "4/3")

        calc = file.add_calculation("INT x:[0, pi]. sqrt(1 + cos(2*x))")
        calc.perform_rule(rules.ApplyIdentity("cos(2*x)", "2 * cos(x)^2 - 1"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.SplitRegion(parser.parse_expr("pi / 2")))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "2 * sqrt(2)")

        calc = file.add_calculation("INT x:[0, 1].x * exp(-x)")
        calc.perform_rule(rules.IntegrationByParts(parser.parse_expr("x"), parser.parse_expr("-exp(-x)")))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "-(2 * exp(-1)) + 1")

        calc = file.add_calculation("INT x:[1, exp(1)]. x * log(x)")
        calc.perform_rule(rules.IntegrationByParts(parser.parse_expr("log(x)/2"), parser.parse_expr("x^2")))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "exp(2) / 4 + 1/4")

        calc = file.add_calculation("INT x:[pi/4, pi/3]. x / sin(x)^2")
        calc.perform_rule(rules.ApplyIdentity("sin(x)^2", "csc(x)^(-2)"))
        calc.perform_rule(rules.IntegrationByParts(parser.parse_expr("x"), parser.parse_expr("-cot(x)")))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.ApplyIdentity("cot(x)", "cos(x) / sin(x)"))
        calc.perform_rule(rules.Substitution("u", parser.parse_expr("sin(x)")))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "-(sqrt(3) * pi / 9) - log(2) / 2 + log(3) / 2 + pi / 4")

        calc = file.add_calculation("INT x:[1, 4]. log(x) / sqrt(x)")
        calc.perform_rule(rules.IntegrationByParts(parser.parse_expr("2*log(x)"), parser.parse_expr("sqrt(x)")))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "8 * log(2) - 4")

        calc = file.add_calculation("INT x:[0, 1]. x * atan(x)")
        calc.perform_rule(rules.IntegrationByParts(parser.parse_expr("atan(x)/2"), parser.parse_expr("x^2")))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("x^2 / (2 * x^2 + 2)", "(1 - 1 / (x^2 + 1)) / 2"))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "pi / 4 - 1/2")

        calc = file.add_calculation("INT x:[0, pi/2]. exp(2*x)*cos(x)")
        calc.perform_rule(rules.IntegrationByParts(parser.parse_expr("exp(2*x)"), parser.parse_expr("sin(x)")))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IntegrationByParts(parser.parse_expr("exp(2*x)"), parser.parse_expr("-cos(x)")))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IntegrateByEquation(parser.parse_expr("INT x:[0, pi/2]. exp(2*x)*cos(x)")))
        self.assertEqual(str(calc.last_expr), "(exp(pi) - 2) / 5")

        calc = file.add_calculation("INT x:[0,pi]. (x * sin(x))^2")
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.ApplyIdentity("sin(x)^2", "(1 - cos(2*x)) / 2"))
        calc.perform_rule(rules.ExpandPolynomial())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IntegrationByParts(parser.parse_expr("x^2/2"), parser.parse_expr("sin(2*x)")))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IntegrationByParts(parser.parse_expr("x/2"), parser.parse_expr("-cos(2*x)")))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "pi ^ 3 / 6 - pi / 4")

        calc = file.add_calculation("INT x:[1, exp(1)]. sin(log(x))")
        calc.perform_rule(rules.Substitution("u", parser.parse_expr("log(x)")))
        calc.perform_rule(rules.IntegrationByParts(parser.parse_expr("-exp(u)"), parser.parse_expr("cos(u)")))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IntegrationByParts(parser.parse_expr("exp(u)"), parser.parse_expr("sin(u)")))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IntegrateByEquation(parser.parse_expr("INT u:[0,1]. exp(u) * sin(u)")))
        self.assertEqual(str(calc.last_expr), "(-(cos(1) * exp(1)) + exp(1) * sin(1) + 1) / 2")

        calc = file.add_calculation("INT x:[1/exp(1), exp(1)]. abs(log(x))")
        calc.perform_rule(rules.SplitRegion(parser.parse_expr("1")))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IntegrationByParts(parser.parse_expr("log(x)"), parser.parse_expr("x")))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IntegrationByParts(parser.parse_expr("log(x)"), parser.parse_expr("x")))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "-(2 * exp(-1)) + 2")

        # self.checkAndOutput(file)

    def testMIT2019(self):
        file = compstate.CompFile("MIT", "MIT2019")

        calc = file.add_calculation("INT x:[0,pi/100]. (sin(20*x)+sin(19*x)) / (cos(20*x)+cos(19*x))")
        calc.perform_rule(rules.ApplyIdentity("sin(20*x)+sin(19*x)", "2*cos(1/2*x)*sin(39/2*x)"))
        calc.perform_rule(rules.ApplyIdentity("cos(20*x)+cos(19*x)", "2*cos(1/2*x)*cos(39/2*x)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Substitution("u", "cos(39/2*x)"))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "-(2/39 * log(cos(39 * pi / 200)))")

        # self.checkAndOutput(file)

    def testLHopital(self):
        file = compstate.CompFile("UCDavis", "LHopital")

        calc = file.add_calculation("LIM {x -> 1}. (x^2 - 1) / (x^2 + 3*x - 4)")
        calc.perform_rule(rules.LHopital())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "2/5")

        calc = file.add_calculation("LIM {x -> 4}. (x - 4) / (sqrt(x) - 2)")
        calc.perform_rule(rules.LHopital())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "4")

        calc = file.add_calculation("LIM {x -> 0}. sin(x) / x")
        calc.perform_rule(rules.LHopital())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "1")

        calc = file.add_calculation("LIM {x -> 0}. (3^x - 2^x) / (x^2 - x)")
        calc.perform_rule(rules.LHopital())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "log(2) - log(3)")

        calc = file.add_calculation("LIM {x -> 3}. (1/x - 1/3) / (x^2 - 9)")
        calc.perform_rule(rules.LHopital())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "-1/54")

        calc = file.add_calculation("LIM {x -> 0}. x * tan(x) / sin(3*x)")
        calc.perform_rule(rules.LHopital())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "0")

        # self.checkAndOutput(file)

    def testExponential(self):
        file = compstate.CompFile("UCDavis", 'Exponential')

        calc = file.add_calculation("INT x. 5*exp(x)")
        calc.perform_rule(rules.IndefiniteIntegralIdentity())

        calc = file.add_calculation("INT x. 2 - 3*exp(x)")
        calc.perform_rule(rules.IndefiniteIntegralIdentity())

        calc = file.add_calculation("INT x:[0,log(2)/7]. 14*exp(7*x)")
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        calc = file.add_calculation("INT x. 7 ^ (2*x+3)")
        calc.perform_rule(rules.Substitution("u", "2 * x"))
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.ReplaceSubstitution())
        calc.perform_rule(rules.Simplify())

        calc = file.add_calculation("INT x. exp(5*x) * (exp(2*x) / 7 + 3 / exp(3*x))")
        calc.perform_rule(rules.ExpandPolynomial())
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        calc = file.add_calculation("INT x. exp(x) * (1 + 2*exp(x)) ^ 4")
        calc.perform_rule(rules.Substitution("u", "2 * exp(x) + 1"))
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.ReplaceSubstitution())
        calc.perform_rule(rules.Simplify())

        calc = file.add_calculation("INT x. (exp(4*x) - exp(-4*x))^2")
        calc.perform_rule(rules.ExpandPolynomial())
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        calc = file.add_calculation("INT x. exp(x) * (1 - exp(x)) * (1 + exp(x)) ^ 10")
        calc.perform_rule(rules.Substitution("u", "exp(x) + 1"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.ExpandPolynomial())
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        goal01 = file.add_goal("(INT x:[0,1]. (3^x + 4^x) / 5^x) = (-2/5) / (log(3) - log(5)) - (1/5) / (2 * log(2) - log(5))")
        proof01 = goal01.proof_by_calculation()
        calc = proof01.lhs_calc
        calc.perform_rule(rules.Equation("(3^x + 4^x) / 5^x", "(3^x/5^x) + (4^x/5^x)"))
        calc.perform_rule(rules.ApplyIdentity("3^x/5^x", "(3/5)^x"))
        calc.perform_rule(rules.ApplyIdentity("4^x/5^x", "(4/5)^x"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation(None, "-2/5 / (log(3) - log(5)) - 1/5 / (2 * log(2) - log(5))"))

        goal02 = file.add_goal("(INT x. 30 * exp(-3*x) * (1 + 3 * exp(-x)) ^ 5) = "
                             "(-5/36)*(3*exp(-x) + 1)^8 + (20/63)*(3*exp(-x) + 1)^7 + (-5/27)*(3*exp(-x) + 1)^6 + SKOLEM_CONST(C)")
        proof02 = goal02.proof_by_calculation()
        calc = proof02.lhs_calc
        calc.perform_rule(rules.Substitution("u", "3 * exp(-x) + 1"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("u ^ 5 * (u - 1) ^ 2", "u ^ 7 - 2 * u ^ 6 + u ^ 5"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.ReplaceSubstitution())
        calc.perform_rule(rules.Simplify())
        calc = proof02.rhs_calc
        calc.perform_rule(rules.Simplify())

        goal03 = file.add_goal("(INT x. 8*exp(x)*(3+exp(x)) / sqrt(exp(2*x) + 6*exp(x) + 1)) = "
                             "8 * sqrt(exp(2 * x) + 6 * exp(x) + 1) + SKOLEM_CONST(C)")
        proof03 = goal03.proof_by_calculation()
        calc = proof03.lhs_calc
        calc.perform_rule(rules.Equation("8*exp(x)*(3+exp(x))", "4*(2 * exp(2 * x) + 6 * exp(x))"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Substitution("u", "exp(2*x) + 6*exp(x) + 1"))
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.ReplaceSubstitution())

        goal04 = file.add_goal("(INT x. (27*exp(9*x) + exp(12*x)) ^ (1/3)) = "
                             "(exp(3 * x) + 27) ^ (4/3) / 4 + SKOLEM_CONST(C)")
        proof04 = goal04.proof_by_calculation()
        calc = proof04.lhs_calc
        calc.perform_rule(rules.Equation("27*exp(9*x) + exp(12*x)", "exp(9*x) * (27 + exp(3*x))"))
        calc.perform_rule(rules.Equation("(exp(9*x) * (27 + exp(3*x))) ^ (1/3)",
                                         "exp(3*x) * (27 + exp(3*x)) ^ (1/3)"))
        calc.perform_rule(rules.Substitution("u", "exp(3*x) + 27"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.ReplaceSubstitution())
        goal04.print_entry()

        # self.checkAndOutput(file)

    def testTrigonometric(self):
        file = compstate.CompFile("UCDavis", 'Trigonometric')

        calc = file.add_calculation("INT x. sin(3*x)")
        calc.perform_rule(rules.IndefiniteIntegralIdentity())

        calc = file.add_calculation("INT x. tan(5*x)", conds=["x > 0", "x < pi/10"])
        calc.perform_rule(rules.ApplyIdentity("tan(5*x)", "sin(5*x) / cos(5*x)"))
        calc.perform_rule(rules.Substitution("u", "cos(5*x)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.ReplaceSubstitution())
        calc.perform_rule(rules.Simplify())

        calc = file.add_calculation("INT x. 5 * sec(4*x) * tan(4*x)", conds=["x > 0", "x < pi/8"])
        calc.perform_rule(rules.ApplyIdentity("sec(4*x)", "1 / cos(4*x)"))
        calc.perform_rule(rules.ApplyIdentity("tan(4*x)", "sin(4*x) / cos(4*x)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Substitution("u", "cos(4*x)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.ReplaceSubstitution())
        calc.perform_rule(rules.Simplify())

        calc = file.add_calculation("INT x. (sin(x) + cos(x)) ^ 2", conds=["x > 0", "x < pi/2"])
        calc.perform_rule(rules.ExpandPolynomial())
        calc.perform_rule(rules.Equation("2 * cos(x) * sin(x) + cos(x) ^ 2 + sin(x) ^ 2", "2 * cos(x) * sin(x) + (sin(x) ^ 2 + cos(x) ^ 2)"))
        calc.perform_rule(rules.ApplyIdentity("sin(x) ^ 2 + cos(x) ^ 2", "1"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Substitution("u", "sin(x)"))
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.ReplaceSubstitution())
        calc.perform_rule(rules.Simplify())

        calc = file.add_calculation("INT x. 3 * cos(5*x) ^ 2")
        calc.perform_rule(rules.Substitution("u", "5 * x"))
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.ReplaceSubstitution())
        calc.perform_rule(rules.Simplify())

        calc = file.add_calculation("INT x. (2 + tan(x)) ^ 2")
        calc.perform_rule(rules.ExpandPolynomial())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.ApplyIdentity("tan(x) ^ 2", "sec(x) ^ 2 - 1"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        calc = file.add_calculation("INT x. sin(x) ^ 3", conds=["x > 0", "x < pi/2"])
        calc.perform_rule(rules.Equation("sin(x)^3", "sin(x) * sin(x)^2"))
        calc.perform_rule(rules.ApplyIdentity("sin(x)^2", "1 - cos(x)^2"))
        calc.perform_rule(rules.ExpandPolynomial())
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.Substitution("u", "cos(x)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.ReplaceSubstitution())
        calc.perform_rule(rules.Simplify())

        calc = file.add_calculation("INT x. cos(5*x) / (3 + sin(5*x))", conds=["x > 0", "x < pi/10"])
        calc.perform_rule(rules.Substitution("u", "sin(5*x)"))
        calc.perform_rule(rules.Substitution("v", "u + 3"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.ReplaceSubstitution())
        calc.perform_rule(rules.ReplaceSubstitution())
        calc.perform_rule(rules.Simplify())

        calc = file.add_calculation("INT x. cos(x)^2 / (1 + sin(x))", conds=["x > 0", "x < pi/2"])
        calc.perform_rule(rules.ApplyIdentity("cos(x)^2", "1 - sin(x)^2"))
        calc.perform_rule(rules.Equation("1 - sin(x)^2", "(1 + sin(x)) * (1 - sin(x))"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        calc = file.add_calculation("INT x. sin(x) / (1 + sin(x))", conds=['x>-pi/2', 'x<pi/2'])
        calc.perform_rule(rules.Equation("sin(x) / (1 + sin(x))", "sin(x) * (1 - sin(x)) / (1 - sin(x)^2)"))
        calc.perform_rule(rules.ApplyIdentity("1 - sin(x)^2", "cos(x)^2"))
        calc.perform_rule(rules.Equation("sin(x) * (1 - sin(x)) / cos(x) ^ 2", "(1/cos(x)) * (sin(x)/cos(x)) - ((sin(x)/cos(x))^2)"))
        calc.perform_rule(rules.ApplyIdentity("1 / cos(x)", "sec(x)"))
        calc.perform_rule(rules.ApplyIdentity("sin(x) / cos(x)", "tan(x)"))
        calc.perform_rule(rules.ApplyIdentity("sin(x) / cos(x)", "tan(x)"))
        calc.perform_rule(rules.ApplyIdentity("tan(x)^2", "sec(x)^2 - 1"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        calc = file.add_calculation("INT x. (csc(3*x) + cot(3*x)) ^ 2")
        calc.perform_rule(rules.ExpandPolynomial())
        calc.perform_rule(rules.ApplyIdentity("cot(3*x)^2", "csc(3*x)^2 - 1"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Substitution("u", "3 * x"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.ReplaceSubstitution())
        calc.perform_rule(rules.Substitution("u", "3 * x"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.ReplaceSubstitution())
        calc.perform_rule(rules.Simplify())

        calc = file.add_calculation("INT x. (sec(x))^2 * sqrt(5 + tan(x))")
        calc.perform_rule(rules.Substitution("u", "5 + tan(x)"))
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.ReplaceSubstitution())
        calc.perform_rule(rules.Simplify())

        calc = file.add_calculation("INT x. tan(x)^5")
        calc.perform_rule(rules.Equation("tan(x)^5", "tan(x)^3 * tan(x)^2"))
        calc.perform_rule(rules.ApplyIdentity("tan(x)^2", "sec(x)^2 - 1"))
        calc.perform_rule(rules.ExpandPolynomial())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Substitution("u", "tan(x)"))
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.ReplaceSubstitution())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("tan(x)^3", "tan(x) * tan(x)^2"))
        calc.perform_rule(rules.ApplyIdentity("tan(x)^2", "sec(x)^2 - 1"))
        calc.perform_rule(rules.Equation("tan(x) * (sec(x) ^ 2 - 1)", "tan(x) * sec(x) ^ 2 - tan(x)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Substitution("u", "tan(x)"))
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.ReplaceSubstitution())
        calc.perform_rule(rules.Simplify())

        calc = file.add_calculation("INT x. (sin(x) - cos(x)) * (sin(x) + cos(x)) ^ 5")
        calc.perform_rule(rules.Equation("(sin(x) - cos(x)) * (sin(x) + cos(x)) ^ 5", "-(sin(x) + cos(x)) ^ 5 * (cos(x) - sin(x))"))
        calc.perform_rule(rules.Substitution("u", "sin(x) + cos(x)"))
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.ReplaceSubstitution())
        calc.perform_rule(rules.Simplify())

        calc = file.add_calculation("INT x. cos(x) * exp(4 + sin(x))")
        calc.perform_rule(rules.Equation("cos(x) * exp(4 + sin(x))", "exp(4 + sin(x)) * cos(x)"))
        calc.perform_rule(rules.Substitution("u", "4 + sin(x)"))
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.ReplaceSubstitution())
        calc.perform_rule(rules.Simplify())

        calc = file.add_calculation("INT x. sin(3*x) * sin(cos(3*x))")
        calc.perform_rule(rules.Equation("sin(3*x) * sin(cos(3*x))", "(-sin(cos(3*x))/3) * (-(3 * sin(3*x)))"))
        calc.perform_rule(rules.Substitution("u", "cos(3*x)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.ReplaceSubstitution())
        calc.perform_rule(rules.Simplify())

        calc = file.add_calculation("INT x. cos(x) * log(sin(x)) / sin(x)")
        calc.perform_rule(rules.Equation("cos(x) * log(sin(x)) / sin(x)", "log(sin(x)) * (cos(x) / sin(x))"))
        calc.perform_rule(rules.Substitution("u", "log(sin(x))"))
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.ReplaceSubstitution())
        calc.perform_rule(rules.Simplify())

        calc = file.add_calculation("INT x. sec(x) * tan(x) * sqrt(4 + 3 * sec(x))")
        calc.perform_rule(rules.Equation("sec(x) * tan(x) * sqrt(4 + 3 * sec(x))", "(sqrt(4 + 3 * sec(x)) / 3) * (3 * sec(x) * tan(x))"))
        calc.perform_rule(rules.Substitution("u", "4 + 3 * sec(x)"))
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.ReplaceSubstitution())
        calc.perform_rule(rules.Simplify())

        calc = file.add_calculation("INT x. exp(x) * cos(exp(x))")
        calc.perform_rule(rules.Substitution("u", "exp(x)"))
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.ReplaceSubstitution())
        calc.perform_rule(rules.Simplify())

        calc = file.add_calculation("INT x. x * sin(3*x)")
        calc.perform_rule(rules.IntegrationByParts("x", "(-1/3) * cos(3*x)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        calc = file.add_calculation("INT x. x ^ 2 * cos(x)")
        calc.perform_rule(rules.IntegrationByParts("x^2", "sin(x)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IntegrationByParts("x", "-cos(x)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        calc = file.add_calculation("INT x. sin(x) * cos(x) * exp(sin(x))")
        calc.perform_rule(rules.Equation("sin(x) * cos(x) * exp(sin(x))", "sin(x) * exp(sin(x)) * cos(x)"))
        calc.perform_rule(rules.Substitution("u", "sin(x)"))
        calc.perform_rule(rules.IntegrationByParts("u", "exp(u)"))
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.ReplaceSubstitution())
        calc.perform_rule(rules.Simplify())

        calc = file.add_calculation("INT x. exp(x) * sin(x)")
        calc.perform_rule(rules.IntegrationByParts("exp(x)", "-cos(x)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IntegrationByParts("exp(x)", "sin(x)"))
        calc.perform_rule(rules.IntegrateByEquation("INT x. exp(x) * sin(x)"))

        calc = file.add_calculation("INT x. sin(3*x) * cos(4*x)")
        calc.perform_rule(rules.IntegrationByParts("sin(3*x)", "sin(4*x)/4"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IntegrationByParts("cos(3*x)", "-cos(4*x)/4"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IntegrateByEquation("INT x. sin(3*x) * cos(4*x)"))
        calc.print_entry()

        calc = file.add_calculation("INT x. sec(x) * sqrt(sec(x) + tan(x))")
        calc.perform_rule(rules.Equation("sec(x) * sqrt(sec(x) + tan(x))", "2 * ((sec(x) * tan(x) + sec(x) ^ 2) / (2 * sqrt(sec(x) + tan(x))))"))
        calc.perform_rule(rules.Substitution("u", "sqrt(sec(x) + tan(x))"))
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.ReplaceSubstitution())

        calc = file.add_calculation("INT x. (sin(2*x) - cos(2*x)) / (sin(2*x) + cos(2*x))")
        calc.perform_rule(rules.Equation("(sin(2*x) - cos(2*x)) / (sin(2*x) + cos(2*x))", "(-1/(2 * (sin(2*x) + cos(2*x)))) * (2 * cos(2 * x) - 2 * sin(2 * x))"))
        calc.perform_rule(rules.Substitution("u", "sin(2*x) + cos(2*x)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.ReplaceSubstitution())
        calc.perform_rule(rules.Simplify())

        calc = file.add_calculation("INT x. (sin(x) + cos(x)) / (exp(-x) + sin(x))")
        calc.perform_rule(rules.Equation("(sin(x) + cos(x)) / (exp(-x) + sin(x))", "(1 / (1 + exp(x) * sin(x))) * (cos(x) * exp(x) + exp(x) * sin(x))"))
        calc.perform_rule(rules.Substitution("u", "1 + exp(x) * sin(x)"))
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.ReplaceSubstitution())
        calc.perform_rule(rules.Simplify())

        # self.checkAndOutput(file)

    def testWallis(self):
        # Reference:
        # Irresistable Integrals, Section 2.3
        file = compstate.CompFile("base", 'wallis')

        goal = file.add_goal("(INT x:[0,oo]. 1/(x^2+b)^(m+1)) = pi / 2^(2*m+1) * binom(2*m, m) * (1/(b^((2*m+1)/2)))", conds=["b > 0", "m >= 0"])

        # Make definition
        goal.add_definition("I(m,b) = (INT x:[0,oo]. 1/(x^2+b)^(m+1))", conds=["b > 0", "m >= 0"])

        # Prove the following equality
        Eq1 = goal.add_subgoal("1", "(D b. I(m,b)) = -(m+1) * I(m+1, b)", conds=["b > 0", "m >= 0"])
        proof = Eq1.proof_by_calculation()

        calc = proof.lhs_calc
        calc.perform_rule(rules.OnSubterm(rules.ExpandDefinition("I")))
        calc.perform_rule(rules.DerivIntExchange())
        calc.perform_rule(rules.Simplify())

        calc = proof.rhs_calc
        calc.perform_rule(rules.OnSubterm(rules.ExpandDefinition("I")))
        calc.perform_rule(rules.Simplify())
        assert Eq1.is_finished()

        # Prove the following by induction
        Eq2 = goal.add_subgoal("2", "I(m, b) = pi / 2^(2*m+1) * binom(2*m, m) * (1/(b^((2*m+1)/2)))", conds=["b > 0", "m >= 0"])
        proof = Eq2.proof_by_induction("m")
        proof_base = proof.base_case.proof_by_calculation()
        proof_induct = proof.induct_case.proof_by_calculation()

        # Base case
        calc = proof_base.lhs_calc
        calc.perform_rule(rules.ExpandDefinition("I"))
        calc.perform_rule(rules.SubstitutionInverse("u", "x", parser.parse_expr("sqrt(b) * u")))
        calc.perform_rule(rules.Simplify())
        s1 = parser.parse_expr("1 / (b * u^2 + b)")
        s2 = parser.parse_expr("(1/b) * (1 / (1^2 + u^2))")
        calc.perform_rule(rules.Equation(s1, s2))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        # Induction case, LHS
        calc = proof_induct.lhs_calc
        s = calc.parse_expr("I(m + 1,b)")
        calc.perform_rule(rules.ApplyEquation("1", s))
        calc.perform_rule(rules.OnSubterm(rules.ApplyInductHyp()))
        calc.perform_rule(rules.Simplify())
        s1 = parser.parse_expr("-((2 * m + 1) / 2) - 1")
        s2 = parser.parse_expr("-m - 3/2")
        calc.perform_rule(rules.Equation(s1, s2))
        s = parser.parse_expr("b ^ (-m - 3/2) * 2 ^ -(2 * m) * pi * (2 * m + 1) / (4 * m + 4) * binom(2 * m,m)")
        calc.perform_rule(rules.Equation(None, s))

        # Induction step, RHS
        calc = proof_induct.rhs_calc
        s1 = parser.parse_expr("binom(2*m+2, m+1)")
        s2 = parser.parse_expr("2 * binom(2*m, m) * ((2*m+1) / (m+1))")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        s1 = parser.parse_expr("-((2 * m + 3) / 2)")
        s2 = parser.parse_expr("-m - 3/2")
        calc.perform_rule(rules.Equation(s1, s2))
        calc.perform_rule(rules.Simplify())
        assert Eq2.is_finished()

        proof = goal.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.OnSubterm(rules.FoldDefinition("I")))
        calc.perform_rule(rules.ApplyEquation("2", calc.parse_expr("I(m, b)")))
        goal.print_entry(is_toplevel=True)
        assert goal.is_finished()

    def testGammaFunction(self):
        # Reference:
        # Inside interesting integrals, Section 4.1
        file = compstate.CompFile("interesting", "GammaFunction")

        # Definition of Gamma function
        file.add_definition("Gamma(n) = (INT x:[0,oo]. exp(-x) * x^(n-1))", conds=["n > 0"])

        goal1 = file.add_goal("Gamma(n) = (n - 1) * Gamma(n - 1)", conds=["n > 1"])

        # Recursive equation for gamma function
        proof = goal1.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.ExpandDefinition("Gamma"))
        u = calc.parse_expr("x ^ (n - 1)")
        v = calc.parse_expr("-exp(-x)")
        calc.perform_rule(rules.IntegrationByParts(u, v))
        calc.perform_rule(rules.Simplify())
        calc = proof.rhs_calc
        calc.perform_rule(rules.OnSubterm(rules.ExpandDefinition("Gamma")))
        assert goal1.is_finished()

        # Gamma function and factorial
        goal2 = file.add_goal("Gamma(n) = factorial(n - 1)", conds=["n >= 1"])

        proof = goal2.proof_by_induction("n", 1)
        proof_base = proof.base_case.proof_by_calculation()
        proof_induct = proof.induct_case.proof_by_calculation()

        calc = proof_base.lhs_calc
        calc.perform_rule(rules.ExpandDefinition("Gamma"))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        calc = proof_induct.lhs_calc
        s = calc.parse_expr("Gamma(n + 1)")
        calc.perform_rule(rules.ApplyEquation(goal1.goal, s))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.OnSubterm(rules.ApplyInductHyp()))
        source = calc.parse_expr("n * factorial(n - 1)")
        target = calc.parse_expr("factorial(n)")
        calc.perform_rule(rules.ApplyIdentity(source, target))

        assert goal2.is_finished()

        # Application
        calc = file.add_calculation("INT x:[0,oo]. exp(-x^3)")
        calc.perform_rule(rules.Substitution('y', parser.parse_expr('x^3')))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("exp(-y) / y ^ (2/3)", "exp(-y) * y ^ (1/3 - 1)"))
        calc.perform_rule(rules.OnSubterm(rules.FoldDefinition("Gamma")))
        calc.perform_rule(rules.Equation(None, "(4/3 - 1) * Gamma(4/3 - 1)"))
        calc.perform_rule(rules.ApplyEquation(goal1.goal, calc.parse_expr("(4/3 - 1) * Gamma(4/3 - 1)")))
        goal1.print_entry(is_toplevel=True)
        goal2.print_entry(is_toplevel=True)
        calc.print_entry()

    def testChapter1Section5(self):
        # Reference:
        # Inside interesting integrals, Section 1.5
        file = compstate.CompFile("interesting", "chapter1section5")

        goal = file.add_goal("(INT x:[0,oo]. log(x) / (x ^ 2 + 1)) = 0")
        proof = goal.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.SplitRegion("1"))
        calc.perform_rule(rules.SubstitutionInverse("u", "x", "1 / u"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("u ^ 2 * (1 / u ^ 2 + 1)", "u ^ 2 + 1"))
        calc.perform_rule(rules.Simplify())

        # self.checkAndOutput(file)

    def testChapter1Section7(self):
        # Reference:
        # Inside interesting integrals, Section 1.7
        file = compstate.CompFile("interesting", "chapter1section7")

        goal = file.add_goal("(INT x:[0,1]. (x^4*(1-x)^4)/(1+x^2)) = 22/7 - pi")
        proof = goal.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Equation("(x^4*(1-x)^4)/(1+x^2)", "(x^6-4*x^5+5*x^4-4*x^2+4)-4/(1+x^2)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        # self.checkAndOutput(file)

    def testEasy01(self):
        # Reference:
        # Inside interesting integrals, Section 2.1.a
        file = compstate.CompFile("interesting", "easy01")

        goal = file.add_goal("(INT x:[1,oo]. 1 / ((x+a)*sqrt(x-1))) = pi / sqrt(a+1)", conds=["a > -1"])
        proof = goal.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Substitution(var_name='t', var_subst=parser.parse_expr("sqrt(x-1)")))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Substitution(var_name='y', var_subst=parser.parse_expr("t / sqrt(a + 1)")))
        calc.perform_rule(rules.Equation("y ^ 2 * (a + 1) + a + 1",
                                         "(a + 1) * (y^2 + 1)"))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        # self.checkAndOutput(file)

    def testEasy02(self):
        # Reference:
        # Inside interesting integrals, Section 2.1.b
        file = compstate.CompFile("interesting", "easy02")

        goal = file.add_goal("(INT x:[0, oo]. log(1 + a^2 / x^2)) = a * pi", conds=["a > 0"])
        proof = goal.proof_by_calculation()
        calc = proof.lhs_calc
        u = parser.parse_expr("log(1+a^2/x^2)")
        v = parser.parse_expr("x")
        calc.perform_rule(rules.IntegrationByParts(u = u, v = v))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("x ^ 2 * (a ^ 2 / x ^ 2 + 1)", "(a^2 + x^2)"))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        # self.checkAndOutput(file)

    def testEasy03(self):
        # Reference:
        # Inside interesting integrals, Section 2.1.c
        file = compstate.CompFile("interesting", "easy03")

        goal = file.add_goal("(INT x:[0, oo]. log(x) / (x^2+b^2)) = pi * log(b) / (2*b)", conds=["b > 0"])
        proof = goal.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.SubstitutionInverse("t", "x", parser.parse_expr("1/t")))
        old_expr = parser.parse_expr("log(1 / t)")
        new_expr = parser.parse_expr("-log(t)")
        calc.perform_rule(rules.Equation(old_expr=old_expr, new_expr=new_expr))
        # assoc rewrite because Equation Rule can't find (b ^ 2 + (1 / t) ^ 2) ^ (-1) * -(t ^ (-2))
        old_expr = parser.parse_expr("-log(t) / ((1 / t) ^ 2 + b ^ 2) * -(1 / t ^ 2)")
        new_expr = parser.parse_expr("log(t) / (1 + b^2*t^2)")
        calc.perform_rule(rules.Equation(old_expr=old_expr, new_expr=new_expr))

        calc.perform_rule(rules.SubstitutionInverse('s', "t", parser.parse_expr("s/b")))
        calc.perform_rule(rules.Simplify())

        old_expr = parser.parse_expr("log(s/b)")
        new_expr = parser.parse_expr("log(s) - log(b)")
        calc.perform_rule(rules.Equation(old_expr=old_expr, new_expr=new_expr))

        calc.perform_rule(rules.ExpandPolynomial())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        # self.checkAndOutput(file)

    def testEasy04(self):
        # Reference:
        # Inside interesting integrals, Section 2.1.d
        file = compstate.CompFile("interesting", "easy04")

        goal = file.add_goal("(INT x:[1, oo]. log(x)/(x+1)^2) = log(2)")

        goal01 = goal.add_subgoal("1", "(INT x:[0, oo]. 1/(1 + exp(a*x))) = (log(2)/a)", conds=["a > 0"])
        proof = goal01.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Substitution(var_name="u", var_subst=parser.parse_expr("exp(a*x)")))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("1 / (u * (u + 1))", "1/u - 1/(u+1)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Substitution(var_name="y", var_subst=parser.parse_expr("u+1")))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        goal02 = goal.add_subgoal("2", "(INT x:[1,oo]. log(x) / (a ^ 2 * (x + 1) ^ 2))=log(2)/a^2", conds=["a > 0"])
        proof = goal02.proof_by_rewrite_goal(begin="1")
        calc = proof.begin
        calc.perform_rule(rules.DerivEquation("a"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Substitution(var_name="y", var_subst="exp(a*x)"))
        calc.perform_rule(rules.SolveEquation("INT y:[1,oo]. log(y) / (a ^ 2 * (y + 1) ^ 2)"))

        proof = goal.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Equation("(x+1)^2", "1^2 * (x+1)^2"))
        s = calc.parse_expr("INT x:[1,oo]. log(x) / (1 ^ 2 * (x + 1) ^ 2)")
        calc.perform_rule(rules.ApplyEquation("2", s))
        calc.perform_rule(rules.Simplify())
        assert goal.is_finished()

    def testEasy05(self):
        # Reference:
        # Inside interesting integrals, Section 2.1.e
        file = compstate.CompFile("interesting", "easy05")

        goal = file.add_goal("(INT x:[sqrt(2), oo]. 1/(x+x^sqrt(2)))=(1+sqrt(2))*log(1+2^(1/2*(1-sqrt(2))))")
        proof = goal.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Equation("1/(x+x^sqrt(2))", "x^(-sqrt(2))/(x^(1-sqrt(2))+1)"))
        calc.perform_rule(rules.Substitution(var_name="u", var_subst="log(x^(1-sqrt(2))+1)"))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.OnCount(rules.Equation("-sqrt(2) + 1", "-1/(1+sqrt(2))"), 2))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.OnCount(rules.Equation("sqrt(2)", "2^(1/2)"), 2))
        calc.perform_rule(rules.ApplyIdentity("2 ^ (1/2) ^ (-sqrt(2) + 1)", "2^(1/2*(-sqrt(2)+1))"))
        calc.perform_rule(rules.Equation("(sqrt(2) + 1) * log(2 ^ (1/2 * (-sqrt(2) + 1)) + 1)",
                                         "(1+sqrt(2))*log(1+2^(1/2*(1-sqrt(2))))"))
        assert goal.is_finished()

    def testEasy06(self):
        # Reference:
        # Inside interesting integrals, Section 2.1.f
        file = compstate.CompFile("interesting", "easy06")
        goal = file.add_goal("(INT x:[-oo, oo]. 1/cosh(x)) = pi")
        proof = goal.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.OnSubterm(rules.ExpandDefinition("cosh")))
        calc.perform_rule(rules.Substitution("t", parser.parse_expr("exp(x)")))
        calc.perform_rule(rules.Equation("-log(t)", "log(1/t)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("1 / (t * (1 / t + t))", "1 / (1 + t ^ 2)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        # self.checkAndOutput(file)

    def testTrick2a(self):
        # Reference:
        # Inside interesting integrals, Section 2.2, example 1
        file = compstate.CompFile("interesting", "Trick2a")

        calc = file.add_calculation("INT x:[0,pi/2]. sqrt(sin(x)) / (sqrt(sin(x)) + sqrt(cos(x)))")
        calc.perform_rule(rules.Substitution("y", parser.parse_expr("pi / 2 - x")))
        calc.perform_rule(rules.Equation(
            "sqrt(cos(y)) / (sqrt(cos(y)) + sqrt(sin(y)))",
            "1 - sqrt(sin(y)) / (sqrt(cos(y)) + sqrt(sin(y)))"))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.IntegrateByEquation(
            parser.parse_expr("INT x:[0,pi/2]. sqrt(sin(x)) / (sqrt(sin(x)) + sqrt(cos(x)))")))
        self.assertEqual(str(calc.last_expr), "pi / 4")

        # self.checkAndOutput(file)

    def testTrick2b(self):
        # Reference:
        # Inside interesting integrals, Section 2.2, example 2
        file = compstate.CompFile("interesting", "Trick2b")

        calc = file.add_calculation("INT x:[0,pi]. x * sin(x) / (1 + cos(x)^2)")
        calc.perform_rule(rules.Substitution("y", parser.parse_expr("pi - x")))
        calc.perform_rule(rules.ExpandPolynomial())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IntegrateByEquation(parser.parse_expr("INT x:[0,pi]. x * sin(x) / (1 + cos(x)^2)")))
        calc.perform_rule(rules.Substitution("u", parser.parse_expr("cos(y)")))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertEqual(str(calc.last_expr), "pi ^ 2 / 4")

        # self.checkAndOutput(file)

    def testTrick2c(self):
        # Reference:
        # Inside interesting integrals, Section 2.2, example 3
        file = compstate.CompFile("interesting", "Trick2c")

        goal = file.add_goal("(INT x:[0, pi/2]. sin(x) ^ 2 / (sin(x) + cos(x))) = " + \
                             "sqrt(2) / 4 * log(3 + 2*sqrt(2))")

        goal01 = goal.add_subgoal("1", "(INT x:[0, pi/2]. sin(x) ^ 2 / (sin(x) + cos(x))) = " + \
                                  "(INT x:[0, pi/2]. cos(x) ^ 2 / (sin(x) + cos(x)))")
        proof = goal01.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Substitution(var_name='y', \
                          var_subst=parser.parse_expr("pi/2 - x")))

        goal02 = goal.add_subgoal("2", "(INT x:[0, pi/2]. sin(x) ^ 2 / (sin(x) + cos(x))) " + \
                               "= (1/2 * INT x:[0, pi/2]. 1 / (sin(x) + cos(x)))")
        proof = goal02.proof_by_calculation()
        calc = proof.rhs_calc
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("1", "sin(x) ^ 2 + cos(x) ^ 2"))
        calc.perform_rule(rules.ExpandPolynomial())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.OnCount(rules.Equation("cos(x) + sin(x)", "sin(x) + cos(x)"), 1))
        s = calc.parse_expr("INT x:[0,pi / 2]. cos(x) ^ 2 / (sin(x) + cos(x))")
        calc.perform_rule(rules.ApplyEquation("1", s))
        calc.perform_rule(rules.Simplify())

        proof = goal.proof_by_calculation()
        calc = proof.lhs_calc
        s = calc.parse_expr("INT x:[0,pi / 2]. sin(x) ^ 2 / (sin(x) + cos(x))")
        calc.perform_rule(rules.ApplyEquation("2", s))
        calc.perform_rule(rules.Substitution('z', 'tan(x/2)'))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("(-(z ^ 2) + 1) / (z ^ 2 + 1) + 2 * z / (z ^ 2 + 1)",
                                         "(2 - (z-1)^2) / (z^2 + 1)"))
        calc.perform_rule(rules.Equation("(z ^ 2 + 1) * ((2 - (z - 1) ^ 2) / (z ^ 2 + 1))",
                                         "2 - (z - 1) ^ 2"))
        calc.perform_rule(rules.Equation("2 - (z - 1) ^ 2",
                                         "(sqrt(2) + (z-1))*(sqrt(2) - (z-1))"))
        calc.perform_rule(rules.Equation("1 / ((sqrt(2) + (z - 1)) * (sqrt(2) - (z - 1)))",
                                         "sqrt(2) / 4 * (1/(sqrt(2) + (z - 1)) + 1/(sqrt(2) - (z - 1)))"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.OnCount(rules.Substitution("u", "sqrt(2) + 1 - z"), 1))
        calc.perform_rule(rules.OnCount(rules.Substitution("u", "sqrt(2) - 1 + z"), 2))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("sqrt(2) * (log(sqrt(2) + 1) - log(sqrt(2) - 1)) / 4",
                                         "1/4 * sqrt(2) * (log(sqrt(2) + 1) - log(sqrt(2) - 1))"))
        calc.perform_rule(rules.ApplyIdentity("log(sqrt(2) + 1) - log(sqrt(2) - 1)",
                                              "log((sqrt(2) + 1) / (sqrt(2) - 1))"))
        calc.perform_rule(rules.Equation("(sqrt(2) + 1) / (sqrt(2) - 1)", \
                                         "3 + 2 * sqrt(2)"))
        calc.perform_rule(rules.Simplify())
        assert goal.is_finished()

    def testTrick2d(self):
        # Reference:
        # Inside interesting integrals, Section 2.2, example 4
        file = compstate.CompFile("interesting", "Trick2d")

        goal = file.add_goal("(INT x:[0,1]. log(x+1) / (x^2 + 1)) = pi / 8 * log(2)")

        goal01 = goal.add_subgoal("1", "(INT x:[0,1]. log(x+1) / (x^2 + 1)) = (INT x:[0,pi / 4]. log(tan(x) + 1))")
        proof = goal01.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.SubstitutionInverse("u", "x", "tan(u)"))
        calc.perform_rule(rules.ApplyIdentity("sec(u) ^ 2" , "tan(u)^2 + 1"))
        calc.perform_rule(rules.Simplify())

        goal02 = goal.add_subgoal("2", "(INT x:[0,1]. log(x+1) / (x^2 + 1)) "\
                                       +"= (pi / 4 * log(2) - (INT x:[0,1]. log(x+1) / (x^2 + 1)))")
        proof = goal02.proof_by_calculation()
        calc = proof.lhs_calc
        s = calc.parse_expr("INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)")
        calc.perform_rule(rules.ApplyEquation("1", s))
        calc.perform_rule(rules.SubstitutionInverse("y", "x", "pi/4 - y"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.ApplyIdentity("tan(pi / 4 - y)",\
                                              "(tan(pi / 4) - tan(y)) / (1 + tan(pi / 4) * tan(y))"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("(-tan(y) + 1) / (tan(y) + 1) + 1",\
                                         "2 / (1 + tan(y))"))
        calc.perform_rule(rules.ApplyIdentity("log(2 / (1 + tan(y)))", \
                                              "log(2) - log(1 + tan(y))"))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        s = calc.parse_expr("INT x:[0,pi / 4]. log(tan(x) + 1)")
        calc.perform_rule(rules.ApplyEquation("1", s))

        proof = goal.proof_by_rewrite_goal(begin="2")
        calc = proof.begin
        calc.perform_rule(rules.SolveEquation("(INT x:[0,1]. log(x+1) / (x^2 + 1))"))
        assert goal.is_finished()

    def testTrick2e(self):
        # Reference:
        # Inside interesting integrals, Section 2.2, example 5
        file = compstate.CompFile("interesting", "Trick2e")

        goal = file.add_goal("(INT t:[0,a]. log(t+a) / (t^2 + a^2)) = pi /(8*a) * log(2*a^2)",
                               conds=["a > 0"])

        goal01 = goal.add_subgoal("1", "(INT x:[0,1]. log(x+1) / (x^2+1)) = (a * INT t:[0,a]. log(t+a) / (t^2 + a^2)) - pi / 4 * log(a)",
                                  conds=["a > 0"])
        proof = goal01.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.SubstitutionInverse("t", "x", "t/a"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("1 / (t ^ 2 / a ^ 2 + 1) * log(t / a + 1)",
                                         "log(t / a + 1) * a^2 / (t ^ 2  + a ^ 2)"))
        calc.perform_rule(rules.Equation("t / a + 1", "(t + a) / a"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.ApplyIdentity("log((a + t) / a)", "log(a+t) - log(a)"))
        calc.perform_rule(rules.Equation("1 / (a ^ 2 + t ^ 2) * (log(a + t) - log(a))",
                                         "log(a + t) / (a ^ 2 + t ^ 2) - log(a) / (a ^ 2 + t ^ 2)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.ExpandPolynomial())

        goal02 = goal.add_subgoal("2", "((a * INT t:[0,a]. log(t+a) / (t^2 + a^2)) - pi / 4 * log(a)) = pi * log(2) / 8",
                               conds=["a > 0"])
        proof = goal02.proof_by_calculation()
        calc = proof.lhs_calc
        s = calc.parse_expr("(a * INT t:[0,a]. log(t+a) / (t^2 + a^2)) - pi / 4 * log(a)")
        calc.perform_rule(rules.ApplyEquation("1", s))
        calc.perform_rule(rules.DefiniteIntegralIdentity())

        proof = goal.proof_by_rewrite_goal(begin="2")
        calc = proof.begin
        calc.perform_rule(rules.SolveEquation("INT t:[0,a]. log(t+a) / (t^2 + a^2)"))
        calc.perform_rule(rules.Equation("pi * log(a) / 4", "1/8 * pi * (2 * log(a))"))
        calc.perform_rule(rules.Equation("2 * log(a)", "log(a^2)"))
        calc.perform_rule(rules.Equation("1/8 * pi * log(a ^ 2) + pi * log(2) / 8", "1/8 * pi * (log(2) + log(a^2))"))
        calc.perform_rule(rules.Equation("(log(2) + log(a ^ 2))", "log(2 * a^2)"))
        calc.perform_rule(rules.Equation("1 / a * (1/8 * pi * log(2 * a ^ 2))", "pi / (8 * a) * log(2 * a ^ 2)"))
        assert goal.is_finished()

    def testPartialFraction(self):
        # Reference
        # Inside interesting integrals, Section 2.3, example 2
        file = compstate.CompFile("interesting", 'partialFraction')

        goal = file.add_goal("(INT x:[0,oo]. 1 / (x^4 + 2*x^2*cosh(2*a) + 1)) = pi / (4 * cosh(a))")
        proof = goal.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.OnSubterm(rules.ExpandDefinition("cosh")))
        calc.perform_rule(rules.Equation(
            "x ^ 4 + 2 * x ^ 2 * ((exp(-(2 * a)) + exp(2 * a)) / 2) + 1",
            "(x ^ 2 + exp(2 * a)) * (x ^ 2 + exp(-(2 * a)))"))
        calc.perform_rule(rules.Equation(
            "1 / ((x ^ 2 + exp(2 * a)) * (x ^ 2 + exp(-(2 * a))))",
            "1 / (exp(2*a) - exp(-(2*a))) * (1 / (x^2 + exp(-(2 * a))) - 1 / (x^2 + exp(2*a)))"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("exp(-(2 * a))", "exp(-a) ^ 2"))
        calc.perform_rule(rules.Equation("exp(-(2 * a))", "exp(-a) ^ 2"))
        calc.perform_rule(rules.Equation("exp(2 * a)", "exp(a) ^ 2"))
        calc.perform_rule(rules.Equation("exp(2 * a)", "exp(a) ^ 2"))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation(None, "pi / (4 * ((exp(a) + exp(-a)) / 2))"))
        calc.perform_rule(rules.OnSubterm(rules.FoldDefinition("cosh")))

        # self.checkAndOutput(file)

    def testPartialFraction03(self):
        # Reference
        # Inside interesting integrals, Section 2.3, example 3
        file = compstate.CompFile("interesting", 'partialFraction03')

        file.add_definition("I(a) = (INT x:[0,oo]. 1 / (x^4 + 2*x^2*cos(2*a) + 1))")

        goal = file.add_goal("2*x^2*cos(2*a) + x^4 + 1 != 0", conds=["cos(a) != 0"])
        cond = parser.parse_expr("x != 0")
        proof = goal.proof_by_case(cond)
        proofa = proof.cases[0].proof_by_calculation()
        calc = proofa.lhs_calc
        calc.perform_rule(rules.Equation(None, "(x^2 - 1) ^ 2 + 2*x^2*(1+cos(2*a))"))
        calc.perform_rule(rules.ApplyIdentity("cos(2*a)", "2 * cos(a)^2 - 1"))
        calc.perform_rule(rules.Simplify())
        proofb = proof.cases[1].proof_by_calculation()
        calc = proofb.lhs_calc
        calc.perform_rule(rules.Simplify())
        self.assertTrue(proof.is_finished())

        goal = file.add_goal("x ^ 4 + 2 * x ^ 2 * cos(2 * a) + 1 != 0", conds=["cos(a) != 0"])
        proof = goal.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Simplify())
        self.assertTrue(proof.is_finished())

        goal = file.add_goal("(x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1) != 0", conds=["cos(a) != 0"])
        cond = parser.parse_expr("x != 0")
        proof = goal.proof_by_case(cond)
        proofa = proof.cases[0].proof_by_calculation()
        calc = proofa.lhs_calc
        calc.perform_rule(rules.ExpandPolynomial())
        calc.perform_rule(rules.ApplyIdentity("sin(a)^2", "1 - cos(a)^2"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation(None, "(x^2 - 1) ^ 2 + 4*x^2*(cos(a)^2)"))
        proofb = proof.cases[1].proof_by_calculation()
        calc = proofb.lhs_calc
        calc.perform_rule(rules.Simplify())
        self.assertTrue(proof.is_finished())

        goal = file.add_goal("-(2 * x * sin(a)) + x^2 + 1 > 0", conds=["cos(a) != 0"])
        cond = parser.parse_expr("x")
        proof = goal.proof_by_case(cond)
        proofa = proof.cases[0].proof_by_calculation()
        calc = proofa.lhs_calc
        calc.perform_rule(rules.Equation(None, "(x + 1) ^ 2 - 2 * x * (1 + sin(a))"))
        proofb = proof.cases[1].proof_by_calculation()
        calc = proofb.lhs_calc
        calc.perform_rule(rules.Simplify())
        proofc = proof.cases[2].proof_by_calculation()
        calc = proofc.lhs_calc
        calc.perform_rule(rules.Equation(None, "(x - 1) ^ 2 + 2 * x * (1 - sin(a))"))
        self.assertTrue(proof.is_finished())

        goal01 = file.add_goal("I(a) = (INT x:[0,oo]. x^2 / (x^4 + 2*x^2*cos(2*a) + 1))", conds=["cos(a) != 0"])
        proof = goal01.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.ExpandDefinition("I"))
        calc.perform_rule(rules.SubstitutionInverse("y", "x", "1/y"))
        calc.perform_rule(rules.Equation("1 / (2 * (1 / y) ^ 2 * cos(2 * a) + (1 / y) ^ 4 + 1) * - (1 / y ^ 2)", \
                                         "-y^2 / (y^4 + 2*y^2*cos(2*a) + 1)"))
        calc.perform_rule(rules.Simplify())
        calc = proof.rhs_calc
        calc.perform_rule(rules.Simplify())

        goal02 = file.add_goal("2*I(a) = (INT x:[0,oo]. (1+x^2) / (x^4 + 2*x^2*cos(2*a) + 1))", conds=["cos(a) != 0"])
        proof = goal02.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Equation("2*I(a)", "I(a) + I(a)"))
        calc.perform_rule(rules.OnLocation(rules.ExpandDefinition("I"), "0"))
        s = calc.parse_expr("I(a)")
        calc.perform_rule(rules.ApplyEquation(goal01.goal, s))
        calc.perform_rule(rules.Equation("(INT x:[0,oo]. 1 / (2 * x ^ 2 * cos(2 * a) + x ^ 4 + 1)) + (INT x:[0,oo]. x ^ 2 / (x ^ 4 + 2 * x ^ 2 * cos(2 * a) + 1))", \
                                         "(INT x:[0,oo]. 1 / (2 * x ^ 2 * cos(2 * a) + x ^ 4 + 1) + x ^ 2 / (x ^ 4 + 2 * x ^ 2 * cos(2 * a) + 1))"))
        calc.perform_rule(
            rules.Equation("1 / (2 * x ^ 2 * cos(2 * a) + x ^ 4 + 1) + x ^ 2 / (x ^ 4 + 2 * x ^ 2 * cos(2 * a) + 1)",
                           "(1+x^2) / (2 * x ^ 2 * cos(2 * a) + x ^ 4 + 1)"))
        calc.perform_rule(rules.Simplify())
        calc = proof.rhs_calc
        calc.perform_rule(rules.Simplify())

        goal03a = file.add_goal("(INT x:[-oo,oo]. (x ^ 2 + 1) / (2 * x ^ 2 * cos(2 * a) + x ^ 4 + 1)) = " \
                                "2 * (INT x:[0,oo]. (x ^ 2 + 1) / (2 * x ^ 2 * cos(2 * a) + x ^ 4 + 1))", conds=["cos(a) != 0"])
        proof = goal03a.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.SplitRegion("0"))
        calc.perform_rule(rules.Substitution("x", "-x"))
        calc.perform_rule(rules.Simplify())

        goal03b = file.add_goal("(INT x:[-oo,oo]. (-2*x*sin(a)) / ((x^2-2*x*sin(a)+1)*(x^2+2*x*sin(a)+1))) = 0", conds=["cos(a) != 0"])
        proof = goal03b.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.SplitRegion("0"))
        calc.perform_rule(rules.Substitution("x", "-x"))
        calc.perform_rule(rules.Simplify())

        goal03 = file.add_goal("I(a) = (1/4 * (INT x:[-oo,oo]. 1 / (cos(a) ^ 2 + x ^ 2)))", conds=["cos(a) != 0"])
        proof = goal03.proof_by_rewrite_goal(begin=goal02)
        calc = proof.begin
        calc.perform_rule(rules.SolveEquation("I(a)"))
        s = calc.parse_expr("INT x:[0,oo]. (x ^ 2 + 1) / (2 * x ^ 2 * cos(2 * a) + x ^ 4 + 1)")
        calc.perform_rule(rules.ApplyEquation(goal03a.goal, s))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.ApplyIdentity("cos(2*a)", "1-2*sin(a)^2"))
        calc.perform_rule(rules.Equation("2*x^2*(1-2*sin(a)^2) + x^4 + 1",
                                         "(x^2-2*x*sin(a)+1)*(x^2+2*x*sin(a)+1)"))
        calc.perform_rule(rules.Equation("1/4 * (INT x:[-oo,oo]. (x ^ 2 + 1) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 "
                                         "+ 2 * x * sin(a) + 1)))",
                                         "1/4 * (INT x:[-oo,oo]. (x ^ 2 + 1) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 "
                                         "+ 2 * x * sin(a) + 1))) + 1/4 * 0"))
        s = calc.parse_expr("0")
        calc.perform_rule(rules.ApplyEquation(goal03b.goal, s))
        calc.perform_rule(rules.Equation("1/4 * (INT x:[-oo,oo]. (x ^ 2 + 1) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1))) + 1/4 * (INT x:[-oo,oo]. -2 * x * sin(a) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1)))",
                                         "1/4 * ((INT x:[-oo,oo]. (x ^ 2 + 1) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1))) + (INT x:[-oo,oo]. -2 * x * sin(a) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1))))"))
        calc.perform_rule(rules.Equation("(INT x:[-oo,oo]. (x ^ 2 + 1) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1))) + (INT x:[-oo,oo]. -2 * x * sin(a) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1)))",
                                         "(INT x:[-oo,oo]. (x ^ 2 + 1) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1)) - 2 * x * sin(a) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1)))"))
        calc.perform_rule(rules.Equation("(x ^ 2 + 1) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1)) - 2 * x * sin(a) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1))",
                                         "(x^2+1-2*x*sin(a)) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1))"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.OnLocation(rules.Equation("1", "sin(a)^2 +cos(a)^2"), "1.1.0.1.1"))
        calc.perform_rule(rules.Equation("(2 * x * sin(a) + x ^ 2 + (sin(a) ^ 2 + cos(a) ^ 2))",
                                         "(x+sin(a))^2 + cos(a) ^ 2"))
        calc.perform_rule(rules.Substitution("u", "x+sin(a)"))

        goal04 = file.add_goal("I(a) = pi / (4 * cos(a))", conds=["cos(a)>0"])
        proof = goal04.proof_by_calculation()
        calc = proof.lhs_calc
        s = calc.parse_expr("I(a)")
        calc.perform_rule(rules.ApplyEquation(goal03.goal, s))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        goal05 = file.add_goal("I(a) = -(pi / (4 * cos(a)))", conds=["cos(a)<0"])
        proof = goal05.proof_by_calculation()
        calc = proof.lhs_calc
        s = calc.parse_expr("I(a)")
        calc.perform_rule(rules.ApplyEquation(goal03.goal, s))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        goal06 = file.add_goal("(INT x:[0,oo]. 1/(x^4+1))=(pi*sqrt(2))/4")
        proof = goal06.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Equation("x^4+1", "x^4 + 2*x^2*cos(2*(pi/4))+1"))
        calc.perform_rule(rules.FoldDefinition("I"))
        s = calc.parse_expr("I(pi / 4)")
        calc.perform_rule(rules.ApplyEquation(goal04.goal, s))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("sqrt(2) * pi / 4", "pi * sqrt(2) / 4"))

        goal07 = file.add_goal("(INT x:[0,oo]. 1/(x^4+1)) = (INT x:[0,oo]. x^2/(x^4+1))")
        proof = goal07.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Substitution(var_name='u', var_subst="1/x"))
        calc.perform_rule(rules.Equation("1 / (u ^ 2 * (1 / u ^ 4 + 1))", "u^2/(u^4+1)"))

        goal08 = file.add_goal("(INT x:[0,oo]. 1/(x^4+x^2+1))=pi/(2*sqrt(3))")
        proof = goal08.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Equation("x^4+x^2+1", "x^4 + 2*x^2*cos(2*(pi/6))+1"))
        calc.perform_rule(rules.FoldDefinition("I"))
        s = calc.parse_expr("I(pi / 6)")
        calc.perform_rule(rules.ApplyEquation(goal04.goal, s))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("sqrt(3) * pi / 6", "pi/(2*sqrt(3))"))

        goal09 = file.add_goal("(INT x:[0,oo]. 1/(x^4-x^2+1))=pi/2")
        proof = goal09.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Equation("x^4-x^2+1", "x^4 + 2*x^2*cos(2*(pi/3))+1"))
        calc.perform_rule(rules.FoldDefinition("I"))
        s = calc.parse_expr("I(pi / 3)")
        calc.perform_rule(rules.ApplyEquation(goal04.goal, s))
        calc.perform_rule(rules.Simplify())

        goal10 = file.add_goal("(INT x:[0,oo]. 1/(x^4+2*x^2+1))=pi/4")
        proof = goal10.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Equation("x^4+2*x^2+1", "x^4 + 2*x^2*cos(2*0)+1"))
        calc.perform_rule(rules.FoldDefinition("I"))
        s = calc.parse_expr("I(0)")
        calc.perform_rule(rules.ApplyEquation(goal04.goal, s))
        calc.perform_rule(rules.Simplify())

        # self.checkAndOutput(file)

    def testLeibniz01(self):
        # Reference
        # Inside interesting integrals, Section 3.1, example 1
        file = compstate.CompFile("interesting", "Leibniz01")

        goal = file.add_goal("(INT x:[0,oo]. 1 / (x^2 + a^2)^3) = 3*pi / (16 * a^5)", conds=["a > 0"])

        # Basic result: integral of 1 / (x^2 + a^2)
        goal1 = goal.add_subgoal("1", "(INT x:[0,oo]. 1 / (x^2 + a^2)) = pi / (2 * a)", conds=["a > 0"])
        proof = goal1.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.SubstitutionInverse("u", "x", parser.parse_expr("a * u")))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("1 / (a ^ 2 * u ^ 2 + a ^ 2)", "1 / ((a ^ 2) * (u^2 + 1))"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        # Derivate to get integral of 1 / (x^2 + a^2)^2
        goal2 = goal.add_subgoal("2", "(INT x:[0,oo]. 1 / (x^2 + a^2)^2) = pi / (4 * a^3)", conds=["a > 0"])
        proof = goal2.proof_by_rewrite_goal(begin="1")
        calc = proof.begin
        calc.perform_rule(rules.DerivEquation('a'))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.SolveEquation(parser.parse_expr("INT x:[0,oo]. 1 / (a ^ 2 + x ^ 2) ^ 2")))

        # Derivate again:
        proof = goal.proof_by_rewrite_goal(begin="2")
        calc = proof.begin
        calc.perform_rule(rules.DerivEquation('a'))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.SolveEquation(parser.parse_expr("INT x:[0,oo]. 1 / (a ^ 2 + x ^ 2) ^ 3")))
        assert goal.is_finished()

    def testLeibniz02(self):
        # Reference
        # Inside interesting integrals, Section 3.1, example 2
        file = compstate.CompFile("interesting", 'leibniz02')

        # Overall goal: (INT x:[-oo,oo]. exp(-(x^2)/2)) = sqrt(2*pi)
        goal = file.add_goal("(INT x:[0,1]. 1/sqrt(-log(x))) = sqrt(pi)")

        # Make definition
        goal.add_definition("g(t) = (INT x:[0,t].exp(-(x^2)/2))^2")

        Eq1 = goal.add_subgoal("1", "(INT x:[-oo,oo]. exp(-x^2/2)) = 2 * LIM {t->oo}. sqrt(g(t))")
        Eq1_proof = Eq1.proof_by_calculation()
        calc = Eq1_proof.lhs_calc
        calc.perform_rule(rules.SplitRegion(expr.Const(0)))
        calc.perform_rule(rules.Substitution('y', parser.parse_expr("-x")))
        calc.perform_rule(rules.Substitution('x', parser.parse_expr("y")))
        calc.perform_rule(rules.Simplify())
        calc = Eq1_proof.rhs_calc
        calc.perform_rule(rules.OnSubterm(rules.ExpandDefinition("g")))
        calc.perform_rule(rules.Simplify())

        Eq2 = goal.add_subgoal("2", "(D t. g(t) + 2 * INT y:[0, 1].exp(-(1+y^2)*t^2/2)/(1+y^2)) = 0", conds=["t > 0"])
        Eq2_proof = Eq2.proof_by_calculation()
        calc = Eq2_proof.lhs_calc
        calc.perform_rule(rules.OnSubterm(rules.ExpandDefinition("g")))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.OnCount(rules.Substitution('y', parser.parse_expr('x/t')), 2))
        calc.perform_rule(rules.Equation("exp(t ^ 2 * (-(y ^ 2) - 1) / 2)", "exp(1/2 * t ^ 2 * (-(y ^ 2) - 1))"))
        calc.perform_rule(rules.Equation("1/2 * t ^ 2 * (-(y ^ 2) - 1)", "-1/2 * t ^ 2 * y ^ 2 + 1/2 * t ^ 2 * (-1)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.ApplyIdentity(
            "exp(-(t ^ 2 * y ^ 2 / 2) - t ^ 2 / 2)",
            "exp(-1/2 * t ^ 2 * y ^ 2) * exp(-1/2 * t ^ 2)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("(-(y ^ 2) - 1) / (y ^ 2 + 1)", "-1"))
        calc.perform_rule(rules.Simplify())
        assert Eq2.is_finished()

        Eq3 = goal.add_subgoal("3", "2 * (INT y:[0,1]. exp(1/2 * t ^ 2 * (-(y ^ 2) - 1)) * (y ^ 2 + 1) ^ -1) + g(t) = SKOLEM_CONST(C)", conds=['t>0'])
        Eq3_proof = Eq3.proof_by_rewrite_goal(begin="2")
        calc = Eq3_proof.begin
        calc.perform_rule(rules.IntegralEquation())
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        assert Eq3.is_finished()

        Eq4 = goal.add_subgoal("4", "pi/2 = SKOLEM_CONST(C)")
        proof_of_Eq4 = Eq4.proof_by_rewrite_goal(begin="3")
        calc = proof_of_Eq4.begin
        calc.perform_rule(rules.LimitEquation('t', expr.Const(0)))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.OnSubterm(rules.ExpandDefinition("g")))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        assert Eq4.is_finished()

        Eq5 = goal.add_subgoal("5", "g(t) = -(2 * (INT y:[0,1]. 1 / (y ^ 2 + 1) * exp(t ^ 2 * (-(y ^ 2) - 1) / 2))) + pi / 2", conds=['t>0'])
        proof_of_Eq5 = Eq5.proof_by_calculation()
        calc = proof_of_Eq5.lhs_calc
        s = calc.parse_expr("g(t)")
        calc.perform_rule(rules.ApplyEquation("3", s))
        s = calc.parse_expr("SKOLEM_CONST(C)")
        calc.perform_rule(rules.ApplyEquation("4", s))
        calc.perform_rule(rules.Simplify())
        assert Eq5.is_finished()

        Eq6 = goal.add_subgoal("6", "(INT x:[-oo,oo]. exp(-x^2/2)) = sqrt(2*pi)")
        proof_of_Eq6 = Eq6.proof_by_calculation()
        calc = proof_of_Eq6.lhs_calc
        s = calc.parse_expr("INT x:[-oo,oo]. exp(-(x ^ 2) / 2)")
        calc.perform_rule(rules.ApplyEquation("1", s))
        s = calc.parse_expr("g(t)")
        calc.perform_rule(rules.ApplyEquation("5", s))
        calc.perform_rule(rules.Simplify())
        calc = proof_of_Eq6.rhs_calc
        calc.perform_rule(rules.Simplify())
        assert Eq6.is_finished()

        Eq7 = goal.add_subgoal("7", "(INT x:[0,oo]. exp(-x^2/2)) = sqrt(2) * sqrt(pi) / 2")
        proof_of_Eq7 = Eq7.proof_by_rewrite_goal(begin="6")
        calc = proof_of_Eq7.begin
        calc.perform_rule(rules.SplitRegion(expr.Const(0)))
        calc.perform_rule(rules.Substitution('y', parser.parse_expr("-x")))
        calc.perform_rule(rules.Substitution('x', parser.parse_expr("y")))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.SolveEquation(parser.parse_expr("INT x:[0,oo]. exp(-(x ^ 2 / 2))")))
        assert Eq7.is_finished()

        Eq8 = goal.add_subgoal("8", "(INT x:[-oo,oo]. exp(-(a*x^2))) = sqrt(pi / a)", conds=["a > 0"])
        proof_of_Eq8 = Eq8.proof_by_calculation()
        calc = proof_of_Eq8.lhs_calc
        calc.perform_rule(rules.Substitution("u", parser.parse_expr("sqrt(2*a) * x")))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Substitution('x', parser.parse_expr("u")))
        calc.perform_rule(rules.Equation("-(x ^ 2 / 2)", "-x^2 / 2"))
        s = calc.parse_expr("INT x:[-oo,oo]. exp(-(x ^ 2) / 2)")
        calc.perform_rule(rules.ApplyEquation("6", s))
        calc.perform_rule(rules.Simplify())

        calc = proof_of_Eq8.rhs_calc
        calc.perform_rule(rules.Simplify())
        assert Eq8.is_finished()

        Eq9 = goal.add_subgoal("9", "(INT x:[0,oo]. exp(-x^2)) = sqrt(pi)/2")
        proof_of_Eq9 = Eq9.proof_by_rewrite_goal(begin="7")
        calc = proof_of_Eq9.begin
        calc.perform_rule(rules.Substitution(var_name="x", var_subst="x/sqrt(2)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.SolveEquation("(INT x:[0,oo]. exp(-x^2))"))
        assert Eq9.is_finished()

        proof_of_Eq10 = goal.proof_by_rewrite_goal(begin="9")
        calc = proof_of_Eq10.begin
        calc.perform_rule(rules.Substitution(var_name="t", var_subst="exp(-x^2)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.SolveEquation("(INT x:[0,1]. 1/sqrt(-log(x)))"))
        goal.print_entry()
        assert goal.is_finished()

    def testLeibniz03(self):
        # Reference:
        # Inside interesting integrals, Section 3.1, example #3

        # Overall goal: INT x:[0,oo]. cos(tx)*exp(-(x^2)/2) = sqrt(pi/2)*exp(-(t^2)/2)
        # TODO: remove conditions I(t) > 0

        # Initial state
        file = compstate.CompFile("interesting", 'leibniz03')

        # Make definition
        file.add_definition("I(t) = INT x:[0,oo]. cos(t*x)*exp(-(x^2)/2)")

        Eq0 = file.add_goal("I(0) = sqrt(pi/2)")
        Eq0_proof = Eq0.proof_by_calculation()
        calc = Eq0_proof.lhs_calc
        calc.perform_rule(rules.ExpandDefinition("I"))
        calc.perform_rule(rules.Equation("-(x ^ 2 / 2)", "-(x^2)/2"))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc = Eq0_proof.rhs_calc
        calc.perform_rule(rules.Simplify())
        assert Eq0_proof.is_finished()
        # Prove the following equality
        Eq1 = file.add_goal("(D t. I(t)) = -t*I(t)")
        Eq1_proof = Eq1.proof_by_calculation()
        calc = Eq1_proof.lhs_calc
        calc.perform_rule(rules.OnSubterm(rules.ExpandDefinition("I")))
        calc.perform_rule(rules.Simplify())
        u = parser.parse_expr('sin(t*x)')
        v = parser.parse_expr('-exp(-x^2/2)')
        calc.perform_rule(rules.IntegrationByParts(u, v))
        calc.perform_rule(rules.Simplify())
        calc = Eq1_proof.rhs_calc
        calc.perform_rule(rules.OnSubterm(rules.ExpandDefinition("I")))
        calc.perform_rule(rules.Simplify())
        assert Eq1.is_finished()
        Eq2 = file.add_goal("(D t. log(I(t)) + t^2/2) = 0", conds=["I(t) > 0"])
        Eq2_proof = Eq2.proof_by_calculation()
        calc = Eq2_proof.lhs_calc
        calc.perform_rule(rules.Simplify())
        s = calc.parse_expr("D t. I(t)")
        calc.perform_rule(rules.ApplyEquation(Eq1.goal,s))
        calc.perform_rule(rules.Simplify())
        assert Eq2.is_finished()
        Eq3 = file.add_goal("1/2 * t ^ 2 + log(I(t)) = SKOLEM_CONST(C)", conds=["I(t) > 0"])
        Eq3_proof = Eq3.proof_by_rewrite_goal(begin = Eq2)
        calc = Eq3_proof.begin
        calc.perform_rule(rules.IntegralEquation())
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        assert Eq3.is_finished()
        Eq4 = file.add_goal("log(sqrt(pi / 2)) = SKOLEM_CONST(C)")
        Eq4_proof = Eq4.proof_by_rewrite_goal(begin = Eq3)
        calc = Eq4_proof.begin
        calc.perform_rule(rules.LimitEquation('t', expr.Const(0)))
        calc.perform_rule(rules.Simplify())
        s = calc.parse_expr("I(0)")
        calc.perform_rule(rules.ApplyEquation(Eq0.goal,s))
        assert Eq4.is_finished()
        Eq5 = file.add_goal("log(I(t)) = -t ^ 2 / 2 + log(sqrt(pi / 2))", conds=["I(t) > 0"])
        Eq5_proof = Eq5.proof_by_calculation()
        calc = Eq5_proof.lhs_calc
        s = calc.parse_expr("log(I(t))")
        calc.perform_rule(rules.ApplyEquation(Eq3.goal, s))
        s = calc.parse_expr("SKOLEM_CONST(C)")
        calc.perform_rule(rules.ApplyEquation(Eq4.goal, s))
        calc.perform_rule(rules.Simplify())
        calc = Eq5_proof.rhs_calc
        calc.perform_rule(rules.Simplify())
        assert Eq5.is_finished()
        Eq6 = file.add_goal("I(t) = sqrt(pi/2) * exp(-t^2/2)", conds=["I(t) > 0"])
        Eq6_proof = Eq6.proof_by_rewrite_goal(begin = Eq5)
        calc = Eq6_proof.begin
        calc.perform_rule(rules.SolveEquation(parser.parse_expr("I(t)")))
        calc.perform_rule(rules.Equation(
            "exp(-(t ^ 2 / 2) - log(2) / 2 + log(pi) / 2)",
            "2 ^ (1/2) ^ (-1) * pi ^ (1/2) / exp(1/2 * t ^ 2)"))
        calc.perform_rule(rules.Simplify())
        assert Eq6.is_finished()
        Eq7 = file.add_goal("(INT x:[-oo, oo]. exp(-x^2/2) * cos(s + t*x)) = sqrt(2*pi)*exp(-t^2/2)*cos(s)", conds=['I(t)>0'])
        Eq7_proof = Eq7.proof_by_calculation()
        calc = Eq7_proof.lhs_calc
        calc.perform_rule(rules.ApplyIdentity("cos(s + t*x)", "cos(s)*cos(t*x)-sin(s)*sin(t*x)"))
        calc.perform_rule(rules.Equation("exp(-(x ^ 2) / 2) * (cos(s) * cos(t * x) - sin(s) * sin(t * x))",
                                         "exp(-(x ^ 2) / 2) * cos(s) * cos(t * x) - exp(-(x^2)/2)*sin(s) * sin(t * x)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.OnLocation(rules.SplitRegion("0"), "1.1"))
        calc.perform_rule(rules.OnLocation(rules.Substitution(var_name="x", var_subst="-x"), "1.1.0"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.SplitRegion("0"))
        calc.perform_rule(rules.OnLocation(rules.Substitution(var_name="x", var_subst="-x"), "1.0"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("-(x^2/2)", "-x^2/2"))
        calc.perform_rule(rules.OnSubterm(rules.FoldDefinition("I")))
        s = calc.parse_expr("I(t)")
        calc.perform_rule(rules.ApplyEquation(Eq6.goal, s))
        calc.perform_rule(rules.Equation("2 * cos(s) * (sqrt(pi / 2) * exp(-(t ^ 2) / 2))",
                                         "sqrt(2*pi)*exp(-t^2/2)*cos(s)"))
        assert Eq7.is_finished()
        # self.checkAndOutput(file)

    def testLeibniz03New(self):
        # Reference:
        # Inside interesting integrals, Section 3.1, example #3

        # Overall goal: INT x:[0,oo]. cos(tx)*exp(-(x^2)/2) = sqrt(pi/2)*exp(-(t^2)/2)

        # Initial state
        file = compstate.CompFile("interesting", 'leibniz03_new')
        # Make definition
        file.add_definition("I(t) = INT x:[0,oo]. cos(t*x)*exp(-(x^2)/2)")
        goal = file.add_goal("Gamma(n+1/2) = sqrt(pi) * factorial(2*n) / (4^n * factorial(n))", conds=["n>=0", 'isInt(n)'])
        proof = goal.proof_by_induction(induct_var='n')
        base_proof = proof.base_case.proof_by_calculation()
        induct_proof = proof.induct_case.proof_by_calculation()
        calc = base_proof.lhs_calc
        # calc.perform_rule(rules.Simplify())
        s = calc.parse_expr("Gamma(1/2)")
        calc.perform_rule(rules.ApplyEquation("Gamma(1/2) = sqrt(pi)", s))
        calc = induct_proof.lhs_calc
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("n+3/2")
        s2 = calc.parse_expr("(n+1/2)+1")
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = parser.parse_expr("Gamma(n + 1/2 + 1)")
        s2 = parser.parse_expr("(n+1/2) * Gamma(n+1/2)")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        calc.perform_rule(rules.OnLocation(rules.ApplyInductHyp(), '1'))
        s1 = parser.parse_expr("n+1/2")
        s2 = parser.parse_expr("(2*n+1)*(2*n+2) / (4 * (n+1))")
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = parser.parse_expr("(2 * n + 1) * (2 * n + 2) / (4 * (n + 1)) * (sqrt(pi) * factorial(2 * n) / (4 ^ n * factorial(n)))")
        s2 = parser.parse_expr("sqrt(pi) * ((2*n+1+1) * ((2*n+1) * factorial(2*n)))/ (4^1*4^n*((n+1)*factorial(n)))")
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = parser.parse_expr("4^1 * 4^n")
        s2 = parser.parse_expr("4^(1+n)")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        s1 = parser.parse_expr("(2*n+1)*factorial(2*n)")
        s2 = parser.parse_expr("factorial(2*n+1)")
        calc.perform_rule(rules.ApplyIdentity(s1,s2))
        s1 = parser.parse_expr("(2*n+1+1)*factorial(2*n+1)")
        s2 = parser.parse_expr("factorial(2*n+2)")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        s1 = parser.parse_expr("(n+1)*factorial(n)")
        s2 = parser.parse_expr("factorial(n+1)")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        calc.perform_rule(rules.Simplify())
        self.assertTrue(goal.is_finished())
        goal1 = file.add_goal("converges(SUM(n, 0, oo, INT x:[0,oo]. (-(-1)) ^ n * (abs(t) * abs(x)) ^ (2 * n) / factorial(2 * n) * exp(-(x ^ 2 / 2))))")
        cond = parser.parse_expr("t")
        proof = goal1.proof_by_case(cond)
        proofa = proof.cases[2].proof_by_calculation()
        calc = proofa.arg_calc
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("(t*x)^(2*n)")
        s2 = calc.parse_expr("t^(2*n) * x^(2*n)")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.OnLocation(rules.Substitution('u', 'x^2/2'), '0.1'))
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("(sqrt(u) * sqrt(2)) ^ (2 * n)")
        s2 = calc.parse_expr("((sqrt(u)*sqrt(2))^2)^n")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("(2*u)^n")
        s2 = calc.parse_expr("2^n * u^n")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("u ^ (n - 1/2) * exp(-u)")
        s2 = calc.parse_expr("exp(-u) * u^(n+1/2 - 1)")
        calc.perform_rule(rules.Equation(s1, s2))
        calc.perform_rule(rules.OnLocation(rules.FoldDefinition("Gamma"), "1.0.1"))
        s = calc.parse_expr("Gamma(n+1/2)")
        calc.perform_rule(rules.ApplyEquation(goal.goal, s))
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("t^(2*n)")
        s2 = calc.parse_expr("(t^2)^n")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        s1 = calc.parse_expr("(t ^ 2 ) ^ n * 2 ^ n")
        s2 = calc.parse_expr("(t^2 * 2) ^ n")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        s1 = calc.parse_expr("4^-n")
        s2 = calc.parse_expr("(1/4)^n")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        s1 = calc.parse_expr("(t ^ 2 * 2) ^ n * (1/4) ^ n")
        s2 = calc.parse_expr("(t^2/2)^n")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        calc.perform_rule(rules.OnLocation(rules.SeriesEvaluationIdentity(), "1"))
        proofb = proof.cases[0].proof_by_calculation()
        calc = proofb.arg_calc
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("(t*x)^(2*n)")
        s2 = calc.parse_expr("t^(2*n) * x^(2*n)")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.OnLocation(rules.Substitution('u', 'x^2/2'), '0.1'))
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("(sqrt(u) * sqrt(2)) ^ (2 * n)",)
        s2 = calc.parse_expr("((sqrt(u)*sqrt(2))^2)^n")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("(2*u)^n")
        s2 = calc.parse_expr("2^n * u^n")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("u ^ (n - 1/2) * exp(-u)")
        s2 = calc.parse_expr("exp(-u) * u^(n+1/2 - 1)")
        calc.perform_rule(rules.Equation(s1, s2))
        calc.perform_rule(rules.OnLocation(rules.FoldDefinition("Gamma"), "1.0.1"))
        # eq = calc.parse_expr("Gamma(n+1/2) = sqrt(pi) * factorial(2*n) / (4^n * factorial(n))")
        s = calc.parse_expr("Gamma(n+1/2)")
        calc.perform_rule(rules.ApplyEquation(goal.goal, s))
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("t^(2*n)")
        s2 = calc.parse_expr("(t^2)^n")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        s1 = calc.parse_expr("(t ^ 2 ) ^ n * 2 ^ n")
        s2 = calc.parse_expr("(t^2 * 2) ^ n")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        s1 = calc.parse_expr("4^-n")
        s2 = calc.parse_expr("(1/4)^n")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        s1 = calc.parse_expr("(t ^ 2 * 2) ^ n * (1/4) ^ n")
        s2 = calc.parse_expr("(t^2/2)^n")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        calc.perform_rule(rules.OnLocation(rules.SeriesEvaluationIdentity(), "1"))
        proofc = proof.cases[1].proof_by_calculation()
        calc = proofc.arg_calc
        calc.perform_rule(rules.Simplify())
        self.assertTrue(goal1.is_finished())

        goal2 = file.add_goal("I(t) = sqrt(pi/2) * exp(-t^2/2)", conds=['t<0'])
        proof = goal2.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.ExpandDefinition("I"))
        calc.perform_rule(rules.OnLocation(rules.SeriesExpansionIdentity(), "0.0"))
        s1 = parser.parse_expr("SUM(n, 0, oo, (-1) ^ n * (t * x) ^ (2 * n) / factorial(2 * n)) * exp(-(x ^ 2 / 2))")
        s2 = parser.parse_expr("SUM(n, 0, oo, (-1) ^ n * (t * x) ^ (2 * n) / factorial(2 * n) * exp(-(x^2/2)))")
        calc.perform_rule(rules.Equation(s1,s2))
        calc.perform_rule(rules.IntSumExchange())
        s1 = calc.parse_expr("(t*x)^(2*n)")
        s2 = calc.parse_expr("t^(2*n) * x^(2*n)")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.OnLocation(rules.Substitution('u', 'x^2/2'), '0.1'))
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("(sqrt(u) * sqrt(2)) ^ (2 * n)")
        s2 = calc.parse_expr("((sqrt(u)*sqrt(2))^2)^n")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("(2*u)^n")
        s2 = calc.parse_expr("2^n * u^n")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("u ^ (n - 1/2) * exp(-u)")
        s2 = calc.parse_expr("exp(-u) * u^(n+1/2 - 1)")
        calc.perform_rule(rules.Equation(s1, s2))
        calc.perform_rule(rules.OnLocation(rules.FoldDefinition("Gamma"), "1.0.1"))
        eq = calc.parse_expr("Gamma(n+1/2) = sqrt(pi) * factorial(2*n) / (4^n * factorial(n))")
        s = calc.parse_expr("Gamma(n+1/2)")
        calc.perform_rule(rules.ApplyEquation(goal.goal, s))
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("t^(2*n)")
        s2 = calc.parse_expr("(t^2)^n")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        s1 = calc.parse_expr("(t^2)^n * (-1)^n")
        s2 = calc.parse_expr("(t^2*-1)^n")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        s1 = calc.parse_expr("(t ^ 2 * -1) ^ n * 2 ^ n")
        s2 = calc.parse_expr("(t^2 * -2) ^ n")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        s1 = calc.parse_expr("4^-n")
        s2 = calc.parse_expr("(1/4)^n")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        s1 = calc.parse_expr("(t ^ 2 * -2) ^ n * (1/4) ^ n")
        s2 = calc.parse_expr("(-t^2/2)^n")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        calc.perform_rule(rules.OnLocation(rules.SeriesEvaluationIdentity(), '1'))
        calc.perform_rule(rules.Simplify())
        calc = proof.rhs_calc
        calc.perform_rule(rules.Simplify())
        self.assertTrue(goal2.is_finished())
        goal3 = file.add_goal("I(t) = sqrt(pi/2) * exp(-t^2/2)", conds=['t>0'])
        proof = goal3.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.ExpandDefinition("I"))
        calc.perform_rule(rules.OnLocation(rules.SeriesExpansionIdentity(), "0.0"))
        s1 = parser.parse_expr("SUM(n, 0, oo, (-1) ^ n * (t * x) ^ (2 * n) / factorial(2 * n)) * exp(-(x ^ 2 / 2))")
        s2 = parser.parse_expr("SUM(n, 0, oo, (-1) ^ n * (t * x) ^ (2 * n) / factorial(2 * n) * exp(-(x^2/2)))")
        calc.perform_rule(rules.Equation(s1, s2))
        calc.perform_rule(rules.IntSumExchange())
        s1 = calc.parse_expr("(t*x)^(2*n)")
        s2 = calc.parse_expr("t^(2*n) * x^(2*n)")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.OnLocation(rules.Substitution('u', 'x^2/2'), '0.1'))
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("(sqrt(u) * sqrt(2)) ^ (2 * n)")
        s2 = calc.parse_expr("((sqrt(u)*sqrt(2))^2)^n")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("(2*u)^n")
        s2 = calc.parse_expr("2^n * u^n")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("u ^ (n - 1/2) * exp(-u)")
        s2 = calc.parse_expr("exp(-u) * u^(n+1/2 - 1)")
        calc.perform_rule(rules.Equation(s1, s2))
        calc.perform_rule(rules.OnLocation(rules.FoldDefinition("Gamma"), "1.0.1"))
        eq = calc.parse_expr("Gamma(n+1/2) = sqrt(pi) * factorial(2*n) / (4^n * factorial(n))")
        s = calc.parse_expr("Gamma(n+1/2)")
        calc.perform_rule(rules.ApplyEquation(goal.goal, s))
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("t^(2*n)")
        s2 = calc.parse_expr("(t^2)^n")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        s1 = calc.parse_expr("(t^2)^n * (-1)^n")
        s2 = calc.parse_expr("(t^2*-1)^n")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        s1 = calc.parse_expr("(t ^ 2 * -1) ^ n * 2 ^ n")
        s2 = calc.parse_expr("(t^2 * -2) ^ n")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        s1 = calc.parse_expr("4^-n")
        s2 = calc.parse_expr("(1/4)^n")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        s1 = calc.parse_expr("(t ^ 2 * -2) ^ n * (1/4) ^ n")
        s2 = calc.parse_expr("(-t^2/2)^n")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        calc.perform_rule(rules.OnLocation(rules.SeriesEvaluationIdentity(), '1'))
        calc.perform_rule(rules.Simplify())
        calc = proof.rhs_calc
        calc.perform_rule(rules.Simplify())
        self.assertTrue(goal3.is_finished())
        # self.checkAndOutput(file)

    def testGaussianPowerExp(self):
        # Reference:
        # Inside interesting integrals, Section 2.3
        file = compstate.CompFile("interesting", 'gaussianPowerExp')
        file.add_definition("I(n) = (INT x:[0, oo]. x^(2*n) * exp(-x^2))", conds=["n>=0", "isInt(n)"])

        goal01 = file.add_goal("(INT x:[0, oo]. (D x. x^(2*n-1)*exp(-x^2))) = 0", conds=["n>=1", "isInt(n)"])
        proof = goal01.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Simplify())
        self.assertTrue(goal01.is_finished())
        goal02 = file.add_goal("(INT x:[0, oo]. (D x. x^(2*n-1)*exp(-x^2))) = (2*n-1)*I(n-1) - 2 * I(n)",  conds=["n>=1", "isInt(n)"])
        proof = goal02.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.OnSubterm(rules.ExpandDefinition("I")))
        calc.perform_rule(rules.OnLocation(rules.Simplify(), "0"))
        calc.perform_rule(rules.Simplify())
        e1 = parser.parse_expr("2*n-2")
        e2 = parser.parse_expr("2*(n-1)")
        calc.perform_rule(rules.Equation(e1, e2))
        calc.perform_rule(rules.OnSubterm(rules.FoldDefinition("I")))
        self.assertTrue(goal02.is_finished())
        goal03 = file.add_goal("I(n) = (2 * n - 1) / 2 * I(n - 1)",  conds=["n>=1", "isInt(n)"])
        proof = goal03.proof_by_rewrite_goal(begin = goal01)
        calc = proof.begin
        s = calc.parse_expr("INT x:[0,oo]. D x. x ^ (2 * n - 1) * exp(-(x ^ 2))")
        calc.perform_rule(rules.ApplyEquation(goal02.goal,s))
        x = calc.parse_expr("I(n)")
        calc.perform_rule(rules.SolveEquation(x))
        e1 = calc.parse_expr("I(n - 1) * (2 * n - 1) / 2")
        e2 = calc.parse_expr("(2 * n - 1) / 2 * I(n - 1)")
        calc.perform_rule(rules.Equation(e1, e2))
        self.assertTrue(goal03.is_finished())
        goal04 = file.add_goal("I(n) = factorial(2*n)/(4^n*factorial(n))*(1/2)*sqrt(pi)", conds=["n>=0", "isInt(n)"])
        proof = goal04.proof_by_induction("n", 0)
        proof_base = proof.base_case.proof_by_calculation()
        proof_induct = proof.induct_case.proof_by_calculation()
        calc = proof_base.lhs_calc
        calc.perform_rule(rules.ExpandDefinition("I"))
        calc.perform_rule(rules.Substitution(var_name="x", var_subst="sqrt(2)*x"))
        e1 = parser.parse_expr("-(x ^ 2 / 2)")
        e2 = parser.parse_expr("-(x ^ 2) / 2")
        calc.perform_rule(rules.Equation(e1, e2))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        calc = proof_induct.lhs_calc
        calc.perform_rule(rules.ApplyEquation(goal03.goal, calc.parse_expr("I(n+1)")))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.OnSubterm(rules.ApplyInductHyp()))
        e1 = calc.parse_expr("factorial(2 * n) / (4 ^ n * factorial(n)) * (1/2) * sqrt(pi) * (2 * n + 1) / 2")
        e2 = calc.parse_expr("4 ^ -n * sqrt(pi) * (((2*n+1)*factorial(2*n))/(4 * factorial(n)))")
        calc.perform_rule(rules.Equation(e1, e2))
        e1 = parser.parse_expr("(2*n+1)*factorial(2*n)")
        e2 = parser.parse_expr("factorial(2*n+1)")
        calc.perform_rule(rules.ApplyIdentity(e1, e2))
        e1 = parser.parse_expr("factorial(2 * n + 1) / (4 * factorial(n))")
        e2 = parser.parse_expr("((2*n+1 + 1)*factorial(2 * n + 1)) / (8 * ((n+1)*factorial(n)))")
        calc.perform_rule(rules.Equation(e1, e2))
        e1 = parser.parse_expr("(2*n+1 + 1)*factorial(2*n+1)")
        e2 = parser.parse_expr("factorial(2*n+2)")
        calc.perform_rule(rules.ApplyIdentity(e1, e2))
        e1 = parser.parse_expr("(n+1)*factorial(n)")
        e2 = parser.parse_expr("factorial(n+1)")
        calc.perform_rule(rules.ApplyIdentity(e1, e2))
        e1 = parser.parse_expr("4 ^ -n * sqrt(pi) * (factorial(2 * n + 2) / (8 * factorial(n + 1)))")
        e2 = parser.parse_expr("4 ^ -n * sqrt(pi) * factorial(2 * n + 2) / (8 * factorial(n + 1))")
        calc.perform_rule(rules.Equation(e1, e2))

        self.assertTrue(goal04.is_finished())
        # self.checkAndOutput(file)

    def testEulerLogSineIntegral(self):
        # Reference:
        # Inside interesting integrals, Section 2.4
        file = compstate.CompFile("interesting", "euler_log_sin")

        goal = file.add_goal("(INT x:[0,pi/2]. log(a * sin(x))) = pi/2 * log(a/2)", conds=["a > 0"])

        # Define I(a)
        goal.add_definition("I(a) = INT x:[0,pi/2]. log(a * sin(x))", conds=["a > 0"])

        # Define J(a)
        goal.add_definition("J(a) = INT x:[0,pi/2]. log(a * sin(2*x))", conds=["a > 0"])

        # Prove J(a) = I(a)
        goal1 = goal.add_subgoal("1", "J(a) = I(a)", conds=['a > 0'])
        proof = goal1.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.ExpandDefinition("J"))
        calc.perform_rule(rules.Substitution("t", parser.parse_expr("2*x")))
        calc.perform_rule(rules.SplitRegion(parser.parse_expr('pi/2')))
        calc.perform_rule(rules.OnCount(rules.Substitution('x', parser.parse_expr('pi - t')), 2))
        calc.perform_rule(rules.Simplify())

        calc = proof.rhs_calc
        calc.perform_rule(rules.ExpandDefinition("I"))

        # Prove J(a) = pi/2 * log(2/a) + 2 * I(a)
        goal2 = goal.add_subgoal("2", "J(a) = pi/2 * log(2/a) + 2 * I(a)", conds=["a > 0"])
        proof = goal2.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.ExpandDefinition("J"))
        calc.perform_rule(rules.ApplyIdentity("sin(2*x)", "2 * sin(x) * cos(x)"))
        calc.perform_rule(rules.Equation("a * (2 * sin(x) * cos(x))",
                                         "(2/a) * (a*sin(x)) * (a*cos(x))"))
        calc.perform_rule(rules.ApplyIdentity(
            "log(2/a * (a*sin(x)) * (a*cos(x)))",
            "log(2/a * (a*sin(x))) + log(a*cos(x))"))
        calc.perform_rule(rules.ApplyIdentity(
            "log(2/a * (a*sin(x)))",
            "log(2/a) + log(a*sin(x))"))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        calc.perform_rule(rules.Substitution('t', parser.parse_expr("pi/2 - x")))
        calc.perform_rule(rules.Substitution('x', parser.parse_expr("t")))
        calc.perform_rule(rules.Simplify())

        calc = proof.rhs_calc
        calc.perform_rule(rules.OnSubterm(rules.ExpandDefinition("I")))
        calc.perform_rule(rules.Simplify())

        # Finally show I(a) = pi/2 * log(a/2)
        proof = goal.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.OnSubterm(rules.FoldDefinition("I")))
        s = calc.parse_expr("I(a)")
        calc.perform_rule(rules.ApplyEquation("1", s))
        s = calc.parse_expr("J(a)")
        calc.perform_rule(rules.ApplyEquation("2", s))
        calc.perform_rule(rules.IntegrateByEquation(parser.parse_expr("I(a)")))
        calc.perform_rule(rules.ApplyIdentity("log(2 / a)", "log(2) + log(1 / a)"))
        calc.perform_rule(rules.ExpandPolynomial())
        calc.perform_rule(rules.Simplify())
        calc = proof.rhs_calc
        calc.perform_rule(rules.ApplyIdentity("log(a / 2)", "log(a) - log(2)"))
        calc.perform_rule(rules.ExpandPolynomial())
        assert goal.is_finished()

    def testEulerLogSineIntegral02(self):
        # Reference:
        # Inside interesting integrals, Section 2.4 (2.4.2)
        file = compstate.CompFile("interesting", "euler_log_sin02")
        goal = file.add_goal("(INT x:[0, pi/2]. log(sin(x) / x)) = pi/2 * (1-log(pi))")
        proof = goal.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.ApplyIdentity("log(sin(x) / x)", "log(sin(x)) - log(x)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("log(sin(x))", "log(1*sin(x))"))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.IntegrationByParts("log(x)", "x"))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.ExpandPolynomial())
        calc.perform_rule(rules.Simplify())

        calc = proof.rhs_calc
        calc.perform_rule(rules.ExpandPolynomial())
        # self.checkAndOutput(file)

    def testEulerLogSineIntegral0304(self):
        # Reference:
        # Inside interesting integrals, Section 2.4 (2.4.3 & 2.4.4)

        file = compstate.CompFile("interesting", "euler_log_sin0304")

        goal = file.add_goal("(INT x:[0, 1]. log(x+1/x)/(x^2+1)) = pi/2 * log(2)")

        goal01 = goal.add_subgoal("1", "(INT x:[0, oo]. log(x^2+1)/(x^2+1)) = pi * log(2)")
        proof = goal01.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.SubstitutionInverse("u", "x", "tan(u)"))
        calc.perform_rule(rules.ApplyIdentity("sec(u)^2", "tan(u)^2 + 1"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.ApplyIdentity("tan(u)^2+1", "sec(u)^2"))
        calc.perform_rule(rules.ApplyIdentity("sec(u)", "cos(u)^(-1)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Substitution(var_name="x", var_subst="pi/2 - u"))
        calc.perform_rule(rules.Equation("sin(x)", "1 * sin(x)"))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        proof = goal.proof_by_rewrite_goal(begin="1")
        calc = proof.begin
        calc.perform_rule(rules.SplitRegion("1"))

        # introduce new variable to void infinite recursion when checking conditions
        calc.perform_rule(rules.OnCount(rules.Substitution(var_name="y", var_subst="1/x"), 2))
        calc.perform_rule(rules.Equation("y ^ 2 * (1 / y ^ 2 + 1)", "y^2+1"))
        calc.perform_rule(rules.Equation("1 / (y ^ 2 + 1) * log(1 / y ^ 2 + 1)",
                                         "log(1/y^2 + 1)/(y^2+1)"))
        calc.perform_rule(rules.Equation("(INT y:[0,1]. log(y ^ 2 + 1) / (y ^ 2 + 1)) + (INT y:[0,1]. log(1 / y ^ 2 + 1) / (y ^ 2 + 1))",
                                         "(INT y:[0,1]. log(y ^ 2 + 1) / (y ^ 2 + 1) + log(1 / y ^ 2 + 1) / (y ^ 2 + 1))"))
        calc.perform_rule(rules.Equation("log(y ^ 2 + 1) / (y ^ 2 + 1) + log(1 / y ^ 2 + 1) / (y ^ 2 + 1)",
                                         "(log(y ^ 2 + 1) + log(1 / y ^ 2 + 1)) / (y ^ 2 + 1)"))
        calc.perform_rule(rules.Equation("log(y ^ 2 + 1) + log(1 / y ^ 2 + 1)", "log((y^2+1)*(1/y^2+1))"))
        calc.perform_rule(rules.Equation("(y ^ 2 + 1) * (1 / y ^ 2 + 1)", "(y + 1/y)^2"))
        calc.perform_rule(rules.ApplyIdentity("log((y + 1 / y) ^ 2)", "2*log(y+1/y)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("1 / (y ^ 2 + 1) * log(1 / y + y)", "log(y+1/y)/(y^2+1)"))
        calc.perform_rule(rules.SolveEquation("(INT y:[0, 1]. log(y+1/y)/(y^2+1))"))
        assert goal.is_finished()

    def testEulerLogSineIntegral05(self):
        # Reference:
        # Inside interesting integrals, Section 2.4 (2.4.5)
        file = compstate.CompFile("interesting", "euler_log_sin05")

        goal = file.add_goal("(INT x:[0,oo]. log(x) / (x^2-b*x+1)) = 0", conds=["b > -2", "b < 2"])

        sub_goal1 = goal.add_subgoal("1", "x ^ 2 - b * x + 1 != 0")
        proof = sub_goal1.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Equation("x ^ 2 - b * x + 1", "(x - 1/2 * b) ^ 2 + 1 - 1/4 * b^2"))
        assert sub_goal1.is_finished()

        sub_goal2 = goal.add_subgoal("2", "(INT x:[0, oo]. log(x^a+1) / (x^2 - b*x + 1)) = \
        (INT x:[0, oo]. log(x^a+1) / (x^2 - b*x + 1)) - a * INT x:[0,oo]. log(x) / (x^2-b*x+1)", conds=["a>0"])
        proof = sub_goal2.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.SubstitutionInverse("u", "x", "1/u"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.ExpandPolynomial())
        calc.perform_rule(rules.ApplyIdentity("(1/u)^a", "1^a / u^a"))
        calc.perform_rule(rules.Equation("1 ^ a / u ^ a + 1", "(1+u^a) / u^a"))
        calc.perform_rule(rules.Equation("log((1 + u ^ a) / u ^ a)", "log(1+u^a) - log(u^a)"))
        calc.perform_rule(rules.ExpandPolynomial())
        calc.perform_rule(rules.Simplify())
        assert sub_goal2.is_finished()
        proof = goal.proof_by_rewrite_goal(begin="2")
        calc = proof.begin
        calc.perform_rule(rules.SolveEquation("INT x:[0,oo]. log(x) / (x^2-b*x+1)"))
        assert goal.is_finished()

    def testEulerLogSineIntegral06(self):
        # Reference:
        # Inside interesting integrals, Section 2.4 (2.4.6)
        file = compstate.CompFile("interesting", "euler_log_sin06")
        goal = file.add_goal("(INT x:[0,1]. (1-x) / (1+x+x^2)) = sqrt(3) * pi / 6 - log(3) / 2")
        proof = goal.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Equation("(1+x+x^2)", "(x+1/2)^2 + 3/4"))
        calc.perform_rule(rules.Substitution("u", "2*(x+1/2)/sqrt(3)"))
        calc.perform_rule(rules.Equation("3 * u ^ 2 / 2 + 3/2", "3/2*(u^2+1)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("1 / (u ^ 2 + 1) * (-(u * sqrt(3) / 2) + 3/2)", "-sqrt(3) / 2 * (u / (u ^ 2 + 1)) + 3 / 2 * (1 / (u ^ 2 + 1))"))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Substitution("t", "u^2 + 1"))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.ExpandPolynomial())
        calc.perform_rule(rules.Simplify())
        # self.checkAndOutput(file)

    def testDirichletIntegral(self):
        # Reference:
        # Inside interesting integrals, Section 3.2
        file = compstate.CompFile("interesting", "dirichletIntegral")

        # Define g(y)
        file.add_definition("g(y, a) = INT x:[0,oo]. exp(-x * y) * sin(a * x) / x", conds=["y > 0"])

        # Differentiate g(y)
        goal2 = file.add_goal("(D y. g(y, a)) = - a / (a ^ 2 + y ^ 2)", conds=["y > 0", "a!=0"])
        proof = goal2.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.OnSubterm(rules.ExpandDefinition("g")))
        calc.perform_rule(rules.DerivIntExchange())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc = proof.rhs_calc
        calc.perform_rule(rules.Simplify())
        assert goal2.is_finished()

        # Integrate the previous equation on both sides
        goal3 = file.add_goal("g(y, a) = -atan(y / a) + SKOLEM_FUNC(C(a))", conds=["y > 0", "a > 0"])
        proof = goal3.proof_by_rewrite_goal(begin = goal2)
        calc = proof.begin
        calc.perform_rule(rules.IntegralEquation())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        assert goal3.is_finished()

        goal3_1 = file.add_goal("g(y, a) = -atan(y / a) + SKOLEM_FUNC(C(a))", conds=["y > 0", "a < 0"])
        proof = goal3_1.proof_by_rewrite_goal(begin=goal2)
        calc = proof.begin
        calc.perform_rule(rules.IntegralEquation())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        assert goal3_1.is_finished()

        # Evaluate the case y = oo
        goal4 = file.add_goal("(LIM {y -> oo}. g(y, a)) = 0", conds=["y > 0"])
        proof = goal4.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.OnSubterm(rules.ExpandDefinition("g")))
        calc.perform_rule(rules.Simplify())
        assert goal4.is_finished()
        # Evaluate C(a) for a > 0
        goal5 = file.add_goal("SKOLEM_FUNC(C(a)) = pi / 2", conds=["a > 0"])
        proof = goal5.proof_by_rewrite_goal(begin = goal3)
        calc = proof.begin
        calc.perform_rule(rules.LimitEquation("y", parser.parse_expr("oo")))
        s = calc.parse_expr("LIM {y -> oo}. g(y,a)")
        calc.perform_rule(rules.ApplyEquation(goal4.goal, s))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.SolveEquation(parser.parse_expr("SKOLEM_FUNC(C(a))")))
        assert goal5.is_finished()
        # Evaluate C(a) for a < 0
        goal6 = file.add_goal("SKOLEM_FUNC(C(a)) = -(pi / 2)", conds=["a < 0"])
        proof = goal6.proof_by_rewrite_goal(begin = goal3_1)
        calc = proof.begin
        calc.perform_rule(rules.LimitEquation("y", parser.parse_expr("oo")))
        s = calc.parse_expr("LIM {y -> oo}. g(y,a)")
        calc.perform_rule(rules.ApplyEquation(goal4.goal, s))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.SolveEquation(parser.parse_expr("SKOLEM_FUNC(C(a))")))
        assert goal6.is_finished()

        goal11 = file.add_goal("g(y,a) = pi / 2 -atan(y / a)", conds=["a>0", "y>0"])
        proof = goal11.proof_by_rewrite_goal(begin = goal3)
        calc = proof.begin
        s = calc.parse_expr("SKOLEM_FUNC(C(a))")
        calc.perform_rule(rules.ApplyEquation(goal5.goal, s))
        assert goal11.is_finished()

        goal12 = file.add_goal("g(y,a) = -pi / 2 -atan(y / a)", conds=["a<0", "y>0"])
        proof = goal12.proof_by_rewrite_goal(begin=goal3_1)
        calc = proof.begin
        s = calc.parse_expr("SKOLEM_FUNC(C(a))")
        calc.perform_rule(rules.ApplyEquation(goal6.goal, s))
        assert goal12.is_finished()
        # TODO: This proof is not valid
        # Case y = 0: g(0) = INT x:[0, oo]. sin(a * x) / x
        # goal7 = file.add_goal("g(0, a) = INT x:[0,oo]. sin(a * x) / x")
        # proof = goal7.proof_by_calculation()
        # calc = proof.lhs_calc
        # calc.perform_rule(rules.ExpandDefinition("g"))  # check 0 satisfies condition y >= 0
        # assert goal7.is_finished()
        # Final result: case a > 0
        # goal8 = file.add_goal("(INT x:[0,oo]. sin(a * x) / x) = pi / 2", conds=["a > 0"])
        # proof = goal8.proof_by_calculation()
        # calc = proof.lhs_calc
        # calc.perform_rule(rules.ApplyEquation(goal7.goal))
        # calc.perform_rule(rules.ApplyEquation(goal3.goal))
        # calc.perform_rule(rules.OnLocation(rules.ApplyEquation(goal5.goal), "1"))
        # calc.perform_rule(rules.Simplify())
        # print(goal8)
        # assert goal8.is_finished()
        # # Final result: case a < 0
        # goal9 = file.add_goal("(INT x:[0,oo]. sin(a * x) / x) = -(pi / 2)", conds=["a < 0"])
        # proof = goal9.proof_by_calculation()
        # calc = proof.lhs_calc
        # calc.perform_rule(rules.ApplyEquation(goal7.goal))
        # calc.perform_rule(rules.ApplyEquation(goal3.goal))
        # calc.perform_rule(rules.OnLocation(rules.ApplyEquation(goal6.goal), "1"))
        # calc.perform_rule(rules.Simplify())
        # assert goal9.is_finished()
        # # Final result: case a = 0
        # goal10 = file.add_goal("(INT x:[0,oo]. sin(a * x) / x) = 0", conds=["a = 0"])
        # proof = goal10.proof_by_calculation()
        # calc = proof.lhs_calc
        # calc.perform_rule(rules.Simplify())
        # assert goal10.is_finished()
        # self.checkAndOutput(file)

    # TODO: Change the order of integration.
    # def testFlipSide02(self):
    #     # Reference:
    #     # Inside interesting integrals, Section 3.4 example #2
    #     file = compstate.CompFile("interesting", "filpside02")
    #
    #     goal01 = file.add_goal("(INT t:[0,oo]. (exp(-p*t^2)-exp(-q*t^2))/t^2) = (INT t:[0,oo]. (INT a:[p,q]. exp(-a*t^2)))", conds=["p>0","q>0"])
    #     proof_of_goal01 = goal01.proof_by_calculation()
    #     calc = proof_of_goal01.rhs_calc
    #     calc.perform_rule(rules.OnLocation(rules.Substitution(var_name="x", var_subst="-a*t"), "0"))
    #     calc.perform_rule(rules.OnLocation(rules.DefiniteIntegralIdentity(), "0"))
    #     calc.perform_rule(rules.Simplify())
    #     calc.perform_rule(rules.ReplaceSubstitution())
    #     calc.perform_rule(rules.Equation("1 / t * (-(exp(-(p * t ^ 2)) / t) + exp(-(q * t ^ 2)) / t)",
    #                                      "-exp(-p * t ^ 2) / t^2 + exp(-q * t ^ 2) / t^2"))
    #     calc.perform_rule(rules.Equation("-(INT t:[0,oo]. -exp(-p * t ^ 2) / t ^ 2 + exp(-q * t ^ 2) / t ^ 2)",
    #                                      "(INT t:[0,oo]. exp(-p * t ^ 2) / t ^ 2 - exp(-q * t ^ 2) / t ^ 2)"))
    #     calc.perform_rule(rules.Equation("exp(-p * t ^ 2) / t ^ 2 - exp(-q * t ^ 2) / t ^ 2",
    #                                      "(exp(-p*t^2)-exp(-q*t^2))/t^2"))
    #
    #     goal02 = file.add_goal("(INT t:[0,oo]. (INT a:[p,q]. exp(-a*t^2))) = sqrt(pi)*(sqrt(q) - sqrt(p))", conds=["p>0","q>0"])
    #     proof_of_goal02 = goal02.proof_by_calculation()
    #     calc = proof_of_goal02.lhs_calc
    #
    #     # self.checkAndOutput(file)
    def testFlipside03(self):
        # Reference:
        # Inside interesting integrals, Section 3.4, example #3

        # goal: INT x:[0, 1]. (x ^ a - 1) / log(x) = log(a + 1)
        file = compstate.CompFile("interesting", "flipside03")

        # introduce definition
        goal = file.add_goal("(INT x:[0, 1]. (x ^ a - 1) / log(x)) = log(a + 1)", conds=["a > -1"])

        goal.add_definition("I(a) = INT x:[0, 1]. (x ^ a - 1) / log(x)", conds=["a > -1"])

        # verify the following equation: D a. I(a) = 1/(a+1)
        goal1 = goal.add_subgoal("1", "(D a. I(a)) = 1/(a+1)", conds=["a > -1"])
        proof_of_goal1 = goal1.proof_by_calculation()
        calc = proof_of_goal1.lhs_calc
        calc.perform_rule(rules.OnSubterm(rules.ExpandDefinition("I")))
        calc.perform_rule(rules.DerivIntExchange())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        goal3 = goal.add_subgoal("2", "I(a) = log(a+1) + SKOLEM_CONST(C)", conds=["a > -1"])
        proof_of_goal3 = goal3.proof_by_rewrite_goal(begin="1")
        calc = proof_of_goal3.begin
        calc.perform_rule(rules.IntegralEquation())
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        assert goal3.is_finished()

        goal4 = goal.add_subgoal("3", "SKOLEM_CONST(C) = 0")
        proof_of_goal4 = goal4.proof_by_rewrite_goal(begin="2")
        calc = proof_of_goal4.begin
        calc.perform_rule(rules.VarSubsOfEquation([{'var': 'a', 'expr': "0"}]))
        calc.perform_rule(rules.OnSubterm(rules.ExpandDefinition("I")))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.SolveEquation(parser.parse_expr("SKOLEM_CONST(C)")))
        assert goal4.is_finished()

        proof_of_goal5 = goal.proof_by_calculation()
        calc = proof_of_goal5.lhs_calc
        calc.perform_rule(rules.OnSubterm(rules.FoldDefinition("I")))
        s = calc.parse_expr("I(a)")
        calc.perform_rule(rules.ApplyEquation("2", s))
        s = calc.parse_expr("SKOLEM_CONST(C)")
        calc.perform_rule(rules.ApplyEquation("3", s))
        calc.perform_rule(rules.Simplify())

        assert goal.is_finished()

    def testFlipside04(self):
        # Reference:
        # Inside interesting integrals, Section 3.4, example #4
        file = compstate.CompFile("interesting", "flipside04")

        goal = file.add_goal("(INT x:[0,1]. (x^a-x^b)/log(x)) = log((a+1)/(b+1))", conds=["a>-1", "b>-1"])
        proof = goal.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Equation("x^a-x^b", "(x^a-1)-(x^b-1)"))
        calc.perform_rule(rules.Equation("(x ^ a - 1 - (x ^ b - 1)) / log(x)", "(x ^ a - 1)/log(x) - (x ^ b - 1) / log(x)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.ApplyIdentity("log(a + 1) - log(b + 1)", "log((a+1)/(b+1))"))
        # self.checkAndOutput(file)

    def testFlipside07(self):
        # Reference:
        # Inside interesting integrals, Section 3.4, example #7
        file = compstate.CompFile("interesting", "flipside07")

        goal = file.add_goal("(INT x:[0,1]. x^a * (log(x))^2) = 2/(a+1)^3", conds=["a > -1"])

        goal01 = goal.add_subgoal("1", "(D a. (D a. (INT x:[0,1]. x^a))) = 2/(a+1)^3", conds=["a>-1"])
        proof = goal01.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        proof = goal.proof_by_rewrite_goal(begin="1")
        calc = proof.begin
        calc.perform_rule(rules.Simplify())
        assert goal.is_finished()

    # TODO: Solve LIM {z -> oo}. atan(z * sqrt(a - b) / sqrt(a + b))
    # def testFlipSide08(self):
    #     # Reference:
    #     # Inside interesting integrals, Section 3.4, example #7 and #8
    #     file = compstate.CompFile("interesting", "flipside08")
    #
    #     file.add_definition("I(a, b) = (INT x:[0,pi]. log(a+b*cos(x)))", conds=["a>b", "b>=0"])
    #
    #     goal01 = file.add_goal("(D a. I(a, b)) = (INT x:[0,pi]. 1/(a+b*cos(x)))", conds=["a>b", "b>=0"])
    #     proof = goal01.proof_by_calculation()
    #     calc = proof.lhs_calc
    #     calc.perform_rule(rules.OnSubterm(rules.ExpandDefinition("I")))
    #     calc.perform_rule(rules.Simplify())
    #     calc.perform_rule(rules.Equation("b*cos(x)+a", "a+b*cos(x)"))
    #
    #     goal02 = file.add_goal("(INT x:[0,pi]. 1/(a+b*cos(x))) = pi/(sqrt(a^2-b^2))", conds=["a>b", "b>=0"])
    #     proof = goal02.proof_by_calculation()
    #     calc = proof.lhs_calc
    #     calc.perform_rule(rules.Substitution(var_name="z", var_subst="tan(x/2)"))
    #     calc.perform_rule(rules.Simplify())
    #     calc.perform_rule(rules.Equation("1 / ((z ^ 2 + 1) * (b * (-(z ^ 2) + 1) / (z ^ 2 + 1) + a))",
    #                                      "1/(a*(1+z^2) + b*(1-z^2))"))
    #     calc.perform_rule(rules.Equation("1/(a*(1+z^2) + b*(1-z^2))", "(1/(a-b)) * 1/((a+b)/(a-b)+z^2)"))
    #     calc.perform_rule(rules.Simplify())
    #     calc.perform_rule(rules.Equation("(a+b)/(a-b)+z^2", "z^2+(a+b)/(a-b)"))
    #     calc.perform_rule(rules.DefiniteIntegralIdentity())
    #     calc.perform_rule(rules.Simplify())
    #
    #     # self.checkAndOutput(file)

    def testFrullaniIntegral01(self):
        # Reference:
        # Inside interesting integrals, Section 3.3
        file = compstate.CompFile("interesting", "FrullaniIntegral01")

        goal = file.add_goal("(INT x:[0,oo]. (atan(a*x) - atan(b*x)) / x) = pi * log(a) / 2 - pi * log(b) / 2",
                             conds=["a > 0", "b > 0"])

        # Define I(a, b)
        goal.add_definition("I(a, b) = INT x:[0,oo]. (atan(a*x) - atan(b*x)) / x", conds=["a > 0", "b > 0"])

        # Evalute D a. I(a, b) for a > 0
        goal2 = goal.add_subgoal("1", "(D a. I(a,b)) = pi / (2*a)", conds=["a > 0", "b > 0"])
        proof_of_goal2 = goal2.proof_by_calculation()
        calc = proof_of_goal2.lhs_calc
        calc.perform_rule(rules.OnSubterm(rules.ExpandDefinition("I")))
        calc.perform_rule(rules.DerivIntExchange())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Substitution('u' , parser.parse_expr("a*x")))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        # Integrate the previous result to get formula for I(a,b)
        goal3 = goal.add_subgoal("2", "I(a,b) = pi * log(a) / 2 + SKOLEM_FUNC(C(b))", conds=["a > 0", "b > 0"])
        proof_of_goal3 = goal3.proof_by_rewrite_goal(begin="1")
        calc = proof_of_goal3.begin
        calc.perform_rule(rules.IntegralEquation())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        # Obtain value of Skolem function
        goal4 = goal.add_subgoal("3", "SKOLEM_FUNC(C(a)) = -(pi * log(a) / 2)", conds=["a > 0"])
        proof_of_goal4 = goal4.proof_by_rewrite_goal(begin="2")
        calc = proof_of_goal4.begin
        calc.perform_rule(rules.VarSubsOfEquation([{'var': "b", 'expr': "a"}]))
        calc.perform_rule(rules.SolveEquation(parser.parse_expr("SKOLEM_FUNC(C(a))")))
        calc.perform_rule(rules.OnSubterm(rules.ExpandDefinition("I")))
        calc.perform_rule(rules.Simplify())

        # Final result
        proof_of_goal6 = goal.proof_by_calculation()
        calc = proof_of_goal6.lhs_calc
        calc.perform_rule(rules.OnSubterm(rules.FoldDefinition("I")))
        s = calc.parse_expr("I(a,b)")
        calc.perform_rule(rules.ApplyEquation("2", s))
        s = calc.parse_expr("SKOLEM_FUNC(C(b))")
        calc.perform_rule(rules.ApplyEquation("3", s))
        calc.perform_rule(rules.Simplify())
        goal.print_entry()
        assert goal.is_finished()

    # TODO: Calculating the limit lim{x->oo}. exp(-a*x) = 0
    # def testFrullaniIntegral02(self):
    #     # Reference:
    #     # Inside interesting integrals, Section 3.3
    #     file = compstate.CompFile("interesting", "FrullaniIntegral02")
    #
    #     file.add_definition("I(a) = (INT x:[0, oo]. (exp(-a*x)-1)/x)", conds=["a>0"])
    #
    #     goal01 = file.add_goal("(D a. I(a)) = -1/a")
    #     proof_of_goal01 = goal01.proof_by_calculation()
    #     calc = proof_of_goal01.lhs_calc
    #     calc.perform_rule(rules.OnSubterm(rules.ExpandDefinition("I")))
    #     calc.perform_rule(rules.Simplify())
    #     calc.perform_rule(rules.DefiniteIntegralIdentity())
    #     calc.perform_rule(rules.Simplify())
    #     # self.checkAndOutput(file)

    def testCombiningTwoTricks(self):
        # Reference:
        # Inside interesting integrals, Section 3.5
        file = compstate.CompFile("interesting", 'CombiningTwoTricks')

        file.add_definition("I(a,b,n) = (INT x:[0, pi/2]. 1/(a*cos(x)^2+b*sin(x)^2)^n)", conds=["a > 0", "b > 0", "n > 0", "isInt(n)"])
        file.add_definition("J(a,y,n) = (INT x:[0, y]. 1/(x^2+a^2)^n)", conds=["a>0", "y>=0", "n>0", "isInt(n)"])

        goal01 = file.add_goal("I(a, b, n) = -1/(n-1)*((D a. I(a,b, n-1))+(D b. I(a,b,n-1)))", conds=["a>0","b>0","n>=2","isInt(n)"])
        proof_of_goal01 = goal01.proof_by_calculation()
        calc = proof_of_goal01.rhs_calc
        calc.perform_rule(rules.OnSubterm(rules.ExpandDefinition("I")))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("(-n + 1) * (INT x:[0,pi / 2]. cos(x) ^ 2 * (a * cos(x) ^ 2 + b * sin(x) ^ 2) ^ -n) + (-n + 1) * (INT x:[0,pi / 2]. sin(x) ^ 2 * (a * cos(x) ^ 2 + b * sin(x) ^ 2) ^ -n)",
                                         "(-n+1)*((INT x:[0,pi / 2]. cos(x) ^ 2 * (a * cos(x) ^ 2 + b * sin(x) ^ 2) ^ -n) + (INT x:[0,pi / 2]. sin(x) ^ 2 * (a * cos(x) ^ 2 + b * sin(x) ^ 2) ^ -n))"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("(-n + 1) / (n - 1)", "-1"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("(INT x:[0,pi / 2]. cos(x) ^ 2 * (a * cos(x) ^ 2 + b * sin(x) ^ 2) ^ -n) + (INT x:[0,pi / 2]. sin(x) ^ 2 * (a * cos(x) ^ 2 + b * sin(x) ^ 2) ^ -n)",
                                         "(INT x:[0,pi / 2]. cos(x) ^ 2 * (a * cos(x) ^ 2 + b * sin(x) ^ 2) ^ -n + sin(x) ^ 2 * (a * cos(x) ^ 2 + b * sin(x) ^ 2) ^ -n)"))
        calc.perform_rule(rules.Equation("cos(x) ^ 2 * (a * cos(x) ^ 2 + b * sin(x) ^ 2) ^ -n + sin(x) ^ 2 * (a * cos(x) ^ 2 + b * sin(x) ^ 2) ^ -n",
                                         "(sin(x)^2+cos(x)^2)*(a * cos(x) ^ 2 + b * sin(x) ^ 2) ^ -n"))
        calc.perform_rule(rules.ApplyIdentity("sin(x)^2+cos(x)^2", "1"))
        calc.perform_rule(rules.Equation("1 * (a * cos(x) ^ 2 + b * sin(x) ^ 2) ^ -n",
                                         "1 / (a * cos(x) ^ 2 + b * sin(x) ^ 2) ^ n"))
        calc.perform_rule(rules.FoldDefinition("I"))

        goal02 = file.add_goal("I(a,b,1) = pi/(2*sqrt(a*b))", conds=["a>0", "b>0"])
        proof_of_goal02 = goal02.proof_by_calculation()
        calc = proof_of_goal02.lhs_calc
        calc.perform_rule(rules.ExpandDefinition("I"))
        calc.perform_rule(rules.Equation("1/(a*cos(x)^2+b*sin(x)^2)", "(1/cos(x))^2/(a+b*(sin(x)/cos(x))^2)"))
        calc.perform_rule(rules.ApplyIdentity("sin(x)/cos(x)", "tan(x)"))
        calc.perform_rule(rules.ApplyIdentity("1/cos(x)", "sec(x)"))
        calc.perform_rule(rules.Substitution(var_name="y", var_subst="tan(x)"))
        calc.perform_rule(rules.Equation("1 / (b * y ^ 2 + a)", "(1/b)/(y^2+a/b)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("a / b", "(sqrt(a) / sqrt(b))^2"))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("2 * sqrt(a) * sqrt(b)", "2 * sqrt(a * b)"))

        goal03 = file.add_goal("I(a,b,2) = pi/(4*sqrt(a*b)) * (1/a + 1/b)", conds=["a>0", "b>0"])
        proof_of_goal03 = goal03.proof_by_calculation()
        calc = proof_of_goal03.lhs_calc
        s = calc.parse_expr("I(a,b,2)")
        calc.perform_rule(rules.ApplyEquation(goal01.goal, s))
        calc.perform_rule(rules.Simplify())
        s = calc.parse_expr("I(a,b,1)")
        calc.perform_rule(rules.ApplyEquation(goal02.goal, s))
        s = calc.parse_expr("I(a,b,1)")
        calc.perform_rule(rules.ApplyEquation(goal02.goal, s))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("pi / (4 * a ^ (3/2) * sqrt(b)) + pi / (4 * sqrt(a) * b ^ (3/2))",
                                         "pi/(4*sqrt(a*b)) * (1/a + 1/b)"))

        goal04 = file.add_goal("I(a,b,3) = pi/(16*sqrt(a*b)) * (3/a^2 + 3/b^2 + 2/(a*b))", conds=["a>0", "b>0"])
        proof_of_goal04 = goal04.proof_by_calculation()
        calc = proof_of_goal04.lhs_calc
        s = calc.parse_expr(" I(a,b,3)")
        calc.perform_rule(rules.ApplyEquation(goal01.goal, s))
        calc.perform_rule(rules.Simplify())
        s = calc.parse_expr("I(a,b,2)")
        calc.perform_rule(rules.ApplyEquation(goal03.goal, s))
        calc.perform_rule(rules.ApplyEquation(goal03.goal, s))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("pi / (16 * a ^ (3/2) * sqrt(b)) * (1 / a + 1 / b) + pi / (16 * sqrt(a) * b ^ (3/2)) * (1 / a + 1 / b) + pi / (8 * a ^ (5/2) * sqrt(b)) + pi / (8 * sqrt(a) * b ^ (5/2))",
                                         "pi/(16*sqrt(a*b)) * (3/a^2 + 3/b^2 + 2/(a*b))"))

        goal05 = file.add_goal("J(a,y,n) = y/(y^2+a^2)^n+2*n*J(a,y,n)-(2*n*a^2)*J(a,y,n+1)", conds=["a>0", "y>=0", "n>0", "isInt(n)"])
        proof_of_goal05 = goal05.proof_by_calculation()
        calc = proof_of_goal05.lhs_calc
        calc.perform_rule(rules.ExpandDefinition("J"))
        calc.perform_rule(rules.IntegrationByParts(u="1/(x^2+a^2)^n", v="x"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("2 * n * (INT x:[0,y]. x ^ 2 * (a ^ 2 + x ^ 2) ^ (-n - 1))",
                                         "n * (INT x:[0,y]. 2 * x ^ 2 / (a ^ 2 + x ^ 2) ^ (n + 1))"))
        calc.perform_rule(rules.Equation("2 * x ^ 2 / (a ^ 2 + x ^ 2) ^ (n + 1)",
                                         "2 * (x ^ 2 + a ^ 2) / (a ^ 2 + x ^ 2) ^ (n + 1) - 2 * a ^ 2 / (a ^ 2 + x ^ 2) ^ (n + 1)"))
        calc.perform_rule(rules.Equation("(INT x:[0,y]. 2 * (x ^ 2 + a ^ 2) / (a ^ 2 + x ^ 2) ^ (n + 1) - 2 * a ^ 2 / (a ^ 2 + x ^ 2) ^ (n + 1))",
                                         "2 * (INT x:[0,y]. 1 / (x ^ 2 + a ^ 2) ^ n) - 2 * a ^ 2 *(INT x:[0,y]. 1 / (x ^ 2 + a ^ 2) ^ (n + 1))"))
        calc.perform_rule(rules.OnSubterm(rules.FoldDefinition("J")))
        calc.perform_rule(rules.Equation("n * (2 * J(a,y,n) - 2 * a ^ 2 * J(a,y,n + 1)) + y * (a ^ 2 + y ^ 2) ^ -n",
                                         "y/(y^2+a^2)^n+2*n*J(a,y,n)-(2*n*a^2)*J(a,y,n+1)"))

        goal06 = file.add_goal("J(a,y,n+1) = y/(2*n*a^2*(y^2+a^2)^n) + (2*n-1)/(2*n*a^2) * J(a,y,n)", conds=["a>0", "y>=0", "n>0", "isInt(n)"])
        proof_of_goal06 = goal06.proof_by_rewrite_goal(begin = goal05)
        calc = proof_of_goal06.begin
        calc.perform_rule(rules.SolveEquation("J(a,y,n+1)"))
        calc.perform_rule(rules.Equation("(y * (a ^ 2 + y ^ 2) ^ -n + 2 * n * J(a,y,n) - J(a,y,n)) / (2 * a ^ 2 * n)",
                                         "y/(2*n*a^2*(y^2+a^2)^n) + (2*n-1)/(2*n*a^2) * J(a,y,n)"))

        # self.checkAndOutput(file)

    def testCatalanConstant01(self):
        # Reference:
        # Inside interesting integrals, Section 5.1, example #1
        file = compstate.CompFile("interesting", 'CatalanConstant01')
        # Define Catalan's constant
        file.add_definition("G = SUM(n, 0, oo, (-1)^n / (2*n+1)^2)")

        goal = file.add_goal("(INT x:[0, 1]. atan(x) / x) = G")

        goal1 = goal.add_subgoal("1", "converges(SUM(n, 0, oo, INT x:[0,1]. x ^ (2 * n) / (2 * n + 1)))")
        proof = goal1.proof_by_calculation()
        calc = proof.arg_calc
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertTrue(proof.is_finished())

        # Evaluate integral of atan(x) / x
        proof_of_goal = goal.proof_by_calculation()
        calc = proof_of_goal.lhs_calc
        calc.perform_rule(rules.SeriesExpansionIdentity(old_expr="atan(x)", index_var="n"))
        s1 = calc.parse_expr("x ^ (2 * n + 1)")
        s2 = calc.parse_expr("x^ (2*n) * x")
        calc.perform_rule(rules.Equation(s1, s2))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IntSumExchange())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        calc = proof_of_goal.rhs_calc
        calc.perform_rule(rules.ExpandDefinition("G"))
        goal.print_entry()
        assert goal.is_finished()

    def testCatalanConstant02(self):
        # Reference:
        # Inside interesting integrals, Section 5.1, example #2 and example #3
        file = compstate.CompFile("interesting", 'CatalanConstant02')

        goal = file.add_goal("(INT x:[0, oo]. log(x + 1) / (x ^ 2 + 1)) = pi / 4 * log(2) + G")

        # Evaluate I(k)
        goal1 = goal.add_subgoal("1", "(INT x:[1,oo]. log(x) / (x^k)) = 1/(k-1)^2", conds=["k > 1"])
        proof_of_goal1 = goal1.proof_by_calculation()
        calc = proof_of_goal1.lhs_calc
        u = parser.parse_expr("log(x)")
        v = parser.parse_expr("(x^(1-k)) / (1-k)")
        calc.perform_rule(rules.ElimInfInterval())
        calc.perform_rule(rules.IntegrationByParts(u=u, v=v))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        calc = proof_of_goal1.rhs_calc
        s1 = parser.parse_expr("1 / (k-1)^2")
        s2 = parser.parse_expr("1 / (-k+1)^2")
        calc.perform_rule(rules.Equation(s1, s2))
        self.assertTrue(goal1.is_finished())

        goal2 = goal.add_subgoal("2", "converges(SUM(n, 0, oo, INT x:[1,oo]. x ^ (-(2 * n) - 2) * log(x)))")
        proof = goal2.proof_by_calculation()
        calc = proof.arg_calc
        s1 = calc.parse_expr("x ^ (-(2 * n) - 2) * log(x)")
        s2 = calc.parse_expr("log(x) / x^(2*n+2)")
        calc.perform_rule(rules.Equation(s1, s2))
        source = calc.parse_expr("INT x:[1,oo]. log(x) / x ^ (2 * n + 2)")
        calc.perform_rule(rules.ApplyEquation("1", source))
        calc.perform_rule(rules.Simplify())
        self.assertTrue(proof.is_finished())

        goal5 = goal.add_subgoal("3", "(INT x:[1,oo]. log(x) / (x^2+1)) = G")
        proof_of_goal5 = goal5.proof_by_calculation()
        calc = proof_of_goal5.lhs_calc
        s1 = parser.parse_expr("log(x)/(x^2+1)")
        s2 = parser.parse_expr("log(x) * (x^-2) * (1 + (1/x^2))^-1")
        calc.perform_rule(rules.Equation(s1, s2))
        calc.perform_rule(rules.SeriesExpansionIdentity(old_expr="(1+(1/x^2))^-1"))
        s1 = parser.parse_expr("log(x) * x ^ (-2) * SUM(n, 0, oo, (-1) ^ n * (1 / x ^ 2) ^ n)")
        s2 = parser.parse_expr("SUM(n, 0, oo, (-1) ^ n * ((1 / x ^ 2) ^ n) * log(x) * x ^ -2)")
        calc.perform_rule(rules.Equation(s1, s2))
        calc.perform_rule(rules.IntSumExchange())
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("x ^ (-(2 * n) - 2) * log(x)")
        s2 = calc.parse_expr("log(x) / x^(2*n+2)")
        calc.perform_rule(rules.Equation(s1, s2))
        source = calc.parse_expr("(INT x:[1,oo]. log(x) / x ^ (2 * n + 2))")
        calc.perform_rule(rules.ApplyEquation("1", source))
        calc.perform_rule(rules.Simplify())
        calc = proof_of_goal5.rhs_calc
        calc.perform_rule(rules.ExpandDefinition("G"))
        self.assertTrue(goal5.is_finished())

        proof_of_goal6 = goal.proof_by_calculation()
        calc = proof_of_goal6.lhs_calc
        calc.perform_rule(rules.SplitRegion("1"))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        s1 = parser.parse_expr("x + 1")
        s2 = parser.parse_expr("x * (1 + 1 / x)")
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = parser.parse_expr("log(x * (1 + 1 / x))")
        s2 = parser.parse_expr("log(x) + log(1 + 1 / x)")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        s1 = parser.parse_expr("(log(x) + log(1 + 1 / x)) / (x ^ 2 + 1)")
        s2 = parser.parse_expr("log(x) / (x ^ 2 + 1) + log(1 + 1 / x) / (x ^ 2 + 1)")
        calc.perform_rule(rules.Equation(s1, s2))
        calc.perform_rule(rules.Simplify())
        source = calc.parse_expr("INT x:[1,oo]. log(x) / (x ^ 2 + 1)")
        calc.perform_rule(rules.ApplyEquation("3", source))
        calc.perform_rule(rules.Substitution(var_name="u", var_subst="1 / x"))
        s1 = calc.parse_expr("u ^ 2 * (1 / u ^ 2 + 1)")
        s2 = calc.parse_expr("u ^ 2 + 1")
        calc.perform_rule(rules.Equation(s1, s2))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("pi * log(2) / 4")
        s2 = calc.parse_expr("pi / 4 * log(2)")
        calc.perform_rule(rules.Equation(s1, s2))

        goal.print_entry()
        self.assertTrue(goal.is_finished())


    def testCatalanConstant03(self):
        # Reference:
        # Inside interesting integrals, Section 5.1, example #4
        # Notice: This proof is different from the book's.
        # If we use the series method mentioned in the textbook, the theorem is only valid when a>b>0;
        # In fact, the theorem can be proved to hold for a>0 and b>0 using the method below.

        file = compstate.CompFile("interesting", 'CatalanConstant03')

        goal = file.add_goal("(INT x:[0,pi]. x*sin(x)/(a+b*(cos(x))^2)) = pi/sqrt(a*b) *atan(sqrt(b/a))", conds=["a>0", "b>0"])

        goal.add_definition("I(a, b) = (INT x:[0,pi]. x*sin(x)/(a+b*(cos(x))^2))", conds=["a>0", "b>0"])
        goal01 = goal.add_subgoal("1", "I(a,b) = (INT x:[0,pi]. (pi-x)*sin(x)/(a+b*(cos(x))^2))", conds=["a>0","b>0"])
        proof = goal01.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.ExpandDefinition("I"))
        calc.perform_rule(rules.Substitution(var_name="x", var_subst="pi-x"))
        calc.perform_rule(rules.Equation("sin(x) * (-x + pi) / (b * cos(x) ^ 2 + a)", "(pi-x)*sin(x)/(a+b*(cos(x))^2)"))
        assert goal01.is_finished()
        proof = goal.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.OnSubterm(rules.FoldDefinition("I")))
        print(proof)
        calc.perform_rule(rules.Equation("I(a,b)", "1/2 * (I(a,b) + I(a,b))"))
        calc.perform_rule(rules.OnCount(rules.ExpandDefinition("I"), 1))
        source = calc.parse_expr("I(a,b)")
        calc.perform_rule(rules.ApplyEquation("1", source))
        calc.perform_rule(rules.Equation("(INT x:[0,pi]. x * sin(x) / (b * cos(x) ^ 2 + a)) + (INT x:[0,pi]. (pi - x) * sin(x) / (a + b * cos(x) ^ 2))",
                                         "(INT x:[0,pi]. x * sin(x)/(a+b*cos(x)^2) + (pi-x) * sin(x)/(a+b*cos(x)^2))"))
        calc.perform_rule(rules.Equation("x * sin(x)/(a+b*cos(x)^2) + (pi-x) * sin(x)/(a+b*cos(x)^2)",
                                         "pi * sin(x)/(a+b *cos(x)^2)"))
        calc.perform_rule(rules.Substitution(var_name="u", var_subst="cos(x)"))
        calc.perform_rule(rules.Substitution(var_name="x", var_subst="sqrt(b/a) * u"))
        calc.perform_rule(rules.Equation("a * x ^ 2 + a", "a*(x^2+1)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.ApplyIdentity("atan(-(sqrt(b) / sqrt(a)))", "-atan(sqrt(b) / sqrt(a))"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("sqrt(a) * sqrt(b)", "sqrt(a*b)"))
        calc.perform_rule(rules.Equation("sqrt(b) / sqrt(a)", "sqrt(b/a)"))
        goal.print_entry()
        assert goal.is_finished()

    def testLogFunction01(self):
        # Reference:
        # Inside interesting integrals, Section 5.2, example #1
        file = compstate.CompFile("interesting", "LogFunction01")

        goal = file.add_goal("(INT x:[0,1]. log(1 + x) / x) = (pi^2) / 12")

        # Convergence
        goal01 = goal.add_subgoal("1", "converges(SUM(n,0,oo, (INT x:[0,1]. x ^ n / (n + 1))))")
        proof = goal01.proof_by_calculation()
        calc = proof.arg_calc
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertTrue(proof.is_finished())

        # Main result
        proof_of_goal = goal.proof_by_calculation()
        calc = proof_of_goal.lhs_calc
        calc.perform_rule(rules.SeriesExpansionIdentity(old_expr="log(1+x)"))
        calc.perform_rule(rules.Equation(
            "SUM(n,0,oo,(-1) ^ n * x ^ (n + 1) / (n + 1)) / x",
            "SUM(n,0,oo,(-1) ^ n * x ^ (n + 1) / (n + 1) * (1/x))"))
        calc.perform_rule(rules.IntSumExchange())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.SeriesEvaluationIdentity())

        goal.print_entry()
        assert goal.is_finished()

    def testLogFunction02(self):
        # Reference:
        # Inside interesting integrals, Section 5.2, example #2 (5.2.4)
        file = compstate.CompFile("interesting", 'LogFunction02')
        goal01 = file.add_goal("-log(1-x) - log(1+x) = \
                -SUM(k,0,oo,(-1)^k*(-x)^(k+1) / (k+1))-SUM(k,0,oo,(-1)^k*x^(k+1)/(k+1))", conds=["abs(x) < 1"])
        proof_of_goal01 = goal01.proof_by_calculation()
        calc = proof_of_goal01.lhs_calc
        calc.perform_rule(rules.SeriesExpansionIdentity(old_expr="log(1-x)", index_var='k'))
        calc.perform_rule(rules.SeriesExpansionIdentity(old_expr="log(1+x)", index_var='k'))
        assert goal01.is_finished()

        goal02 = file.add_goal("x / (-(x ^ 2) + 1) = \
                1/2 * SUM(k, 0, oo, x ^ k) - 1/2 * SUM(k, 0, oo, x ^ k * (-1) ^ k)",conds = ["abs(x) < 1"])
        proof_of_goal02 = goal02.proof_by_rewrite_goal(begin = goal01)
        calc = proof_of_goal02.begin
        calc.perform_rule(rules.DerivEquation('x'))
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("1 / (-x + 1) - 1 / (x + 1)")
        s2 = calc.parse_expr("2 * (x / (1-x^2))")
        calc.perform_rule(rules.Equation(s1, s2))
        calc.perform_rule(rules.SolveEquation(calc.parse_expr("x / (1-x^2)")))
        s1 = calc.parse_expr("(-1) ^ k * (-x) ^ k")
        s2 = calc.parse_expr("x ^ k")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        calc.perform_rule(rules.OnLocation(rules.ExpandPolynomial(), "1"))
        s1 = calc.parse_expr("-(SUM(k, 0, oo, x ^ k * (-1) ^ k) / 2) + SUM(k, 0, oo, x ^ k) / 2")
        s2 = calc.parse_expr("1/2 * SUM(k, 0, oo, x ^ k) - 1/2 * SUM(k, 0, oo, x ^ k * (-1) ^ k)")
        calc.perform_rule(rules.Equation(s1, s2))
        assert goal02.is_finished()

        goal = file.add_goal("converges(SUM(k, 0, oo, INT y:[0,1]. -(y ^ k * log(y))))")
        proof = goal.proof_by_calculation()
        calc = proof.arg_calc
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertTrue(goal.is_finished())

        goal03 = file.add_goal("(INT x:[0, pi/2]. cos(x)/sin(x) * log(1/cos(x))) = pi^2/24")
        proof_of_goal03 = goal03.proof_by_calculation()
        calc = proof_of_goal03.lhs_calc
        e = calc.parse_expr("cos(x)")
        calc.perform_rule(rules.Substitution(var_name='t', var_subst=e))
        calc.perform_rule(rules.Simplify())
        e = calc.parse_expr("t")
        calc.perform_rule(rules.Substitution(var_name='y', var_subst=e))
        s1 = calc.parse_expr("y * log(y) / (-(y ^ 2) + 1)")
        s2 = calc.parse_expr("log(y) * (y / (-(y ^ 2) + 1))")
        calc.perform_rule(rules.Equation(s1, s2))
        s = calc.parse_expr("y / (-(y ^ 2) + 1)")
        calc.perform_rule(rules.ApplyEquation(goal02.goal, s))
        calc.perform_rule(rules.OnLocation(rules.ExpandPolynomial(), "0"))
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("log(y) * SUM(k, 0, oo, y ^ k)")
        s2 = calc.parse_expr("SUM(k, 0, oo, log(y) * y ^ k)")
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = calc.parse_expr("log(y) * SUM(k, 0, oo, y ^ k * (-1) ^ k)")
        s2 = calc.parse_expr("SUM(k, 0, oo, log(y) * y ^ k * (-1) ^ k)")
        calc.perform_rule(rules.Equation(s1, s2))
        calc.perform_rule(rules.OnSubterm(rules.IntSumExchange()))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.OnSubterm(rules.SeriesEvaluationIdentity()))
        calc.perform_rule(rules.Simplify())
        assert goal03.is_finished()
        # self.checkAndOutput(file)


    def testBernoulliIntegral(self):
        # Reference:
        # Inside interesting integrals, Section 6.1
        file = compstate.CompFile("interesting", "BernoulliIntegral")

        goal1 = file.add_goal("(INT x:[0,1]. x^(c*x^a)) = SUM(k,0,oo,(-c)^k / (k*a+1)^(k+1))", conds=["a > 0", "c != 0"])

        goal11 = goal1.add_subgoal("1", "converges(SUM(k, 0, oo, abs(INT x:[0,1]. (c * x ^ a * log(x)) ^ k / factorial(k))))")
        proof = goal11.proof_by_calculation()
        calc = proof.arg_calc
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("(c*x^a * log(x))^k")
        s2 = calc.parse_expr("(c*x^a)^k * log(x)^k")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        s1 = calc.parse_expr("(c*x^a)^k")
        s2 = calc.parse_expr("c^k * x^a^k")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        assert goal11.is_finished()

        proof_of_goal1 = goal1.proof_by_calculation()
        calc = proof_of_goal1.lhs_calc
        s1 = calc.parse_expr("x^(c*x^a)")
        s2 = calc.parse_expr("exp(log(x^(c*x^a)))")
        calc.perform_rule(rules.Equation(s1, s2))
        calc.perform_rule(rules.SeriesExpansionIdentity(old_expr="exp(log(x^(c*x^a)))", index_var='k'))
        calc.perform_rule(rules.IntSumExchange())
        s1 = calc.parse_expr("log(x^(c*x^a))")
        s2 = calc.parse_expr("(c*x^a) * log(x)")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        s1 = calc.parse_expr("(c*x^a * log(x))^k")
        s2 = calc.parse_expr("(c*x^a)^k * log(x)^k")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        s1 = calc.parse_expr("(c*x^a)^k")
        s2 = calc.parse_expr("c^k * x^a^k")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("c^k * (-1)^k")
        s2 = calc.parse_expr("(-c)^k")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        calc = proof_of_goal1.rhs_calc
        calc.perform_rule(rules.Simplify())
        assert goal1.is_finished()

        goal2 = file.add_goal("(INT x:[0,1]. x^x) = SUM(k, 0, oo, (-1) ^ k * (k + 1) ^ (-k - 1))")
        proof_of_goal2 = goal2.proof_by_calculation()
        calc = proof_of_goal2.lhs_calc
        s1 = parser.parse_expr("x^x")
        s2 = parser.parse_expr("x^(1*x^1)")
        calc.perform_rule(rules.Equation(s1, s2))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        # source = parser.parse_expr("INT x:[0,1]. x ^ (1 * x ^ 1)")
        # calc.perform_rule(rules.ApplyEquation(goal1.goal, source))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        assert goal2.is_finished()

        goal3 = file.add_goal("(INT x:[0,1]. x^(-x)) = SUM(k, 0, oo, (k + 1) ^ (-k - 1))")
        proof_of_goal3 = goal3.proof_by_calculation()
        calc = proof_of_goal3.lhs_calc
        s1 = parser.parse_expr("x^(-x)")
        s2 = parser.parse_expr("x^(-1*x^1)")
        calc.perform_rule(rules.Equation(s1, s2))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        assert goal3.is_finished()

        goal4 = file.add_goal("(INT x:[0,1]. x^(x^2)) = SUM(k, 0, oo, (-1) ^ k * (2 * k + 1) ^ (-k - 1))")
        proof_of_goal4 = goal4.proof_by_calculation()
        calc = proof_of_goal4.lhs_calc
        s1 = parser.parse_expr("x^(x^2)")
        s2 = parser.parse_expr("x^(1*x^2)")
        calc.perform_rule(rules.Equation(s1, s2))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        assert goal4.is_finished()

        goal5 = file.add_goal("(INT x:[0,1]. x^(sqrt(x))) = SUM(k, 0, oo, (-1) ^ k * (2 / (k + 2)) ^ (k+1))")
        proof_of_goal5 = goal5.proof_by_calculation()
        calc = proof_of_goal5.lhs_calc
        s1 = parser.parse_expr("x^(sqrt(x))")
        s2 = parser.parse_expr("x^(1*x^(1/2))")
        calc.perform_rule(rules.Equation(s1, s2))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("(k / 2 + 1)")
        s2 = calc.parse_expr("(2/(k+2)) ^ (-1)")
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = calc.parse_expr("(2 / (k + 2)) ^ (-1) ^ (-k - 1)")
        s2 = calc.parse_expr("(2/(k+2))^(k+1)")
        calc.perform_rule(rules.OnLocation(rules.ApplyIdentity(s1, s2), "0.1"))
        assert goal5.is_finished()

    def testAhmedIntegral(self):
        # Reference:
        # Inside interesting integrals, Section 6.2
        file = compstate.CompFile("interesting", "AhmedIntegral")

        goal = file.add_goal("(INT x:[0,1]. atan(sqrt(2+x^2)) / ((1+x^2)*sqrt(2+x^2))) = 5*pi^2/96")

        # Define I
        goal.add_definition("I(u) = (INT x:[0,1]. atan(u * sqrt(2+x^2)) / ((1+x^2)*sqrt(2+x^2)))", conds=["u > 0"])

        goal01 = goal.add_subgoal("1", "I(1) = INT x:[0,1]. atan(sqrt(x ^ 2 + 2)) / ((x ^ 2 + 1) * sqrt(x ^ 2 + 2))")
        proof = goal01.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.ExpandDefinition("I"))
        assert goal01.is_finished()

        goal02 = goal.add_subgoal("2", "(D u. I(u)) = \
                1 / (1 + u^2) * (pi/4 - u / sqrt(1+2*u^2) * atan(u/sqrt(1+2*u^2)))", conds=["u > 0"])
        proof = goal02.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.OnSubterm(rules.ExpandDefinition("I")))
        calc.perform_rule(rules.DerivIntExchange())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation(
            "1 / ((x ^ 2 + 1) * (u ^ 2 * (x ^ 2 + 2) + 1))",
            "1 / (u ^ 2 + 1) * (1 / (1 + x ^ 2) - u ^ 2 / (1 + 2*u^2 + u^2*x^2))"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation(
            "1 / (u ^ 2 * x ^ 2 + 2 * u ^ 2 + 1)",
            "u^(-2) * (x ^ 2 + (2 * u ^ 2 + 1)/u^2) ^ (-1)"))
        calc.perform_rule(rules.Simplify())
        e = parser.parse_expr("y * sqrt(u ^ (-2) * (2 * u ^ 2 + 1))")
        calc.perform_rule(rules.SubstitutionInverse('y', 'x', e))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation(
            "1 / (y ^ 2 * (2 * u ^ 2 + 1) / u ^ 2 + (2 * u ^ 2 + 1) / u ^ 2)",
            "(1 / (y ^ 2 + 1)) * (u ^ 2 / (2 * u ^ 2 + 1))"))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        assert goal02.is_finished()

        goal03 = goal.add_subgoal("3", "(INT u:[1, oo]. D u. I(u)) = pi^2/12 - I(1)")
        proof = goal03.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.OnCount(rules.ExpandDefinition("I"), 1))
        calc.perform_rule(rules.Simplify())
        u = expr.Const(1)
        v = parser.parse_expr("atan(x/sqrt(2+x^2))")
        calc.perform_rule(rules.IntegrationByParts(u=u, v=v))
        calc.perform_rule(rules.Simplify())
        assert goal03.is_finished()

        goal04 = goal.add_subgoal("4", "(INT u:[1,oo]. D u. I(u)) = - (pi^2 / 48) + I(1)")
        proof_of_goal04 = goal04.proof_by_calculation()
        calc = proof_of_goal04.lhs_calc
        source = calc.parse_expr("D u. I(u)")
        calc.perform_rule(rules.ApplyEquation("2", source))
        calc.perform_rule(rules.ExpandPolynomial())
        calc.perform_rule(rules.Simplify())
        e = parser.parse_expr("1/x")
        calc.perform_rule(rules.SubstitutionInverse('x', 'u', e))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("x ^ 3 * (1 / x ^ 2 + 1) * sqrt(2 / x ^ 2 + 1)",
                                         "sqrt((1 + x^2) ^ 2 * (2 + x^2))"))
        calc.perform_rule(rules.Equation("x * sqrt(2 / x ^ 2 + 1)",
                                         "sqrt(x ^ 2 + 2)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("1 / sqrt(x ^ 2 + 2)", "sqrt(x^2+2) ^ (-1)"))
        calc.perform_rule(rules.ApplyIdentity("atan((sqrt(x^2 + 2))^(-1))", "pi/2 - atan(sqrt(x^2 + 2))"))
        calc.perform_rule(rules.ExpandPolynomial())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("atan(sqrt(x ^ 2 + 2)) / (x ^ 2 * sqrt(x ^ 2 + 2) + sqrt(x ^ 2 + 2))",\
                                         "atan(sqrt(x ^ 2 + 2)) / ((x ^ 2 + 1) * sqrt(x ^ 2 + 2))"))
        source = calc.parse_expr("INT x:[0,1]. atan(sqrt(x ^ 2 + 2)) / ((x ^ 2 + 1) * sqrt(x ^ 2 + 2))")
        calc.perform_rule(rules.ApplyEquation("1", source))
        u = expr.Const(1)
        v = parser.parse_expr("atan(x/sqrt(2+x^2))")
        calc.perform_rule(rules.IntegrationByParts(u=u, v=v))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        assert goal04.is_finished()

        proof_of_goal05 = goal.proof_by_rewrite_goal(begin="3")
        calc = proof_of_goal05.begin
        source = calc.parse_expr("(INT u:[1, oo]. D u. I(u))")
        calc.perform_rule(rules.ApplyEquation("4", source))
        calc.perform_rule(rules.SolveEquation(calc.parse_expr("I(1)")))
        calc.perform_rule(rules.OnSubterm(rules.ExpandDefinition("I")))

        assert goal.is_finished()

    def testEulerConstant01(self):
        # Reference:
        # Inside interesting integrals, Section 5.4.1 & 5.4.3
        # Notice that the proof below is different from the book's.
        file = compstate.CompFile("base", "EulerConstant01")

        # Define Euler's Constant
        file.add_definition("G = - (INT x:[0, oo]. exp(-x) * log(x)) ")

        # Prove these two integrals equal
        goal = file.add_goal("(INT x:[0, 1]. (1 - exp(-x)) / x) - (INT x:[1, oo]. exp(-x) / x) = G")
        proof = goal.proof_by_calculation()
        calc = proof.lhs_calc
        u = parser.parse_expr("exp(-x)")
        v = parser.parse_expr("log(x)")
        calc.perform_rule(rules.OnCount(rules.IntegrationByParts(u=u, v=v), 2))
        u = parser.parse_expr("1 - exp(-x)")
        v = parser.parse_expr("log(x)")
        calc.perform_rule(rules.IntegrationByParts(u=u, v=v))
        calc.perform_rule(rules.Simplify())
        calc = proof.rhs_calc
        calc.perform_rule(rules.ExpandDefinition("G"))
        calc.perform_rule(rules.SplitRegion(parser.parse_expr("1")))
        calc.perform_rule(rules.Simplify())
        assert goal.is_finished()

    # TODO: Substitute s for t - sqrt(a*b)/t
    # def testProbability(self):
    #     # Reference:
    #     # Inside interesting integrals, Section 3.7
    #     file = compstate.CompFile("interesting", "Probability")
    #
    #     file.add_definition("I(a, b) = (INT x:[0, oo]. exp(-a*x^2 - b/x^2))", conds=["a>0", "b>=0"])
    #     file.add_definition("J(a, b) = (INT t:[0, oo]. exp(-t^2 - a*b/t^2))", conds=["a>0", "b>=0"])
    #
    #     goal01 = file.add_goal("I(a, b) = 1/sqrt(a) * J(a, b)", conds=["a>0", "b>0"])
    #     proof01 = goal01.proof_by_calculation()
    #     calc = proof01.lhs_calc
    #     calc.perform_rule(rules.ExpandDefinition("I"))
    #     calc.perform_rule(rules.Substitution(var_name="t", var_subst="x*sqrt(a)"))
    #     calc.perform_rule(rules.Simplify())
    #     calc.perform_rule(rules.Equation("-(a * b / t ^ 2) - t ^ 2", "-t^2 - a*b/t^2"))
    #     calc.perform_rule(rules.OnSubterm((rules.FoldDefinition("J"))))
    #
    #     goal02 = file.add_goal("J(a, b) = sqrt(a*b) * (INT t:[0,oo]. exp(-t^2-a*b/t^2)/t^2)", conds=["a>0","b>0"])
    #     proof02 = goal02.proof_by_calculation()
    #     calc = proof02.lhs_calc
    #     calc.perform_rule(rules.ExpandDefinition("J"))
    #     calc.perform_rule(rules.Substitution(var_name="t", var_subst="sqrt(a*b)/t"))
    #     calc.perform_rule(rules.Simplify())
    #     calc.perform_rule(rules.Equation("1 / t ^ 2 * exp(-(a * b / t ^ 2) - t ^ 2)",
    #                                      "exp(-t^2-a*b/t^2)/t^2"))
    #     calc.perform_rule(rules.Equation("sqrt(a)*sqrt(b)", "sqrt(a*b)"))
    #
    #     goal03 = file.add_goal("2*J(a,b) = exp(-2*sqrt(a*b)) * (INT s:[0,oo]. exp(-s^2))", conds=["a>0","b>0"])
    #     proof03 = goal03.proof_by_calculation()
    #     calc = proof03.lhs_calc
    #     calc.perform_rule(rules.Equation("2*J(a, b)", "J(a,b) + J(a, b)"))
    #     calc.perform_rule(rules.OnLocation(rules.ExpandDefinition("J"), "0"))
    #     calc.perform_rule(rules.OnSubterm(rules.ApplyEquation(goal02.goal)))
    #     calc.perform_rule(rules.Equation("(INT t:[0,oo]. exp(-(a * b / t ^ 2) - t ^ 2)) + sqrt(a * b) * (INT t:[0,oo]. exp(-(t ^ 2) - a * b / t ^ 2) / t ^ 2)",
    #                                      "(INT t:[0,oo]. exp(-(a * b / t ^ 2) - t ^ 2) + sqrt(a * b) * exp(-(t ^ 2) - a * b / t ^ 2) / t ^ 2)"))
    #     calc.perform_rule(rules.Substitution(var_name="s", var_subst="t-sqrt(a*b)/t"))
    #     # self.checkAndOutput(file)

    # TODO: Show I(1/a) = 0 for a>1.
    # def testDini(self):
    #     # Reference:
    #     # Inside interesting integrals, Section 3.8
    #     file = compstate.CompFile("interesting", "dini")
    #
    #     file.add_definition("I(a) = (INT x:[0,pi]. log(1-2*a*cos(x)+a^2))", conds=["a>0", "a!=1"])
    #
    #     goal01 = file.add_goal("(D a. I(a)) = "
    #                            "1 / a * (-(abs(a + 1) * (-(2 * a ^ 2) + 2) / ((a + 1) ^ 2 * abs(-a + 1)) * (LIM {z -> oo}. atan(z * abs(a + 1) / abs(-a + 1)))) + pi)")
    #     proof_of_goal01 = goal01.proof_by_calculation()
    #     calc = proof_of_goal01.lhs_calc
    #     calc.perform_rule(rules.OnSubterm(rules.ExpandDefinition("I")))
    #     calc.perform_rule(rules.Simplify())
    #     calc.perform_rule(rules.Equation("(-(2 * cos(x)) + 2 * a)", "1/a*(2*a^2-2*a*cos(x))"))
    #     calc.perform_rule(rules.Simplify())
    #     calc.perform_rule(rules.Equation("(-(2 * a * cos(x)) + 2 * a ^ 2) / (-(2 * a * cos(x)) + a ^ 2 + 1)",
    #                                      "1 - (1-a^2) / (-(2 * a * cos(x)) + a ^ 2 + 1)"))
    #     calc.perform_rule(rules.Simplify())
    #     calc.perform_rule(rules.DefiniteIntegralIdentity())
    #     calc.perform_rule(rules.Simplify())
    #     calc.perform_rule(rules.Substitution(var_name="z", var_subst="tan(x/2)"))
    #     calc.perform_rule(rules.Equation("(z ^ 2 + 1) * (-(2 * a * (-(z ^ 2) + 1) / (z ^ 2 + 1)) + a ^ 2 + 1)",
    #                                      "(1+a)^2*z^2+(1-a)^2"))
    #     calc.perform_rule(rules.Equation("2/((1+a)^2*z^2+(1-a)^2)", "2/(1+a)^2 * 1/(((1-a)/(1+a))^2+z^2)"))
    #     calc.perform_rule(rules.Simplify())
    #     calc.perform_rule(rules.Equation("(-a + 1) ^ 2 / (a + 1) ^ 2 + z ^ 2", "z ^ 2 + ((1-a) / (1+a)) ^ 2"))
    #     calc.perform_rule(rules.DefiniteIntegralIdentity())
    #     calc.perform_rule(rules.Simplify())
    #
    #
    #     goal02 = file.add_goal("I(a)=SKOLEM_CONST(C)", conds=["a>=0", "a<1"])
    #     proof_of_goal02 = goal02.proof_by_rewrite_goal(begin=goal01)
    #     calc = proof_of_goal02.begin
    #     calc.perform_rule(rules.Simplify())
    #     calc.perform_rule(rules.Equation("(-a + 1) * (2 * a + 2)", "-(2 * a ^ 2) + 2"))
    #     calc.perform_rule(rules.Simplify())
    #     calc.perform_rule(rules.IntegralEquation())
    #     calc.perform_rule(rules.IndefiniteIntegralIdentity())
    #     calc.perform_rule(rules.Simplify())
    #
    #     goal03 = file.add_goal("SKOLEM_CONST(C) = 0", conds=["a>=0", "a<1"])
    #     proof_of_goal03 = goal03.proof_by_rewrite_goal(begin=goal02)
    #     calc = proof_of_goal03.begin
    #     calc.perform_rule(rules.VarSubsOfEquation([{'var': 'a', 'expr': "0"}]))
    #     calc.perform_rule(rules.OnSubterm(rules.ExpandDefinition("I")))
    #     calc.perform_rule(rules.Simplify())
    #     calc.perform_rule(rules.SolveEquation(parser.parse_expr("SKOLEM_CONST(C)")))
    #
    #     goal04 = file.add_goal("I(a) = 0", conds=["a>=0", "a<1"])
    #     proof_of_goal04 = goal04.proof_by_rewrite_goal(begin=goal02)
    #     calc = proof_of_goal04.begin
    #     calc.perform_rule(rules.OnSubterm(rules.ApplyEquation(goal03.goal)))
    #
    #     goal05 = file.add_goal("I(1/a) = I(a) - 2*pi*log(a)", conds=["a>1"])
    #     proof_of_goal05 = goal05.proof_by_calculation()
    #     calc = proof_of_goal05.lhs_calc
    #     calc.perform_rule(rules.ExpandDefinition("I"))
    #     calc.perform_rule(rules.Equation("-(2 * cos(x) / a) + 1 / a ^ 2 + 1", "(a^2-2*a*cos(x)+1)/a^2"))
    #     calc.perform_rule(rules.ApplyIdentity("log((a^2-2*a*cos(x)+1)/a^2)", "log(a^2-2*a*cos(x)+1) - log(a^2)"))
    #     calc.perform_rule(rules.Simplify())
    #     calc.perform_rule(rules.DefiniteIntegralIdentity())
    #     calc.perform_rule(rules.Simplify())
    #     calc.perform_rule(rules.Equation("-(2 * a * cos(x)) + a ^ 2 + 1", "1-2*a*cos(x)+a^2"))
    #     calc.perform_rule(rules.OnSubterm(rules.FoldDefinition("I")))
    #
    #     goal06 = file.add_goal("I(a) = 2*pi*log(a)", conds=["a>1"])
    #     proof_of_goal06 = goal06.proof_by_rewrite_goal(begin=goal05)
    #     calc = proof_of_goal06.begin
    #     calc.perform_rule(rules.OnLocation(rules.ApplyEquation(goal04.goal), "0"))
    #
    #     # self.checkAndOutput(file)

    def testChapter3Practice01(self):
        # Reference:
        # Inside interesting integrals, Section 3.10, C3.1
        file = compstate.CompFile("interesting", "Chapter3Practice01")

        goal = file.add_goal("(INT x:[0, oo]. log(1 + a ^ 2 * x ^ 2) / (b ^ 2 + x ^ 2)) = (pi / b) * log(1 + a * b)",
                             conds=["a > 0", "b > 0"])

        goal.add_definition("I(a, b) = (INT x:[0, oo]. log(1 + a ^ 2 * x ^ 2) / (b ^ 2 + x ^ 2))",
                             conds=["a >= 0", "b > 0"])

        goal01 = goal.add_subgoal("1", "(D a. I(a,b)) = pi / (1 + a * b)", conds=["a > 0", "b > 0"])
        proof_of_goal01 = goal01.proof_by_calculation()
        calc = proof_of_goal01.lhs_calc
        calc.perform_rule(rules.OnSubterm(rules.ExpandDefinition("I")))
        calc.perform_rule(rules.DerivIntExchange())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("x ^ 2 / ((b ^ 2 + x ^ 2) * (a ^ 2 * x ^ 2 + 1))",
                                         "1 / (1 - a ^ 2 * b ^ 2) * (1 / (1 + a ^ 2 * x ^ 2) - b ^ 2 / (b ^ 2 + x ^ 2))"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Substitution("t", parser.parse_expr("a * x")))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation(None, "pi / (1 + a * b)"))
        assert goal01.is_finished()

        goal02 = goal.add_subgoal("2", "I(a, b) = (pi / b) * log(1 + a * b) + SKOLEM_FUNC(C(b))", conds=["a > 0", "b > 0"])
        proof_of_goal02 = goal02.proof_by_rewrite_goal(begin="1")
        calc = proof_of_goal02.begin
        calc.perform_rule(rules.IntegralEquation())
        calc.perform_rule(rules.Substitution(var_name="u", var_subst="1 + a * b"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.IndefiniteIntegralIdentity())
        calc.perform_rule(rules.ReplaceSubstitution())
        calc.perform_rule(rules.Equation("abs(1 + a * b)", "(1 + a * b)"))
        assert goal02.is_finished()

        goal03 = goal.add_subgoal("3", "I(0, b) = 0", conds=['b > 0'])
        proof_of_goal03 = goal03.proof_by_calculation()
        calc = proof_of_goal03.lhs_calc
        calc.perform_rule(rules.ExpandDefinition("I"))
        calc.perform_rule(rules.Simplify())
        assert goal03.is_finished()

        goal04 = goal.add_subgoal("4", "SKOLEM_FUNC(C(b)) = 0", conds=['b > 0'])
        proof_of_goal04 = goal04.proof_by_rewrite_goal(begin="2")
        calc = proof_of_goal04.begin
        calc.perform_rule(rules.LimitEquation("a", parser.parse_expr("0")))
        calc.perform_rule(rules.Simplify())
        source = calc.parse_expr("I(0,b)")
        calc.perform_rule(rules.ApplyEquation("3", source))
        calc.perform_rule(rules.SolveEquation("SKOLEM_FUNC(C(b))"))
        assert goal04.is_finished()

        proof_of_goal05 = goal.proof_by_rewrite_goal(begin="2")
        calc = proof_of_goal05.begin
        s = calc.parse_expr("SKOLEM_FUNC(C(b))")
        calc.perform_rule(rules.ApplyEquation("4", s))
        calc.perform_rule(rules.OnSubterm(rules.ExpandDefinition("I")))
        calc.perform_rule(rules.Simplify())

        assert goal.is_finished()

    def testChapter3Practice02(self):
        # Reference:
        # Inside interesting integrals, Section C3.10, C3.2
        file = compstate.CompFile("interesting", "Chapter3Practice02")

        goal = file.add_goal("(INT x:[-oo, oo]. cos(a * x) / (b ^ 2 - x ^ 2)) = pi * sin(a * b) / b",
                             conds=["a > 0", "b > 0", "b!=x"])
        proof = goal.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Equation("b ^ 2 - x ^ 2",
                                         "(b + x) * (b - x)"))
        calc.perform_rule(rules.Equation("cos(a * x) / ((b + x) * (b - x))",
                                         "(1 / (2 * b)) * (cos(a * x) / (b + x) + cos(a * x) / (b - x))"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Substitution(var_name="u", var_subst="b + x"))
        calc.perform_rule(rules.OnCount(rules.Substitution(var_name="u", var_subst="b - x"), 2))
        calc.perform_rule(rules.Equation("a * (-b + u)", "-(a * (b - u))"))
        calc.perform_rule(rules.ApplyIdentity("cos(-(a * (b - u)))", "cos(a * (b - u))"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("cos(a * (b - u))", "cos(a * b - a * u)"))
        calc.perform_rule(
            rules.ApplyIdentity("cos(a * b - a * u)", "cos(a * b) * cos(a * u) + sin(a * b) * sin(a * u)"))
        calc.perform_rule(rules.Equation("(cos(a * b) * cos(a * u) + sin(a * b) * sin(a * u)) / u",
                                         "cos(a * b) * cos(a * u) / u + sin(a * b) * sin(a * u) / u"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("(INT u:[-oo,oo]. cos(a * u) / u)", "0"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.SplitRegion("0"))
        calc.perform_rule(rules.Substitution(var_name="u", var_subst="-u"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        assert goal.is_finished()

    def testChapter3Practice03(self):
        # Reference:
        # Inside interesting integrals, Section C3.10, C3.3
        file = compstate.CompFile("interesting", "Chapter3Practice03")

        goal = file.add_goal(
            "(INT x:[-oo, oo]. cos(a * x) / (b ^ 4 - x ^ 4)) = pi * (exp(-(a * b)) + sin(a * b)) / (2 * b ^ 3)",
            conds=["a > 0", "b > 0", "b != x"])
        proof_of_goal = goal.proof_by_calculation()
        calc = proof_of_goal.lhs_calc
        calc.perform_rule(rules.Equation("b ^ 4 - x ^ 4", "(b ^ 2 + x ^ 2) * (b ^ 2- x ^ 2)"))
        calc.perform_rule(rules.Equation("cos(a * x) / ((b ^ 2 + x ^ 2) * (b ^ 2 - x ^ 2))",
                                         "(1 / (2 * b ^ 2)) * (cos(a * x) / (b ^ 2 + x ^ 2) + cos(a * x) / (b ^ 2 - x ^ 2))"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.SplitRegion("0"))
        calc.perform_rule(rules.Substitution(var_name="x", var_subst="-x"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("b ^ 2 + x ^ 2", "x ^ 2 + b ^ 2"))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation(None, "pi * (exp(-(a * b)) + sin(a * b)) / (2 * b ^ 3)"))

        # self.checkAndOutput(file)

    def testChapter3Practice04(self):
        # Reference:
        # Inside interesting integrals, Section 3.10, C3.4
        file = compstate.CompFile("interesting", "Chapter3Practice04")

        goal = file.add_goal("(INT x:[0, oo]. x * sin(a * x) / (x ^ 2 - b ^ 2)) = pi / 2 * cos(a * b)",
                             conds=["a > 0", "b > 0", "b != x"])

        goal01 = goal.add_subgoal("1",
            "(INT x:[0, oo]. x * sin(a * x) / (x ^ 2 - b ^ 2)) = \
            1/2 * (INT x:[-oo, oo]. x * sin(a * x) / (x ^ 2 - b ^ 2))",\
            conds=["x != b", "a > 0", "b > 0"])
        proof_of_goal01 = goal01.proof_by_calculation()
        calc = proof_of_goal01.lhs_calc
        calc.perform_rule(rules.Simplify())
        calc = proof_of_goal01.rhs_calc
        calc.perform_rule(rules.SplitRegion("0"))
        calc.perform_rule(rules.Substitution(var_name="x", var_subst="-x"))
        calc.perform_rule(rules.Simplify())
        assert goal01.is_finished()

        proof_of_goal02 = goal.proof_by_rewrite_goal(begin="1")
        calc = proof_of_goal02.begin
        calc.perform_rule(rules.OnCount(rules.Equation("x ^ 2 - b ^ 2", "(x + b) * (x - b)"), 2))
        calc.perform_rule(rules.Equation("x * sin(a * x) / ((x + b) * (x - b))",
                                         "-x * sin(a * x) / ((b - x) * (b + x))"))
        calc.perform_rule(rules.Equation("-x * sin(a * x) / ((b - x) * (b + x))",
                                         "-1/(2 * b) * (x * sin(a * x) / (b + x) + x * sin(a * x) / (b - x))"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.OnCount(rules.Substitution(var_name="u", var_subst="b + x"), 2))
        calc.perform_rule(rules.OnCount(rules.Substitution(var_name="u", var_subst="b - x"), 3))
        calc.perform_rule(rules.Equation("sin(a * (-b + u))", "sin(-(a * (b - u)))"))
        calc.perform_rule(rules.ApplyIdentity("sin(-(a * (b - u)))", "-sin(a * (b - u))"))
        calc.perform_rule(rules.Equation("(-b + u) * -sin(a * (b - u))", "(b - u) * sin(a * (b - u))"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("(INT u:[-oo,oo]. (b - u) * sin(a * (b - u)) / u)",
                                         "(INT u:[-oo,oo]. ((b - u) / u) * sin(a * (b - u)))"))
        calc.perform_rule(rules.Equation("((b - u) / u) * sin(a * (b - u))",
                                         "(b / u - 1) * sin(a * b - a * u)"))
        calc.perform_rule(rules.Equation("(b / u - 1) * sin(a * b - a * u)",
                                         "(b / u) * sin(a * b - a * u) - sin(a * b - a * u)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.OnCount(rules.Substitution(var_name="s", var_subst="a * b - a * u"), 3))
        calc.perform_rule(rules.OnCount(rules.SplitRegion("0"), 3))
        calc.perform_rule(rules.OnCount(rules.Substitution(var_name="s", var_subst="-s"), 3))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(
            rules.ApplyIdentity("sin(a * b - a * u)", "sin(a * b) * cos(a * u) - cos(a * b) * sin(a * u)"))
        calc.perform_rule(rules.Equation("(sin(a * b) * cos(a * u) - cos(a * b) * sin(a * u)) / u",
                                         "sin(a * b) * cos(a * u) / u - cos(a * b) * sin(a * u) / u"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.OnCount(rules.SplitRegion("0"), 3))
        calc.perform_rule(rules.OnCount(rules.Substitution(var_name="u", var_subst="-u"), 3))
        calc.perform_rule(rules.OnCount(rules.SplitRegion("0"), 2))
        calc.perform_rule(rules.OnCount(rules.Substitution(var_name="u", var_subst="-u"), 2))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())

        assert goal.is_finished()

    def testChapter3Practice05(self):
        # Reference:
        # Inside interesting integrals, Section 3.10, C3.5
        file = compstate.CompFile("interesting", "Chapter3Practice05")

        file.add_definition("I(a, b) = (INT x:[0, oo]. cos(a * x) * sin(b * x) / x)", conds=["a > 0", "b > 0"])

        goal01 = file.add_goal("I(a, b) = 1/2 * (INT x:[0, oo]. sin((b + a) * x) / x) + 1/2 * (INT x:[0, oo]. sin((b - a) * x) / x)",
                               conds=["a > 0", "b > 0"])
        proof_of_goal01 = goal01.proof_by_calculation()
        calc = proof_of_goal01.lhs_calc
        calc.perform_rule(rules.ExpandDefinition("I"))
        calc.perform_rule(rules.ApplyIdentity("cos(a * x) * sin(b * x)", "1/2 * (sin(b * x + a * x) - sin(a * x - b * x))"))
        calc.perform_rule(rules.Equation("1/2 * (sin(b * x + a * x) - sin(a * x - b * x)) / x",
                                         "1/2 * sin((b + a) * x) / x - 1/2 * sin(-((b - a) * x)) / x"))
        calc.perform_rule(rules.ApplyIdentity("sin(-((b - a) * x))", "-sin((b - a) * x)"))
        calc.perform_rule(rules.Simplify())
        calc = proof_of_goal01.rhs_calc
        calc.perform_rule(rules.Simplify())

        # Case a < b
        goal02 = file.add_goal("I(a, b) = pi / 2", conds=["b - a > 0", "a > 0", "b > 0"])
        proof_of_goal02 = goal02.proof_by_calculation()
        calc = proof_of_goal02.lhs_calc
        s = calc.parse_expr("I(a,b)")
        calc.perform_rule(rules.ApplyEquation(goal01.goal, s))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        # Case a > b
        goal03 = file.add_goal("I(a, b) = 0", conds=["b - a < 0", "a > 0", "b > 0"])
        proof_of_goal03 = goal03.proof_by_calculation()
        calc = proof_of_goal03.lhs_calc
        s = calc.parse_expr("I(a,b)")
        calc.perform_rule(rules.ApplyEquation(goal01.goal, s))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        # Case a = b
        goal04 = file.add_goal("I(a, b) = pi / 4", conds=["b - a = 0", "a > 0", "b > 0"])
        proof_of_goal04 = goal04.proof_by_calculation()
        calc = proof_of_goal04.lhs_calc
        s = calc.parse_expr("I(a,b)")
        calc.perform_rule(rules.ApplyEquation(goal01.goal, s))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        # self.checkAndOutput(file)

    def testChapter3Practice06(self):
        # Reference:
        # Inside interesting integrals, Section 3.10, C3.6
        file = compstate.CompFile("interesting", "Chapter3Practice06")

        goal = file.add_goal("(INT x:[-1, 1]. ((1 + x) / (1 - x)) ^ (1/2)) = pi")
        proof_of_goal = goal.proof_by_calculation()
        calc = proof_of_goal.lhs_calc
        calc.perform_rule(rules.SubstitutionInverse("u", "x", "cos(2 * u)"))
        calc.perform_rule(rules.OnCount(rules.ApplyIdentity("cos(2 * u)", "2 * cos(u) ^ 2 - 1"), 1))
        calc.perform_rule(rules.ApplyIdentity("cos(2 * u)", "1 - 2 * sin(u) ^ 2"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.ApplyIdentity("sin(2 * u)", "2 * sin(u) * cos(u)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        assert goal.is_finished()

    def testChapter3Practice07a(self):
        # Reference:
        # Inside interesting integrals, Section 3.10, C3.8
        file = compstate.CompFile("interesting", "Chapter3Practice07")

        goal = file.add_goal("(INT x:[-oo, oo]. x * exp(-x ^ 2 - x)) = -(1/2) * sqrt(pi * sqrt(exp(1)))")

        goal.add_definition("I(a, b) = (INT x:[-oo, oo]. exp(-a * x ^ 2 + b * x))", conds=["a > 0"])

        goal01 = goal.add_subgoal("1", "I(a, b) = exp(b ^ 2 / (4 * a)) * sqrt(pi / a)", conds=["a > 0"])
        proof_of_goal01 = goal01.proof_by_calculation()
        calc = proof_of_goal01.lhs_calc
        calc.perform_rule(rules.ExpandDefinition("I"))
        calc.perform_rule(rules.Equation("-(a * x ^ 2) + b * x",
                                         "b ^ 2 / (4 * a) - a * (x - b / (2 * a)) ^ 2"))
        calc.perform_rule(rules.ApplyIdentity("exp(b ^ 2 / (4 * a) - a * (x - b / (2 * a)) ^ 2)",
                                              "exp(b ^ 2 / (4 * a)) * exp(-a * (x - b / (2 * a)) ^ 2)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Substitution(var_name="y", var_subst="x - b / (2 * a)"))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        goal03 = goal.add_subgoal("2", "(D b. I(a, b)) = (b / (2 * a)) * exp(b ^ 2 / (4 * a)) * sqrt(pi / a)", conds=["a > 0"])
        proof_of_goal03 = goal03.proof_by_calculation()
        calc = proof_of_goal03.lhs_calc
        s = calc.parse_expr("I(a,b)")
        calc.perform_rule(rules.ApplyEquation("1", s))
        calc.perform_rule(rules.Simplify())

        goal04 = goal.add_subgoal("3", "(INT x:[-oo, oo]. x * exp(-(a * x ^ 2) + b * x)) = (b / (2 * a)) * exp(b ^ 2 / (4 * a)) * sqrt(pi / a)",
                                  conds=["a > 0"])
        proof_of_goal04 = goal04.proof_by_rewrite_goal(begin="2")
        calc = proof_of_goal04.begin
        calc.perform_rule(rules.OnSubterm(rules.ExpandDefinition("I")))
        calc.perform_rule(rules.Simplify())

        proof_of_goal06 = goal.proof_by_calculation()
        calc = proof_of_goal06.lhs_calc
        calc.perform_rule(rules.Equation("x * exp(-x ^ 2 - x)", "x * exp(-(1 * x ^ 2) + (-1) * x)"))
        s = calc.parse_expr("INT x:[-oo,oo]. x * exp(-(1 * x ^ 2) + -1 * x)")
        calc.perform_rule(rules.ApplyEquation("3", s))
        calc.perform_rule(rules.Equation(None, "-(1/2) * sqrt(pi * sqrt(exp(1)))"))

        assert goal.is_finished()

    def testChapter3Practice07b(self):
        # Reference:
        # Inside interesting integrals, Section 3.10, C3.8
        file = compstate.CompFile("interesting", "Chapter3Practice07")

        goal = file.add_goal("(INT x:[-oo, oo]. x ^ 2 * exp(-x ^ 2 - x)) = 3/4 * sqrt(pi * sqrt(exp(1)))")

        goal.add_definition("I(a, b) = (INT x:[-oo, oo]. exp(-a * x ^ 2 + b * x))", conds=["a > 0"])

        goal01 = goal.add_subgoal("1", "I(a, b) = exp(b ^ 2 / (4 * a)) * sqrt(pi / a)", conds=["a > 0"])
        proof_of_goal01 = goal01.proof_by_calculation()
        calc = proof_of_goal01.lhs_calc
        calc.perform_rule(rules.ExpandDefinition("I"))
        calc.perform_rule(rules.Equation("-(a * x ^ 2) + b * x",
                                         "b ^ 2 / (4 * a) - a * (x - b / (2 * a)) ^ 2"))
        calc.perform_rule(rules.ApplyIdentity("exp(b ^ 2 / (4 * a) - a * (x - b / (2 * a)) ^ 2)",
                                              "exp(b ^ 2 / (4 * a)) * exp(-a * (x - b / (2 * a)) ^ 2)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Substitution(var_name="y", var_subst="x - b / (2 * a)"))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        goal02 = goal.add_subgoal("2", "(D a. I(a, b)) = -(b ^ 2 / (4 * a ^ 2)) * exp(b ^ 2 / (4 * a)) * sqrt(pi / a) - 1 / (2 * a) * exp(b ^ 2 / (4 * a)) * sqrt(pi / a)", conds=["a > 0"])
        proof_of_goal02 = goal02.proof_by_calculation()
        calc = proof_of_goal02.lhs_calc
        s = calc.parse_expr("I(a,b)")
        calc.perform_rule(rules.ApplyEquation("1", s))
        calc.perform_rule(rules.Simplify())
        calc = proof_of_goal02.rhs_calc
        calc.perform_rule(rules.Simplify())
        assert goal02.is_finished()

        goal05 = goal.add_subgoal("3", "(INT x:[-oo, oo]. x ^ 2 * exp(-(a * x ^ 2) + b * x)) = (b ^ 2 / (4 * a ^ 2)) * exp(b ^ 2 / (4 * a)) * sqrt(pi / a) + 1 / (2 * a) * exp(b ^ 2 / (4 * a)) * sqrt(pi / a)",
                                  conds=["a > 0"])
        proof_of_goal05 = goal05.proof_by_rewrite_goal(begin="2")
        calc = proof_of_goal05.begin
        calc.perform_rule(rules.OnSubterm(rules.ExpandDefinition("I")))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.SolveEquation("(INT x:[-oo,oo]. x ^ 2 * exp(-(a * x ^ 2) + b * x))"))

        proof_of_goal07 = goal.proof_by_calculation()
        calc = proof_of_goal07.lhs_calc
        calc.perform_rule(rules.Equation("x ^ 2 * exp(-x ^ 2 - x)", "x ^ 2 * exp(-(1 * x ^ 2) + (-1) * x)"))
        s = calc.parse_expr("INT x:[-oo,oo]. x ^ 2 * exp(-(1 * x ^ 2) + -1 * x)")
        calc.perform_rule(rules.ApplyEquation("3", s))
        calc.perform_rule(rules.Equation(None, "3/4 * sqrt(pi * sqrt(exp(1)))"))

        assert goal.is_finished()

    def testChapter3Practice08(self):
        # Reference:
        # Inside interesting integrals, Section 3.10, C3.9
        file = compstate.CompFile("interesting", "Chapter3Practice08")

        goal = file.add_goal("(INT x:[0, oo]. sin(m * x) / (x * (a ^ 2 + x ^ 2) ^ 2)) = (pi / (2 * a ^ 4)) * (1 - ((2 + m * a) / 2) * exp(- a * m))",
                             conds=["a > 0", "m > 0"])

        goal01 = goal.add_subgoal("1", "(INT x:[0, oo]. sin(m * x) / (x * (a ^ 2 + x ^ 2))) = (pi * (1 - exp(-a * m))) / (2 * a ^ 2)", conds=["a > 0", "m > 0"])
        proof_of_goal01 = goal01.proof_by_calculation()
        calc = proof_of_goal01.lhs_calc
        calc.perform_rule(rules.DefiniteIntegralIdentity())

        proof_of_goal02 = goal.proof_by_rewrite_goal(begin="1")
        calc = proof_of_goal02.begin
        calc.perform_rule(rules.DerivEquation("a"))
        calc.perform_rule(rules.OnSubterm(rules.DerivIntExchange()))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.SolveEquation("(INT x:[0,oo]. sin(m * x) / (x * (a ^ 2 + x ^ 2) ^ 2))"))
        calc.perform_rule(rules.Equation("-((2 * a ^ 2 * m * pi * exp(-(a * m)) - 4 * a * pi * (-exp(-(a * m)) + 1)) / (8 * a ^ 5))",
                                         "(pi / (2 * a ^ 4)) * (1 - ((2 + m * a) / 2) * exp(- a * m))"))

        assert goal.is_finished()

    def testChapter3Practice09(self):
        # Reference:
        # Inside interesting integrals, Section 3.10, C3.10
        file = compstate.CompFile("interesting", "Chapter3Practice09")

        goal = file.add_goal("(INT x:[0, 1]. x / (a * x + b * (1 - x)) ^ 3) = 1 / (2 * a ^ 2 * b)",
                             conds=["a > 0", "b > 0", "a > b"])

        goal01 = goal.add_subgoal("1", "(INT x:[0, 1]. 1 / (a * x + b * (1 - x)) ^ 2) = 1 / (a * b)",
                                  conds=["a > 0", "b > 0", "a > b"])
        proof_of_goal01 = goal01.proof_by_calculation()
        calc = proof_of_goal01.lhs_calc
        calc.perform_rule(rules.Substitution(var_name="u", var_subst="(a - b) * x + b"))
        calc.perform_rule(rules.Equation("1 / ((a - b) * (b * (-((-b + u) / (a - b)) + 1) + a * (-b + u) / (a - b)) ^ 2)",
                                         "1 / (u ^ 2 * (a - b))"))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("1 / (a - b) * (-(1 / a) + 1 / b)", "1 / (a * b)"))

        proof_of_goal02 = goal.proof_by_rewrite_goal(begin="1")
        calc = proof_of_goal02.begin
        calc.perform_rule(rules.DerivEquation("a"))
        calc.perform_rule(rules.OnSubterm(rules.DerivIntExchange()))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("(b * (-x + 1) + a * x) ^ 3", "(a * x + b * (1 - x)) ^ 3"))
        calc.perform_rule(rules.SolveEquation("(INT x:[0,1]. x / (a * x + b * (1 - x)) ^ 3)"))
        assert goal01.is_finished()
        assert goal.is_finished()

    def testBetaWallis(self):
        # Reference:
        # Inside interesting integrals, Section 4.2
        file = compstate.CompFile("interesting", "BetaWallis")

        # Definition of Beta function
        file.add_definition("B(m, n) = INT x:[0,1]. x^(m-1) * (1-x)^(n-1)", conds=["m > 0", "n > 0"])

        goal01 = file.add_goal("B(m, n) = 2 * (INT x:[0, pi / 2]. cos(x) ^ (2 * m - 1) * sin(x) ^ (2 * n - 1))", conds=["m > 0", "n > 0"])
        proof_of_goal01 = goal01.proof_by_calculation()
        calc = proof_of_goal01.lhs_calc
        calc.perform_rule(rules.ExpandDefinition("B"))
        calc.perform_rule(rules.SubstitutionInverse("t", "x", "cos(t) ^ 2"))
        calc.perform_rule(rules.OnLocation(rules.ApplyIdentity("cos(t) ^ 2", "1 - sin(t) ^ 2"), "0.0.0.1.0.0.0"))
        calc.perform_rule(rules.Simplify())
        assert goal01.is_finished()

        goal = file.add_goal("x - x ^ 2 > 0", conds=["x > 0", "x < 1"])
        proof = goal.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Equation(None, "x * (1 - x)"))
        assert goal.is_finished()

        goal02 = file.add_goal("(INT x:[0, 1]. (x - x ^ 2) ^ n) = factorial(n) ^ 2 / factorial(2 * n + 1)", conds=["n > -1"])
        proof_of_goal02 = goal02.proof_by_calculation()
        calc = proof_of_goal02.lhs_calc
        calc.perform_rule(rules.Equation("x - x ^ 2", "x * (1 - x)"))
        s1 = calc.parse_expr("(x * (1 - x)) ^ n")
        s2 = calc.parse_expr("x ^ n * (1 - x) ^ n")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        s1 = calc.parse_expr("n")
        s2 = calc.parse_expr("n + 1 - 1")
        calc.perform_rule(rules.OnLocation(rules.Equation(s1, s2), "0.0.1"))
        s1 = calc.parse_expr("n")
        s2 = calc.parse_expr("n + 1 - 1")
        calc.perform_rule(rules.OnLocation(rules.Equation(s1, s2), "0.1.1"))
        calc.perform_rule(rules.FoldDefinition("B"))
        s1 = calc.parse_expr("B(n + 1, n + 1)")
        s2 = calc.parse_expr("Gamma(n + 1) * Gamma(n + 1) / Gamma(2 * n + 2)")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("Gamma(n + 1)")
        s2 = calc.parse_expr("factorial(n)")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        s1 = calc.parse_expr("Gamma(2 * n + 2)")
        s2 = calc.parse_expr("factorial(2 * n + 1)")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        assert goal02.is_finished()

        goal03 = file.add_goal("(INT x:[0, 1]. (x - x ^ 2) ^ (1/2)) = pi / 8")
        proof_of_goal03 = goal03.proof_by_calculation()
        calc = proof_of_goal03.lhs_calc
        calc.perform_rule(rules.SubstitutionInverse("u", "x", "1/2 + 1/2 * sin(u)"))
        calc.perform_rule(rules.Equation("1/2 + 1/2 * sin(u) - (1/2 + 1/2 * sin(u)) ^ 2", "1/4 * (1 - sin(u) ^ 2)"))
        calc.perform_rule(rules.ApplyIdentity("1 - sin(u) ^ 2", "cos(u) ^ 2"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        assert goal03.is_finished()

        goal04 = file.add_goal("(INT x:[0, 1]. (x - x ^ 2) ^ (1/2)) = (factorial(1/2) ^ 2) / 2")
        proof_of_goal04 = goal04.proof_by_calculation()
        calc = proof_of_goal04.lhs_calc
        source = calc.parse_expr("INT x:[0,1]. (x - x ^ 2) ^ (1/2)")
        calc.perform_rule(rules.ApplyEquation(goal02.goal, source))
        calc.perform_rule(rules.Equation("2 * (1/2) + 1", "2"))
        calc.perform_rule(rules.Simplify())
        assert goal04.is_finished()

        goal05 = file.add_goal("factorial(1/2) = sqrt(pi) / 2")
        proof_of_goal05 = goal05.proof_by_rewrite_goal(begin = goal03)
        calc = proof_of_goal05.begin
        source = calc.parse_expr("INT x:[0,1]. (x - x ^ 2) ^ (1/2)")
        calc.perform_rule(rules.ApplyEquation(goal04.goal, source))
        calc.perform_rule(rules.SolveEquation("factorial(1/2)"))
        assert goal05.is_finished()

        goal06 = file.add_goal("(INT x:[0, oo]. exp(-x) * sqrt(x)) = sqrt(pi) / 2")
        proof_of_goal06 = goal06.proof_by_calculation()
        calc = proof_of_goal06.lhs_calc
        calc.perform_rule(rules.Equation("sqrt(x)", "x ^ (1/2)"))
        calc.perform_rule(rules.Equation("1/2", "3/2 - 1"))
        calc.perform_rule(rules.FoldDefinition("Gamma"))
        calc.perform_rule(rules.ApplyIdentity("Gamma(3 / 2)", "factorial(1 / 2)"))
        source = calc.parse_expr("factorial(1 / 2)")
        calc.perform_rule(rules.ApplyEquation(goal05.goal, source))
        assert goal06.is_finished()

        goal07 = file.add_goal("(INT x:[0, 1]. sqrt(-log(x))) = sqrt(pi) / 2")
        proof_of_goal07 = goal07.proof_by_calculation()
        calc = proof_of_goal07.lhs_calc
        calc.perform_rule(rules.SubstitutionInverse("y", "x", "exp(-y)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("sqrt(y) * exp(-y)", "exp(-y) * sqrt(y)"))
        source = calc.parse_expr("INT y:[0,oo]. exp(-y) * sqrt(y)")
        calc.perform_rule(rules.ApplyEquation(goal06.goal, source))
        assert goal07.is_finished()

        goal08 = file.add_goal("(INT x:[0, oo]. exp(-x) / sqrt(x)) = (INT x:[-oo, oo]. exp(-(x^2)))")
        proof_of_goal08 = goal08.proof_by_calculation()
        calc = proof_of_goal08.lhs_calc
        calc.perform_rule(rules.Substitution(var_name="y", var_subst="sqrt(x)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("2 * (INT y:[0,oo]. exp(-(y ^ 2)))",
                                         "(INT y:[0,oo]. exp(-(y ^ 2))) + (INT y:[0,oo]. exp(-(y ^ 2)))"))
        calc.perform_rule(rules.OnLocation(rules.Substitution(var_name="z", var_subst="-y"), "0"))
        calc = proof_of_goal08.rhs_calc
        calc.perform_rule(rules.SplitRegion("0"))
        assert goal08.is_finished()

        goal09 = file.add_goal("(INT x:[0, oo]. exp(-x) / sqrt(x)) = sqrt(pi)")
        proof_of_goal09 = goal09.proof_by_rewrite_goal(begin = goal08)
        calc = proof_of_goal09.begin
        calc.perform_rule(rules.Equation("x ^ 2", "1 * x ^ 2"))
        calc.perform_rule(rules.OnLocation(rules.DefiniteIntegralIdentity(), "1"))
        assert goal09.is_finished()

        goal10 = file.add_goal("Gamma(1/2) = sqrt(pi)")
        proof_of_goal10 = goal10.proof_by_calculation()
        calc = proof_of_goal10.lhs_calc
        calc.perform_rule(rules.ExpandDefinition("Gamma"))
        source = calc.parse_expr("INT x:[0,oo]. exp(-x) / sqrt(x)")
        calc.perform_rule(rules.ApplyEquation(goal09.goal, source))
        assert goal10.is_finished()

        goal11 = file.add_goal("factorial(-1/2) = sqrt(pi)")
        proof_of_goal11 = goal11.proof_by_rewrite_goal(begin = goal10)
        calc = proof_of_goal11.begin
        calc.perform_rule(rules.ApplyIdentity("Gamma(1/2)", "factorial(-1/2)"))
        assert goal11.is_finished()

        goal12 = file.add_goal("(INT x:[0, pi / 2]. sqrt(sin(x))) = Gamma(3/4) * Gamma(1/2) / (2 * Gamma(5/4))")
        proof_of_goal12 = goal12.proof_by_calculation()
        calc = proof_of_goal12.lhs_calc
        calc.perform_rule(rules.Substitution(var_name="u", var_subst="sin(x) ^ 2"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("1 / (u ^ (1/4) * sqrt(-u + 1))", "u ^ (3/4 - 1) * (1 - u) ^ (1/2 - 1)"))
        calc.perform_rule(rules.OnLocation(rules.FoldDefinition("B"), "1"))
        calc.perform_rule(rules.ApplyIdentity("B(3/4, 1/2)", "Gamma(3/4) * Gamma(1/2) / Gamma(5/4)"))
        calc.perform_rule(rules.Equation("1/2 * (Gamma(3/4) * Gamma(1/2) / Gamma(5/4))",
                                         "Gamma(3/4) * Gamma(1/2) / (2 * Gamma(5/4))"))
        assert goal12.is_finished()
        goal13 = file.add_goal("(INT x:[0, pi / 2]. sqrt(cos(x))) = Gamma(3/4) * Gamma(1/2) / (2 * Gamma(5/4))")
        proof_of_goal13 = goal13.proof_by_calculation()
        calc = proof_of_goal13.lhs_calc
        calc.perform_rule(rules.Substitution(var_name="u", var_subst="cos(x) ^ 2"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("1 / (u ^ (1/4) * sqrt(-u + 1))", "u ^ (3/4 - 1) * (1 - u) ^ (1/2 - 1)"))
        calc.perform_rule(rules.OnLocation(rules.FoldDefinition("B"), "1"))
        calc.perform_rule(rules.ApplyIdentity("B(3/4, 1/2)", "Gamma(3/4) * Gamma(1/2) / Gamma(5/4)"))
        calc.perform_rule(rules.Equation("1/2 * (Gamma(3/4) * Gamma(1/2) / Gamma(5/4))",
                                         "Gamma(3/4) * Gamma(1/2) / (2 * Gamma(5/4))"))
        assert goal13.is_finished()
        goal14 = file.add_goal("(INT x:[0, pi / 2]. 1 / sqrt(sin(x) * cos(x))) = 1/2 * B(1/4, 1/4)")
        proof_of_goal14 = goal14.proof_by_calculation()
        calc = proof_of_goal14.rhs_calc
        source = calc.parse_expr("B(1/4,1/4)")
        calc.perform_rule(rules.ApplyEquation(goal01.goal, source))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("sqrt(cos(x)) * sqrt(sin(x))", "sqrt(sin(x) * cos(x))"))
        assert goal14.is_finished()
        goal15 = file.add_goal("(INT x:[0, pi / 2]. 1 / sqrt(sin(x) * cos(x))) = Gamma(1/4) ^ 2 / (2 * sqrt(pi))")
        proof_of_goal15 = goal15.proof_by_rewrite_goal(begin = goal14)
        calc = proof_of_goal15.begin
        calc.perform_rule(rules.OnLocation(rules.ApplyIdentity("B(1/4, 1/4)", "Gamma(1/4) * Gamma(1/4) / Gamma(1/2)"), "1"))
        source = calc.parse_expr("Gamma(1/2)")
        calc.perform_rule(rules.ApplyEquation(goal10.goal, source))
        assert goal15.is_finished()
        goal16 = file.add_goal("(INT x:[0, pi / 2]. 1 / sqrt(sin(x))) = Gamma(1/4) ^ 2 / (2 * sqrt(2 * pi))")
        proof_of_goal16 = goal16.proof_by_rewrite_goal(begin = goal15)
        calc = proof_of_goal16.begin
        calc.perform_rule(rules.OnLocation(rules.Substitution(var_name="x", var_subst="2 * x"), "0"))
        calc.perform_rule(rules.Equation("2 * sqrt(cos(x / 2)) * sqrt(sin(x / 2))", "2 * sqrt(sin(1/2 * x) * cos(1/2 * x))"))
        calc.perform_rule(rules.ApplyIdentity("sin(1/2 * x) * cos(1/2 * x)", "(1 / 2) * (sin(1 / 2 * x + 1 / 2 * x) + sin(1 / 2 * x - 1 / 2 * x))"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.OnLocation(rules.SplitRegion("pi / 2"), "0.1"))
        calc.perform_rule(rules.OnLocation(rules.Substitution(var_name="x", var_subst="pi - x"), "0.1.1"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.SolveEquation("(INT x:[0, pi / 2]. 1 / sqrt(sin(x)))"))
        assert goal16.is_finished()
        goal17 = file.add_goal("(INT x:[0, pi / 2]. 1 / sqrt(cos(x))) = Gamma(1/4) ^ 2 / (2 * sqrt(2 * pi))")
        proof_of_goal17 = goal17.proof_by_rewrite_goal(begin = goal16)
        calc = proof_of_goal17.begin
        calc.perform_rule(rules.Substitution(var_name="x", var_subst="pi / 2 - x"))
        calc.perform_rule(rules.Simplify())
        assert goal17.is_finished()
        goal18 = file.add_goal("B(m, n) = (INT x:[0, oo]. x ^ (m - 1) / (x + 1) ^ (m + n))", conds=["m > 0", "n > 0"])
        proof_of_goal18 = goal18.proof_by_calculation()
        calc = proof_of_goal18.lhs_calc
        calc.perform_rule(rules.ExpandDefinition("B"))
        calc.perform_rule(rules.Substitution(var_name="y", var_subst="1 / (1 - x) - 1"))
        calc.perform_rule(rules.Equation("-(1 / (y + 1)) + 1", "y / (y + 1)"))
        calc.perform_rule(rules.ApplyIdentity("(1 / (y + 1)) ^ (n - 1)", "1 / (y + 1) ^ (n - 1)"))
        calc.perform_rule(rules.ApplyIdentity("(y / (y + 1)) ^ (m - 1)", "y ^ (m - 1) / (y + 1) ^ (m - 1)"))
        calc.perform_rule(rules.Equation("1 / (y + 1) ^ 2 * (1 / (y + 1) ^ (n - 1)) * (y ^ (m - 1) / (y + 1) ^ (m - 1))",
                                         "y ^ (m - 1) / (y + 1) ^ (m + n)"))
        assert goal18.is_finished()
        goal19 = file.add_goal("(INT x:[0, oo]. x ^ (m - 1) / (x + 1)) = B(m, 1 - m)", conds=["m > 0", "m < 1"])
        proof_of_goal19 = goal19.proof_by_calculation()
        calc = proof_of_goal19.rhs_calc
        source = calc.parse_expr("B(m,1 - m)")
        calc.perform_rule(rules.ApplyEquation(goal18.goal, source))
        calc.perform_rule(rules.Simplify())
        assert goal19.is_finished()
        goal20 = file.add_goal("B(m, 1 - m) = Gamma(m) * Gamma(1 - m)", conds=["m > 0", "m < 1"])
        proof_of_goal20 = goal20.proof_by_calculation()
        calc = proof_of_goal20.lhs_calc
        calc.perform_rule(rules.ApplyIdentity("B(m, 1 - m)", "Gamma(m) * Gamma(1 - m) / Gamma(1)"))
        calc.perform_rule(rules.ApplyIdentity("Gamma(1)", "1"))
        calc.perform_rule(rules.Equation("Gamma(m) * Gamma(1 - m) / 1", "Gamma(m) * Gamma(1 - m)"))
        assert goal20.is_finished()
        goal21 = file.add_goal("Gamma(m) * Gamma(1 - m) = pi / sin(m * pi)", conds=["m > 0", "m < 1"])
        proof_of_goal21 = goal21.proof_by_rewrite_goal(begin = goal19)
        calc = proof_of_goal21.begin
        source = calc.parse_expr("B(m,1 - m)")
        calc.perform_rule(rules.ApplyEquation(goal20.goal, source))
        calc.perform_rule(rules.OnLocation(rules.DefiniteIntegralIdentity(), "0"))
        calc.perform_rule(rules.SolveEquation("Gamma(m) * Gamma(1 - m)"))
        assert goal21.is_finished()
        goal22 = file.add_goal("factorial(z) * factorial(z) / factorial(2 * z + 1) = (INT x:[0, 1]. x ^ z * (1 - x) ^ z)",
                               conds=["z > -1"])
        proof_of_goal22 = goal22.proof_by_calculation()
        calc = proof_of_goal22.rhs_calc
        calc.perform_rule(rules.OnLocation(rules.Equation("z", "(z + 1) - 1"), "0.0"))
        calc.perform_rule(rules.OnLocation(rules.Equation("z", "(z + 1) - 1"), "0.1"))
        calc.perform_rule(rules.FoldDefinition("B"))
        calc.perform_rule(rules.ApplyIdentity("B(z + 1, z + 1)", "Gamma(z + 1) * Gamma(z + 1) / Gamma(2 * z + 2)"))
        calc.perform_rule(rules.ApplyIdentity("Gamma(z + 1)", "factorial(z)"))
        calc.perform_rule(rules.ApplyIdentity("Gamma(z + 1)", "factorial(z)"))
        calc.perform_rule(rules.ApplyIdentity("Gamma(2 * z + 2)", "factorial(2 * z + 1)"))
        assert goal22.is_finished()
        goal23 = file.add_goal("factorial(z) * factorial(z + 1/2) = 2 ^ (-2 * z - 1) * sqrt(pi) * factorial(2 * z + 1)",
                               conds=["z > -1"])
        proof_of_goal23 = goal23.proof_by_rewrite_goal(begin = goal22)
        calc = proof_of_goal23.begin
        calc.perform_rule(rules.OnLocation(rules.Substitution(var_name="s", var_subst="2 * x - 1"), "1"))
        calc.perform_rule(rules.Equation("-((s + 1) / 2) + 1", "1/2 * (1 - s)"))
        calc.perform_rule(rules.Equation("(s + 1) / 2", "1/2 * (1 + s)"))
        calc.perform_rule(rules.ApplyIdentity("(1/2 * (1 - s)) ^ z", "(1/2) ^ z * (1 - s) ^ z"))
        calc.perform_rule(rules.ApplyIdentity("(1/2 * (1 + s)) ^ z", "(1/2) ^ z * (1 + s) ^ z"))
        calc.perform_rule(rules.OnLocation(rules.Simplify(), "1"))
        calc.perform_rule(rules.ApplyIdentity("(s + 1) ^ z * (-s + 1) ^ z", "((s + 1) * (-s + 1)) ^ z"))
        calc.perform_rule(rules.Equation("(s + 1) * (-s + 1)", "(1 - s ^ 2)"))
        calc.perform_rule(rules.SplitRegion("0"))
        calc.perform_rule(rules.OnLocation(rules.Substitution(var_name="s", var_subst="-s"), "1.1.0"))
        calc.perform_rule(rules.OnLocation(rules.Simplify(), "1"))
        calc.perform_rule(rules.Substitution(var_name="u", var_subst="s ^ 2"))
        calc.perform_rule(rules.OnLocation(rules.Simplify(), "1"))
        calc.perform_rule(rules.Equation("(-u + 1) ^ z / sqrt(u)", "u ^ (1/2 - 1) * (1 - u) ^ (z + 1 - 1)"))
        calc.perform_rule(rules.OnLocation(rules.FoldDefinition("B"), "1.1"))
        calc.perform_rule(rules.ApplyIdentity("B(1/2,z + 1)", "Gamma(1/2) * Gamma(z + 1) / Gamma(z + 3/2)"))
        source = calc.parse_expr("Gamma(1/2)")
        calc.perform_rule(rules.ApplyEquation(goal10.goal, source))
        calc.perform_rule(rules.ApplyIdentity("Gamma(z + 1)", "factorial(z)"))
        calc.perform_rule(rules.ApplyIdentity("Gamma(z + 3/2)", "factorial(z + 1/2)"))
        calc.perform_rule(rules.SolveEquation("factorial(2 * z + 1)"))
        calc.perform_rule(rules.Equation("2 * (1/2) ^ -(2 * z) * factorial(z) * factorial(z + 1/2) / sqrt(pi)",
                                         "2 / sqrt(pi) * (1/2) ^ (-2 * z) * (factorial(z) * factorial(z + 1/2))"))
        calc.perform_rule(rules.SolveEquation("factorial(z) * factorial(z + 1/2)"))
        calc.perform_rule(rules.Equation("(1/2) ^ (2 * z) * sqrt(pi) * factorial(2 * z + 1) / 2", "sqrt(pi) * 2 ^ (-1) * (2 ^ (-1)) ^ (2 * z) * factorial(2 * z + 1)"))
        calc.perform_rule(rules.ApplyIdentity("(2 ^ (-1)) ^ (2 * z)", "2 ^ (-2 * z)"))
        calc.perform_rule(rules.Equation("sqrt(pi) * 2 ^ (-1) * 2 ^ (-2 * z)", "2 ^ (-2 * z - 1) * sqrt(pi)"))
        assert goal23.is_finished()
        # self.checkAndOutput(file)

    # def testChapter1Practice0101(self):
    #     # Reference:
    #     # Inside interesting integrals, C1.1
    #     file = compstate.CompFile("interesting", "chapter1_practice01_01")
    #     goal = file.add_goal("(INT x:[0,8]. 1/(x-2)) = log(3)")
    #     proof = goal.proof_by_calculation()
    #     calc = proof.lhs_calc
    #     calc.perform_rule(rules.Substitution("u", "x-2"))
    #     calc.perform_rule(rules.SplitRegion("0"))
    #     calc.perform_rule(rules.DefiniteIntegralIdentity())
    #     calc.perform_rule(rules.Simplify())
    #     # self.checkAndOutput(file)

    # def testChapter1Practice0102(self):
    #     # Reference:
    #     # Inside interesting integrals, C1.1
    #     file = compstate.CompFile("interesting", "chapter1_practice01_02")
    #     goal = file.add_goal("(INT x:[0,3]. 1/(x-1)^(2/3)) = 3 * (1+2^(1/3))")
    #     proof = goal.proof_by_calculation()
    #     calc = proof.lhs_calc
    #     calc.perform_rule(rules.Substitution("u", "x-1"))
    #     calc.perform_rule(rules.SplitRegion("0"))
    #     calc.perform_rule(rules.DefiniteIntegralIdentity())
    #     calc.perform_rule(rules.Simplify())
    #     calc = proof.rhs_calc
    #     calc.perform_rule(rules.ExpandPolynomial())
    #     # self.checkAndOutput(file)

    def testChapter1Practice0104(self):
        # Reference:
        # Inside interesting integrals, C1.5
        file = compstate.CompFile("interesting", "chapter1_practice0104")

        goal = file.add_goal("(INT x:[0,pi/3]. 1/cos(x))=log(2+sqrt(3))")
        proof = goal.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Equation("1/cos(x)", "cos(x)/(cos(x))^2"))
        calc.perform_rule(rules.ApplyIdentity("(cos(x))^2", "1-(sin(x))^2"))
        calc.perform_rule(rules.Substitution(var_name="u", var_subst="sin(x)"))
        calc.perform_rule(rules.Equation("1/(-(u^2)+1)", "1/2 * (1/(1-u)+1/(1+u))"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Substitution(var_name="t", var_subst="-u+1"))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("-(1/2 * log(-(sqrt(3) / 2) + 1)) + 1/2 * log(sqrt(3) / 2 + 1)",
                                         "1/2 * (log(sqrt(3) / 2 + 1) - log(-(sqrt(3) / 2) + 1))"))
        calc.perform_rule(rules.ApplyIdentity("log(sqrt(3) / 2 + 1) - log(-(sqrt(3) / 2) + 1)",
                                              "log((sqrt(3) / 2 + 1)/(-(sqrt(3) / 2) + 1))"))
        calc.perform_rule(rules.Equation("(sqrt(3) / 2 + 1)/(-(sqrt(3) / 2) + 1)", "(2+sqrt(3))^2"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("sqrt(3)+2", "2+sqrt(3)"))

        # self.checkAndOutput(file)

    def testChapter2Practice01(self):
        # Reference:
        # Inside interesting integrals, C2.1
        file = compstate.CompFile("interesting", "chapter2_practice01")

        goal = file.add_goal("(INT x:[0,4]. log(x)/sqrt(4*x-x^2)) = 0")

        goal01 = goal.add_subgoal("1", "(INT y:[0,1]. 1/ (sqrt(y) * sqrt(1-y))) = pi")
        proof = goal01.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.SubstitutionInverse("x", "y", "sin(x)^2"))
        calc.perform_rule(rules.OnCount(rules.ApplyIdentity("sin(x)^2", "1-cos(x)^2"), 2))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        goal02 = goal.add_subgoal("2", "(INT y:[0,1]. log(y)/ (sqrt(y) * sqrt(1-y))) = -(2*pi*log(2))")
        proof = goal02.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.SubstitutionInverse("x", "y", "sin(x)^2"))
        calc.perform_rule(rules.Equation("log(sin(x)^2)", "2*log(sin(x))"))
        calc.perform_rule(rules.OnCount(rules.ApplyIdentity("sin(x)^2", "1-cos(x)^2"), 2))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("sin(x)", "1*sin(x)"))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        sub_goal = goal.add_subgoal("3", "4 * x - x ^ 2 >= 0", conds=["x > 0", "x < 4"])
        proof = sub_goal.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Equation("4 * x - x ^ 2", "x*(4 - x)"))

        sub_goal2 = goal.add_subgoal("4", "sqrt(4 * x - x ^ 2) != 0", conds=["x > 0", "x < 4"])
        proof = sub_goal2.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Equation("4 * x - x ^ 2", "x*(4 - x)"))

        proof = goal.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Substitution("y", "x/4"))
        calc.perform_rule(rules.Equation("log(4*y)", "log(4)+log(y)"))
        calc.perform_rule(rules.Equation("sqrt(-(16 * y ^ 2) + 16 * y)", "4 * sqrt(-y ^ 2 + y)"))
        calc.perform_rule(rules.Equation("sqrt(-y^2 + y)", "sqrt(y) * sqrt(1-y)"))
        calc.perform_rule(rules.ExpandPolynomial())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("-y+1", "1-y"))
        calc.perform_rule(rules.Equation("-y+1", "1-y"))
        source = calc.parse_expr("INT y:[0,1]. 1/ (sqrt(y) * sqrt(1-y))")
        calc.perform_rule(rules.ApplyEquation("1", source))
        source = calc.parse_expr("INT y:[0,1]. log(y)/ (sqrt(y) * sqrt(1-y))")
        calc.perform_rule(rules.ApplyEquation("2", source))
        calc.perform_rule(rules.Simplify())
        assert goal.is_finished()

    def testChapter2Practice02(self):
        # Reference:
        # Inside interesting integrals, C2.2
        file = compstate.CompFile("interesting", "chapter2_practice02")
        goal01 = file.add_goal("(INT x:[0,1]. (x - 2) / (x ^ 2 - x + 1)) = -pi/sqrt(3)")
        proof = goal01.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Equation("(x ^ 2 - x + 1)", "(x-1/2)^2 + 3/4"))
        calc.perform_rule(rules.Substitution("u", "x-1/2"))
        calc.perform_rule(rules.OnLocation(rules.ExpandPolynomial(), "0"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("INT u:[-1/2,1/2]. u / (u ^ 2 + 3/4)", "0"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Substitution("x", "2*u/sqrt(3)"))
        calc.perform_rule(rules.Equation("3 * x ^ 2 / 2 + 3/2", "3/2*(x^2+1)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        calc = proof.rhs_calc
        calc.perform_rule(rules.Simplify())

        goal02 = file.add_goal("(INT x:[0,1]. 1/(x^3+1)) = 1/3 *log(2) + pi / (3*sqrt(3))")
        proof = goal02.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Equation("1/(x^3 + 1)", "1/((x+1)*(x^2-x+1))"))
        calc.perform_rule(rules.Equation("1/((x+1)*(x^2-x+1))", "1/3*(1/(x+1)-(x-2)/(x^2-x+1))"))
        calc.perform_rule(rules.OnLocation(rules.ExpandPolynomial(), "0"))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        calc = proof.rhs_calc
        calc.perform_rule(rules.Simplify())
        # self.checkAndOutput(file)


    def testChapter2Practice03(self):
        # Reference:
        # Inside interesting integrals, C2.3
        file = compstate.CompFile("interesting", "chapter2_practice03")

        goal = file.add_goal("(INT x:[0, oo]. 1/(x^4+1)^(m+1)) = (4*m-1)/(4*m) * (INT x:[0, oo]. 1/(x^4+1)^m)",
                             conds=["m>=1", "isInt(m)"])

        goal01 = goal.add_subgoal("1", "(INT x:[0, oo]. 1/(x^4+1)^m) = 4*m*((INT x:[0, oo]. 1 / (x^4+1)^m)-(INT x:[0, oo]. 1 / (x^4+1)^(m+1)))",
                                  conds=["m>=1", "isInt(m)"])
        proof = goal01.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.IntegrationByParts(u="1/(x^4+1)^m", v="x"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("x ^ 4 * (x ^ 4 + 1) ^ (-m - 1)", "(x^4+1)/(x^4+1)^(m+1)-1/(x^4+1)^(m+1)"))
        calc.perform_rule(rules.Equation("(INT x:[0,oo]. (x ^ 4 + 1) / (x ^ 4 + 1) ^ (m + 1) - 1 / (x ^ 4 + 1) ^ (m + 1))",
                                         "(INT x:[0, oo]. 1 / (x^4+1)^m)-(INT x:[0, oo]. 1 / (x^4+1)^(m+1))"))
        assert goal01.is_finished()

        proof = goal.proof_by_rewrite_goal(begin="1")
        calc = proof.begin
        calc.perform_rule(rules.SolveEquation("(INT x:[0, oo]. 1/(x^4+1)^(m+1))"))
        calc.perform_rule(rules.Equation("-(1 / (4 * m) * (INT x:[0,oo]. (x ^ 4 + 1) ^ -m)) + (INT x:[0,oo]. (x ^ 4 + 1) ^ -m)",
                                         "(4*m-1)/(4*m) * (INT x:[0, oo]. 1/(x^4+1)^m)"))
        assert goal.is_finished()

    def testChapter2Practice05(self):
        # Reference:
        # Inside interesting integrals, C2.5
        file = compstate.CompFile("interesting", "chapter2_practice05")
        goal = file.add_goal("(INT x:[0,oo]. log(1+x) / (x*sqrt(x))) = 2 * pi")
        proof = goal.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.IntegrationByParts("log(1+x)", "-2/sqrt(x)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Substitution("t", "sqrt(x)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        # self.checkAndOutput(file)

    def testChapter4Practice01(self):
        # Reference:
        # Inside interesting integrals, C4.1
        file = compstate.CompFile("interesting", "chapter4_practice01")
        goal01 = file.add_goal("B(2, n+1) = INT u:[0,1]. u * (1-u)^n", conds=['n>-1'])
        proof = goal01.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.ExpandDefinition("B"))
        calc = proof.rhs_calc
        calc.perform_rule(rules.Simplify())
        assert goal01.is_finished()

        goal02 = file.add_goal("(INT x:[0,1]. (1-sqrt(x))^n) = 2 / ((n+1)*(n+2))", conds=["n>-1"])
        proof = goal02.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Substitution("u", "sqrt(x)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("-u+1","1-u"))
        s = calc.parse_expr("INT u:[0,1]. u * (1 - u) ^ n")
        l = calc.parse_expr("B(2, n+1) = INT u:[0,1]. u * (1-u)^n")
        calc.perform_rule(rules.ApplyEquation(l, s))
        s = calc.parse_expr("B(2,n+1)")
        t = calc.parse_expr("Gamma(2)*Gamma(n+1) / Gamma(n+3)")
        calc.perform_rule(rules.ApplyIdentity(s, t))
        s = calc.parse_expr("Gamma(2)")
        t = calc.parse_expr("factorial(1)")
        calc.perform_rule(rules.ApplyIdentity(s, t))
        s = calc.parse_expr("Gamma(n+1)")
        t = calc.parse_expr("factorial(n)")
        calc.perform_rule(rules.ApplyIdentity(s, t))
        s = calc.parse_expr("Gamma(n+3)")
        t = calc.parse_expr("factorial(n+2)")
        calc.perform_rule(rules.ApplyIdentity(s, t))
        calc.perform_rule(rules.Simplify())
        s = calc.parse_expr("factorial(n+2)")
        t = calc.parse_expr("(n+2)*factorial(n+1)")
        calc.perform_rule(rules.ApplyIdentity(s, t))
        s = calc.parse_expr("factorial(n+1)")
        t = calc.parse_expr("(n+1)*factorial(n)")
        calc.perform_rule(rules.ApplyIdentity(s, t))
        calc.perform_rule(rules.Simplify())
        assert goal02.is_finished()

        goal03 = file.add_goal("(INT x:[0,1]. (1-sqrt(x))^9) = 1 / 55")
        proof = goal03.proof_by_rewrite_goal(begin = goal02)
        calc = proof.begin
        calc.perform_rule(rules.VarSubsOfEquation([{'var': 'n', 'expr': "9"}]))
        calc.perform_rule(rules.Simplify())
        assert goal03.is_finished()

        self.checkAndOutput(file)

    # The condition of goal02 can not be weakened until complex integration is supported.
    def testChapter4Practice02(self):
        # Reference:
        # Inside interesting integrals, C4.2
        file = compstate.CompFile("interesting", "chapter4_practice02")

        goal01 = file.add_goal("Gamma(n+1) = INT t:[0,oo]. t^n * exp(-t)")
        proof = goal01.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.ExpandDefinition("Gamma"))

        goal02 = file.add_goal("(INT x:[0,1]. x^m * log(x)^n) = (-1)^n * factorial(n) / (m+1)^(n+1)", conds=["m > -1", "n >= 0", "isInt(n)"])
        proof = goal02.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.SubstitutionInverse("u", "x", "exp(-u)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("(-u) ^ n * exp(-u) * exp(-u) ^ m", "(-u) ^ n * (exp(-u) * exp(-u) ^ m)"))
        calc.perform_rule(rules.Equation("exp(-u) * exp(-u) ^ m", "exp(-u)^1 * exp(-u)^m"))
        calc.perform_rule(rules.OnLocation(rules.ApplyIdentity("exp(-u)^1 * exp(-u)^m","exp(-u)^(m+1)"), "0.1"))
        calc.perform_rule(rules.Equation("-u", "-1*u"))
        calc.perform_rule(rules.OnLocation(rules.ApplyIdentity("(-1*u)^n", "(-1)^n * u^n"), "0.0"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Substitution("t", "(m+1)*u"))
        calc.perform_rule(rules.ApplyIdentity("(t / (m + 1)) ^ n", "t ^ n / (m+1)^n"))
        calc.perform_rule(rules.ApplyIdentity("exp(-(t / (m + 1)))^(m + 1)", "exp(-(t/(m+1)) * (m+1))"))
        calc.perform_rule(rules.Simplify())
        l = parser.parse_expr("Gamma(n+1) = INT t:[0,oo]. t ^ n * exp(-t)")
        s = calc.parse_expr("(INT t:[0,oo]. t ^ n * exp(-t))")
        calc.perform_rule(rules.ApplyEquation(l, s))
        calc.perform_rule(rules.ApplyIdentity("Gamma(n+1)", "factorial(n)"))
        calc.perform_rule(rules.Simplify())
        calc = proof.rhs_calc
        calc.perform_rule(rules.Simplify())
        # self.checkAndOutput(file)

    def testChapter4Practice03(self):
        # Reference:
        # Inside interesting integrals, C4.3
        file = compstate.CompFile("interesting", "chapter4_practice03")
        goal01 = file.add_goal("B(a+1,b+2) = (INT x:[0,1]. x ^ a * (-x + 1) ^ (b + 1))", conds=['a>-1', 'b>-2'])
        proof = goal01.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.ExpandDefinition("B"))
        assert goal01.is_finished()
        goal02 = file.add_goal("(INT x:[0,1]. x^a * (INT y:[0 ,1-x]. y^b)) = factorial(a) * factorial(b) / factorial(a+b+2)", conds=["b >= 0", 'a>=0'])
        proof = goal02.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.OnLocation(rules.DefiniteIntegralIdentity(), "0.1"))
        calc.perform_rule(rules.Simplify())
        l = parser.parse_expr("B(a+1,b+2) = (INT x:[0,1]. x ^ a * (-x + 1) ^ (b + 1))")
        s = calc.parse_expr("(INT x:[0,1]. x ^ a * (-x + 1) ^ (b + 1))")
        calc.perform_rule(rules.ApplyEquation(l, s))
        calc.perform_rule(rules.ApplyIdentity("B(a + 1,b + 2)", "Gamma(a+1) * Gamma(b+2)/Gamma(a+b+3)"))
        calc.perform_rule(rules.ApplyIdentity("Gamma(a+1)", "factorial(a)"))
        calc.perform_rule(rules.ApplyIdentity("Gamma(b+2)", "factorial(b+1)"))
        calc.perform_rule(rules.ApplyIdentity("Gamma(a+b+3)", "factorial(a+b+2)"))
        calc.perform_rule(rules.ApplyIdentity("factorial(b+1)","(b+1)*factorial(b)"))
        calc.perform_rule(rules.Simplify())
        assert goal02.is_finished()
        # self.checkAndOutput(file)

    def testExpSinh(self):
        file = compstate.CompFile("interesting", "exp_sinh")

        sub_goal = file.add_goal("9 - 10 * p ^ 2 + p ^ 4 != 0", conds=["p > 3"])
        proof = sub_goal.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Equation("9 - 10 * p ^ 2 + p ^ 4", "(p ^ 2 - 5) ^ 2 - 16"))
        self.assertTrue(proof.is_finished())

        goal = file.add_goal("(INT t:[0, oo]. exp(-p*t) * sinh(t)^3) = 6 / (9 - 10 * p ^ 2 + p ^ 4)", conds=["p > 3"])
        proof = goal.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.OnSubterm(rules.ExpandDefinition("sinh")))
        calc.perform_rule(rules.ExpandPolynomial())
        calc.perform_rule(rules.Equation("-(p * t) - 3 * t", "(-p-3)*t"))
        calc.perform_rule(rules.Equation("-(p * t) + 3 * t", "(3-p)*t"))
        calc.perform_rule(rules.Equation("-(p * t) - t", "(-p-1)*t"))
        calc.perform_rule(rules.Equation("-(p * t) + t", "(1-p)*t"))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("3 / (-(8 * p) + 8) - 1 / (-(8 * p) + 24) - 3 / (-(8 * p) - 8) + 1 / (-(8 * p) - 24)",
                                         "6 / (9 - 10 * p^2 + p^4)"))
        # self.checkAndOutput(file)

    def testZetaFunction(self):
        # Reference:
        # Inside interesting integrals, section 5.3
        file = compstate.CompFile("interesting", "zeta_function")
        file.add_definition("zeta(s) = SUM(k, 0, oo, 1/(k+1)^s)", conds=["type(s,1)"])
        s = "(INT y:[0,1]. (INT x:[0,1]. x^a * y^a / (1-x*y))) = SUM(n, 0, oo, 1/(n+1+a)^2)"
        goal01 = file.add_goal(s, conds=["a>-1"])
        proof = goal01.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("x ^ a / (-(x * y) + 1)")
        s2 = calc.parse_expr("x ^ a * (1-x*y)^(-1)")
        calc.perform_rule(rules.Equation(s1, s2))
        calc.perform_rule(rules.SeriesExpansionIdentity(old_expr="(1 - x * y) ^ (-1)", index_var='k'))
        s1 = calc.parse_expr("x ^ a * SUM(k, 0, oo, (x * y) ^ k)")
        s2 = calc.parse_expr("SUM(k, 0, oo, x ^ a * (x * y) ^ k)")
        calc.perform_rule(rules.Equation(s1, s2))
        calc.perform_rule(rules.OnLocation(rules.IntSumExchange(), "0.1"))
        s1 = calc.parse_expr(" (x * y) ^ k")
        s2 = calc.parse_expr(" x ^ k * y ^ k")
        calc.perform_rule(rules.OnLocation(rules.ApplyIdentity(s1, s2), "0.1.0.0.1"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.OnLocation(rules.DefiniteIntegralIdentity(), "0"))
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("y ^ a * SUM(k, 0, oo, y ^ k / (a + k + 1))")
        s2 = calc.parse_expr("SUM(k, 0, oo, y ^ a * (y ^ k / (a + k + 1)))")
        calc.perform_rule(rules.Equation(s1, s2))
        calc.perform_rule(rules.IntSumExchange())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        calc = proof.rhs_calc
        calc.perform_rule(rules.Simplify())
        assert goal01.is_finished()

        #application
        s = "(INT y:[0,1]. (INT x:[0,1]. 1 / (1-x*y))) = zeta(2)"
        goal02 = file.add_goal(s)
        proof = goal02.proof_by_rewrite_goal(begin=goal01)
        calc = proof.begin
        calc.perform_rule(rules.VarSubsOfEquation([{'var':'a', 'expr':'0'}]))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("-(x * y) + 1", "1-x*y"))
        calc.perform_rule(rules.OnLocation(rules.FoldDefinition('zeta'), '1'))
        assert goal02.is_finished()
        s1 = "(INT y:[0,1]. (INT x:[0,1]. log(x*y)^(s-2)*(x^a * y^a) / (1-x*y)))"
        s2 = "(-1)^s * factorial(s-1) * SUM(n, 0, oo, 1/(n+a+1)^s)"
        goal03 = file.add_goal(s1 + "=" + s2, conds=["a > -1", "s >= 2", "isInt(s)"])
        proof = goal03.proof_by_induction('s', 2)
        proof_base = proof.base_case.proof_by_calculation()
        proof_induct = proof.induct_case.proof_by_rewrite_goal(begin = goal03)
        calc = proof_base.lhs_calc
        calc.perform_rule(rules.OnLocation(rules.Simplify(), "0.0"))
        calc.perform_rule(rules.Equation("-(x * y) + 1", "1-x*y"))
        s = calc.parse_expr("INT y:[0,1]. INT x:[0,1]. x ^ a * y ^ a / (1 - x * y)")
        calc.perform_rule(rules.ApplyEquation(goal01.goal, s))
        calc.perform_rule(rules.Simplify())
        calc = proof_induct.begin
        calc.perform_rule(rules.DerivEquation('a'))
        calc.perform_rule(rules.OnLocation(rules.DerivIntExchange(), "0"))
        calc.perform_rule(rules.OnLocation(rules.DerivIntExchange(), "0.0"))
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("x ^ a * y ^ a * log(x) + x ^ a * y ^ a * log(y)")
        s2 = calc.parse_expr("x ^ a * y ^ a * (log(x) + log(y))")
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = calc.parse_expr("-(x * y) + 1")
        s2 = calc.parse_expr("1 - x*y")
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = calc.parse_expr("log(x) + log(y)")
        s2 = calc.parse_expr("log(x*y)")
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = calc.parse_expr("log(x * y) ^ (s - 2) * (x ^ a * y ^ a * log(x * y))")
        s2 = calc.parse_expr("x ^ a * y ^ a * log(x * y) ^ (s - 1)")
        calc.perform_rule(rules.Equation(s1, s2))

        s1 = calc.parse_expr("-(s * (-1) ^ s * factorial(s - 1) * SUM(n, 0, oo, (a + n + 1) ^ (-s - 1)))")
        s2 = calc.parse_expr("(-1) ^ (s+1) * (s* factorial(s - 1)) * SUM(n, 0, oo, (a + n + 1) ^ (-s - 1))")
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = calc.parse_expr("s* factorial(s - 1)")
        s2 = calc.parse_expr("factorial(s)")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        assert goal03.is_finished()
        s1 = 'zeta(s)'
        s2 = "(-1)^s / factorial(s-1) * \
            (INT y:[0,1]. (INT x:[0,1]. log(x*y)^(s-2)/ (1-x*y)))"
        goal04 = file.add_goal(s1 + '=' + s2, conds=["s >= 2", "isInt(s)"])
        proof = goal04.proof_by_rewrite_goal(begin = goal03)
        calc = proof.begin
        calc.perform_rule(rules.VarSubsOfEquation([{'var':'a', 'expr':'0'}]))
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("(n + 1) ^ -s")
        s2 = calc.parse_expr("1/(n+1)^s")
        calc.perform_rule(rules.Equation(s1, s2))
        calc.perform_rule(rules.OnLocation(rules.FoldDefinition('zeta'), "1.1"))
        s1 = calc.parse_expr("-(x * y) + 1")
        s2 = calc.parse_expr("1 - x*y")
        calc.perform_rule(rules.Equation(s1, s2))
        s = calc.parse_expr("zeta(s)")
        calc.perform_rule(rules.SolveEquation(s))
        s1 = calc.parse_expr("-s")
        s2 = calc.parse_expr("-1*s")
        calc.perform_rule(rules.Equation(s1,s2))
        s1 = calc.parse_expr("(-1) ^ (-1 * s)")
        s2 = calc.parse_expr("((-1)^-1^s)")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        s1 = calc.parse_expr("(-1)^-1")
        s2 = calc.parse_expr("-1")
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = calc.parse_expr("-(x * y) + 1")
        s2 = calc.parse_expr("1 - x * y")
        calc.perform_rule(rules.Equation(s1, s2))

        s1 = "(INT x:[0, oo]. exp(-k*x) * x^(s-1))"
        s2 = "Gamma(s) / k^s"
        goal05 = file.add_goal(s1+'='+s2, conds=["k>0"])
        proof = goal05.proof_by_calculation()
        calc = proof.lhs_calc
        s1 = "u"
        s2 = calc.parse_expr("k*x")
        calc.perform_rule(rules.Substitution(s1, s2))
        s1 = calc.parse_expr("(u / k) ^ (s - 1)")
        s2 = calc.parse_expr("u ^ (s - 1) / k ^ (s - 1)")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("u ^ (s - 1) * exp(-u)")
        s2 = calc.parse_expr("exp(-u) * u ^ (s - 1)")
        calc.perform_rule(rules.Equation(s1, s2))
        calc.perform_rule(rules.OnLocation(rules.FoldDefinition("Gamma"), "1"))
        calc = proof.rhs_calc
        calc.perform_rule(rules.Simplify())

        s1 = "(INT x:[0,oo]. x^(s-1)/(exp(x) - 1))"
        s2 = "Gamma(s) * zeta(s)"
        goal06 = file.add_goal(s1+'='+s2)
        proof = goal06.proof_by_rewrite_goal(begin = goal05)
        calc = proof.begin
        calc.perform_rule(rules.SummationEquation('k', '1', 'oo'))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.OnLocation(rules.ChangeSummationIndex('0'), "1.1"))
        s1 = calc.parse_expr("(k + 1) ^ -s")
        s2 = calc.parse_expr("1/(k+1)^s")
        calc.perform_rule(rules.Equation(s1, s2))
        calc.perform_rule(rules.OnLocation(rules.FoldDefinition("zeta"),"1.1"))
        calc.perform_rule(rules.OnLocation(rules.IntSumExchange(), "0"))
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("exp(-(k * x))")
        s2 = calc.parse_expr("exp(-x*k)")
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = calc.parse_expr("exp(-x*k)")
        s2 = calc.parse_expr("exp(-x)^k")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        calc.perform_rule(rules.OnLocation(rules.SeriesEvaluationIdentity(), "0.0.1"))
        s1 = calc.parse_expr("x ^ (s - 1) * (exp(-x) ^ 1 / (1 - exp(-x)))")
        s2 = calc.parse_expr("x^(s-1) / (exp(x) - 1)")
        calc.perform_rule(rules.Equation(s1, s2))
        # assert goal03.is_finished()
        assert goal01.is_finished()
        assert goal02.is_finished()
        assert goal03.is_finished()
        assert goal04.is_finished()
        assert goal05.is_finished()
        assert goal06.is_finished()


        # self.checkAndOutput(file)

    def testCoxeterIntegral(self):
        # TODO: find some problems about condition inheritance
        file = compstate.CompFile("interesting", "coxeter")
        lemma = file.add_goal("cos(2*x) = 2*cos(x)^2 - 1", conds=["x>=0", "x<=pi/2"])
        proof = lemma.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.ApplyIdentity("cos(2*x)", "2*cos(x)^2 - 1"))

        goal01 = file.add_goal("acos(a) = 2*acos(sqrt((1+a) / 2))", conds=["abs(a) <= 1"])
        proof = goal01.proof_by_rewrite_goal(begin = lemma)
        calc = proof.begin
        calc.perform_rule(rules.VarSubsOfEquation([{"var":"x", "expr":"acos(u)"}]))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.FunEquation("acos"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.VarSubsOfEquation([{"var": "u", "expr": "sqrt((1+a)/2)"}]))
        calc.perform_rule(rules.Simplify())
        assert goal01.is_finished()
        # self.checkAndOutput(file)

    def testPowerfulElementaryIntegral(self):
        # Reference: Impossible, Integrals, Sums, and Series
        # section 1.1
        file = compstate.CompFile("impossible", "powerful_elementry_integral")

        s1 = "(INT x:[0, 1]. 1/ ((1 + y*x)*sqrt(1-x^2)))"
        s2 = "2 / sqrt(-(y ^ 2) + 1) * (atan((y + 1) / sqrt(-(y ^ 2) + 1)) - atan(y / sqrt(-(y ^ 2) + 1))) "
        goal01 = file.add_goal(s1+"="+s2, conds=["abs(y) < 1"])
        proof = goal01.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.SubstitutionInverse("t", "x", "sin(t)"))
        calc.perform_rule(rules.OnLocation(rules.Equation("1", "sin(t)^2 + cos(t)^2"), "0.0.1.1"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.OnLocation(rules.Equation("1", "sin(t/2)^2 + cos(t/2)^2"), "0.1"))
        calc.perform_rule(rules.OnLocation(rules.Equation("t", "2*(t/2)"), "0.1.0"))
        calc.perform_rule(rules.OnLocation(rules.ApplyIdentity("sin(2*(t/2))","2*sin(t/2)*cos(t/2)"), "0.1.0.1"))
        s1 = "1 / (y * (2 * sin(t / 2) * cos(t / 2)) + (sin(t / 2) ^ 2 + cos(t / 2) ^ 2))"
        s2 = "(1/cos(t/2)^2) / ((y * (2 * sin(t / 2) * cos(t / 2)) + (sin(t / 2) ^ 2 + cos(t / 2) ^ 2))/cos(t/2)^2)"
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = "1/cos(t/2)^2"
        s2 = "(cos(t/2)^-1)^2"
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = "cos(t/2)^-1"
        s2 = "sec(t/2)"
        calc.perform_rule(rules.OnLocation(rules.ApplyIdentity(s1, s2), "0.0.0"))
        calc.perform_rule(rules.OnLocation(rules.ExpandPolynomial(), "0.1"))
        s1 = "1 / cos(t / 2) ^ 2 * sin(t / 2) ^ 2"
        s2 = "(sin(1/2 * t) / cos(1/2*t)) ^ 2"
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = "sin(1/2 * t) / cos(1/2*t)"
        s2 = "tan(1/2*t)"
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        s1 = "2 * y / cos(t / 2) * sin(t / 2)"
        s2 = "2 * y * (sin(1/2 * t) / cos(1/2 * t))"
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = "sin(1/2 * t) / cos(1/2 * t)"
        s2 = "tan(1/2 * t)"
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        calc.perform_rule(rules.Substitution("u", "tan(t/2)"))
        calc.perform_rule(rules.Simplify())
        s1 = "2 * u * y + u ^ 2 + 1"
        s2 = "(u+y)^2 + (sqrt(1-y^2)^2)"
        calc.perform_rule(rules.Equation(s1, s2))
        calc.perform_rule(rules.Substitution("s", "u+y"))
        s1 = "(s ^ 2 - y ^ 2 + 1)"
        s2 = "s^2 + (sqrt(1-y^2)^2)"
        calc.perform_rule(rules.Equation(s1, s2))
        calc.perform_rule(rules.Substitution("v", "s/sqrt(1-y^2)"))
        calc.perform_rule(rules.Simplify())
        s1 = "v ^ 2 * (-(y ^ 2) + 1) - y ^ 2 + 1"
        s2 = "(v^2+1)*(1-y^2)"
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = "1 / ((v ^ 2 + 1) * (1 - y ^ 2))"
        s2 = "1/(v^2+1) * (1-y^2)^-1"
        calc.perform_rule(rules.Equation(s1, s2))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        # application
        s1 = "(INT y:[-1, 1]. INT x:[0,1]. 1/((1+y*x)*sqrt(1-x^2)))"
        s2 = "pi^2 / 2"
        goal02 = file.add_goal(s1+"="+s2)
        proof = goal02.proof_by_calculation()
        calc = proof.lhs_calc
        s = calc.parse_expr("INT x:[0,1]. 1 / ((1 + y * x) * sqrt(1 - x ^ 2))")
        calc.perform_rule(rules.ApplyEquation(goal01.goal, s))
        s1 = "atan((y + 1) / sqrt(-(y ^ 2) + 1)) - atan(y / sqrt(-(y ^ 2) + 1))"
        s2 = "atan(((y + 1) / sqrt(-(y ^ 2) + 1) - (y / sqrt(-(y ^ 2) + 1))) / \
        (1 + (y + 1) / sqrt(-(y ^ 2) + 1) * (y / sqrt(-(y ^ 2) + 1))))"
        calc.perform_rule(rules.ApplyIdentity(s1,s2))
        s1 = "(y + 1) / sqrt(-(y ^ 2) + 1) - y / sqrt(-(y ^ 2) + 1)"
        s2 = "1 / sqrt(-(y ^ 2) + 1)"
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = "(y + 1) / sqrt(-(y ^ 2) + 1) * (y / sqrt(-(y ^ 2) + 1))"
        s2 = "y*(y+1) / (1-y^2)"
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = "y * (y + 1) / (1 - y ^ 2)"
        s2 = "y / (1-y)"
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = "1 + y / (1 - y)"
        s2 = "1 / (1-y)"
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = "1 / sqrt(-(y ^ 2) + 1) / (1 / (1 - y))"
        s2 = "(1-y) / sqrt((1-y)*(1+y))"
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = "(1-y) / sqrt((1-y)*(1+y))"
        s2 = "sqrt((1-y)^2/((1-y)*(1+y)))"
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = "(1-y)^2/((1-y)*(1+y))"
        s2 = "(1-y) / (1+y)"
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = "atan(sqrt((1 - y) / (1 + y)))"
        s2 = "acos(y) / 2"
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Substitution("x", "acos(y)"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        # self.checkAndOutput(file)

    def testElementaryLogIntegral(self):
        # Reference: Impossible, Integrals, Sums, and Series
        # section 1.2
        file = compstate.CompFile("impossible", "elementary_log")
        file.add_definition("I(m,n) = (INT x:[0,1]. x^m * log(x)^n)")
        s1 = "I(m,n)"
        s2 = "-(n / (m + 1)) * I(m, n - 1)"
        goal01 = file.add_goal(s1 + "=" + s2, conds=["m >= 0", "n >= 0", "isInt(n)", "isInt(m)"])
        proof = goal01.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.ExpandDefinition("I"))
        u = "log(x)^n"
        v = "x^(m+1) / (m+1)"
        calc.perform_rule(rules.IntegrationByParts(u, v))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.OnLocation(rules.FoldDefinition("I"), "0.1"))
        calc.perform_rule(rules.Equation("-(n / (m + 1) * I(m,n - 1))", "-(n / (m + 1)) * I(m, n - 1)"))
        calc = proof.rhs_calc
        assert goal01.is_finished()
        s1 = "I(m,n)"
        s2 = "(-1)^n * (factorial(n) / (m+1)^(n+1))"
        goal02 = file.add_goal(s1 + "=" + s2, conds=["m >= 0", "n >= 0", "isInt(n)", "isInt(m)"])
        proof = goal02.proof_by_induction("n")
        base_proof = proof.base_case.proof_by_calculation()
        calc = base_proof.lhs_calc
        calc.perform_rule(rules.ExpandDefinition("I"))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())

        induct_proof = proof.induct_case.proof_by_calculation()
        calc = induct_proof.lhs_calc
        s = calc.parse_expr("I(m,n + 1)")
        calc.perform_rule(rules.ApplyEquation(goal01.goal, s))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.OnLocation(rules.ApplyInductHyp(), "0.0.1"))
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("-((-1) ^ n * factorial(n) * (m + 1) ^ (-n - 2) * (n + 1))")
        s2 = calc.parse_expr("(-1)^(n+1) * ((n+1) * factorial(n)) / (m+1)^(n+2)")
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = calc.parse_expr("(n+1) * factorial(n)")
        s2 = calc.parse_expr("factorial(n+1)")
        calc.perform_rule(rules.OnLocation(rules.ApplyIdentity(s1, s2), "0.1"))
        calc.perform_rule(rules.Simplify())
        assert goal02.is_finished()

        s1 = "(INT x:[0,a]. x^m *log(x)^n)"
        s2 = "a^(m+1) * SUM(k, 0, n, (-1)^k*binom(n, k)*factorial(k)*log(a)^(n-k)/(m+1)^(k+1))"
        s = parser.parse_expr(s1+"="+s2)
        goal03 = file.add_goal(s, conds=["m>=0", "n>=0", "isInt(n)", "isInt(m)", "a>0"])
        proof = goal03.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.SubstitutionInverse("y", "x", "a*y"))
        source = calc.parse_expr("(a * y) ^ m")
        target = calc.parse_expr("a^m * y^m")
        calc.perform_rule(rules.ApplyIdentity(source,target))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.Equation("a*y", "y*a"))
        s1 = "log(y*a)"
        s2 = "log(y) + log(a)"
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        calc.perform_rule(rules.OnLocation(rules.SeriesExpansionIdentity(index_var = "k"), "1.0.1"))
        s1 = calc.parse_expr("y ^ m * SUM(k, 0, n, binom(n,k) * log(y) ^ k * log(a) ^ (n - k))")
        s2 = calc.parse_expr("SUM(k, 0, n, y^m * binom(n,k) * log(y) ^ k * log(a) ^ (n - k))")
        calc.perform_rule(rules.Equation(s1, s2))
        calc.perform_rule(rules.OnLocation(rules.IntSumExchange(), "1"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.OnLocation(rules.FoldDefinition("I"), "1.0.1"))
        s = calc.parse_expr("I(m,k)")
        calc.perform_rule(rules.ApplyEquation(goal02.goal, s))
        calc.perform_rule(rules.Simplify())
        calc = proof.rhs_calc
        calc.perform_rule(rules.Simplify())
        assert goal03.is_finished()
        # self.checkAndOutput(file)

    def testHarmonicSeries(self):
        # Reference: Impossible, Integrals, Sums, and Series
        # section 1.3
        file = compstate.CompFile("impossible", "harmonic_series")
        file.add_definition("H(n) = SUM(k, 1, n, 1/k)")
        file.add_definition("I(n) = (INT x:[0,1]. x^(n-1) * log(1-x))")

        goal01 = file.add_goal("I(n) = -(H(n)/n)", conds=["n>0", "isInt(n)"])
        proof = goal01.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.ExpandDefinition("I"))
        u = "log(1-x)"
        v = "(x^n - 1)/n"
        calc.perform_rule(rules.IntegrationByParts(u, v))
        calc.perform_rule(rules.Simplify())
        s1 = "(x ^ n - 1) / (-x + 1)"
        s2 = "-((1-x^n)/(1-x))"
        calc.perform_rule(rules.Equation(s1, s2))
        calc.perform_rule(rules.OnLocation(rules.SeriesExpansionIdentity(index_var="k"), "1.0.0"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.OnLocation(rules.IntSumExchange(), "0.1"))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.OnLocation(rules.FoldDefinition("H"), "0.1"))
        calc.perform_rule(rules.Simplify())

        # self.checkAndOutput(file)

    def testUsefulLogIntegral(self):
        # Reference: Impossible, Integrals, Sums, and Series
        # section 1.4
        file = compstate.CompFile("impossible", "useful_log")

        file.add_definition("Li(s, x) = SUM(k, 1, oo, x^k /k^s)")
        goal = file.add_goal("x/(1-x) = SUM(k, 1, oo, x^k)", conds=["abs(x) < 1"])
        split_cond = parser.parse_expr("x != 0")
        proof = goal.proof_by_case(split_cond)
        proofa = proof.cases[0].proof_by_calculation()
        calc = proofa.lhs_calc
        calc.perform_rule(rules.Equation("x/(1-x)", "x * (1-x)^(-1)"))
        calc.perform_rule(rules.OnLocation(rules.SeriesExpansionIdentity(), "1"))
        s1 = calc.parse_expr("x * SUM(n, 0, oo, x ^ n)")
        s2 = calc.parse_expr("SUM(n, 0, oo, x * x^n)")
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = calc.parse_expr("x*x^n")
        s2 = calc.parse_expr("x^1 * x^n")
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = calc.parse_expr("x^1 * x^n")
        s2 = calc.parse_expr("x^(n+1)")
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        calc.perform_rule(rules.ChangeSummationIndex("1"))
        calc.perform_rule(rules.Simplify())
        proofb = proof.cases[1].proof_by_calculation()
        calc = proofb.lhs_calc
        calc.perform_rule(rules.Simplify())
        calc = proofb.rhs_calc
        calc.perform_rule(rules.Simplify())
        self.assertTrue(goal.is_finished())
        goal = file.add_goal("Li(0, x) = x/(1-x)", conds=["abs(x) < 1"])
        proof = goal.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.ExpandDefinition("Li"))
        s = calc.parse_expr("SUM(k, 1, oo, x ^ k)")
        calc.perform_rule(rules.ApplyEquation("x/(1-x) = SUM(k, 1, oo, x^k)", s))
        self.assertTrue(goal.is_finished())

        goal = file.add_goal("Li(s+1, x) = (INT t:[0, x]. Li(s, t) / t)", conds=["abs(x)<1", "s>=0", "isInt(s)"])
        proof = goal.proof_by_induction("s")
        proof_base = proof.base_case.proof_by_calculation()
        calc = proof_base.lhs_calc
        calc.perform_rule(rules.ExpandDefinition("Li"))
        calc = proof_base.rhs_calc
        eq = parser.parse_expr("Li(0,x)=x/(1-x)")
        s = calc.parse_expr("Li(0,t)")
        calc.perform_rule(rules.ApplyEquation(eq, s))
        eq = parser.parse_expr("x/(1-x) = SUM(k, 1, oo, x^k)")
        s = calc.parse_expr("t / (1 - t)")
        calc.perform_rule(rules.ApplyEquation(eq, s))
        s1 = calc.parse_expr("SUM(k, 1, oo, t ^ k) / t")
        s2 = calc.parse_expr("1/t * SUM(k, 1, oo, t ^ k)")
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = calc.parse_expr("1/t * SUM(k, 1, oo, t ^ k)")
        s2 = calc.parse_expr("SUM(k, 1, oo, 1/t * t ^ k)")
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = calc.parse_expr("1/t * t ^ k")
        s2 = calc.parse_expr("t^(-1) * t^k")
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = calc.parse_expr("t^(-1)*t^k")
        s2 = calc.parse_expr("t^(-1+k)")
        calc.perform_rule(rules.OnLocation(rules.ApplyIdentity(s1, s2),"0.0"))
        calc.perform_rule(rules.IntSumExchange())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertTrue(proof_base.is_finished())
        proof_induct = proof.induct_case.proof_by_calculation()
        calc = proof_induct.lhs_calc
        calc.perform_rule(rules.ExpandDefinition("Li"))
        calc = proof_induct.rhs_calc
        calc.perform_rule(rules.OnLocation(rules.ExpandDefinition("Li"), "0.0"))
        s1 = calc.parse_expr("SUM(k, 1, oo, k ^ (-s - 1) * t ^ k) / t")
        s2 = calc.parse_expr("t^-1 * SUM(k, 1, oo, k ^ (-s - 1) * t ^ k)")
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = calc.parse_expr("t^-1 * SUM(k, 1, oo, k ^ (-s - 1) * t ^ k)")
        s2 = calc.parse_expr("SUM(k, 1, oo, t^-1 * (k ^ (-s - 1) * t ^ k))")
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = calc.parse_expr("t^-1 * (k ^ (-s - 1) * t ^ k)")
        s2 = calc.parse_expr("t^-1 * t^k * k^(-s-1)")
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = calc.parse_expr("t^-1 * t^k")
        s2 = calc.parse_expr("t^(-1 + k)")
        calc.perform_rule(rules.OnLocation(rules.ApplyIdentity(s1, s2),"0.0.0"))
        calc.perform_rule(rules.IntSumExchange())
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        self.assertTrue(goal.is_finished())

        goal = file.add_goal("(D x. Li(2, 1-x)) = log(x)/(1-x)", conds=["x < 1", "x>0"])
        proof = goal.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.OnLocation(rules.ExpandDefinition("Li"),"0"))
        calc.perform_rule(rules.Simplify())
        calc = proof.rhs_calc
        s1 = parser.parse_expr("log(x)")
        s2 = parser.parse_expr("log(1-(1-x))")
        calc.perform_rule(rules.Equation(s1, s2))
        calc.perform_rule(rules.OnLocation(rules.SeriesExpansionIdentity(), "0"))
        s1 = calc.parse_expr("SUM(n, 0, oo, (-1) ^ n * (-(1 - x)) ^ (n + 1) / (n + 1)) / (1 - x)")
        s2 = calc.parse_expr("-(-(1-x))^(-1) * SUM(n, 0, oo, (-1) ^ n * (-(1 - x)) ^ (n + 1) / (n + 1))")
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = calc.parse_expr("-(-(1-x))^(-1) * SUM(n, 0, oo, (-1) ^ n * (-(1 - x)) ^ (n + 1) / (n + 1))")
        s2 = calc.parse_expr("SUM(n, 0, oo, -(-(1-x))^(-1) * ((-1) ^ n * (-(1 - x)) ^ (n + 1) / (n + 1)))")
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = calc.parse_expr(" -(-(1-x))^(-1) * ((-1) ^ n * (-(1 - x)) ^ (n + 1) / (n + 1))")
        s2 = calc.parse_expr("-1 * ((-(1-x))^(-1) * (-(1 - x)) ^ (n + 1)) * (-1) ^ n  / (n + 1)")
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = calc.parse_expr("(-(1 - x)) ^ (-1) * (-(1 - x)) ^ (n + 1)")
        s2 = calc.parse_expr('(-(1-x))^(-1 + (n+1))')
        calc.perform_rule(rules.ApplyIdentity(s1, s2))
        calc.perform_rule(rules.ChangeSummationIndex("1"))
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("x-1")
        s2 = calc.parse_expr("-1 * (-x+1)")
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = calc.parse_expr("(-1 * (-x + 1)) ^ (n - 1)")
        s2 = calc.parse_expr("(-1)^(n-1) * (-x+1)^(n-1)")
        calc.perform_rule(rules.OnLocation(rules.ApplyIdentity(s1, s2), "0.0.1"))
        calc.perform_rule(rules.Simplify())
        self.assertTrue(goal.is_finished())

        goal = file.add_goal("(D x. -Li(3, 1-x)) = Li(2, -x+1) / (-x+1)", conds = ["x>0","x<1"])
        proof = goal.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.Equation("3", "2+1"))
        eq = parser.parse_expr("Li(s+1, x) = (INT t:[0, x]. Li(s, t) / t)")
        s = calc.parse_expr("Li(2 + 1,1 - x)")
        calc.perform_rule(rules.ApplyEquation(eq, s))
        calc.perform_rule(rules.Simplify())
        self.assertTrue(goal.is_finished())

        goal = file.add_goal("Li(s,1) = zeta(s)", conds=["isInt(s)", "s>1"])
        proof = goal.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.ExpandDefinition("Li"))
        calc = proof.rhs_calc
        calc.perform_rule(rules.ExpandDefinition("zeta"))
        calc.perform_rule(rules.ChangeSummationIndex("1"))
        calc.perform_rule(rules.Simplify())
        self.assertTrue(goal.is_finished())

        s1 = "(INT t:[0,x]. log(1-t)^2 / t)"
        s2 = "log(x) * log(1-x)^2 + 2 * log(1-x) * \
                Li(2, 1-x) - 2*Li(3,1-x) + 2 * zeta(3)"
        goal = file.add_goal(s1+'='+s2, conds=["x>0", "x<1"])
        proof = goal.proof_by_calculation()
        calc = proof.lhs_calc
        u = parser.parse_expr("log(1-t)^2")
        v = parser.parse_expr("log(t)")
        calc.perform_rule(rules.IntegrationByParts(u, v))
        calc.perform_rule(rules.Simplify())
        # (D x. Li(2, 1-x)) = log(x)/(1-x)
        s1 = parser.parse_expr(" log(t) * log(-t + 1) / (-t + 1)")
        s2 = parser.parse_expr("log(t) / (1-t) * log(-t+1)")
        calc.perform_rule(rules.Equation(s1, s2))
        eq = parser.parse_expr("(D x. Li(2, 1-x)) = log(x)/(1-x)")
        s = calc.parse_expr("log(t) / (1 - t)")
        calc.perform_rule(rules.ApplyEquation(eq, s))
        u = parser.parse_expr("log(-t+1)")
        v = parser.parse_expr("Li(2, 1-t)")
        calc.perform_rule(rules.OnLocation(rules.IntegrationByParts(u, v), "0.1"))
        calc.perform_rule(rules.Simplify())
        eq = parser.parse_expr("(D x. -Li(3, 1-x)) = Li(2, -x+1) / (-x+1)")
        s = calc.parse_expr("Li(2,-t + 1) / (-t + 1)")
        calc.perform_rule(rules.ApplyEquation(eq, s))
        calc.perform_rule(rules.Simplify())
        eq = parser.parse_expr("Li(s,1) = zeta(s)")
        s = calc.parse_expr("Li(3,1)")
        calc.perform_rule(rules.OnLocation(rules.ApplyEquation(eq, s), "1.1"))
        calc = proof.rhs_calc
        calc.perform_rule(rules.Simplify())
        self.assertTrue(goal.is_finished())

        s1 = "(D x. Li(2, 1/(1+x)))"
        s2 = "log(x/(1+x)) / (1+x)"
        goal = file.add_goal(s1+"="+s2, conds=['x>0', 'x<1'])
        proof = goal.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.OnLocation(rules.ExpandDefinition("Li"), "0"))
        calc.perform_rule(rules.Simplify())
        calc = proof.rhs_calc
        s1 = calc.parse_expr("x/(1+x)")
        s2 = calc.parse_expr("1 - 1/(1+x)")
        calc.perform_rule(rules.Equation(s1, s2))
        calc.perform_rule(rules.OnLocation(rules.SeriesExpansionIdentity(), "0"))
        s1 = calc.parse_expr("-(1/(1+x))")
        s2 = calc.parse_expr("(-1) * (1/(1+x))")
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = calc.parse_expr("((-1) * (1/(1+x)))^(n+1)")
        s2 = calc.parse_expr("(-1)^(n+1) * (1/(1+x))^(n+1)")
        calc.perform_rule(rules.OnLocation(rules.ApplyIdentity(s1, s2), "0.0.0.1"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.OnLocation(rules.ChangeSummationIndex('1'), '0.1'))
        calc.perform_rule(rules.Simplify())
        s1 = calc.parse_expr("(1 / (x + 1)) ^ n")
        s2 = calc.parse_expr("(1 / (x + 1)) ^ (n-1 + 1)")
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = calc.parse_expr("(1 / (x + 1)) ^ (n-1 + 1)")
        s2 = calc.parse_expr("(1 / (x + 1)) ^ (n-1) * (1/(x+1))^1")
        calc.perform_rule(rules.OnLocation(rules.ApplyIdentity(s1, s2), "0.1.0.1"))
        calc.perform_rule(rules.Simplify())
        self.assertTrue(goal.is_finished())


        s1 = "1 / (x + 1) * Li(2,1 / (x + 1))"
        s2 = "(D x. -Li(3, 1/(x+1)))"
        goal = file.add_goal(s1+"="+s2, conds=["x>0", "x<1"])
        proof = goal.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.OnLocation(rules.ExpandDefinition("Li"), "1"))
        s1 = calc.parse_expr("(1 / (x + 1)) ^ k")
        s2 = calc.parse_expr("(1 / (x + 1)) ^ (k-1 + 1)")
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = calc.parse_expr("(1 / (x + 1)) ^ (k-1 + 1)")
        s2 = calc.parse_expr("(1 / (x + 1)) ^ (k-1) * (1/(x+1))^1")
        calc.perform_rule(rules.OnLocation(rules.ApplyIdentity(s1, s2), "1.0.1"))
        calc.perform_rule(rules.Simplify())
        calc = proof.rhs_calc
        calc.perform_rule(rules.OnLocation(rules.ExpandDefinition("Li"), "0.0"))
        calc.perform_rule(rules.Simplify())
        self.assertTrue(goal.is_finished())

        s1 = parser.parse_expr("(INT t:[0, x]. log(1+t)^2 / t)")
        s2 = parser.parse_expr("log(x)*log(1+x)^2-(2/3)*log(1+x)^3-\
              2*log(1+x)*Li(2,1/(1+x))-2*Li(3,1/(1+x))+2*zeta(3)")
        goal = file.add_goal(expr.Op('=', s1, s2), conds=["x>0", "x<1"])
        proof = goal.proof_by_calculation()
        calc = proof.lhs_calc
        u = parser.parse_expr("log(1+t)^2")
        v = parser.parse_expr("log(t)")
        calc.perform_rule(rules.IntegrationByParts(u, v))
        calc.perform_rule(rules.Simplify())
        s1 = parser.parse_expr("log(t)")
        s2 = parser.parse_expr("log((1+t) * (t/(1+t)))")
        calc.perform_rule(rules.Equation(s1, s2))
        s1 = parser.parse_expr("log((1+t) * (t/(1+t)))")
        s2 = parser.parse_expr("log(1+t) + log(t/(1+t))")
        calc.perform_rule(rules.OnLocation(rules.ApplyIdentity(s1, s2), "0.0.1.0.0.0"))
        calc.perform_rule(rules.OnLocation(rules.ExpandPolynomial(), "0.0.1.0.0"))
        calc.perform_rule(rules.OnLocation(rules.ExpandPolynomial(), "0.0.1.0"))
        calc.perform_rule(rules.Simplify())
        calc.perform_rule(rules.OnLocation(rules.Substitution('u', 'log(t+1)'), '0.1.1'))
        calc.perform_rule(rules.DefiniteIntegralIdentity())
        calc.perform_rule(rules.Simplify())
        s1 = parser.parse_expr("log(t + 1) / (t + 1) * log(t / (t + 1))")
        s2 = parser.parse_expr("log(t/(1+t)) / (1+t) * log(t+1)")
        calc.perform_rule(rules.Equation(s1, s2))
        eq = parser.parse_expr("(D x. Li(2, 1/(1+x))) = log(x/(1+x)) / (1+x)")
        s = calc.parse_expr("log(t / (1 + t)) / (1 + t)")
        calc.perform_rule(rules.ApplyEquation(eq, s))
        u = parser.parse_expr("log(t+1)")
        v = parser.parse_expr("Li(2, 1/(1+t))")
        calc.perform_rule(rules.OnLocation(rules.IntegrationByParts(u, v), "0.0.0.1"))
        calc.perform_rule(rules.Simplify())
        eq = parser.parse_expr("1 / (x + 1) * Li(2,1 / (x + 1)) = (D x. -Li(3, 1/(x+1)))")
        s = calc.parse_expr("1 / (t + 1) * Li(2,1 / (t + 1))")
        calc.perform_rule(rules.ApplyEquation(eq, s))
        calc.perform_rule(rules.Simplify())
        eq = parser.parse_expr("Li(s,1) = zeta(s)")
        s = calc.parse_expr("Li(3,1)")
        calc.perform_rule(rules.ApplyEquation(eq, s))
        calc = proof.rhs_calc
        calc.perform_rule(rules.Simplify())
        self.assertTrue(goal.is_finished())
        # self.checkAndOutput(file)


if __name__ == "__main__":
    unittest.main()
