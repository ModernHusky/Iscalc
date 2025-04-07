"""Overall test for integrals."""

import unittest
import json

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
