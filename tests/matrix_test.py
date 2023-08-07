import json
import unittest
from typing import List

from integral import rules, context, parser, compstate
from integral.expr import Op, Var, Const, Matrix, Vector, Expr, Fun


class MatrixTest(unittest.TestCase):
    def checkAndOutput(self, file: compstate.CompFile, omit_finish: bool = False):
        # Test parsing of json file
        json_file = file.export()
        for i, item in enumerate(json_file['content']):
            self.assertEqual(compstate.parse_item(file.content[i].parent, item).export(), file.content[i].export(), "%d:%sasdfasdf%s"%(i, compstate.parse_item(file.content[i].parent, item).export(), file.content[i].export()))

        # Output to file
        with open('../examples/' + file.name + '.json', 'w', encoding='utf-8') as f:
            json.dump(file.export(), f, indent=4, ensure_ascii=False, sort_keys=True)

        # Test goals are finished
        if not omit_finish:
            for content in file.content:
                self.assertTrue(content.is_finished())
    r = rules.FullSimplify()
    ctx = context.Context()
    def testTranspose1(self):
        m = Matrix((2,3), [[1,2,3],[4,5,6]])
        res = "[1, 4]\n[2, 5]\n[3, 6]"
        assert str(m.t) == res

    def testTranspose2(self):
        m = Matrix((3, 3), [[1, 2, 3], [4, 5, 6], [7,8,9]])
        res = "[1, 4, 7]\n[2, 5, 8]\n[3, 6, 9]"
        assert str(m.t) == res

    def testMultiplication1(self):
        m1 = Matrix((2,3), [[Const(1), Const(3), Const(5)],
                            [Const(2), Const(4), Const(7)]])
        m2 = Matrix((3,3), [[Const(-5), Const(8), Const(11)],
                            [Const(3), Const(9), Const(21)],
                            [Const(4), Const(0), Const(8)]])
        res = "[1 * -5 + 3 * 3 + 5 * 4, 1 * 8 + 3 * 9 + 5 * 0, 1 * 11 + 3 * 21 + 5 * 8]\n\
[2 * -5 + 4 * 3 + 7 * 4, 2 * 8 + 4 * 9 + 7 * 0, 2 * 11 + 4 * 21 + 7 * 8]"
        assert res == str(m1 * m2)
        res = self.r.eval(m1 * m2, self.ctx)
        s = "[24, 35, 114]\n[30, 52, 162]"
        assert str(res) == s



    def testMultiplication2(self):
        v1 = Vector([Const(1), Const(1), Const(0), Const(0)], is_column=False)
        v2 = Vector([Const(1), Const(2), Const(3), Const(4)])
        res = "1 * 1 + 1 * 2 + 0 * 3 + 0 * 4"
        self.assertEqual(str(v1 * v2), res)
        self.assertIsInstance(v1*v2, Expr)

    def testMultiplication3(self):
        m = Matrix((2,3), [[Const(1), Const(3), Const(5)],
                            [Const(2), Const(4), Const(7)]])
        v = Vector([Const(1), Const(2), Const(1)])
        res = m * v
        s = "[1 * 1 + 3 * 2 + 5 * 1]\n\
[2 * 1 + 4 * 2 + 7 * 1]"
        self.assertIsInstance(res, Vector)
        self.assertTrue(res.is_column)
        self.assertEqual(str(res), s)
        res = self.r.eval(res, self.ctx)
        e = parser.parse_expr("1 * 1 + 3 * 2 + 5 * 1")
        self.assertEqual(str(res), "[12]\n[17]")


    def testMultiplication4(self):
        v = Vector([Const(1), Const(2)], is_column=False)
        m = Matrix((2, 3), [[Const(1), Const(3), Const(5)],
                            [Const(2), Const(4), Const(7)]])

        res = v * m
        s = "[1 * 1 + 2 * 2, 1 * 3 + 2 * 4, 1 * 5 + 2 * 7]"
        self.assertIsInstance(res, Vector)
        self.assertTrue(res.is_row)
        self.assertEqual(str(res), s)

    def testHatAndVee1(self):
        v = Vector([Var('a1'), Var('a2'), Var('a3')])
        s = "[0, -a3, a2]\n\
[a3, 0, -a1]\n\
[-a2, a1, 0]"
        self.assertEqual(str(v.hat), s)
        self.assertEqual(str(v.hat.vee), str(v))

    def testHatAndVee2(self):
        v = Vector([Var('v1'), Var('v2'), Var('v3'), Var('w1'), Var('w2'), Var('w3')])
        s = "[0, -w3, w2, v1]\n\
[w3, 0, -w1, v2]\n\
[-w2, w1, 0, v3]\n\
[0, 0, 0, 0]"
        self.assertEqual(str(v.hat), s)
        self.assertEqual(str(v.hat.vee), str(v))

    def testExample02(self):
        file = compstate.CompFile("matrix", "matrix_example02")
        goal = file.add_goal('hat({a1,a2,a3})^2 = \
        T({a1,a2,a3}) * {a1,a2,a3} - norm({a1,a2,a3})^2 * unit_matrix(3)')
        proof = goal.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.OnLocation(rules.ExpandMatVecFunc(), "0"))
        olde = "{{0, -a3, a2}, {a3, 0, -a1}, {-a2, a1, 0}}^2"
        newe = "{{0, -a3, a2}, {a3, 0, -a1}, {-a2, a1, 0}}^1 * {{0, -a3, a2},{a3, 0, -a1},{-a2, a1, 0}}^(2-1)"
        calc.perform_rule(rules.ApplyIdentity(olde, newe))
        calc.perform_rule(rules.FullSimplify())
        calc = proof.rhs_calc
        calc.perform_rule(rules.FullSimplify())
        calc.perform_rule(rules.OnLocation(rules.MatrixRewrite(), '1'))
        calc.perform_rule(rules.FullSimplify())
        calc.perform_rule(rules.OnLocation(rules.MatrixRewrite(), '0'))
        calc.perform_rule(rules.FullSimplify())
        self.checkAndOutput(file)

    def testExample03(self):
        file = compstate.CompFile("matrix", "matrix_example03")
        goal = file.add_goal("T({a1,a2,a3})*{a1,a2,a3}*hat({a1,a2,a3}) = zero_matrix(3,3)")
        proof = goal.proof_by_calculation()
        calc = proof.lhs_calc
        calc.perform_rule(rules.FullSimplify())
        calc = proof.rhs_calc
        calc.perform_rule(rules.FullSimplify())
        self.checkAndOutput(file)

    def testExample04(self):
        file = compstate.CompFile("MIRM", "matrix_example04")
        file.add_assumption("norm({w1,w2,w3}) = 1")
        file.add_assumption("isInt(n)")
        file.add_assumption("n>=0")
        goal = file.add_goal("hat({w1,w2,w3})^(2*n+1) = (-1)^n * hat({w1,w2,w3})")
        proof = goal.proof_by_induction('n', 0)
        base_proof = proof.base_case.proof_by_calculation()
        induct_proof = proof.induct_case.proof_by_calculation()
        calc = induct_proof.lhs_calc
        calc.perform_rule(rules.Equation("(2 * n + 3)", "2 + (2 * n + 1)"))
        calc.perform_rule(rules.ApplyIdentity("hat({w1,w2,w3})^(2 + (2 * n + 1))", "hat({w1,w2,w3})^2* hat({w1,w2,w3})^(2 * n + 1)"))
        calc.perform_rule(rules.OnLocation(rules.ApplyInductHyp(),'1'))
        eq = "hat({a1, a2, a3}) ^ 2 = T({a1, a2, a3}) * {a1, a2,a3} - norm({a1, a2, a3}) ^ 2 * unit_matrix(3)"
        calc.perform_rule(rules.OnLocation(rules.ApplyEquation(eq), "0"))
        eq = "norm({w1,w2,w3}) = 1"
        calc.perform_rule(rules.OnLocation(rules.ApplyEquation(eq), "0.1.0.0"))
        calc.perform_rule(rules.ExpandPolynomial())
        old_e = "(-1) ^ n * T({w1, w2, w3}) * {w1, w2, w3} * hat({w1, w2, w3})"
        new_e = "(-1) ^ n * (T({w1, w2, w3}) * {w1, w2, w3} * hat({w1, w2, w3}))"
        calc.perform_rule(rules.Equation(old_e, new_e))
        eq = "T({a1, a2, a3}) * {a1, a2, a3} * hat({a1, a2, a3}) = zero_matrix(3,3)"
        calc.perform_rule(rules.OnLocation(rules.ApplyEquation(eq), "0.1"))
        calc.perform_rule(rules.FullSimplify())
        old_e = "(-1) ^ n * {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}} * {{0, -w3, w2}, {w3, 0, -w1}, {-w2, w1, 0}}"
        new_e = "(-1) ^ n * ({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}} * {{0, -w3, w2}, {w3, 0, -w1}, {-w2, w1, 0}})"
        calc.perform_rule(rules.Equation(old_e, new_e))
        calc.perform_rule(rules.OnLocation(rules.FullSimplify(), "0.0.1"))
        calc.perform_rule(rules.OnLocation(rules.MatrixRewrite(), "1"))
        old_e = "-((-1) ^ n * {{0, -w3, w2}, {w3, 0, -w1}, {-w2, w1, 0}})"
        new_e = "(-1) ^(n+1) * {{0, -w3, w2}, {w3, 0, -w1}, {-w2, w1, 0}}"
        calc.perform_rule(rules.Equation(old_e, new_e))
        calc.perform_rule(rules.OnLocation(rules.MatrixRewrite(), "0"))
        calc.perform_rule(rules.FullSimplify())
        calc = induct_proof.rhs_calc
        calc.perform_rule(rules.FullSimplify())
        old_e = "-((-1) ^ n * {{0, -w3, w2}, {w3, 0, -w1}, {-w2, w1, 0}})"
        new_e = "(-1) ^(n+1) * {{0, -w3, w2}, {w3, 0, -w1}, {-w2, w1, 0}}"
        calc.perform_rule(rules.Equation(old_e, new_e))
        calc.perform_rule(rules.MatrixRewrite())
        calc.perform_rule(rules.FullSimplify())
        assert proof.is_finished()
        self.checkAndOutput(file)

    def testExample01(self):
        file = compstate.CompFile("MIRM", "matrix_example01")
        file.add_definition("matrix P[n][n]")
        file.add_definition("matrix A[n][n]")
        file.add_definition("int n")
        file.add_assumption("invertible(P)")
        file.add_assumption("n > 0")
        # find same name variable definition
        # and then use that definition
        goal = file.add_goal("(inv(P)*A*P)^n = inv(P)*A^n*P")
        proof = goal.proof_by_induction(induct_var='n', start=0)
        base_proof = proof.base_case.proof_by_calculation()
        induct_proof = proof.induct_case.proof_by_calculation()
        calc = base_proof.lhs_calc
        calc = base_proof.rhs_calc
        calc = induct_proof.lhs_calc
        old_expr = "(inv(matrix P[n][n]) * matrix A[n][n] * matrix P[n][n]) ^ (int n + 1)"
        new_expr = "(inv(matrix P[n][n]) * matrix A[n][n] * matrix P[n][n]) ^ (int n) * (inv(matrix P[n][n]) * matrix A[n][n] * matrix P[n][n])"
        calc.perform_rule(rules.ApplyIdentity(old_expr, new_expr))
        calc.perform_rule(rules.OnSubterm(rules.ApplyInductHyp()))
        calc.perform_rule(rules.FullSimplify())
        old_expr = "inv(matrix P[n][n]) * matrix A[n][n] ^ int n * matrix A[n][n] * matrix P[n][n]"
        new_expr = "inv(matrix P[n][n]) * (matrix A[n][n] ^ int n * matrix A[n][n]) * matrix P[n][n]"
        calc.perform_rule(rules.Equation(old_expr, new_expr))
        old_expr = "matrix A[n][n] ^ int n * matrix A[n][n]"
        new_expr = "matrix A[n][n] ^ (int n + 1)"
        calc.perform_rule(rules.ApplyIdentity(old_expr, new_expr))
        self.checkAndOutput(file)

    def testNormalize(self):
        e = parser.parse_expr("matrix A[n][n]")
        e = e ^ Const(0)
        from integral import poly
        from integral.poly import normalize
        ctx = context.Context()
        e = normalize(e, ctx)
        self.assertEqual(str(e), "unit_matrix(n)")

    def testGetType(self):
        e = parser.parse_expr("inv(matrix P[n][n]) * matrix A[n][n] * matrix P[n][n]")
        ty = e.get_type()
        print(ty)
        e = e ^ Const(0)
        from integral import poly
        from integral.poly import normalize
        ctx = context.Context()
        e = normalize(e, ctx)
        print(e)


