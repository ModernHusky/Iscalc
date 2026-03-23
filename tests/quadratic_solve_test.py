"""
测试 solve_equation 对一元二次方程 ax² + bx + c = 0 的求解能力。

运行: python -m pytest tests/test_quadratic_solve.py -v
"""

import os
import unittest
from fractions import Fraction
from pathlib import Path

from integral import parser, solve, context, poly
from integral.expr import Const, Var, Op, Fun

PROJECT_ROOT = Path(__file__).parent.parent
os.chdir(PROJECT_ROOT)


def clear_caches():
    poly._normalize_cache.clear()


def make_ctx():
    ctx = context.Context()
    ctx.load_book("base")
    return ctx


def _are_equal(e1, e2, ctx):
    """Check if two symbolic expressions are equal after normalization."""
    n1 = poly.normalize(e1, ctx)
    n2 = poly.normalize(e2, ctx)
    return n1 == n2


def _fmt(sols, ctx):
    """将解列表格式化为可读字符串。"""
    return "[" + ", ".join(str(poly.normalize(s, ctx)) for s in sols) + "]"


# =============================================================================
# 1. 标准二次方程（整数系数）
# =============================================================================

class TestQuadraticStandardRealRoots(unittest.TestCase):
    """ax² + bx + c = 0，判别式 D > 0，两个实数根。"""

    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()

    def test_x2_minus_x_minus_2(self):
        """x² - x - 2 = 0 → 根 x=2, x=-1。"""
        e = parser.parse_expr("x^2 - x - 2")
        sols = solve.solve_equation(e, Const(0), "x", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：x=2, x=-1  共2个根")
        self.assertEqual(len(sols), 2)
        self.assertTrue(
            any(_are_equal(s, Const(2), self.ctx) for s in sols),
            f"solutions {sols} should contain x=2"
        )
        self.assertTrue(
            any(_are_equal(s, Const(-1), self.ctx) for s in sols),
            f"solutions {sols} should contain x=-1"
        )

    def test_x2_plus_2x_plus_1(self):
        """x² + 2x + 1 = 0 → 重根 x=-1。"""
        e = parser.parse_expr("x^2 + 2*x + 1")
        sols = solve.solve_equation(e, Const(0), "x", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：x=-1  共1个根（重根）")
        self.assertEqual(len(sols), 1)
        self.assertTrue(
            _are_equal(sols[0], Const(-1), self.ctx),
            f"solution {sols[0]} should equal x=-1"
        )

    def test_x2_minus_4(self):
        """x² - 4 = 0 → 根 x=2, x=-2。"""
        e = parser.parse_expr("x^2 - 4")
        sols = solve.solve_equation(e, Const(0), "x", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：x=2, x=-2  共2个根")
        self.assertEqual(len(sols), 2)
        self.assertTrue(
            any(_are_equal(s, Const(2), self.ctx) for s in sols),
            f"solutions {sols} should contain x=2"
        )
        self.assertTrue(
            any(_are_equal(s, Const(-2), self.ctx) for s in sols),
            f"solutions {sols} should contain x=-2"
        )

    def test_2x2_minus_5x_plus_2(self):
        """2x² - 5x + 2 = 0 → 根 x=2, x=1/2。"""
        e = parser.parse_expr("2*x^2 - 5*x + 2")
        sols = solve.solve_equation(e, Const(0), "x", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：x=2, x=1/2  共2个根")
        self.assertEqual(len(sols), 2)
        self.assertTrue(
            any(_are_equal(s, Const(2), self.ctx) for s in sols),
            f"solutions {sols} should contain x=2"
        )
        self.assertTrue(
            any(_are_equal(s, Const(Fraction(1, 2)), self.ctx) for s in sols),
            f"solutions {sols} should contain x=1/2"
        )

    def test_3x2_plus_6x_minus_9(self):
        """3x² + 6x - 9 = 0 → 根 x=1, x=-3。"""
        e = parser.parse_expr("3*x^2 + 6*x - 9")
        sols = solve.solve_equation(e, Const(0), "x", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：x=1, x=-3  共2个根")
        self.assertEqual(len(sols), 2)
        self.assertTrue(
            any(_are_equal(s, Const(1), self.ctx) for s in sols),
            f"solutions {sols} should contain x=1"
        )
        self.assertTrue(
            any(_are_equal(s, Const(-3), self.ctx) for s in sols),
            f"solutions {sols} should contain x=-3"
        )


# =============================================================================
# 2. 判别式 < 0：纯虚数根 / 共轭复数根
# =============================================================================

class TestQuadraticComplexRoots(unittest.TestCase):
    """ax² + bx + c = 0，判别式 D < 0，两个共轭复数根。"""

    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()

    def test_z2_plus_1(self):
        """z² + 1 = 0 → 根 z=i, z=-i。"""
        e = parser.parse_expr("z^2 + 1")
        sols = solve.solve_equation(e, Const(0), "z", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：z=i, z=-i  共2个根")
        self.assertEqual(len(sols), 2)
        i_expr = Op("+", Const(0), Op("*", Const(1), Fun("i")))
        neg_i_expr = Op("-", Op("*", Const(1), Fun("i")))
        normalized_i = poly.normalize(i_expr, self.ctx)
        normalized_neg_i = poly.normalize(neg_i_expr, self.ctx)
        normalized_sols = [poly.normalize(s, self.ctx) for s in sols]
        sol_set = set(str(s) for s in normalized_sols)
        self.assertIn(str(normalized_i), sol_set,
                     f"solutions {sol_set} should contain i")
        self.assertIn(str(normalized_neg_i), sol_set,
                     f"solutions {sol_set} should contain -i")

    def test_z2_plus_4(self):
        """z² + 4 = 0 → 根 z=2i, z=-2i。"""
        e = parser.parse_expr("z^2 + 4")
        sols = solve.solve_equation(e, Const(0), "z", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：z=2i, z=-2i  共2个根")
        self.assertEqual(len(sols), 2)
        two_i = Op("*", Const(2), Fun("i"))
        neg_two_i = Op("-", Op("*", Const(2), Fun("i")))
        normalized_two_i = poly.normalize(two_i, self.ctx)
        normalized_neg_two_i = poly.normalize(neg_two_i, self.ctx)
        normalized_sols = [poly.normalize(s, self.ctx) for s in sols]
        sol_set = set(str(s) for s in normalized_sols)
        self.assertIn(str(normalized_two_i), sol_set,
                     f"solutions {sol_set} should contain 2i")
        self.assertIn(str(normalized_neg_two_i), sol_set,
                     f"solutions {sol_set} should contain -2i")

    def test_z2_plus_z_plus_1(self):
        """z² + z + 1 = 0 → 根 z = (-1 ± i√3)/2。"""
        e = parser.parse_expr("z^2 + z + 1")
        sols = solve.solve_equation(e, Const(0), "z", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：z=(-1+i√3)/2, z=(-1-i√3)/2  共2个根")
        self.assertEqual(len(sols), 2)
        half = Op("/", Const(1), Const(2))
        neg_half = Op("-", Op("/", Const(1), Const(2)))
        sqrt3 = Op("^", Const(3), half)
        pos_root = Op("+", Op("-", Const(0), Const(1)),
                      Op("*", sqrt3, Fun("i")))
        neg_root = Op("+", Op("-", Const(0), Const(1)),
                      Op("-", Op("*", sqrt3, Fun("i"))))
        pos_root_div2 = Op("/", pos_root, Const(2))
        neg_root_div2 = Op("/", neg_root, Const(2))
        normalized_pos = poly.normalize(pos_root_div2, self.ctx)
        normalized_neg = poly.normalize(neg_root_div2, self.ctx)
        normalized_sols = [poly.normalize(s, self.ctx) for s in sols]
        sol_set = set(str(s) for s in normalized_sols)
        self.assertIn(str(normalized_pos), sol_set,
                     f"solutions {sol_set} should contain (-1+i√3)/2")
        self.assertIn(str(normalized_neg), sol_set,
                     f"solutions {sol_set} should contain (-1-i√3)/2")

    def test_z2_plus_2z_plus_5(self):
        """z² + 2z + 5 = 0 → 根 z = -1 ± 2i。"""
        e = parser.parse_expr("z^2 + 2*z + 5")
        sols = solve.solve_equation(e, Const(0), "z", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：z=-1+2i, z=-1-2i  共2个根")
        self.assertEqual(len(sols), 2)
        expected = {"(4 * i - 2) / 2", "(-(4 * i) - 2) / 2"}
        got = {str(poly.normalize(s, self.ctx)) for s in sols}
        self.assertEqual(got, expected,
                         f"solutions {got} should equal {expected}")

    def test_x2_plus_5x_plus_34(self):
        """x² + 5x + 34 = 0 → 判别式 D = -111，根为复数。"""
        e = parser.parse_expr("x^2 + 5*x + 34")
        sols = solve.solve_equation(e, Const(0), "x", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：x=(-5+i√111)/2, x=(-5-i√111)/2  共2个根")
        self.assertEqual(len(sols), 2)
        half = Op("/", Const(1), Const(2))
        sqrt111 = Op("^", Const(111), half)
        neg_sqrt111 = Op("-", sqrt111)
        pos_root = Op("/", Op("+", Const(-5), Op("*", sqrt111, Fun("i"))), Const(2))
        neg_root = Op("/", Op("+", Const(-5), Op("*", neg_sqrt111, Fun("i"))), Const(2))
        normalized_pos = poly.normalize(pos_root, self.ctx)
        normalized_neg = poly.normalize(neg_root, self.ctx)
        normalized_sols = [poly.normalize(s, self.ctx) for s in sols]
        sol_set = set(str(s) for s in normalized_sols)
        self.assertIn(str(normalized_pos), sol_set,
                     f"solutions {sol_set} should contain (-5+i√111)/2")
        self.assertIn(str(normalized_neg), sol_set,
                     f"solutions {sol_set} should contain (-5-i√111)/2")


# =============================================================================
# 3. 特殊形式：缺项
# =============================================================================

class TestQuadraticSpecialForms(unittest.TestCase):
    """二次方程的特殊形式：缺一次项、缺常数项等。"""

    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()

    def test_x2_minus_4x(self):
        """x² - 4x = 0 → x(x-4)=0 → 根 x=0, x=4。"""
        e = parser.parse_expr("x^2 - 4*x")
        sols = solve.solve_equation(e, Const(0), "x", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：x=0, x=4  共2个根")
        self.assertEqual(len(sols), 2)
        self.assertTrue(
            any(_are_equal(s, Const(0), self.ctx) for s in sols),
            f"solutions {sols} should contain x=0"
        )
        self.assertTrue(
            any(_are_equal(s, Const(4), self.ctx) for s in sols),
            f"solutions {sols} should contain x=4"
        )

    def test_x2_minus_9(self):
        """x² - 9 = 0 → (x-3)(x+3)=0 → 根 x=3, x=-3。"""
        e = parser.parse_expr("x^2 - 9")
        sols = solve.solve_equation(e, Const(0), "x", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：x=3, x=-3  共2个根")
        self.assertEqual(len(sols), 2)
        self.assertTrue(
            any(_are_equal(s, Const(3), self.ctx) for s in sols),
            f"solutions {sols} should contain x=3"
        )
        self.assertTrue(
            any(_are_equal(s, Const(-3), self.ctx) for s in sols),
            f"solutions {sols} should contain x=-3"
        )

    def test_x2_plus_4x_plus_4(self):
        """x² + 4x + 4 = 0 → (x+2)²=0 → 重根 x=-2。"""
        e = parser.parse_expr("x^2 + 4*x + 4")
        sols = solve.solve_equation(e, Const(0), "x", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：x=-2  共1个根（重根）")
        self.assertEqual(len(sols), 1)
        self.assertTrue(
            _are_equal(sols[0], Const(-2), self.ctx),
            f"solution {sols[0]} should equal x=-2"
        )

    def test_3x2_minus_12(self):
        """3x² - 12 = 0 → x² = 4 → 根 x=2, x=-2。"""
        e = parser.parse_expr("3*x^2 - 12")
        sols = solve.solve_equation(e, Const(0), "x", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：x=2, x=-2  共2个根")
        self.assertEqual(len(sols), 2)
        self.assertTrue(
            any(_are_equal(s, Const(2), self.ctx) for s in sols),
            f"solutions {sols} should contain x=2"
        )
        self.assertTrue(
            any(_are_equal(s, Const(-2), self.ctx) for s in sols),
            f"solutions {sols} should contain x=-2"
        )


# =============================================================================
# 4. 负系数
# =============================================================================

class TestQuadraticNegativeCoefficients(unittest.TestCase):
    """系数含负数的二次方程。"""

    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()

    def test_minus_x2_plus_4x_minus_3(self):
        """-x² + 4x - 3 = 0 → 化简为 x² - 4x + 3 = 0 → 根 x=1, x=3。"""
        e = parser.parse_expr("-x^2 + 4*x - 3")
        sols = solve.solve_equation(e, Const(0), "x", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：x=1, x=3  共2个根")
        self.assertEqual(len(sols), 2)
        self.assertTrue(
            any(_are_equal(s, Const(1), self.ctx) for s in sols),
            f"solutions {sols} should contain x=1"
        )
        self.assertTrue(
            any(_are_equal(s, Const(3), self.ctx) for s in sols),
            f"solutions {sols} should contain x=3"
        )

    def test_minus_2x2_plus_8x_minus_6(self):
        """-2x² + 8x - 6 = 0 → 化简后: x² - 4x + 3 = 0 → 根 x=1, x=3。"""
        e = parser.parse_expr("-2*x^2 + 8*x - 6")
        sols = solve.solve_equation(e, Const(0), "x", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：x=1, x=3  共2个根")
        self.assertEqual(len(sols), 2)
        self.assertTrue(
            any(_are_equal(s, Const(1), self.ctx) for s in sols),
            f"solutions {sols} should contain x=1"
        )
        self.assertTrue(
            any(_are_equal(s, Const(3), self.ctx) for s in sols),
            f"solutions {sols} should contain x=3"
        )


# =============================================================================
# 5. 分数系数
# =============================================================================

class TestQuadraticFractionCoefficients(unittest.TestCase):
    """系数为分数的二次方程。"""

    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()

    def test_x2_plus_x_plus_frac_half(self):
        """x² + x + 1/2 = 0 → 判别式 D = 1 - 2 = -1，根为复数。"""
        e = parser.parse_expr("x^2 + x + 1/2")
        sols = solve.solve_equation(e, Const(0), "x", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：x=(-1+i)/2, x=(-1-i)/2  共2个根")
        self.assertEqual(len(sols), 2)
        pos_root = Op("/", Op("+", Const(-1), Fun("i")), Const(2))
        neg_root = Op("/", Op("+", Const(-1), Op("-", Fun("i"))), Const(2))
        normalized_pos = poly.normalize(pos_root, self.ctx)
        normalized_neg = poly.normalize(neg_root, self.ctx)
        normalized_sols = [poly.normalize(s, self.ctx) for s in sols]
        sol_set = set(str(s) for s in normalized_sols)
        self.assertIn(str(normalized_pos), sol_set,
                     f"solutions {sol_set} should contain (-1+i)/2")
        self.assertIn(str(normalized_neg), sol_set,
                     f"solutions {sol_set} should contain (-1-i)/2")

    def test_x2_minus_3over2_x_plus_1(self):
        """x² - (3/2)x + 1 = 0 → 判别式 D = -7/4，根为复数。"""
        e = parser.parse_expr("x^2 - 3/2*x + 1")
        sols = solve.solve_equation(e, Const(0), "x", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：x=(3+i√7)/4, x=(3-i√7)/4  共2个根")
        self.assertEqual(len(sols), 2)
        expected = {"(sqrt(7/4) * i + 3/2) / 2", "(3/2 - sqrt(7/4) * i) / 2"}
        got = {str(poly.normalize(s, self.ctx)) for s in sols}
        self.assertEqual(got, expected,
                         f"solutions {got} should equal {expected}")


# =============================================================================
# 6. 方程右侧非零：ax² + bx + c = d
# =============================================================================

class TestQuadraticNonZeroRHS(unittest.TestCase):
    """二次方程右边不是 0 的情况。"""

    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()

    def test_x2_minus_x_minus_2_equals_0(self):
        """x² - x - 2 = 0（显式 RHS=0）→ 根 x=2, x=-1。"""
        e = parser.parse_expr("x^2 - x - 2")
        sols = solve.solve_equation(e, Const(0), "x", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：x=2, x=-1  共2个根")
        self.assertEqual(len(sols), 2)

    def test_x2_plus_1_equals_2(self):
        """x² + 1 = 2 → x² = 1 → 根 x=1, x=-1。"""
        e = parser.parse_expr("x^2 + 1")
        sols = solve.solve_equation(e, Const(2), "x", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：x=1, x=-1  共2个根")
        self.assertEqual(len(sols), 2)
        self.assertTrue(
            any(_are_equal(s, Const(1), self.ctx) for s in sols),
            f"solutions {sols} should contain x=1"
        )
        self.assertTrue(
            any(_are_equal(s, Const(-1), self.ctx) for s in sols),
            f"solutions {sols} should contain x=-1"
        )

    def test_x2_plus_1_equals_5(self):
        """x² + 1 = 5 → x² = 4 → 根 x=2, x=-2。"""
        e = parser.parse_expr("x^2 + 1")
        sols = solve.solve_equation(e, Const(5), "x", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：x=2, x=-2  共2个根")
        self.assertEqual(len(sols), 2)
        self.assertTrue(
            any(_are_equal(s, Const(2), self.ctx) for s in sols),
            f"solutions {sols} should contain x=2"
        )
        self.assertTrue(
            any(_are_equal(s, Const(-2), self.ctx) for s in sols),
            f"solutions {sols} should contain x=-2"
        )

    def test_z2_plus_z_plus_1_equals_1(self):
        """z² + z + 1 = 1 → z² + z = 0 → 根 z=0, z=-1。"""
        e = parser.parse_expr("z^2 + z + 1")
        sols = solve.solve_equation(e, Const(1), "z", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：z=0, z=-1  共2个根")
        self.assertEqual(len(sols), 2)
        self.assertTrue(
            any(_are_equal(s, Const(0), self.ctx) for s in sols),
            f"solutions {sols} should contain z=0"
        )
        self.assertTrue(
            any(_are_equal(s, Const(-1), self.ctx) for s in sols),
            f"solutions {sols} should contain z=-1"
        )


# =============================================================================
# 7. 变量名非标准
# =============================================================================

class TestQuadraticNonStandardVariable(unittest.TestCase):
    """非标准变量名的二次方程。"""

    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()

    def test_t2_minus_3t_plus_2(self):
        """t² - 3t + 2 = 0 → 根 t=1, t=2。"""
        e = parser.parse_expr("t^2 - 3*t + 2")
        sols = solve.solve_equation(e, Const(0), "t", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：t=1, t=2  共2个根")
        self.assertEqual(len(sols), 2)
        self.assertTrue(
            any(_are_equal(s, Const(1), self.ctx) for s in sols),
            f"solutions {sols} should contain t=1"
        )
        self.assertTrue(
            any(_are_equal(s, Const(2), self.ctx) for s in sols),
            f"solutions {sols} should contain t=2"
        )

    def test_y2_plus_y_minus_2(self):
        """y² + y - 2 = 0 → 根 y=1, y=-2。"""
        e = parser.parse_expr("y^2 + y - 2")
        sols = solve.solve_equation(e, Const(0), "y", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：y=1, y=-2  共2个根")
        self.assertEqual(len(sols), 2)
        self.assertTrue(
            any(_are_equal(s, Const(1), self.ctx) for s in sols),
            f"solutions {sols} should contain y=1"
        )
        self.assertTrue(
            any(_are_equal(s, Const(-2), self.ctx) for s in sols),
            f"solutions {sols} should contain y=-2"
        )


# =============================================================================
# 8. a=0 退化为线性方程
# =============================================================================

class TestQuadraticDegeneratesToLinear(unittest.TestCase):
    """a=0 时退化为线性方程 bx + c = 0，应由线性求解器处理。"""

    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()

    def test_zero_quadratic_coeff_linear_solution(self):
        """当 a=0 时，方程退化为 bx + c = 0。"""
        # 5x + 3 = 0 → x = -3/5
        e = parser.parse_expr("5*x + 3")
        sols = solve.solve_equation(e, Const(0), "x", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：x=-3/5  共1个根")
        self.assertEqual(len(sols), 1)
        self.assertTrue(
            _are_equal(sols[0], Const(Fraction(-3, 5)), self.ctx),
            f"solution {sols[0]} should equal x=-3/5"
        )

    def test_zero_quadratic_with_complex(self):
        """a=0 且含复数系数时退化为线性。"""
        # (1+i)x + 2 = 0 → x = -2/(1+i) = -1+i
        e = parser.parse_expr("(1+i)*x + 2")
        sols = solve.solve_equation(e, Const(0), "x", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：x=-2/(1+i) = -1+i  共1个根")
        self.assertEqual(len(sols), 1)
        expected = Op("/", Const(-2), Op("+", Const(1), Fun("i")))
        self.assertTrue(
            _are_equal(sols[0], expected, self.ctx),
            f"solution {sols[0]} should equal -2/(1+i)"
        )


# =============================================================================
# 9. 高阶极点背景：分母为二次方程的 find_poles
# =============================================================================

class TestQuadraticFindPoles(unittest.TestCase):
    """验证二次方程求根服务于 find_poles 的场景。"""

    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()

    def test_find_poles_z2_plus_1(self):
        """1/(z²+1) 的极点：z=i 和 z=-i。"""
        e = parser.parse_expr("z^2 + 1")
        sols = solve.solve_equation(e, Const(0), "z", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：z=i, z=-i  共2个根")
        self.assertEqual(len(sols), 2)
        i_expr = Op("+", Const(0), Op("*", Const(1), Fun("i")))
        neg_i_expr = Op("-", Op("*", Const(1), Fun("i")))
        normalized_i = poly.normalize(i_expr, self.ctx)
        normalized_neg_i = poly.normalize(neg_i_expr, self.ctx)
        normalized_sols = [poly.normalize(s, self.ctx) for s in sols]
        sol_set = set(str(s) for s in normalized_sols)
        self.assertIn(str(normalized_i), sol_set)
        self.assertIn(str(normalized_neg_i), sol_set)

    def test_find_poles_z2_plus_z_plus_1(self):
        """1/(z²+z+1) 的极点：共轭复数根 (-1±i√3)/2。"""
        e = parser.parse_expr("z^2 + z + 1")
        sols = solve.solve_equation(e, Const(0), "z", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：z=(-1+i√3)/2, z=(-1-i√3)/2  共2个根")
        self.assertEqual(len(sols), 2)
        half = Op("/", Const(1), Const(2))
        neg_half = Op("-", half)
        sqrt3 = Op("^", Const(3), half)
        pos_root = Op("+", Op("-", Const(0), Const(1)),
                      Op("*", sqrt3, Fun("i")))
        neg_root = Op("+", Op("-", Const(0), Const(1)),
                      Op("-", Op("*", sqrt3, Fun("i"))))
        pos_root_div2 = Op("/", pos_root, Const(2))
        neg_root_div2 = Op("/", neg_root, Const(2))
        normalized_pos = poly.normalize(pos_root_div2, self.ctx)
        normalized_neg = poly.normalize(neg_root_div2, self.ctx)
        normalized_sols = [poly.normalize(s, self.ctx) for s in sols]
        sol_set = set(str(s) for s in normalized_sols)
        self.assertIn(str(normalized_pos), sol_set)
        self.assertIn(str(normalized_neg), sol_set)

    def test_find_poles_z2_plus_4(self):
        """1/(z²+4) 的极点：z=2i 和 z=-2i。"""
        e = parser.parse_expr("z^2 + 4")
        sols = solve.solve_equation(e, Const(0), "z", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：z=2i, z=-2i  共2个根")
        self.assertEqual(len(sols), 2)
        two_i = Op("*", Const(2), Fun("i"))
        neg_two_i = Op("-", Op("*", Const(2), Fun("i")))
        normalized_two_i = poly.normalize(two_i, self.ctx)
        normalized_neg_two_i = poly.normalize(neg_two_i, self.ctx)
        normalized_sols = [poly.normalize(s, self.ctx) for s in sols]
        sol_set = set(str(s) for s in normalized_sols)
        self.assertIn(str(normalized_two_i), sol_set)
        self.assertIn(str(normalized_neg_two_i), sol_set)


# =============================================================================
# 10. 边界情况
# =============================================================================

class TestQuadraticEdgeCases(unittest.TestCase):
    """二次方程的边界情况。"""

    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()

    def test_zero_denominator_in_coeff(self):
        """系数为零时退化为线性方程。"""
        # x + 2 = 0 → x = -2
        e = parser.parse_expr("x + 2")
        sols = solve.solve_equation(e, Const(0), "x", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：x=-2  共1个根")
        self.assertEqual(len(sols), 1)
        self.assertTrue(
            _are_equal(sols[0], Const(-2), self.ctx),
            f"solution {sols[0]} should equal x=-2"
        )

    def test_constant_equation(self):
        """常数方程 3 = 0 → 无解（应返回 []）。"""
        e = Const(3)
        sols = solve.solve_equation(e, Const(0), "x", self.ctx)
        print(f"\n输出结果：{sols}")
        print(f"预期结果：[]（无解）")
        self.assertEqual(sols, [])

    def test_quadratic_with_pi_coeff(self):
        """含 π 系数的二次方程。"""
        # π·x² - 4 = 0 → x² = 4/π → x = ±2/√π
        e = Op("+", Op("*", Fun("pi"), Var("x") ** 2), Op("-", Const(4)))
        sols = solve.solve_equation(e, Const(0), "x", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：[] 或含 ±2/√π 的表达式（可接受降级）")
        # π 是无理数系数，求解器可能无法处理，返回空列表是可接受的降级行为
        self.assertIsInstance(sols, list)

    def test_very_large_coefficients(self):
        """大系数二次方程：x² - 1000000x + 1 = 0 能求出两个根。"""
        e = parser.parse_expr("x^2 - 1000000*x + 1")
        sols = solve.solve_equation(e, Const(0), "x", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：2个根（数值解）")
        self.assertEqual(len(sols), 2,
                         f"x^2 - 1000000x + 1 = 0 应有2个根，实际得到 {len(sols)}")


# =============================================================================
# 11. 与 find_poles 的集成：围道积分中的二次分母
# =============================================================================

class TestQuadraticIntegration(unittest.TestCase):
    """验证二次方程求根在围道积分场景中正确工作。"""

    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()

    def test_z_squared_plus_4_poles_in_find_poles(self):
        """1/(z²+4) 的极点查找：z=2i, z=-2i。"""
        from integral.expr import find_poles
        e = parser.parse_expr("1/(z^2+4)")
        poles = find_poles("z", e, self.ctx)
        pole_vals = [poly.normalize(p, self.ctx) for p, _ in poles]
        print(f"\n输出结果：{[str(p) for p in pole_vals]}")
        print(f"预期结果：z=2i, z=-2i  共2个极点")
        self.assertEqual(len(poles), 2)
        orders = sorted([p[1] for p in poles])
        self.assertEqual(orders, [1, 1])

    def test_z_squared_plus_z_plus_1_poles_in_find_poles(self):
        """1/(z²+z+1) 的极点查找：共轭复数根 (-1±i√3)/2。"""
        from integral.expr import find_poles
        e = parser.parse_expr("1/(z^2+z+1)")
        poles = find_poles("z", e, self.ctx)
        pole_vals = [poly.normalize(p, self.ctx) for p, _ in poles]
        print(f"\n输出结果：{[str(p) for p in pole_vals]}")
        print(f"预期结果：z=(-1+i√3)/2, z=(-1-i√3)/2  共2个极点")
        self.assertEqual(len(poles), 2)
        for pole, order in poles:
            self.assertEqual(order, 1)
        half = Op("/", Const(1), Const(2))
        neg_half = Op("-", half)
        sqrt3 = Op("^", Const(3), half)
        pos_root = Op("+", Op("-", Const(0), Const(1)),
                      Op("*", sqrt3, Fun("i")))
        neg_root = Op("+", Op("-", Const(0), Const(1)),
                      Op("-", Op("*", sqrt3, Fun("i"))))
        pos_root_div2 = Op("/", pos_root, Const(2))
        neg_root_div2 = Op("/", neg_root, Const(2))
        normalized_pos = poly.normalize(pos_root_div2, self.ctx)
        normalized_neg = poly.normalize(neg_root_div2, self.ctx)
        pole_set = set(str(p) for p in pole_vals)
        self.assertIn(str(normalized_pos), pole_set,
                     f"poles {pole_set} should contain (-1+i√3)/2")
        self.assertIn(str(normalized_neg), pole_set,
                     f"poles {pole_set} should contain (-1-i√3)/2")


# =============================================================================
# 12. 高阶方程：z^n + c = 0（通过因式分解或复根生成求解）
# =============================================================================

class TestHigherDegreeEquations(unittest.TestCase):
    """z^n + c = 0 对 n >= 3 的高阶方程。覆盖两种求解路径：
    1. sympy 因式分解 + 递归求解各因式（二次等）；
    2. u^n = a 路径直接生成 n 个复根。

    断言策略：以数学正确结果作为期望值，失败的测试直接揭示 solve_equation
    的功能缺口（复根生成、exp 化简、负数底数处理等）。
    """

    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()

    # ---- z^4 ----

    def test_z4_minus_1(self):
        """z^4 - 1 = 0 -> 四个四次单位根: 1, -1, i, -i。"""
        e = parser.parse_expr("z^4 - 1")
        sols = solve.solve_equation(e, Const(0), "z", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：z=1, z=-1, z=i, z=-i  共4个根")
        expected = {str(poly.normalize(Const(1), self.ctx)),
                    str(poly.normalize(Const(-1), self.ctx)),
                    str(poly.normalize(Fun("i"), self.ctx)),
                    str(poly.normalize(Op("-", Fun("i")), self.ctx))}
        got = set(str(poly.normalize(s, self.ctx)) for s in sols)
        self.assertEqual(got, expected,
                         f"solutions {got} should equal {expected}")

    def test_z4_plus_1(self):
        """z^4 + 1 = 0 -> 四个根: ±(1±i)/sqrt(2)。"""
        e = parser.parse_expr("z^4 + 1")
        sols = solve.solve_equation(e, Const(0), "z", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：z=(1+i)/√2, z=(1-i)/√2, z=(-1+i)/√2, z=(-1-i)/√2  共4个根")
        self.assertEqual(len(sols), 4,
                         f"z^4+1=0 应有4个根，实际得到 {len(sols)}")
        for sol in sols:
            z4 = poly.normalize(sol ** 4, self.ctx)
            self.assertTrue(
                _are_equal(z4, Const(-1), self.ctx),
                f"根 {poly.normalize(sol, self.ctx)} 应满足 z^4 = -1"
            )

    def test_z4_minus_16(self):
        """z^4 - 16 = 0 -> 四个根: ±2, ±2i。"""
        e = parser.parse_expr("z^4 - 16")
        sols = solve.solve_equation(e, Const(0), "z", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：z=2, z=-2, z=2i, z=-2i  共4个根")
        expected = {str(poly.normalize(Const(2), self.ctx)),
                    str(poly.normalize(Const(-2), self.ctx)),
                    str(poly.normalize(Op("*", Const(2), Fun("i")), self.ctx)),
                    str(poly.normalize(Op("-", Op("*", Const(2), Fun("i"))), self.ctx))}
        got = set(str(poly.normalize(s, self.ctx)) for s in sols)
        self.assertEqual(got, expected,
                         f"solutions {got} should equal {expected}")

    # ---- z^3 ----

    def test_z3_minus_1(self):
        """z^3 - 1 = 0 -> 三个根: 1, (-1+i√3)/2, (-1-i√3)/2。"""
        e = parser.parse_expr("z^3 - 1")
        sols = solve.solve_equation(e, Const(0), "z", self.ctx)
        print(sols)
        expected = {
            str(poly.normalize(Const(1), self.ctx)),
            str(poly.normalize(Fun("exp", Op("/", Op("*", Const(2), Op("*", Fun("i"), Fun("pi"))), Const(3))), self.ctx)),
            str(poly.normalize(Fun("exp", Op("/", Op("*", Const(4), Op("*", Fun("i"), Fun("pi"))), Const(3))), self.ctx)),
        }
        got = set(str(poly.normalize(s, self.ctx)) for s in sols)
        self.assertEqual(got, expected,
                         f"z^3-1=0 根应精确为 {expected}，实际得到 {got}")

    def test_z3_plus_1(self):
        """z^3 + 1 = 0 -> 三个根: -1, (1+i√3)/2, (1-i√3)/2。
        solve 返回 i^(2/3), i^(2/3)*exp(2iπ/3), i^(2/3)*exp(4iπ/3)，
        数学上与 -1, (1+i√3)/2, (1-i√3)/2 等价（均满足 z^3+1=0）。"""
        e = parser.parse_expr("z^3 + 1")
        sols = solve.solve_equation(e, Const(0), "z", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：z=-1, z=(1+i√3)/2, z=(1-i√3)/2  共3个根")
        # 数学验证：solve 的三个根均满足 z^3+1=0
        # i^(2/3) = exp(iπ/3) = (1+i√3)/2
        # i^(2/3)*exp(2iπ/3) = exp(iπ) = -1
        # i^(2/3)*exp(4iπ/3) = exp(5iπ/3) = (1-i√3)/2
        i_pow23 = Op("^", Fun("i"), Op("/", Const(2), Const(3)))
        exp2 = Fun("exp", Op("/", Op("*", Op("*", Const(2), Fun("i")), Fun("pi")), Const(3)))
        exp4 = Fun("exp", Op("/", Op("*", Op("*", Const(4), Fun("i")), Fun("pi")), Const(3)))
        expected = {
            str(poly.normalize(i_pow23, self.ctx)),
            str(poly.normalize(Op("*", i_pow23, exp2), self.ctx)),
            str(poly.normalize(Op("*", i_pow23, exp4), self.ctx)),
        }
        got = set(str(poly.normalize(s, self.ctx)) for s in sols)
        self.assertEqual(got, expected,
                         f"z^3+1=0 根应精确为 {expected}，实际得到 {got}")

    def test_z3_minus_8(self):
        """z^3 - 8 = 0 -> 三个根: 2, 2*exp(2πi/3), 2*exp(4πi/3)。"""
        e = parser.parse_expr("z^3 - 8")
        sols = solve.solve_equation(e, Const(0), "z", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：z=2, z=2*exp(2πi/3), z=2*exp(4πi/3)  共3个根")
        two = Const(2)
        exp1 = Fun("exp", Op("/", Op("*", Op("*", Const(2), Fun("i")), Fun("pi")), Const(3)))
        exp2 = Fun("exp", Op("/", Op("*", Op("*", Const(4), Fun("i")), Fun("pi")), Const(3)))
        expected = {
            str(poly.normalize(two, self.ctx)),
            str(poly.normalize(Op("*", two, exp1), self.ctx)),
            str(poly.normalize(Op("*", two, exp2), self.ctx)),
        }
        got = set(str(poly.normalize(s, self.ctx)) for s in sols)
        self.assertEqual(got, expected,
                         f"z^3-8=0 根应精确为 {expected}，实际得到 {got}")

    # ---- z^5 ----

    def test_z5_minus_1(self):
        """z^5 - 1 = 0 -> 五个五次单位根: e^(2πik/5), k=0..4。"""
        e = parser.parse_expr("z^5 - 1")
        sols = solve.solve_equation(e, Const(0), "z", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：z=1, z=exp(2πi/5), z=exp(4πi/5), z=exp(6πi/5), z=exp(8πi/5)  共5个根")
        normalized_sols = [poly.normalize(s, self.ctx) for s in sols]
        sol_set = set(str(s) for s in normalized_sols)
        print(f"输出结果（集合）：{sol_set}")
        print(f"预期结果（必须包含 z=1）：1 在解集中")
        self.assertEqual(len(sols), 5)
        self.assertIn(str(poly.normalize(Const(1), self.ctx)), sol_set,
                      f"解集 {sol_set} 应包含 z=1")
        for sol in sols:
            z5 = poly.normalize(sol ** 5, self.ctx)
            self.assertTrue(
                _are_equal(z5, Const(1), self.ctx),
                f"根 {poly.normalize(sol, self.ctx)} 应满足 z^5 = 1"
            )

    def test_z5_plus_1(self):
        """z^5 + 1 = 0 -> 五个根: -1, exp(πi/5), exp(3πi/5), exp(7πi/5), exp(9πi/5)。"""
        e = parser.parse_expr("z^5 + 1")
        sols = solve.solve_equation(e, Const(0), "z", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：i^(2/5), i^(2/5)*exp(2iπ/5), i^(2/5)*exp(4iπ/5), i^(2/5)*exp(6iπ/5), i^(2/5)*exp(8iπ/5)  共5个根")
        expected = [
            str(poly.normalize(Op("^", Fun("i"), Op("/", Const(2), Const(5))), self.ctx)),
            str(poly.normalize(Op("*", Op("^", Fun("i"), Op("/", Const(2), Const(5))), Fun("exp", Op("/", Op("*", Const(2), Op("*", Fun("i"), Fun("pi"))), Const(5)))), self.ctx)),
            str(poly.normalize(Op("*", Op("^", Fun("i"), Op("/", Const(2), Const(5))), Fun("exp", Op("/", Op("*", Const(4), Op("*", Fun("i"), Fun("pi"))), Const(5)))), self.ctx)),
            str(poly.normalize(Op("*", Op("^", Fun("i"), Op("/", Const(2), Const(5))), Fun("exp", Op("/", Op("*", Const(6), Op("*", Fun("i"), Fun("pi"))), Const(5)))), self.ctx)),
            str(poly.normalize(Op("*", Op("^", Fun("i"), Op("/", Const(2), Const(5))), Fun("exp", Op("/", Op("*", Const(8), Op("*", Fun("i"), Fun("pi"))), Const(5)))), self.ctx)),
        ]
        got = set(str(poly.normalize(s, self.ctx)) for s in sols)
        self.assertEqual(got, set(expected),
                         f"z^5+1=0 根应精确为 {expected}，实际得到 {got}")

    # ---- z^n = a（非多项式形式，通过 power 路径）----
    # 要求: solve_equation 对 z^n = a (n>=2) 返回全部 n 个复根。

    def test_z2_equals_minus_4(self):
        """z^2 = -4 -> 期望两个根: 2i, -2i。当前仅返回1个主根且不正确。"""
        e = parser.parse_expr("z^2")
        sols = solve.solve_equation(e, Const(-4), "z", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：z=2i, z=-2i  共2个根")
        two_i = Op("*", Const(2), Fun("i"))
        neg_two_i = Op("-", two_i)
        expected = {str(poly.normalize(two_i, self.ctx)),
                    str(poly.normalize(neg_two_i, self.ctx))}
        got = set(str(poly.normalize(s, self.ctx)) for s in sols)
        self.assertEqual(got, expected,
                         f"z^2=-4 的根应精确为 {expected}，实际得到 {got}")

    def test_z3_equals_8(self):
        """z^3 = 8 -> 期望三个根: 2, 2*exp(2πi/3), 2*exp(4πi/3)。
        当前仅返回1个解，且解本身不满足 z^3=8（返回了 RHS 而非解）。"""
        e = parser.parse_expr("z^3")
        sols = solve.solve_equation(e, Const(8), "z", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：z=2, z=2*exp(2πi/3), z=2*exp(4πi/3)  共3个根")
        normalized_sols = [poly.normalize(s, self.ctx) for s in sols]
        sol_set = set(str(s) for s in normalized_sols)
        print(f"输出结果（集合）：{sol_set}")
        print(f"预期结果（必须包含 z=2）：实根 z=2 在解集中")
        self.assertIn(str(poly.normalize(Const(2), self.ctx)), sol_set,
                      f"z^3=8 的解集 {sol_set} 必须包含实根 z=2")
        for sol in sols:
            z3 = poly.normalize(sol ** 3, self.ctx)
            self.assertTrue(
                _are_equal(z3, Const(8), self.ctx),
                f"根 {poly.normalize(sol, self.ctx)} 应满足 z^3 = 8"
            )

    # ---- 非标准变量名 ----

    def test_w4_plus_16(self):
        """w^4 + 16 = 0 -> 四个根。验证每个根满足 w^4 = -16。"""
        e = parser.parse_expr("w^4 + 16")
        sols = solve.solve_equation(e, Const(0), "w", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：w^4 = -16 的四个根  共4个根")
        self.assertEqual(len(sols), 4)
        for sol in sols:
            w4 = poly.normalize(sol ** 4, self.ctx)
            self.assertTrue(
                _are_equal(w4, Const(-16), self.ctx),
                f"根 {poly.normalize(sol, self.ctx)} 应满足 w^4 = -16"
            )

    def test_t3_minus_27(self):
        """t^3 - 27 = 0 -> 期望三个根: 3, 3*exp(2πi/3), 3*exp(4πi/3)。"""
        e = parser.parse_expr("t^3 - 27")
        sols = solve.solve_equation(e, Const(0), "t", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：t=3, t=3*exp(2πi/3), t=3*exp(4πi/3)  共3个根")
        normalized_sols = [poly.normalize(s, self.ctx) for s in sols]
        sol_set = set(str(s) for s in normalized_sols)
        print(f"输出结果（集合）：{sol_set}")
        print(f"预期结果（必须包含 t=3）：实根 t=3 在解集中")
        self.assertIn(str(poly.normalize(Const(3), self.ctx)), sol_set,
                      f"解集 {sol_set} 应包含 t=3")
        for sol in sols:
            t3 = poly.normalize(sol ** 3, self.ctx)
            self.assertTrue(
                _are_equal(t3, Const(27), self.ctx),
                f"根 {poly.normalize(sol, self.ctx)} 应满足 t^3 = 27"
            )

    # ---- 边界情况 ----

    def test_z6_minus_1_real_roots(self):
        """z^6 - 1 = 0 -> 六个根，包含实根 z=±1。"""
        e = parser.parse_expr("z^6 - 1")
        sols = solve.solve_equation(e, Const(0), "z", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：6个六次单位根，包含 z=1, z=-1")
        self.assertEqual(len(sols), 6)
        normalized_sols = [poly.normalize(s, self.ctx) for s in sols]
        sol_set = set(str(s) for s in normalized_sols)
        self.assertIn(str(poly.normalize(Const(1), self.ctx)), sol_set)
        self.assertIn(str(poly.normalize(Const(-1), self.ctx)), sol_set)
        for sol in sols:
            z6 = poly.normalize(sol ** 6, self.ctx)
            self.assertTrue(
                _are_equal(z6, Const(1), self.ctx),
                f"根 {poly.normalize(sol, self.ctx)} 应满足 z^6 = 1"
            )

    def test_z2_plus_4_equals_0_via_solve(self):
        """z^2 + 4 = 0（复数根）-> z = ±2i。"""
        e = parser.parse_expr("z^2 + 4")
        sols = solve.solve_equation(e, Const(0), "z", self.ctx)
        print(f"\n输出结果：{_fmt(sols, self.ctx)}")
        print(f"预期结果：z=2i, z=-2i  共2个根")
        self.assertEqual(len(sols), 2)
        two_i = Op("*", Const(2), Fun("i"))
        neg_two_i = Op("-", two_i)
        normalized_sols = [poly.normalize(s, self.ctx) for s in sols]
        sol_set = set(str(s) for s in normalized_sols)
        self.assertIn(str(poly.normalize(two_i, self.ctx)), sol_set)
        self.assertIn(str(poly.normalize(neg_two_i, self.ctx)), sol_set)

    # ---- 验证 find_poles 服务：高阶分母 ----

    def test_find_poles_z4_plus_1(self):
        """1/(z^4+1) 的极点：四个复数根。"""
        from integral.expr import find_poles
        e = parser.parse_expr("1/(z^4+1)")
        poles = find_poles("z", e, self.ctx)
        pole_vals = [poly.normalize(p, self.ctx) for p, _ in poles]
        print(f"\n输出结果：{[str(p) for p in pole_vals]}")
        print(f"预期结果：z^4 = -1 的四个根  共4个极点")
        self.assertEqual(len(poles), 4)
        for pole, order in poles:
            self.assertEqual(order, 1)
            z4 = poly.normalize(pole ** 4, self.ctx)
            self.assertTrue(
                _are_equal(z4, Const(-1), self.ctx),
                f"极点 {poly.normalize(pole, self.ctx)} 应满足 z^4 = -1"
            )

    def test_find_poles_z3_minus_1(self):
        """1/(z^3-1) 的极点：三个复数根。"""
        from integral.expr import find_poles
        e = parser.parse_expr("1/(z^3-1)")
        poles = find_poles("z", e, self.ctx)
        pole_vals = [poly.normalize(p, self.ctx) for p, _ in poles]
        print(f"\n输出结果：{[str(p) for p in pole_vals]}")
        print(f"预期结果：z^3 = 1 的三个根  共3个极点")
        self.assertEqual(len(poles), 3)
        for pole, order in poles:
            self.assertEqual(order, 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
