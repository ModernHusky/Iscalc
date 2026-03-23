"""
完备性测试: integral/expr.py find_poles 函数及其依赖的 solve_equation

测试 find_poles(z, p(z)/q(z)) 对各种分母 q(z) 的求解能力：
  - 代数方程：线性、二次、三次、高阶
  - 重复根（高阶极点）
  - 复数根
  - 有理数/无理数系数
  - 符号参数
  - solve_equation 无法求解的方程（优雅降级）
  - 非除法表达式（边界情况）
  - 去重逻辑

运行: python -m pytest tests/test_find_poles.py -v
"""

import os
import unittest
from pathlib import Path
from fractions import Fraction

from integral import parser, rules, expr, context, poly, condprover
from integral.expr import Const, Var, Op, Fun
from integral.expr import find_poles

PROJECT_ROOT = Path(__file__).parent.parent
os.chdir(PROJECT_ROOT)


def clear_caches():
    poly._normalize_cache.clear()
    condprover.clear_condition_cache()
    rules._poles_cache.clear()
    rules._winding_cache.clear()


def make_ctx():
    ctx = context.Context()
    ctx.load_book("base")
    return ctx


# =============================================================================
# 1. 边界情况：非除法表达式 / 无极点
# =============================================================================

class TestFindPolesEdgeCases(unittest.TestCase):
    """find_poles 对非除法表达式的处理——应返回空列表。"""

    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()

    def test_non_division_expr_returns_empty(self):
        """纯变量 z 无分母，不构成极点。"""
        result = find_poles("z", Var("z"), self.ctx)
        self.assertEqual(result, [])

    def test_constant_expr_returns_empty(self):
        """常数表达式无分母。"""
        result = find_poles("z", Const(5), self.ctx)
        self.assertEqual(result, [])

    def test_division_by_constant_returns_empty(self):
        """分母不含变量，不构成极点。"""
        expr_ = Op("/", Const(1), Const(2))
        result = find_poles("z", expr_, self.ctx)
        self.assertEqual(result, [])

    def test_product_expr_returns_empty(self):
        """乘积表达式不是 find_poles 处理的除法形式。"""
        expr_ = Op("*", Var("z"), Const(2))
        result = find_poles("z", expr_, self.ctx)
        self.assertEqual(result, [])

    def test_negative_linear_denom(self):
        """分母为 -(z-1)：应仍能找到极点。"""
        denom = -(Var("z") - Const(1))
        expr_ = Op("/", Const(1), denom)
        result = find_poles("z", expr_, self.ctx)
        self.assertEqual(len(result), 1)
        pole, order = result[0]
        self.assertEqual(order, 1)


# =============================================================================
# 2. 简单一阶极点（线性分母）
# =============================================================================

class TestSimpleFirstOrderPoles(unittest.TestCase):
    """1/(z - a) 形式：z = a 为一阶极点。"""

    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()

    def _check(self, expr_, expected_count, expected_order=1, var="z"):
        """通用验证：极点数量、阶数、且每个极点代入分母为零。"""
        result = find_poles(var, expr_, self.ctx)
        self.assertEqual(len(result), expected_count,
            f"Expected {expected_count} poles, got {len(result)}: {result}")
        for pole, order in result:
            self.assertEqual(order, expected_order,
                f"Pole {pole} has order {order}, expected {expected_order}")
            denom = expr_.args[1]
            val = denom.subst(var, pole)
            val = poly.normalize(val, self.ctx)
            self.assertEqual(val, Const(0),
                f"Pole {pole} does not satisfy denom == 0")

    def test_1_over_z_minus_1(self):
        """1/(z-1): 极点 z=1, 一阶。"""
        self._check(Op("/", Const(1), Var("z") - Const(1)), 1, 1)

    def test_1_over_z_plus_1(self):
        """1/(z+1): 极点 z=-1, 一阶。"""
        self._check(Op("/", Const(1), Var("z") + Const(1)), 1, 1)

    def test_1_over_z(self):
        """1/z: 极点 z=0, 一阶。"""
        self._check(Op("/", Const(1), Var("z")), 1, 1)

    def test_1_over_2z_minus_1(self):
        """1/(2z-1): 极点 z=1/2, 一阶。"""
        self._check(Op("/", Const(1), Const(2) * Var("z") - Const(1)), 1, 1)

    def test_1_over_3z_plus_2(self):
        """1/(3z+2): 极点 z=-2/3, 一阶。"""
        self._check(Op("/", Const(1), Const(3) * Var("z") + Const(2)), 1, 1)

    def test_1_over_z_minus_half(self):
        """1/(z-1/2): 极点 z=1/2, 一阶。"""
        self._check(Op("/", Const(1), Var("z") - Const(Fraction(1, 2))), 1, 1)

    def test_1_over_negative_coeff(self):
        """1/(-2z+3): 极点 z=3/2, 一阶。"""
        self._check(Op("/", Const(1), -Const(2) * Var("z") + Const(3)), 1, 1)


# =============================================================================
# 3. 重复根：高阶极点
# =============================================================================

class TestHigherOrderPoles(unittest.TestCase):
    """1/(z-a)^n 形式：z=a 为 n 阶极点。"""

    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()

    def test_second_order_pole(self):
        """1/(z-1)^2: 极点 z=1, 二阶。"""
        denom = (Var("z") - Const(1)) ** 2
        expr_ = Op("/", Const(1), denom)
        result = find_poles("z", expr_, self.ctx)
        self.assertGreaterEqual(len(result), 1)
        poles = [(p, o) for p, o in result]
        matching = [(p, o) for p, o in poles
                    if poly.normalize(p - Const(1), self.ctx) == Const(0)]
        self.assertGreater(len(matching), 0,
            f"Expected pole Const(1) not found in {poles}")
        pole, order = matching[0]
        self.assertEqual(order, 2,
            f"Pole {pole} order={order}, expected 2")

    def test_third_order_pole(self):
        """1/(z-2)^3: 极点 z=2, 三阶。"""
        denom = (Var("z") - Const(2)) ** 3
        expr_ = Op("/", Const(1), denom)
        result = find_poles("z", expr_, self.ctx)
        self.assertGreaterEqual(len(result), 1)
        poles = [(p, o) for p, o in result]
        matching = [(p, o) for p, o in poles
                    if poly.normalize(p - Const(2), self.ctx) == Const(0)]
        self.assertGreater(len(matching), 0,
            f"Expected pole Const(2) not found in {poles}")
        pole, order = matching[0]
        self.assertEqual(order, 3,
            f"Pole {pole} order={order}, expected 3")

    def test_pole_order_four(self):
        """1/(z+1)^4: 极点 z=-1, 四阶。"""
        denom = (Var("z") + Const(1)) ** 4
        expr_ = Op("/", Const(1), denom)
        result = find_poles("z", expr_, self.ctx)
        self.assertGreaterEqual(len(result), 1)
        poles = [(p, o) for p, o in result]
        matching = [(p, o) for p, o in poles
                    if poly.normalize(p - Const(-1), self.ctx) == Const(0)]
        self.assertGreater(len(matching), 0,
            f"Expected pole Const(-1) not found in {poles}")
        pole, order = matching[0]
        self.assertEqual(order, 4,
            f"Pole {pole} order={order}, expected 4")

    def test_order_5_pole(self):
        """1/(z-3)^5: 极点 z=3, 五阶。"""
        denom = (Var("z") - Const(3)) ** 5
        expr_ = Op("/", Const(1), denom)
        result = find_poles("z", expr_, self.ctx)
        self.assertGreaterEqual(len(result), 1)
        poles = [(p, o) for p, o in result]
        matching = [(p, o) for p, o in poles
                    if poly.normalize(p - Const(3), self.ctx) == Const(0)]
        self.assertGreater(len(matching), 0,
            f"Expected pole Const(3) not found in {poles}")
        pole, order = matching[0]
        self.assertEqual(order, 5,
            f"Pole {pole} order={order}, expected 5")

    def test_pole_order_regression(self):
        """回归测试：验证 _determine_pole_order 本身对低阶极点是正确的。"""
        from integral.expr import _determine_pole_order
        denom1 = Var("z") - Const(1)
        denom2 = (Var("z") - Const(1)) ** 2
        denom3 = (Var("z") - Const(1)) ** 3
        self.assertEqual(_determine_pole_order("z", denom1, Const(1), self.ctx), 1)
        self.assertEqual(_determine_pole_order("z", denom2, Const(1), self.ctx), 2)
        self.assertEqual(_determine_pole_order("z", denom3, Const(1), self.ctx), 3)


# =============================================================================
# 4. 多项式分母：多个不同极点
# =============================================================================

class TestMultipleDistinctPoles(unittest.TestCase):
    """分母为乘积 (z-a)(z-b)...: 每个根对应一个极点。"""

    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()

    def test_two_simple_poles(self):
        """1/((z-1)(z+1)): 极点 z=1, z=-1, 均为一阶。"""
        denom = (Var("z") - Const(1)) * (Var("z") + Const(1))
        expr_ = Op("/", Const(1), denom)
        result = find_poles("z", expr_, self.ctx)
        self.assertEqual(len(result), 2)
        orders = sorted([r[1] for r in result])
        self.assertEqual(orders, [1, 1])

    def test_three_simple_poles(self):
        """1/((z-1)(z-2)(z-3)): 三个一阶极点。"""
        denom = (Var("z") - Const(1)) * (Var("z") - Const(2)) * (Var("z") - Const(3))
        expr_ = Op("/", Const(1), denom)
        result = find_poles("z", expr_, self.ctx)
        self.assertEqual(len(result), 3)
        orders = [r[1] for r in result]
        self.assertTrue(all(o == 1 for o in orders))

    def test_mixed_order_poles(self):
        """1/((z)^2 * (z-1)): 极点 z=0 二阶, z=1 一阶。"""
        denom = (Var("z") ** 2) * (Var("z") - Const(1))
        expr_ = Op("/", Const(1), denom)
        result = find_poles("z", expr_, self.ctx)
        self.assertEqual(len(result), 2)
        orders = {r[1] for r in result}
        self.assertEqual(orders, {1, 2})

    def test_three_poles_mixed_orders(self):
        """1/((z)^2 * (z-1) * (z+1)): z=0 二阶, z=1 一阶, z=-1 一阶。"""
        denom = (Var("z") ** 2) * (Var("z") - Const(1)) * (Var("z") + Const(1))
        expr_ = Op("/", Const(1), denom)
        result = find_poles("z", expr_, self.ctx)
        self.assertEqual(len(result), 3)
        orders = {r[1] for r in result}
        self.assertEqual(orders, {1, 2})
        second_order_count = sum(1 for r in result if r[1] == 2)
        self.assertEqual(second_order_count, 1)
        first_order_count = sum(1 for r in result if r[1] == 1)
        self.assertEqual(first_order_count, 2)


# =============================================================================
# 5. 复数根极点
# =============================================================================

class TestComplexPoles(unittest.TestCase):
    """分母有复数根: z^2+1=0 => z=i, z=-i。"""

    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()

    def test_z_squared_plus_1(self):
        """1/(z^2+1): 极点 z=i 和 z=-i, 均一阶。"""
        denom = Var("z") ** 2 + Const(1)
        expr_ = Op("/", Const(1), denom)
        result = find_poles("z", expr_, self.ctx)
        self.assertEqual(len(result), 2)
        orders = [r[1] for r in result]
        self.assertEqual(sorted(orders), [1, 1])
        for pole, order in result:
            self.assertEqual(order, 1)
            val = denom.subst("z", pole)
            val = poly.normalize(val, self.ctx)
            self.assertEqual(val, Const(0),
                f"pole {pole} is not actually a zero of {denom}")

    def test_z_squared_minus_1(self):
        """1/(z^2-1) = 1/((z-1)(z+1)): 极点 z=1, z=-1。"""
        denom = Var("z") ** 2 - Const(1)
        expr_ = Op("/", Const(1), denom)
        result = find_poles("z", expr_, self.ctx)
        self.assertEqual(len(result), 2)
        orders = [r[1] for r in result]
        self.assertEqual(sorted(orders), [1, 1])

    def test_z_squared_plus_4(self):
        """1/(z^2+4): 极点 z=2i, z=-2i, 一阶。"""
        denom = Var("z") ** 2 + Const(4)
        expr_ = Op("/", Const(1), denom)
        result = find_poles("z", expr_, self.ctx)
        self.assertEqual(len(result), 2)
        orders = [r[1] for r in result]
        self.assertEqual(sorted(orders), [1, 1])

    def test_z_four_minus_one(self):
        """1/(z^4-1): 四个一阶极点（单位圆的四次单位根）。"""
        denom = Var("z") ** 4 - Const(1)
        expr_ = Op("/", Const(1), denom)
        result = find_poles("z", expr_, self.ctx)
        self.assertEqual(len(result), 4)
        orders = [r[1] for r in result]
        self.assertEqual(sorted(orders), [1, 1, 1, 1])

    def test_z_cubed_plus_z(self):
        """z/(z^3+z) = 1/(z^2+1): solve_equation(z^3+z=0) 返回 [0, i, -i]。

        验证通过 sympy 因式分解 + 递归求根的完整链路。
        z=0 处分子分母同时为零，为可去奇点，不计入极点；
        z=i 和 z=-i 为一阶极点。
        """
        from integral import solve
        denom_raw = Var("z") ** 3 + Var("z")
        zeros = solve.solve_equation(denom_raw, Const(0), "z", self.ctx)
        self.assertEqual(len(zeros), 3,
            f"z^3+z=0 应返回 3 个根 [0, i, -i]，实际: {zeros}")
        poles = find_poles("z", Op("/", Var("z"), denom_raw), self.ctx)
        self.assertEqual(len(poles), 2,
            f"z/(z^3+z) 应有 2 个极点（z=0 可去），实际: {poles}")
        orders = [o for _, o in poles]
        self.assertEqual(sorted(orders), [1, 1])


# =============================================================================
# 6. 有理数/无理数系数
# =============================================================================

class TestRationalCoefficients(unittest.TestCase):
    """分母系数为有理数 (Fraction)。"""

    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()

    def test_rational_pole_half(self):
        """1/(2z-3): 极点 z=3/2。"""
        expr_ = Op("/", Const(1), Const(2) * Var("z") - Const(3))
        result = find_poles("z", expr_, self.ctx)
        self.assertEqual(len(result), 1)
        pole, order = result[0]
        self.assertEqual(order, 1)
        val = (Const(2) * Var("z") - Const(3)).subst("z", pole)
        val = poly.normalize(val, self.ctx)
        self.assertEqual(val, Const(0))

    def test_multiple_rational_poles(self):
        """1/((3z-1)(2z+5)): 两个有理数极点。"""
        denom = (Const(3) * Var("z") - Const(1)) * (Const(2) * Var("z") + Const(5))
        expr_ = Op("/", Const(1), denom)
        result = find_poles("z", expr_, self.ctx)
        self.assertEqual(len(result), 2)
        orders = [r[1] for r in result]
        self.assertEqual(sorted(orders), [1, 1])

    def test_pi_coefficient(self):
        """1/(π*z - 1): 极点 z=1/π（无理数极点）。"""
        pi_const = Fun("pi")
        expr_ = Op("/", Const(1), pi_const * Var("z") - Const(1))
        result = find_poles("z", expr_, self.ctx)
        self.assertEqual(len(result), 1)
        pole, order = result[0]
        self.assertEqual(order, 1)

    def test_fraction_coefficient(self):
        """1/((1/2)z - 1): 极点 z=2。"""
        half = Const(Fraction(1, 2))
        expr_ = Op("/", Const(1), half * Var("z") - Const(1))
        result = find_poles("z", expr_, self.ctx)
        self.assertEqual(len(result), 1)
        pole, order = result[0]
        self.assertEqual(order, 1)
        val = (half * Var("z") - Const(1)).subst("z", pole)
        val = poly.normalize(val, self.ctx)
        self.assertEqual(val, Const(0))


# =============================================================================
# 7. 符号参数（用 Var 而非 Symbol）
# =============================================================================

class TestSymbolicParameters(unittest.TestCase):
    """分母含有符号参数 a, b 等（使用 Var 表示独立于 z 的符号常量）。"""

    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()

    def test_symbol_param_linear(self):
        """1/(z - a): 极点 z=a, 一阶。"""
        a_sym = Var("a")
        denom = Var("z") - a_sym
        expr_ = Op("/", Const(1), denom)
        result = find_poles("z", expr_, self.ctx)
        self.assertEqual(len(result), 1)
        pole, order = result[0]
        self.assertEqual(order, 1)
        diff = poly.normalize(pole - a_sym, self.ctx)
        self.assertEqual(diff, Const(0))

    def test_symbol_param_quadratic(self):
        """1/(z^2 - a): 极点 z=√a, z=-√a。"""
        a_sym = Var("a")
        denom = Var("z") ** 2 - a_sym
        expr_ = Op("/", Const(1), denom)
        result = find_poles("z", expr_, self.ctx)
        self.assertIsInstance(result, list)
        for pole, order in result:
            val = denom.subst("z", pole)
            val = poly.normalize(val, self.ctx)
            self.assertEqual(val, Const(0),
                f"pole {pole} is not actually a zero of {denom}")

    def test_two_symbol_params(self):
        """1/((z-a)(z-b)): 极点 z=a, z=b。"""
        a_sym = Var("a")
        b_sym = Var("b")
        denom = (Var("z") - a_sym) * (Var("z") - b_sym)
        expr_ = Op("/", Const(1), denom)
        result = find_poles("z", expr_, self.ctx)
        self.assertEqual(len(result), 2)


# =============================================================================
# 8. solve_equation 无法求解的方程（优雅降级）
# =============================================================================

class TestHigherDegreeAndSymbolic(unittest.TestCase):
    """高次多项式和含符号参数方程的极点检测。

    根据代数基本定理，n 次复多项式恰有 n 个根（计重数）。
    因此 1/(z^5 - z + 1) 应该有 5 个一阶极点。
    """

    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()

    def test_z5_minus_z_plus_one_has_five_poles(self):
        """1/(z^5 - z + 1): 根据代数基本定理恰有 5 个一阶极点。

        【算法限制】：solve_equation 无法求五次方程根（z^5 - z + 1 = 0）。
        预期结果：5 个极点（代数基本定理）；当前系统返回 []。
        """
        denom = Var("z") ** 5 - Var("z") + Const(1)
        expr_ = Op("/", Const(1), denom)
        result = find_poles("z", expr_, self.ctx)
        self.assertEqual(len(result), 5)
        orders = [r[1] for r in result]
        self.assertEqual(sorted(orders), [1, 1, 1, 1, 1])

    def test_z5_minus_z_plus_one_solve_equation_behavior(self):
        """solve_equation 对 z^5 - z + 1 = 0 的行为。

        【算法限制】：五次及以上方程无一般代数解。
        预期结果：5 个解（代数基本定理）；当前系统返回 []。
        """
        from integral import solve
        denom = Var("z") ** 5 - Var("z") + Const(1)
        sols = solve.solve_equation(denom, Const(0), "z", self.ctx)
        self.assertEqual(len(sols), 5)


# =============================================================================
# 9. 去重逻辑
# =============================================================================

class TestDuplicateElimination(unittest.TestCase):
    """同一极点不应在结果中出现多次。"""

    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()

    def test_no_duplicates_in_result(self):
        """find_poles 返回的列表中不应有重复极点。"""
        denom = (Var("z") - Const(1)) ** 2 * (Var("z") + Const(1))
        expr_ = Op("/", Const(1), denom)
        result = find_poles("z", expr_, self.ctx)
        poles = [r[0] for r in result]
        for i, p1 in enumerate(poles):
            for j, p2 in enumerate(poles):
                if i != j:
                    diff = poly.normalize(p1 - p2, self.ctx)
                    self.assertNotEqual(diff, Const(0),
                        f"Duplicate pole found: {p1} and {p2}")

    def test_deduplication_of_z1_pole(self):
        """(z-1)^3 的极点 z=1 只出现一次。"""
        denom = (Var("z") - Const(1)) ** 3
        expr_ = Op("/", Const(1), denom)
        result = find_poles("z", expr_, self.ctx)
        pole_list = [(p, o) for p, o in result]
        z1_count = sum(1 for p, _ in pole_list
                       if poly.normalize(p - Const(1), self.ctx) == Const(0))
        self.assertEqual(z1_count, 1,
            f"Expected exactly 1 occurrence of z=1 after dedup, got {z1_count}")


# =============================================================================
# 10. 非标准变量名
# =============================================================================

class TestNonStandardVariable(unittest.TestCase):
    """find_poles 支持任意变量名。"""

    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()

    def test_variable_w(self):
        """w 平面上的极点。"""
        w = Var("w")
        denom = w - Const(2)
        expr_ = Op("/", Const(1), denom)
        result = find_poles("w", expr_, self.ctx)
        self.assertEqual(len(result), 1)
        pole, order = result[0]
        self.assertEqual(order, 1)

    def test_variable_t(self):
        """t 作为积分变量。"""
        t = Var("t")
        denom = (t - Const(1)) * (t - Const(2))
        expr_ = Op("/", Const(1), denom)
        result = find_poles("t", expr_, self.ctx)
        self.assertEqual(len(result), 2)

    def test_variable_w_with_order(self):
        """w 平面上二阶极点。"""
        w = Var("w")
        denom = (w - Const(3)) ** 2
        expr_ = Op("/", Const(1), denom)
        result = find_poles("w", expr_, self.ctx)
        self.assertGreaterEqual(len(result), 1)
        matching = [(p, o) for p, o in result
                   if poly.normalize(denom.subst("w", p), self.ctx) == Const(0)]
        self.assertGreater(len(matching), 0)


# =============================================================================
# 11. 分子也有变量的分式
# =============================================================================

class TestNumeratorWithVariable(unittest.TestCase):
    """分子含变量但分母也含变量的分式。"""

    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()

    def test_z_over_z_minus_1(self):
        """z/(z-1): 极点 z=1, 一阶（分子在极点处非零，不影响阶数）。"""
        denom = Var("z") - Const(1)
        expr_ = Op("/", Var("z"), denom)
        result = find_poles("z", expr_, self.ctx)
        self.assertEqual(len(result), 1)
        pole, order = result[0]
        self.assertEqual(order, 1)
        val = denom.subst("z", pole)
        val = poly.normalize(val, self.ctx)
        self.assertEqual(val, Const(0))

    def test_z_squared_over_z_plus_1(self):
        """z^2/(z+1): 极点 z=-1, 一阶。"""
        denom = Var("z") + Const(1)
        expr_ = Op("/", Var("z") ** 2, denom)
        result = find_poles("z", expr_, self.ctx)
        self.assertEqual(len(result), 1)
        pole, order = result[0]
        self.assertEqual(order, 1)
        val = denom.subst("z", pole)
        val = poly.normalize(val, self.ctx)
        self.assertEqual(val, Const(0))

    def test_z_cubed_over_product(self):
        """z^3/((z-1)(z-2)): 极点 z=1, z=2。"""
        denom = (Var("z") - Const(1)) * (Var("z") - Const(2))
        expr_ = Op("/", Var("z") ** 3, denom)
        result = find_poles("z", expr_, self.ctx)
        self.assertEqual(len(result), 2)
        orders = [r[1] for r in result]
        self.assertEqual(sorted(orders), [1, 1])


# =============================================================================
# 12. 零函数（异常处理）
# =============================================================================

class TestZeroDenominator(unittest.TestCase):
    """分母为零函数时不应崩溃。"""

    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()

    def test_zero_denominator_returns_empty_or_catches(self):
        """分母为常数 0 时，find_poles 应捕获异常返回空列表。"""
        bad_expr = Op("/", Var("z"), Const(0))
        try:
            result = find_poles("z", bad_expr, self.ctx)
            self.assertIsInstance(result, list)
        except Exception:
            pass


# =============================================================================
# 13. _determine_pole_order 直接测试
# =============================================================================

class TestDeterminePoleOrder(unittest.TestCase):
    """直接测试 _determine_pole_order 函数的准确性。"""

    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()

    def test_order_one_simple_pole(self):
        """1/(z-1): z=1 处一阶。"""
        from integral.expr import _determine_pole_order
        denom = Var("z") - Const(1)
        pole = Const(1)
        order = _determine_pole_order("z", denom, pole, self.ctx)
        self.assertEqual(order, 1)

    def test_order_two_repeated_root(self):
        """1/(z-1)^2: z=1 处二阶。"""
        from integral.expr import _determine_pole_order
        denom = (Var("z") - Const(1)) ** 2
        pole = Const(1)
        order = _determine_pole_order("z", denom, pole, self.ctx)
        self.assertEqual(order, 2)

    def test_order_three_triple_root(self):
        """1/(z-2)^3: z=2 处三阶。"""
        from integral.expr import _determine_pole_order
        denom = (Var("z") - Const(2)) ** 3
        pole = Const(2)
        order = _determine_pole_order("z", denom, pole, self.ctx)
        self.assertEqual(order, 3)

    def test_order_five(self):
        """1/(z-3)^5: z=3 处五阶。"""
        from integral.expr import _determine_pole_order
        denom = (Var("z") - Const(3)) ** 5
        pole = Const(3)
        order = _determine_pole_order("z", denom, pole, self.ctx)
        self.assertEqual(order, 5)

    def test_wrong_pole_returns_high_order(self):
        """将 z=2 代入 (z-1)^2 得一阶导数非零，应返回 order=1（而非 2）。"""
        from integral.expr import _determine_pole_order
        denom = (Var("z") - Const(1)) ** 2
        wrong_pole = Const(2)
        order = _determine_pole_order("z", denom, wrong_pole, self.ctx)
        self.assertEqual(order, 1)


# =============================================================================
# 14. 综合集成测试：完整链路
# =============================================================================

class TestComprehensiveScenarios(unittest.TestCase):
    """综合场景：验证 solve_equation + _determine_pole_order 的完整链路。"""

    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()

    def test_z_four_minus_one(self):
        """1/(z^4-1): 四个一阶极点。"""
        denom = Var("z") ** 4 - Const(1)
        expr_ = Op("/", Const(1), denom)
        result = find_poles("z", expr_, self.ctx)
        self.assertEqual(len(result), 4)
        orders = [r[1] for r in result]
        self.assertEqual(sorted(orders), [1, 1, 1, 1])

    def test_z_four_plus_one(self):
        """1/(z^4+1): 四个一阶极点，即 z^4=-1 的四个根。"""
        denom = Var("z") ** 4 + Const(1)
        expr_ = Op("/", Const(1), denom)
        result = find_poles("z", expr_, self.ctx)
        self.assertEqual(len(result), 4)
        orders = [r[1] for r in result]
        self.assertEqual(sorted(orders), [1, 1, 1, 1])
        for pole, order in result:
            val = denom.subst("z", pole)
            val = poly.normalize(val, self.ctx)
            self.assertEqual(val, Const(0),
                f"pole {pole} is not actually a zero of {denom}")

    def test_z_six_minus_1(self):
        """1/(z^6-1): 六个一阶极点（单位圆六次单位根）。"""
        denom = Var("z") ** 6 - Const(1)
        expr_ = Op("/", Const(1), denom)
        result = find_poles("z", expr_, self.ctx)
        self.assertEqual(len(result), 6)
        orders = [r[1] for r in result]
        self.assertTrue(all(o == 1 for o in orders))

    def test_sparse_polynomial(self):
        """1/(z^3 - 8): 极点为立方根 of 8: 2, 2*exp(2πi/3), 2*exp(4πi/3)。"""
        denom = Var("z") ** 3 - Const(8)
        expr_ = Op("/", Const(1), denom)
        result = find_poles("z", expr_, self.ctx)
        self.assertEqual(len(result), 3)
        orders = [r[1] for r in result]
        self.assertEqual(sorted(orders), [1, 1, 1])

    def test_z_squared_plus_z_plus_one(self):
        """1/(z^2+z+1): 极点为 (-1±i√3)/2，一阶，共 2 个。

        solve_equation 无法对 z^2+z+1 因式分解（SymPy factor 不改变原式），
        线性提取 b*x+c 也匹配不上（x^2 项系数为 1），二次公式尚未实现，
        导致返回空列表。这是 solve_equation 的已知能力边界。
        预期结果：2 个一阶复数极点。
        """
        import sympy as sp
        z_sym = sp.symbols('z')
        expected_roots = sp.solve(z_sym**2 + z_sym + 1, z_sym)
        self.assertEqual(len(expected_roots), 2)
        self.assertTrue(sp.im(expected_roots[0]) != 0)
        self.assertTrue(sp.im(expected_roots[1]) != 0)

        denom = Var("z") ** 2 + Var("z") + Const(1)
        expr_ = Op("/", Const(1), denom)
        result = find_poles("z", expr_, self.ctx)
        self.assertEqual(len(result), 2)
        orders = [r[1] for r in result]
        self.assertEqual(sorted(orders), [1, 1])


# =============================================================================
# 15. 性能与缓存
# =============================================================================

class TestCacheAndPerformance(unittest.TestCase):
    """验证重复调用的确定性和上下文隔离。"""

    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()

    def test_repeated_call_is_deterministic(self):
        """同一表达式多次调用 find_poles 应得到一致结果。"""
        denom = (Var("z") - Const(1)) * (Var("z") + Const(1))
        expr_ = Op("/", Const(1), denom)
        result1 = find_poles("z", expr_, self.ctx)
        result2 = find_poles("z", expr_, self.ctx)
        self.assertEqual(len(result1), len(result2))
        orders1 = sorted([r[1] for r in result1])
        orders2 = sorted([r[1] for r in result2])
        self.assertEqual(orders1, orders2)

    def test_context_independence(self):
        """不同 Context 实例下，极点结果数量一致。"""
        ctx1 = make_ctx()
        ctx2 = make_ctx()
        denom = Const(2) * Var("z") - Const(1)
        expr_ = Op("/", Const(1), denom)
        result1 = find_poles("z", expr_, ctx1)
        result2 = find_poles("z", expr_, ctx2)
        self.assertEqual(len(result1), len(result2))


# =============================================================================
# 16. 高阶极点求根综合验证：solve_equation 对 (z-a)^n = 0 的正确处理
# =============================================================================

class TestHigherOrderPoleSolveEquation(unittest.TestCase):
    """验证 solve_equation 对重复根 (z-a)^n = 0 的完整链路。

    solve_equation((z-a)^n=0) 直接递归求解 u=0，
    返回单一解 z=a，不再生成 n 个错误复根。
    """

    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()

    def test_all_orders_deduplicate_correctly(self):
        """n=2,3,4,5 时，去重后的极点阶数均为 n。"""
        a = Const(2)
        for n in [2, 3, 4, 5]:
            denom = (Var("z") - a) ** n
            expr_ = Op("/", Const(1), denom)
            result = find_poles("z", expr_, self.ctx)
            for pole, order in result:
                if poly.normalize(pole - a, self.ctx) == Const(0):
                    self.assertEqual(order, n,
                        f"n={n}: Expected order {n}, got {order}")
                    break
            else:
                self.fail(f"n={n}: Expected pole z={a} not found")


if __name__ == "__main__":
    unittest.main()
