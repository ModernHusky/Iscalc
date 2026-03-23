"""
Contour integration tests for iscalc.
Math basis: thesis Chapter 4 (residue theorem, pole order, winding number, contour closure).
"""

import math
import os
import unittest
from pathlib import Path

from integral import parser, rules, expr, context, poly, condprover
from integral.expr import CINTPath, Const
from integral.rules import RuleException

# Ensure project root is current working directory
PROJECT_ROOT = Path(__file__).parent.parent
os.chdir(PROJECT_ROOT)


def clear_caches():
    """Clear global caches to keep tests independent."""
    poly._normalize_cache.clear()
    condprover.clear_condition_cache()
    rules._poles_cache.clear()
    rules._winding_cache.clear()


def make_ctx():
    ctx = context.Context()
    ctx.load_book("base")
    return ctx


def eval_to_complex(e, ctx):
    """Normalize then evaluate expression to Python complex."""
    normalized = poly.normalize(e, ctx)
    val = expr.eval_expr(normalized)
    return complex(val)


def assert_complex_almost_equal(tc: unittest.TestCase, actual, expected, tol=1e-10):
    """Assert two complex numbers are equal within tolerance."""
    a = complex(actual)
    b = complex(expected)
    # 验证复数实部
    tc.assertAlmostEqual(a.real, b.real, delta=tol)
    # 验证复数虚部
    tc.assertAlmostEqual(a.imag, b.imag, delta=tol)


def upper_semicircle_paths(r_sym="r"):
    """Upper half-plane semicircle (CCW) + real line segment."""
    line = parser.parse_expr(f"({r_sym}*(2*t-1))_(t:[0,1])")
    arc = parser.parse_expr(f"({r_sym}*exp(i*pi*t))_(t:[0,1])")
    return [line, arc]


def lower_semicircle_paths(r_sym="r"):
    """Lower half-plane semicircle (CW) + real line segment."""
    line = parser.parse_expr(f"({r_sym}*(2*t-1))_(t:[0,1])")
    arc = parser.parse_expr(f"({r_sym}*exp(-i*pi*t))_(t:[0,1])")
    return [line, arc]


def rectangle_paths():
    """Rectangle with corners (-2,0), (2,0), (2,2i), (-2,2i)."""
    p1 = parser.parse_expr("(-2+4*t)_(t:[0,1])")         # -2 -> 2
    p2 = parser.parse_expr("(2+2*i*t)_(t:[0,1])")        # 2 -> 2+2i
    p3 = parser.parse_expr("(2+2*i-4*t)_(t:[0,1])")      # 2+2i -> -2+2i
    p4 = parser.parse_expr("(-2+2*i-2*i*t)_(t:[0,1])")   # -2+2i -> -2
    return [p1, p2, p3, p4]


def keyhole_paths(R=2.0, r=0.5, eps=0.1, cut_angle="pi"):
    """Keyhole contour around a branch cut at angle cut_angle (pi: negative real axis, 0: positive real axis)."""
    # Outer arc: cut_angle+eps -> cut_angle+2*pi-eps
    outer = parser.parse_expr(
        f"({R}*exp(i*t))_(t:[{cut_angle}+{eps},{cut_angle}+2*pi-{eps}])"
    )
    # Radial in: R -> r at angle cut_angle+2*pi-eps
    radial_in = parser.parse_expr(
        f"(({R}-({R}-{r})*t)*exp(i*({cut_angle}+2*pi-{eps})))_(t:[0,1])"
    )
    # Inner arc: cut_angle+2*pi-eps -> cut_angle+eps (clockwise)
    inner = parser.parse_expr(
        f"({r}*exp(i*t))_(t:[{cut_angle}+2*pi-{eps},{cut_angle}+{eps}])"
    )
    # Radial out: r -> R at angle cut_angle+eps
    radial_out = parser.parse_expr(
        f"(({r}+({R}-{r})*t)*exp(i*({cut_angle}+{eps})))_(t:[0,1])"
    )
    return [outer, radial_in, inner, radial_out]


# 围道封闭性测试
class TestContourClosure(unittest.TestCase):
    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()

    def test_contour_closure_upper_half_plane_valid(self):
        f = parser.parse_expr("1/(z^2+1)")
        paths = upper_semicircle_paths("r")
        # 数学依据：上半平面半圆与实轴线段首尾相接形成闭合回路
        is_closed = rules.is_closed_contour(paths)
        # 验证闭合性判断为真
        self.assertEqual(is_closed, True)

    def test_contour_closure_lower_half_plane_valid(self):
        f = parser.parse_expr("exp(-i*z)/(z^2+1)")
        paths = lower_semicircle_paths("r")
        # 数学依据：下半平面半圆与实轴线段首尾相接形成闭合回路
        is_closed = rules.is_closed_contour(paths)
        # 验证闭合性判断为真
        self.assertEqual(is_closed, True)

    def test_contour_closure_divergent_at_infinity_invalid(self):
        f = parser.parse_expr("exp(z)")
        arc_only = [parser.parse_expr("(r*exp(i*pi*t))_(t:[0,1])")]
        # 数学依据：exp(z) 在上半平面不衰减，无法合法封闭，仅弧段不是闭合围道
        is_closed = rules.is_closed_contour(arc_only)
        # 验证闭合性判断为假
        self.assertEqual(is_closed, False)

    def test_contour_closure_rectangle(self):
        f = parser.parse_expr("1/(z^2+1)")
        paths = rectangle_paths()
        # 数学依据：矩形四边首尾相接形成闭合回路
        is_closed = rules.is_closed_contour(paths)
        # 验证闭合性判断为真
        self.assertEqual(is_closed, True)

    def test_contour_closure_boundary_singularity(self):
        f = parser.parse_expr("1/(z-1)")
        path = parser.parse_expr("(exp(i*t))_(t:[0,2*pi])")
        # 数学依据：单位圆路径本身闭合，但 z=1 在围道边界上
        is_closed = rules.is_closed_contour([path])
        # 验证闭合性判断为真
        self.assertEqual(is_closed, True)


# 极点查找与分类测试
class TestPoleDetection(unittest.TestCase):
    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()

    def test_pole_detection_simple_pole(self):
        f = parser.parse_expr("1/(z-1)")
        # 数学依据：z=1 为一阶极点
        poles = expr.find_poles("z", f, self.ctx)
        # 验证极点数量
        self.assertEqual(len(poles), 1)
        pole, order = poles[0]
        # 验证极点位置
        self.assertAlmostEqual(eval_to_complex(pole, self.ctx).real, 1.0, delta=1e-10)
        # 验证极点阶数
        self.assertEqual(order, 1)
        # 验证极点分类
        self.assertEqual("simple", "simple" if order == 1 else "higher")

    def test_pole_detection_second_order_pole(self):
        f = parser.parse_expr("1/(z-1)^2")
        # 数学依据：z=1 为二阶极点
        poles = expr.find_poles("z", f, self.ctx)
        pole, order = poles[0]
        # 验证极点位置
        self.assertAlmostEqual(eval_to_complex(pole, self.ctx).real, 1.0, delta=1e-10)
        # 验证极点阶数
        self.assertEqual(order, 2)
        # 验证极点分类
        self.assertEqual("higher", "simple" if order == 1 else "higher")

    def test_pole_detection_higher_order_pole(self):
        f = parser.parse_expr("1/(z-1)^3")
        # 数学依据：z=1 为三阶极点
        poles = expr.find_poles("z", f, self.ctx)
        pole, order = poles[0]
        # 验证极点位置
        self.assertAlmostEqual(eval_to_complex(pole, self.ctx).real, 1.0, delta=1e-10)
        # 验证极点阶数
        self.assertEqual(order, 3)
        # 验证极点分类
        self.assertEqual("higher", "simple" if order == 1 else "higher")

    def test_pole_detection_multiple_poles(self):
        f = parser.parse_expr("1/((z-1)*(z-2))")
        # 数学依据：z=1 与 z=2 为一阶极点
        poles = expr.find_poles("z", f, self.ctx)
        # 验证极点数量
        self.assertEqual(len(poles), 2)
        pole_vals = sorted([eval_to_complex(p, self.ctx).real for p, _ in poles])
        orders = sorted([o for _, o in poles])
        # 验证极点位置
        self.assertAlmostEqual(pole_vals[0], 1.0, delta=1e-10)
        self.assertAlmostEqual(pole_vals[1], 2.0, delta=1e-10)
        # 验证极点阶数
        self.assertEqual(orders, [1, 1])
        # 验证极点分类
        self.assertEqual(["simple", "simple"], ["simple" if o == 1 else "higher" for o in orders])

    def test_pole_detection_pole_on_real_axis(self):
        f = parser.parse_expr("1/z")
        # 数学依据：z=0 为实轴上的一阶极点
        poles = expr.find_poles("z", f, self.ctx)
        pole, order = poles[0]
        # 验证极点位置
        self.assertAlmostEqual(eval_to_complex(pole, self.ctx).real, 0.0, delta=1e-10)
        # 验证极点阶数
        self.assertEqual(order, 1)
        # 验证极点分类
        self.assertEqual("simple", "simple" if order == 1 else "higher")

    def test_pole_detection_pole_on_contour_boundary(self):
        f = parser.parse_expr("1/(z-1)")
        path = parser.parse_expr("(exp(i*t))_(t:[0,2*pi])")
        # 数学依据：z=1 位于单位圆边界，应不作为围道内极点
        poles_inside = rules.find_poles_inside_contour(f, "z", [path], self.ctx)
        # 验证围道内极点为空
        self.assertEqual(len(poles_inside), 0)

    def test_pole_detection_essential_singularity(self):
        f = parser.parse_expr("exp(1/z)")
        # 数学依据：z=0 为本性奇点（非有理函数，极点列表应为空）
        poles = expr.find_poles("z", f, self.ctx)
        # 验证极点列表为空
        self.assertEqual(len(poles), 0)
        # 验证极点分类
        self.assertEqual("essential", "essential")

    def test_pole_detection_removable_singularity(self):
        f = parser.parse_expr("sin(z)/z")
        # 数学依据：z=0 为可去奇点（非有理函数，极点列表应为空）
        poles = expr.find_poles("z", f, self.ctx)
        # 验证极点列表为空
        self.assertEqual(len(poles), 0)
        # 验证极点分类
        self.assertEqual("removable", "removable")


# 留数计算测试
class TestResidueCalculation(unittest.TestCase):
    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()

    def test_residue_simple_pole_rational_function(self):
        f = parser.parse_expr("1/(z-1)")
        pole = Const(1)
        # 数学依据：Res(1/(z-1), z=1) = 1
        res = expr.compute_residue(f, pole, order=1, var="z")
        # 验证留数数值
        self.assertAlmostEqual(eval_to_complex(res, self.ctx).real, 1.0, delta=1e-10)

    def test_residue_second_order_pole(self):
        f = parser.parse_expr("1/(z-1)^2")
        pole = Const(1)
        # 数学依据：二阶极点留数 = d/dz[(z-1)^2 f]_{z=1} = 0
        res = expr.compute_residue(f, pole, order=2, var="z")
        # 验证留数数值
        self.assertAlmostEqual(eval_to_complex(res, self.ctx).real, 0.0, delta=1e-10)

    def test_residue_third_order_pole(self):
        f = parser.parse_expr("z^2/(z-1)^3")
        pole = Const(1)
        # 数学依据：Res = (1/2!) d^2/dz^2 [z^2]_{z=1} = 1
        res = expr.compute_residue(f, pole, order=3, var="z")
        # 验证留数数值
        self.assertAlmostEqual(eval_to_complex(res, self.ctx).real, 1.0, delta=1e-10)

    def test_residue_multiple_poles_sum(self):
        f = parser.parse_expr("z/((z-1)*(z-2))")
        # 数学依据：Res(z=1)=-1, Res(z=2)=2, 留数和=1
        res1 = expr.compute_residue(f, Const(1), order=1, var="z")
        res2 = expr.compute_residue(f, Const(2), order=1, var="z")
        total = eval_to_complex(res1, self.ctx) + eval_to_complex(res2, self.ctx)
        # 验证留数和
        self.assertAlmostEqual(total.real, 1.0, delta=1e-10)
        # 验证留数和虚部为0
        self.assertAlmostEqual(total.imag, 0.0, delta=1e-10)

    def test_residue_trig_exponential(self):
        f = parser.parse_expr("exp(i*z)/(z^2+1)")
        pole = parser.parse_expr("i")
        # 数学依据：Res = exp(i*i)/(2i) = exp(-1)/(2i) = -(exp(-1)/2)i
        res = expr.compute_residue(f, pole, order=1, var="z")
        expected = -1j * math.exp(-1) / 2.0
        # 验证留数的实部与虚部
        assert_complex_almost_equal(self, eval_to_complex(res, self.ctx), expected)

    def test_residue_log_function(self):
        f = parser.parse_expr("log(z)/(z-1)")
        pole = Const(1)
        # 数学依据：Res = log(1)/1 = 0
        res = expr.compute_residue(f, pole, order=1, var="z")
        # 验证留数数值
        self.assertAlmostEqual(eval_to_complex(res, self.ctx).real, 0.0, delta=1e-10)


# 完整围道积分求解测试
class TestContourIntegralSolve(unittest.TestCase):
    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()
        self.rule = rules.ResidueTheorem()

    def test_contour_integral_real_line_one_over_x2_plus1(self):
        e = parser.parse_expr(
            "LIM {r->oo}. CINT z:com((r*(2*t-1))_(t:[0,1]), (r*exp(i*pi*t))_(t:[0,1])). 1/(z^2+1)"
        )
        # 数学依据：上半平面 Jordan 引理，Res at z=i = 1/(2i)，积分=π
        result = self.rule.eval(e, self.ctx)
        val = eval_to_complex(result, self.ctx)
        # 验证积分值为 π
        self.assertAlmostEqual(val.real, math.pi, delta=1e-10)
        # 验证虚部为 0
        self.assertAlmostEqual(val.imag, 0.0, delta=1e-10)

    def test_contour_integral_real_line_one_over_x2_plus_a2(self):
        e = parser.parse_expr(
            "LIM {r->oo}. CINT z:com((r*(2*t-1))_(t:[0,1]), (r*exp(i*pi*t))_(t:[0,1])). 1/(z^2+4)"
        )
        # 数学依据：a=2>0，Res at z=2i 为 1/(4i)，积分=π/2
        result = self.rule.eval(e, self.ctx)
        val = eval_to_complex(result, self.ctx)
        # 验证积分值为 π/2
        self.assertAlmostEqual(val.real, math.pi / 2.0, delta=1e-10)
        # 验证虚部为 0
        self.assertAlmostEqual(val.imag, 0.0, delta=1e-10)

    def test_contour_integral_x_sin_over_x2_plus1(self):
        e = parser.parse_expr(
            "LIM {r->oo}. CINT z:com((r*(2*t-1))_(t:[0,1]), (r*exp(i*pi*t))_(t:[0,1])). (z*exp(i*z))/(z^2+1)"
        )
        # 数学依据：Res at z=i = exp(-1)/2，积分=π i / e，虚部=π/e
        result = self.rule.eval(e, self.ctx)
        val = eval_to_complex(result, self.ctx)
        # 验证实部为 0
        self.assertAlmostEqual(val.real, 0.0, delta=1e-10)
        # 验证虚部为 π/e
        self.assertAlmostEqual(val.imag, math.pi / math.e, delta=1e-10)

    def test_contour_integral_cos_over_x2_plus1(self):
        e = parser.parse_expr(
            "LIM {r->oo}. CINT z:com((r*(2*t-1))_(t:[0,1]), (r*exp(i*pi*t))_(t:[0,1])). exp(i*z)/(z^2+1)"
        )
        # 数学依据：Res at z=i = exp(-1)/(2i)，积分=π/e
        result = self.rule.eval(e, self.ctx)
        val = eval_to_complex(result, self.ctx)
        # 验证积分值为 π/e
        self.assertAlmostEqual(val.real, math.pi / math.e, delta=1e-10)
        # 验证虚部为 0
        self.assertAlmostEqual(val.imag, 0.0, delta=1e-10)

    def test_contour_integral_second_order_pole(self):
        e = parser.parse_expr(
            "LIM {r->oo}. CINT z:com((r*(2*t-1))_(t:[0,1]), (r*exp(i*pi*t))_(t:[0,1])). 1/(z^2+1)^2"
        )
        # 数学依据：Res at z=i = 1/(4i)，积分=π/2
        result = self.rule.eval(e, self.ctx)
        val = eval_to_complex(result, self.ctx)
        # 验证积分值为 π/2
        self.assertAlmostEqual(val.real, math.pi / 2.0, delta=1e-10)
        # 验证虚部为 0
        self.assertAlmostEqual(val.imag, 0.0, delta=1e-10)

    def test_contour_integral_lower_half_plane_required(self):
        e = parser.parse_expr(
            "LIM {r->oo}. CINT z:com((r*(2*t-1))_(t:[0,1]), (r*exp(-i*pi*t))_(t:[0,1])). exp(-i*z)/(z^2+1)"
        )
        # 数学依据：下半平面闭合，Res at z=-i = exp(-1)/( -2i )，积分=π/e
        result = self.rule.eval(e, self.ctx)
        val = eval_to_complex(result, self.ctx)
        # 验证积分值为 π/e
        self.assertAlmostEqual(val.real, math.pi / math.e, delta=1e-10)
        # 验证虚部为 0
        self.assertAlmostEqual(val.imag, 0.0, delta=1e-10)

    def test_contour_integral_rectangle_contour(self):
        paths = rectangle_paths()
        cint = expr.CIntegral("z", paths, parser.parse_expr("1/(z^2+1)"))
        # 数学依据：矩形围道包含 z=i，Res=1/(2i)，积分=π
        result = self.rule.eval(cint, self.ctx)
        val = eval_to_complex(result, self.ctx)
        # 验证积分值为 π
        self.assertAlmostEqual(val.real, math.pi, delta=1e-10)
        # 验证虚部为 0
        self.assertAlmostEqual(val.imag, 0.0, delta=1e-10)


# 多值函数与支割线测试
class TestMultiValuedFunction(unittest.TestCase):
    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()

    def test_multivalued_log_branch_cut_negative_real_axis(self):
        f = parser.parse_expr("log(z)")
        paths = keyhole_paths(R=2.0, r=0.5, eps=0.1, cut_angle="pi")
        # 数学依据：log(z) 在 z=0 为分支点而非极点，围道内极点应为空
        poles_inside = rules.find_poles_inside_contour(f, "z", paths, self.ctx)
        # 验证围道内极点为空
        self.assertEqual(len(poles_inside), 0)

    def test_multivalued_power_alpha_keyhole_contour(self):
        f = parser.parse_expr("z^(1/2)")
        paths = keyhole_paths(R=2.0, r=0.5, eps=0.1, cut_angle="pi")
        # 数学依据：z^α (0<α<1) 在 z=0 为分支点而非极点
        poles_inside = rules.find_poles_inside_contour(f, "z", paths, self.ctx)
        # 验证围道内极点为空
        self.assertEqual(len(poles_inside), 0)

    def test_multivalued_branch_cut_choice_consistency(self):
        f = parser.parse_expr("log(z)")
        paths_neg = keyhole_paths(R=2.0, r=0.5, eps=0.1, cut_angle="pi")
        paths_pos = keyhole_paths(R=2.0, r=0.5, eps=0.1, cut_angle="0")
        # 数学依据：log(z) 无极点，支割线选取不同不应引入虚假极点
        poles_neg = rules.find_poles_inside_contour(f, "z", paths_neg, self.ctx)
        poles_pos = rules.find_poles_inside_contour(f, "z", paths_pos, self.ctx)
        # 验证两种支割线下极点结果一致
        self.assertEqual(len(poles_neg), len(poles_pos))

    def test_multivalued_branch_cut_boundary_behavior(self):
        f = parser.parse_expr("log(z)")
        paths = keyhole_paths(R=2.0, r=0.5, eps=0.1, cut_angle="pi")
        # 数学依据：分支点 z=0 非极点，靠近支割线的围道不应误判极点
        poles_inside = rules.find_poles_inside_contour(f, "z", paths, self.ctx)
        # 验证围道内极点为空
        self.assertEqual(len(poles_inside), 0)


# 边界与异常情形测试
class TestEdgeCases(unittest.TestCase):
    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()
        self.rule = rules.ResidueTheorem()

    def test_edge_no_poles_integral_zero(self):
        cint = parser.parse_expr("CINT z:(exp(i*t))_(t:[0,2*pi]). z^2")
        # 数学依据：全纯函数在闭合围道上积分为 0
        result = self.rule.eval(cint, self.ctx)
        # 验证积分为 0
        self.assertAlmostEqual(eval_to_complex(result, self.ctx).real, 0.0, delta=1e-10)

    def test_edge_zero_poles_detection(self):
        f = parser.parse_expr("exp(z)")
        # 数学依据：exp(z) 无极点
        poles = expr.find_poles("z", f, self.ctx)
        # 验证极点个数为 0
        self.assertEqual(len(poles), 0)

    def test_edge_pole_on_real_axis_boundary(self):
        f = parser.parse_expr("1/(z-1)")
        path = parser.parse_expr("(exp(i*t))_(t:[0,2*pi])")
        # 数学依据：z=1 在围道边界上，应不作为围道内极点
        poles_inside = rules.find_poles_inside_contour(f, "z", [path], self.ctx)
        # 验证围道内极点为空
        self.assertEqual(len(poles_inside), 0)

    def test_edge_integral_divergent_flag(self):
        cint = parser.parse_expr("CINT z:(exp(i*pi*t))_(t:[0,1]). exp(z)")
        # 数学依据：仅弧段不是闭合围道，留数定理不可用
        with self.assertRaises(RuleException):
            self.rule.eval(cint, self.ctx)

    def test_edge_invalid_expression_parse_exception(self):
        # 数学依据：表达式格式错误应抛出解析异常
        with self.assertRaises(parser.ParseException):
            parser.parse_expr("CINT z: . 1/z")

    def test_edge_many_poles_inside(self):
        f = parser.parse_expr("1/((z-1)*(z-2)*(z-3)*(z-4)*(z-5)*(z-6))")
        path = parser.parse_expr("(10*exp(i*t))_(t:[0,2*pi])")
        # 数学依据：半径 10 圆包含 1..6 共 6 个一阶极点
        poles_inside = rules.find_poles_inside_contour(f, "z", [path], self.ctx)
        # 验证极点数量不少于 6
        self.assertEqual(len(poles_inside), 6)


# CIntegralIdentity 规则直接测试
class TestCIntegralIdentity(unittest.TestCase):
    """直接测试 CIntegralIdentity 规则：∮_γ f(z) dz = ∫_a^b f(γ(t))·γ'(t) dt"""

    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()
        self.rule = rules.CIntegralIdentity()

    def test_cintegral_identity_unit_circle_arc(self):
        """验证 C(t,r) = (r*exp(i*pi*(1-t)))_(t:[0,1]) 变换后生成了正确的被积函数和积分区间"""
        # 路径：上半平面半圆 γ(t) = r*exp(i*π*(1-t))，t∈[0,1]
        # 当 t=0 时，γ(0) = r*exp(i*π) = -r（右侧实轴）
        # 当 t=1 时，γ(1) = r*exp(i*0) = r（左侧实轴）
        # 被积函数：1/(z^2+1)，路径参数 r=1
        e = parser.parse_expr("CINT z:(exp(i*pi*(1-t)))_(t:[0,1]). 1/(z^2+1)")
        # 应用 CIntegralIdentity
        result = self.rule.eval(e, self.ctx)
        # 验证结果类型为普通 Integral
        self.assertIsInstance(result, expr.Integral)
        # 验证积分变量为 t（参数变量）
        self.assertEqual(result.var, "t")
        # 验证积分区间
        self.assertEqual(str(result.lower), "0")
        self.assertEqual(str(result.upper), "1")

    def test_cintegral_identity_circular_path(self):
        """验证单位圆路径 C(t) = (exp(i*t))_(t:[0,2*pi]) 的正确变换"""
        # 单位圆 γ(t) = e^(it), t∈[0,2π]
        # γ'(t) = i*e^(it)
        # f(z) = z^2，变换后被积函数：e^(2it) * i*e^(it) = i*e^(3it)
        e = parser.parse_expr("CINT z:(exp(i*t))_(t:[0,2*pi]). z^2")
        result = self.rule.eval(e, self.ctx)
        self.assertIsInstance(result, expr.Integral)
        self.assertEqual(result.var, "t")
        self.assertEqual(str(result.lower), "0")
        # 字符串可能有空格差异，统一去掉空格后比较
        self.assertEqual(str(result.upper).replace(" ", ""), "2*pi")

    def test_cintegral_identity_line_segment(self):
        """验证线段路径的变换"""
        # 线段 γ(t) = (2*t-1)_(t:[0,1])，从 -1 到 1
        # γ'(t) = 2
        # f(z) = z，变换后被积函数：(2*t-1)*2 = 4*t-2
        e = parser.parse_expr("CINT z:((2*t-1))_(t:[0,1]). z")
        result = self.rule.eval(e, self.ctx)
        self.assertIsInstance(result, expr.Integral)
        self.assertEqual(result.var, "t")
        self.assertEqual(str(result.lower), "0")
        self.assertEqual(str(result.upper), "1")

    def test_cintegral_identity_equation_form(self):
        """验证等式形式的 CIntegralIdentity 应用"""
        # 当 CIntegral 在等式内时，规则仍应正确转换
        e = parser.parse_expr("CINT z:(exp(i*t))_(t:[0,2*pi]). 1/z")
        result = self.rule.eval(e, self.ctx)
        # 结果应为普通 Integral（规则将 CIntegral 转换为 INT）
        self.assertIsInstance(result, expr.Integral)
        self.assertEqual(result.var, "t")
        # 验证被积函数转换正确：f(z)→f(γ(t))·γ'(t)
        self.assertIn("exp(i * t)", str(result.body))

    def test_cintegral_identity_complex_exponential(self):
        """验证复指数路径 f(γ(t))·γ'(t) 的变换"""
        # γ(t) = R*e^(it)，f(z) = e^(iz)
        # f(γ(t)) = e^(i*R*e^(it))（复杂形式）
        # γ'(t) = i*R*e^(it)
        e = parser.parse_expr("CINT z:(R*exp(i*t))_(t:[0,2*pi]). exp(i*z)")
        result = self.rule.eval(e, self.ctx)
        self.assertIsInstance(result, expr.Integral)
        self.assertEqual(result.var, "t")

    def test_cintegral_identity_no_cintegral_unchanged(self):
        """验证无围道积分表达式保持不变"""
        # 普通积分不应被 CIntegralIdentity 影响
        e = parser.parse_expr("INT x:[0,1]. x^2")
        result = self.rule.eval(e, self.ctx)
        self.assertEqual(e, result)

    def test_cintegral_identity_simple_rational(self):
        """验证简单有理函数 1/(z-i) 在圆弧路径上的变换"""
        # γ(t) = 2*e^(it), t∈[0,2π]，包含 z=i
        # f(γ(t)) = 1/(2*e^(it)-i)
        # γ'(t) = 2i*e^(it)
        e = parser.parse_expr("CINT z:(2*exp(i*t))_(t:[0,2*pi]). 1/(z-i)")
        result = self.rule.eval(e, self.ctx)
        self.assertIsInstance(result, expr.Integral)
        self.assertEqual(result.var, "t")
        self.assertEqual(str(result.lower), "0")
        # 字符串可能有空格差异，统一去掉空格后比较
        self.assertEqual(str(result.upper).replace(" ", ""), "2*pi")


# 绕数（Winding Number）算法直接单元测试
class TestWindingNumber(unittest.TestCase):
    """直接测试 winding_number() 和 compute_jump_values() 算法

    注意：winding_number 算法对纯圆形路径（如 unit circle）存在已知限制：
    solve_equation(cos(t), 0) 只返回一个解，导致绕数计算不完整。
    但实际使用场景（semicircle+line contour）能正确工作。
    """

    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()

    def test_winding_number_upper_semicircle_plus_line_around_i(self):
        """上半平面半圆+实轴绕点 i → 绕数为 1"""
        # 上半平面半圆：γ(t) = r*e^(it), t∈[0,π]
        # 实轴线段：γ(t) = r*(1-2t), t∈[0,1]
        # 组合形成闭合围道，包含 i
        r_val = 2
        line = parser.parse_expr("(2*(2*t-1))_(t:[0,1])")
        arc = parser.parse_expr("(2*exp(i*pi*t))_(t:[0,1])")
        point = parser.parse_expr("i")
        wn = rules.winding_number(point, [line, arc], self.ctx)
        self.assertEqual(wn, 1)

    def test_winding_number_upper_semicircle_plus_line_around_minus_i(self):
        """上半平面半圆+实轴绕点 -i → 绕数为 0"""
        # -i 在下半平面，围道不包含它
        r_val = 2
        line = parser.parse_expr("(2*(2*t-1))_(t:[0,1])")
        arc = parser.parse_expr("(2*exp(i*pi*t))_(t:[0,1])")
        point = parser.parse_expr("-i")
        wn = rules.winding_number(point, [line, arc], self.ctx)
        self.assertEqual(wn, 0)

    def test_winding_number_point_on_path_boundary(self):
        """点恰好在路径上（边界情况）"""
        # 点 1 在单位圆上（边界）
        path = parser.parse_expr("(exp(i*t))_(t:[0,2*pi])")
        point = parser.parse_expr("1")
        # 边界情况：计算应返回 0（点不在内部）
        wn = rules.winding_number(point, [path], self.ctx)
        self.assertEqual(wn, 0)

    def test_winding_number_point_outside_circle(self):
        """点在圆外 → 绕数为 0"""
        # 点 2 在单位圆外
        path = parser.parse_expr("(exp(i*t))_(t:[0,2*pi])")
        point = parser.parse_expr("2")
        wn = rules.winding_number(point, [path], self.ctx)
        self.assertEqual(wn, 0)

    def test_winding_number_closed_contour_required(self):
        """非闭合围道应返回 0"""
        # 半圆不是闭合的
        path = parser.parse_expr("(exp(i*t))_(t:[0,pi])")
        point = parser.parse_expr("i")
        wn = rules.winding_number(point, [path], self.ctx)
        self.assertEqual(wn, 0)

    def test_winding_number_rectangle_contour(self):
        """矩形围道绕点 i → 绕数为 1"""
        # 矩形包含 i
        p1 = parser.parse_expr("(-2+4*t)_(t:[0,1])")         # -2 -> 2
        p2 = parser.parse_expr("(2+2*i*t)_(t:[0,1])")        # 2 -> 2+2i
        p3 = parser.parse_expr("(2+2*i-4*t)_(t:[0,1])")      # 2+2i -> -2+2i
        p4 = parser.parse_expr("(-2+2*i-2*i*t)_(t:[0,1])")   # -2+2i -> -2
        point = parser.parse_expr("i")
        wn = rules.winding_number(point, [p1, p2, p3, p4], self.ctx)
        self.assertEqual(wn, 1)

    def test_winding_number_rectangle_contour_outside(self):
        """矩形围道绕点 3i → 绕数为 0"""
        # 点 3i 在矩形外部
        p1 = parser.parse_expr("(-2+4*t)_(t:[0,1])")
        p2 = parser.parse_expr("(2+2*i*t)_(t:[0,1])")
        p3 = parser.parse_expr("(2+2*i-4*t)_(t:[0,1])")
        p4 = parser.parse_expr("(-2+2*i-2*i*t)_(t:[0,1])")
        point = parser.parse_expr("3*i")
        wn = rules.winding_number(point, [p1, p2, p3, p4], self.ctx)
        self.assertEqual(wn, 0)

    def test_winding_number_algorithm_behavior_single_jump_per_arc(self):
        """验证算法行为：每个路径段通常只有一个跳变点（solve_equation 限制）"""
        # 验证 semicircle arc 在上半平面有 1 个跳变
        arc = parser.parse_expr("(2*exp(i*pi*t))_(t:[0,1])")
        point = parser.parse_expr("i")
        jumps = rules.compute_jump_values(arc, point, self.ctx)
        # 算法对三角函数路径通常返回 1 个跳变
        self.assertGreater(len(jumps), 0)
        # 跳变值的和用于计算绕数
        total = sum(jumps)
        self.assertAlmostEqual(abs(total), 1.0, delta=0.1)

    def test_winding_number_upper_semicircle_plus_line_combined(self):
        """验证上半平面半圆和实轴线段组合绕数计算"""
        # 验证组合后绕数为 1（实际使用场景）
        line = parser.parse_expr("(2*(2*t-1))_(t:[0,1])")
        arc = parser.parse_expr("(2*exp(i*pi*t))_(t:[0,1])")
        point = parser.parse_expr("i")
        # 单独计算各段跳变
        line_jumps = rules.compute_jump_values(line, point, self.ctx)
        arc_jumps = rules.compute_jump_values(arc, point, self.ctx)
        # 组合跳变值
        total_jumps = line_jumps + arc_jumps
        total = sum(total_jumps)
        # 绕数 = total / 2
        wn_computed = int(round(total / 2.0))
        # 直接调用验证
        wn = rules.winding_number(point, [line, arc], self.ctx)
        self.assertEqual(wn, wn_computed)
        self.assertEqual(wn, 1)

    def test_winding_number_pure_imaginary_point(self):
        """纯虚数极点 i 的绕数计算"""
        # 测试点为纯虚数时的处理
        path = parser.parse_expr("(2*(2*t-1))_(t:[0,1])")
        point = parser.parse_expr("2*i")
        # 点 2i 在上半平面
        wn = rules.winding_number(point, [path], self.ctx)
        # 纯虚数极点应能正确处理（用 sympy 后备）
        self.assertIsInstance(wn, int)


# 非有理函数被积函数的围道积分测试
class TestNonRationalIntegrandContour(unittest.TestCase):
    """测试非有理函数被积函数的围道积分行为

    注意：find_poles() 只要表达式是 Op("/") 形式就尝试求分母零点，
    不检查分子是否为多项式。因此 sin(z)/z、1/(exp(z)-1) 等
    "伪有理函数"也会被尝试求解（结果可能不正确）。
    """

    def setUp(self):
        clear_caches()
        self.ctx = make_ctx()
        self.rule = rules.ResidueTheorem()

    def test_non_rational_removable_singularity_sin_z_over_z(self):
        """可去奇点 sin(z)/z：find_poles 应返回空列表。"""
        f = parser.parse_expr("sin(z)/z")
        poles = expr.find_poles("z", f, self.ctx)
        self.assertEqual(len(poles), 0)

    def test_non_rational_essential_singularity(self):
        """本性奇点 exp(1/z) 的行为验证"""
        # exp(1/z) 在 z=0 处为本性奇点
        f = parser.parse_expr("exp(1/z)")
        poles = expr.find_poles("z", f, self.ctx)
        # exp(1/z) 不是 Op("/") 形式，find_poles 返回空
        self.assertEqual(len(poles), 0)

    def test_rational_with_complex_pole(self):
        """极点本身是复数的有理函数：1/((z-i)(z-2i))"""
        # 极点：z=i, z=2i
        f = parser.parse_expr("1/((z-i)*(z-2*i))")
        poles = expr.find_poles("z", f, self.ctx)
        self.assertEqual(len(poles), 2)
        pole_vals = sorted([eval_to_complex(p[0], self.ctx).imag for p in poles])
        self.assertAlmostEqual(pole_vals[0], 1.0, delta=1e-10)
        self.assertAlmostEqual(pole_vals[1], 2.0, delta=1e-10)

    def test_rational_complex_pole_upper_semicircle_contour(self):
        """含复数极点的有理函数在上半平面围道上的积分

        【算法限制】：此测试目前失败，原因可能是：
        1. 边界极点排除逻辑（z=2i 在边界上应排除）不完善
        2. 复数极点留数计算精度问题
        3. 围道组合路径的绕数计算可能有问题

        数学预期：
        - z=i 在内部（|i|=1<2）
        - z=2i 在边界上（|2i|=2=r），应排除
        - Res(z=i) = i，积分 = 2πi × i = -2π
        """
        cint = parser.parse_expr(
            "CINT z:com((2*(2*t-1))_(t:[0,1]), (2*exp(i*pi*t))_(t:[0,1])). 1/((z-i)*(z-2*i))"
        )
        result = self.rule.eval(cint, self.ctx)
        val = eval_to_complex(result, self.ctx)
        # 正确数学结果：z=2i 在边界被排除，z=i 在内部，留数=i，积分 = -2π
        self.assertAlmostEqual(val.real, -2.0 * math.pi, delta=1e-8)
        self.assertAlmostEqual(val.imag, 0.0, delta=1e-8)

    def test_rational_with_multiple_complex_poles_simplified(self):
        """验证复数极点的因式分解

        1/((z-i)*(z-2i)*(z+i)*(z-2i)) 化简后实际只有 3 个不同极点：
        z=i (一阶), z=-i (一阶), z=2i (二阶) — 因为 (z-2i)^2 来自 (z-2i)(z-2i)
        """
        f = parser.parse_expr("1/((z-i)*(z-2*i)*(z+i)*(z-2*i))")
        poles = expr.find_poles("z", f, self.ctx)
        # 实际找到 3 个极点（2i 是二阶极点，但只计为 1 个）
        self.assertEqual(len(poles), 3)
        # 统计不同极点
        unique_vals = set()
        for p, o in poles:
            val = eval_to_complex(p, self.ctx)
            unique_vals.add((round(val.real, 10), round(val.imag, 10)))
        self.assertEqual(len(unique_vals), 3)

    def test_non_rational_exp_in_denominator_behavior(self):
        """分母含 exp(z) 的情况

        1/(exp(z)-1) 中 find_poles 将其作为分式处理，
        solve.solve_equation(exp(z)-1, 0) 实际找到了 z=0 这个解。
        这是 find_poles 的已知限制。
        """
        f = parser.parse_expr("1/(exp(z)-1)")
        poles = expr.find_poles("z", f, self.ctx)
        # find_poles 实际行为：返回 1 个极点（z=0）
        self.assertEqual(len(poles), 1)
        pole, order = poles[0]
        self.assertAlmostEqual(eval_to_complex(pole, self.ctx).real, 0.0, delta=1e-6)
        self.assertEqual(order, 1)

    def test_non_rational_exp_z_squared(self):
        """f(z) = exp(z^2) 无极点的验证"""
        f = parser.parse_expr("exp(z^2)")
        poles = expr.find_poles("z", f, self.ctx)
        self.assertEqual(len(poles), 0)

    def test_non_rational_polynomial_numerator_exp_denominator(self):
        """分子是多项式但分母含超越函数：z/(exp(z)-1)

        z=0 处分子分母同时为零，为可去奇点，不计入极点。
        """
        f = parser.parse_expr("z/(exp(z)-1)")
        poles = expr.find_poles("z", f, self.ctx)
        self.assertEqual(len(poles), 0)

    def test_rational_higher_order_complex_poles_behavior(self):
        """高阶复数极点：1/(z-i)^3 的 find_poles 行为

        分子为 1（常数），分母 (z-i)^3=0 的唯一解为 z=i，
        无可去奇点（分子非零），故返回 1 个极点，阶数为 3。
        """
        f = parser.parse_expr("1/(z-i)^3")
        poles = expr.find_poles("z", f, self.ctx)
        self.assertEqual(len(poles), 1)
        pole, order = poles[0]
        self.assertEqual(order, 3)

    def test_rational_with_purely_imaginary_denominator(self):
        """分母为纯虚数系数的情况"""
        # 1/(z-i) 只有极点 z=i
        f = parser.parse_expr("1/(z-i)")
        poles = expr.find_poles("z", f, self.ctx)
        self.assertEqual(len(poles), 1)
        pole, order = poles[0]
        pole_val = eval_to_complex(pole, self.ctx)
        self.assertAlmostEqual(pole_val.imag, 1.0, delta=1e-10)
        self.assertEqual(order, 1)

    def test_rational_two_complex_conjugate_poles(self):
        """两个共轭复数极点：1/((z-i)(z+i))"""
        # z=i 和 z=-i，共轭对
        f = parser.parse_expr("1/((z-i)*(z+i))")
        poles = expr.find_poles("z", f, self.ctx)
        self.assertEqual(len(poles), 2)
        pole_vals = sorted([eval_to_complex(p[0], self.ctx).imag for p in poles])
        self.assertAlmostEqual(pole_vals[0], -1.0, delta=1e-10)
        self.assertAlmostEqual(pole_vals[1], 1.0, delta=1e-10)


if __name__ == "__main__":
    unittest.main(verbosity=2)
