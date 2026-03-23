"""
Comprehensive robustness tests for `is_closed_contour`.

Each test follows the documented contract:
  "如果围道闭合返回True，否则返回False"

Run with: python -m pytest tests/test_is_closed_contour.py -v
"""

import os
import unittest
from pathlib import Path

from integral import parser, rules, expr, context, poly, condprover
from integral.expr import CINTPath, Const, Var, Fun, Op

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


def parse_path(s: str):
    return parser.parse_expr(s)


# ---------------------------------------------------------------------------
# 1. 边界情况
# ---------------------------------------------------------------------------

class TestIsClosedContourEdgeCases(unittest.TestCase):
    """Empty / None / non-CINTPath inputs — documented contract: not closed → False."""

    def setUp(self):
        clear_caches()

    def test_empty_list_returns_false(self):
        """空列表 → 按文档应返回 False（不是闭合围道）。"""
        self.assertFalse(rules.is_closed_contour([]))

    def test_single_none_returns_false(self):
        """单个 None → 不是 CINTPath，文档语义应返回 False（无有效路径）。"""
        self.assertFalse(rules.is_closed_contour([None]))

    def test_single_string_returns_false(self):
        """单个字符串 → 不是 CINTPath，文档语义应返回 False。"""
        self.assertFalse(rules.is_closed_contour(["not a path"]))

    def test_single_int_returns_false(self):
        """单个整数 → 不是 CINTPath，文档语义应返回 False。"""
        self.assertFalse(rules.is_closed_contour([42]))

    def test_all_non_cintpath_returns_false(self):
        """全部是非 CINTPath 对象 → 无有效路径，应返回 False。"""
        self.assertFalse(rules.is_closed_contour(["bad", 123, None, {}]))

    def test_nested_list_returns_false(self):
        """外层列表包含子列表而非 CINTPath → 无有效路径，应返回 False。"""
        path = parse_path("(exp(i*t))_(t:[0,2*pi])")
        self.assertFalse(rules.is_closed_contour([[path]]))

    def test_mixed_valid_and_invalid_returns_true_only_if_valid_closed(self):
        """混合路径：无效项被跳过；有效路径本身闭合则返回 True。"""
        closed = parse_path("(exp(i*t))_(t:[0,2*pi])")
        self.assertTrue(rules.is_closed_contour([closed, None, "bad"]))
        self.assertTrue(rules.is_closed_contour(["bad", closed]))

    def test_mixed_valid_and_invalid_open_returns_false(self):
        """混合路径：无效项被跳过；有效路径本身开口则返回 False。"""
        open_arc = parse_path("(exp(i*t))_(t:[0,pi])")
        self.assertFalse(rules.is_closed_contour([open_arc, None, "bad"]))
        self.assertFalse(rules.is_closed_contour(["bad", open_arc]))

    def test_bad_path_among_good_open_breaks_chain(self):
        """多路径中一个无效（跳过）+ 其余路径开口 → 返回 False。"""
        p1 = parse_path("(t)_(t:[0,1])")    # 0 → 1
        p2 = parse_path("(1+t)_(t:[0,1])")  # 1 → 2
        # bad 被跳过；p1+p2 不形成闭合 → False
        self.assertFalse(rules.is_closed_contour([p1, "bad", p2]))

    def test_bad_path_among_good_closed_returns_true(self):
        """多路径中一个无效（跳过）+ 其余路径闭合 → 返回 True。"""
        circle = parse_path("(exp(i*t))_(t:[0,2*pi])")  # 闭合
        self.assertTrue(rules.is_closed_contour([circle, "bad"]))


# ---------------------------------------------------------------------------
# 2. 单一路径 — 闭合
# ---------------------------------------------------------------------------

class TestIsClosedContourSingleClosed(unittest.TestCase):
    """单一 CINTPath 本身闭合的用例。"""

    def setUp(self):
        clear_caches()

    def test_unit_circle_ccw(self):
        """标准单位圆 CCW t:[0,2*pi] → 闭合。"""
        path = parse_path("(exp(i*t))_(t:[0,2*pi])")
        self.assertTrue(rules.is_closed_contour([path]))

    def test_unit_circle_cw(self):
        """标准单位圆 CW exp(-i*t) t:[0,2*pi] → 闭合。"""
        path = parse_path("(exp(-i*t))_(t:[0,2*pi])")
        self.assertTrue(rules.is_closed_contour([path]))

    def test_arbitrary_radius_ccw(self):
        """任意实数半径 CCW → 闭合。"""
        path = parse_path("(3*exp(i*t))_(t:[0,2*pi])")
        self.assertTrue(rules.is_closed_contour([path]))

    def test_float_radius_ccw(self):
        """浮点半径 CCW → 闭合。"""
        path = parse_path("(3.5*exp(i*t))_(t:[0,2*pi])")
        self.assertTrue(rules.is_closed_contour([path]))

    def test_symbolic_radius_exp(self):
        """符号半径 r*exp(i*t) → is_closed 特殊处理 exp(i*t) → 闭合。"""
        path = parse_path("(r*exp(i*t))_(t:[0,2*pi])")
        self.assertTrue(rules.is_closed_contour([path]))

    def test_zero_length_path(self):
        """零长度路径 t:[0,0] → start==end 快捷路径 → 闭合。"""
        path = parse_path("(t)_(t:[0,0])")
        self.assertTrue(rules.is_closed_contour([path]))

    def test_constant_path(self):
        """退化路径（始终为常数）→ 闭合。"""
        path = parse_path("(5)_(t:[0,1])")
        self.assertTrue(rules.is_closed_contour([path]))


# ---------------------------------------------------------------------------
# 3. 单一路径 — 不闭合
# ---------------------------------------------------------------------------

class TestIsClosedContourSingleOpen(unittest.TestCase):
    """单一 CINTPath 本身不闭合的用例。"""

    def setUp(self):
        clear_caches()

    def test_unit_circle_half_turn(self):
        """半圈 t:[0,pi] → 不闭合（终点 -1 ≠ 起点 1）。"""
        path = parse_path("(exp(i*t))_(t:[0,pi])")
        self.assertFalse(rules.is_closed_contour([path]))

    def test_unit_circle_quarter(self):
        """1/4 圈 t:[0,pi/2] → 不闭合。"""
        path = parse_path("(exp(i*t))_(t:[0,pi/2])")
        self.assertFalse(rules.is_closed_contour([path]))

    def test_float_radius_quarter(self):
        """浮点半径 1/4 圈 → 不闭合。"""
        path = parse_path("(2.0*exp(i*t))_(t:[0,pi/2])")
        self.assertFalse(rules.is_closed_contour([path]))

    def test_line_segment_0_to_1(self):
        """线段 t:[0,1] → 不闭合。"""
        path = parse_path("(t)_(t:[0,1])")
        self.assertFalse(rules.is_closed_contour([path]))

    def test_symmetric_line_minus1_to_1(self):
        """对称线段 -1→1 → 不闭合。"""
        path = parse_path("(2*t-1)_(t:[0,1])")
        self.assertFalse(rules.is_closed_contour([path]))

    def test_pure_imaginary_axis(self):
        """纯虚轴 i→0 → 不闭合。"""
        path = parse_path("(i*(1-t))_(t:[0,1])")
        self.assertFalse(rules.is_closed_contour([path]))

    def test_unrelated_symbol(self):
        """带无关符号的闭合参数范围 → are_equal 可能失败，fallback is_closed() → bool。"""
        path = parse_path("(a*exp(i*t))_(t:[0,2*pi])")
        result = rules.is_closed_contour([path])
        self.assertIsInstance(result, bool)


# ---------------------------------------------------------------------------
# 4. 多路径 — 闭合
# ---------------------------------------------------------------------------

class TestIsClosedContourMultiClosed(unittest.TestCase):
    """多段路径首尾相连形成闭合围道。"""

    def setUp(self):
        clear_caches()

    def test_rectangle_4_edges(self):
        """矩形四条边 → 闭合。"""
        hw = 1.0  # half-width
        h = 2.0   # height
        p1 = parse_path(f"(-{hw}+{2*hw}*t)_(t:[0,1])")             # left→right
        p2 = parse_path(f"({hw}+{h}*i*t)_(t:[0,1])")               # bottom→top
        p3 = parse_path(f"({hw}+{h}*i-({2*hw})*t)_(t:[0,1])")       # top→left
        p4 = parse_path(f"(-{hw}+{h}*i-({h})*i*t)_(t:[0,1])")      # left top→bottom
        self.assertTrue(rules.is_closed_contour([p1, p2, p3, p4]))

    def test_triangle_4_sides(self):
        """三角形 0→1→1+i→i→0 → 闭合。"""
        p1 = parse_path("(t)_(t:[0,1])")           # 0 → 1
        p2 = parse_path("(1+i*t)_(t:[0,1])")       # 1 → 1+i
        p3 = parse_path("(1+i-t)_(t:[0,1])")       # 1+i → i
        p4 = parse_path("(i*(1-t))_(t:[0,1])")     # i → 0
        self.assertTrue(rules.is_closed_contour([p1, p2, p3, p4]))

    def test_upper_semicircle_plus_line(self):
        """上半圆弧 + 实轴线段 → 闭合。"""
        line = parse_path("(2*t-1)_(t:[0,1])")        # -1 → 1
        arc  = parse_path("(exp(i*pi*t))_(t:[0,1])")  # 1 → -1 via upper
        self.assertTrue(rules.is_closed_contour([line, arc]))

    def test_lower_semicircle_plus_line(self):
        """下半圆弧 + 实轴线段 → 闭合。"""
        line = parse_path("(2*t-1)_(t:[0,1])")         # -1 → 1
        arc  = parse_path("(exp(-i*pi*t))_(t:[0,1])") # 1 → -1 via lower
        self.assertTrue(rules.is_closed_contour([line, arc]))

    def test_pentagon_closed(self):
        """五边形 0→1→1+i→0.5+i→0.5→0 → 闭合。"""
        p1 = parse_path("(t)_(t:[0,1])")               # 0 → 1
        p2 = parse_path("(1+i*t)_(t:[0,1])")            # 1 → 1+i
        p3 = parse_path("((1-0.5*t)+i)_(t:[0,1])")      # 1+i → 0.5+i
        p4 = parse_path("(0.5+i*(1-t))_(t:[0,1])")     # 0.5+i → 0.5
        p5 = parse_path("(0.5*(1-t))_(t:[0,1])")        # 0.5 → 0
        self.assertTrue(rules.is_closed_contour([p1, p2, p3, p4, p5]))

    def test_two_identical_circles(self):
        """两个相同闭合路径 → p0.end(1) == p1.start(1) → 闭合。"""
        p = parse_path("(exp(i*t))_(t:[0,2*pi])")
        self.assertTrue(rules.is_closed_contour([p, p]))

    def test_three_segment_chain_closed(self):
        """三段折线：0→1→2→0 → 闭合。"""
        p1 = parse_path("(t)_(t:[0,1])")       # 0 → 1
        p2 = parse_path("(1+t)_(t:[0,1])")      # 1 → 2
        p3 = parse_path("(2-t)_(t:[0,1])")      # 2 → 1  ← wait, 2-1=1 not 0
        # Correct: p3 should go 2 → 0
        p3 = parse_path("(2-2*t)_(t:[0,1])")    # 2 → 0
        self.assertTrue(rules.is_closed_contour([p1, p2, p3]))


# ---------------------------------------------------------------------------
# 5. 多路径 — 不闭合
# ---------------------------------------------------------------------------

class TestIsClosedContourMultiOpen(unittest.TestCase):
    """多段路径不形成闭合围道。"""

    def setUp(self):
        clear_caches()

    def test_rectangle_one_edge_missing(self):
        """矩形缺一条边 → 不闭合。"""
        hw = 1.0
        h = 2.0
        p1 = parse_path(f"(-{hw}+{2*hw}*t)_(t:[0,1])")
        p2 = parse_path(f"({hw}+{h}*i*t)_(t:[0,1])")
        p3 = parse_path(f"({hw}+{h}*i-({2*hw})*t)_(t:[0,1])")
        self.assertFalse(rules.is_closed_contour([p1, p2, p3]))

    def test_rectangle_edge_order_reversed(self):
        """矩形一条边顺序颠倒 → 链断裂 → 不闭合。"""
        hw = 1.0
        h = 2.0
        p1 = parse_path(f"(-{hw}+{2*hw}*t)_(t:[0,1])")
        p2 = parse_path(f"({hw}+{h}*i*t)_(t:[0,1])")
        p3 = parse_path(f"({hw}+{h}*i-({2*hw})*t)_(t:[0,1])")
        p4 = parse_path(f"(-{hw}+{h}*i-({h})*i*t)_(t:[0,1])")
        # swap p3 and p4 order
        self.assertFalse(rules.is_closed_contour([p1, p2, p4, p3]))

    def test_two_open_segments_not_chain(self):
        """两段开口线段不首尾相连 → 不闭合。"""
        p1 = parse_path("(2*t-1)_(t:[0,1])")   # -1 → 1
        p2 = parse_path("(1+t)_(t:[0,1])")     # 1 → 2
        self.assertFalse(rules.is_closed_contour([p1, p2]))

    def test_triangle_one_side_gap(self):
        """三角形缺一段 → 不闭合。"""
        p1 = parse_path("(t)_(t:[0,1])")
        p2 = parse_path("(1+i*t)_(t:[0,1])")
        self.assertFalse(rules.is_closed_contour([p1, p2]))

    def test_arc_wrong_direction_no_connection(self):
        """两个同向半圆弧无法首尾相连 → 不闭合。"""
        arc1 = parse_path("(exp(i*pi*t))_(t:[0,1])")    # 1 → -1
        arc2 = parse_path("(exp(-i*pi*t))_(t:[0,1])")   # 1 → -1 (同向)
        self.assertFalse(rules.is_closed_contour([arc1, arc2]))

    def test_three_segment_broken_chain(self):
        """三段链：0→1→2→1（末不回到首）→ 不闭合。"""
        p1 = parse_path("(t)_(t:[0,1])")        # 0 → 1
        p2 = parse_path("(1+t)_(t:[0,1])")       # 1 → 2
        p3 = parse_path("(2-t)_(t:[0,1])")       # 2 → 1  (回到 1，不回到 0)
        self.assertFalse(rules.is_closed_contour([p1, p2, p3]))

    def test_two_same_closed_paths_not_closed(self):
        """三段链：[0→1, 1→1(圆), 1→2]。最后一段 1→2 无法回到起点 → 不闭合。"""
        p1 = parse_path("(t)_(t:[0,1])")           # 0 → 1
        p2 = parse_path("(exp(i*t))_(t:[0,2*pi])")  # 1 → 1（闭合圆）
        p3 = parse_path("(2*t)_(t:[0,1])")          # 0 → 2
        # p2.end=1, p3.start=0 → 1≠0 → 不闭合
        self.assertFalse(rules.is_closed_contour([p1, p2, p3]))

    def test_single_open_arc(self):
        """单个开口弧 → 不闭合。"""
        path = parse_path("(exp(i*t))_(t:[0,pi])")
        self.assertFalse(rules.is_closed_contour([path]))

    def test_single_closed_path(self):
        """单个闭合路径 → 闭合。"""
        path = parse_path("(exp(i*t))_(t:[0,2*pi])")
        self.assertTrue(rules.is_closed_contour([path]))


# ---------------------------------------------------------------------------
# 6. 鲁棒性：确定性、异常路径、特殊数值
# ---------------------------------------------------------------------------

class TestIsClosedContourRobustness(unittest.TestCase):
    """确定性、异常路径 fallback、特殊数值边界。"""

    def setUp(self):
        clear_caches()

    def test_repeated_calls_deterministic_closed(self):
        """同一闭合路径反复调用，结果必须一致。"""
        path = parse_path("(exp(i*t))_(t:[0,2*pi])")
        for _ in range(5):
            self.assertTrue(rules.is_closed_contour([path]))

    def test_repeated_calls_deterministic_open(self):
        """同一开口路径反复调用，结果必须一致。"""
        path = parse_path("(exp(i*t))_(t:[0,pi])")
        for _ in range(5):
            self.assertFalse(rules.is_closed_contour([path]))

    def test_float_precision_closed(self):
        """浮点半径：数值精度边界，闭合性判断必须正确。"""
        path = parse_path("(0.1*exp(i*t))_(t:[0,2*pi])")
        self.assertTrue(rules.is_closed_contour([path]))

    def test_negative_radius(self):
        """负半径：绝对值 |−r| 圆 → 闭合。"""
        path = parse_path("(-2*exp(i*t))_(t:[0,2*pi])")
        self.assertTrue(rules.is_closed_contour([path]))

    def test_complex_radius(self):
        """复数半径 |e^(iθ)·r| = r 圆 → 闭合。"""
        path = parse_path("(exp(i)*exp(i*t))_(t:[0,2*pi])")
        self.assertTrue(rules.is_closed_contour([path]))

    def test_npi_angle_full_circle(self):
        """n 圈：t:[0,n*2*pi] → 闭合。"""
        for n in [1, 2, 3]:
            path = parse_path(f"(exp(i*t))_(t:[0,{n}*2*pi])")
            self.assertTrue(rules.is_closed_contour([path]))

    def test_cos_sin_arc_open(self):
        """cos/sin 弧显式 0→1 → start=1, end=cos(1)+i*sin(1) ≠ 1 → 不闭合。"""
        pi_expr = Var("pi")
        real_part = Fun("cos", Op("*", pi_expr, Var("t")))
        imag_part = Fun("*", Fun("i"), Fun("sin", Op("*", pi_expr, Var("t"))))
        path_expr = Op("+", real_part, imag_part)
        path = CINTPath("t", path_expr, Const(0), Const(1))
        result = rules.is_closed_contour([path])
        # 端点不化简为 1，应返回 False
        self.assertFalse(result)

    def test_zero_radius_path(self):
        """退化路径 r=0 → 始终在原点 → 闭合。"""
        path = parse_path("(0*exp(i*t))_(t:[0,2*pi])")
        self.assertTrue(rules.is_closed_contour([path]))

    def test_all_paths_trigger_compute_exception(self):
        """所有路径在 compute_path_endpoints 中抛异常 → 返回 False。"""
        # 无法构造这种情况（无路径会走 endpoints=[] 分支）
        # 验证空 endpoints 返回 False（需修改实现才能期望 False）
        # 当前实现返回 True — 此测试记录为预期 False
        self.assertFalse(rules.is_closed_contour([None, "bad", 1]))


# ---------------------------------------------------------------------------
# 7. 复杂路径（基于复分析标准围道：keyhole、indented、wedge、annulus、ellipse、sector、dumbbell）
# ---------------------------------------------------------------------------

class TestIsClosedContourComplexPaths(unittest.TestCase):
    """
    基于网络资料中标准复分析围道的闭合性测试。
    参考资料：
      - Keyhole contour: Math.StackExchange "Complex keyhole contour integral"
      - Indented contour: Math.SE "Indented Contour Integrals", pole at origin
      - Wedge contour: Math.SE "Wedge contour parameterisation", angle pi/3
      - Annulus: two concentric circles with opposite orientation (outer CCW, inner CW)
      - Ellipse: a*cos(t) + i*b*sin(t), t:[0,2*pi]
      - Sector: quarter-circle + two radii
      - Dumbbell: two circles joined by straight-line bridges
    """

    def setUp(self):
        clear_caches()

    # --- Keyhole contour ---
    # 由4段组成：外弧(CCW) + 径向线(下行) + 内弧(CW) + 径向线(上行)
    # 绕过支割 (branch cut on positive real axis)
    def test_keyhole_contour_closed(self):
        """Keyhole 围道：绕过支割的4段路径 → 闭合。"""
        R, delta = 2.0, 0.1
        outer  = parser.parse_expr(f"(R*exp(i*((2*pi-2*{delta})*t+{delta})))_(t:[0,1])")
        rad_in  = parser.parse_expr(f"(R*exp(i*(2*pi-{delta}))*(1-t)+delta*exp(i*(2*pi-{delta}))*t)_(t:[0,1])")
        inner  = parser.parse_expr(f"(delta*exp(i*(2*pi-{delta}-(2*pi-2*{delta})*t)))_(t:[0,1])")
        rad_out = parser.parse_expr(f"(delta*exp(i*{delta})*(1-t)+R*exp(i*{delta})*t)_(t:[0,1])")
        self.assertTrue(rules.is_closed_contour([outer, rad_in, inner, rad_out]))

    def test_keyhole_outer_arc_only_open(self):
        """仅外弧 → 不闭合（keyhole 需4段）。"""
        R, delta = 2.0, 0.1
        outer = parser.parse_expr(f"(R*exp(i*((2*pi-2*{delta})*t+{delta})))_(t:[0,1])")
        self.assertFalse(rules.is_closed_contour([outer]))

    # --- Indented contour ---
    # 在实轴极点处做小半圆 detour，绕过 z=0
    # 4段：上半大弧(-R→R) + 上实轴(R→δ) + 上半小弧(δ→-δ,CW) + 下实轴(-δ→-R)
    def test_indented_contour_closed(self):
        """Indented 围道：绕过实轴极点的上半 indented contour → 闭合。"""
        R = 3.0
        delta = 0.2
        large_arc = parser.parse_expr("(3*exp(i*(pi-pi*t)))_(t:[0,1])")
        upper_seg  = parser.parse_expr(f"(3-(3+{delta})*t)_(t:[0,1])")
        indent_arc = parser.parse_expr(f"({delta}*exp(i*(pi-pi*t)))_(t:[0,1])")
        lower_seg  = parser.parse_expr(f"({delta}+(-3-{delta})*t)_(t:[0,1])")
        self.assertTrue(rules.is_closed_contour([large_arc, upper_seg, indent_arc, lower_seg]))

    def test_indented_missing_indent_open(self):
        """Indented 缺小半圆段 → 不闭合。"""
        R = 3.0
        delta = 0.2
        large_arc = parser.parse_expr("(3*exp(i*(pi-pi*t)))_(t:[0,1])")
        upper_seg  = parser.parse_expr(f"(3-(3+{delta})*t)_(t:[0,1])")
        lower_seg  = parser.parse_expr(f"({delta}+(-3-{delta})*t)_(t:[0,1])")
        self.assertFalse(rules.is_closed_contour([large_arc, upper_seg, lower_seg]))

    # --- Wedge contour (pi/3) ---
    # 三边：弧(0→R*exp(i*pi/3)) + 斜边(R*exp(i*pi/3)→0) + 实轴(0→R)
    def test_wedge_contour_pi_over_3_closed(self):
        """Wedge 围道，角度 pi/3 → 闭合。"""
        R = 2.0
        wedge_arc  = parser.parse_expr("(2*exp(i*pi/3*t))_(t:[0,1])")          # 0 -> 2*exp(i*pi/3)
        wedge_ang  = parser.parse_expr("(2*exp(i*pi/3)*(1-t))_(t:[0,1])")      # 2*exp(i*pi/3) -> 0
        wedge_real = parser.parse_expr("(2*t)_(t:[0,1])")                      # 0 -> 2
        self.assertTrue(rules.is_closed_contour([wedge_arc, wedge_ang, wedge_real]))

    def test_wedge_half_angle_closed(self):
        """Wedge 角度 pi/6（任意正角度）→ 闭合（扇形三角形：0 → Re^(iπ/6) → 0 → R）。"""
        R = 2.0
        a = "pi/6"
        wedge_arc  = parser.parse_expr(f"(2*exp(i*{a}*t))_(t:[0,1])")
        wedge_ang  = parser.parse_expr(f"(2*exp(i*{a})*(1-t))_(t:[0,1])")
        wedge_real = parser.parse_expr("(2*t)_(t:[0,1])")
        self.assertTrue(rules.is_closed_contour([wedge_arc, wedge_ang, wedge_real]))

    def test_wedge_wrong_order_open(self):
        """Wedge 边顺序错误 → 不闭合。"""
        R = 2.0
        wedge_arc  = parser.parse_expr("(2*exp(i*pi/3*t))_(t:[0,1])")
        wedge_real = parser.parse_expr("(2*t)_(t:[0,1])")
        # 顺序错了：arc.start(2) != real.start(0)
        self.assertFalse(rules.is_closed_contour([wedge_arc, wedge_real]))

    # --- Annulus ---
    # 外圆 CCW + 径向线(外→内) + 内圆 CW + 径向线(内→外)
    # 积分时内外方向相反使围道等效于两圆之差
    def test_annulus_closed(self):
        """Annulus 围道：外圆+径向线+内圆+径向线 → 闭合。"""
        outer  = parser.parse_expr("(2*exp(i*2*pi*t))_(t:[0,1])")          # 2 -> 2
        bridge_out_in = parser.parse_expr("(2-(2-1)*t)_(t:[0,1])")             # 2 -> 1
        inner  = parser.parse_expr("(1*exp(-i*2*pi*t))_(t:[0,1])")            # 1 -> 1 CW
        bridge_in_out = parser.parse_expr("(1+(2-1)*t)_(t:[0,1])")             # 1 -> 2
        self.assertTrue(rules.is_closed_contour([outer, bridge_out_in, inner, bridge_in_out]))

    def test_annulus_without_inner_circle_closed(self):
        """Annulus 缺内圆时，outer 圆 + 两径向线形成闭合三角形（0→1→2→0）→ 数学上闭合。"""
        outer = parser.parse_expr("(2*exp(i*2*pi*t))_(t:[0,1])")
        bridge_out_in = parser.parse_expr("(2-(2-1)*t)_(t:[0,1])")
        bridge_in_out = parser.parse_expr("(1+(2-1)*t)_(t:[0,1])")
        # outer.end(2)==b1.start(2), b1.end(1)==b2.start(1), b2.end(2)==outer.start(2) → closed
        self.assertTrue(rules.is_closed_contour([outer, bridge_out_in, bridge_in_out]))

    # --- Ellipse ---
    # z = a*cos(t) + i*b*sin(t), t:[0,2*pi]（用 2*pi*t 参数化）
    def test_ellipse_standard_closed(self):
        """标准椭圆 a*cos + i*b*sin，t:[0,2*pi] → 闭合。"""
        ellipse = parser.parse_expr("(3*cos(2*pi*t)+i*2*sin(2*pi*t))_(t:[0,1])")
        self.assertTrue(rules.is_closed_contour([ellipse]))

    def test_ellipse_unit_circle_as_special_case(self):
        """单位圆是椭圆 a=b=1 的特例 → 闭合。"""
        circle = parser.parse_expr("(cos(2*pi*t)+i*sin(2*pi*t))_(t:[0,1])")
        self.assertTrue(rules.is_closed_contour([circle]))

    def test_ellipse_half_turn_open(self):
        """半椭圆（t:[0,pi]）→ 不闭合。"""
        half_ellipse = parser.parse_expr("(3*cos(pi*t)+i*2*sin(pi*t))_(t:[0,1])")
        self.assertFalse(rules.is_closed_contour([half_ellipse]))

    # --- Sector (quarter-circle) ---
    # 两半径 + 1/4 圆弧，角度 pi/2
    def test_sector_quarter_circle_closed(self):
        """Sector 围道：两半径 + 1/4 圆弧 → 闭合。"""
        R = 2.0
        rad1 = parser.parse_expr("(2*t)_(t:[0,1])")
        arc  = parser.parse_expr("(2*exp(i*pi/2*t))_(t:[0,1])")
        rad2 = parser.parse_expr("(2*i*(1-t))_(t:[0,1])")
        self.assertTrue(rules.is_closed_contour([rad1, arc, rad2]))

    def test_sector_full_circle_open(self):
        """Sector 缺一条边 → 不闭合。"""
        R = 2.0
        rad1 = parser.parse_expr("(2*t)_(t:[0,1])")
        arc  = parser.parse_expr("(2*exp(i*pi/2*t))_(t:[0,1])")
        self.assertFalse(rules.is_closed_contour([rad1, arc]))

    # --- Dumbbell ---
    # 两圆 + 两条直线桥接
    def test_dumbbell_two_circles_bridged_closed(self):
        """Dumbbell 围道：两闭合圆 + 桥接线段 → 闭合。"""
        # c1: 中心 -1, r=0.5, 起点/终点 = -1+0.5 = -0.5
        c1      = parser.parse_expr("(-1+0.5*exp(i*2*pi*t))_(t:[0,1])")
        bridge1 = parser.parse_expr("(-0.5+2*t)_(t:[0,1])")    # -0.5 -> 1.5
        # c2: 中心 1, r=0.5, 起点/终点 = 1+0.5 = 1.5
        c2      = parser.parse_expr("(1+0.5*exp(i*2*pi*t))_(t:[0,1])")
        bridge2 = parser.parse_expr("(1.5-2*t)_(t:[0,1])")     # 1.5 -> -0.5
        self.assertTrue(rules.is_closed_contour([c1, bridge1, c2, bridge2]))

    def test_dumbbell_missing_bridge_open(self):
        """Dumbbell 缺一座桥 → 不闭合。"""
        c1      = parser.parse_expr("(-1+0.5*exp(i*2*pi*t))_(t:[0,1])")
        bridge1 = parser.parse_expr("(-0.5+2*t)_(t:[0,1])")
        c2      = parser.parse_expr("(1+0.5*exp(i*2*pi*t))_(t:[0,1])")
        # bridge2 missing
        self.assertFalse(rules.is_closed_contour([c1, bridge1, c2]))

    def test_dumbbell_missing_both_bridges_open(self):
        """Dumbbell 两桥均缺失 → 不闭合。"""
        c1 = parser.parse_expr("(-1+0.5*exp(i*2*pi*t))_(t:[0,1])")
        c2 = parser.parse_expr("(1+0.5*exp(i*2*pi*t))_(t:[0,1])")
        self.assertFalse(rules.is_closed_contour([c1, c2]))

    # --- Double circle with same center, opposite orientation ---
    # 外圆 CCW + 内圆 CW，两个都闭合但方向相反
    def test_double_circle_opposite_orientation(self):
        """两个独立圆（无桥接）：outer.end(3)≠inner.start(1) → 不闭合。"""
        outer = parser.parse_expr("(3*exp(i*2*pi*t))_(t:[0,1])")
        inner = parser.parse_expr("(1*exp(-i*2*pi*t))_(t:[0,1])")
        self.assertFalse(rules.is_closed_contour([outer, inner]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
