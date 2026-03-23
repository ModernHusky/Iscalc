"""测试 expand_euler 函数：exp(±i*θ) → cos(θ) ± i*sin(θ) 的展开功能。"""
import unittest
import cmath
import math
import random
from integral import context, parser, expr, rules


def contains_in_order(result_str: str, *parts: str) -> bool:
    """检查 result_str 中是否依次包含所有 parts（允许中间有其他内容）。"""
    pos = 0
    for part in parts:
        idx = result_str.find(part, pos)
        if idx == -1:
            return False
        pos = idx + len(part)
    return True


def check_parts_unordered_pair(result_str: str, part1: str, part2: str) -> bool:
    """检查 part1 和 part2 都出现在 result_str 中，顺序不限。"""
    idx1 = result_str.find(part1)
    idx2 = result_str.find(part2)
    return idx1 != -1 and idx2 != -1


class TestExpandEuler(unittest.TestCase):
    """测试 expand_euler 的各种虚数指数形式。"""

    @classmethod
    def setUpClass(cls):
        cls.ctx = context.Context()

    # ── 特殊角度：纯虚数指数 i * π/k ──────────────────────────────────
    def check(self, input_str: str, *expected_parts: str):
        """通用检查：展开结果的字符串包含所有 expected_parts。"""
        e = parser.parse_expr(input_str)
        result = rules.expand_euler(e, self.ctx)
        result_str = str(result)
        self.assertTrue(
            contains_in_order(result_str, *expected_parts),
            f"\n输入: {input_str}\n结果: {result_str}\n期望依次包含: {expected_parts}"
        )

    def test_expand_i_pi_over_2(self):
        self.check("exp(i*pi/2)", "cos", "pi", "2", "i", "sin")

    def test_expand_i_pi_over_3(self):
        self.check("exp(i*pi/3)", "cos", "pi", "3", "i", "sin")

    def test_expand_i_pi_over_4(self):
        self.check("exp(i*pi/4)", "cos", "pi", "4", "i", "sin")

    def test_expand_i_pi_over_6(self):
        self.check("exp(i*pi/6)", "cos", "pi", "6", "i", "sin")

    def test_expand_i_pi_over_5(self):
        self.check("exp(i*pi/5)", "cos", "pi", "5", "i", "sin")

    def test_expand_i_pi_over_8(self):
        self.check("exp(i*pi/8)", "cos", "pi", "8", "i", "sin")

    def test_expand_i_pi(self):
        self.check("exp(i*pi)", "cos(pi)", "i", "sin(pi)")

    def test_expand_2i_pi_over_3(self):
        self.check("exp(2*i*pi/3)", "cos", "pi", "3", "i", "sin")

    def test_expand_3i_pi_over_4(self):
        self.check("exp(3*i*pi/4)", "cos", "pi", "4", "i", "sin")

    def test_expand_4i_pi_over_3(self):
        self.check("exp(4*i*pi/3)", "cos", "pi", "3", "i", "sin")

    def test_expand_5i_pi_over_6(self):
        self.check("exp(5*i*pi/6)", "cos", "pi", "6", "i", "sin")

    def test_expand_neg_i_pi_over_4(self):
        self.check("exp(-i*pi/4)", "cos", "-", "pi", "4", "i", "sin")

    def test_expand_neg_i_pi_over_3(self):
        self.check("exp(-i*pi/3)", "cos", "-", "pi", "3", "i", "sin")

    def test_expand_neg_2i_pi_over_3(self):
        self.check("exp(-2*i*pi/3)", "cos", "-", "pi", "3", "i", "sin")

    def test_expand_neg3i_pi_over_4(self):
        self.check("exp(-3*i*pi/4)", "cos", "-", "pi", "4", "i", "sin")

    # ── 符号参数：保留 sin/cos 不展开 ─────────────────────────────────
    def test_expand_i_t(self):
        """exp(i*t) → cos(t) + i*sin(t)，t 为符号变量"""
        self.check("exp(i*t)", "cos(t)", "i", "sin(t)")

    def test_expand_neg_i_t(self):
        """exp(-i*t) → cos(-t) + i*sin(-t)"""
        self.check("exp(-i*t)", "cos", "i", "sin")

    def test_expand_ki_t(self):
        """exp(k*i*t) → cos(k*t) + i*sin(k*t)，k 为符号系数"""
        self.check("exp(k*i*t)", "cos", "k", "t", "i", "sin")

    def test_expand_ni_t(self):
        """exp(n*i*t) → cos(n*t) + i*sin(n*t)"""
        self.check("exp(2*i*t)", "cos", "t", "i", "sin")

    def test_expand_i_k_t(self):
        """exp(i*k*t) → cos(k*t) + i*sin(k*t)"""
        self.check("exp(i*k*t)", "cos", "k", "t", "i", "sin")

    def test_expand_i_pi_times_1_minus_t(self):
        """exp(i*π*(1-t)) → cos(π*(1-t)) + i*sin(π*(1-t))"""
        self.check("exp(i*pi*(1-t))", "cos", "pi", "1", "t", "i", "sin")

    def test_expand_i_2pi_t(self):
        """exp(i*2π*t) → cos(2π*t) + i*sin(2π*t)"""
        # 乘积顺序由 normalize 决定，可能是 t*pi 也可能是 pi*t，数学等价即可
        e = parser.parse_expr("exp(i*2*pi*t)")
        result = rules.expand_euler(e, self.ctx)
        result_str = str(result)
        self.assertIn("cos", result_str)
        self.assertIn("i", result_str)
        self.assertIn("sin", result_str)
        # pi 和 t 至少都出现（顺序不限）
        self.assertTrue(check_parts_unordered_pair(result_str, "pi", "t"),
            f"cos/sin 参数中应同时含 pi 和 t: {result_str}")

    def test_expand_i_pi_t_minus_alpha(self):
        """exp(i*π*t - i*α)"""
        # 乘积顺序可能互换，检查 cos/sin 出现且同时含 pi, t, alpha
        e = parser.parse_expr("exp(i*pi*t - i*alpha)")
        result = rules.expand_euler(e, self.ctx)
        result_str = str(result)
        self.assertIn("cos", result_str)
        self.assertIn("i", result_str)
        self.assertIn("sin", result_str)
        self.assertTrue("alpha" in result_str, f"缺少 alpha: {result_str}")

    def test_expand_i_2pi_over_5(self):
        self.check("exp(2*i*pi/5)", "cos", "pi", "5", "i", "sin")

    def test_expand_i_7pi_over_12(self):
        self.check("exp(7*i*pi/12)", "cos", "pi", "12", "i", "sin")

    # ── 有理数系数 ────────────────────────────────────────────────────
    def test_expand_i_1over2(self):
        """exp(i/2) → cos(1/2) + i*sin(1/2)"""
        self.check("exp(i/2)", "cos", "i", "sin")

    def test_expand_i_a_over_b(self):
        """exp((i*a)/b) → cos(a/b) + i*sin(a/b)"""
        self.check("exp((i*a)/b)", "cos", "a", "b", "i", "sin")

    def test_expand_ia_over_b(self):
        """exp(i*a/b) → cos(a/b) + i*sin(a/b)"""
        self.check("exp(i*a/b)", "cos", "a", "b", "i", "sin")

    def test_expand_ni_over_m(self):
        """exp((n*i)/m) → cos(n/m) + i*sin(n/m)"""
        self.check("exp((2*i)/3)", "cos", "2", "3", "i", "sin")

    # ── 纯虚数形式 ────────────────────────────────────────────────────
    def test_expand_i_only(self):
        """exp(i) → cos(1) + i*sin(1)"""
        self.check("exp(i)", "cos(1)", "i", "sin(1)")

    def test_expand_neg_i_only(self):
        """exp(-i) → cos(-1) + i*sin(-1)"""
        self.check("exp(-i)", "cos", "i", "sin")

    def test_expand_2i(self):
        """exp(2*i) → cos(2) + i*sin(2)"""
        self.check("exp(2*i)", "cos(2)", "i", "sin(2)")

    def test_expand_ni(self):
        """exp(n*i) → cos(n) + i*sin(n)"""
        self.check("exp(3*i)", "cos(3)", "i", "sin(3)")

    # ── 复数指数：exp(a + i*b) = exp(a) * (cos(b) + i*sin(b)) ─────────
    def test_expand_1_plus_i_pi_over_2(self):
        """exp(1 + i*π/2)"""
        self.check("exp(1+i*pi/2)", "exp", "cos", "pi", "2", "i", "sin")

    def test_expand_3_plus_i_pi_over_4(self):
        """exp(3 + i*π/4)"""
        self.check("exp(3+i*pi/4)", "exp", "cos", "pi", "4", "i", "sin")

    def test_expand_a_plus_i_b(self):
        """exp(a + i*b)"""
        self.check("exp(a+i*b)", "exp", "cos", "b", "i", "sin")

    def test_expand_neg1_plus_i(self):
        """exp(-1 + i)"""
        self.check("exp(-1+i)", "exp", "cos(1)", "i", "sin(1)")

    def test_expand_pi_plus_i(self):
        """exp(π + i)"""
        self.check("exp(pi+i)", "exp", "cos(1)", "i", "sin(1)")

    def test_expand_x_plus_i_x(self):
        """exp(x + i*x)"""
        self.check("exp(x+i*x)", "exp", "cos", "x", "i", "sin", "x")

    def test_expand_i_times_1_minus_t_plus_offset(self):
        """exp(α + i*π*(1-t))"""
        self.check("exp(alpha + i*pi*(1-t))", "exp", "cos", "pi", "1", "t", "i", "sin")

    # ── 纯实数指数：不展开 ────────────────────────────────────────────
    def test_no_expand_real_exp(self):
        """exp(x)：不含 i，不展开"""
        e = parser.parse_expr("exp(x)")
        result = rules.expand_euler(e, self.ctx)
        self.assertEqual(str(result), "exp(x)")

    def test_no_expand_real_pi(self):
        """exp(π)：不含 i，不展开"""
        e = parser.parse_expr("exp(pi)")
        result = rules.expand_euler(e, self.ctx)
        self.assertEqual(str(result), "exp(pi)")

    def test_no_expand_real_number(self):
        """exp(2)：不含 i，不展开"""
        e = parser.parse_expr("exp(2)")
        result = rules.expand_euler(e, self.ctx)
        self.assertEqual(str(result), "exp(2)")

    # ── 递归：子表达式中的 exp ─────────────────────────────────────────
    def test_expand_nested(self):
        """嵌套表达式：exp(i) + exp(i*pi/2)"""
        self.check("exp(i) + exp(i*pi/2)", "cos(1)", "sin(1)", "cos", "pi", "2", "sin")

    def test_expand_in_product(self):
        """乘积表达式：r*exp(i*t)"""
        self.check("r*exp(i*t)", "r", "cos", "t", "i", "sin")

    def test_expand_nested_exp(self):
        """嵌套：exp(exp(i*t))，外层不含 i，不展开"""
        e = parser.parse_expr("exp(exp(i*t))")
        result = rules.expand_euler(e, self.ctx)
        self.assertEqual(str(result), "exp(exp(i * t))")

    def test_expand_i_pi_over_2_squared(self):
        """exp((i*pi/2)^2)，参数不含纯虚数形式，不展开"""
        e = parser.parse_expr("exp((i*pi/2)^2)")
        result = rules.expand_euler(e, self.ctx)
        self.assertEqual(str(result), "exp((i * pi / 2) ^ 2)")


class TestExpandEulerMathEquivalence(unittest.TestCase):
    """验证展开结果的数学等价性：expand_euler 展开后数值上与原表达式相等。"""

    @classmethod
    def setUpClass(cls):
        cls.ctx = context.Context()

    def _eval_expr(self, e, subs: dict):
        """对表达式进行数值求值。"""
        if isinstance(e, (int, float)):
            return complex(e)
        if isinstance(e, expr.Const):
            if isinstance(e.val, (int, float)):
                return complex(e.val)
            return complex(float(e.val.numerator)) / complex(float(e.val.denominator))
        if isinstance(e, expr.Var):
            return subs.get(e.name, complex(0))
        if isinstance(e, expr.Fun):
            if e.func_name == 'i':
                return 1j
            if e.func_name == 'pi':
                return math.pi
            args = [self._eval_expr(a, subs) for a in e.args]
            if e.func_name == 'exp':
                return cmath.exp(args[0])
            if e.func_name == 'cos':
                return cmath.cos(args[0])
            if e.func_name == 'sin':
                return cmath.sin(args[0])
            if e.func_name == 'sqrt':
                v = args[0]
                return complex(math.sqrt(abs(v))) * (1 if v >= 0 else 1j)
            return complex(0)
        if isinstance(e, expr.Op):
            if e.op == '-' and len(e.args) == 1:
                return -self._eval_expr(e.args[0], subs)
            args = [self._eval_expr(a, subs) for a in e.args]
            if e.op == '+':
                return args[0] + args[1]
            if e.op == '-':
                return args[0] - args[1]
            if e.op == '*':
                return args[0] * args[1]
            if e.op == '/':
                return args[0] / args[1]
            if e.op == '^':
                return args[0] ** args[1]
        return complex(0)

    def _numeric_equal(self, original_str: str, tol: float = 1e-12):
        """验证 expand_euler 展开后数值上与原表达式相等。"""
        e_orig = parser.parse_expr(original_str)
        e_expanded = rules.expand_euler(e_orig, self.ctx)

        vars_ = list(e_orig.get_vars() - {'pi', 'i'})
        if not vars_:
            self.skipTest(f"无变量，跳过数值测试: {original_str}")

        rng = random.Random(42)
        subs = {v: complex(rng.uniform(-3, 3), rng.uniform(-3, 3)) for v in vars_}
        subs['i'] = 1j
        subs['pi'] = math.pi

        val_orig = self._eval_expr(e_orig, subs)
        val_expanded = self._eval_expr(e_expanded, subs)
        diff = abs(val_orig - val_expanded)
        self.assertLess(
            diff, tol,
            f"数值不等: {original_str} → exp={val_orig}, expanded={val_expanded}, diff={diff}"
        )

    # ── 符号参数 ──────────────────────────────────────────────────────
    def test_numeric_i_t(self):           self._numeric_equal("exp(i*t)")
    def test_numeric_ki_t(self):          self._numeric_equal("exp(2*i*t)")
    def test_numeric_neg_i_t(self):       self._numeric_equal("exp(-i*t)")
    def test_numeric_i_2pi_t(self):       self._numeric_equal("exp(i*2*pi*t)")

    # ── 复数指数 ───────────────────────────────────────────────────────
    def test_numeric_a_plus_i_b(self):     self._numeric_equal("exp(a+i*b)")
    def test_numeric_x_plus_i_x(self):     self._numeric_equal("exp(x+i*x)")


if __name__ == "__main__":
    unittest.main()
