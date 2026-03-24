"""Unit test for cintegrals using internal language."""

import unittest
import sys
import cProfile
import pstats
import os
from pathlib import Path

from integral import compstate
from integral import state
from integral import parser
from integral import context

# 获取项目根目录（tests 目录的父目录）
PROJECT_ROOT = Path(__file__).parent.parent
os.chdir(PROJECT_ROOT)


class ActionTest(unittest.TestCase):
    def setUp(self):
        """Clear all global caches before each test to prevent interference"""
        # Clear normalize cache
        from integral import poly
        poly._normalize_cache.clear()

        # Clear condition check cache
        from integral import condprover
        condprover.clear_condition_cache()

        # Clear poles and winding caches
        from integral import rules
        rules._poles_cache.clear()
        rules._winding_cache.clear()

    def check_file(self, filename: str, *, print_lines=False, print_state=False,
                   write_stats=False):
        with open(f'theories/{filename}.thy', 'r', encoding='utf-8') as f:
            content = f.read()
        state.check_actions(content, print_lines=print_lines, print_state=print_state,
                            write_stats=write_stats, filename=filename)

    def testCalculationFinished2(self):
        ctx = context.Context()
        ctx.load_book("base")
        st = state.InitialState(ctx)

        actions = """
            prove (INT x:[-oo,oo]. 1/(x^2+1)) = pi
            define C(t,r) = (r*exp(i*pi*(1-t)))_(t:[0,1])
            define L(t,r) = (r*(1-2*t))_(t:[0,1])
            subgoal 1: (LIM {r->oo}. CINT z:com(C(t,r),L(t,r)). 1/(z^2+1))=-pi
            lhs:
                apply residue theorem
            done
            subgoal 2: (LIM {r->oo}. CINT z:C(t,r). 1/(z^2+1))=0
            lhs:
                apply cintegral identity
                simplify
            done
            subgoal 3: (LIM {r->oo}. ((CINT z:com(C(t,r),L(t,r)). 1/(z^2+1)) - (CINT z:C(t,r). 1/(z^2+1)))) = -(INT x:[-oo,oo]. 1/(x^2+1))
            lhs:
                simplify
                rewrite (CINT z:com(C(t,r),L(t,r)). 1 / (z ^ 2 + 1)) to (CINT z:C(t,r). 1 / (z ^ 2 + 1)) + (CINT z:L(t,r). 1 / (z ^ 2 + 1))
                simplify
                rewrite (CINT z:L(t,r). 1 / (z ^ 2 + 1)) to (-(INT x:[-r,r]. 1/(x^2+1)))
                rewrite to -(INT x:[-oo,oo]. 1/(x^2+1))
            done
            lhs:
                apply 3 on (INT x:[-oo,oo]. 1/(x^2+1))
                apply 1 on (LIM {r -> oo}. CINT z:com(C(t,r),L(t,r)). 1 / (z ^ 2 + 1))
                apply 2 on (LIM {r -> oo}. CINT z:C(t,r). 1 / (z ^ 2 + 1))
                simplify
            done
        """
        actions = [s for s in actions.split('\n') if s.strip()]
        for act in actions:
            a = parser.parse_action(act)
            st = st.process_action(a)
        assert st.is_finished()


# ============ A类：标准有理函数 ============

def test_contour_single_pole_upper_half():
    """测试：∫_{-∞}^{+∞} dx/(x²+4) = π/2
    极点：z = 2i (一阶), z = -2i (一阶)
    围道：上半平面半圆
    留数：Res[f, z=2i] = lim_{z→2i} (z-2i)/(z²+4) = 1/(4i)
    预期结果：2πi × 1/(4i) = π/2
    """
    ctx = context.Context()
    ctx.load_book("base")
    st = state.InitialState(ctx)

    actions = """
        prove (INT x:[-oo,oo]. 1/(x^2+4)) = pi/2
        define C(t,r) = (r*exp(i*pi*(1-t)))_(t:[0,1])
        define L(t,r) = (r*(1-2*t))_(t:[0,1])
        subgoal 1: (LIM {r->oo}. CINT z:com(C(t,r),L(t,r)). 1/(z^2+4)) = -pi/2
        lhs:
            apply residue theorem
        done
        subgoal 2: (LIM {r->oo}. CINT z:C(t,r). 1/(z^2+4)) = 0
        lhs:
            apply cintegral identity
            simplify
        done
        subgoal 3: (LIM {r->oo}. ((CINT z:com(C(t,r),L(t,r)). 1/(z^2+4)) - (CINT z:C(t,r). 1/(z^2+4)))) = -(INT x:[-oo,oo]. 1/(x^2+4))
        lhs:
            simplify
            rewrite (CINT z:com(C(t,r),L(t,r)). 1/(z^2+4)) to (CINT z:C(t,r). 1/(z^2+4)) + (CINT z:L(t,r). 1/(z^2+4))
            simplify
            rewrite (CINT z:L(t,r). 1/(z^2+4)) to (-(INT x:[-r,r]. 1/(x^2+4)))
            rewrite to -(INT x:[-oo,oo]. 1/(x^2+4))
        done
        lhs:
            apply 3 on (INT x:[-oo,oo]. 1/(x^2+4))
            apply 1 on (LIM {r->oo}. CINT z:com(C(t,r),L(t,r)). 1/(z^2+4))
            apply 2 on (LIM {r->oo}. CINT z:C(t,r). 1/(z^2+4))
            simplify
        done
    """
    actions = [s for s in actions.split('\n') if s.strip()]
    for act in actions:
        a = parser.parse_action(act)
        st = st.process_action(a)
    assert st.is_finished()


def test_contour_two_poles_upper_half():
    """测试：∫_{-∞}^{+∞} dx/((x²+1)(x²+9)) = π/12
    极点：z = i, 3i, -i, -3i (全部一阶)
    上半平面：i, 3i
    留数：Res[f,z=i] = -i/16
          Res[f,z=3i] = i/48
    和：-i/16 + i/48 = -i/24
    预期结果：2πi × (-i/24) = 2π/24 = π/12
    """
    ctx = context.Context()
    ctx.load_book("base")
    st = state.InitialState(ctx)

    actions = """
        prove (INT x:[-oo,oo]. 1/((x^2+1)*(x^2+9))) = pi/12
        define C(t,r) = (r*exp(i*pi*(1-t)))_(t:[0,1])
        define L(t,r) = (r*(1-2*t))_(t:[0,1])
        subgoal 1: (LIM {r->oo}. CINT z:com(C(t,r),L(t,r)). 1/((z^2+1)*(z^2+9))) = -pi/12
        lhs:
            apply residue theorem
        done
        subgoal 2: (LIM {r->oo}. CINT z:C(t,r). 1/((z^2+1)*(z^2+9))) = 0
        lhs:
            apply cintegral identity
            simplify
        done
        subgoal 3: (LIM {r->oo}. ((CINT z:com(C(t,r),L(t,r)). 1/((z^2+1)*(z^2+9))) - (CINT z:C(t,r). 1/((z^2+1)*(z^2+9))))) = -(INT x:[-oo,oo]. 1/((x^2+1)*(x^2+9)))
        lhs:
            simplify
            rewrite (CINT z:com(C(t,r),L(t,r)). 1/((z^2+1)*(z^2+9))) to (CINT z:C(t,r). 1/((z^2+1)*(z^2+9))) + (CINT z:L(t,r). 1/((z^2+1)*(z^2+9)))
            simplify
            rewrite (CINT z:L(t,r). 1/((z^2+1)*(z^2+9))) to (-(INT x:[-r,r]. 1/((x^2+1)*(x^2+9))))
            rewrite to -(INT x:[-oo,oo]. 1/((x^2+1)*(x^2+9)))
        done
        lhs:
            apply 3 on (INT x:[-oo,oo]. 1/((x^2+1)*(x^2+9)))
            apply 1 on (LIM {r->oo}. CINT z:com(C(t,r),L(t,r)). 1/((z^2+1)*(z^2+9)))
            apply 2 on (LIM {r->oo}. CINT z:C(t,r). 1/((z^2+1)*(z^2+9)))
            simplify
        done
    """
    actions = [s for s in actions.split('\n') if s.strip()]
    for act in actions:
        a = parser.parse_action(act)
        st = st.process_action(a)
    assert st.is_finished()


def test_contour_conjugate_poles():
    """测试：∫_{-∞}^{+∞} dx/(x⁴+1) = π/√2
    极点：四次单位根 e^(iπ/4), e^(i3π/4), e^(i5π/4), e^(i7π/4)
    上半平面：e^(iπ/4) = (1+i)/√2, e^(i3π/4) = (-1+i)/√2
    留数计算较复杂，标准结果为 π/√2
    """
    ctx = context.Context()
    ctx.load_book("base")
    st = state.InitialState(ctx)

    actions = """
        prove (INT x:[-oo,oo]. 1/(x^4+1)) = pi/sqrt(2)
        define C(t,r) = (r*exp(i*pi*(1-t)))_(t:[0,1])
        define L(t,r) = (r*(1-2*t))_(t:[0,1])
        subgoal 1: (LIM {r->oo}. CINT z:com(C(t,r),L(t,r)). 1/(z^4+1)) = -pi/sqrt(2)
        lhs:
            apply residue theorem
            expand polynomial
            rewrite pi / (2 * sqrt(i)) to pi*sqrt(i)/(2*i)
            rewrite pi * sqrt(i) / (2 * i) to pi*sqrt(i)*i/(-2)
            simplify
            rewrite i to exp(i*pi/2)
            rewrite i to exp(i*pi/2) (at 2)
            simplify
            rewrite to pi/2*(exp(3 * i * pi / 4) - exp(i * pi / 4))
            rewrite exp(3 * i * pi / 4) to -sqrt(2)/2 + i*sqrt(2)/2
            rewrite exp(i * pi / 4) to sqrt(2)/2 + i*sqrt(2)/2
            simplify
            rewrite to -pi/sqrt(2)
        done
        subgoal 2: (LIM {r->oo}. CINT z:C(t,r). 1/(z^4+1)) = 0
        lhs:
            apply cintegral identity
            simplify
        done
        subgoal 3: (LIM {r->oo}. ((CINT z:com(C(t,r),L(t,r)). 1/(z^4+1)) - (CINT z:C(t,r). 1/(z^4+1)))) = -(INT x:[-oo,oo]. 1/(x^4+1))
        lhs:
            simplify
            rewrite (CINT z:com(C(t,r),L(t,r)). 1/(z^4+1)) to (CINT z:C(t,r). 1/(z^4+1)) + (CINT z:L(t,r). 1/(z^4+1))
            simplify
            rewrite (CINT z:L(t,r). 1/(z^4+1)) to (-(INT x:[-r,r]. 1/(x^4+1)))
            rewrite to -(INT x:[-oo,oo]. 1/(x^4+1))
        done
        lhs:
            apply 3 on (INT x:[-oo,oo]. 1/(x^4+1))
            apply 1 on (LIM {r->oo}. CINT z:com(C(t,r),L(t,r)). 1/(z^4+1))
            apply 2 on (LIM {r->oo}. CINT z:C(t,r). 1/(z^4+1))
            simplify
        done
    """
    actions = [s for s in actions.split('\n') if s.strip()]
    for act in actions:
        a = parser.parse_action(act)
        st = st.process_action(a)
    assert st.is_finished()


def test_contour_quadratic_discriminant():
    """测试：∫_{-∞}^{+∞} dx/(x²+x+1) = 2π/√3
    极点：z = (-1 ± i√3)/2，上半平面为 (-1+i√3)/2
    分母判别式：1-4 = -3 < 0
    留数：Res = 1/(2z+1) 在 z = (-1+i√3)/2 处
    """
    ctx = context.Context()
    ctx.load_book("base")
    st = state.InitialState(ctx)

    actions = """
        prove (INT x:[-oo,oo]. 1/(x^2+x+1)) = 2*pi/sqrt(3) for x^2+x+1!=0
        define C(t,r) = (r*exp(i*pi*(1-t)))_(t:[0,1])
        define L(t,r) = (r*(1-2*t))_(t:[0,1])
        subgoal 1: (LIM {r->oo}. CINT z:com(C(t,r),L(t,r)). 1/(z^2+z+1)) = -2*pi/sqrt(3)
        lhs:
            apply residue theorem
            rewrite to -2*pi/sqrt(3)
        done
        subgoal 2: (LIM {r->oo}. CINT z:C(t,r). 1/(z^2+z+1)) = 0
        lhs:
            apply cintegral identity
            simplify
        done
        subgoal 3: (LIM {r->oo}. ((CINT z:com(C(t,r),L(t,r)). 1/(z^2+z+1)) - (CINT z:C(t,r). 1/(z^2+z+1)))) = -(INT x:[-oo,oo]. 1/(x^2+x+1))
        lhs:
            simplify
            rewrite (CINT z:com(C(t,r),L(t,r)). 1/(z^2+z+1)) to (CINT z:C(t,r). 1/(z^2+z+1)) + (CINT z:L(t,r). 1/(z^2+z+1))
            simplify
            rewrite (CINT z:L(t,r). 1/(z^2+z+1)) to (-(INT x:[-r,r]. 1/(x^2+x+1)))
            rewrite to -(INT x:[-oo,oo]. 1/(x^2+x+1))
        done
        lhs:
            apply 3 on (INT x:[-oo,oo]. 1/(x^2+x+1))
            apply 1 on (LIM {r->oo}. CINT z:com(C(t,r),L(t,r)). 1/(z^2+z+1))
            apply 2 on (LIM {r->oo}. CINT z:C(t,r). 1/(z^2+z+1))
            simplify
        done
    """
    actions = [s for s in actions.split('\n') if s.strip()]
    for act in actions:
        a = parser.parse_action(act)
        st = st.process_action(a)
    assert st.is_finished()


# ============ B类：有理函数乘三角函数 ============

def test_contour_rational_cos_a1():
    """测试：∫_{-∞}^{+∞} cos(x)/(x²+1) dx = π/e
    构造：∫ e^(ix)/(x²+1) dx，取实部
    极点：z = i (上半平面), z = -i (下半平面)
    留数：Res[e^(iz)/(z²+1), z=i] = e^(i·i)/(2i) = e^(-1)/(2i)
    结果：2πi × e^(-1)/(2i) = π/e
    """
    ctx = context.Context()
    ctx.load_book("base")
    st = state.InitialState(ctx)

    actions = """
        prove (INT x:[-oo,oo]. exp(i*x)/(x^2+1)) = pi/exp(1)
        define C(t,r) = (r*exp(i*pi*(1-t)))_(t:[0,1])
        define L(t,r) = (r*(1-2*t))_(t:[0,1])
        subgoal 1: (LIM {r->oo}. CINT z:com(C(t,r),L(t,r)). exp(i*z)/(z^2+1)) = -pi/exp(1)
        lhs:
            apply residue theorem
            rewrite to -pi/exp(1)
        done
        subgoal 2: (LIM {r->oo}. CINT z:C(t,r). exp(i*z)/(z^2+1)) = 0
        lhs:
            apply cintegral identity
            substitute z for r * exp(i * pi * (1 - t))
            simplify
        done
        subgoal 3: (LIM {r->oo}. ((CINT z:com(C(t,r),L(t,r)). exp(i*z)/(z^2+1)) - (CINT z:C(t,r). exp(i*z)/(z^2+1)))) = -(INT x:[-oo,oo]. exp(i*x)/(x^2+1))
        lhs:
            simplify
            rewrite (CINT z:com(C(t,r),L(t,r)). exp(i*z)/(z^2+1)) to (CINT z:C(t,r). exp(i*z)/(z^2+1)) + (CINT z:L(t,r). exp(i*z)/(z^2+1))
            simplify
            rewrite (CINT z:L(t,r). exp(i*z)/(z^2+1)) to (-(INT x:[-r,r]. exp(i*x)/(x^2+1)))
            rewrite to -(INT x:[-oo,oo]. exp(i*x)/(x^2+1))
        done
        lhs:
            apply 3 on (INT x:[-oo,oo]. exp(i*x)/(x^2+1))
            apply 1 on (LIM {r->oo}. CINT z:com(C(t,r),L(t,r)). exp(i*z)/(z^2+1))
            apply 2 on (LIM {r->oo}. CINT z:C(t,r). exp(i*z)/(z^2+1))
            simplify
        done
    """
    actions = [s for s in actions.split('\n') if s.strip()]
    for act in actions:
        a = parser.parse_action(act)
        st = st.process_action(a)
    assert st.is_finished()

def test_contour_rational_cos_a2():
    """测试：∫_{-∞}^{+∞} cos(2x)/(x²+4) dx = π/(2e²)
    构造：∫ e^(2ix)/(x²+4) dx，取实部
    极点：z = 2i (上半平面), z = -2i (下半平面)
    留数：Res[e^(2iz)/(z²+4), z=2i] = e^(2i·2i)/(4i) = e^(-4)/(4i)
    结果：2πi × e^(-4)/(4i) = πe^(-4) = π/(2e²)
    """
    ctx = context.Context()
    ctx.load_book("base")
    st = state.InitialState(ctx)

    actions = """
        prove (INT x:[-oo,oo]. exp(2*i*x)/(x^2+4)) = pi/(2*exp(4))
        define C(t,r) = (r*exp(i*pi*(1-t)))_(t:[0,1])
        define L(t,r) = (r*(1-2*t))_(t:[0,1])
        subgoal 1: (LIM {r->oo}. CINT z:com(C(t,r),L(t,r)). exp(2*i*z)/(z^2+4)) = -pi/(2*exp(4))
        lhs:
            apply residue theorem
            rewrite to -pi/(2*exp(4))
        done
        subgoal 2: (LIM {r->oo}. CINT z:C(t,r). exp(2*i*z)/(z^2+4)) = 0
        lhs:
            apply cintegral identity
            substitute z for r * exp(i * pi * (1 - t))
            simplify
        done
        subgoal 3: (LIM {r->oo}. ((CINT z:com(C(t,r),L(t,r)). exp(2*i*z)/(z^2+4)) - (CINT z:C(t,r). exp(2*i*z)/(z^2+4)))) = -(INT x:[-oo,oo]. exp(2*i*x)/(x^2+4))
        lhs:
            simplify
            rewrite (CINT z:com(C(t,r),L(t,r)). exp(2*i*z)/(z^2+4)) to (CINT z:C(t,r). exp(2*i*z)/(z^2+4)) + (CINT z:L(t,r). exp(2*i*z)/(z^2+4))
            simplify
            rewrite (CINT z:L(t,r). exp(2*i*z)/(z^2+4)) to (-(INT x:[-r,r]. exp(2*i*x)/(x^2+4)))
            rewrite to -(INT x:[-oo,oo]. exp(2*i*x)/(x^2+4))
        done
        lhs:
            apply 3 on (INT x:[-oo,oo]. exp(2*i*x)/(x^2+4))
            apply 1 on (LIM {r->oo}. CINT z:com(C(t,r),L(t,r)). exp(2*i*z)/(z^2+4))
            apply 2 on (LIM {r->oo}. CINT z:C(t,r). exp(2*i*z)/(z^2+4))
            simplify
        done
    """
    actions = [s for s in actions.split('\n') if s.strip()]
    for act in actions:
        a = parser.parse_action(act)
        st = st.process_action(a)
    assert st.is_finished()

# ============ C类：高阶极点 ============

def test_contour_second_order_pole():
    """测试：∫_{-∞}^{+∞} dx/(x²+1)² = π/4
    极点：z = i (二阶), z = -i (二阶)
    二阶极点留数公式：Res[f,z₀] = lim_{z→z₀} d/dz[(z-z₀)²f(z)]
    """
    ctx = context.Context()
    ctx.load_book("base")
    st = state.InitialState(ctx)

    actions = """
        prove (INT x:[-oo,oo]. 1/(x^2+1)^2) = pi/2
        define C(t,r) = (r*exp(i*pi*(1-t)))_(t:[0,1])
        define L(t,r) = (r*(1-2*t))_(t:[0,1])
        subgoal 1: (LIM {r->oo}. CINT z:com(C(t,r),L(t,r)). 1/(z^2+1)^2) = -pi/2
        lhs:
            apply residue theorem
        done
        subgoal 2: (LIM {r->oo}. CINT z:C(t,r). 1/(z^2+1)^2) = 0
        lhs:
            apply cintegral identity
            substitute z for r * exp(i * pi * (1 - t))
            simplify
        done
        subgoal 3:(LIM {r->oo}. CINT z:com(C(t,r),L(t,r)). 1/(z^2+1)^2) - (LIM {r->oo}. CINT z:C(t,r). 1/(z^2+1)^2) = -(INT x:[-oo,oo]. 1/(x^2+1)^2)
        lhs:
            rewrite CINT z:com(C(t,r),L(t,r)). 1/(z^2+1)^2 to (CINT z:C(t,r). 1/(z^2+1)^2) + (CINT z:L(t,r). 1/(z^2+1)^2)
            simplify
            rewrite CINT z:L(t,r). 1 / (z ^ 2 + 1) ^ 2 to -(INT x:[-r,r]. 1/(x^2+1)^2)
            rewrite to -(INT x:[-oo,oo]. 1/(x^2+1)^2)
        done
        lhs:
            apply 3 on (INT x:[-oo,oo]. 1/(x^2+1)^2)
            apply 1 on (LIM {r -> oo}. CINT z:com(C(t,r),L(t,r)). 1 / (z ^ 2 + 1) ^ 2)
            apply 2 on (LIM {r->oo}. CINT z:C(t,r). 1/(z^2+1)^2)
            simplify
        done
    """
    actions = [s for s in actions.split('\n') if s.strip()]
    for act in actions:
        a = parser.parse_action(act)
        st = st.process_action(a)
    assert st.is_finished()


def test_contour_mixed_order_poles():
    """测试：∫_{-∞}^{+∞} dx/((x²+1)²(x²+4))
    同时存在二阶极点 (z=i) 和一阶极点 (z=2i)
    """
    ctx = context.Context()
    ctx.load_book("base")
    st = state.InitialState(ctx)

    actions = """
        prove (INT x:[-oo,oo]. 1/((x^2+1)^2*(x^2+4))) = pi/3
        define C(t,r) = (r*exp(i*pi*(1-t)))_(t:[0,1])
        define L(t,r) = (r*(1-2*t))_(t:[0,1])
        lhs:
            simplify
        done
    """
    actions = [s for s in actions.split('\n') if s.strip()]
    for act in actions:
        a = parser.parse_action(act)
        st = st.process_action(a)
    assert st.is_finished()


def test_contour_third_order_pole():
    """测试：∫_{-∞}^{+∞} dx/(x²+1)³ = 3π/8
    极点：z = i (三阶), z = -i (三阶)
    三阶极点留数公式需要二阶导数
    """
    ctx = context.Context()
    ctx.load_book("base")
    st = state.InitialState(ctx)

    actions = """
        prove (INT x:[-oo,oo]. 1/(x^2+1)^3) = 3*pi/8
        define C(t,r) = (r*exp(i*pi*(1-t)))_(t:[0,1])
        define L(t,r) = (r*(1-2*t))_(t:[0,1])
        lhs:
            simplify
        done
    """
    actions = [s for s in actions.split('\n') if s.strip()]
    for act in actions:
        a = parser.parse_action(act)
        st = st.process_action(a)
    assert st.is_finished()


# ============ D类：需要下半平面围道 ============

def test_contour_lower_half_exp_neg():
    """测试：∫_{-∞}^{+∞} e^(-ix)/(x²+1) dx = -πi/e
    由于 e^(-iz) 在上半平面指数增长，需选取下半平面围道
    极点：z = -i (下半平面)
    留数：Res[e^(-iz)/(z²+1), z=-i] = e^(-i·(-i))/(-2i) = e^(-1)/(-2i) = -e^(-1)/(2i)
    方向：顺时针（负向），需加负号
    结果：-2πi × (-e^(-1)/(2i)) = π/e，但需考虑方向为 -πi/e
    """
    ctx = context.Context()
    ctx.load_book("base")
    st = state.InitialState(ctx)

    actions = """
        prove (INT x:[-oo,oo]. exp(-i*x)/(x^2+1)) = -pi*i/exp(1)
        define C(t,r) = (r*exp(-i*pi*t))_(t:[0,1])
        define L(t,r) = (r*(1-2*t))_(t:[0,1])
        lhs:
            simplify
        done
    """
    actions = [s for s in actions.split('\n') if s.strip()]
    for act in actions:
        a = parser.parse_action(act)
        st = st.process_action(a)
    assert st.is_finished()


def test_contour_lower_half_two_poles():
    """测试：∫_{-∞}^{+∞} e^(-2ix)/((x²+1)(x²+4)) dx
    下半平面极点：z = -i, z = -2i
    """
    ctx = context.Context()
    ctx.load_book("base")
    st = state.InitialState(ctx)

    actions = """
        prove (INT x:[-oo,oo]. exp(-2*i*x)/((x^2+1)*(x^2+4))) = -pi*i/3*(exp(2)-exp(1))
        define C(t,r) = (r*exp(-i*pi*t))_(t:[0,1])
        define L(t,r) = (r*(1-2*t))_(t:[0,1])
        lhs:
            simplify
        done
    """
    actions = [s for s in actions.split('\n') if s.strip()]
    for act in actions:
        a = parser.parse_action(act)
        st = st.process_action(a)
    assert st.is_finished()


def test_contour_lower_half_opposite():
    """测试：∫_{-∞}^{+∞} e^(-ix)/(x²+4) dx = -πi/(2e)
    下半平面极点：z = -2i
    """
    ctx = context.Context()
    ctx.load_book("base")
    st = state.InitialState(ctx)

    actions = """
        prove (INT x:[-oo,oo]. exp(-i*x)/(x^2+4)) = -pi*i/(2*exp(2))
        define C(t,r) = (r*exp(-i*pi*t))_(t:[0,1])
        define L(t,r) = (r*(1-2*t))_(t:[0,1])
        lhs:
            simplify
        done
    """
    actions = [s for s in actions.split('\n') if s.strip()]
    for act in actions:
        a = parser.parse_action(act)
        st = st.process_action(a)
    assert st.is_finished()


# ============ E类：含对数或根式的多值函数 ============

def test_contour_log_rational():
    """测试：∫₀^∞ log(x)/(x²+1) dx = 0
    使用钥匙孔围道，支点 z=0
    结果为 0（对称性）
    """
    ctx = context.Context()
    ctx.load_book("base")
    st = state.InitialState(ctx)

    actions = """
        prove (INT x:[0,oo]. log(x)/(x^2+1)) = 0
        define C(t,r) = (r*exp(i*pi*(1-t)))_(t:[0,1])
        define L(t,r) = (r*(1-2*t))_(t:[0,1])
        lhs:
            simplify
        done
    """
    actions = [s for s in actions.split('\n') if s.strip()]
    for act in actions:
        a = parser.parse_action(act)
        st = st.process_action(a)
    assert st.is_finished()


def test_contour_sqrt_rational():
    """测试：∫₀^∞ √x/(x²+1) dx = π/2
    使用钥匙孔围道，√x = e^(1/2 log z)
    """
    ctx = context.Context()
    ctx.load_book("base")
    st = state.InitialState(ctx)

    actions = """
        prove (INT x:[0,oo]. sqrt(x)/(x^2+1)) = pi/2
        define C(t,r) = (r*exp(i*pi*(1-t)))_(t:[0,1])
        define L(t,r) = (r*(1-2*t))_(t:[0,1])
        lhs:
            simplify
        done
    """
    actions = [s for s in actions.split('\n') if s.strip()]
    for act in actions:
        a = parser.parse_action(act)
        st = st.process_action(a)
    assert st.is_finished()


def test_contour_fractional_power():
    """测试：∫₀^∞ x^(1/3)/(x²+1) dx
    非整数指数 α = 1/3，使用钥匙孔围道
    """
    ctx = context.Context()
    ctx.load_book("base")
    st = state.InitialState(ctx)

    actions = """
        prove (INT x:[0,oo]. x^(1/3)/(x^2+1)) = pi/2*sec(pi/6)
        define C(t,r) = (r*exp(i*pi*(1-t)))_(t:[0,1])
        define L(t,r) = (r*(1-2*t))_(t:[0,1])
        lhs:
            simplify
        done
    """
    actions = [s for s in actions.split('\n') if s.strip()]
    for act in actions:
        a = parser.parse_action(act)
        st = st.process_action(a)
    assert st.is_finished()


# ============ F类：边界与验证性用例 ============

def test_contour_no_poles():
    """测试：∫₀^R dx/(x²+4) 在有限区间上，极点在虚轴上
    验证系统对有限区间积分的处理
    """
    ctx = context.Context()
    ctx.load_book("base")
    st = state.InitialState(ctx)

    actions = """
        prove (INT x:[0,2]. 1/(x^2+4)) = pi/4
        lhs:
            simplify
        done
    """
    actions = [s for s in actions.split('\n') if s.strip()]
    for act in actions:
        a = parser.parse_action(act)
        st = st.process_action(a)
    assert st.is_finished()


def test_contour_real_axis_pole():
    """测试：主值积分 ∫_{-1}^{1} dx/x³
    极点在实轴上，需要处理主值
    """
    ctx = context.Context()
    ctx.load_book("base")
    st = state.InitialState(ctx)

    actions = """
        prove (INT x:[-1,1]. 1/x^3) = 0
        lhs:
            simplify
        done
    """
    actions = [s for s in actions.split('\n') if s.strip()]
    for act in actions:
        a = parser.parse_action(act)
        st = st.process_action(a)
    assert st.is_finished()


def test_contour_parameter_result():
    """测试：∫_{-∞}^{+∞} dx/(x²+a²) = π/a (a > 0)
    含参数 a 的积分，验证系统对参数的处理
    """
    ctx = context.Context()
    ctx.load_book("base")
    st = state.InitialState(ctx)

    actions = """
        let a = ?a
        prove (INT x:[-oo,oo]. 1/(x^2+a^2)) = pi/a
        define C(t,r) = (r*exp(i*pi*(1-t)))_(t:[0,1])
        define L(t,r) = (r*(1-2*t))_(t:[0,1])
        lhs:
            simplify
        done
    """
    actions = [s for s in actions.split('\n') if s.strip()]
    for act in actions:
        a = parser.parse_action(act)
        st = st.process_action(a)
    assert st.is_finished()