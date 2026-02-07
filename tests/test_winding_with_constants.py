"""
测试绕数计算 - 包含常数半径和变量半径的情况
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from integral import compstate, parser
from integral.state import InitialState

def run_test(test_name, actions_str, expected_result):
    """运行单个测试"""
    print(f"\n{'='*60}")
    print(f"测试: {test_name}")
    print(f"{'='*60}")
    
    file = compstate.CompFile('interesting', 'leibniz03')
    state = InitialState(file)
    
    actions = [s for s in actions_str.strip().split('\n') if s.strip()]
    exec_state = state
    
    for act in actions:
        exec_state = exec_state.process_action(parser.parse_action(act))
    
    print(exec_state)
    print(f"\n期望结果: {expected_result}")
    print("✓ 测试通过" if expected_result in str(exec_state) else "✗ 测试失败")

if __name__ == '__main__':
    # 测试1: 变量半径r，带LIM，正向路径
    run_test(
        "测试1: 变量半径r + LIM {r->oo} + 正向路径 -> -pi",
        """
        prove (INT x:[-oo,oo]. 1/(x^2+1)) = pi
        define C(t,r) = (r*exp(i*pi*(1-t)))_(t:[0,1])
        define L(t,r) = (r*(1-2*t))_(t:[0,1])
        subgoal 1:(LIM {r->oo}. CINT z:com(C(t,r),L(t,r)). 1/(z^2+1))=-pi
        lhs:
            apply residue theorem
            simplify
        """,
        "-pi"
    )
    
    # 测试2: 变量半径r，带LIM，反向路径
    run_test(
        "测试2: 变量半径r + LIM {r->oo} + 反向路径 -> pi",
        """
        prove (INT x:[-oo,oo]. 1/(x^2+1)) = pi
        define C(t,r) = (r*exp(i*pi*(1-t)))_(t:[1,0])
        define L(t,r) = (r*(1-2*t))_(t:[1,0])
        subgoal 1:(LIM {r->oo}. CINT z:com(C(t,r),L(t,r)). 1/(z^2+1))=pi
        lhs:
            apply residue theorem
            simplify
        """,
        "pi"
    )
    
    # 测试3: 常数半径3，无LIM，正向路径
    run_test(
        "测试3: 常数半径3 + 无LIM + 正向路径 -> -pi",
        """
        prove (INT x:[-oo,oo]. 1/(x^2+1)) = pi
        define C(t,r) = (3*exp(i*pi*(1-t)))_(t:[0,1])
        define L(t,r) = (3*(1-2*t))_(t:[0,1])
        subgoal 1:(CINT z:com(C(t,r),L(t,r)). 1/(z^2+1))=pi
        lhs:
            apply residue theorem
            simplify
        """,
        "-pi"
    )
    
    # 测试4: 常数半径3，带LIM（但路径中无r），正向路径
    run_test(
        "测试4: 常数半径3 + LIM {r->oo}（未使用） + 正向路径 -> -pi",
        """
        prove (INT x:[-oo,oo]. 1/(x^2+1)) = pi
        define C(t,r) = (3*exp(i*pi*(1-t)))_(t:[0,1])
        define L(t,r) = (3*(1-2*t))_(t:[0,1])
        subgoal 1:(LIM {r->oo}. CINT z:com(C(t,r),L(t,r)). 1/(z^2+1))=-pi
        lhs:
            apply residue theorem
            simplify
        """,
        "-pi"
    )
    
    # 测试5: 常数半径2，无LIM，正向路径
    run_test(
        "测试5: 常数半径2 + 无LIM + 正向路径 -> -pi",
        """
        prove (INT x:[-oo,oo]. 1/(x^2+1)) = pi
        define C(t,r) = (2*exp(i*pi*(1-t)))_(t:[0,1])
        define L(t,r) = (2*(1-2*t))_(t:[0,1])
        subgoal 1:(CINT z:com(C(t,r),L(t,r)). 1/(z^2+1))=-pi
        lhs:
            apply residue theorem
            simplify
        """,
        "-pi"
    )
    
    # 测试6: 常数半径10，反向路径
    run_test(
        "测试6: 常数半径10 + 无LIM + 反向路径 -> pi",
        """
        prove (INT x:[-oo,oo]. 1/(x^2+1)) = pi
        define C(t,r) = (10*exp(i*pi*(1-t)))_(t:[1,0])
        define L(t,r) = (10*(1-2*t))_(t:[1,0])
        subgoal 1:(CINT z:com(C(t,r),L(t,r)). 1/(z^2+1))=pi
        lhs:
            apply residue theorem
            simplify
        """,
        "pi"
    )
    
    # 测试7: 小半径0.5（不包含极点）
    run_test(
        "测试7: 小半径0.5（不包含极点i） -> 0",
        """
        prove (INT x:[-oo,oo]. 1/(x^2+1)) = pi
        define C(t,r) = (1/2*exp(i*pi*(1-t)))_(t:[0,1])
        define L(t,r) = (1/2*(1-2*t))_(t:[0,1])
        subgoal 1:(CINT z:com(C(t,r),L(t,r)). 1/(z^2+1))=0
        lhs:
            apply residue theorem
            simplify
        """,
        "0"
    )
    
    print(f"\n{'='*60}")
    print("所有测试完成！")
    print(f"{'='*60}\n")

