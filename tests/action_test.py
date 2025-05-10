"""Unit test for integrals using internal language."""

import unittest
import time
import sys
import cProfile
import pstats

from integral import compstate
from integral import state
from integral import parser
from integral import expr
from integral import context


class ActionTest(unittest.TestCase):
    def check_actions(self, content: str, *, print_lines=False, print_state=False,
                      write_stats=False, filename=""):
        actions = content.split('\n')
        ctx = context.Context()
        st = state.InitialState(ctx)
        start_time = None
        cur_goal = None
        for i, act in enumerate(actions, 1):
            if print_lines:
                print(act)
            if not act.strip():
                # empty line
                continue
            if act.lstrip().startswith('#') or act.lstrip().startswith('//'):
                # title or comment
                continue
            a = parser.parse_action(act)
            if isinstance(a, state.ImportsAction):
                for thy_name in a.theories:
                    ctx.load_book(thy_name)
            if isinstance(a, (state.ProveAction, state.CalculateAction)):
                cur_goal = a
                if write_stats:
                    start_time = time.time()
                    with open("stats.txt", "a", encoding='utf-8') as stats_file:
                        stats_file.write(f"{filename} {i} {cur_goal}\n")
            try:
                st = st.process_action(a)
                if isinstance(st, state.InitialState):
                    if isinstance(cur_goal, state.ProveAction):
                        if cur_goal.expr.is_equals() and expr.is_indefinite_integral(cur_goal.expr.lhs):
                            ctx.add_indefinite_integral(cur_goal.expr, cur_goal.conditions, cur_goal.attrs)
                        elif cur_goal.expr.is_equals() and expr.is_integral(cur_goal.expr.lhs):
                            ctx.add_definite_integral(cur_goal.expr, cur_goal.conditions, cur_goal.attrs)
                        else:
                            ctx.add_other_identities(cur_goal.expr, cur_goal.conditions, cur_goal.attrs)
                    if cur_goal and write_stats:
                        elapsed_time = time.time() - start_time
                        with open("stats.txt", "a", encoding='utf-8') as stats_file:
                            stats_file.write(f"{elapsed_time:.2f} seconds\n")
                    cur_goal = None
            except Exception as e:
                print(cur_goal)
                print(st)
                raise e
        if print_state:
            print(st)
        if not print_state and not isinstance(st, state.InitialState):
            raise AssertionError("Does not end in initial state (add print_state=True to debug)")
                
    def check_file(self, filename: str, *, print_lines=False, print_state=False,
                   write_stats=False):
        with open(f'theories/{filename}.thy', 'r', encoding='utf-8') as f:
            content = f.read()
        self.check_actions(content, print_lines=print_lines, print_state=print_state,
                           write_stats=write_stats, filename=filename)

    def testCalculationFinished(self):
        ctx = context.Context()
        ctx.load_book("base")
        st = state.InitialState(ctx)

        actions = """
            calculate INT x. (3 - 2*x)^3
                substitute u for 3 - 2*x
                apply integral identity
                simplify
        """
        actions = [s for s in actions.split('\n') if s.strip()]
        for act in actions:
            a = parser.parse_action(act)
            st = st.process_action(a)
        self.assertFalse(st.is_finished())

    def testCalculationFinished2(self):
        ctx = context.Context()
        ctx.load_book("base")
        st = state.InitialState(ctx)

        actions = """
            calculate INT x. (3 - 2*x)^3
                substitute u for 3 - 2*x
                apply integral identity
                simplify
                replace substitution
                simplify
        """
        actions = [s for s in actions.split('\n') if s.strip()]
        for act in actions:
            a = parser.parse_action(act)
            st = st.process_action(a)
        self.assertTrue(st.is_finished())

    def testStandard(self):
        self.check_file("standard")

    def testStandard2(self):
        self.check_file("standard2")

    def testStandard3(self):
        self.check_file("standard3")

    def testStandard4(self):
        self.check_file("standard4")

    def testMIT2019(self):
        self.check_file("mit2019")

    def testLHopital(self):
        self.check_file("lhopital")

    def testTongji(self):
        self.check_file("tongji05")

    def testTongjiIndefSubstitution(self):
        self.check_file("tongji0402")

    def testTongjiIndefByParts(self):
        self.check_file("tongji0403")

    def testTongjiIndefRational(self):
        self.check_file("tongji0404")

    def testUSubstitution(self):
        self.check_file("ucdavisUSubst")

    def testUCDavisPartialFraction(self):
        self.check_file("ucdavisPartial")

    def testIntegrateByParts(self):
        self.check_file("ucdavisByParts")

    def testExponential(self):
        self.check_file("ucdavisExponential")

    def testTrigonometric(self):
        self.check_file("ucdavisTrigonometric")

    def testLogAndArcTangent(self):
        self.check_file("ucdavisLogArctan")

    def testPowerSubstitution(self):
        self.check_file("ucdavisPowerSubst")

    def testTrigSubstitution(self):
        self.check_file("ucdavisTrigSubst")

    def testWallis(self):
        self.check_file("irresistable")

    def testInteresting1(self):
        self.check_file("interesting1")

    def testInteresting2(self):
        self.check_file("interesting2")

    def testInteresting3(self):
        self.check_file("interesting3")

    def testInteresting4(self):
        self.check_file("interesting4")

    def testInteresting5(self):
        self.check_file("interesting5")

    def testInteresting6(self):
        self.check_file("interesting6")

    def testLeibniz03(self):
        # Inside interesting integrals, Section 3.1
        # TODO: Still cannot remove the condition I(t) > 0
        actions = """
            imports interesting3

            prove (INT x:[0,oo]. cos(t*x)*exp(-(x^2)/2)) = sqrt(pi/2)*exp(-(t^2)/2)
            let I(t) = INT x:[0,oo]. cos(t*x)*exp(-(x^2)/2)
            subgoal 1: I(0) = sqrt(pi/2)
            lhs:
                expand definition for I
                rewrite -(x^2/2) to -(x^2)/2
                apply integral identity
            rhs:
                simplify
            done

            subgoal 2: (D t. I(t)) = -t*I(t)
            lhs:
                expand definition for I (all)
                simplify
                integrate by parts with u = sin(t*x), v = -exp(-x^2/2)
                simplify
            rhs:
                expand definition for I (all)
                simplify
            done

            subgoal 3: (D t. log(I(t)) + t^2/2) = 0
            lhs:
                simplify
                apply 2 on D t. I(t)
                simplify
            done

            subgoal 4: 1/2 * t ^ 2 + log(I(t)) = SKOLEM_CONST(C) for I(t) > 0
            from 3:
                integrate both sides
                apply integral identity
            done

            subgoal 5: log(sqrt(pi / 2)) = SKOLEM_CONST(C) for I(t) > 0
            from 4:
                apply limit t -> 0 both sides
                simplify
                apply 1 on I(0)
            done

            subgoal 6: log(I(t)) = -t^2 / 2 + log(sqrt(pi / 2)) for I(t) > 0
            lhs:
                apply 4 on log(I(t))
                apply 5 on SKOLEM_CONST(C)
            done

            subgoal 7: I(t) = sqrt(pi / 2) * exp(-t^2/2) for I(t) > 0
            from 6:
                solve equation for I(t)
                rewrite exp(-(t ^ 2 / 2) - log(2) / 2 + log(pi) / 2) to exp(-(t^2/2)) / exp(log(2)/2) * exp(log(pi)/2)
                rewrite exp(log(2) / 2) to sqrt(2)
                rewrite exp(log(pi) / 2) to sqrt(pi)
            done
        """
        try:
            self.check_actions(actions)
        except compstate.CheckFinishedException as e:
            ()

    def testGaussianPowerExp(self):
        # Inside interesting integrals, Section 2.3
        actions = """
            imports interesting3

            prove (INT x:[0, oo]. x^(2*n) * exp(-x^2)) = factorial(2*n)/(4^n*factorial(n))*(1/2)*sqrt(pi) for n: int, n >= 0
            let I(n) = (INT x:[0, oo]. x^(2*n) * exp(-x^2))
            subgoal 1: (INT x:[0, oo]. (D x. x^(2*n-1)*exp(-x^2))) = 0 for n>=1
            lhs:
                simplify
            done
            subgoal 2: (D x. x^(2*n-1)*exp(-x^2)) = (2*n-1)*x^(2*n-2)*exp(-x^2)-2*x^(2*n)*exp(-x^2) for x > 0
            lhs:
                simplify
                rewrite x ^ (2 * n - 2) * exp(-(x ^ 2)) * (2 * n - 1) to (2*n-1)*x^(2*n-2)*exp(-x^2) 
            done
            subgoal 3: (INT x:[0, oo]. x^(2*n) * exp(-x^2)) = I(n)
            rhs:
                expand definition for I
            done
            subgoal 4: (INT x:[0, oo]. (D x. x^(2*n-1)*exp(-x^2))) = (2*n-1)*I(n-1) - 2 * I(n) for n>=1
            lhs:
                apply 2 on (D x. x^(2*n-1)*exp(-x^2))
                simplify
                apply 3 on (INT x:[0,oo]. x ^ (2 * n) * exp(-(x ^ 2)))
                rewrite (INT x:[0,oo]. x ^ (2 * n - 2) * exp(-(x ^ 2))) to (INT x:[0,oo]. x ^ (2 * (n - 1)) * exp(-(x ^ 2)))
                apply 3 on (INT x:[0,oo]. x ^ (2 * (n - 1)) * exp(-(x ^ 2)))
            done
            subgoal 5: I(n) = I(n-1)*(2*n-1)/2
            from 4:
                apply 1 on (INT x:[0, oo]. (D x. x^(2*n-1)*exp(-x^2)))
                solve equation for I(n)
                rewrite I(n - 1) * (2 * n - 1) / 2 to I(n - 1) * ((2*n)*(2 * n - 1)) / (2*(2*n))
                rewrite (2 * (2 * n)) to (4 * n)
            done
            subgoal 6: I(n) = I(0)*factorial(2*n)/(4^n*factorial(n))
            induction on n
                base:
                lhs:
                done
                induct:
                lhs:
                    apply 5 on I(n+1)
                    simplify
                    apply induction hypothesis (all)
                rhs:
                    rewrite factorial(n+1) to factorial(n) * (n+1)
                    rewrite factorial(2*n+2) to factorial(2*n+1) * (2*n+2)
                    rewrite factorial(2*n+1) to factorial(2*n) * (2*n+1)
                    simplify
                    rewrite to I(0) * factorial(2 * n) / (4 ^ n * factorial(n)) * (2 * n + 1) / 2
                done
            done
            subgoal 7: I(0) = 1/2*sqrt(pi)
            lhs:
                expand definition for I
                substitute y for x*sqrt(2)
                rewrite -(y^2/2) to -(y^2)/2
                apply integral identity
            done

            lhs:
                fold definition for I
                apply 6 on I(n)
                apply 7 on I(0)
            done
            """
        # requires evaluation of probability integral
        self.check_actions(actions)

    # def testFlipside08(self):
    #     actions = """
    #         prove (INT x:[0,pi]. 1/(a+b*cos(x))) = pi/(sqrt(a^2-b^2)) for a > b, b >= 0
    #         define I(a, b) = (INT x:[0,pi]. log(a+b*cos(x)))
    #         subgoal 1: (D a. I(a, b)) = (INT x:[0,pi]. 1/(a+b*cos(x)))
    #         lhs:
    #             expand definition for I (all)
    #             simplify
    #         done
    #     """
    #     try:
    #         self.check_actions(actions, "interesting")
    #     except compstate.CheckFinishedException as e:
    #         ()

    def testEulerFormula1(self):
        # TODO ([log(x)]_x=1,oo) - 1/2 * ([log(x - i)]_x=1,oo) - 1/2 * ([log(x + i)]_x=1,oo) ->
        #                                        [(log(x)) - 1/2 * (log(x - i)) - 1/2 * (log(x + i))]_x=1,oo
        actions = """
            imports standard
            prove (INT x:[1,oo]. 1/(x*(x^2+1))) = log(2)/2 for x:real, x!=0
            lhs:
                rewrite 1/(x*(x^2+1)) to 1/x - 1/(2*(x-i)) - 1/(2*(x+i))
                apply integral identity
                simplify
                rewrite to (log(-i + 1) + log(i + 1)) / 2
                rewrite log(-i + 1) + log(i + 1) to log((-i + 1)*(i + 1))
                rewrite (-i + 1)*(i + 1) to (-i*i-i+i+1)
                simplify
            done
        """
        self.check_actions(actions)

    def testEulerFormula2(self):
        actions = """
            imports standard
            prove (INT x:[0,oo]. sin(b*x)*exp(-x*y)) = b/(y^2+b^2) for b: real, y > 0
            lhs:
                rewrite sin(b*x) to (exp(i*(b*x)) - exp(-i*(b*x))) / (2*i)
                expand polynomial
                simplify
                rewrite -(b * x * i) - x * y to (-y - b*i) * x
                rewrite b * x * i - x * y to (-y + b*i) * x
                apply integral identity
                rewrite x * (-(b*i) - y) to  -x * (b*i) - x * y
                rewrite x * (b * i - y) to x * b * i - x * y
                simplify
                rewrite to b / (y^2 + b^2)
            done
        """
        self.check_actions(actions)

    def testPostgraduateIndefinitePart1SectionA(self):
        self.check_file("postgradIndef1a")

    def testPostgraduateIndefinitePart1SectionB(self):
        self.check_file("postgradIndef1b")

    def testPostgraduateIndefinitePart2SectionA(self):
        self.check_file("postgradIndef2a")

    def testPostgraduateIndefinitePart2SectionB(self):
        self.check_file("postgradIndef2b")

    def testPostgraduateIndefinitePart3SectionA(self):
        self.check_file("postgradIndef3a")

    def testPostgraduateIndefinitePart4SectionA(self):
        self.check_file("postgradIndef4a")

    def testPostgraduateIndefinitePart4SectionB(self):
        self.check_file("postgradIndef4b")

    def testPostgraduateIndefinitePart5SectionA(self):
        self.check_file("postgradIndef5a")

    def testPostgraduateIndefinitePart5SectionB(self):
        self.check_file("postgradIndef5b")

    def testPostgraduateIndefinitePart6SectionA(self):
        self.check_file("postgradIndef6a")

    def testPostgraduateIndefinitePart6SectionB(self):
        self.check_file("postgradIndef6b")

    def testPostgraduateDefinitePart1SectionA(self):
        self.check_file("postgradDef1a")

    def testPostgraduateDefinitePart1SectionB(self):
        self.check_file("postgradDef1b")


if __name__ == "__main__":
    if "--profile" in sys.argv:
        sys.argv.remove("--profile")
        profiler = cProfile.Profile()
        profiler.enable()
        
        # Run tests
        unittest.TestProgram(exit=False)
        
        profiler.disable()
        print("\n\n--- Profiling Results ---")
        stats = pstats.Stats(profiler)
        stats.sort_stats(pstats.SortKey.CUMULATIVE)
        stats.print_stats(50)  # Show top 50 functions by cumulative time
    else:
        unittest.main()
