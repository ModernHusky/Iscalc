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
    def check_actions(self, base_file: str, actions: str,
                      *, print_lines=False, print_state=False, write_stats=True):
        ctx = context.Context()
        ctx.load_book(base_file)
        st = state.InitialState(ctx)
        actions = [s for s in actions.split('\n') if s.strip()]
        start_time = None
        cur_goal = None
        for act in actions:
            if print_lines:
                print(act)
            if act.lstrip().startswith('#') or act.lstrip().startswith('//'):
                # title or comment
                continue
            a = parser.parse_action(act)
            if isinstance(a, (state.ProveAction, state.CalculateAction)):
                cur_goal = a
                if write_stats:
                    start_time = time.time()
                    with open("stats.txt", "a", encoding='utf-8') as stats_file:
                        stats_file.write(f"{cur_goal}\n")
            try:
                st = st.process_action(a)
                if isinstance(st, state.InitialState):
                    if isinstance(cur_goal, state.ProveAction):
                        if cur_goal.expr.is_equals() and expr.is_indefinite_integral(cur_goal.expr.lhs):
                            ctx.add_indefinite_integral(cur_goal.expr, cur_goal.conditions)
                        elif cur_goal.expr.is_equals() and expr.is_integral(cur_goal.expr.lhs):
                            ctx.add_definite_integral(cur_goal.expr, cur_goal.conditions)
                        else:
                            ctx.add_other_identities(cur_goal.expr, cur_goal.attrs, cur_goal.conditions)
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
        with open('theories/standard.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("base", actions)

    def testStandard2(self):
        with open('theories/standard2.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("base", actions)

    def testStandard3(self):
        with open('theories/standard3.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("base", actions)

    def testStandard4(self):
        with open('theories/standard4.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("base", actions)

    def testMIT2019(self):
        actions = """
            calculate INT x:[0,pi / 100]. (sin(20 * x) + sin(19 * x)) / (cos(20 * x) + cos(19 * x))
                rewrite sin(20 * x) + sin(19 * x) to 2 * cos(1/2 * x) * sin(39/2 * x)
                rewrite cos(20 * x) + cos(19 * x) to 2 * cos(1/2 * x) * cos(39/2 * x)
                simplify
                substitute u for cos(39/2 * x)
                apply integral identity
                simplify
            done
        """        
        self.check_actions("standard", actions)

    def testLHopital(self):
        with open('theories/lhopital.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("standard", actions)

    def testTongji(self):
        with open('theories/tongji05.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("standard", actions)

    def testTongjiIndefSubstitution(self):
        with open('theories/tongji0402.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("standard", actions)

    def testTongjiIndefByParts(self):
        with open('theories/tongji0403.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("standard", actions)

    def testTongjiIndefRational(self):
        with open('theories/tongji0404.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("standard", actions)

    def testUSubstitution(self):
        with open('theories/ucdavisUSubst.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("standard", actions)

    def testActions2(self):
        actions = """
        prove (INT x. 1/sqrt(-(x^2)+a)) = arcsin(x/sqrt(a))+ SKOLEM_CONST(C) for a > 0, -x^2 + a > 0, x / sqrt(a) <= 1, x / sqrt(a) >= -1
        lhs:
            rewrite sqrt(-(x^2)+a) to sqrt(a - x^2)
            rewrite sqrt(a - x^2) to sqrt(a*(1 - x^2/a))
            rewrite sqrt(a*(1 - x^2/a)) to sqrt(a)*sqrt(1 - (x/sqrt(a))^2)
            rewrite 1/(sqrt(a)*sqrt(1 - (x/sqrt(a))^2)) to (1/sqrt(a))*(1/sqrt(1 - (x/sqrt(a))^2))
            substitute u for x/sqrt(a)
            simplify
            apply integral identity
            replace substitution
        done
        """
        self.check_actions("standard", actions)

    def testUCDavisPartialFraction(self):
        with open('theories/ucdavisPartial.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("standard", actions)

    def testIntegrateByParts(self):
        with open('theories/ucdavisByParts.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("standard", actions)

    def testExponential(self):
        with open('theories/ucdavisExponential.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("standard", actions)

    def testTrigonometric(self):
        with open('theories/ucdavisTrigonometric.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("standard", actions)

    def testLogAndArcTangent(self):
        with open('theories/ucdavisLogArctan.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("standard", actions)

    def testPowerSubstitution(self):
        with open('theories/ucdavisPowerSubst.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("standard", actions)

    def testTrigSubstitution(self):
        with open('theories/ucdavisTrigSubst.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("standard", actions)

    def testWallis(self):
        # Reference:
        # Irresistable Integrals, Section 2.3
        actions = """
            prove (INT x:[0,oo]. 1 / (x ^ 2 + b) ^ (m + 1)) = pi / 2 ^ (2 * m + 1) * binom(2 * m,m) * (1 / b ^ ((2 * m + 1) / 2)) for m: int, b: real, b > 0, m >= 0
            let I(m,b) = (INT x:[0,oo]. 1 / (x ^ 2 + b) ^ (m + 1)) for b > 0, m >= 0
            subgoal 1: (D b. I(m,b)) = -(m + 1) * I(m + 1,b)
            lhs:
                expand definition for I (all)
                exchange derivative and integral
                simplify
            rhs:
                expand definition for I (all)
                simplify
            done

            subgoal 2: I(m,b) = pi / 2 ^ (2 * m + 1) * binom(2 * m,m) * (1 / b ^ ((2 * m + 1) / 2)) for m: int
            induction on m
                base:
                    lhs:
                        expand definition for I
                        substitute sqrt(b) * u for x
                        simplify
                        rewrite 1 / (b * u ^ 2 + b) to 1 / b * (1 / (1 ^ 2 + u ^ 2))
                        apply integral identity
                        simplify
                done
                induct:
                    lhs:
                        apply 1 on I(m + 1,b)
                        apply induction hypothesis (all)
                        simplify
                        rewrite -((2 * m + 1) / 2) - 1 to -m - 3/2
                        rewrite to b ^ (-m - 3/2) * 2 ^ -(2 * m) * pi * (2 * m + 1) / (4 * m + 4) * binom(2 * m,m)
                    rhs:
                        rewrite binom(2 * m + 2,m + 1) to 2 * binom(2 * m,m) * ((2 * m + 1) / (m + 1))
                        rewrite -((2 * m + 3) / 2) to -m - 3/2
                        simplify
                done
            done

            lhs:
                fold definition for I (all)
                apply 2 on I(m,b)
            done
        """
        self.check_actions("standard", actions)

    def testInteresting1(self):
        with open("theories/interesting1.thy", 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("standard", actions)

    def testInteresting2(self):
        with open("theories/interesting2.thy", 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("interesting1", actions)

    def testInteresting3(self):
        with open("theories/interesting3.thy", 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("interesting2", actions)

    def testInteresting4(self):
        with open("theories/interesting4.thy", 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("interesting3", actions)

    def testInteresting5(self):
        with open("theories/interesting5.thy", 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("interesting4", actions)

    def testInteresting6(self):
        with open("theories/interesting6.thy", 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("interesting5", actions)

    def testLeibniz03(self):
        # Inside interesting integrals, Section 3.1
        # TODO: Still cannot remove the condition I(t) > 0
        actions = """
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
            self.check_actions("interesting3", actions)
        except compstate.CheckFinishedException as e:
            ()

    def testGaussianPowerExp(self):
        # Inside interesting integrals, Section 2.3
        actions = """
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
                rhs:
                    simplify
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
        self.check_actions("interesting3", actions)

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
    #         self.check_actions("interesting", "flipside08", actions)
    #     except compstate.CheckFinishedException as e:
    #         ()

    def testEulerFormula1(self):
        # TODO ([log(x)]_x=1,oo) - 1/2 * ([log(x - i)]_x=1,oo) - 1/2 * ([log(x + i)]_x=1,oo) ->
        #                                        [(log(x)) - 1/2 * (log(x - i)) - 1/2 * (log(x + i))]_x=1,oo
        actions = """
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
        self.check_actions("standard", actions)

    def testEulerFormula2(self):
        actions = """
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
        self.check_actions("standard", actions)

    def testPostgraduateIndefinitePart1SectionA(self):
        with open('theories/postgradIndef1a.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("standard", actions)

    def testPostgraduateIndefinitePart1SectionB(self):
        with open('theories/postgradIndef1b.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("standard", actions)

    def testPostgraduateIndefinitePart2SectionA(self):
        with open('theories/postgradIndef2a.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("standard", actions)

    def testPostgraduateIndefinitePart2SectionB(self):
        with open('theories/postgradIndef2b.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("standard", actions)

    def testPostgraduateIndefinitePart3SectionA(self):
        with open('theories/postgradIndef3a.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("standard", actions)

    def testPostgraduateIndefinitePart4SectionA(self):
        with open('theories/postgradIndef4a.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("standard", actions)

    def testPostgraduateIndefinitePart4SectionB(self):
        with open('theories/postgradIndef4b.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("standard", actions)

    def testPostgraduateIndefinitePart5SectionA(self):
        with open('theories/postgradIndef5a.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("standard", actions)

    def testPostgraduateIndefinitePart5SectionB(self):
        with open('theories/postgradIndef5b.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("standard", actions)

    def testPostgraduateIndefinitePart6SectionA(self):
        with open('theories/postgradIndef6a.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("standard", actions)

    def testPostgraduateIndefinitePart6SectionB(self):
        with open('theories/postgradIndef6b.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("standard", actions)

    def testPostgraduateDefinitePart1SectionA(self):
        with open('theories/postgradDef1a.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("standard", actions)

    def testPostgraduateDefinitePart1SectionB(self):
        with open('theories/postgradDef1b.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("standard", actions)


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
