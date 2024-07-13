"""Unit test for integrals using internal language."""

import unittest
import lark

from integral import compstate
from integral import action
from integral import parser


class ActionTest(unittest.TestCase):
    def check_actions(self, base_file, current_file, actions, print_state=False):
        file = compstate.CompFile(base_file, current_file)
        state = action.InitialState(file)
        actions = [s for s in actions.split('\n') if s.strip()]
        for act in actions:
            a = parser.parse_action(act)
            state = state.process_action(a)
        while not isinstance(state.past, action.InitialState):
            state = state.past
        if print_state:
            print(state)
        self.assertTrue(state.is_finished())

    def testTongji(self):
        actions = """
            calculate INT x:[2,3]. 2 * x + x ^ 2
            apply integral identity
            simplify
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT x:[0,1]. (3 * x + 1) ^ (-2)
            substitute u for 3 * x + 1
            apply integral identity
            simplify
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT x:[0,1]. exp(6 * x)
            apply integral identity
            simplify
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT x:[-1,2]. x * exp(x)
            integrate by parts with u = x, v = exp(x)
            apply integral identity
            simplify
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT x:[0,pi/4]. sin(x)
            apply integral identity
            simplify
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT x:[0,1]. 3*x^2 - x + 1
            apply integral identity
            simplify
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT x:[1,2]. x^2 + 1/x^4
            apply integral identity
            simplify
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT x:[pi/3, pi]. sin(2*x + pi/3)
            substitute u for 2*x + pi/3
            apply integral identity
            simplify
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT x:[4, 9]. x ^ (1 / 3) * (x ^ (1 / 2) + 1)
            expand polynomial
            apply integral identity
            simplify
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT x:[-1, 0]. (3 * x ^ 4 + 3 * x ^ 2 + 1) / (x ^ 2 + 1)
            partial fraction decomposition
            apply integral identity
            simplify
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT x:[4, exp(1) + 3]. (x ^ 3 - 12 * x ^ 2 - 42) / (x - 3)
            partial fraction decomposition
            apply integral identity
            simplify
            substitute u for x - 3
            apply integral identity
            expand polynomial
            simplify
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT x:[0, pi / 2]. sin(x) * cos(x) ^ 3
            substitute u for cos(x)
            apply integral identity
            simplify
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT x:[0, pi]. 1 - sin(x) ^ 3
            simplify
            rewrite sin(x) ^ 3 to sin(x) * sin(x) ^ 2
            rewrite sin(x) ^ 2 to 1 - cos(x) ^ 2 using identity
            substitute u for cos(x)
            apply integral identity
            simplify
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT x:[pi/6, pi/2]. cos(x) ^ 2
            apply integral identity
            simplify
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT x:[0, 1]. (1 - x^2) ^ (1/2)
            inverse substitute sin(u) for x creating u
            rewrite 1 - sin(u) ^ 2 to cos(u) ^ 2 using identity
            simplify
            apply integral identity
            simplify
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT x:[0, sqrt(2)]. sqrt(2 - x^2)
            inverse substitute sqrt(2) * sin(u) for x creating u
            simplify
            rewrite sin(u) ^ 2 to 1 - cos(u) ^ 2 using identity
            simplify
            apply integral identity
            simplify
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT y:[-sqrt(2), sqrt(2)]. sqrt(8 - 2*y^2)
            inverse substitute 2 * sin(u) for y creating u
            simplify
            rewrite sin(u) ^ 2 to 1 - cos(u) ^ 2 using identity
            simplify
            apply integral identity
            expand polynomial
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT x:[1/sqrt(2), 1]. sqrt(1 - x^2) / x ^ 2
            inverse substitute sin(u) for x creating u
            simplify
            rewrite sin(u) ^ 2 to 1 - cos(u) ^ 2 using identity
            simplify
            rewrite cos(u) ^ 2 to 1 - sin(u) ^ 2 using identity
            expand polynomial
            simplify
            rewrite 1 / sin(u) ^ 2 to csc(u) ^ 2 using identity
            apply integral identity
            simplify
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT x:[-1, 1]. x / sqrt(5 - 4 * x)
            substitute u for 5 - 4 * x
            expand polynomial
            simplify
            apply integral identity
            simplify
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT x:[1,4]. 1 / (1 + sqrt(x))
            substitute u for sqrt(x)
            substitute v for u + 1
            expand polynomial
            simplify
            apply integral identity
            simplify
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT x:[3/4, 1]. 1 / (sqrt(1-x) - 1)
            substitute u for sqrt(1 - x)
            substitute v for u - 1
            expand polynomial
            simplify
            apply integral identity
            simplify
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT t:[0, 1]. t * exp(-t ^ 2 / 2)
            substitute u for t ^ 2 / 2
            apply integral identity
            simplify
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT x:[1, exp(2)]. 1 / (x * sqrt(1 + log(x)))
            substitute u for 1 + log(x)
            apply integral identity
            simplify
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT x:[-2, 0]. (x + 2) / (x^2 + 2*x + 2)
            rewrite x^2 + 2*x + 2 to (x + 1) ^ 2 + 1
            substitute u for x + 1
            expand polynomial
            simplify
            split region at 0
            substitute v for u ^ 2 + 1
            apply integral identity
            substitute v for u ^ 2 + 1
            apply integral identity
            simplify
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT x:[-pi/2, pi/2]. cos(x) ^ 4
            rewrite cos(x) ^ 4 to (cos(x) ^ 2) ^ 2
            rewrite cos(x) ^ 2 to (1 + cos(2*x)) / 2 using identity
            expand polynomial
            substitute u for 2 * x
            simplify
            apply integral identity
            simplify
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT x:[-pi/2, pi/2]. sqrt(cos(x) - cos(x)^3)
            rewrite cos(x) - cos(x)^3 to cos(x) * (1 - cos(x)^2)
            rewrite 1 - cos(x)^2 to sin(x)^2 using identity
            simplify
            split region at 0
            simplify
            substitute u for cos(x)
            apply integral identity
            substitute u for cos(x)
            apply integral identity
            simplify
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT x:[0, pi]. sqrt(1 + cos(2*x))
            rewrite cos(2 * x) to 2 * cos(x) ^ 2 - 1 using identity
            simplify
            split region at pi / 2
            simplify
            apply integral identity
            simplify
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT x:[0, 1]. x * exp(-x)
            integrate by parts with u = x, v = -exp(-x)
            simplify
            apply integral identity
            simplify
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT x:[1, exp(1)]. x * log(x)
            integrate by parts with u = log(x) / 2, v = x ^ 2
            apply integral identity
            simplify
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT x:[pi/4, pi/3]. x / sin(x)^2
            rewrite sin(x) ^ 2 to csc(x) ^ -2 using identity
            integrate by parts with u = x, v = -cot(x)
            simplify
            rewrite cot(x) to cos(x) / sin(x) using identity
            substitute u for sin(x)
            apply integral identity
            simplify
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT x:[1, 4]. log(x) / sqrt(x)
            integrate by parts with u = 2 * log(x), v = sqrt(x)
            apply integral identity
            simplify
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT x:[0, 1]. x * atan(x)
            integrate by parts with u = atan(x) / 2, v = x ^ 2
            simplify
            rewrite x^2 / (2 * x^2 + 2) to (1 - 1 / (x^2 + 1)) / 2
            apply integral identity
            simplify
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT x:[0, pi/2]. exp(2*x) * cos(x)
            integrate by parts with u = exp(2*x), v = sin(x)
            simplify
            integrate by parts with u = exp(2*x), v = -cos(x)
            simplify
            solve integral INT x:[0, pi/2]. exp(2*x)*cos(x)
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT x:[0,pi]. (x * sin(x)) ^ 2
            simplify
            rewrite sin(x) ^ 2 to (1 - cos(2*x)) / 2 using identity
            expand polynomial
            simplify
            integrate by parts with u = x^2 / 2, v = sin(2*x)
            simplify
            integrate by parts with u = x / 2, v = -cos(2*x)
            simplify
            apply integral identity
            simplify
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT x:[1, exp(1)]. sin(log(x))
            substitute u for log(x)
            integrate by parts with u = -exp(u), v = cos(u)
            simplify
            integrate by parts with u = exp(u), v = sin(u)
            simplify
            solve integral INT u:[0,1]. exp(u) * sin(u)
        """
        self.check_actions("base", "tongji", actions)

        actions = """
            calculate INT x:[1/exp(1), exp(1)]. abs(log(x))
            split region at 1
            simplify
            integrate by parts with u = log(x), v = x
            simplify
            integrate by parts with u = log(x), v = x
            apply integral identity
            simplify
        """
        self.check_actions("base", "tongji", actions)

    def testChapter1Section5(self):
        actions = """
            prove (INT x:[0,oo]. log(x) / (x ^ 2 + 1)) = 0
            lhs:
                split region at 1
                inverse substitute 1 / u for x creating u
                simplify
                rewrite u ^ 2 * (1 / u ^ 2 + 1) to u ^ 2 + 1
                simplify
        """
        self.check_actions("interesting", "chapter1section5", actions)

    def testChapter1Section7(self):
        actions = """
            prove (INT x:[0,1]. (x^4*(1-x)^4)/(1+x^2)) = 22/7 - pi
            lhs:
                rewrite (x^4*(1-x)^4)/(1+x^2) to (x^6-4*x^5+5*x^4-4*x^2+4)-4/(1+x^2)
                simplify
                apply integral identity
                simplify
        """
        self.check_actions("interesting", "chapter1section7", actions)

    def testEasy01(self):
        actions = """
            prove (INT x:[1,oo]. 1 / ((x+a)*sqrt(x-1))) = pi / sqrt(a+1) for a > -1
            lhs:
                substitute t for sqrt(x - 1)
                simplify
                substitute y for t / sqrt(a + 1)
                rewrite y ^ 2 * (a + 1) + a + 1 to (a + 1) * (y^2 + 1)
                apply integral identity
                simplify
        """
        self.check_actions("interesting", "easy01", actions)

    def testEasy02(self):
        actions = """
            prove (INT x:[0, oo]. log(1 + a^2 / x^2)) = a * pi for a > 0
            lhs:
                integrate by parts with u = log(1+a^2/x^2), v = x
                simplify
                rewrite x^2 * (a^2 / x^2 + 1) to a^2 + x^2
                apply integral identity
                simplify
        """
        self.check_actions("interesting", "easy02", actions)

    def testEasy03(self):
        actions = """
            prove (INT x:[0, oo]. log(x) / (x^2+b^2)) = pi * log(b) / (2*b) for b > 0
            lhs:
                inverse substitute 1/t for x creating t
                rewrite log(1/t) to -log(t)
                rewrite -log(t) / ((1/t)^2 + b^2) * -(1/t^2) to log(t) / (1 + b^2*t^2)
                inverse substitute s/b for t creating s
                rewrite log(s/b) to log(s) - log(b)
                expand polynomial
                apply integral identity
                simplify
        """
        self.check_actions("interesting", "easy03", actions)

    def testEasy04(self):
        actions = """
            prove (INT x:[1,oo]. log(x) / (x+1)^2) = log(2)
            subgoal 1: (INT x:[0,oo]. 1 / (1 + exp(a*x))) = log(2) / a for a > 0
            lhs:
                substitute u for exp(a * x)
                simplify
                rewrite 1 / (u * (u+1)) to 1/u - 1/(u+1)
                simplify
                substitute y for u + 1
                apply integral identity
                simplify
            done
            subgoal 2: (INT x:[1,oo]. log(x) / (a ^ 2 * (x + 1) ^ 2)) = log(2) / a^2 for a > 0
            from 1:
                differentiate both sides at a
                simplify
                substitute y for exp(a * x)
                solve equation for INT y:[1,oo]. log(y) / (a ^ 2 * (y + 1) ^ 2)
            done
            lhs:
                rewrite (x+1) ^ 2 to 1^2 * (x+1)^2
                apply 2 on INT x:[1,oo]. log(x) / (1 ^ 2 * (x + 1) ^ 2)
                simplify
        """
        self.check_actions("interesting", "easy04", actions)

    def testEasy05(self):
        actions = """
            prove (INT x:[sqrt(2),oo]. 1 / (x + x ^ sqrt(2))) = (1 + sqrt(2)) * log(1 + 2 ^ (1/2 * (1 - sqrt(2))))
            lhs:
                rewrite 1 / (x + x ^ sqrt(2)) to x ^ -sqrt(2) / (x ^ (1 - sqrt(2)) + 1)
                substitute u for log(x ^ (1 - sqrt(2)) + 1)
                apply integral identity
                simplify
                rewrite -sqrt(2) + 1 to -1 / (1 + sqrt(2)) (at 2)
                simplify
                rewrite sqrt(2) to 2 ^ (1/2) (at 2)
                rewrite 2 ^ (1/2) ^ (-sqrt(2) + 1) to 2 ^ (1/2 * (-sqrt(2) + 1)) using identity
                rewrite (sqrt(2) + 1) * log(2 ^ (1/2 * (-sqrt(2) + 1)) + 1) to (1 + sqrt(2)) * log(1 + 2 ^ (1/2 * (1 - sqrt(2))))
        """
        self.check_actions("interesting", "easy05", actions)

    def testEasy06(self):
        actions = """
            prove (INT x:[-oo,oo]. 1 / cosh(x)) = pi
            lhs:
                expand definition for cosh (all)
                substitute t for exp(x)
                rewrite -log(t) to log(1 / t)
                simplify
                rewrite 1 / (t * (1 / t + t)) to 1 / (1 + t ^ 2)
                simplify
                apply integral identity
                simplify
        """
        self.check_actions("interesting", "easy06", actions)

    def testTrick2a(self):
        actions = """
            calculate INT x:[0,pi / 2]. sqrt(sin(x)) / (sqrt(sin(x)) + sqrt(cos(x)))
            substitute y for pi / 2 - x
            rewrite sqrt(cos(y)) / (sqrt(cos(y)) + sqrt(sin(y))) to 1 - sqrt(sin(y)) / (sqrt(cos(y)) + sqrt(sin(y)))
            apply integral identity
            solve integral INT x:[0,pi / 2]. sqrt(sin(x)) / (sqrt(sin(x)) + sqrt(cos(x)))
        """
        self.check_actions("interesting", "Trick2a", actions)

    def testTrick2b(self):
        actions = """
            calculate INT x:[0,pi]. x * sin(x) / (1 + cos(x) ^ 2)
            substitute y for pi - x
            expand polynomial
            simplify
            solve integral INT x:[0,pi]. x * sin(x) / (1 + cos(x) ^ 2)
            substitute u for cos(y)
            apply integral identity
            simplify
        """
        self.check_actions("interesting", "Trick2b", actions)

    def testTrick2c(self):
        actions = """
            prove (INT x:[0,pi / 2]. sin(x) ^ 2 / (sin(x) + cos(x))) = sqrt(2) / 4 * log(3 + 2 * sqrt(2))

            subgoal 1: (INT x:[0,pi / 2]. sin(x) ^ 2 / (sin(x) + cos(x))) = (INT x:[0,pi / 2]. cos(x) ^ 2 / (sin(x) + cos(x)))
            lhs:
                substitute y for pi / 2 - x
            done

            subgoal 2: (INT x:[0,pi / 2]. sin(x) ^ 2 / (sin(x) + cos(x))) = 1/2 * (INT x:[0,pi / 2]. 1 / (sin(x) + cos(x)))
            rhs:
                simplify
                rewrite 1 to sin(x) ^ 2 + cos(x) ^ 2
                expand polynomial
                simplify
                rewrite cos(x) + sin(x) to sin(x) + cos(x) (at 1)
                apply 1 on INT x:[0,pi / 2]. cos(x) ^ 2 / (sin(x) + cos(x))
                simplify
            done

            lhs:
                apply 2 on INT x:[0,pi / 2]. sin(x) ^ 2 / (sin(x) + cos(x))
                substitute z for tan(x / 2)
                simplify
                rewrite (-(z ^ 2) + 1) / (z ^ 2 + 1) + 2 * z / (z ^ 2 + 1) to (2 - (z - 1) ^ 2) / (z ^ 2 + 1)
                rewrite (z ^ 2 + 1) * ((2 - (z - 1) ^ 2) / (z ^ 2 + 1)) to 2 - (z - 1) ^ 2
                rewrite 2 - (z - 1) ^ 2 to (sqrt(2) + (z - 1)) * (sqrt(2) - (z - 1))
                rewrite 1 / ((sqrt(2) + (z - 1)) * (sqrt(2) - (z - 1))) to sqrt(2) / 4 * (1 / (sqrt(2) + (z - 1)) + 1 / (sqrt(2) - (z - 1)))
                simplify
                substitute u for sqrt(2) + 1 - z (at 1)
                substitute u for sqrt(2) - 1 + z (at 2)
                apply integral identity
                simplify
                rewrite sqrt(2) * (log(sqrt(2) + 1) - log(sqrt(2) - 1)) / 4 to 1/4 * sqrt(2) * (log(sqrt(2) + 1) - log(sqrt(2) - 1))
                rewrite log(sqrt(2) + 1) - log(sqrt(2) - 1) to log((sqrt(2) + 1) / (sqrt(2) - 1)) using identity
                rewrite (sqrt(2) + 1) / (sqrt(2) - 1) to 3 + 2 * sqrt(2)
                simplify
        """
        self.check_actions("interesting", "Trick2c", actions)

    def testTrick2d(self):
        actions = """
            prove (INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)) = pi / 8 * log(2)
            subgoal 1: (INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)) = (INT x:[0,pi / 4]. log(tan(x) + 1))
            lhs:
                inverse substitute tan(u) for x creating u
                rewrite sec(u) ^ 2 to tan(u) ^ 2 + 1 using identity
                simplify
            done
            subgoal 2: (INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)) = pi / 4 * log(2) - (INT x:[0,1]. log(x + 1) / (x ^ 2 + 1))
            lhs:
                apply 1 on INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)
                inverse substitute pi / 4 - y for x creating y
                simplify
                rewrite tan(pi / 4 - y) to (tan(pi / 4) - tan(y)) / (1 + tan(pi / 4) * tan(y)) using identity
                simplify
                rewrite (-tan(y) + 1) / (tan(y) + 1) + 1 to 2 / (1 + tan(y))
                rewrite log(2 / (1 + tan(y))) to log(2) - log(1 + tan(y)) using identity
                apply integral identity
                simplify
                apply 1 on INT x:[0,pi / 4]. log(tan(x) + 1)
            done
            from 2:
                solve equation for INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)
        """
        self.check_actions("interesting", "Trick2d", actions)

    def testTrick2e(self):
        actions = """
            prove (INT t:[0,a]. log(t + a) / (t ^ 2 + a ^ 2)) = pi / (8 * a) * log(2 * a ^ 2) for a > 0
            subgoal 1: (INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)) = a * (INT t:[0,a]. log(t + a) / (t ^ 2 + a ^ 2)) - pi / 4 * log(a) for a > 0
            lhs:
                inverse substitute t / a for x creating t
                simplify
                rewrite 1 / (t ^ 2 / a ^ 2 + 1) * log(t / a + 1) to log(t / a + 1) * a ^ 2 / (t ^ 2 + a ^ 2)
                rewrite t / a + 1 to (t + a) / a
                simplify
                rewrite log((a + t) / a) to log(a + t) - log(a) using identity
                rewrite 1 / (a ^ 2 + t ^ 2) * (log(a + t) - log(a)) to log(a + t) / (a ^ 2 + t ^ 2) - log(a) / (a ^ 2 + t ^ 2)
                simplify
                apply integral identity
                simplify
                expand polynomial
            done
            subgoal 2: a * (INT t:[0,a]. log(t + a) / (t ^ 2 + a ^ 2)) - pi / 4 * log(a) = pi * log(2) / 8 for a > 0
            lhs:
                apply 1 on a * (INT t:[0,a]. log(t + a) / (t ^ 2 + a ^ 2)) - pi / 4 * log(a)
                apply integral identity
            done
            from 2:
                solve equation for INT t:[0,a]. log(t + a) / (t ^ 2 + a ^ 2)
                rewrite pi * log(a) / 4 to 1/8 * pi * (2 * log(a))
                rewrite 2 * log(a) to log(a ^ 2)
                rewrite 1/8 * pi * log(a ^ 2) + pi * log(2) / 8 to 1/8 * pi * (log(2) + log(a ^ 2))
                rewrite log(2) + log(a ^ 2) to log(2 * a ^ 2)
                rewrite 1 / a * (1/8 * pi * log(2 * a ^ 2)) to pi / (8 * a) * log(2 * a ^ 2)
        """
        self.check_actions("interesting", "Trick2e", actions)

    def testPartialFraction(self):
        actions = """
            prove (INT x:[0,oo]. 1 / (x ^ 4 + 2 * x ^ 2 * cosh(2 * a) + 1)) = pi / (4 * cosh(a))
            lhs:
                expand definition for cosh (all)
                rewrite x ^ 4 + 2 * x ^ 2 * ((exp(-(2 * a)) + exp(2 * a)) / 2) + 1 to (x ^ 2 + exp(2 * a)) * (x ^ 2 + exp(-(2 * a)))
                rewrite 1 / ((x ^ 2 + exp(2 * a)) * (x ^ 2 + exp(-(2 * a)))) to 1 / (exp(2 * a) - exp(-(2 * a))) * (1 / (x ^ 2 + exp(-(2 * a))) - 1 / (x ^ 2 + exp(2 * a)))
                simplify
                rewrite exp(-(2 * a)) to exp(-a) ^ 2
                rewrite exp(-(2 * a)) to exp(-a) ^ 2
                rewrite exp(2 * a) to exp(a) ^ 2
                rewrite exp(2 * a) to exp(a) ^ 2
                apply integral identity
                simplify
                rewrite to pi / (4 * ((exp(a) + exp(-a)) / 2))
                fold definition for cosh (all)
        """
        self.check_actions("interesting", "partialFraction", actions)

    def testLeibniz01(self):
        actions = """
            prove (INT x:[0,oo]. 1 / (x ^ 2 + a ^ 2) ^ 3) = 3 * pi / (16 * a ^ 5) for a > 0
            subgoal 1: (INT x:[0,oo]. 1 / (x ^ 2 + a ^ 2)) = pi / (2 * a) for a > 0
            lhs:
                inverse substitute a * u for x creating u
                simplify
                rewrite 1 / (a ^ 2 * u ^ 2 + a ^ 2) to 1 / (a ^ 2 * (u ^ 2 + 1))
                simplify
                apply integral identity
                simplify
            done

            subgoal 2: (INT x:[0,oo]. 1 / (x ^ 2 + a ^ 2) ^ 2) = pi / (4 * a ^ 3) for a > 0
            from 1:
                differentiate both sides at a
                simplify
                solve equation for INT x:[0,oo]. 1 / (a ^ 2 + x ^ 2) ^ 2
            done

            from 2:
                differentiate both sides at a
                simplify
                solve equation for INT x:[0,oo]. 1 / (a ^ 2 + x ^ 2) ^ 3
        """
        self.check_actions("interesting", "Leibniz01", actions)

    def testLeibniz02(self):
        actions = """
            prove (INT x:[0,1]. 1 / sqrt(-log(x))) = sqrt(pi)
            define g(t) = (INT x:[0,t]. exp(-(x ^ 2) / 2)) ^ 2
            subgoal 1: (INT x:[-oo,oo]. exp(-(x ^ 2) / 2)) = 2 * (LIM {t -> oo}. sqrt(g(t)))
            lhs:
                split region at 0
                substitute y for -x
                substitute x for y
                simplify
            rhs:
                expand definition for g (all)
                simplify
            done
            subgoal 2: (D t. g(t) + 2 * (INT y:[0,1]. exp(-(1 + y ^ 2) * t ^ 2 / 2) / (1 + y ^ 2))) = 0 for t > 0
            lhs:
                expand definition for g (all)
                simplify
                substitute y for x / t (at 2)
                rewrite exp(t ^ 2 * (-(y ^ 2) - 1) / 2) to exp(1/2 * t ^ 2 * (-(y ^ 2) - 1))
                rewrite 1/2 * t ^ 2 * (-(y ^ 2) - 1) to -1/2 * t ^ 2 * y ^ 2 + 1/2 * t ^ 2 * -1
                simplify
                rewrite exp(-(t ^ 2 * y ^ 2 / 2) - t ^ 2 / 2) to exp(-1/2 * t ^ 2 * y ^ 2) * exp(-1/2 * t ^ 2) using identity
                simplify
                rewrite (-(y ^ 2) - 1) / (y ^ 2 + 1) to -1
                simplify
            done
            subgoal 3: 2 * (INT y:[0,1]. exp(1/2 * t ^ 2 * (-(y ^ 2) - 1)) * (y ^ 2 + 1) ^ (-1)) + g(t) = SKOLEM_CONST(C) for t > 0
            from 2:
                integrate both sides
                apply indefinite integral
                simplify
            done
            subgoal 4: pi / 2 = SKOLEM_CONST(C)
            from 3:
                apply limit t -> 0 both sides
                simplify
                expand definition for g (all)
                apply integral identity
                simplify
            done
            subgoal 5: g(t) = -(2 * (INT y:[0,1]. 1 / (y ^ 2 + 1) * exp(t ^ 2 * (-(y ^ 2) - 1) / 2))) + pi / 2 for t > 0
            lhs:
                apply 3 on g(t)
                apply 4 on SKOLEM_CONST(C)
                simplify
            done
            subgoal 6: (INT x:[-oo,oo]. exp(-(x ^ 2) / 2)) = sqrt(2 * pi)
            lhs:
                apply 1 on INT x:[-oo,oo]. exp(-(x ^ 2) / 2)
                apply 5 on g(t)
                simplify
            rhs:
                simplify
            done
            subgoal 7: (INT x:[0,oo]. exp(-(x ^ 2) / 2)) = sqrt(2) * sqrt(pi) / 2
            from 6:
                split region at 0
                substitute y for -x
                substitute x for y
                simplify
                solve equation for INT x:[0,oo]. exp(-(x ^ 2 / 2))
            done
            subgoal 8: (INT x:[-oo,oo]. exp(-(a * x ^ 2))) = sqrt(pi / a) for a > 0
            lhs:
                substitute u for sqrt(2 * a) * x
                simplify
                substitute x for u
                rewrite -(x ^ 2 / 2) to -(x ^ 2) / 2
                apply 6 on INT x:[-oo,oo]. exp(-(x ^ 2) / 2)
                simplify
            rhs:
                simplify
            done
            subgoal 9: (INT x:[0,oo]. exp(-(x ^ 2))) = sqrt(pi) / 2
            from 7:
                substitute x for x / sqrt(2)
                simplify
                solve equation for INT x:[0,oo]. exp(-(x ^ 2))
            done
            from 9:
                substitute t for exp(-(x ^ 2))
                simplify
                solve equation for INT x:[0,1]. 1 / sqrt(-log(x))
        """
        self.check_actions("interesting", "Leibniz02", actions)

    def testEulerLogSineIntegral(self):
        actions = """
            prove (INT x:[0,pi / 2]. log(a * sin(x))) = pi / 2 * log(a / 2) for a > 0
            define I(a) = INT x:[0,pi/2]. log(a * sin(x)) for a > 0
            define J(a) = INT x:[0,pi/2]. log(a * sin(2*x)) for a > 0
            subgoal 1: J(a) = I(a) for a > 0
            lhs:
                expand definition for J
                substitute t for 2 * x
                split region at pi / 2
                substitute x for pi - t (at 2)
                simplify
            rhs:
                expand definition for I
            done
            subgoal 2: J(a) = pi / 2 * log(2 / a) + 2 * I(a) for a > 0
            lhs:
                expand definition for J
                rewrite sin(2 * x) to 2 * sin(x) * cos(x) using identity
                rewrite a * (2 * sin(x) * cos(x)) to 2 / a * (a * sin(x)) * (a * cos(x))
                rewrite log(2 / a * (a * sin(x)) * (a * cos(x))) to log(2 / a * (a * sin(x))) + log(a * cos(x)) using identity
                rewrite log(2 / a * (a * sin(x))) to log(2 / a) + log(a * sin(x)) using identity
                apply integral identity
                simplify
                substitute t for pi / 2 - x
                substitute x for t
                simplify
            rhs:
                expand definition for I (all)
                simplify
            done
            lhs:
                fold definition for I (all)
                apply 1 on I(a)
                apply 2 on J(a)
                solve integral I(a)
                rewrite log(2 / a) to log(2) + log(1 / a) using identity
                expand polynomial
                simplify
            rhs:
                rewrite log(a / 2) to log(a) - log(2) using identity
                expand polynomial
        """
        self.check_actions("interesting", "euler_log_sin", actions)

    def testEulerLogSineIntegral02(self):
        actions = """
            prove (INT x:[0,pi / 2]. log(sin(x) / x)) = pi / 2 * (1 - log(pi))
            lhs:
                rewrite log(sin(x) / x) to log(sin(x)) - log(x) using identity
                simplify
                rewrite log(sin(x)) to log(1 * sin(x))
                apply integral identity
                integrate by parts with u = log(x), v = x
                apply integral identity
                simplify
                expand polynomial
                simplify
            rhs:
                expand polynomial
        """
        self.check_actions("interesting", "euler_log_sin02", actions)

    def testEulerLogSineIntegral0304(self):
        actions = """
            prove (INT x:[0,1]. log(x + 1 / x) / (x ^ 2 + 1)) = pi / 2 * log(2)
            subgoal 1: (INT x:[0,oo]. log(x ^ 2 + 1) / (x ^ 2 + 1)) = pi * log(2)
            lhs:
                inverse substitute tan(u) for x creating u
                rewrite sec(u) ^ 2 to tan(u) ^ 2 + 1 using identity
                simplify
                rewrite tan(u) ^ 2 + 1 to sec(u) ^ 2 using identity
                rewrite sec(u) to cos(u) ^ (-1) using identity
                simplify
                substitute x for pi / 2 - u
                rewrite sin(x) to 1 * sin(x)
                apply integral identity
                simplify
            done

            from 1:
                split region at 1
                substitute y for 1 / x (at 2)
                rewrite y ^ 2 * (1 / y ^ 2 + 1) to y ^ 2 + 1
                rewrite 1 / (y ^ 2 + 1) * log(1 / y ^ 2 + 1) to log(1 / y ^ 2 + 1) / (y ^ 2 + 1)
                rewrite (INT y:[0,1]. log(y ^ 2 + 1) / (y ^ 2 + 1)) + (INT y:[0,1]. log(1 / y ^ 2 + 1) / (y ^ 2 + 1)) to INT y:[0,1]. log(y ^ 2 + 1) / (y ^ 2 + 1) + log(1 / y ^ 2 + 1) / (y ^ 2 + 1)
                rewrite log(y ^ 2 + 1) / (y ^ 2 + 1) + log(1 / y ^ 2 + 1) / (y ^ 2 + 1) to (log(y ^ 2 + 1) + log(1 / y ^ 2 + 1)) / (y ^ 2 + 1)
                rewrite log(y ^ 2 + 1) + log(1 / y ^ 2 + 1) to log((y ^ 2 + 1) * (1 / y ^ 2 + 1))
                rewrite (y ^ 2 + 1) * (1 / y ^ 2 + 1) to (y + 1 / y) ^ 2
                rewrite log((y + 1 / y) ^ 2) to 2 * log(y + 1 / y) using identity
                simplify
                rewrite 1 / (y ^ 2 + 1) * log(1 / y + y) to log(y + 1 / y) / (y ^ 2 + 1)
                solve equation for INT y:[0,1]. log(y + 1 / y) / (y ^ 2 + 1)
        """
        self.check_actions("interesting", "euler_log_sin0304", actions)

    def testEulerLogSineIntegral05(self):
        actions = """
            prove (INT x:[0,oo]. log(x) / (x ^ 2 - b * x + 1)) = 0 for b > -2, b < 2
            subgoal 1: x ^ 2 - b * x + 1 != 0 for b > -2, b < 2
            lhs:
                rewrite x ^ 2 - b * x + 1 to (x - 1/2 * b) ^ 2 + 1 - 1/4 * b ^ 2
            done
            subgoal 2: (INT x:[0,oo]. log(x ^ a + 1) / (x ^ 2 - b * x + 1)) = (INT x:[0,oo]. log(x ^ a + 1) / (x ^ 2 - b * x + 1)) - a * (INT x:[0,oo]. log(x) / (x ^ 2 - b * x + 1)) for a > 0, b > -2, b < 2
            lhs:
                inverse substitute 1 / u for x creating u
                simplify
                expand polynomial
                rewrite (1 / u) ^ a to 1 ^ a / u ^ a using identity
                rewrite 1 ^ a / u ^ a + 1 to (1 + u ^ a) / u ^ a
                rewrite log((1 + u ^ a) / u ^ a) to log(1 + u ^ a) - log(u ^ a)
                expand polynomial
                simplify
            done
            from 2:
                solve equation for INT x:[0,oo]. log(x) / (x ^ 2 - b * x + 1)
        """
        self.check_actions("interesting", "euler_log_sin05", actions)

    def testEulerLogSineIntegral06(self):
        actions = """
            prove (INT x:[0,1]. (1 - x) / (1 + x + x ^ 2)) = sqrt(3) * pi / 6 - log(3) / 2
            lhs:
                rewrite 1 + x + x ^ 2 to (x + 1/2) ^ 2 + 3/4
                substitute u for 2 * (x + 1/2) / sqrt(3)
                rewrite 3 * u ^ 2 / 2 + 3/2 to 3/2 * (u ^ 2 + 1)
                simplify
                rewrite 1 / (u ^ 2 + 1) * (-(u * sqrt(3) / 2) + 3/2) to -sqrt(3) / 2 * (u / (u ^ 2 + 1)) + 3/2 * (1 / (u ^ 2 + 1))
                apply integral identity
                simplify
                substitute t for u ^ 2 + 1
                apply integral identity
                simplify
                expand polynomial
                simplify
        """
        self.check_actions("interesting", "euler_log_sin06", actions)


if __name__ == "__main__":
    unittest.main()
