"""Unit test for integrals using internal language."""

import unittest
import lark

from integral import compstate
from integral import action
from integral import parser

import os
os.chdir('E:\\=graduatelife======\\learn-git\\iscalc')

class ActionTest(unittest.TestCase):
    def check_actions(self, base_file, current_file, actions: str,
                      *, print_lines=False, print_state=False):
        file = compstate.CompFile(base_file, current_file)
        state = action.InitialState(file)
        actions = [s for s in actions.split('\n') if s.strip()]
        for act in actions:
            if print_lines:
                print(act)
            if act.lstrip().startswith('#') or act.lstrip().startswith('//'):
                # title or comment
                continue
            a = parser.parse_action(act)
            state = state.process_action(a)
        if print_state:
            print(state)
        if not print_state and not isinstance(state, action.InitialState):
            raise AssertionError("Does not end in initial state")
        
    def testCalculationFinished(self):
        file = compstate.CompFile("base", "standard")
        state = action.InitialState(file)

        actions = """
            calculate INT x. (3 - 2*x)^3
                substitute u for 3 - 2*x
                apply integral identity
                simplify
        """
        actions = [s for s in actions.split('\n') if s.strip()]
        for act in actions:
            a = parser.parse_action(act)
            state = state.process_action(a)
        self.assertFalse(state.is_finished())

    def testCalculationFinished2(self):
        file = compstate.CompFile("base", "standard")
        state = action.InitialState(file)

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
            state = state.process_action(a)
        self.assertTrue(state.is_finished())

    def testStandard(self):
        actions = """
            prove (INT x. 1 / (x + a)) = log(abs(x + a)) + SKOLEM_CONST(C) for x + a != 0
            lhs:
                substitute u for x + a
                apply integral identity
                replace substitution
            done

            prove (INT x. exp(a * x)) = exp(a * x) / a + SKOLEM_CONST(C) for a != 0
            lhs:
                substitute u for a * x
                simplify
                apply integral identity
                replace substitution
            done

            prove (INT x. sin(a * x)) = -(cos(a * x) / a) + SKOLEM_CONST(C) for a != 0
            lhs:
                substitute u for a * x
                simplify
                apply integral identity
                replace substitution
                simplify
            done

            prove (INT x. cos(a * x)) = sin(a * x) / a + SKOLEM_CONST(C) for a != 0
            lhs:
                substitute u for a * x
                simplify
                apply integral identity
                replace substitution
            done

            prove (INT x. 1 / (a ^ 2 + x ^ 2)) = 1 / a * arctan(x / a) + SKOLEM_CONST(C) for a != 0
            lhs:
                substitute u for x / a
                rewrite a ^ 2 * u ^ 2 + a ^ 2 to a ^ 2 * (u ^ 2 + 1)
                simplify
                apply integral identity
                replace substitution
                simplify
            done

            prove (INT x. x ^ k * log(x)) = x ^ (k + 1) * log(x) / (k + 1) - x ^ (k + 1) / (k + 1) ^ 2 + SKOLEM_CONST(C) for x > 0, k != -1
            lhs:
                integrate by parts with u = log(x), v = x ^ (k + 1) / (k + 1)
                simplify
                apply integral identity
                simplify
            done

            prove (INT x:[0,1]. x ^ m * log(x) ^ n) = (-1) ^ n * factorial(n) / (m + 1) ^ (n + 1) for m >= 0, n >= 0, isInt(n)
            induction on n
            base:
                lhs:
                    apply integral identity
                    simplify
                rhs:
                    simplify
                done
            induct:
                lhs:
                    integrate by parts with u = log(x) ^ (n + 1), v = x ^ (m + 1) / (m + 1)
                    simplify
                    apply induction hypothesis (all)
                    simplify
                    rewrite to (-1) ^ (n + 1) * (m + 1) ^ (-n - 2) * ((n + 1) * factorial(n))
                    rewrite (n + 1) * factorial(n) to factorial(n + 1)
                    simplify
                rhs:
                    simplify
                done
            done

            prove (INT x:[0,oo]. exp(-(x * y)) * sin(a * x)) = a / (a ^ 2 + y ^ 2) for y > 0
            lhs:
                integrate by parts with u = exp(-(x * y)), v = -cos(a * x) / a
                simplify
                integrate by parts with u = exp(-(x * y)), v = sin(a * x) / a
                simplify
                solve integral INT x:[0,oo]. exp(-(x * y)) * sin(a * x)
                rewrite to a / (a ^ 2 + y ^ 2)
            done

            prove (INT x. a ^ x) = a ^ x / log(a) + SKOLEM_CONST(C) for a > 0, a != 1
            lhs:
                rewrite a ^ x to exp(log(a) * x)
                apply integral identity
                simplify
            done

            prove (INT x. cos(x) ^ 2) = 1/2 * (sin(2 * x) / 2 + x) + SKOLEM_CONST(C)
            lhs:
                rewrite cos(x) ^ 2 to (1 + cos(2 * x)) / 2
                apply integral identity
                simplify
            done
        """
        self.check_actions("base", "standard", actions)

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
        self.check_actions("MIT", "MIT2019", actions)

    def testLHopital(self):
        actions = """
            calculate LIM {x -> 1 }. (x ^ 2 - 1) / (x ^ 2 + 3 * x - 4)
                l'Hopital's rule
                simplify
            done

            calculate LIM {x -> 4 }. (x - 4) / (sqrt(x) - 2)
                l'Hopital's rule
                simplify
            done

            calculate LIM {x -> 0 }. sin(x) / x
                l'Hopital's rule
                simplify
            done
            
            calculate LIM {x -> 0 }. (3 ^ x - 2 ^ x) / (x ^ 2 - x)
                l'Hopital's rule
                simplify
            done

            calculate LIM {x -> 3 }. (1 / x - 1/3) / (x ^ 2 - 9)
                l'Hopital's rule
                simplify
            done

            calculate LIM {x -> 0 }. x * tan(x) / sin(3 * x)
                l'Hopital's rule
                simplify
            done
        """
        self.check_actions("UCDavis", "LHopital", actions)

    def testTongji(self):
        with open('theories/tongji05.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("base", "tongji", actions)

    def testTongjiIndefSubstitution(self):
        with open('theories/tongji0402.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("base", "tongji0402", actions)

    def testTongjiIndefByParts(self):
        with open('theories/tongji0403.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("base", "tongji0403", actions)

    def testTongjiIndefRational(self):
        with open('theories/tongji0404.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("base", "tongji0404", actions)

    def testUSubstitution(self):
        with open('theories/ucdavisUSubst.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("UCDavis", "USubstitution", actions)

    def testUCDavisPartialFraction(self):
        with open('theories/ucdavisPartial.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("UCDavis", "PartialFraction", actions)

    def testIntegrateByParts(self):
        with open('theories/ucdavisByParts.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("UCDavis", "IntegrateByParts", actions)

    def testExponential(self):
        with open('theories/ucdavisExponential.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("UCDavis", "Exponential", actions)

    def testTrigonometric(self):
        with open('theories/ucdavisTrigonometric.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("UCDavis", "Trigonometric", actions)

    def testLogAndArcTangent(self):
        with open('theories/ucdavisLogArctan.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("UCDavis", "LogAndArcTangent", actions)

    def testPowerSubstitution(self):
        with open('theories/ucdavisPowerSubst.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("UCDavis", "PowerSubstitution", actions)

    def testTrigSubstitution(self):
        with open('theories/ucdavisTrigSubst.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("UCDavis", "TrigSubstitution", actions)

    def testWallis(self):
        actions = """
            prove (INT x:[0,oo]. 1 / (x ^ 2 + b) ^ (m + 1)) = pi / 2 ^ (2 * m + 1) * binom(2 * m,m) * (1 / b ^ ((2 * m + 1) / 2)) for b > 0, m >= 0
            define I(m,b) = (INT x:[0,oo]. 1 / (x ^ 2 + b) ^ (m + 1)) for b > 0, m >= 0
            subgoal 1: (D b. I(m,b)) = -(m + 1) * I(m + 1,b) for b > 0, m >= 0
            lhs:
                expand definition for I (all)
                exchange derivative and integral
                simplify
            rhs:
                expand definition for I (all)
                simplify
            done

            subgoal 2: I(m,b) = pi / 2 ^ (2 * m + 1) * binom(2 * m,m) * (1 / b ^ ((2 * m + 1) / 2)) for b > 0, m >= 0
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
        self.check_actions("base", "wallis", actions)

    def testGammaFunction(self):
        actions = """
            define Gamma(n) = (INT x:[0,oo]. exp(-x) * x^(n-1)) for n > 0
            prove Gamma(n) = (n - 1) * Gamma(n - 1) for n > 1
            lhs:
                expand definition for Gamma
                integrate by parts with u = x ^ (n - 1), v = -exp(-x)
                simplify
            rhs:
                expand definition for Gamma (all)
            done

            prove Gamma(n) = factorial(n - 1) for n >= 1
            induction on n starting from 1
                base:
                lhs:
                    expand definition for Gamma
                    apply integral identity
                    simplify
                done
                induct:
                lhs:
                    apply Gamma(n) = (n - 1) * Gamma(n - 1) on Gamma(n + 1)
                    simplify
                    apply induction hypothesis (all)
                    rewrite n * factorial(n - 1) to factorial(n)
                done
            done

            calculate INT x:[0,oo]. exp(-(x ^ 3))
                substitute y for x ^ 3
                simplify
                rewrite exp(-y) / y ^ (2/3) to exp(-y) * y ^ (1/3 - 1)
                fold definition for Gamma (all)
                rewrite to (4/3 - 1) * Gamma(4/3 - 1)
                apply Gamma(n) = (n - 1) * Gamma(n - 1) on (4/3 - 1) * Gamma(4/3 - 1)
            done
        """
        self.check_actions("interesting", "GammaFunction", actions)

    def testChapter1Section5(self):
        actions = """
            prove (INT x:[0,oo]. log(x) / (x ^ 2 + 1)) = 0
            lhs:
                split region at 1
                substitute 1 / u for x
                simplify
                rewrite u ^ 2 * (1 / u ^ 2 + 1) to u ^ 2 + 1
                simplify
            done
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
            done
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
            done
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
            done
        """
        self.check_actions("interesting", "easy02", actions)

    def testEasy03(self):
        actions = """
            prove (INT x:[0, oo]. log(x) / (x^2+b^2)) = pi * log(b) / (2*b) for b > 0
            lhs:
                substitute 1/t for x
                rewrite log(1/t) to -log(t)
                rewrite -log(t) / ((1/t)^2 + b^2) * -(1/t^2) to log(t) / (1 + b^2*t^2)
                substitute s/b for t
                rewrite log(s/b) to log(s) - log(b)
                expand polynomial
                apply integral identity
                simplify
            done
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
            done
        """
        self.check_actions("interesting", "easy04", actions)

    def testEasy05(self):
        actions = """
            prove (INT x:[sqrt(2),oo]. 1 / (x + x ^ sqrt(2))) = (1 + sqrt(2)) * log(1 + 2 ^ (1/2 * (1 - sqrt(2))))
            lhs:
                rewrite 1 / (x + x ^ sqrt(2)) to x ^ -sqrt(2) / (x ^ (-sqrt(2) + 1) + 1)
                substitute u for log(x ^ (1 - sqrt(2)) + 1)
                apply integral identity
                simplify
                rewrite -sqrt(2) + 1 to -1 / (1 + sqrt(2)) (at 2)
                simplify
                rewrite sqrt(2) to 2 ^ (1/2) (at 2)
                rewrite 2 ^ (1/2) ^ (-sqrt(2) + 1) to 2 ^ (1/2 * (-sqrt(2) + 1))
                rewrite (sqrt(2) + 1) * log(2 ^ (1/2 * (-sqrt(2) + 1)) + 1) to (1 + sqrt(2)) * log(1 + 2 ^ (1/2 * (1 - sqrt(2))))
            done
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
            done
        """
        self.check_actions("interesting", "easy06", actions)

    def testTrick2a(self):
        actions = """
            calculate INT x:[0,pi / 2]. sqrt(sin(x)) / (sqrt(sin(x)) + sqrt(cos(x)))
                substitute y for pi / 2 - x
                rewrite sqrt(cos(y)) / (sqrt(cos(y)) + sqrt(sin(y))) to 1 - sqrt(sin(y)) / (sqrt(cos(y)) + sqrt(sin(y)))
                apply integral identity
                solve integral INT x:[0,pi / 2]. sqrt(sin(x)) / (sqrt(sin(x)) + sqrt(cos(x)))
            done
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
            done
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
                rewrite log(sqrt(2) + 1) - log(sqrt(2) - 1) to log((sqrt(2) + 1) / (sqrt(2) - 1))
                rewrite (sqrt(2) + 1) / (sqrt(2) - 1) to 3 + 2 * sqrt(2)
                simplify
            done
        """
        self.check_actions("interesting", "Trick2c", actions)

    def testTrick2d(self):
        actions = """
            prove (INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)) = pi / 8 * log(2)
            subgoal 1: (INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)) = (INT x:[0,pi / 4]. log(tan(x) + 1))
            lhs:
                substitute tan(u) for x
                rewrite sec(u) ^ 2 to tan(u) ^ 2 + 1
                simplify
            done
            subgoal 2: (INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)) = pi / 4 * log(2) - (INT x:[0,1]. log(x + 1) / (x ^ 2 + 1))
            lhs:
                apply 1 on INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)
                substitute pi / 4 - y for x
                simplify
                rewrite tan(pi / 4 - y) to (tan(pi / 4) - tan(y)) / (1 + tan(pi / 4) * tan(y))
                simplify
                rewrite (-tan(y) + 1) / (tan(y) + 1) + 1 to 2 / (1 + tan(y))
                rewrite log(2 / (1 + tan(y))) to log(2) - log(1 + tan(y))
                apply integral identity
                simplify
                apply 1 on INT x:[0,pi / 4]. log(tan(x) + 1)
            done
            from 2:
                solve equation for INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)
            done
        """
        self.check_actions("interesting", "Trick2d", actions)

    def testTrick2e(self):
        actions = """
            prove (INT t:[0,a]. log(t + a) / (t ^ 2 + a ^ 2)) = pi / (8 * a) * log(2 * a ^ 2) for a > 0
            subgoal 1: (INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)) = a * (INT t:[0,a]. log(t + a) / (t ^ 2 + a ^ 2)) - pi / 4 * log(a) for a > 0
            lhs:
                substitute t / a for x
                simplify
                rewrite 1 / (t ^ 2 / a ^ 2 + 1) * log(t / a + 1) to log(t / a + 1) * a ^ 2 / (t ^ 2 + a ^ 2)
                rewrite t / a + 1 to (t + a) / a
                simplify
                rewrite log((a + t) / a) to log(a + t) - log(a)
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
            done
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
            done
        """
        self.check_actions("interesting", "partialFraction", actions)

    def testPartialFraction03(self):#?
        # 输出<integral.action.CaseAnalysisState object at 0x0000021C0BB3A3B0>
        # Inside interesting integrals, Section 2.3, example 3
        actions = """
            prove (INT x:[0,oo]. 1/(x^4+2*x^2*cos(2*a)+1)) = pi/abs((4*cos(a))) for cos(a)!=0
            subgoal 1:(INT x:[0,oo]. x^2 / (x ^ 4 + 2 * x^2* cos(2 * a) + 1)) = (INT x:[0,oo]. 1/(x^4+2*x^2*cos(2*a)+1))
            rhs:
                substitute y for 1/x
                rewrite 1 / (y ^ 2 * (2 * cos(2 * a) / y ^ 2 + 1 / y ^ 4 + 1)) to (1/y ^ 2)/ (2 * cos(2 * a) / y ^ 2 + 1 / y ^ 4 + 1)
                rewrite 1 / y ^ 2 / (2 * cos(2 * a) / y ^ 2 + 1 / y ^ 4 + 1) to (y^4*(1 / y ^ 2)) / (y^4*(2 * cos(2 * a) / y ^ 2 + 1 / y ^ 4 + 1))
                rewrite y ^ 4 * (1 / y ^ 2) / (y ^ 4 * (2 * cos(2 * a) / y ^ 2 + 1 / y ^ 4 + 1)) to y^2/(y^4+2*y^2*cos(2*a)+1)
                substitute x for y
            done

            subgoal 2:(INT x:[0,oo]. 1/(x^4+2*x^2*cos(2*a)+1)) = 1/2*(INT x:[0,oo]. (1 + x^2)/(x^4+2*x^2*cos(2*a)+1))
            rhs:
                rewrite (1 + x^2)/(x^4+2*x^2*cos(2*a)+1) to (1/(x^4+2*x^2*cos(2*a)+1) + x^2/(x^4+2*x^2*cos(2*a)+1))
                simplify
                rewrite (2 * x ^ 2 * cos(2 * a) + x ^ 4 + 1) to (x ^ 4 + 2 * x^2* cos(2 * a) + 1)
                apply 1 on (INT x:[0,oo]. x ^ 2 / (x ^ 4 + 2 * x^2* cos(2 * a) + 1))
                simplify
            done


            subgoal 3:(INT x:[0,oo]. 1/(x^4+2*x^2*cos(2*a)+1)) = 1/4*(INT x:[-oo,oo]. (1 + x^2)/(x^4+2*x^2*cos(2*a)+1))
            rhs:
                split region at 0
                substitute u for -x
                substitute x for u
                simplify
                rewrite (INT x:[0,oo]. (x ^ 2 + 1) / (2 * x ^ 2 * cos(2 * a) + x ^ 4 + 1)) to (INT x:[0,oo]. (1 + x^2)/(x^4+2*x^2*cos(2*a)+1))
                apply 2 on (INT x:[0,oo]. (1 + x ^ 2) / (x ^ 4 + 2 * x ^ 2 * cos(2 * a) + 1))
                simplify
                rewrite to (INT x:[0,oo]. 1/(x^4+2*x^2*cos(2*a)+1))
            done
            subgoal 4:(INT x:[-oo,oo]. 2*x*sin(a) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1))) = -(INT x:[-oo,oo]. 2*x*sin(a) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1)))
            lhs:
                substitute u for -x
                substitute x for u
                rewrite (INT x:[-oo,oo]. -(2 * x * sin(a) / ((2 * x * sin(a) + x ^ 2 + 1) * (-(2 * x * sin(a)) + x ^ 2 + 1)))) to -(INT x:[-oo,oo]. 2*x*sin(a) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1)))
            done
            subgoal 5:(INT x:[-oo,oo]. 2*x*sin(a) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1))) = 0
            lhs:
                rewrite to 1/2*(INT x:[-oo,oo]. 2*x*sin(a) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1)))+1/2*(INT x:[-oo,oo]. 2*x*sin(a) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1)))
                apply 4 on (INT x:[-oo,oo]. 2*x*sin(a) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1)))
                rewrite to 1/2*((INT x:[-oo,oo]. 2*x*sin(a) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1))) - (INT x:[-oo,oo]. 2*x*sin(a) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1))))
                simplify
            done
            subgoal 6:(INT x:[-oo,oo]. (1 + x ^ 2) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1))) = (INT x:[-oo,oo]. (1 + 2*x*sin(a) + x ^ 2) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1)))
            rhs:
                expand polynomial
                rewrite (-(4 * x ^ 2 * sin(a) ^ 2) + 2 * x ^ 2 + x ^ 4 + 1) to ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1))
                simplify
                rewrite 2 * sin(a) * (INT x:[-oo,oo]. x / ((2 * x * sin(a) + x ^ 2 + 1) * (-(2 * x * sin(a)) + x ^ 2 + 1))) to (INT x:[-oo,oo]. (2*x*sin(a)) / ((2 * x * sin(a) + x ^ 2 + 1) * (-(2 * x * sin(a)) + x ^ 2 + 1)))
                rewrite ((2 * x * sin(a) + x ^ 2 + 1) * (-(2 * x * sin(a)) + x ^ 2 + 1)) to ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1))
                apply 5 on (INT x:[-oo,oo]. 2 * x * sin(a) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1)))
                rewrite to (INT x:[-oo,oo]. (x ^ 2 / (-(4 * x ^ 2 * sin(a) ^ 2) + 2 * x ^ 2 + x ^ 4 + 1) + 1 / (-(4 * x ^ 2 * sin(a) ^ 2) + 2 * x ^ 2 + x ^ 4 + 1)))
                rewrite x ^ 2 / (-(4 * x ^ 2 * sin(a) ^ 2) + 2 * x ^ 2 + x ^ 4 + 1) + 1 / (-(4 * x ^ 2 * sin(a) ^ 2) + 2 * x ^ 2 + x ^ 4 + 1) to (1 + x ^ 2) / (-(4 * x ^ 2 * sin(a) ^ 2) + 2 * x ^ 2 + x ^ 4 + 1)
                rewrite (-(4 * x ^ 2 * sin(a) ^ 2) + 2 * x ^ 2 + x ^ 4 + 1) to ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1))
            done
            subgoal 7:(INT x:[0,oo]. 1/(x^4+2*x^2*cos(2*a)+1)) = pi/(4*cos(a)) for cos(a)>0
            lhs:
                apply 3 on (INT x:[0,oo]. 1/(x^4+2*x^2*cos(2*a)+1))
                rewrite cos(2*a) to 1 - 2*(sin(a))^2
                rewrite 2 * x ^ 2 * (1 - 2 * sin(a) ^ 2) to 2*x^2 - 4*x^2*sin(a)^2
                rewrite (x ^ 4 + (2 * x ^ 2 - 4 * x ^ 2 * sin(a) ^ 2) + 1) to (x^2 - 2*x*sin(a)+1)*(x^2+2*x*sin(a)+1)
                apply 6 on (INT x:[-oo,oo]. (1 + x ^ 2) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1)))
                rewrite (1 + 2 * x * sin(a) + x ^ 2) to (x ^ 2 + 2 * x * sin(a) + 1)
                rewrite (x ^ 2 + 2 * x * sin(a) + 1) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1)) to 1 / (x ^ 2 - 2 * x * sin(a) + 1)
                rewrite 1 to sin(a)^2 + cos(a)^2
                rewrite 1 to sin(a)^2 + cos(a)^2
                rewrite (x ^ 2 - 2 * x * sin(a) + (sin(a) ^ 2 + cos(a) ^ 2)) to (x ^ 2 - 2 * x * sin(a) + sin(a) ^ 2 + cos(a) ^ 2)
                rewrite x ^ 2 - 2 * x * sin(a) + sin(a) ^ 2 to (x-sin(a))^2
                rewrite sin(a)^2 + cos(a)^2 to 1
                substitute u for (x - sin(a))
                apply integral identity
                simplify
                rewrite to 1 / (4 * cos(a))*((LIM {u -> oo}. arctan(u / cos(a)))-(LIM {u -> oo}. arctan(-(u / cos(a)))))
                rewrite arctan(-(u / cos(a))) to -arctan((u / cos(a)))
                simplify
            done
            subgoal 8:(INT x:[0,oo]. 1/(x^4+2*x^2*cos(2*a)+1)) = pi/(4*cos(a)) for cos(a)<0
            lhs:
                apply 3 on (INT x:[0,oo]. 1/(x^4+2*x^2*cos(2*a)+1))
                rewrite cos(2*a) to 1 - 2*(sin(a))^2
                rewrite 2 * x ^ 2 * (1 - 2 * sin(a) ^ 2) to 2*x^2 - 4*x^2*sin(a)^2
                rewrite (x ^ 4 + (2 * x ^ 2 - 4 * x ^ 2 * sin(a) ^ 2) + 1) to (x^2 - 2*x*sin(a)+1)*(x^2+2*x*sin(a)+1)
                apply 6 on (INT x:[-oo,oo]. (1 + x ^ 2) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1)))
                rewrite (1 + 2 * x * sin(a) + x ^ 2) to (x ^ 2 + 2 * x * sin(a) + 1)
                rewrite (x ^ 2 + 2 * x * sin(a) + 1) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1)) to 1 / (x ^ 2 - 2 * x * sin(a) + 1)
                rewrite 1 to sin(a)^2 + cos(a)^2
                rewrite 1 to sin(a)^2 + cos(a)^2
                rewrite (x ^ 2 - 2 * x * sin(a) + (sin(a) ^ 2 + cos(a) ^ 2)) to (x ^ 2 - 2 * x * sin(a) + sin(a) ^ 2 + cos(a) ^ 2)
                rewrite x ^ 2 - 2 * x * sin(a) + sin(a) ^ 2 to (x-sin(a))^2
                rewrite sin(a)^2 + cos(a)^2 to 1
                substitute u for (x - sin(a))
                apply integral identity
                simplify
                rewrite to 1 / (4 * cos(a))*((LIM {u -> oo}. arctan(u / cos(a)))-(LIM {u -> oo}. arctan(-(u / cos(a)))))
                rewrite arctan(-(u / cos(a))) to -arctan((u / cos(a)))
                simplify
            done
            case analysis on cos(a)
            case negative:
            lhs:
                apply 7 on (INT x:[0,oo]. 1/(x^4+2*x^2*cos(2*a)+1))
            done
            case positive:
            lhs:
                apply 8 on (INT x:[0,oo]. 1/(x^4+2*x^2*cos(2*a)+1))
            done
            """
        self.check_actions("interesting", "partialFraction03", actions)

    def testLeibniz01(self):
        actions = """
            prove (INT x:[0,oo]. 1 / (x ^ 2 + a ^ 2) ^ 3) = 3 * pi / (16 * a ^ 5) for a > 0
            subgoal 1: (INT x:[0,oo]. 1 / (x ^ 2 + a ^ 2)) = pi / (2 * a) for a > 0
            lhs:
                substitute a * u for x
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
            done
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
                rewrite exp(-(t ^ 2 * y ^ 2 / 2) - t ^ 2 / 2) to exp(-1/2 * t ^ 2 * y ^ 2) * exp(-1/2 * t ^ 2)
                simplify
                rewrite (-(y ^ 2) - 1) / (y ^ 2 + 1) to -1
                simplify
            done
            subgoal 3: 2 * (INT y:[0,1]. exp(1/2 * t ^ 2 * (-(y ^ 2) - 1)) * (y ^ 2 + 1) ^ (-1)) + g(t) = SKOLEM_CONST(C) for t > 0
            from 2:
                integrate both sides
                apply integral identity
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
            done
        """
        self.check_actions("interesting", "Leibniz02", actions)

    def testLeibniz03New(self):#?
        # apply integral identity无法计算 (INT x:[0,oo]. x ^ (-1/2) * exp(-x))
        # Overall goal: INT x:[0,oo]. cos(tx)*exp(-(x^2)/2) = sqrt(pi/2)*exp(-(t^2)/2)
        
        actions = """
            prove (INT x:[0,oo]. cos(t*x)*exp(-(x^2)/2)) = sqrt(pi/2)*exp(-(t^2)/2)
            define I(t) = (INT x:[0,oo]. cos(t*x)*exp(-(x^2)/2))
            subgoal 1:I(0) = sqrt(2*pi)
            lhs:
                expand definition for I
                substitute u for -x^2/2
                simplify
                substitute x for -u
                apply integral identity
            """
        self.check_actions("interesting", "leibniz03_new", actions)

    def testGaussianPowerExp(self):# ?
        # Inside interesting integrals, Section 2.3
        # apply integral identity无法计算 (INT x:[0,oo]. x ^ (-1/2) * exp(-x))
        actions = """
            prove (INT x:[0, oo]. x^(2*n) * exp(-x^2)) = factorial(2*n)/(4^n*factorial(n))*(1/2)*sqrt(pi) for isInt(n)
            define I(n) = (INT x:[0, oo]. x^(2*n) * exp(-x^2)) for n>=0,isInt(n)
            subgoal 1:(INT x:[0, oo]. (D x. x^(2*n-1)*exp(-x^2))) = 0 for n>=1,isInt(n)
            lhs:
                simplify
            done
            subgoal 2:(D x. x^(2*n-1)*exp(-x^2)) = (2*n-1)*x^(2*n-2)*exp(-x^2)-2*x^(2*n)*exp(-x^2)
            lhs:
                simplify
                rewrite x ^ (2 * n - 2) * exp(-(x ^ 2)) * (2 * n - 1) to (2*n-1)*x^(2*n-2)*exp(-x^2) 
                rewrite x ^ (2 * n - 1) to x^(2*n)/x
                rewrite 2 * (x ^ (2 * n) / x) * x * exp(-(x ^ 2)) to 2 *x ^ (2 * n) / x * x * exp(-(x ^ 2))
                rewrite 2 * x ^ (2 * n) / x * x * exp(-(x ^ 2)) to 2 * x ^ (2 * n) * exp(-(x ^ 2))
            done
            subgoal 3:(INT x:[0, oo]. x^(2*n) * exp(-x^2)) = I(n) for n>=0,isInt(n)
            rhs:
                expand definition for I
            done
            subgoal 4:(INT x:[0, oo]. (D x. x^(2*n-1)*exp(-x^2))) = (2*n-1)*I(n-1) - 2 * I(n) for n>=1,isInt(n)
            lhs:
                apply 2 on (D x. x^(2*n-1)*exp(-x^2))
                simplify
                apply 3 on (INT x:[0,oo]. x ^ (2 * n) * exp(-(x ^ 2)))
                rewrite (INT x:[0,oo]. x ^ (2 * n - 2) * exp(-(x ^ 2))) to (INT x:[0,oo]. x ^ (2 * (n - 1)) * exp(-(x ^ 2)))
                apply 3 on (INT x:[0,oo]. x ^ (2 * (n - 1)) * exp(-(x ^ 2)))
            done
            subgoal 5:I(n) = I(n-1)*(2*n*(2*n-1))/(4*n) for n>=0,isInt(n)
            from 4:
                apply 1 on (INT x:[0, oo]. (D x. x^(2*n-1)*exp(-x^2)))
                solve equation for I(n)
                rewrite I(n - 1) * (2 * n - 1) / 2 to I(n - 1) * ((2*n)*(2 * n - 1)) / (2*(2*n))
                rewrite (2 * (2 * n)) to (4 * n)
            done
            subgoal 6 :I(n) = I(0)*factorial(2*n)/(4^n*factorial(n)) for n>=1,isInt(n)
            induction on n starting from 1
                base:
                lhs:
                    apply 5 on I(1)
                    simplify
                done
                induct:
                lhs:
                    apply 5 on I(n+1)
                    simplify
                    apply induction hypothesis(all)
                    rewrite (4 * n + 4) to 2*(2*n+2)
                    rewrite I(0) * factorial(2 * n) / (4 ^ n * factorial(n)) * (2 * n + 1) * (2 * n + 2) / (2 * (2 * n + 2)) to I(0) * factorial(2 * n) / (4 ^ n * factorial(n)) * (2 * n + 1)/ 2
                done
            done
            subgoal 7:I(0) = 1/2*sqrt(pi)
            lhs:
                expand definition for I(all)
                substitute u for -(x^2)
                apply integral identity
                substitute -x for u
                simplify
                rewrite exp(-x) / sqrt(x) to x^(-1/2)*exp(-x)
                apply integral identity
            """
        self.check_actions("interesting", "gaussianPowerExp", actions)

    def testEulerLogSineIntegral(self):#?+
        # Inside interesting integrals, Section 2.4
        actions = """
            prove (INT x:[0,pi/2]. log(a * sin(x))) = pi/2 * log(a/2) for a>0
            subgoal 1:(INT x:[0,pi/2]. log(a * sin(x))) = (INT x:[0,pi/2]. log(a * cos(x)))
            lhs:
                substitute y for pi/2-x
            done
            subgoal 2:(INT x:[0,pi/2]. log(a * sin(2*x))) = (INT x:[0,pi/2]. log(a * sin(x)))
            lhs:
                substitute t for 2*x
                simplify
                split region at pi/2
                simplify
                substitute u for pi-t
                simplify
                substitute x for pi-u
            done
            subgoal 3:2*cos(x)*sin(x) = sin(2*x)
            rhs:
                rewrite to 2*cos(x)*sin(x)
            done
            subgoal 4:(INT x:[0,pi/2]. log(a * sin(x)))=1/2 * (INT x:[0,pi / 2]. log(a * sin(x))) + pi * log(a) / 4 - pi * log(2) / 4
            lhs:
                rewrite to 1/2*((INT x:[0,pi/2]. log(a * sin(x)))+(INT x:[0,pi/2]. log(a * sin(x))))
                apply 1 on (INT x:[0,pi/2]. log(a * sin(x)))
                rewrite to 1/2*(INT x:[0,pi/2]. (log(a * sin(x))+log(a*cos(x))))
                rewrite log(a*cos(x)) to log(a )+ log(cos(x))
                rewrite to 1/2 * (INT x:[0,pi / 2]. log(a * sin(x)) + log(a) + log(cos(x)))
                rewrite log(a * sin(x)) + log(a) to log(a * sin(x)*a)
                rewrite log(a * sin(x) * a) + log(cos(x)) to log(a * sin(x) * a*cos(x))
                rewrite to 1/2 * (INT x:[0,pi / 2]. log(a ^ 2 *1/2*(2 * cos(x) * sin(x))))
                apply 3 on (2 * cos(x) * sin(x))
                rewrite log(a ^ 2 * 1 / 2 * sin(2 * x)) to log(a*1/2*a*sin(2*x))
                rewrite log(a*1/2*a*sin(2*x)) to log(a*sin(2*x)*a*1/2)
                rewrite log(a*sin(2*x)*a*1/2) to log(a*sin(2*x)*a)+log(1/2)
                rewrite log(a * sin(2 * x) * a) to log(a * sin(2 * x))+log(a)
                apply integral identity
                simplify
                apply 2 on (INT x:[0,pi / 2]. log(a * sin(2 * x)))
            done
            subgoal 5:(INT x:[0,pi / 2]. log(a * sin(x))) = pi * log(a) / 2 - pi * log(2) / 2 
            from 4:
                solve equation for INT x:[0,pi / 2]. log(a * sin(x))
            done
            lhs:
                apply 5 on (INT x:[0,pi / 2]. log(a * sin(x)))
                rewrite to pi/2*(log(a)-log(2))
                rewrite log(a) - log(2) to log(a/2)
            done
            """
        self.check_actions("interesting", "euler_log_sin", actions)

    def testEulerLogSineIntegral02(self):
        actions = """
            prove (INT x:[0,pi / 2]. log(sin(x) / x)) = pi / 2 * (1 - log(pi))
            lhs:
                rewrite log(sin(x) / x) to log(sin(x)) - log(x)
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
            done
        """
        self.check_actions("interesting", "euler_log_sin02", actions)

    def testEulerLogSineIntegral0304(self):
        actions = """
            prove (INT x:[0,1]. log(x + 1 / x) / (x ^ 2 + 1)) = pi / 2 * log(2)
            subgoal 1: (INT x:[0,oo]. log(x ^ 2 + 1) / (x ^ 2 + 1)) = pi * log(2)
            lhs:
                substitute tan(u) for x
                rewrite sec(u) ^ 2 to tan(u) ^ 2 + 1
                simplify
                rewrite tan(u) ^ 2 + 1 to sec(u) ^ 2
                rewrite sec(u) to cos(u) ^ (-1)
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
                rewrite log((y + 1 / y) ^ 2) to 2 * log(y + 1 / y)
                simplify
                rewrite 1 / (y ^ 2 + 1) * log(1 / y + y) to log(y + 1 / y) / (y ^ 2 + 1)
                solve equation for INT y:[0,1]. log(y + 1 / y) / (y ^ 2 + 1)
            done
        """
        self.check_actions("interesting", "euler_log_sin0304", actions)

    def testEulerLogSineIntegral05(self):
        actions = """
            prove (INT x:[0,oo]. log(x) / (x ^ 2 - b * x + 1)) = 0 for b > -2, b < 2
            subgoal 1: x ^ 2 - b * x + 1 != 0
            lhs:
                rewrite x ^ 2 - b * x + 1 to (x - 1/2 * b) ^ 2 + 1 - 1/4 * b ^ 2
            done
            subgoal 2: (INT x:[0,oo]. log(x ^ a + 1) / (x ^ 2 - b * x + 1)) = (INT x:[0,oo]. log(x ^ a + 1) / (x ^ 2 - b * x + 1)) - a * (INT x:[0,oo]. log(x) / (x ^ 2 - b * x + 1)) for a > 0
            lhs:
                substitute 1 / u for x
                simplify
                expand polynomial
                rewrite (1 / u) ^ a to 1 ^ a / u ^ a
                rewrite 1 ^ a / u ^ a + 1 to (1 + u ^ a) / u ^ a
                rewrite log((1 + u ^ a) / u ^ a) to log(1 + u ^ a) - log(u ^ a)
                expand polynomial
                simplify
            done
            from 2:
                solve equation for INT x:[0,oo]. log(x) / (x ^ 2 - b * x + 1)
            done
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
            done
        """
        self.check_actions("interesting", "euler_log_sin06", actions)

    def testDirichletIntegral(self):#?
        # Inside interesting integrals, Section 3.2
        # 未实现sgn功能
        actions = """
            prove (INT x:[0,oo]. sin(a*x)/x) = pi/2 * sgn(a)
            define g(y,a) = INT x:[0,oo]. exp(-x * y) * sin(a * x) / x for y>0
            subgoal 1: (D y. g(y, a)) = - a / (a ^ 2 + y ^ 2) for y>0,a!=0
            lhs:
                expand definition for g(all)
                exchange derivative and integral
                simplify
                apply integral identity
            rhs:
                simplify
            done
            subgoal 2:g(y, a) = -arctan(y / a) + SKOLEM_FUNC(C(a)) for y>0,a>0
            from 1:
                integrate both sides
                apply integral identity
                simplify
            done
            subgoal 3:g(y, a) = -arctan(y / a) + SKOLEM_FUNC(C(a)) for y>0,a<0
            from 1:
                integrate both sides
                apply integral identity
                simplify
            done
            subgoal 4:(LIM {y -> oo}. g(y, a)) = 0 for y>0
            lhs:
                expand definition for g(all)
                simplify
            done
            subgoal 5:SKOLEM_FUNC(C(a)) = pi / 2 for a>0
            from 2:
                apply limit y -> oo both sides
                apply 4 on LIM {y -> oo}. g(y,a)
                simplify
                solve equation for SKOLEM_FUNC(C(a))
            done
            subgoal 6:SKOLEM_FUNC(C(a)) = -pi / 2 for a<0
            from 3:
                apply limit y -> oo both sides
                apply 4 on LIM {y -> oo}. g(y,a)
                simplify
                solve equation for SKOLEM_FUNC(C(a))
            done
            subgoal 7:g(y,a) = pi / 2 for a>0,y=0
            from 2:
                apply 5 on SKOLEM_FUNC(C(a)) 
                simplify
            done
            subgoal 8:g(y,a) = -pi / 2 for a<0,y=0
            from 3:
                apply 6 on SKOLEM_FUNC(C(a)) 
                simplify
            done
            define g(y,a) = INT x:[0,oo]. exp(-x * y) * sin(a * x) / x for y=0
            subgoal 9:g(y,a) = 0 for a=0,y=0
            lhs:
                expand definition for g(all)
                simplify
            """
        self.check_actions("interesting", "dirichletIntegral", actions)

    def testFlipside03(self):
        actions = """
            prove (INT x:[0,1]. (x ^ a - 1) / log(x)) = log(a + 1) for a > -1
            define I(a) = INT x:[0, 1]. (x ^ a - 1) / log(x) for a > -1
            subgoal 1: (D a. I(a)) = 1 / (a + 1) for a > -1
            lhs:
                expand definition for I (all)
                exchange derivative and integral
                simplify
                apply integral identity
                simplify
            done
            subgoal 2: I(a) = log(a + 1) + SKOLEM_CONST(C) for a > -1
            from 1:
                integrate both sides
                apply integral identity
                simplify
            done
            subgoal 3: SKOLEM_CONST(C) = 0
            from 2:
                substitute a for 0 in equation
                expand definition for I (all)
                simplify
                solve equation for SKOLEM_CONST(C)
            done
            lhs:
                fold definition for I (all)
                apply 2 on I(a)
                apply 3 on SKOLEM_CONST(C)
                simplify
            done
        """
        self.check_actions("interesting", "flipside03", actions)

    def testFlipside04(self):
        actions = """
            prove (INT x:[0,1]. (x ^ a - x ^ b) / log(x)) = log((a + 1) / (b + 1)) for a > -1, b > -1
            lhs:
                rewrite x ^ a - x ^ b to x ^ a - 1 - (x ^ b - 1)
                rewrite (x ^ a - 1 - (x ^ b - 1)) / log(x) to (x ^ a - 1) / log(x) - (x ^ b - 1) / log(x)
                simplify
                apply integral identity
                rewrite log(a + 1) - log(b + 1) to log((a + 1) / (b + 1))
            done
        """
        self.check_actions("interesting", "flipside04", actions)

    def testFlipside05(self):#? new_done
        # 报错StateException: Use done when goal is not finished
        actions = """
            prove (INT x:[0, oo]. exp(- (t * x))*(( cos(a * x) - cos(b * x) ) / x)) = log(sqrt((t^2+b^2)/(t^2+a^2))) for (t^2+b^2)/(t^2+a^2) > 0
            subgoal 1: (INT s:[a,b]. sin(x*s)) = (cos(a*x)-cos(b*x))/x
            lhs:
                apply integral identity
                simplify
                rewrite cos(a * x) / x - cos(b * x) / x to (cos(a*x)-cos(b*x))/x
            done
            subgoal 2: (INT x:[0,oo]. exp(-(t * x)) * sin(s * x)) = -(t ^ 2 / s ^ 2 * (INT x:[0,oo]. exp(-(t * x)) * sin(s * x))) + 1 / s for t>0,s>0
            lhs:
                integrate by parts with u=exp(-(t * x)),v=-1/s*cos(s*x)
                simplify
                integrate by parts with u=exp(-(t * x)),v=1/s*sin(s*x)
                simplify
            done
            subgoal 3:(INT x:[0,oo]. exp(-(t * x)) * sin(s * x)) = s / (t ^ 2 + s ^ 2)
            from 2:
                solve equation for INT x:[0,oo]. exp(-(t * x)) * sin(s * x)
                rewrite 1 / (s * (t ^ 2 / s ^ 2 + 1)) to s/(t^2+s^2)
            done
            subgoal 4:log(sqrt((b ^ 2 + t ^ 2) / (a ^ 2 + t ^ 2))) = 1/2*log((b ^ 2 + t ^ 2) / (a ^ 2 + t ^ 2))
            lhs:
                rewrite sqrt((b ^ 2 + t ^ 2) / (a ^ 2 + t ^ 2)) to ((b ^ 2 + t ^ 2) / (a ^ 2 + t ^ 2))^(1/2)
                simplify
            done
            lhs:
                apply 1 on (cos(a*x)-cos(b*x))/x
                rewrite INT x:[0,oo]. exp(-(t * x)) * (INT s:[a,b]. sin(x * s)) to INT x:[0,oo]. (INT s:[a,b]. exp(-(t * x)) * sin(s * x))
                exchange integral and integral
                apply 3 on INT x:[0,oo]. exp(-(t * x)) * sin(s * x)
                substitute u for t^2+s^2
                apply integral identity
                simplify
                rewrite -(log(a ^ 2 + t ^ 2) / 2) + log(b ^ 2 + t ^ 2) / 2 to -1/2*log(a ^ 2 + t ^ 2) + 1/2*log(b ^ 2 + t ^ 2)
                rewrite to 1/2*(log(b ^ 2 + t ^ 2) - log(a ^ 2 + t ^ 2))
                rewrite log(b ^ 2 + t ^ 2) - log(a ^ 2 + t ^ 2) to log((b ^ 2 + t ^ 2)/(a ^ 2 + t ^ 2))
                apply 4 on log((b ^ 2 + t ^ 2) / (a ^ 2 + t ^ 2))
                simplify
                rewrite to log(sqrt((t^2+b^2)/(t^2+a^2)))
            done
            """
        self.check_actions("interesting", "flipside05", actions)

    def testFlipside06(self):#? new_done
        # 报错StateException: Use done when goal is not finished
        actions = """
            prove (INT x:[0,oo]. (cos(a*x)-cos(b*x)) / x) = log(abs(b/a)) for a!=0,b!=0
            subgoal 1: (INT s:[a,b]. sin(x*s)) = (cos(a*x)-cos(b*x))/x for a!=0,b!=0
            lhs:
                apply integral identity
                simplify
                rewrite cos(a * x) / x - cos(b * x) / x to (cos(a*x)-cos(b*x))/x
            done
            subgoal 2: (INT x:[0,oo]. exp(-(t * x)) * sin(s * x)) = -(t ^ 2 / s ^ 2 * (INT x:[0,oo]. exp(-(t * x)) * sin(s * x))) + 1 / s for t>0,s>0
            lhs:
                integrate by parts with u=exp(-(t * x)),v=-1/s*cos(s*x)
                simplify
                integrate by parts with u=exp(-(t * x)),v=1/s*sin(s*x)
                simplify
            done
            subgoal 3:(INT x:[0,oo]. exp(-(t * x)) * sin(s * x)) = s / (t ^ 2 + s ^ 2)
            from 2:
                solve equation for INT x:[0,oo]. exp(-(t * x)) * sin(s * x)
                rewrite 1 / (s * (t ^ 2 / s ^ 2 + 1)) to s/(t^2+s^2)
            done
            subgoal 4:log(sqrt((b ^ 2 + t ^ 2) / (a ^ 2 + t ^ 2))) = 1/2*log((b ^ 2 + t ^ 2) / (a ^ 2 + t ^ 2))
            lhs:
                rewrite sqrt((b ^ 2 + t ^ 2) / (a ^ 2 + t ^ 2)) to ((b ^ 2 + t ^ 2) / (a ^ 2 + t ^ 2))^(1/2)
                simplify
            done

            subgoal 5: (INT x:[0, oo]. exp(-(t * x))*((cos(a*x)-cos(b*x))/x)) = log(sqrt((t^2+b^2)/(t^2+a^2))) for (t^2+b^2)/(t^2+a^2) > 0
            lhs:
                apply 1 on (cos(a*x)-cos(b*x))/x
                rewrite INT x:[0,oo]. exp(-(t * x)) * (INT s:[a,b]. sin(x * s)) to INT x:[0,oo]. (INT s:[a,b]. exp(-(t * x)) * sin(s * x))
                exchange integral and integral
                apply 3 on INT x:[0,oo]. exp(-(t * x)) * sin(s * x)
                substitute u for t^2+s^2
                apply integral identity
                simplify
                rewrite -(log(a ^ 2 + t ^ 2) / 2) + log(b ^ 2 + t ^ 2) / 2 to -1/2*log(a ^ 2 + t ^ 2) + 1/2*log(b ^ 2 + t ^ 2)
                rewrite to 1/2*(log(b ^ 2 + t ^ 2) - log(a ^ 2 + t ^ 2))
                rewrite log(b ^ 2 + t ^ 2) - log(a ^ 2 + t ^ 2) to log((b ^ 2 + t ^ 2)/(a ^ 2 + t ^ 2))
                apply 4 on log((b ^ 2 + t ^ 2) / (a ^ 2 + t ^ 2))
                simplify
                rewrite to log(sqrt((t^2+b^2)/(t^2+a^2)))
            done
            lhs:
                rewrite (cos(a*x)-cos(b*x)) / x to exp(-(0 * x))*((cos(a*x)-cos(b*x))/x)
                apply 5 on INT x:[0,oo]. exp(-(0 * x)) * ((cos(a * x) - cos(b * x)) / x)
                simplify
                rewrite log(abs(b) / abs(a)) to log(abs(b / a))
            done
            """
        self.check_actions("interesting", "flipside06", actions)

    def testFrullaniIntegral01(self):
        actions = """
            prove (INT x:[0,oo]. (arctan(a * x) - arctan(b * x)) / x) = pi * log(a) / 2 - pi * log(b) / 2 for a > 0, b > 0
            define I(a,b) = (INT x:[0,oo]. (arctan(a * x) - arctan(b * x)) / x) for a > 0, b > 0
            subgoal 1: (D a. I(a,b)) = pi / (2 * a) for a > 0, b > 0
            lhs:
                expand definition for I (all)
                exchange derivative and integral
                simplify
                substitute u for a * x
                apply integral identity
                simplify
            done
            subgoal 2: I(a,b) = pi * log(a) / 2 + SKOLEM_FUNC(C(b)) for a > 0, b > 0
            from 1:
                integrate both sides
                simplify
                apply integral identity
                simplify
            done
            subgoal 3: SKOLEM_FUNC(C(a)) = -(pi * log(a) / 2) for a > 0
            from 2:
                substitute b for a in equation
                solve equation for SKOLEM_FUNC(C(a))
                expand definition for I (all)
                simplify
            done
            lhs:
                fold definition for I (all)
                apply 2 on I(a,b)
                apply 3 on SKOLEM_FUNC(C(b))
                simplify
            done
        """
        self.check_actions("interesting", "FrullaniIntegral01", actions)

    def testCatalanConstant01(self):
        actions = """
            define G = SUM(n, 0, oo, (-1)^n / (2*n+1)^2)
            prove (INT x:[0,1]. arctan(x) / x) = G
            subgoal 1: converges(SUM(n, 0, oo, INT x:[0,1]. x ^ (2 * n) / (2 * n + 1)))
            arg:
                simplify
                apply integral identity
                simplify
            done
            lhs:
                apply series expansion on arctan(x) index n
                rewrite x ^ (2 * n + 1) to x ^ (2 * n) * x
                simplify
                exchange integral and sum
                apply integral identity
                simplify
            rhs:
                expand definition for G
            done
        """
        self.check_actions("interesting", "CatalanConstant01", actions)

    def testCatalanConstant02(self):
        actions = """
            prove (INT x:[0,oo]. log(x + 1) / (x ^ 2 + 1)) = pi / 4 * log(2) + G
            subgoal 1: (INT x:[1,oo]. log(x) / x ^ k) = 1 / (k - 1) ^ 2 for k > 1
            lhs:
                improper integral to limit creating t
                integrate by parts with u = log(x), v = x ^ (1 - k) / (1 - k)
                simplify
                apply integral identity
                simplify
            rhs:
                rewrite 1 / (k - 1) ^ 2 to 1 / (-k + 1) ^ 2
            done
            subgoal 2: converges(SUM(n, 0, oo, INT x:[1,oo]. x ^ (-(2 * n) - 2) * log(x)))
            arg:
                rewrite x ^ (-(2 * n) - 2) * log(x) to log(x) / x ^ (2 * n + 2)
                apply 1 on INT x:[1,oo]. log(x) / x ^ (2 * n + 2)
                simplify
            done
            subgoal 3: (INT x:[1,oo]. log(x) / (x ^ 2 + 1)) = G
            lhs:
                rewrite log(x) / (x ^ 2 + 1) to log(x) * x ^ (-2) * (1 + 1 / x ^ 2) ^ (-1)
                apply series expansion on (1 + 1 / x ^ 2) ^ (-1) index n
                rewrite log(x) * x ^ (-2) * SUM(n, 0, oo, (-1) ^ n * (1 / x ^ 2) ^ n) to SUM(n, 0, oo, (-1) ^ n * (1 / x ^ 2) ^ n * log(x) * x ^ (-2))
                exchange integral and sum
                simplify
                rewrite x ^ (-(2 * n) - 2) * log(x) to log(x) / x ^ (2 * n + 2)
                apply 1 on INT x:[1,oo]. log(x) / x ^ (2 * n + 2)
                simplify
            rhs:
                expand definition for G
            done
            lhs:
                split region at 1
                apply integral identity
                rewrite x + 1 to x * (1 + 1 / x)
                rewrite log(x * (1 + 1 / x)) to log(x) + log(1 + 1 / x)
                rewrite (log(x) + log(1 + 1 / x)) / (x ^ 2 + 1) to log(x) / (x ^ 2 + 1) + log(1 + 1 / x) / (x ^ 2 + 1)
                simplify
                apply 3 on INT x:[1,oo]. log(x) / (x ^ 2 + 1)
                substitute u for 1 / x
                rewrite u ^ 2 * (1 / u ^ 2 + 1) to u ^ 2 + 1
                apply integral identity
                simplify
                rewrite pi * log(2) / 4 to pi / 4 * log(2)
            done
        """
        self.check_actions("interesting", "CatalanConstant02", actions)

    def testCatalanConstant03(self):
        actions = """
            prove (INT x:[0,pi]. x * sin(x) / (a + b * cos(x) ^ 2)) = pi / sqrt(a * b) * arctan(sqrt(b / a)) for a > 0, b > 0
            define I(a,b) = (INT x:[0,pi]. x * sin(x) / (a + b * cos(x) ^ 2)) for a > 0, b > 0
            subgoal 1: I(a,b) = (INT x:[0,pi]. (pi - x) * sin(x) / (a + b * cos(x) ^ 2)) for a > 0, b > 0
            lhs:
                expand definition for I
                substitute x for pi - x
                rewrite sin(x) * (-x + pi) / (b * cos(x) ^ 2 + a) to (pi - x) * sin(x) / (a + b * cos(x) ^ 2)
            done
            lhs:
                fold definition for I (all)
                rewrite I(a,b) to 1/2 * (I(a,b) + I(a,b))
                expand definition for I (at 1)
                apply 1 on I(a,b)
                rewrite (INT x:[0,pi]. x * sin(x) / (b * cos(x) ^ 2 + a)) + (INT x:[0,pi]. (pi - x) * sin(x) / (a + b * cos(x) ^ 2)) to INT x:[0,pi]. x * sin(x) / (a + b * cos(x) ^ 2) + (pi - x) * sin(x) / (a + b * cos(x) ^ 2)
                rewrite x * sin(x) / (a + b * cos(x) ^ 2) + (pi - x) * sin(x) / (a + b * cos(x) ^ 2) to pi * sin(x) / (a + b * cos(x) ^ 2)
                substitute u for cos(x)
                substitute x for sqrt(b / a) * u
                apply integral identity
                simplify
                rewrite arctan(-(sqrt(b) / sqrt(a))) to -arctan(sqrt(b) / sqrt(a))
                simplify
                rewrite sqrt(a) * sqrt(b) to sqrt(a * b)
                rewrite sqrt(b) / sqrt(a) to sqrt(b / a)
            done
        """
        self.check_actions("interesting", "CatalanConstant03", actions)

    def testLogFunction01(self):
        actions = """
            prove (INT x:[0,1]. log(1 + x) / x) = pi ^ 2 / 12
            subgoal 1: converges(SUM(n, 0, oo, INT x:[0,1]. x ^ n / (n + 1)))
            arg:
                simplify
                apply integral identity
                simplify
            done
            lhs:
                apply series expansion on log(1 + x) index n
                rewrite SUM(n, 0, oo, (-1) ^ n * x ^ (n + 1) / (n + 1)) / x to SUM(n, 0, oo, (-1) ^ n * x ^ (n + 1) / (n + 1) * (1 / x))
                exchange integral and sum
                simplify
                apply integral identity
                simplify
                apply series evaluation
            done
        """
        self.check_actions("interesting", "LogFunction01", actions)

    def testLogFunction02(self):#?
        # Inside interesting integrals, Section 5.2, example #2 (5.2.4)
        # 从subgoal 3开始原本的intsumexchange代码无法实现交换
        # SUM(k, 0, oo, (-1) ^ k / (k + 1) ^ 2)无法simplify

        actions = """
            prove (INT x:[0, pi/2]. cos(x)/sin(x) * log(1/cos(x))) = pi^2/24
            subgoal 1: (-log(1-x) - log(1+x)) = -SUM(k,0,oo,(-1)^k*(-x)^(k+1) / (k+1))-SUM(k,0,oo,(-1)^k*x^(k+1)/(k+1)) for abs(x) < 1
            lhs:
                apply series expansion on log(1-x) index k
                apply series expansion on log(1+x) index k
            done
            subgoal 2:x / (-(x ^ 2) + 1) = 1/2 * SUM(k, 0, oo, x ^ k) - 1/2 * SUM(k, 0, oo, x ^ k * (-1) ^ k) for abs(x) < 1
            from 1:
                differentiate both sides at x
                simplify
                rewrite 1 / (-x + 1) - 1 / (x + 1) to 2 * (x / (1-x^2))
                solve equation for x / (1-x^2)
                rewrite (-1) ^ k * (-x) ^ k to x ^ k
                expand polynomial
                rewrite (SUM(k, 0, oo, x ^ k) - SUM(k, 0, oo, x ^ k * (-1) ^ k)) / 2 to 1/2 * SUM(k, 0, oo, x ^ k) - 1/2 * SUM(k, 0, oo, x ^ k * (-1) ^ k)
            done
            subgoal 3:(INT y:[0,1]. (SUM(k, 0, oo, log(y) * y ^ k * (-1) ^ k))) = -SUM(k, 0, oo, (-1) ^ k / (k + 1) ^ 2)
            lhs:
                exchange integral and sum
                apply integral identity
                integrate by parts with u=log(y),v=y^(k+1)/(k+1)
                simplify
                apply integral identity
                simplify
                apply series evaluation
            done
            subgoal 4:(INT y:[0,1]. SUM(k, 0, oo, log(y) * y ^ k)) = -SUM(k, 0, oo, 1 / (k + 1) ^ 2)
            lhs:
                exchange integral and sum
                apply integral identity
                integrate by parts with u=log(y),v=y^(k+1)/(k+1)
                simplify
                apply integral identity
                simplify
                apply series evaluation
            done
            lhs:
                substitute t for cos(x)
                simplify
                substitute y for t
                rewrite y * log(y) / (-(y ^ 2) + 1) to log(y) * (y / (-(y ^ 2) + 1))
                apply 2 on y / (-(y ^ 2) + 1)
                rewrite log(y) * (1/2 * SUM(k, 0, oo, y ^ k) - 1/2 * SUM(k, 0, oo, y ^ k * (-1) ^ k)) to 1/2 * log(y) * SUM(k, 0, oo, y ^ k) - 1/2 * log(y) * SUM(k, 0, oo, y ^ k * (-1) ^ k)
                expand polynomial 
                simplify
                rewrite (log(y) * SUM(k, 0, oo, y ^ k)) to SUM(k, 0, oo, log(y) * y ^ k)
                rewrite log(y) * SUM(k, 0, oo, y ^ k * (-1) ^ k) to SUM(k, 0, oo, log(y) * y ^ k * (-1) ^ k)
                apply 3 on (INT y:[0,1]. (SUM(k, 0, oo, log(y) * y ^ k * (-1) ^ k)))
                apply 4 on (INT y:[0,1]. SUM(k, 0, oo, log(y) * y ^ k))
                apply series evaluation
                simplify
            """
        self.check_actions("interesting", "LogFunction02", actions)

    def testLogFunction03(self):#? new_done
        # 报错StateException: Use done when goal is not finished
        # Inside interesting integrals, Section 5.2, example #3 (5.2.2)
        actions = """
            prove (INT x:[0, 1]. log(1 - x) / x) = -(pi ^ 2 / 6)
            subgoal 1:converges(-SUM(n,1,oo,INT x:[0,1]. x^n/(x*n)))
            arg:
                simplify
                apply integral identity
                simplify
            done
            subgoal 2:(INT x:[0,1]. -(x ^ n / (n + 1))) = -(1/(n+1)^2) for n>=0
            lhs:
                apply integral identity
                simplify
            done
            subgoal 3:SUM(n, 0, oo, 1 / (n + 1) ^ 2) = pi ^ 2 / 6
            lhs:
                apply series evaluation
            done
            lhs:
                apply series expansion on log(1-x) index n
                rewrite SUM(n, 0, oo, (-1) ^ n * (-x) ^ (n + 1) / (n + 1)) / x to SUM(n, 0, oo, (-1) ^ n * (-x) ^ (n + 1) / (n + 1))*(1/x)
                rewrite (-x) ^ (n + 1) to (-1) ^ (n + 1) * x ^ (n + 1)
                rewrite SUM(n, 0, oo, (-1) ^ n * ((-1) ^ (n + 1) * x ^ (n + 1)) / (n + 1)) * (1 / x) to SUM(n, 0, oo, (-1) ^ n * ((-1) ^ (n + 1) * x ^ (n + 1)) / (n + 1) * (1 / x))
                rewrite (-1) ^ n * ((-1) ^ (n + 1) * x ^ (n + 1)) / (n + 1) * (1 / x) to (-1) ^ n * ((-1) ^ (n + 1) * x ^ n) / (n + 1)
                rewrite SUM(n, 0, oo, (-1) ^ n * ((-1) ^ (n + 1) * x ^ n) / (n + 1)) to SUM(n, 0, oo, (-1) ^ (2*n+1) * (x ^ n) / (n + 1))
                simplify
                rewrite -(INT x:[0,1]. SUM(n, 0, oo, x ^ n / (n + 1))) to INT x:[0,1]. SUM(n, 0, oo, -(x ^ n / (n + 1)))
                exchange integral and sum
                apply 2 on (INT x:[0,1]. -(x ^ n / (n + 1)))
                simplify
                apply 3 on SUM(n, 0, oo, 1 / (n + 1) ^ 2)
            done
            """
        self.check_actions("interesting", "LogFunction03(不存在)", actions)

    def testBernoulliIntegral(self):
        actions = """
            prove (INT x:[0,1]. x ^ (c * x ^ a)) = SUM(k, 0, oo, (-c) ^ k / (k * a + 1) ^ (k + 1)) for a > 0, c != 0
            subgoal 1: converges(SUM(k, 0, oo, abs(INT x:[0,1]. (c * x ^ a * log(x)) ^ k / factorial(k))))
            arg:
                simplify
                rewrite (c * x ^ a * log(x)) ^ k to (c * x ^ a) ^ k * log(x) ^ k
                rewrite (c * x ^ a) ^ k to c ^ k * x ^ a ^ k
                simplify
                apply integral identity
                simplify
            done
            lhs:
                rewrite x ^ (c * x ^ a) to exp(log(x ^ (c * x ^ a)))
                apply series expansion on exp(log(x ^ (c * x ^ a))) index k
                exchange integral and sum
                rewrite log(x ^ (c * x ^ a)) to c * x ^ a * log(x)
                rewrite (c * x ^ a * log(x)) ^ k to (c * x ^ a) ^ k * log(x) ^ k
                rewrite (c * x ^ a) ^ k to c ^ k * x ^ a ^ k
                simplify
                apply integral identity
                simplify
                rewrite c ^ k * (-1) ^ k to (-c) ^ k
            done
            prove (INT x:[0,1]. x ^ x) = SUM(k, 0, oo, (-1) ^ k * (k + 1) ^ (-k - 1))
            lhs:
                rewrite x ^ x to x ^ (1 * x ^ 1)
                apply integral identity
                simplify
            done
            prove (INT x:[0,1]. x ^ -x) = SUM(k, 0, oo, (k + 1) ^ (-k - 1))
            lhs:
                rewrite x ^ -x to x ^ (-1 * x ^ 1)
                apply integral identity
                simplify
            done
            prove (INT x:[0,1]. x ^ (x ^ 2)) = SUM(k, 0, oo, (-1) ^ k * (2 * k + 1) ^ (-k - 1))
            lhs:
                rewrite x ^ (x ^ 2) to x ^ (1 * x ^ 2)
                apply integral identity
                simplify
            done
            prove (INT x:[0,1]. x ^ sqrt(x)) = SUM(k, 0, oo, (-1) ^ k * (2 / (k + 2)) ^ (k + 1))
            lhs:
                rewrite x ^ sqrt(x) to x ^ (1 * x ^ (1/2))
                apply integral identity
                simplify
                rewrite k / 2 + 1 to (2 / (k + 2)) ^ (-1)
                rewrite (2 / (k + 2)) ^ (-1) ^ (-k - 1) to (2 / (k + 2)) ^ (k + 1)
            done
        """
        self.check_actions("interesting", "BernoulliIntegral", actions)

    def testAhmedIntegral(self):
        actions = """
            prove (INT x:[0,1]. arctan(sqrt(2 + x ^ 2)) / ((1 + x ^ 2) * sqrt(2 + x ^ 2))) = 5 * pi ^ 2 / 96
            define I(u) = (INT x:[0,1]. arctan(u * sqrt(2 + x ^ 2)) / ((1 + x ^ 2) * sqrt(2 + x ^ 2))) for u > 0
            subgoal 1: I(1) = (INT x:[0,1]. arctan(sqrt(x ^ 2 + 2)) / ((x ^ 2 + 1) * sqrt(x ^ 2 + 2)))
            lhs:
                expand definition for I
            done
            subgoal 2: (D u. I(u)) = 1 / (1 + u ^ 2) * (pi / 4 - u / sqrt(1 + 2 * u ^ 2) * arctan(u / sqrt(1 + 2 * u ^ 2))) for u > 0
            lhs:
                expand definition for I (all)
                exchange derivative and integral
                simplify
                rewrite 1 / ((x ^ 2 + 1) * (u ^ 2 * (x ^ 2 + 2) + 1)) to 1 / (u ^ 2 + 1) * (1 / (1 + x ^ 2) - u ^ 2 / (1 + 2 * u ^ 2 + u ^ 2 * x ^ 2))
                simplify
                rewrite 1 / (u ^ 2 * x ^ 2 + 2 * u ^ 2 + 1) to u ^ (-2) * (x ^ 2 + (2 * u ^ 2 + 1) / u ^ 2) ^ (-1)
                simplify
                substitute y * sqrt(u ^ (-2) * (2 * u ^ 2 + 1)) for x
                simplify
                rewrite 1 / (y ^ 2 * (2 * u ^ 2 + 1) / u ^ 2 + (2 * u ^ 2 + 1) / u ^ 2) to 1 / (y ^ 2 + 1) * (u ^ 2 / (2 * u ^ 2 + 1))
                apply integral identity
                simplify
            done
            subgoal 3: (INT u:[1,oo]. D u. I(u)) = pi ^ 2 / 12 - I(1)
            lhs:
                simplify
                expand definition for I (at 1)
                simplify
                integrate by parts with u = 1, v = arctan(x / sqrt(2 + x ^ 2)) / 2
                simplify
            done
            subgoal 4: (INT u:[1,oo]. D u. I(u)) = -(pi ^ 2 / 48) + I(1)
            lhs:
                apply 2 on D u. I(u)
                expand polynomial
                simplify
                substitute 1 / x for u
                simplify
                rewrite x ^ 3 * (1 / x ^ 2 + 1) * sqrt(2 / x ^ 2 + 1) to sqrt((1 + x ^ 2) ^ 2 * (2 + x ^ 2))
                rewrite x * sqrt(2 / x ^ 2 + 1) to sqrt(x ^ 2 + 2)
                simplify
                rewrite 1 / sqrt(x ^ 2 + 2) to sqrt(x ^ 2 + 2) ^ (-1)
                rewrite arctan(sqrt(x ^ 2 + 2) ^ (-1)) to pi / 2 - arctan(sqrt(x ^ 2 + 2))
                expand polynomial
                simplify
                rewrite arctan(sqrt(x ^ 2 + 2)) / (x ^ 2 * sqrt(x ^ 2 + 2) + sqrt(x ^ 2 + 2)) to arctan(sqrt(x ^ 2 + 2)) / ((x ^ 2 + 1) * sqrt(x ^ 2 + 2))
                apply 1 on INT x:[0,1]. arctan(sqrt(x ^ 2 + 2)) / ((x ^ 2 + 1) * sqrt(x ^ 2 + 2))
                integrate by parts with u = 1, v = arctan(x / sqrt(2 + x ^ 2))
                apply integral identity
                simplify
            done
            from 3:
                apply 4 on INT u:[1,oo]. D u. I(u)
                solve equation for I(1)
                expand definition for I (all)
            done
        """
        self.check_actions("interesting", "AhmedIntegral", actions)

    def testEulerConstant01(self):
        actions = """
            define EulerConstant = -(INT x:[0, oo]. exp(-x) * log(x))
            prove (INT x:[0,1]. (-exp(-x) + 1) / x) - (INT x:[1,oo]. exp(-x) / x) = EulerConstant
            lhs:
                integrate by parts with u = exp(-x), v = log(x) (at 2)
                integrate by parts with u = 1 - exp(-x), v = log(x)
                simplify
            rhs:
                expand definition for EulerConstant
                split region at 1
                simplify
            done
        """
        self.check_actions("interesting", "EulerConstant01", actions)

    def testChapter3Practice01(self):
        actions = """
            prove (INT x:[0,oo]. log(1 + a ^ 2 * x ^ 2) / (b ^ 2 + x ^ 2)) = pi / b * log(1 + a * b) for a > 0, b > 0
            define I(a,b) = (INT x:[0,oo]. log(1 + a ^ 2 * x ^ 2) / (b ^ 2 + x ^ 2)) for a >= 0, b > 0
            subgoal 1: (D a. I(a,b)) = pi / (1 + a * b) for a > 0, b > 0
            lhs:
                expand definition for I (all)
                exchange derivative and integral
                simplify
                rewrite x ^ 2 / ((b ^ 2 + x ^ 2) * (a ^ 2 * x ^ 2 + 1)) to 1 / (1 - a ^ 2 * b ^ 2) * (1 / (1 + a ^ 2 * x ^ 2) - b ^ 2 / (b ^ 2 + x ^ 2))
                simplify
                apply integral identity
                simplify
                rewrite to pi / (1 + a * b)
            done
            subgoal 2: I(a,b) = pi / b * log(1 + a * b) + SKOLEM_FUNC(C(b)) for a > 0, b > 0
            from 1:
                integrate both sides
                substitute u for 1 + a * b
                simplify
                apply integral identity
                replace substitution
                rewrite abs(1 + a * b) to 1 + a * b
            done
            subgoal 3: I(0,b) = 0 for b > 0
            lhs:
                expand definition for I
                simplify
            done
            subgoal 4: SKOLEM_FUNC(C(b)) = 0 for b > 0
            from 2:
                apply limit a -> 0 both sides
                simplify
                apply 3 on I(0,b)
                solve equation for SKOLEM_FUNC(C(b))
            done
            from 2:
                apply 4 on SKOLEM_FUNC(C(b))
                expand definition for I (all)
                simplify
            done
        """
        self.check_actions("interesting", "Chapter3Practice01", actions)

    def testChapter3Practice02(self):
        actions = """
            prove (INT x:[-oo,oo]. cos(a * x) / (b ^ 2 - x ^ 2)) = pi * sin(a * b) / b for a > 0, b > 0, b != x
            lhs:
                rewrite b ^ 2 - x ^ 2 to (b + x) * (b - x)
                rewrite cos(a * x) / ((b + x) * (b - x)) to 1 / (2 * b) * (cos(a * x) / (b + x) + cos(a * x) / (b - x))
                simplify
                substitute u for b + x
                substitute u for b - x (at 2)
                rewrite a * (-b + u) to -(a * (b - u))
                rewrite cos(-(a * (b - u))) to cos(a * (b - u))
                simplify
                rewrite cos(a * (b - u)) to cos(a * b - a * u)
                rewrite cos(a * b - a * u) to cos(a * b) * cos(a * u) + sin(a * b) * sin(a * u)
                rewrite (cos(a * b) * cos(a * u) + sin(a * b) * sin(a * u)) / u to cos(a * b) * cos(a * u) / u + sin(a * b) * sin(a * u) / u
                simplify
                rewrite INT u:[-oo,oo]. cos(a * u) / u to 0
                simplify
                split region at 0
                substitute u for -u
                simplify
                apply integral identity
                simplify
            done
        """
        self.check_actions("interesting", "Chapter3Practice02", actions)

    def testChapter3Practice03(self):
        actions = """
            prove (INT x:[-oo,oo]. cos(a * x) / (b ^ 4 - x ^ 4)) = pi * (exp(-(a * b)) + sin(a * b)) / (2 * b ^ 3) for a > 0, b > 0, b != x
            lhs:
                rewrite b ^ 4 - x ^ 4 to (b ^ 2 + x ^ 2) * (b ^ 2 - x ^ 2)
                rewrite cos(a * x) / ((b ^ 2 + x ^ 2) * (b ^ 2 - x ^ 2)) to 1 / (2 * b ^ 2) * (cos(a * x) / (b ^ 2 + x ^ 2) + cos(a * x) / (b ^ 2 - x ^ 2))
                simplify
                split region at 0
                substitute x for -x
                simplify
                rewrite b ^ 2 + x ^ 2 to x ^ 2 + b ^ 2
                apply integral identity
                simplify
                rewrite to pi * (exp(-(a * b)) + sin(a * b)) / (2 * b ^ 3)
            done
        """
        self.check_actions("interesting", "Chapter3Practice03", actions)

    def testChapter3Practice04(self):
        actions = """
            prove (INT x:[0,oo]. x * sin(a * x) / (x ^ 2 - b ^ 2)) = pi / 2 * cos(a * b) for a > 0, b > 0, b != x
            subgoal 1: (INT x:[0,oo]. x * sin(a * x) / (x ^ 2 - b ^ 2)) = 1/2 * (INT x:[-oo,oo]. x * sin(a * x) / (x ^ 2 - b ^ 2)) for x != b, a > 0, b > 0
            lhs:
                simplify
            rhs:
                split region at 0
                substitute x for -x
                simplify
            done
            from 1:
                rewrite x ^ 2 - b ^ 2 to (x + b) * (x - b) (at 2)
                rewrite x * sin(a * x) / ((x + b) * (x - b)) to -x * sin(a * x) / ((b - x) * (b + x))
                rewrite -x * sin(a * x) / ((b - x) * (b + x)) to -1 / (2 * b) * (x * sin(a * x) / (b + x) + x * sin(a * x) / (b - x))
                simplify
                substitute u for b + x (at 2)
                substitute u for b - x (at 3)
                rewrite sin(a * (-b + u)) to sin(-(a * (b - u)))
                rewrite sin(-(a * (b - u))) to -sin(a * (b - u))
                rewrite (-b + u) * -sin(a * (b - u)) to (b - u) * sin(a * (b - u))
                simplify
                rewrite INT u:[-oo,oo]. (b - u) * sin(a * (b - u)) / u to INT u:[-oo,oo]. (b - u) / u * sin(a * (b - u))
                rewrite (b - u) / u * sin(a * (b - u)) to (b / u - 1) * sin(a * b - a * u)
                rewrite (b / u - 1) * sin(a * b - a * u) to b / u * sin(a * b - a * u) - sin(a * b - a * u)
                simplify
                substitute s for a * b - a * u (at 3)
                split region at 0 (at 3)
                substitute s for -s (at 3)
                simplify
                rewrite sin(a * b - a * u) to sin(a * b) * cos(a * u) - cos(a * b) * sin(a * u)
                rewrite (sin(a * b) * cos(a * u) - cos(a * b) * sin(a * u)) / u to sin(a * b) * cos(a * u) / u - cos(a * b) * sin(a * u) / u
                simplify
                split region at 0 (at 3)
                substitute u for -u (at 3)
                split region at 0 (at 2)
                substitute u for -u (at 2)
                simplify
                apply integral identity
            done
        """
        self.check_actions("interesting", "Chapter3Practice04", actions)

    def testChapter3Practice05(self):#?
        # Inside interesting integrals, Section 3.10, C3.5
        # 输出<integral.action.CaseAnalysisState object at 0x0000021C0ABE5D50>
        actions = """
            prove (INT x:[0, oo]. cos(a * x) * sin(b * x) / x) = pi/2 for a > 0,b > 0
            subgoal 1:(INT x:[0, oo]. cos(a * x) * sin(b * x) / x) = 1/2 * (INT x:[0, oo]. sin((b + a) * x) / x) + 1/2 * (INT x:[0, oo]. sin((b - a) * x) / x) for a > 0,b > 0
            lhs:

                rewrite cos(a * x) * sin(b * x) to 1/2 * (sin(b * x + a * x) - sin(a * x - b * x))
                rewrite 1/2 * (sin(b * x + a * x) - sin(a * x - b * x)) / x to 1/2 * sin((b + a) * x) / x - 1/2 * sin(-((b - a) * x)) / x
                rewrite sin(-((b - a) * x)) to -sin((b - a) * x)
                simplify
            rhs:
                simplify
            done
            define g(y,a) = INT x:[0,oo]. exp(-x * y) * sin(a * x) / x for y>=0
            subgoal 2: (D y. g(y, a)) = - a / (a ^ 2 + y ^ 2) for y>=0,a!=0
            lhs:
                expand definition for g(all)
                exchange derivative and integral
                simplify
                apply integral identity
            rhs:
                simplify
            done
            subgoal 3:g(y, a) = -arctan(y / a) + SKOLEM_FUNC(C(a)) for y>=0,a>0
            from 2:
                integrate both sides
                apply integral identity
                simplify
            done

            subgoal 4:(LIM {y -> oo}. g(y, a)) = 0 for y>=0
            lhs:
                expand definition for g(all)
                simplify
            done
            subgoal 5:SKOLEM_FUNC(C(a)) = pi / 2 for a>0
            from 3:
                apply limit y -> oo both sides
                apply 4 on LIM {y -> oo}. g(y,a)
                simplify
                solve equation for SKOLEM_FUNC(C(a))
            done
            subgoal 6:g(0,a) = pi / 2 for a>0
            from 3:
                apply 5 on SKOLEM_FUNC(C(a)) 
                simplify
            done
            define g(y,a) = INT x:[0,oo]. exp(-x * y) * sin(a * x) / x for y=0
            subgoal 7:(INT x:[0,oo]. sin(a*x) / x) = g(0,a) for a!=0
            rhs:
                expand definition for g(all)
            done
            subgoal 8:(INT x:[0,oo]. sin(a*x) / x) = pi/2 for a>0
            lhs:
                apply 7 on (INT x:[0,oo]. sin(a*x) / x)
                apply 6 on g(0,a)
            done
            subgoal 811:g(0,0) = 0
            lhs:
                expand definition for g(all)
                simplify
            done
            subgoal 81:(INT x:[0,oo]. sin(a*x) / x) = 0 for a=0
            lhs:
                apply 7 on (INT x:[0,oo]. sin(a*x) / x)
                simplify
                apply 811 on g(0,0)
            done
            subgoal 82:(INT x:[0,oo]. sin(a*x) / x) = -pi/2 for a<0
            lhs:
                simplify
            done
            subgoal 9:(INT x:[0, oo]. cos(a * x) * sin(b * x) / x) = pi/2 for b-a>0,a > 0,b > 0
            lhs:
                apply 1 on (INT x:[0, oo]. cos(a * x) * sin(b * x) / x)
                substitute u for (b+a)*x
                rewrite (INT u:[0,oo]. sin(u) / u) to (INT u:[0,oo]. sin(1*u) / u)
                apply 8 on (INT u:[0,oo]. sin(1*u) / u)
                apply 8 on (INT x:[0,oo]. sin((b - a) * x) / x)
                simplify
            done
            subgoal 10:(INT x:[0, oo]. cos(a * x) * sin(b * x) / x) = pi/2 for b-a=0,a > 0,b > 0
            lhs:
                apply 1 on (INT x:[0, oo]. cos(a * x) * sin(b * x) / x)
                substitute u for (b+a)*x
                rewrite (INT u:[0,oo]. sin(u) / u) to (INT u:[0,oo]. sin(1*u) / u)
                apply 8 on (INT u:[0,oo]. sin(1*u) / u)
                apply 81 on (INT x:[0,oo]. sin((b - a) * x) / x)
                simplify
            done
            subgoal 11:(INT x:[0, oo]. cos(a * x) * sin(b * x) / x) = pi/2 for b-a<0,a > 0,b > 0
            lhs:
                apply 1 on (INT x:[0, oo]. cos(a * x) * sin(b * x) / x)
                substitute u for (b+a)*x
                rewrite (INT u:[0,oo]. sin(u) / u) to (INT u:[0,oo]. sin(1*u) / u)
                apply 8 on (INT u:[0,oo]. sin(1*u) / u)
                apply 82 on (INT x:[0,oo]. sin((b - a) * x) / x)
                simplify
            done
            case analysis on b-a
            case negative :
            lhs:
                apply 11 on (INT x:[0,oo]. cos(a * x) * sin(b * x) / x)
            done
            case zero :
            lhs:
                apply 10 on (INT x:[0,oo]. cos(a * x) * sin(b * x) / x)
            done
            case positive :
            lhs:
                apply 9 on (INT x:[0,oo]. cos(a * x) * sin(b * x) / x)
            done

            """
        self.check_actions("interesting", "Chapter3Practice05", actions)

    def testChapter3Practice06(self):
        actions = """
            prove (INT x:[-1,1]. ((1 + x) / (1 - x)) ^ (1/2)) = pi
            lhs:
                substitute cos(2 * u) for x
                rewrite cos(2 * u) to 2 * cos(u) ^ 2 - 1 (at 1)
                rewrite cos(2 * u) to 1 - 2 * sin(u) ^ 2
                simplify
                rewrite sin(2 * u) to 2 * sin(u) * cos(u)
                simplify
                apply integral identity
                simplify
            done
        """
        self.check_actions("interesting", "Chapter3Practice06", actions)

    def testChapter3Practice07a(self):
        actions = """
            prove (INT x:[-oo,oo]. x * exp(-(x ^ 2) - x)) = -1/2 * sqrt(pi * sqrt(exp(1)))
            define I(a,b) = (INT x:[-oo,oo]. exp(-a * x ^ 2 + b * x)) for a > 0
            subgoal 1: I(a,b) = exp(b ^ 2 / (4 * a)) * sqrt(pi / a) for a > 0
            lhs:
                expand definition for I
                rewrite -(a * x ^ 2) + b * x to b ^ 2 / (4 * a) - a * (x - b / (2 * a)) ^ 2
                rewrite exp(b ^ 2 / (4 * a) - a * (x - b / (2 * a)) ^ 2) to exp(b ^ 2 / (4 * a)) * exp(-a * (x - b / (2 * a)) ^ 2)
                simplify
                substitute y for x - b / (2 * a)
                apply integral identity
                simplify
            done
            subgoal 2: (D b. I(a,b)) = b / (2 * a) * exp(b ^ 2 / (4 * a)) * sqrt(pi / a) for a > 0
            lhs:
                apply 1 on I(a,b)
                simplify
            done
            subgoal 3: (INT x:[-oo,oo]. x * exp(-(a * x ^ 2) + b * x)) = b / (2 * a) * exp(b ^ 2 / (4 * a)) * sqrt(pi / a) for a > 0
            from 2:
                expand definition for I (all)
                simplify
            done
            lhs:
                rewrite x * exp(-(x ^ 2) - x) to x * exp(-(1 * x ^ 2) + -1 * x)
                apply 3 on INT x:[-oo,oo]. x * exp(-(1 * x ^ 2) + -1 * x)
                rewrite to -1/2 * sqrt(pi * sqrt(exp(1)))
            done
        """
        self.check_actions("interesting", "Chapter3Practice07a", actions)

    def testChapter3Practice07b(self):
        actions = """
            prove (INT x:[-oo,oo]. x ^ 2 * exp(-(x ^ 2) - x)) = 3/4 * sqrt(pi * sqrt(exp(1)))
            define I(a,b) = (INT x:[-oo,oo]. exp(-a * x ^ 2 + b * x)) for a > 0
            subgoal 1: I(a,b) = exp(b ^ 2 / (4 * a)) * sqrt(pi / a) for a > 0
            lhs:
                expand definition for I
                rewrite -(a * x ^ 2) + b * x to b ^ 2 / (4 * a) - a * (x - b / (2 * a)) ^ 2
                rewrite exp(b ^ 2 / (4 * a) - a * (x - b / (2 * a)) ^ 2) to exp(b ^ 2 / (4 * a)) * exp(-a * (x - b / (2 * a)) ^ 2)
                simplify
                substitute y for x - b / (2 * a)
                apply integral identity
                simplify
            done
            subgoal 2: (D a. I(a,b)) = -(b ^ 2 / (4 * a ^ 2)) * exp(b ^ 2 / (4 * a)) * sqrt(pi / a) - 1 / (2 * a) * exp(b ^ 2 / (4 * a)) * sqrt(pi / a) for a > 0
            lhs:
                apply 1 on I(a,b)
                simplify
            rhs:
                simplify
            done
            subgoal 3: (INT x:[-oo,oo]. x ^ 2 * exp(-(a * x ^ 2) + b * x)) = b ^ 2 / (4 * a ^ 2) * exp(b ^ 2 / (4 * a)) * sqrt(pi / a) + 1 / (2 * a) * exp(b ^ 2 / (4 * a)) * sqrt(pi / a) for a > 0
            from 2:
                expand definition for I (all)
                simplify
                solve equation for INT x:[-oo,oo]. x ^ 2 * exp(-(a * x ^ 2) + b * x)
            done
            lhs:
                rewrite x ^ 2 * exp(-(x ^ 2) - x) to x ^ 2 * exp(-(1 * x ^ 2) + -1 * x)
                apply 3 on INT x:[-oo,oo]. x ^ 2 * exp(-(1 * x ^ 2) + -1 * x)
                rewrite to 3/4 * sqrt(pi * sqrt(exp(1)))
            done
        """
        self.check_actions("interesting", "Chapter3Practice07b", actions)

    def testChapter3Practice08(self):
        actions = """
            prove (INT x:[0,oo]. sin(m * x) / (x * (a ^ 2 + x ^ 2) ^ 2)) = pi / (2 * a ^ 4) * (1 - (2 + m * a) / 2 * exp(-a * m)) for a > 0, m > 0
            subgoal 1: (INT x:[0,oo]. sin(m * x) / (x * (a ^ 2 + x ^ 2))) = pi * (1 - exp(-a * m)) / (2 * a ^ 2) for a > 0, m > 0
            lhs:
                apply integral identity
            done
            from 1:
                differentiate both sides at a
                exchange derivative and integral (all)
                simplify
                solve equation for INT x:[0,oo]. sin(m * x) / (x * (a ^ 2 + x ^ 2) ^ 2)
                rewrite -((2 * a ^ 2 * m * pi * exp(-(a * m)) - 4 * a * pi * (-exp(-(a * m)) + 1)) / (8 * a ^ 5)) to pi / (2 * a ^ 4) * (1 - (2 + m * a) / 2 * exp(-a * m))
            done
        """
        self.check_actions("interesting", "Chapter3Practice08", actions)

    def testChapter3Practice09(self):
        actions = """
            prove (INT x:[0,1]. x / (a * x + b * (1 - x)) ^ 3) = 1 / (2 * a ^ 2 * b) for a > 0, b > 0, a > b
            subgoal 1: (INT x:[0,1]. 1 / (a * x + b * (1 - x)) ^ 2) = 1 / (a * b) for a > 0, b > 0, a > b
            lhs:
                substitute u for (a - b) * x + b
                rewrite 1 / ((a - b) * (b * (-((-b + u) / (a - b)) + 1) + a * (-b + u) / (a - b)) ^ 2) to 1 / (u ^ 2 * (a - b))
                apply integral identity
                simplify
                rewrite 1 / (a - b) * (-(1 / a) + 1 / b) to 1 / (a * b)
            done
            from 1:
                differentiate both sides at a
                exchange derivative and integral (all)
                simplify
                rewrite (b * (-x + 1) + a * x) ^ 3 to (a * x + b * (1 - x)) ^ 3
                solve equation for INT x:[0,1]. x / (a * x + b * (1 - x)) ^ 3
            done
        """
        self.check_actions("interesting", "Chapter3Practice09", actions)

    def testChapter1Practice0104(self):
        actions = """
            prove (INT x:[0,pi / 3]. 1 / cos(x)) = log(2 + sqrt(3))
            lhs:
                rewrite 1 / cos(x) to cos(x) / cos(x) ^ 2
                rewrite cos(x) ^ 2 to 1 - sin(x) ^ 2
                substitute u for sin(x)
                rewrite 1 / (-(u ^ 2) + 1) to 1/2 * (1 / (1 - u) + 1 / (1 + u))
                simplify
                apply integral identity
                simplify
                rewrite -(1/2 * log(-(sqrt(3) / 2) + 1)) + 1/2 * log(sqrt(3) / 2 + 1) to 1/2 * (log(sqrt(3) / 2 + 1) - log(-(sqrt(3) / 2) + 1))
                rewrite log(sqrt(3) / 2 + 1) - log(-(sqrt(3) / 2) + 1) to log((sqrt(3) / 2 + 1) / (-(sqrt(3) / 2) + 1))
                rewrite (sqrt(3) / 2 + 1) / (-(sqrt(3) / 2) + 1) to (2 + sqrt(3)) ^ 2
                simplify
                rewrite sqrt(3) + 2 to 2 + sqrt(3)
            done
        """
        self.check_actions("interesting", "Chapter1Practice0104", actions)

    def testChapter2Practice01(self):
        actions = """
            prove (INT x:[0,4]. log(x) / sqrt(4 * x - x ^ 2)) = 0
            subgoal 1: (INT y:[0,1]. 1 / (sqrt(y) * sqrt(1 - y))) = pi
            lhs:
                substitute sin(x) ^ 2 for y
                rewrite sin(x) ^ 2 to 1 - cos(x) ^ 2 (at 2)
                simplify
                apply integral identity
                simplify
            done
            subgoal 2: (INT y:[0,1]. log(y) / (sqrt(y) * sqrt(1 - y))) = -(2 * pi * log(2))
            lhs:
                substitute sin(x) ^ 2 for y
                rewrite log(sin(x) ^ 2) to 2 * log(sin(x))
                rewrite sin(x) ^ 2 to 1 - cos(x) ^ 2 (at 2)
                simplify
                rewrite sin(x) to 1 * sin(x)
                apply integral identity
                simplify
            done
            subgoal 3: 4 * x - x ^ 2 >= 0 for x > 0, x < 4
            lhs:
                rewrite 4 * x - x ^ 2 to x * (4 - x)
            done
            subgoal 4: sqrt(4 * x - x ^ 2) != 0 for x > 0, x < 4
            lhs:
                rewrite 4 * x - x ^ 2 to x * (4 - x)
            done
            lhs:
                substitute y for x / 4
                rewrite log(4 * y) to log(4) + log(y)
                rewrite sqrt(-(16 * y ^ 2) + 16 * y) to 4 * sqrt(-(y ^ 2) + y)
                rewrite sqrt(-(y ^ 2) + y) to sqrt(y) * sqrt(1 - y)
                expand polynomial
                simplify
                rewrite -y + 1 to 1 - y
                rewrite -y + 1 to 1 - y
                apply 1 on INT y:[0,1]. 1 / (sqrt(y) * sqrt(1 - y))
                apply 2 on INT y:[0,1]. log(y) / (sqrt(y) * sqrt(1 - y))
                simplify
            done
        """
        self.check_actions("interesting", "Chapter2Practice01", actions)

    def testChapter2Practice02(self):#?+
        # Inside interesting integrals, C2.2
        actions = """
            prove (INT x:[0,1]. (x - 2) / (x ^ 2 - x + 1)) = -pi/sqrt(3)
            subgoal 1:(INT u:[-1/2,0]. u / (u ^ 2 + 3/4)) = -(INT u:[0,1/2]. u / (u ^ 2 + 3/4))
            lhs:
                substitute t for -u
                simplify
            done
            subgoal 2:(INT u:[-1/2,1/2]. u/(u^2+3/4)) = 0
            lhs:
                split region at 0
                apply 1 on INT u:[-1/2,0]. u / (u ^ 2 + 3/4)
                simplify
            done
            subgoal 3:3/2*(INT u:[-1/2,1/2]. 1/(u^2+3/4)) = pi/sqrt(3)
            lhs:
                simplify
                rewrite 1 / (u ^ 2 + 3/4) to (4/3)/((4/3)*u^2+(4/3)*(3/4))
                simplify
                substitute t for (2*u)/sqrt(3)
                rewrite sqrt(3) / (2 * t ^ 2 + 2) to (sqrt(3))/2*(1/(t^2+1))
                simplify
                apply integral identity
                simplify
                rewrite to pi/sqrt(3)
            done
            lhs:
                substitute u for x-1/2
                rewrite (u + 1/2) ^ 2 - u + 1/2 to u^2+3/4
                expand polynomial
                simplify
                apply 2 on INT u:[-1/2,1/2]. u / (u ^ 2 + 3/4)
                apply 3 on 3/2 * (INT u:[-1/2,1/2]. 1 / (u ^ 2 + 3/4))
                rewrite to -pi/sqrt(3)
            done
            """
        self.check_actions("interesting", "chapter2_practice02", actions)

    def testChapter2Practice03(self):
        actions = """
            prove (INT x:[0,oo]. 1 / (x ^ 4 + 1) ^ (m + 1)) = (4 * m - 1) / (4 * m) * (INT x:[0,oo]. 1 / (x ^ 4 + 1) ^ m) for m >= 1, isInt(m)
            subgoal 1: (INT x:[0,oo]. (x ^ 4 + 1) ^ -m) = 4 * m * ((INT x:[0,oo]. 1 / (x ^ 4 + 1) ^ m) - (INT x:[0,oo]. 1 / (x ^ 4 + 1) ^ (m + 1))) for m >= 1, isInt(m)
            lhs:
                integrate by parts with u = 1 / (x ^ 4 + 1) ^ m, v = x
                simplify
                rewrite x ^ 4 * (x ^ 4 + 1) ^ (-m - 1) to (x ^ 4 + 1) / (x ^ 4 + 1) ^ (m + 1) - 1 / (x ^ 4 + 1) ^ (m + 1)
                rewrite INT x:[0,oo]. (x ^ 4 + 1) / (x ^ 4 + 1) ^ (m + 1) - 1 / (x ^ 4 + 1) ^ (m + 1) to (INT x:[0,oo]. 1 / (x ^ 4 + 1) ^ m) - (INT x:[0,oo]. 1 / (x ^ 4 + 1) ^ (m + 1))
            done
            from 1:
                solve equation for INT x:[0,oo]. 1 / (x ^ 4 + 1) ^ (m + 1)
                rewrite -(1 / (4 * m) * (INT x:[0,oo]. (x ^ 4 + 1) ^ -m)) + (INT x:[0,oo]. (x ^ 4 + 1) ^ -m) to (4 * m - 1) / (4 * m) * (INT x:[0,oo]. 1 / (x ^ 4 + 1) ^ m)
            done
        """
        self.check_actions("interesting", "Chapter2Practice03", actions)

    def testChapter2Practice05(self):
        actions = """
            prove (INT x:[0,oo]. log(x + 1) / x ^ (3/2)) = 2 * pi
            lhs:
                integrate by parts with u = log(1 + x), v = -2 / sqrt(x)
                simplify
                substitute t for sqrt(x)
                simplify
                apply integral identity
                simplify
            done
        """
        self.check_actions("interesting", "Chapter2Practice05", actions)

    def testPostgraduateIndefinitePart1SectionA(self):
        with open('theories/postgradIndef1a.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("base", "PostgraduateIndefinitePart1SectionA", actions)

    def testPostgraduateIndefinitePart1SectionB(self):
        with open('theories/postgradIndef1b.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("base", "PostgraduateIndefinitePart1SectionB", actions)

    def testPostgraduateIndefinitePart2SectionA(self):
        with open('theories/postgradIndef2a.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("base", "PostgraduateIndefinitePart2SectionA", actions)

    def testPostgraduateIndefinitePart2SectionB(self):
        with open('theories/postgradIndef2b.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("base", "PostgraduateIndefinitePart2SectionB", actions)

    def testPostgraduateIndefinitePart3SectionA(self):
        with open('theories/postgradIndef3a.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("base", "PostgraduateIndefinitePart3SectionA", actions)

    def testPostgraduateIndefinitePart4SectionA(self):
        with open('theories/postgradIndef4a.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("base", "PostgraduateIndefinitePart4SectionA", actions)

    def testPostgraduateIndefinitePart4SectionB(self):
        with open('theories/postgradIndef4b.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("base", "PostgraduateIndefinitePart4SectionB", actions)

    def testPostgraduateIndefinitePart5SectionA(self):
        with open('theories/postgradIndef5a.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("base", "PostgraduateIndefinitePart5SectionA", actions)

    def testPostgraduateIndefinitePart5SectionB(self):
        with open('theories/postgradIndef5b.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("base", "PostgraduateIndefinitePart5SectionB", actions)

    def testPostgraduateIndefinitePart6SectionA(self):
        with open('theories/postgradIndef6a.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("base", "PostgraduateIndefinitePart6SectionA", actions)

    def testPostgraduateIndefinitePart6SectionB(self):
        with open('theories/postgradIndef6b.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("base", "PostgraduateIndefinitePart6SectionB", actions)

    def testPostgraduateDefinitePart1SectionA(self):
        with open('theories/postgradDef1a.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("base", "PostgraduateDefinitePart1SectionA", actions)


if __name__ == "__main__":
    unittest.main()
