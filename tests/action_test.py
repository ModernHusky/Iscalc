"""Unit test for integrals using internal language."""

import unittest
import lark

from integral import compstate
from integral import action
from integral import parser


class ActionTest(unittest.TestCase):
    def check_actions(self, base_file, current_file, actions: str,
                      *, print_lines=False, print_state=False):
        file = compstate.CompFile(base_file, current_file)
        state = action.InitialState(file)
        actions = [s for s in actions.split('\n') if s.strip()]
        for act in actions:
            if print_lines:
                print(act)
            if act.startswith('#') or act.startswith('//'):
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

            prove (INT x. 1 / (a ^ 2 + x ^ 2)) = 1 / a * atan(x / a) + SKOLEM_CONST(C) for a != 0
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

    def testTongji(self):
        actions = """
            calculate INT x:[2,3]. 2 * x + x ^ 2
                apply integral identity
                simplify
            done

            calculate INT x:[0,1]. (3 * x + 1) ^ (-2)
                substitute u for 3 * x + 1
                apply integral identity
                simplify
            done

            calculate INT x:[0,1]. exp(6 * x)
                apply integral identity
                simplify
            done

            calculate INT x:[-1,2]. x * exp(x)
                integrate by parts with u = x, v = exp(x)
                apply integral identity
                simplify
            done

            calculate INT x:[0,pi/4]. sin(x)
                apply integral identity
                simplify
            done

            calculate INT x:[0,1]. 3*x^2 - x + 1
                apply integral identity
                simplify
            done

            calculate INT x:[1,2]. x^2 + 1/x^4
                apply integral identity
                simplify
            done

            calculate INT x:[pi/3, pi]. sin(2*x + pi/3)
                substitute u for 2*x + pi/3
                apply integral identity
                simplify
            done

            calculate INT x:[4, 9]. x ^ (1 / 3) * (x ^ (1 / 2) + 1)
                expand polynomial
                apply integral identity
                simplify
            done

            calculate INT x:[-1, 0]. (3 * x ^ 4 + 3 * x ^ 2 + 1) / (x ^ 2 + 1)
                partial fraction decomposition
                apply integral identity
                simplify
            done

            calculate INT x:[4, exp(1) + 3]. (x ^ 3 - 12 * x ^ 2 - 42) / (x - 3)
                partial fraction decomposition
                apply integral identity
                simplify
                substitute u for x - 3
                apply integral identity
                expand polynomial
                simplify
            done

            calculate INT x:[0, pi / 2]. sin(x) * cos(x) ^ 3
                substitute u for cos(x)
                apply integral identity
                simplify
            done

            calculate INT x:[0, pi]. 1 - sin(x) ^ 3
                simplify
                rewrite sin(x) ^ 3 to sin(x) * sin(x) ^ 2
                rewrite sin(x) ^ 2 to 1 - cos(x) ^ 2
                substitute u for cos(x)
                apply integral identity
                simplify
            done

            calculate INT x:[pi/6, pi/2]. cos(x) ^ 2
                apply integral identity
                simplify
            done

            calculate INT x:[0, 1]. (1 - x^2) ^ (1/2)
                substitute sin(u) for x creating u
                rewrite 1 - sin(u) ^ 2 to cos(u) ^ 2
                simplify
                apply integral identity
                simplify
            done

            calculate INT x:[0, sqrt(2)]. sqrt(2 - x^2)
                substitute sqrt(2) * sin(u) for x creating u
                simplify
                rewrite sin(u) ^ 2 to 1 - cos(u) ^ 2
                rewrite -(2 * (1 - cos(u) ^ 2)) + 2 to 2 * cos(u)^2
                simplify
                apply integral identity
                simplify
            done

            calculate INT y:[-sqrt(2), sqrt(2)]. sqrt(8 - 2*y^2)
                substitute 2 * sin(u) for y creating u
                simplify
                rewrite sin(u) ^ 2 to 1 - cos(u) ^ 2
                rewrite -(8 * (1 - cos(u) ^ 2)) + 8 to 8*cos(u)^2
                simplify
                apply integral identity
                expand polynomial
            done

            calculate INT x:[1/sqrt(2), 1]. sqrt(1 - x^2) / x ^ 2
                substitute sin(u) for x creating u
                simplify
                rewrite sin(u) ^ 2 to 1 - cos(u) ^ 2
                simplify
                rewrite cos(u) ^ 2 to 1 - sin(u) ^ 2
                expand polynomial
                simplify
                rewrite 1 / sin(u) ^ 2 to csc(u) ^ 2
                apply integral identity
                simplify
            done

            calculate INT x:[-1, 1]. x / sqrt(5 - 4 * x)
                substitute u for 5 - 4 * x
                expand polynomial
                simplify
                apply integral identity
                simplify
            done

            calculate INT x:[1,4]. 1 / (1 + sqrt(x))
                substitute u for sqrt(x)
                substitute v for u + 1
                expand polynomial
                simplify
                apply integral identity
                simplify
            done

            calculate INT x:[3/4, 1]. 1 / (sqrt(1-x) - 1)
                substitute u for sqrt(1 - x)
                substitute v for u - 1
                expand polynomial
                simplify
                apply integral identity
                simplify
            done

            calculate INT t:[0, 1]. t * exp(-t ^ 2 / 2)
                substitute u for t ^ 2 / 2
                apply integral identity
                simplify
            done

            calculate INT x:[1, exp(2)]. 1 / (x * sqrt(1 + log(x)))
                substitute u for 1 + log(x)
                apply integral identity
                simplify
            done

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
            done

            calculate INT x:[-pi/2, pi/2]. cos(x) ^ 4
                rewrite cos(x) ^ 4 to (cos(x) ^ 2) ^ 2
                rewrite cos(x) ^ 2 to (1 + cos(2*x)) / 2
                expand polynomial
                substitute u for 2 * x
                simplify
                apply integral identity
                simplify
            done

            calculate INT x:[-pi/2, pi/2]. sqrt(cos(x) - cos(x)^3)
                rewrite cos(x) - cos(x)^3 to cos(x) * (1 - cos(x)^2)
                rewrite 1 - cos(x)^2 to sin(x)^2
                simplify
                split region at 0
                simplify
                substitute u for cos(x)
                apply integral identity
                substitute u for cos(x)
                apply integral identity
                simplify
            done

            calculate INT x:[0, pi]. sqrt(1 + cos(2*x))
                rewrite cos(2 * x) to 2 * cos(x) ^ 2 - 1
                simplify
                split region at pi / 2
                simplify
                apply integral identity
                simplify
            done

            calculate INT x:[0, 1]. x * exp(-x)
                integrate by parts with u = x, v = -exp(-x)
                simplify
                apply integral identity
                simplify
            done

            calculate INT x:[1, exp(1)]. x * log(x)
                integrate by parts with u = log(x) / 2, v = x ^ 2
                apply integral identity
                simplify
            done

            calculate INT x:[pi/4, pi/3]. x / sin(x)^2
                rewrite sin(x) ^ 2 to csc(x) ^ -2
                integrate by parts with u = x, v = -cot(x)
                simplify
                rewrite cot(x) to cos(x) / sin(x)
                substitute u for sin(x)
                apply integral identity
                simplify
            done

            calculate INT x:[1, 4]. log(x) / sqrt(x)
                integrate by parts with u = 2 * log(x), v = sqrt(x)
                apply integral identity
                simplify
            done

            calculate INT x:[0, 1]. x * atan(x)
                integrate by parts with u = atan(x) / 2, v = x ^ 2
                simplify
                rewrite x^2 / (2 * x^2 + 2) to (1 - 1 / (x^2 + 1)) / 2
                apply integral identity
                simplify
            done

            calculate INT x:[0, pi/2]. exp(2*x) * cos(x)
                integrate by parts with u = exp(2*x), v = sin(x)
                simplify
                integrate by parts with u = exp(2*x), v = -cos(x)
                simplify
                solve integral INT x:[0, pi/2]. exp(2*x)*cos(x)
            done

            calculate INT x:[0,pi]. (x * sin(x)) ^ 2
                simplify
                rewrite sin(x) ^ 2 to (1 - cos(2*x)) / 2
                expand polynomial
                simplify
                integrate by parts with u = x^2 / 2, v = sin(2*x)
                simplify
                integrate by parts with u = x / 2, v = -cos(2*x)
                simplify
                apply integral identity
                simplify
            done

            calculate INT x:[1, exp(1)]. sin(log(x))
                substitute u for log(x)
                integrate by parts with u = -exp(u), v = cos(u)
                simplify
                integrate by parts with u = exp(u), v = sin(u)
                simplify
                solve integral INT u:[0,1]. exp(u) * sin(u)
            done

            calculate INT x:[1/exp(1), exp(1)]. abs(log(x))
                split region at 1
                simplify
                integrate by parts with u = log(x), v = x
                simplify
                integrate by parts with u = log(x), v = x
                apply integral identity
                simplify
            done
        """
        self.check_actions("base", "tongji", actions)

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

    def testUSubstitution(self):
        with open('theories/usubstitution.thy', 'r', encoding='utf-8') as file:
            actions = file.read()
        self.check_actions("UCDavis", "USubstitution", actions)

    def testUCDavisPartialFraction(self):
        actions = """
            calculate INT x:[3, 4]. 1 / (x ^ 2 - 4)
                partial fraction decomposition
                simplify
                substitute u for 4 * x + 8
                substitute u for 4 * x - 8 (at 2)
                apply integral identity
                simplify
            done

            calculate INT x:[1, 2]. (2 * x + 3) / (x ^ 2 - 9)
                partial fraction decomposition
                simplify
                substitute u for 2 * x - 6
                substitute u for 2 * x + 6 (at 2)
                apply integral identity
                simplify
            done

            calculate INT x:[1, 2]. (2 - x) / (x ^ 2 + 5 * x)
                partial fraction decomposition
                apply integral identity
                substitute u for 5 * x + 25
                apply integral identity
                simplify
            done

            calculate INT x:[1, 2]. (x^2 - 1) / (x^2 - 16)
                partial fraction decomposition
                simplify
                substitute u for 8 * x + 32
                apply integral identity
                substitute u for 8 * x - 32
                apply integral identity
                simplify
            done

            calculate INT x:[3, 4]. (x^4 + x^3 + x^2 + 1)/(x^2 + x - 2)
                partial fraction decomposition
                simplify
                substitute u for 3 * x + 6
                apply integral identity
                simplify
                substitute u for 3 * x - 3
                apply integral identity
                simplify
            done

            calculate INT x:[2, 4]. (x^2 + x - 1) / (x * (x^2 - 1))
                partial fraction decomposition
                simplify
                substitute u for 2 * x + 2
                apply integral identity
                simplify
                substitute u for 2 * x - 2
                apply integral identity
                simplify
            done

            calculate INT x:[2, 4]. (x + 7) / (x ^ 2 * (x + 2))
                partial fraction decomposition
                apply integral identity
                simplify
                substitute u for 4 * x + 8
                apply integral identity
                simplify
            done

            calculate INT x:[2, 4]. (x^5 + 1) / (x ^ 3 * (x + 2))
                partial fraction decomposition
                apply integral identity
                simplify
                substitute u for 8 * x + 16
                apply integral identity
                simplify
            done

            calculate INT x:[0, 1]. (x ^ 2 - x + 1) / (x + 1)^3
                partial fraction decomposition
                substitute u for x + 1
                apply integral identity
                simplify
            done

            calculate INT x:[2, 3]. (x^3 + 4) / ((x ^ 2 - 1) * (x ^ 2 + 3 * x + 2))
                partial fraction decomposition
                simplify
                substitute u for 3 * x + 6
                substitute u for 4 * x + 4 (at 2)
                substitute u for 12 * x - 12 (at 3)
                substitute u for x + 1 (at 4)
                apply integral identity
                simplify
            done

            calculate INT x:[3, 4]. (x ^ 3 + 2 * x - 1) / (x ^ 2 - x - 2) ^ 2
                partial fraction decomposition
                simplify
                substitute u for 27 * x + 27
                substitute u for 27 * x - 54 (at 2)
                substitute u for x + 1 (at 3)
                substitute u for x - 2 (at 4)
                apply integral identity
                simplify
            done

            calculate INT x:[pi/12, pi/6]. sec(x) ^ 2 / (tan(x) ^ 3 - tan(x) ^ 2)
                substitute u for tan(x)
                partial fraction decomposition
                apply integral identity
                simplify
                substitute v for u - 1
                apply integral identity
                simplify
            done

            calculate INT x:[3, 4]. (x ^ 3 + 8) / ((x ^ 2 - 1) * (x - 2))
                partial fraction decomposition
                simplify
                substitute u for 6 * x + 6
                substitute u for 2 * x - 2 (at 2)
                substitute u for 3 * x - 6 (at 3)
                apply integral identity
                simplify
            done

            calculate INT x:[1, 2]. exp(x) / ((exp(x) - 1) * (exp(x) + 3))
                substitute u for exp(x)
                partial fraction decomposition
                simplify
                substitute v for 4 * u + 12
                substitute v for 4 * u - 4 (at 2)
                apply integral identity
                simplify
            done

            calculate INT x:[0, 1]. 1 / (exp(x) + 1)
                rewrite 1 / (exp(x) + 1) to exp(x) / (exp(x) * (exp(x) + 1))
                substitute u for exp(x)
                partial fraction decomposition
                apply integral identity
                simplify
            done

            calculate INT x:[1, 2]. (3 - x) / (x * (x ^2 + 1))
                partial fraction decomposition
                apply integral identity
                simplify
                rewrite (3 * x + 1) / (x ^ 2 + 1) to 3 * x / (x ^ 2 + 1) + 1 / (x ^ 2 + 1)
                apply integral identity
                simplify
                substitute u for x ^ 2 + 1
                apply integral identity
                simplify
            done

            calculate INT x:[1, 2]. (3 * x + 1) / (x ^ 2 * (x ^ 2 + 25))
                partial fraction decomposition
                apply integral identity
                simplify
                expand polynomial
                simplify
                substitute u for 25 * x ^ 2 + 625
                apply integral identity
                simplify
                substitute u for x / 5
                rewrite 5 / (625 * u ^ 2 + 625) to 1/125 * (1 / (u ^ 2 + 1))
                apply integral identity
                simplify
            done

            calculate INT x:[0, 1]. 1 / (x ^ 4 - 16)
                partial fraction decomposition
                simplify
                substitute u for x / 2
                rewrite 2 / (32 * u ^ 2 + 32) to 1/16 * (1 / (u ^ 2 + 1))
                apply integral identity
                simplify
            done

            calculate INT x:[pi/6, pi/3]. cos(x) / (sin(x) ^ 3 + sin(x))
                substitute u for sin(x)
                partial fraction decomposition
                simplify
                substitute v for u ^ 2 + 1
                apply integral identity
                simplify
            done

            calculate INT x:[0, 1]. 1 / (x ^ 4 + 4)
                partial fraction decomposition
                simplify
                rewrite 8 * x ^ 2 + 16 * x + 16 to 8 * ((x + 1) ^ 2 + 1)
                substitute u for x + 1
                rewrite 8 * x ^ 2 - 16 * x + 16 to 8 * ((x - 1) ^ 2 + 1)
                substitute u for x - 1 (at 2)
                expand polynomial
                simplify
                substitute v for u ^ 2 + 1
                apply integral identity
                simplify
                substitute v for u ^ 2 + 1
                apply integral identity
                simplify
                rewrite 8 * u ^ 2 + 8 to 8 * (u ^ 2 + 1)
                apply integral identity
                simplify
                rewrite 8 * u ^ 2 + 8 to 8 * (u ^ 2 + 1)
                apply integral identity
                simplify
            done
        """
        self.check_actions("UCDavis", "PartialFraction", actions)

    def testIntegrateByParts(self):
        actions = """
            calculate INT x:[0, 1]. x*exp(x)
                integrate by parts with u = x, v = exp(x)
                apply integral identity
                simplify
            done

            calculate INT x:[0, pi/2]. x * sin(x)
                integrate by parts with u = x, v = -cos(x)
                apply integral identity
                simplify
            done

            calculate INT x:[1, exp(1)]. x * log(x)
                integrate by parts with u = log(x) / 2, v = x ^ 2
                apply integral identity
                simplify
            done

            calculate INT x:[0, pi/4]. x * cos(3*x)
                substitute u for 3 * x
                simplify
                integrate by parts with u = u, v = sin(u)
                apply integral identity
                simplify
            done

            calculate INT x:[1, exp(1)]. log(x)/x^5
                integrate by parts with u = log(x), v = -1/(4 * x^4)
                apply integral identity
                simplify
            done

            calculate INT x:[0, 1/3]. asin(3*x)
                integrate by parts with u = asin(3*x), v = x
                simplify
                substitute u for -(9 * x ^ 2) + 1
                apply integral identity
                simplify
            done

            calculate INT x:[1, exp(1)]. log(x)
                integrate by parts with u = log(x), v = x
                apply integral identity
                simplify
            done

            calculate INT x:[0, 1]. 2*x*atan(x)
                integrate by parts with u = atan(x), v = x^2
                simplify
                rewrite x ^ 2 / (x ^ 2 + 1) to 1 - 1 / (x ^ 2 + 1)
                apply integral identity
                simplify
            done

            calculate INT x:[1, 2]. x^2*exp(3*x)
                integrate by parts with u = x^2/3, v = exp(3*x)
                simplify
                integrate by parts with u = x/3, v = exp(3*x)
                apply integral identity
                simplify
            done

            calculate INT x:[1, exp(1)]. x^3*log(5*x)
                integrate by parts with u = log(5*x) / 4, v = x^4
                apply integral identity
                simplify
            done

            calculate INT x:[1, exp(1)]. log(x) ^ 2
                integrate by parts with u = log(x) ^ 2, v = x
                simplify
                integrate by parts with u = log(x), v = x
                apply integral identity
                simplify
            done

            calculate INT x:[1, exp(1)]. log(x) ^ 3
                integrate by parts with u = log(x) ^ 3, v = x
                simplify
                integrate by parts with u = log(x) ^ 2, v = x
                simplify
                integrate by parts with u = log(x), v = x
                apply integral identity
                simplify
            done

            calculate INT x:[-2, 1]. x * sqrt(x + 3)
                integrate by parts with u = x, v = 2/3 * (x + 3) ^ (3/2)
                simplify
                substitute u for x + 3
                apply integral identity
                simplify
            done

            calculate INT x:[0, pi/4]. x * sin(x) * cos(x)
                integrate by parts with u = x, v = 1/2 * sin(x) ^ 2
                simplify
                rewrite sin(x) ^ 2 to (1 - cos(2*x)) / 2
                apply integral identity
                simplify
            done

            calculate INT x:[1, exp(1)]. (log(x)/x)^2
                integrate by parts with u = log(x) ^ 2, v = -1/x
                simplify
                integrate by parts with u = log(x), v = -1/x
                apply integral identity
                simplify
            done

            calculate INT x:[0, 1]. x ^ 5 * exp(x ^ 3)
                integrate by parts with u = x ^ 3, v = exp(x^3) / 3
                simplify
                substitute u for x ^ 3
                apply integral identity
                simplify
            done
            
            calculate INT x:[0, pi/4]. x ^ 3 * cos(x ^ 2)
                integrate by parts with u = x ^ 2, v = sin(x^2) / 2
                simplify
                substitute u for x ^ 2
                apply integral identity
                simplify
            done

            calculate INT x:[0, 1]. x ^ 7 * sqrt(5 + 3 * x ^ 4)
                integrate by parts with u = x ^ 4, v = 1/18 * (5 + 3*x^4) ^ (3/2)
                substitute u for 3 * x ^ 4 + 5
                apply integral identity
                simplify
            done

            calculate INT x:[1, 2]. x ^ 3 / (x ^ 2 + 5) ^ 2
                integrate by parts with u = x^2 / 2, v = -1/(x^2 + 5)
                simplify
                substitute u for x ^ 2 + 5
                apply integral identity
                simplify
            done

            calculate INT x:[0, 1]. exp(6 * x) * sin(exp(3 * x))
                integrate by parts with u = exp(3*x), v = -cos(exp(3*x)) / 3
                simplify
                substitute u for exp(3 * x)
                apply integral identity
                simplify
            done

            calculate INT x:[0, 1]. (x ^ 3 * exp(x^2))/(x^2 + 1) ^ 2
                integrate by parts with u = x^2 * exp(x^2) / 2, v = -1/(x^2 + 1)
                simplify
                rewrite (2 * x ^ 3 * exp(x ^ 2) + 2 * x * exp(x ^ 2)) / (2 * x ^ 2 + 2) to x * exp(x^2)
                substitute u for x ^ 2
                apply integral identity
                simplify
            done

            calculate INT x:[0, 1]. exp(x) * cos(x)
                integrate by parts with u = exp(x), v = sin(x)
                simplify
                integrate by parts with u = exp(x), v = -cos(x)
                simplify
                solve integral INT x:[0, 1]. exp(x) * cos(x)
            done

            calculate INT x:[0, pi/4]. sin(3*x) * cos(5*x)
                integrate by parts with u = sin(3*x), v = 1/5*sin(5*x)
                simplify
                integrate by parts with u = cos(3*x), v = -1/5*cos(5*x)
                simplify
                solve integral INT x:[0, pi/4]. sin(3*x) * cos(5*x)
            done
        """
        self.check_actions("UCDavis", "IntegrateByParts", actions)

    def testExponential(self):
        actions = """
            calculate INT x. 5 * exp(x)
                apply integral identity
            done

            calculate INT x. 2 - 3 * exp(x)
                apply integral identity
                simplify
            done

            calculate INT x:[0,log(2) / 7]. 14 * exp(7 * x)
                apply integral identity
                simplify
            done

            calculate INT x. 7 ^ (2 * x + 3)
                substitute u for 2 * x
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. exp(5 * x) * (exp(2 * x) / 7 + 3 / exp(3 * x))
                expand polynomial
                apply integral identity
                simplify
            done

            calculate INT x. exp(x) * (1 + 2 * exp(x)) ^ 4
                substitute u for 1 + 2 * exp(x)
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. (exp(4 * x) - exp(-4 * x)) ^ 2
                expand polynomial
                apply integral identity
                simplify
            done

            calculate INT x. exp(x) * (1 - exp(x)) * (1 + exp(x)) ^ 10
                substitute u for exp(x) + 1
                simplify
                expand polynomial
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x:[0,1]. (3 ^ x + 4 ^ x) / 5 ^ x
                rewrite (3 ^ x + 4 ^ x) / 5 ^ x to 3 ^ x / 5 ^ x + 4 ^ x / 5 ^ x
                rewrite 3 ^ x / 5 ^ x to (3/5) ^ x
                rewrite 4 ^ x / 5 ^ x to (4/5) ^ x
                simplify
                apply integral identity
                simplify
            done

            calculate INT x. 30 * exp(-3 * x) * (1 + 3 * exp(-x)) ^ 5
                substitute u for 3 * exp(-x) + 1
                simplify
                rewrite u ^ 5 * (u - 1) ^ 2 to u ^ 7 - 2 * u ^ 6 + u ^ 5
                simplify
                apply integral identity
                replace substitution
                simplify
            done

            prove (INT x. (27 * exp(9 * x) + exp(12 * x)) ^ (1/3)) = (exp(3 * x) + 27) ^ (4/3) / 4 + SKOLEM_CONST(C)
            lhs:
                rewrite 27 * exp(9 * x) + exp(12 * x) to exp(9 * x) * (27 + exp(3 * x))
                rewrite (exp(9 * x) * (27 + exp(3 * x))) ^ (1/3) to exp(3 * x) * (27 + exp(3 * x)) ^ (1/3)
                substitute u for 27 + exp(3 * x)
                simplify
                apply integral identity
                simplify
                replace substitution
            done
        """
        self.check_actions("UCDavis", "Exponential", actions)

    def testTrigonometric(self):
        actions = """
            calculate INT x. tan(5 * x) for x > 0, x < pi/10
                rewrite tan(5 * x) to sin(5 * x) / cos(5 * x)
                substitute u for cos(5 * x)
                simplify
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. 5 * sec(4 * x) * tan(4 * x)
                rewrite sec(4 * x) to 1 / cos(4 * x)
                rewrite tan(4 * x) to sin(4 * x) / cos(4 * x)
                simplify
                substitute u for cos(4 * x)
                simplify
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. (sin(x) + cos(x)) ^ 2
                expand polynomial
                rewrite 2 * cos(x) * sin(x) + cos(x) ^ 2 + sin(x) ^ 2 to 2 * cos(x) * sin(x) + (sin(x) ^ 2 + cos(x) ^ 2)
                rewrite sin(x) ^ 2 + cos(x) ^ 2 to 1
                simplify
                substitute u for sin(x)
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. 3 * cos(5 * x) ^ 2
                substitute u for 5 * x
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. (2 + tan(x)) ^ 2
                expand polynomial
                simplify
                apply integral identity
                rewrite tan(x) ^ 2 to sec(x) ^ 2 - 1
                simplify
                apply integral identity
                simplify
            done

            calculate INT x. sin(x) ^ 3
                rewrite sin(x) ^ 3 to sin(x) * sin(x) ^ 2
                rewrite sin(x) ^ 2 to 1 - cos(x) ^ 2
                expand polynomial
                apply integral identity
                substitute u for cos(x)
                simplify
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. cos(5 * x) / (3 + sin(5 * x)) for x > 0, x < pi/10
                substitute u for sin(5 * x)
                substitute v for u + 3
                simplify
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. cos(x) ^ 2 / (1 + sin(x))
                rewrite cos(x) ^ 2 to 1 - sin(x) ^ 2
                rewrite 1 - sin(x) ^ 2 to (1 + sin(x)) * (1 - sin(x))
                simplify
                apply integral identity
                simplify
            done

            calculate INT x. sin(x) / (1 + sin(x))
                rewrite sin(x) / (1 + sin(x)) to sin(x) * (1 - sin(x)) / (1 - sin(x) ^ 2)
                rewrite 1 - sin(x) ^ 2 to cos(x) ^ 2
                rewrite sin(x) * (1 - sin(x)) / cos(x) ^ 2 to 1 / cos(x) * (sin(x) / cos(x)) - (sin(x) / cos(x)) ^ 2
                rewrite 1 / cos(x) to sec(x)
                rewrite sin(x) / cos(x) to tan(x)
                rewrite sin(x) / cos(x) to tan(x)
                rewrite tan(x) ^ 2 to sec(x) ^ 2 - 1
                simplify
                apply integral identity
                simplify
            done

            calculate INT x. (csc(3 * x) + cot(3 * x)) ^ 2
                expand polynomial
                rewrite cot(3 * x) ^ 2 to csc(3 * x) ^ 2 - 1
                simplify
                substitute u for 3 * x
                simplify
                apply integral identity
                substitute v for 3 * x
                simplify
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. sec(x) ^ 2 * sqrt(5 + tan(x))
                substitute u for 5 + tan(x)
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. tan(x) ^ 5
                rewrite tan(x) ^ 5 to tan(x) ^ 3 * tan(x) ^ 2
                rewrite tan(x) ^ 2 to sec(x) ^ 2 - 1
                expand polynomial
                simplify
                substitute u for tan(x)
                apply integral identity
                simplify
                rewrite tan(x) ^ 3 to tan(x) * tan(x) ^ 2
                rewrite tan(x) ^ 2 to sec(x) ^ 2 - 1
                rewrite tan(x) * (sec(x) ^ 2 - 1) to tan(x) * sec(x) ^ 2 - tan(x)
                simplify
                substitute u for tan(x)
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. (sin(x) - cos(x)) * (sin(x) + cos(x)) ^ 5
                rewrite (sin(x) - cos(x)) * (sin(x) + cos(x)) ^ 5 to -((sin(x) + cos(x)) ^ 5) * (cos(x) - sin(x))
                substitute u for sin(x) + cos(x)
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. cos(x) * exp(4 + sin(x))
                rewrite cos(x) * exp(4 + sin(x)) to exp(4 + sin(x)) * cos(x)
                substitute u for 4 + sin(x)
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. sin(3 * x) * sin(cos(3 * x))
                rewrite sin(3 * x) * sin(cos(3 * x)) to -sin(cos(3 * x)) / 3 * -(3 * sin(3 * x))
                substitute u for cos(3 * x)
                simplify
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. cos(x) * log(sin(x)) / sin(x)
                rewrite cos(x) * log(sin(x)) / sin(x) to log(sin(x)) * (cos(x) / sin(x))
                substitute u for log(sin(x))
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. sec(x) * tan(x) * sqrt(4 + 3 * sec(x))
                rewrite sec(x) * tan(x) * sqrt(4 + 3 * sec(x)) to sqrt(4 + 3 * sec(x)) / 3 * (3 * sec(x) * tan(x))
                substitute u for 4 + 3 * sec(x)
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. exp(x) * cos(exp(x))
                substitute u for exp(x)
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. x * sin(3 * x)
                integrate by parts with u = x, v = -1/3 * cos(3 * x)
                simplify
                apply integral identity
                simplify
            done

            calculate INT x. x ^ 2 * cos(x)
                integrate by parts with u = x ^ 2, v = sin(x)
                simplify
                integrate by parts with u = x, v = -cos(x)
                simplify
                apply integral identity
                simplify
            done

            calculate INT x. sin(x) * cos(x) * exp(sin(x))
                rewrite sin(x) * cos(x) * exp(sin(x)) to sin(x) * exp(sin(x)) * cos(x)
                substitute u for sin(x)
                integrate by parts with u = u, v = exp(u)
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. exp(x) * sin(x)
                integrate by parts with u = exp(x), v = -cos(x)
                simplify
                integrate by parts with u = exp(x), v = sin(x)
                solve integral INT x. exp(x) * sin(x)
            done

            calculate INT x. cos(4 * x) * sin(3 * x)
                integrate by parts with u = sin(3 * x), v = sin(4 * x) / 4
                simplify
                integrate by parts with u = cos(3 * x), v = -cos(4 * x) / 4
                simplify
                solve integral INT x. sin(3 * x) * cos(4 * x)
            done

            calculate INT x. sec(x) * sqrt(sec(x) + tan(x))
                rewrite sec(x) * sqrt(sec(x) + tan(x)) to 2 * ((sec(x) * tan(x) + sec(x) ^ 2) / (2 * sqrt(sec(x) + tan(x))))
                substitute u for sqrt(sec(x) + tan(x))
                apply integral identity
                replace substitution
            done

            calculate INT x. (sin(2 * x) - cos(2 * x)) / (sin(2 * x) + cos(2 * x))
                rewrite (sin(2 * x) - cos(2 * x)) / (sin(2 * x) + cos(2 * x)) to -1 / (2 * (sin(2 * x) + cos(2 * x))) * (2 * cos(2 * x) - 2 * sin(2 * x))
                substitute u for sin(2 * x) + cos(2 * x)
                simplify
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. (sin(x) + cos(x)) / (exp(-x) + sin(x))
                rewrite (sin(x) + cos(x)) / (exp(-x) + sin(x)) to 1 / (1 + exp(x) * sin(x)) * (cos(x) * exp(x) + exp(x) * sin(x))
                substitute u for 1 + exp(x) * sin(x)
                apply integral identity
                replace substitution
                simplify
            done
        """
        self.check_actions("UCDavis", "Trigonometric", actions)

    def testLogAndArcTangent(self):
        actions = """
            calculate INT x:[1, 2]. (3 + x^2)/x
                partial fraction decomposition
                apply integral identity
                simplify
            done

            calculate INT x:[0, 1]. 7/(x + 5)
                substitute u for x + 5
                apply integral identity
                simplify
            done

            calculate INT x:[-1, 1]. 3/(x^2 + 4)
                substitute u for x / 2
                rewrite 4 * u ^ 2 + 4 to 4 * (u ^ 2 + 1)
                apply integral identity
                simplify
            done

            calculate INT x:[-1, 0]. x/(x^2 + 4)
                substitute u for x ^ 2 + 4
                apply integral identity
                simplify
            done

            calculate INT x:[-1, 1].  x^3/(x + 2)
                partial fraction decomposition
                apply integral identity
                simplify
            done

            calculate INT x:[-1, 1]. (x^4 + x^3)/(x^2 + 1)
                partial fraction decomposition
                apply integral identity
                simplify
                expand polynomial
                simplify
                split region at 0
                substitute u for x ^ 2 + 1
                substitute u for x ^ 2 + 1 (at 2)
                apply integral identity
                simplify
            done

            calculate INT x:[1, exp(1)]. 1 / (x * (3 + log(x)))
                substitute u for 3 + log(x)
                apply integral identity
                simplify
            done

            calculate INT x:[1, exp(1)]. (7 - log(x)) / (x * (3 + log(x)))
                substitute u for log(x)
                substitute v for u + 3
                expand polynomial
                apply integral identity
                simplify
            done

            calculate INT x:[0, 1]. 1 / (x^2 + 4 * x + 5)
                rewrite x ^ 2 + 4 * x + 5 to (x + 2) ^ 2 + 1
                substitute u for x + 2
                apply integral identity
                simplify
            done

            calculate INT x:[1, 2]. 1 / (2*x^2 + 12*x + 68)
                rewrite 2*x^2 + 12*x + 68 to 2 * (x+3)^2 + 50
                substitute u for x + 3
                substitute v for u / 5
                rewrite 50 * v ^ 2 + 50 to 50 * (v ^ 2 + 1)
                apply integral identity
                simplify
            done

            calculate INT x:[9, 16]. 1 / (sqrt(x) * (x + 9))
                substitute u for sqrt(x)
                substitute v for u / 3
                rewrite 9 * v^2 + 9 to 9 * (v^2 + 1)
                apply integral identity
                simplify
            done

            calculate INT x:[0, pi/4]. sin(2*x)/(1 + cos(2*x))
                substitute u for 1 + cos(2*x)
                apply integral identity
                simplify
            done

            calculate INT x:[0, pi/2]. cos(x) / (1 + sin(x) ^ 2)
                substitute u for sin(x)
                apply integral identity
                simplify
            done

            calculate INT x:[0, 1]. exp(2*x)/(1 + exp(2*x))
                substitute u for exp(2*x)
                substitute v for 2 * u + 2
                apply integral identity
                simplify
            done

            calculate INT x:[0, 1]. exp(x)/(1 + exp(2*x))
                rewrite exp(2*x) to exp(x) * exp(x)
                substitute u for exp(x)
                apply integral identity
                simplify
            done

            calculate INT x:[0, 1]. exp(4*x)/(1 + exp(2*x))
                substitute u for exp(2*x)
                partial fraction decomposition
                simplify
                substitute v for 2 * u + 2
                apply integral identity
                simplify
            done

            calculate INT x:[0, 1]. x / (x^2 + 2 * x + 1)
                rewrite x^2 + 2*x + 1 to (x+1)^2
                substitute u for x + 1
                expand polynomial
                apply integral identity
                simplify
            done

            calculate INT x:[0, 1]. x / (x ^ 2 + 2 * x + 5)
                rewrite x^2 + 2*x + 5 to (x+1)^2 + 4
                substitute u for x + 1
                expand polynomial
                simplify
                substitute v for u ^ 2 + 4
                apply integral identity
                substitute v for u / 2
                rewrite 4 * v ^ 2 + 4 to 4 * (v ^ 2 + 1)
                apply integral identity
                simplify
            done

            calculate INT x:[-1, 1]. (7 * x - 2) / (2*x^2 - 16 * x + 42)
                rewrite 2*x^2 - 16*x + 42 to 2*(x-4)^2 + 10
                substitute u for x - 4
                expand polynomial
                simplify
                substitute v for 2 * u ^ 2 + 10
                apply integral identity
                simplify
                substitute v for u / sqrt(5)
                rewrite 10 * v^2 + 10 to 10 * (v^2 + 1)
                apply integral identity
                simplify
            done

            calculate INT x:[0, 1]. 1 / (1 + exp(2*x))
                rewrite 1 / (1 + exp(2*x)) to exp(-2*x) / (exp(-2*x) + 1)
                substitute u for exp(-2*x) + 1
                rewrite (u - 1) / (u * (2 * u - 2)) to 1/(2*u)
                apply integral identity
                simplify
            done

            calculate INT x:[0, 1]. exp(3*x)/(1 + exp(2*x))
                rewrite exp(3*x) / (1 + exp(2*x)) to exp(x) * exp(x) * exp(x) / (1 + exp(x) * exp(x))
                substitute u for exp(x)
                partial fraction decomposition
                apply integral identity
                simplify
            done

            calculate INT x:[1, 4]. (9 + 6*sqrt(x) + x)/(4*sqrt(x) + x)
                rewrite 4*sqrt(x) + x to sqrt(x) * (4 + sqrt(x))
                substitute (u - 4)^2 for x creating u
                substitute t for u - 4
                simplify
                partial fraction decomposition
                apply integral identity
                simplify
            done
        """
        self.check_actions("UCDavis", "LogAndArcTangent", actions)

    def testPowerSubstitution(self):
        actions = """
            calculate INT x. 1 / (1 + sqrt(x))
                substitute u for sqrt(x)
                partial fraction decomposition
                simplify
                substitute v for u + 1
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. (2 + sqrt(x)) / (3 - sqrt(x))
                substitute u for sqrt(x)
                partial fraction decomposition
                simplify
                substitute v for u - 3
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. 3 / (4 + x ^ (1/3))
                substitute u for x ^ (1/3)
                simplify
                partial fraction decomposition
                simplify
                substitute v for u + 4
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. x * (x - 1) ^ (1/6) for x > 1
                substitute u for (x - 1) ^ (1/6)
                simplify
                expand polynomial
                simplify
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. (3 * x + 2) / sqrt(x - 9) for x > 9
                substitute u for sqrt(x - 9)
                simplify
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. 1 / (x ^ (2/3) - x ^ (1/3))
                substitute u for x ^ (1/3)
                partial fraction decomposition
                substitute v for u - 1
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. (sqrt(x) + 1) / (sqrt(x) * (x ^ (1/3) + 1)) for x > 0
                substitute u for x ^ (1/6)
                simplify
                partial fraction decomposition
                simplify
                expand polynomial
                simplify
                substitute v for u ^ 2 + 1
                simplify
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. sqrt(5 + sqrt(x)) for x >= 0
                substitute u for sqrt(5 + sqrt(x))
                simplify
                expand polynomial
                simplify
                apply integral identity
                simplify
                replace substitution
                simplify
            done
                
            calculate INT x. (1 + sqrt(x - 3)) ^ (1/3) for x > 3
                substitute u for (1 + sqrt(x - 3)) ^ (1/3)
                simplify
                expand polynomial
                simplify
                apply integral identity
                simplify
                replace substitution
                simplify
            done

            calculate INT x. sqrt(2 + sqrt(4 + sqrt(x))) for x > 0
                substitute u for sqrt(2 + sqrt(4 + sqrt(x)))
                simplify
                expand polynomial
                simplify
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. 1 / sqrt(2 + sqrt(1 + sqrt(x))) for x > 0
                substitute u for sqrt(2 + sqrt(1 + sqrt(x)))
                simplify
                expand polynomial
                simplify
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. sqrt(x) / (x - 1)
                substitute u for sqrt(x)
                simplify
                partial fraction decomposition
                simplify
                substitute v for 2 * u + 2
                substitute w for 2 * u - 2 (at 2)
                simplify
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. sqrt(4 - x) / x ^ 2 for x < 4
                substitute u for sqrt(4 - x)
                simplify
                partial fraction decomposition
                simplify
                substitute v1 for 8 * u + 16
                substitute v2 for 8 * u - 16 (at 2)
                substitute v3 for u + 2 (at 3)
                substitute v4 for u - 2 (at 4)
                simplify
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. (x ^ (1/4) + 5) / (x - 16) for x > 0
                substitute u for x ^ (1/4)
                simplify
                partial fraction decomposition
                simplify
                expand polynomial
                simplify
                substitute v1 for 2 * u ^ 2 + 8
                substitute v2 for u / 2 (at 2)
                rewrite 8 * v2 ^ 2 + 8 to 8 * (v2 ^ 2 + 1)
                substitute v3 for 4 * u + 8 (at 3)
                substitute v4 for 4 * u - 8 (at 4)
                simplify
                apply integral identity
                replace substitution
                simplify
            done
        """
        self.check_actions("UCDavis", "PowerSubstitution", actions)

    def testTrigSubstitution(self):
        actions = """
            calculate INT x. sqrt(1 - x^2) for x > -1, x < 1
                substitute u for asin(x)
                rewrite -(sin(u) ^ 2) + 1 to 1 - sin(u)^2
                rewrite 1 - sin(u)^2 to cos(u)^2
                simplify
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. (x^2 - 1)^(3/2) / x for x > 1
                substitute u for asec(x)
                rewrite sec(u)^2 - 1 to tan(u)^2
                simplify
                rewrite tan(u)^4 to tan(u)^2 * tan(u)^2
                rewrite tan(u)^2 to sec(u)^2 - 1 (at 2)
                expand polynomial
                rewrite tan(u)^2 to sec(u)^2 - 1 (at 2)
                simplify
                apply integral identity
                substitute v for tan(u)
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. 1 / (1 - x^2) ^ (3/2) for x > -1, x < 1
                substitute u for asin(x)
                rewrite -(sin(u)^2) + 1 to 1 - sin(u)^2
                rewrite 1 - sin(u)^2 to cos(u)^2
                simplify
                rewrite 1 / cos(u)^2 to (1/cos(u))^2
                rewrite 1/cos(u) to sec(u)
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. sqrt(x^2 + 1) / x
                substitute u for atan(x)
                rewrite tan(u)^2 + 1 to sec(u)^2
                simplify
                rewrite sec(u)^3 to sec(u) * sec(u)^2
                rewrite sec(u)^2 to 1 + tan(u)^2
                rewrite sec(u) * (1 + tan(u)^2) / tan(u) to (1 + tan(u)^2) / sin(u)
                expand polynomial
                rewrite 1 / sin(u) to csc(u)
                rewrite tan(u) ^ 2 / sin(u) to sec(u) * tan(u)
                simplify
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. x^3 * sqrt(4 - 9*x^2) for x > -2/3, x < 2/3
                substitute u for asin(3*x/2)
                simplify
                rewrite -(4 * sin(u)^2) + 4 to 4 * (1 - sin(u)^2)
                rewrite 1 - sin(u)^2 to cos(u)^2
                simplify
                rewrite cos(u)^2 * sin(u)^3 to sin(u) * sin(u)^2 * cos(u)^2
                rewrite sin(u)^2 to 1 - cos(u)^2
                expand polynomial
                simplify
                substitute v for cos(u)
                substitute v2 for cos(u) (at 2)
                simplify
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. sqrt(1 - x^2) / x for x > 0, x < 1
                substitute u for asin(x)
                rewrite -(sin(u)^2) + 1 to 1 - sin(u)^2
                rewrite 1 - sin(u)^2 to cos(u)^2
                simplify
                rewrite cos(u) ^ 2 to 1 - sin(u)^2
                expand polynomial
                rewrite 1 / sin(u) to csc(u)
                simplify
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. sqrt(x^2 - 9) / x^2 for x > 3
                substitute u for asec(x/3)
                rewrite 9 * sec(u)^2 - 9 to 9 * (sec(u)^2 - 1)
                rewrite sec(u)^2 - 1 to tan(u)^2
                simplify
                rewrite tan(u)^2 to sec(u)^2 - 1
                expand polynomial
                rewrite 1 / sec(u) to cos(u)
                simplify
                apply integral identity
                replace substitution
                simplify
            done

            calculate INT x. sqrt(x^2 + 1) / x^2 for x > 0
                substitute u for atan(x)
                rewrite tan(u)^2 + 1 to sec(u)^2
                simplify
                rewrite sec(u)^3 to sec(u) * sec(u)^2
                rewrite sec(u)^2 to 1 + tan(u)^2
                expand polynomial
                simplify
                rewrite sec(u) / tan(u)^2 to cot(u) * csc(u)
                apply integral identity
                replace substitution
                simplify
            done
            
            calculate INT x. sqrt(x^2+25)
                substitute u for atan(x/5)
                simplify
                rewrite sqrt(25 * tan(u) ^ 2 + 25) to 5 * sqrt(tan(u)^2+1)
                rewrite tan(u)^2+1 to sec(u)^2
                simplify
                integrate by parts with u=sec(u), v=tan(u)
                rewrite tan(u)^2 to sec(u)^2 - 1
                expand polynomial
                apply integral identity
                solve integral 25 * INT u. sec(u)^3
                expand polynomial
                simplify
                replace substitution
                simplify
            done
            
            calculate INT x. sqrt(x^2-4) for x > 2
                substitute u for asec(x/2)
                rewrite sqrt(4 * sec(u) ^ 2 - 4) to 2 * sqrt(sec(u) ^ 2 - 1)
                rewrite sec(u) ^ 2 - 1 to tan(u)^2
                simplify
                rewrite tan(u)^2 to sec(u)^2 - 1
                expand polynomial
                simplify
                integrate by parts with u=sec(u), v=tan(u)
                solve integral 4*(INT u. sec(u)*tan(u)^2)
                expand polynomial
                apply integral identity
                simplify
                replace substitution
                simplify
            done
            
            calculate INT x. x / sqrt(x^4-16) for x > 2
                substitute u for x^2
                simplify
                substitute v for asec(u/4)
                rewrite sqrt(16 * sec(v) ^ 2 - 16) to 4*sqrt(sec(v)^2-1)
                rewrite sec(v)^2-1 to tan(v)^2
                simplify
                apply integral identity
                replace substitution
                simplify
            done
            
            calculate INT x. 1 / sqrt(x^2-4*x) for x > 4
                rewrite x^2-4*x to (x-2)^2 - 4
                substitute u for asec((x-2)/2)
                rewrite sqrt(4 * sec(u) ^ 2 - 4) to 2 * sqrt(sec(u) ^ 2 - 1)
                rewrite sec(u) ^ 2 - 1 to tan(u)^2
                simplify
                apply integral identity
                replace substitution
                simplify
            done
            
            calculate INT x. x/sqrt(x^2 + 4*x + 5)
                rewrite x^2+4*x+5 to (x+2)^2 + 1
                substitute u for (x+2)
                substitute v for atan(u)
                rewrite tan(v)^2+1 to sec(v)^2
                simplify
                expand polynomial
                apply integral identity
                replace substitution
                simplify
            done
            
            calculate INT x. x * sqrt(10*x - x^2) for x > 0, x < 10
                rewrite 10*x - x^2 to 25 - (x-5)^2
                substitute u for (x-5) / 5
                substitute v for asin(u)
                rewrite sqrt(-(25*sin(v)^2)+25) to 5*sqrt(1-sin(v)^2)
                rewrite 1-sin(v)^2 to cos(v)^2
                simplify
                expand polynomial
                simplify
                substitute w for cos(v)
                apply integral identity
                simplify
                replace substitution
                simplify
            done
            
            calculate INT x. sqrt((x-1) / x) for x > 1
                substitute u for sqrt(x)
                simplify
                substitute v for asec(u)
                rewrite sec(v)^2-1 to tan(v)^2
                simplify
                integrate by parts with u=tan(v),v=sec(v)
                rewrite sec(v)^3 to sec(v)*sec(v)^2
                rewrite sec(v)^2 to tan(v)^2 + 1
                expand polynomial
                expand polynomial
                simplify
                solve integral 2 * INT v. sec(v)*tan(v)^2
                expand polynomial
                apply integral identity
                simplify
                replace substitution
                simplify
            done
            
            calculate INT x. sqrt(1-x)*sqrt(x+3) for x < 1, x > -3
                rewrite sqrt(1-x)*sqrt(x+3) to sqrt((1-x)*(x+3))
                rewrite (1-x)*(x+3) to 4 - (x+1)^2
                substitute u for(x+1)/2
                rewrite sqrt(-(4 * u ^ 2) + 4) to 2 * sqrt(1 - u ^ 2)
                substitute v for asin(u)
                rewrite -(sin(v) ^ 2) + 1 to 1-sin(v)^2
                rewrite 1-sin(v) ^ 2 to cos(v)^2
                simplify
                integrate by parts with u=cos(v),v=sin(v)
                simplify
                rewrite sin(v)^2 to 1-cos(v)^2
                simplify
                solve integral 4 * (INT v. cos(v) ^ 2)
                expand polynomial
                apply integral identity
                replace substitution
                simplify
            done
        """
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
                        substitute sqrt(b) * u for x creating u
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
                substitute 1 / u for x creating u
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
                substitute 1/t for x creating t
                rewrite log(1/t) to -log(t)
                rewrite -log(t) / ((1/t)^2 + b^2) * -(1/t^2) to log(t) / (1 + b^2*t^2)
                substitute s/b for t creating s
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
                rewrite -(u ^ 2) - 1 to -(u^2 + 1)
                simplify
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
                substitute tan(u) for x creating u
                rewrite sec(u) ^ 2 to tan(u) ^ 2 + 1
                simplify
            done
            subgoal 2: (INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)) = pi / 4 * log(2) - (INT x:[0,1]. log(x + 1) / (x ^ 2 + 1))
            lhs:
                apply 1 on INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)
                substitute pi / 4 - y for x creating y
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
                substitute t / a for x creating t
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

    def testLeibniz01(self):
        actions = """
            prove (INT x:[0,oo]. 1 / (x ^ 2 + a ^ 2) ^ 3) = 3 * pi / (16 * a ^ 5) for a > 0
            subgoal 1: (INT x:[0,oo]. 1 / (x ^ 2 + a ^ 2)) = pi / (2 * a) for a > 0
            lhs:
                substitute a * u for x creating u
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
                rewrite sin(2 * x) to 2 * sin(x) * cos(x)
                rewrite a * (2 * sin(x) * cos(x)) to 2 / a * (a * sin(x)) * (a * cos(x))
                rewrite log(2 / a * (a * sin(x)) * (a * cos(x))) to log(2 / a * (a * sin(x))) + log(a * cos(x))
                rewrite log(2 / a * (a * sin(x))) to log(2 / a) + log(a * sin(x))
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
                rewrite log(2 / a) to log(2) + log(1 / a)
                expand polynomial
                simplify
            rhs:
                rewrite log(a / 2) to log(a) - log(2)
                expand polynomial
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
                substitute tan(u) for x creating u
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
                substitute 1 / u for x creating u
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

    def testFrullaniIntegral01(self):
        actions = """
            prove (INT x:[0,oo]. (atan(a * x) - atan(b * x)) / x) = pi * log(a) / 2 - pi * log(b) / 2 for a > 0, b > 0
            define I(a,b) = (INT x:[0,oo]. (atan(a * x) - atan(b * x)) / x) for a > 0, b > 0
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
            prove (INT x:[0,1]. atan(x) / x) = G
            subgoal 1: converges(SUM(n, 0, oo, INT x:[0,1]. x ^ (2 * n) / (2 * n + 1)))
            arg:
                simplify
                apply integral identity
                simplify
            done
            lhs:
                apply series expansion on atan(x) index n
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
            prove (INT x:[0,pi]. x * sin(x) / (a + b * cos(x) ^ 2)) = pi / sqrt(a * b) * atan(sqrt(b / a)) for a > 0, b > 0
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
                rewrite -(a * x ^ 2) - a to -a * (x ^ 2 + 1)
                simplify
                apply integral identity
                simplify
                rewrite atan(-(sqrt(b) / sqrt(a))) to -atan(sqrt(b) / sqrt(a))
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
            prove (INT x:[0,1]. atan(sqrt(2 + x ^ 2)) / ((1 + x ^ 2) * sqrt(2 + x ^ 2))) = 5 * pi ^ 2 / 96
            define I(u) = (INT x:[0,1]. atan(u * sqrt(2 + x ^ 2)) / ((1 + x ^ 2) * sqrt(2 + x ^ 2))) for u > 0
            subgoal 1: I(1) = (INT x:[0,1]. atan(sqrt(x ^ 2 + 2)) / ((x ^ 2 + 1) * sqrt(x ^ 2 + 2)))
            lhs:
                expand definition for I
            done
            subgoal 2: (D u. I(u)) = 1 / (1 + u ^ 2) * (pi / 4 - u / sqrt(1 + 2 * u ^ 2) * atan(u / sqrt(1 + 2 * u ^ 2))) for u > 0
            lhs:
                expand definition for I (all)
                exchange derivative and integral
                simplify
                rewrite 1 / ((x ^ 2 + 1) * (u ^ 2 * (x ^ 2 + 2) + 1)) to 1 / (u ^ 2 + 1) * (1 / (1 + x ^ 2) - u ^ 2 / (1 + 2 * u ^ 2 + u ^ 2 * x ^ 2))
                simplify
                rewrite 1 / (u ^ 2 * x ^ 2 + 2 * u ^ 2 + 1) to u ^ (-2) * (x ^ 2 + (2 * u ^ 2 + 1) / u ^ 2) ^ (-1)
                simplify
                substitute y * sqrt(u ^ (-2) * (2 * u ^ 2 + 1)) for x creating y
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
                integrate by parts with u = 1, v = atan(x / sqrt(2 + x ^ 2)) / 2
                simplify
            done
            subgoal 4: (INT u:[1,oo]. D u. I(u)) = -(pi ^ 2 / 48) + I(1)
            lhs:
                apply 2 on D u. I(u)
                expand polynomial
                simplify
                substitute 1 / x for u creating x
                simplify
                rewrite x ^ 3 * (1 / x ^ 2 + 1) * sqrt(2 / x ^ 2 + 1) to sqrt((1 + x ^ 2) ^ 2 * (2 + x ^ 2))
                rewrite x * sqrt(2 / x ^ 2 + 1) to sqrt(x ^ 2 + 2)
                simplify
                rewrite 1 / sqrt(x ^ 2 + 2) to sqrt(x ^ 2 + 2) ^ (-1)
                rewrite atan(sqrt(x ^ 2 + 2) ^ (-1)) to pi / 2 - atan(sqrt(x ^ 2 + 2))
                expand polynomial
                simplify
                rewrite atan(sqrt(x ^ 2 + 2)) / (x ^ 2 * sqrt(x ^ 2 + 2) + sqrt(x ^ 2 + 2)) to atan(sqrt(x ^ 2 + 2)) / ((x ^ 2 + 1) * sqrt(x ^ 2 + 2))
                apply 1 on INT x:[0,1]. atan(sqrt(x ^ 2 + 2)) / ((x ^ 2 + 1) * sqrt(x ^ 2 + 2))
                integrate by parts with u = 1, v = atan(x / sqrt(2 + x ^ 2))
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
                substitute t for a * x
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

    def testChapter3Practice06(self):
        actions = """
            prove (INT x:[-1,1]. ((1 + x) / (1 - x)) ^ (1/2)) = pi
            lhs:
                substitute cos(2 * u) for x creating u
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
                substitute t for -u + 1
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
                substitute sin(x) ^ 2 for y creating x
                rewrite sin(x) ^ 2 to 1 - cos(x) ^ 2 (at 2)
                simplify
                apply integral identity
                simplify
            done
            subgoal 2: (INT y:[0,1]. log(y) / (sqrt(y) * sqrt(1 - y))) = -(2 * pi * log(2))
            lhs:
                substitute sin(x) ^ 2 for y creating x
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


if __name__ == "__main__":
    unittest.main()
