"""Unit test for complex number integration."""

from integral import parser
from integral import compstate
from integral import state
import unittest


class ActionTest(unittest.TestCase):
    def check_actions(self, base_file: str, current_file: str, actions: str,
                      *, print_lines=False, print_state=False):
        file = compstate.CompFile(base_file, current_file)
        st = state.InitialState(file)
        actions = [s for s in actions.split('\n') if s.strip()]
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
            try:
                st = st.process_action(a)
            except Exception as e:
                print(cur_goal)
                print(st)
                raise e
        if print_state:
            print(st)
        if not print_state and not isinstance(st, state.InitialState):
            raise AssertionError("Does not end in initial state (add print_state=True to debug)")

    def testEulerFormula(self):
        actions = """
            prove (INT x:[0,oo]. sin(b*x)*exp(-x*y)) = b/(y^2+b^2) for y>0
            subgoal 1: (exp(i*b*x)-exp(-i*b*x))/(2*i)=sin(b*x)
            lhs:
                rewrite exp(i*b*x) to cos(b*x)+i*sin(b*x)
                rewrite exp(-i*b*x) to cos(b*x)-i*sin(b*x)
                simplify
            done
            subgoal 2: (LIM {n -> oo}. exp(-y*n)*exp(b * i * n)) = 0 for y>0
            lhs:
                rewrite exp(b * i * n) to cos(b*n)+i*sin(b*n)
                rewrite  exp(-y * n) * (cos(b * n) + i * sin(b * n)) to exp(-(n * y)) * i * sin(b * n) + exp(-(n * y)) * cos(b * n)
                simplify
            done
            subgoal 3: (LIM {n -> oo}. exp(-y*n)*exp(-b * i * n)) = 0 for y>0
            lhs:
                rewrite exp(-b * i * n) to cos(b*n)-i*sin(b*n)
                rewrite  exp(-y * n) * (cos(b * n) - i * sin(b * n)) to -exp(-(n * y)) * i * sin(b * n) + exp(-(n * y)) * cos(b * n)
                simplify
            done
            lhs:
                apply 1 on sin(b * x)
                rewrite (exp(i * b * x) - exp(-i * b * x)) / (2 * i) * exp(-x * y) to  exp(-x * y)*(exp(i * b * x) - exp(-i * b * x)) / (2 * i)
                rewrite exp(-x * y) * (exp(i * b * x) - exp(-i * b * x)) to exp(-x * y) * exp(i * b * x) - exp(-x * y) * exp(-i * b * x)
                rewrite exp(-x * y) * exp(i * b * x) to exp(-x*(y-i*b))
                rewrite exp(-x * y) * exp(-i * b * x) to exp(-x*(y+i*b))
                simplify
                substitute u for -(x * (-(b * i) + y))
                apply integral identity
                substitute r for -(x * (b * i + y))
                apply integral identity
                simplify
                rewrite exp(-(oo * (-(b * i) + y))) to exp(-( - oo *(b * i) + oo * y))
                simplify
                rewrite exp(b * i * oo - y * oo) to exp(-y*oo)*exp(b*i*oo)
                rewrite exp(-y * oo) * exp(b * i * oo) to LIM {n -> oo}. exp(-y*n)*exp(b * i * n)
                rewrite exp(-(oo * (b * i + y))) to exp(-(oo *(b * i) + oo * y))
                simplify
                rewrite exp(-(b * i * oo) - y * oo) to exp(-y*oo)*exp(-b*i*oo)
                rewrite exp(-y * oo) * exp(-b * i * oo) to LIM {n -> oo}. exp(-y*n)*exp(-b * i * n)
                rewrite exp(b * n * i - n * y) to exp(-y*n)*exp(b * i * n)
                apply 2 on (LIM {n -> oo}. exp(-y*n)*exp(b * i * n))
                apply 3 on (LIM {n -> oo}. exp(-y*n)*exp(-b * i * n))
                simplify
                rewrite 1 / (2 * i) * (1 / (-(b * i) - y) - 1 / (b * i - y)) to  (1 / (-(b * i) - y) - 1 / (b * i - y)) * (1 / (2 * i))
                rewrite 1 / (-(b * i) - y) to -1/(b*i+y)
                rewrite -1 / (b * i + y) - 1 / (b * i - y) to - 1 / (b * i - y) - 1 / (b * i + y) 
                rewrite -1 / (b * i - y) to 1 / (y-b*i)
                rewrite 1 / (y - b * i) to (y+b*i) / ((y - b * i)*(y+b*i))
                rewrite 1 / (b * i + y) to (y-b*i) / ((b * i + y)*(y-b*i))
                rewrite (b * i + y)*(y-b*i) to (y - b * i)*(y+b*i)
                rewrite (y + b * i) / ((y - b * i) * (y + b * i)) - (y - b * i) / ((y - b * i) * (y + b * i)) to ((y + b * i) - (y - b * i)) / ((y - b * i) * (y + b * i)) 
                rewrite (y - b * i) * (y + b * i) to y*y+y*b*i-b*i*y-b*i*b*i
                simplify
                rewrite to b/(y^2+b^2)
            done
        """
        self.check_actions("interesting", "", actions)

if __name__ == "__main__":
    unittest.main()
