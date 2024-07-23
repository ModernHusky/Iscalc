"""Parsing"""

from typing import Tuple
from lark import Lark, Token, Transformer, v_args, exceptions
from decimal import Decimal
from fractions import Fraction

from integral import expr
from integral.expr import Expr


grammar = r"""
    ?atom: CNAME -> var_expr
        | "?" CNAME  -> symbol_expr
        | INT -> int_expr
        | DECIMAL -> decimal_expr
        | "D" CNAME "." expr -> deriv_expr
        | "pi" -> pi_expr
        | "G" -> g_expr
        | "inf" -> pos_inf_expr
        | "oo" -> pos_inf_expr
        | "-inf" -> neg_inf_expr
        | "-oo" -> neg_inf_expr
        | CNAME "(" expr ("," expr)* ")" -> fun_expr
        | "(" expr ")"
        | "\|" expr "\|" -> abs_expr 
        | "INT" CNAME ":[" expr "," expr "]." expr -> integral_expr
        | "INT" CNAME "." expr -> indefinite_integral_expr
        | "INT" CNAME "[" CNAME ("," CNAME)* "]" "." expr -> indefinite_integral_skolem_expr
        | "[" expr "]_" CNAME "=" expr "," expr -> eval_at_expr
        | "LIM" "{" CNAME "->" expr "}" "." expr -> limit_inf_expr
        | "LIM" "{" CNAME "->" expr "-}" "."  expr -> limit_l_expr
        | "LIM" "{" CNAME "->" expr "+}" "."  expr -> limit_r_expr
        | "[" expr ("," expr)* "]" -> vector_expr

    ?uminus: "-" uminus -> uminus_expr | atom  // priority 80

    ?pow: pow "^" uminus -> pow_expr          // priority 75
        | "-" atom "^" uminus -> uminus_pow_expr
        | uminus

    ?times: times "*" pow -> times_expr        // priority 70
        | times "/" pow -> divides_expr 
        | times "%" pow -> modulo_expr
        | pow

    ?plus: plus "+" times -> plus_expr         // priority 65
        | plus "-" times -> minus_expr | times

    ?compare: plus "=" plus -> eq_expr              // priority 50
        | plus "!=" plus -> not_eq_expr
        | plus ">" plus -> greater_expr
        | plus "<" plus -> less_expr
        | plus ">=" plus -> greater_eq_expr
        | plus "<=" plus -> less_eq_expr
        | plus

    ?expr: compare

    ?conditions: expr ("," expr)* -> conditions

    ?prove_action: "prove" expr -> prove_action
        | "prove" expr "for" conditions -> prove_with_condition_action

    ?define_action: "define" expr -> define_action
        | "define" expr "for" conditions -> define_with_condition_action

    ?calculate_action: "calculate" expr -> calculate_action
        | "calculate" expr "for" conditions -> calculate_with_condition_action

    ?subgoal_action: "subgoal" INT ":" expr -> subgoal_action
        | "subgoal" INT ":" expr "for" conditions -> subgoal_with_condition_action

    ?done_action: "done" -> done_action

    ?sorry_action: "sorry" -> sorry_action

    ?rewrite_goal_action: "from" INT ":" -> rewrite_goal_action

    ?induction_action: "induction" "on" CNAME -> induction_action
        | "induction" "on" CNAME "starting" "from" expr -> induction_starting_action

    ?case_analysis_action: "case" "analysis" "on" expr -> case_analysis_action

    ?lhs_action: "lhs" ":" -> lhs_action
    ?rhs_action: "rhs" ":" -> rhs_action
    ?arg_action: "arg" ":" -> arg_action
    ?base_case_action: "base" ":" -> base_case_action
    ?induct_case_action: "induct" ":" -> induct_case_action
    ?case_action: "case" "true" ":" -> case_true
        | "case" "false" ":" -> case_false
        | "case" "negative" ":" -> case_negative
        | "case" "zero" ":" -> case_zero
        | "case" "positive" ":" -> case_positive

    ?inst_equation: CNAME "for" expr -> inst_equation

    ?atomic_rule: "substitute" CNAME "for" expr -> substitute_rule
        | "substitute" expr "for" CNAME "creating" CNAME -> inverse_substitute_rule
        | "apply" "integral" "identity" -> integral_identity_rule
        | "apply" "indefinite" "integral" -> indefinite_integral_rule
        | "integrate" "by" "parts" "with" "u" "=" expr "," "v" "=" expr -> integrate_by_parts_rule
        | "split" "region" "at" expr -> split_region_rule
        | "rewrite" expr "to" expr -> equation_rule
        | "rewrite" "to" expr -> equation_none_rule
        | "expand" "polynomial" -> expand_polynomial_rule
        | "partial" "fraction" "decomposition" -> partial_fraction_decomposition_rule
        | "rewrite" expr "to" expr "using" "identity" -> apply_identity_rule
        | "solve" "integral" expr -> solve_integral_rule
        | "solve" "equation" "for" expr -> solve_equation_rule
        | "differentiate" "both" "sides" "at" CNAME -> deriv_equation_rule
        | "integrate" "both" "sides" -> integral_equation_rule
        | "apply" INT "on" expr -> apply_equation_rule
        | "apply" expr "on" expr -> apply_equation_expr_rule
        | "apply" "limit" CNAME "->" expr "both" "sides" -> apply_limit_rule
        | "expand" "definition" "for" CNAME -> expand_definition_rule
        | "fold" "definition" "for" CNAME -> fold_definition_rule
        | "exchange" "derivative" "and" "integral" -> exchange_deriv_int_rule
        | "substitute" inst_equation ("," inst_equation)* "in" "equation" -> inst_equation_rule
        | "apply" "series" "expansion" "on" expr "index" CNAME -> apply_series_expansion_rule
        | "apply" "series" "evaluation" -> apply_series_evaluation_rule
        | "exchange" "integral" "and" "sum" -> exchange_integral_sum_rule
        | "apply" "induction" "hypothesis" -> apply_induction_hypothesis_rule
        | "linearity" -> apply_linearity_rule
        | "improper" "integral" "to" "limit" "creating" CNAME -> elim_improper_integral_rule
        | "replace" "substitution" -> replace_substitution_rule
        | "l'Hopital's" "rule" -> lhopitals_rule
        | "simplify" -> full_simplify_rule

    ?rule: atomic_rule
        | atomic_rule "(at" INT ")" -> on_count_rule
        | atomic_rule "(all)" -> on_subterms_rule

    ?action: prove_action
        | subgoal_action
        | done_action
        | sorry_action
        | rewrite_goal_action
        | induction_action
        | case_analysis_action
        | define_action
        | calculate_action
        | lhs_action
        | rhs_action
        | arg_action
        | base_case_action
        | induct_case_action
        | case_action
        | rule -> rule_action

    %import common.CNAME
    %import common.WS
    %import common.INT
    %import common.DECIMAL

    %ignore WS
"""


@v_args(inline=True)
class ExprTransformer(Transformer):
    def __init__(self):
        pass

    def var_expr(self, s):
        return expr.Var(str(s))

    def symbol_expr(self, s):
        return expr.Symbol(str(s), [expr.VAR, expr.CONST, expr.OP, expr.FUN])

    def int_expr(self, n):
        return expr.Const(int(n))

    def decimal_expr(self, n):
        return expr.Const(Decimal(n))

    def plus_expr(self, a, b):
        return expr.Op("+", a, b)

    def minus_expr(self, a, b):
        return expr.Op("-", a, b)

    def times_expr(self, a, b):
        return expr.Op("*", a, b)

    def divides_expr(self, a, b):
        if a.ty == expr.CONST and b.ty == expr.CONST:
            return expr.Const(Fraction(a.val) / Fraction(b.val))
        else:
            return expr.Op("/", a, b)

    def modulo_expr(self, a, b):
        if a.ty == expr.CONST and b.ty == expr.CONST and isinstance(a.val, int) and isinstance(b.val, int):
            return expr.Const(a.val % b.val)
        else:
            return expr.Op("%", a, b)

    def pow_expr(self, a, b):
        return expr.Op("^", a, b)

    def eq_expr(self, a, b):
        return expr.Op("=", a, b)

    def not_eq_expr(self, a, b):
        return expr.Op("!=", a, b)

    def less_expr(self, a, b):
        return expr.Op("<", a, b)

    def greater_expr(self, a, b):
        return expr.Op(">", a, b)

    def less_eq_expr(self, a, b):
        return expr.Op("<=", a, b)

    def greater_eq_expr(self, a, b):
        return expr.Op(">=", a, b)

    def uminus_expr(self, a):
        if expr.is_const(a) and a.val > 0:
            return expr.Const(-a.val)
        else:
            return expr.Op("-", a)

    def uminus_pow_expr(self, a, b):
        return expr.Op("-", expr.Op("^", a, b))

    def pi_expr(self):
        return expr.pi

    def g_expr(self):
        return expr.G

    def pos_inf_expr(self):
        return expr.Inf(Decimal("inf"))

    def neg_inf_expr(self):
        return expr.Inf(Decimal("-inf"))

    def fun_expr(self, func_name, *args):
        if func_name == 'SKOLEM_CONST':
            return expr.SkolemFunc(str(args[0]), tuple())
        elif func_name == 'SKOLEM_FUNC':
            return expr.SkolemFunc(str(args[0].func_name), tuple(arg for arg in args[0].args))
        elif func_name == 'SUM':
            e = expr.Summation(str(args[0]), *args[1:])
            return e
        elif func_name == 'MUL':
            e = expr.Product(str(args[0]), *args[1:])
            return e
        return expr.Fun(str(func_name), *args)

    def abs_expr(self, expr):
        return expr.Fun("abs", expr)

    def deriv_expr(self, var, body):
        return expr.Deriv(str(var), body)

    def integral_expr(self, var, lower, upper, body):
        return expr.Integral(str(var), lower, upper, body)

    def indefinite_integral_expr(self, var, body):
        return expr.IndefiniteIntegral(str(var), body, tuple())

    def indefinite_integral_skolem_expr(self, *args):
        var = args[0]
        skolem_args = tuple(str(arg) for arg in args[1:-1])
        body = args[-1]
        return expr.IndefiniteIntegral(str(var), body, skolem_args)

    def eval_at_expr(self, body, var, lower, upper):
        return expr.EvalAt(var, lower, upper, body)

    def limit_inf_expr(self, var, lim, body):
        return expr.Limit(str(var), lim, body)

    def limit_l_expr(self, var, lim, body):
        return expr.Limit(str(var), lim, body, "-")

    def limit_r_expr(self, var, lim, body):
        return expr.Limit(str(var), lim, body, "+")

    def vector_expr(self, *args):
        data = []
        for arg in args:
            if isinstance(arg, expr.Matrix):
                data.append(arg.data)
            else:
                data.append(arg)
        return expr.Matrix(tuple(data))
    
    def conditions(self, *exprs: Expr) -> Tuple[Expr]:
        return tuple(exprs)

    def prove_action(self, expr: Expr):
        from integral import action
        return action.ProveAction(expr)

    def prove_with_condition_action(self, expr: Expr, conditions: Tuple[Expr]):
        from integral import action
        return action.ProveAction(expr, conditions)

    def define_action(self, expr: Expr):
        from integral import action
        return action.DefineAction(expr)
    
    def define_with_condition_action(self, expr: Expr, conditions: Tuple[Expr]):
        from integral import action
        return action.DefineAction(expr, conditions)

    def calculate_action(self, expr: Expr):
        from integral import action
        return action.CalculateAction(expr)
    
    def calculate_with_condition_action(self, expr: Expr, conditions: Tuple[Expr]):
        from integral import action
        return action.CalculateAction(expr, conditions)

    def subgoal_action(self, name: Token, expr: Expr):
        from integral import action
        return action.SubgoalAction(str(name), expr)
    
    def subgoal_with_condition_action(self, name: Token, expr: Expr, conditions: Tuple[Expr]):
        from integral import action
        return action.SubgoalAction(str(name), expr, conditions)
    
    def done_action(self):
        from integral import action
        return action.DoneAction()

    def sorry_action(self):
        from integral import action
        return action.SorryAction()
 
    def rewrite_goal_action(self, name: Token):
        from integral import action
        return action.RewriteGoalAction(str(name))
    
    def induction_action(self, var_name: Token):
        from integral import action
        return action.InductionAction(str(var_name), expr.Const(0))

    def induction_starting_action(self, var_name: Token, start: Expr):
        from integral import action
        return action.InductionAction(str(var_name), start)

    def case_analysis_action(self, split_cond: Expr):
        from integral import action
        return action.CaseAnalysisAction(split_cond)

    def lhs_action(self):
        from integral import action
        return action.LHSAction()

    def rhs_action(self):
        from integral import action
        return action.RHSAction()

    def arg_action(self):
        from integral import action
        return action.ArgAction()

    def base_case_action(self):
        from integral import action
        return action.BaseCaseAction()
    
    def induct_case_action(self):
        from integral import action
        return action.InductCaseAction()

    def case_true(self):
        from integral import action
        return action.CaseAction("true")

    def case_false(self):
        from integral import action
        return action.CaseAction("false")

    def case_negative(self):
        from integral import action
        return action.CaseAction("negative")

    def case_zero(self):
        from integral import action
        return action.CaseAction("zero")

    def case_positive(self):
        from integral import action
        return action.CaseAction("positive")

    def substitute_rule(self, var_name: Token, expr: Token):
        from integral import rules
        return rules.Substitution(str(var_name), str(expr))

    def inverse_substitute_rule(self, expr: Token, old_var: Token, var_name: Token):
        from integral import rules
        return rules.SubstitutionInverse(str(var_name), str(old_var), str(expr))

    def integral_identity_rule(self):
        from integral import rules
        return rules.DefiniteIntegralIdentity()

    def indefinite_integral_rule(self):
        from integral import rules
        return rules.IndefiniteIntegralIdentity()
    
    def integrate_by_parts_rule(self, u_expr: Expr, v_expr: Expr):
        from integral import rules
        return rules.IntegrationByParts(u_expr, v_expr)

    def split_region_rule(self, expr: Expr):
        from integral import rules
        return rules.SplitRegion(expr)
    
    def equation_rule(self, old_expr: Expr, new_expr: Expr):
        from integral import rules
        return rules.Equation(old_expr, new_expr)

    def equation_none_rule(self, new_expr: Expr):
        from integral import rules
        return rules.Equation(old_expr=None, new_expr=new_expr)

    def expand_polynomial_rule(self):
        from integral import rules
        return rules.ExpandPolynomial()
    
    def partial_fraction_decomposition_rule(self):
        from integral import rules
        return rules.PartialFractionDecomposition()

    def apply_identity_rule(self, old_expr: Expr, new_expr: Expr):
        from integral import rules
        return rules.ApplyIdentity(old_expr, new_expr)
    
    def solve_integral_rule(self, expr: Expr):
        from integral import rules
        return rules.IntegrateByEquation(expr)

    def solve_equation_rule(self, expr: Expr):
        from integral import rules
        return rules.SolveEquation(expr)

    def deriv_equation_rule(self, var: Token):
        from integral import rules
        return rules.DerivEquation(str(var))
    
    def integral_equation_rule(self):
        from integral import rules
        return rules.IntegralEquation()
    
    def apply_equation_rule(self, name: Token, source: Expr):
        from integral import rules
        return rules.ApplyEquation(str(name), source)
    
    def apply_equation_expr_rule(self, eq: Expr, source: Expr):
        from integral import rules
        return rules.ApplyEquation(eq, source)

    def apply_limit_rule(self, var_name: Token, limit: Expr):
        from integral import rules
        return rules.LimitEquation(str(var_name), limit)

    def expand_definition_rule(self, func_name: Token):
        from integral import rules
        return rules.ExpandDefinition(str(func_name))

    def fold_definition_rule(self, func_name):
        from integral import rules
        return rules.FoldDefinition(str(func_name))

    def exchange_deriv_int_rule(self):
        from integral import rules
        return rules.DerivIntExchange()
    
    def inst_equation(self, var_name: Token, expr: Expr):
        return {'var': str(var_name), 'expr': expr}
    
    def inst_equation_rule(self, *insts):
        from integral import rules
        return rules.VarSubsOfEquation(list(insts))

    def apply_series_expansion_rule(self, old_expr: Expr, index_var: Token):
        from integral import rules
        return rules.SeriesExpansionIdentity(old_expr=old_expr, index_var=str(index_var))

    def apply_linearity_rule(self):
        from integral import rules
        return rules.Linearity()
    def apply_series_evaluation_rule(self):
        from integral import rules
        return rules.SeriesEvaluationIdentity()

    def exchange_integral_sum_rule(self):
        from integral import rules
        return rules.IntSumExchange()

    def apply_induction_hypothesis_rule(self):
        from integral import rules
        return rules.ApplyInductHyp()

    def elim_improper_integral_rule(self, var_name: Token):
        from integral import rules
        return rules.ElimInfInterval(new_var=str(var_name))

    def replace_substitution_rule(self):
        from integral import rules
        return rules.ReplaceSubstitution()

    def full_simplify_rule(self):
        from integral import rules
        return rules.Simplify()
    
    def lhopitals_rule(self):
        from integral import rules
        return rules.LHopital()

    def on_count_rule(self, rule, n: Token):
        from integral import rules
        return rules.OnCount(rule, int(str(n)))
    
    def on_subterms_rule(self, rule):
        from integral import rules
        return rules.OnSubterm(rule)

    def rule_action(self, rule):
        from integral import action
        return action.RuleAction(rule)


transformer = ExprTransformer()
expr_parser = Lark(grammar, start="expr", parser="lalr", transformer=transformer)
action_parser = Lark(grammar, start="action", parser="lalr", transformer=transformer)


class ParseException(Exception):
    def __init__(self, s: str, msg: str):
        self.s = s
        self.msg = msg

    def __str__(self):
        return "Error while parsing %s\n%s" % (self.s, self.msg)


def parse_expr(s: str) -> Expr:
    """Parse an integral expression."""
    try:
        res = expr_parser.parse(s)
        return res
    except (exceptions.UnexpectedCharacters, exceptions.UnexpectedToken) as e:
        raise ParseException(s, str(e))

def parse_action(s: str):
    try:
        res = action_parser.parse(s)
        return res
    except (exceptions.UnexpectedCharacters, exceptions.UnexpectedToken) as e:
        raise ParseException(s, str(e))
