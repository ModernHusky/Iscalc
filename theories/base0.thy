# Common integrals

axiom (INT x. c) = c * x + SKOLEM_CONST(C)

axiom (INT x. x) = x ^ 2 / 2 + SKOLEM_CONST(C)

axiom (INT x. -x) = -x ^ 2 / 2 + SKOLEM_CONST(C)

axiom (INT x. a * x) = x ^ 2 / a + SKOLEM_CONST(C) for a != 0

axiom (INT x. a*x + b) = x ^ 2 / a + b * x + SKOLEM_CONST(C) for a != 0

axiom (INT x. a*x - b) = x ^ 2 / a - b * x + SKOLEM_CONST(C) for a != 0

axiom (INT x. 1 / x) = log(abs(x)) + SKOLEM_CONST(C)

axiom (INT x. 1 / (-x)) = -log(abs(x)) + SKOLEM_CONST(C)

axiom (INT x. x ^ n) = x ^ (n + 1) / (n + 1) + SKOLEM_CONST(C) for n != -1

axiom (INT x. 1 / x ^ n) = 1 / (-((n - 1) * (x ^ (n - 1)))) + SKOLEM_CONST(C) for n != 1

axiom (INT x. sqrt(x)) = 2/3 * x ^ (3/2) + SKOLEM_CONST(C)

axiom (INT x. 1 / sqrt(x)) = 2 * sqrt(x) + SKOLEM_CONST(C)

axiom (INT x. exp(x)) = exp(x) + SKOLEM_CONST(C)

axiom (INT x. exp(-x)) = -exp(-x) + SKOLEM_CONST(C)

axiom (INT x. sin(x)) = -cos(x) + SKOLEM_CONST(C)

axiom (INT x. cos(x)) = sin(x) + SKOLEM_CONST(C)

axiom (INT x. 1 / (x^2 + 1)) = arctan(x) + SKOLEM_CONST(C)

axiom (INT x. sec(x)^2) = tan(x) + SKOLEM_CONST(C)

axiom (INT x. 1 / cos(x)^2) = tan(x) + SKOLEM_CONST(C)

axiom (INT x. csc(x)^2) = -cot(x) + SKOLEM_CONST(C)

axiom (INT x. 1 / sin(x)^2) = -cot(x) + SKOLEM_CONST(C)

axiom (INT x. cot(x) * csc(x)) = -csc(x) + SKOLEM_CONST(C)

axiom (INT x. sec(x) * tan(x)) = sec(x) + SKOLEM_CONST(C)

axiom (INT x. sec(x)) = log(abs(sec(x) + tan(x))) + SKOLEM_CONST(C)

axiom (INT x. 1 / cos(x)) = log(abs(sec(x) + tan(x))) + SKOLEM_CONST(C)

axiom (INT x. csc(x)) = log(abs(csc(x) - cot(x))) + SKOLEM_CONST(C)

axiom (INT x. 1 / sin(x)) = log(abs(csc(x) - cot(x))) + SKOLEM_CONST(C)
