imports base

# Lemmas for integrals involving square roots
# These are needed for inverse trigonometric function integrals in standard4

## Basic sqrt integrals

# From standard4.thy line 12-25 pattern
prove (INT x. x / sqrt(a ^ 2 - x ^ 2)) = -sqrt(a ^ 2 - x ^ 2) + SKOLEM_CONST(C) for a ^ 2 - x ^ 2 >= 0, sqrt(a ^ 2 - x ^ 2) != 0
lhs:
    substitute u for a ^ 2 - x ^ 2
    simplify
    rewrite 1 / sqrt(u) to u ^ (-1/2)
    apply integral identity
    replace substitution
    simplify
    rewrite SKOLEM_CONST(C) - sqrt(a ^ 2 - x ^ 2) to -sqrt(a ^ 2 - x ^ 2) + SKOLEM_CONST(C)
done

# Helper lemma: INT sqrt(a^2-x^2) - trigonometric substitution
prove (INT x. sqrt(a ^ 2 - x ^ 2)) = x / 2 * sqrt(a ^ 2 - x ^ 2) + a ^ 2 / 2 * arcsin(x / a) + SKOLEM_CONST(C) for abs(x / a) <= 1, a != 0, a ^ 2 - x ^ 2 >= 0, a > 0
lhs:
    substitute u for arcsin(x / a)
    simplify
    rewrite a ^ 2 - a ^ 2 * sin(u) ^ 2 to a ^ 2 * (1 - sin(u) ^ 2)
    rewrite 1 - sin(u) ^ 2 to cos(u) ^ 2
    rewrite sqrt(a ^ 2 * cos(u) ^ 2) to abs(a) * abs(cos(u))
    simplify
    rewrite cos(u) ^ 2 to (1 + cos(2 * u)) / 2
    simplify
    apply integral identity
    simplify
    substitute v for 2 * u
    simplify
    apply integral identity
    replace substitution
    simplify
    rewrite 1 - x ^ 2 / a ^ 2 to (a ^ 2 - x ^ 2) / a ^ 2
    rewrite sqrt((a ^ 2 - x ^ 2) / a ^ 2) to sqrt(a ^ 2 - x ^ 2) / abs(a)
    simplify
    rewrite a ^ 2 * (x * sqrt(a ^ 2 - x ^ 2) / (2 * a ^ 2) + 1/2 * arcsin(x / a)) to x * sqrt(a ^ 2 - x ^ 2) / 2 + a ^ 2 / 2 * arcsin(x / a)
done

# Needed for standard4.thy lines 53-58 (x^2 * arcsin, x^3 * arcsin)
prove (INT x. x ^ 2 / sqrt(a ^ 2 - x ^ 2)) = -x / 2 * sqrt(a ^ 2 - x ^ 2) + a ^ 2 / 2 * arcsin(x / a) + SKOLEM_CONST(C) for abs(x / a) < 1, a != 0, a ^ 2 - x ^ 2 >= 0, sqrt(a ^ 2 - x ^ 2) != 0, a > 0
lhs:
    integrate by parts with u = x, v = -sqrt(a ^ 2 - x ^ 2)
    simplify
    apply integral identity
    simplify
    rewrite a ^ 2 / 2 * arcsin(x / a) - x * sqrt(a ^ 2 - x ^ 2) / 2 + SKOLEM_CONST(C) to -x / 2 * sqrt(a ^ 2 - x ^ 2) + a ^ 2 / 2 * arcsin(x / a) + SKOLEM_CONST(C)
done

prove (INT x. x ^ 3 / sqrt(a ^ 2 - x ^ 2)) = -(x ^ 2 + 2 * a ^ 2) / 3 * sqrt(a ^ 2 - x ^ 2) + SKOLEM_CONST(C) for a ^ 2 - x ^ 2 >= 0, sqrt(a ^ 2 - x ^ 2) != 0
lhs:
    integrate by parts with u = x ^ 2, v = -sqrt(a ^ 2 - x ^ 2)
    simplify
    rewrite x * sqrt(a ^ 2 - x ^ 2) to x * (a ^ 2 - x ^ 2) / sqrt(a ^ 2 - x ^ 2)
    substitute w for a ^ 2 - x ^ 2
    simplify
    rewrite w / sqrt(w) to sqrt(w)
    rewrite sqrt(w) to w ^ (1/2)
    apply integral identity
    replace substitution
    simplify
    rewrite (a ^ 2 - x ^ 2) ^ (3/2) to (a ^ 2 - x ^ 2) * (a ^ 2 - x ^ 2) ^ (1/2)
    rewrite (a ^ 2 - x ^ 2) ^ (1/2) to sqrt(a ^ 2 - x ^ 2)
    rewrite -(2 * ((a ^ 2 - x ^ 2) * sqrt(a ^ 2 - x ^ 2)) / 3) - x ^ 2 * sqrt(a ^ 2 - x ^ 2) + SKOLEM_CONST(C) to -(2 * (a ^ 2 - x ^ 2) * sqrt(a ^ 2 - x ^ 2) / 3 + x ^ 2 * sqrt(a ^ 2 - x ^ 2)) + SKOLEM_CONST(C)
    rewrite 2 * (a ^ 2 - x ^ 2) * sqrt(a ^ 2 - x ^ 2) / 3 + x ^ 2 * sqrt(a ^ 2 - x ^ 2) to sqrt(a ^ 2 - x ^ 2) * (2 * (a ^ 2 - x ^ 2) / 3 + x ^ 2)
    expand polynomial
    simplify
    rewrite -(2 * a ^ 2 * sqrt(a ^ 2 - x ^ 2) / 3) - x ^ 2 * sqrt(a ^ 2 - x ^ 2) / 3 + SKOLEM_CONST(C) to -((2 * a ^ 2 + x ^ 2) / 3 * sqrt(a ^ 2 - x ^ 2)) + SKOLEM_CONST(C)
    rewrite (2 * a ^ 2 + x ^ 2) / 3 to (x ^ 2 + 2 * a ^ 2) / 3
    rewrite -((x ^ 2 + 2 * a ^ 2) / 3 * sqrt(a ^ 2 - x ^ 2)) + SKOLEM_CONST(C) to -(x ^ 2 + 2 * a ^ 2) / 3 * sqrt(a ^ 2 - x ^ 2) + SKOLEM_CONST(C)
done

# Helper lemma: INT (a^2-x^2)^(3/2) - needed for INT x^2*sqrt
prove (INT x. (a ^ 2 - x ^ 2) ^ (3/2)) = x / 4 * (a ^ 2 - x ^ 2) * sqrt(a ^ 2 - x ^ 2) + 3 * a ^ 2 / 8 * (x * sqrt(a ^ 2 - x ^ 2) + a ^ 2 * arcsin(x / a)) + SKOLEM_CONST(C) for abs(x / a) <= 1, a != 0, a > 0, a ^ 2 - x ^ 2 >= 0, sqrt(a ^ 2 - x ^ 2) != 0
sorry

# Helper: INT x^2*sqrt(a^2-x^2) - now provable using algebraic decomposition
prove (INT x. x ^ 2 * sqrt(a ^ 2 - x ^ 2)) = x ^ 3 / 4 * sqrt(a ^ 2 - x ^ 2) - a ^ 2 * x / 8 * sqrt(a ^ 2 - x ^ 2) + a ^ 4 / 8 * arcsin(x / a) + SKOLEM_CONST(C) for abs(x / a) <= 1, a != 0, a > 0, a ^ 2 - x ^ 2 >= 0, sqrt(a ^ 2 - x ^ 2) != 0
lhs:
    rewrite x ^ 2 * sqrt(a ^ 2 - x ^ 2) to (a ^ 2 - (a ^ 2 - x ^ 2)) * sqrt(a ^ 2 - x ^ 2)
    rewrite (a ^ 2 - (a ^ 2 - x ^ 2)) * sqrt(a ^ 2 - x ^ 2) to a ^ 2 * sqrt(a ^ 2 - x ^ 2) - (a ^ 2 - x ^ 2) * sqrt(a ^ 2 - x ^ 2)
    rewrite (a ^ 2 - x ^ 2) * sqrt(a ^ 2 - x ^ 2) to (a ^ 2 - x ^ 2) ^ (3/2)
    simplify
    apply integral identity
    simplify
    rewrite (a ^ 2 - x ^ 2) ^ (3/2) to (a ^ 2 - x ^ 2) * sqrt(a ^ 2 - x ^ 2)
    expand polynomial
    simplify
done

prove (INT x. x ^ 4 / sqrt(a ^ 2 - x ^ 2)) = -x ^ 3 / 4 * sqrt(a ^ 2 - x ^ 2) - 3 / 8 * a ^ 2 * x * sqrt(a ^ 2 - x ^ 2) + 3 / 8 * a ^ 4 * arcsin(x / a) + SKOLEM_CONST(C) for abs(x / a) <= 1, a != 0, a > 0, a ^ 2 - x ^ 2 >= 0, sqrt(a ^ 2 - x ^ 2) != 0
lhs:
    integrate by parts with u = x ^ 3, v = -sqrt(a ^ 2 - x ^ 2)
    simplify
    apply integral identity
    simplify
    rewrite -(3 * a ^ 2 * x * sqrt(a ^ 2 - x ^ 2) / 8) - x ^ 3 * sqrt(a ^ 2 - x ^ 2) / 4 + 3 * a ^ 4 / 8 * arcsin(x / a) to -x ^ 3 / 4 * sqrt(a ^ 2 - x ^ 2) - 3 / 8 * a ^ 2 * x * sqrt(a ^ 2 - x ^ 2) + 3 / 8 * a ^ 4 * arcsin(x / a)
done

# General formula for x^n / sqrt(a^2 - x^2)
prove (INT x. x ^ n / sqrt(a ^ 2 - x ^ 2)) = -x ^ (n - 1) * sqrt(a ^ 2 - x ^ 2) / n + (n - 1) * a ^ 2 / n * (INT x. x ^ (n - 2) / sqrt(a ^ 2 - x ^ 2)) for n >= 2, isInt(n), a ^ 2 - x ^ 2 > 0
sorry

## Integrals with x in denominator and sqrt

# Needed for standard4.thy lines 68-73 (1/x^2 * arcsin, 1/x^3 * arcsin)
prove (INT x. 1 / (x * sqrt(a ^ 2 - x ^ 2))) = -1 / a * log(abs((a + sqrt(a ^ 2 - x ^ 2)) / x)) + SKOLEM_CONST(C) for x != 0, a != 0, a ^ 2 - x ^ 2 > 0
sorry

prove (INT x. 1 / (x ^ 2 * sqrt(a ^ 2 - x ^ 2))) = -sqrt(a ^ 2 - x ^ 2) / (a ^ 2 * x) + SKOLEM_CONST(C) for x != 0, a != 0, a ^ 2 - x ^ 2 > 0
sorry

# General formula
prove (INT x. 1 / (x ^ n * sqrt(a ^ 2 - x ^ 2))) = -sqrt(a ^ 2 - x ^ 2) / ((n - 1) * a ^ 2 * x ^ (n - 1)) - (n - 2) / ((n - 1) * a ^ 2) * (INT x. 1 / (x ^ (n - 2) * sqrt(a ^ 2 - x ^ 2))) for n >= 2, x != 0, a != 0, a ^ 2 - x ^ 2 > 0
sorry
