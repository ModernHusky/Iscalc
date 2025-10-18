imports standard

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

# Helper lemma: INT (a^2-x^2)^(3/2)
# Use trigonometric substitution to break circular dependency
prove (INT x. (a ^ 2 - x ^ 2) ^ (3/2)) = x / 4 * (a ^ 2 - x ^ 2) * sqrt(a ^ 2 - x ^ 2) + 3 * a ^ 2 / 8 * (x * sqrt(a ^ 2 - x ^ 2) + a ^ 2 * arcsin(x / a)) + SKOLEM_CONST(C) for abs(x / a) <= 1, a > 0, a ^ 2 - x ^ 2 >= 0, sqrt(a ^ 2 - x ^ 2) != 0
lhs:
    substitute u for arcsin(x / a)
    simplify
    rewrite a ^ 2 - a ^ 2 * sin(u) ^ 2 to a ^ 2 * (1 - sin(u) ^ 2)
    rewrite 1 - sin(u) ^ 2 to cos(u) ^ 2
    rewrite (a ^ 2 * cos(u) ^ 2) ^ (3/2) to a ^ 3 * (cos(u) ^ 2) ^ (3/2)
    rewrite (cos(u) ^ 2) ^ (3/2) to abs(cos(u)) ^ 3
    rewrite abs(cos(u)) ^ 3 to abs(cos(u) ^ 3)
    simplify
    rewrite cos(u)^4 to cos(u)^2^2
    rewrite cos(u)^2 to (1+cos(2*u))/2
    expand polynomial
    rewrite cos(2 * u) ^ 2 to (1+cos(4*u))/2
    apply integral identity
    simplify
    replace substitution
    expand polynomial
    simplify
    rewrite sqrt(1 - x ^ 2 / a ^ 2) to sqrt((a^2 - x^2)/a^2)
    rewrite sqrt((a^2 - x^2)/a^2) to sqrt(a^2 - x^2) / a
    simplify
    expand polynomial
    simplify
    rewrite sqrt(1 - x ^ 2 / a ^ 2) to sqrt((a ^ 2 - x ^ 2) / a ^ 2)
    rewrite sqrt((a ^ 2 - x ^ 2) / a ^ 2) to sqrt(a ^ 2 - x ^ 2) / a
    simplify
rhs:
    expand polynomial
done


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

prove (INT x. x ^ 4 / sqrt(a ^ 2 - x ^ 2)) = -x ^ 3 / 4 * sqrt(a ^ 2 - x ^ 2) - 3 / 8 * a ^ 2 * x * sqrt(a ^ 2 - x ^ 2) + 3 / 8 * a ^ 4 * arcsin(x / a) + SKOLEM_CONST(C) for abs(x / a) <= 1, a > 0, a ^ 2 - x ^ 2 >= 0, sqrt(a ^ 2 - x ^ 2) != 0
lhs:
    integrate by parts with u = x ^ 3, v = -sqrt(a ^ 2 - x ^ 2)
    simplify
    apply integral identity
    simplify
    rewrite -(3 * a ^ 2 * x * sqrt(a ^ 2 - x ^ 2) / 8) - x ^ 3 * sqrt(a ^ 2 - x ^ 2) / 4 + 3 * a ^ 4 / 8 * arcsin(x / a) to -x ^ 3 / 4 * sqrt(a ^ 2 - x ^ 2) - 3 / 8 * a ^ 2 * x * sqrt(a ^ 2 - x ^ 2) + 3 / 8 * a ^ 4 * arcsin(x / a)
done

# General formula for x^n / sqrt(a^2 - x^2)
prove (INT x. x ^ n / sqrt(a ^ 2 - x ^ 2)) = -x ^ (n - 1) * sqrt(a ^ 2 - x ^ 2) / n + (n - 1) * a ^ 2 / n * (INT x. x ^ (n - 2) / sqrt(a ^ 2 - x ^ 2)) for n >= 2, isInt(n), x / a < 1, x / a > 0, a > 0, a ^ 2 - x ^ 2 > 0
lhs:
    rewrite x ^ n / sqrt(a ^ 2 - x ^ 2) to x ^ (n - 1) * x / sqrt(a ^ 2 - x ^ 2)
    integrate by parts with u = x ^ (n - 1), v = -sqrt(a ^ 2 - x ^ 2)
    simplify
    rewrite x ^ (n - 2) * sqrt(a ^ 2 - x ^ 2) to x ^ (n - 2) * (a ^ 2 - x ^ 2) / sqrt(a ^ 2 - x ^ 2)
    expand polynomial
    simplify
    rewrite x ^ (n - 2) * x ^ 2 to x ^ n
    solve integral INT x. x ^ n / sqrt(a ^ 2 - x ^ 2)
    expand polynomial
    rewrite x ^ (n - 2) * x ^ 2 / sqrt(a ^ 2 - x ^ 2) to x ^ n / sqrt(a ^ 2 - x ^ 2)
    solve integral INT x. x ^ n / sqrt(a ^ 2 - x ^ 2)
    expand polynomial
    simplify
    rewrite (n + 1) * (1 - 1 / (n + 1)) to n
    rewrite a ^ 2 / ((n + 1) * (1 - 1 / (n + 1))) to a ^ 2 / n
    rewrite x ^ (n - 1) * sqrt(a ^ 2 - x ^ 2) / ((n + 1) * (1 - 1 / (n + 1))) to x ^ (n - 1) * sqrt(a ^ 2 - x ^ 2) / n
    simplify
    rewrite a ^ 2 * (INT x. x ^ (n - 2) / sqrt(a ^ 2 - x ^ 2)) - a ^ 2 / n * (INT x. x ^ (n - 2) / sqrt(a ^ 2 - x ^ 2)) to (a ^ 2 - a ^ 2 / n) * (INT x. x ^ (n - 2) / sqrt(a ^ 2 - x ^ 2))
    rewrite a ^ 2 - a ^ 2 / n to a ^ 2 * (n - 1) / n
    rewrite a ^ 2 * (n - 1) / n to (n - 1) * a ^ 2 / n
done

## Integrals with x in denominator and sqrt

# Needed for standard4.thy lines 68-73 (1/x^2 * arcsin, 1/x^3 * arcsin)
prove (INT x. 1 / (x * sqrt(a ^ 2 - x ^ 2))) = -1 / a * log((a + sqrt(a ^ 2 - x ^ 2)) / x) + SKOLEM_CONST(C) for x != 0, a > 0, a ^ 2 - x ^ 2 > 0, x / a > 0, x / a < 1, (a + sqrt(a ^ 2 - x ^ 2)) / x > 0
lhs:
    substitute u for x / a
    rewrite a ^ 2 - a ^ 2 * u ^ 2 to a ^ 2 * (1 - u ^ 2)
    rewrite sqrt(a ^ 2 * (1 - u ^ 2)) to a * sqrt(1 - u ^ 2)
    simplify
    substitute sin(w) for u
    simplify
    rewrite 1 - sin(w) ^ 2 to cos(w) ^ 2
    rewrite sqrt(cos(w) ^ 2) to abs(cos(w))
    simplify
    apply integral identity
    simplify
    rewrite log((cos(w) + 1) / (1 - cos(w))) to log((cos(w) + 1) ^ 2 / ((cos(w) + 1) * (1 - cos(w))))
    rewrite (cos(w) + 1) * (1 - cos(w)) to 1 - cos(w) ^ 2
    rewrite 1 - cos(w) ^ 2 to sin(w) ^ 2
    replace substitution
    simplify
    rewrite 1 - x ^ 2 / a ^ 2 to (a ^ 2 - x ^ 2) / a ^ 2
    rewrite sqrt((a ^ 2 - x ^ 2) / a ^ 2) to sqrt(a ^ 2 - x ^ 2) / a
    rewrite (sqrt(a ^ 2 - x ^ 2) / a + 1) ^ 2 to (1 + sqrt(a ^ 2 - x ^ 2) / a) ^ 2
    rewrite (1 + sqrt(a ^ 2 - x ^ 2) / a) ^ 2 to ((a + sqrt(a ^ 2 - x ^ 2)) / a) ^ 2
    rewrite a ^ 2 / x ^ 2 * ((a + sqrt(a ^ 2 - x ^ 2)) / a) ^ 2 to ((a + sqrt(a ^ 2 - x ^ 2)) / x) ^ 2
    simplify
rhs:
    rewrite -1 / a * log((a + sqrt(a ^ 2 - x ^ 2)) / x) + SKOLEM_CONST(C) to SKOLEM_CONST(C) - 1 / a * log((a + sqrt(a ^ 2 - x ^ 2)) / x)
    rewrite (a + sqrt(a ^ 2 - x ^ 2)) / x to (sqrt(a ^ 2 - x ^ 2) + a) / x
done

# BLOCKED: abs(cos(u)) simplification issue
# After substitution u = arcsin(x/a), we get abs(cos(u)) which should simplify to cos(u)
# because u is in the range [-pi/2, pi/2] where cos is non-negative.
# However, the system doesn't propagate this information through the substitution context.
# Needs: condition propagation through substitutions or context-aware simplification
prove (INT x. 1 / (x ^ 2 * sqrt(a ^ 2 - x ^ 2))) = -sqrt(a ^ 2 - x ^ 2) / (a ^ 2 * x) + SKOLEM_CONST(C) for x != 0, a > 0, abs(x / a) < 1, a ^ 2 - x ^ 2 > 0
subgoal 1: sin(arcsin(x / a)) != 0
lhs:
    rewrite sin(arcsin(x / a)) to x / a
done
lhs:
    substitute u for arcsin(x / a)
    simplify
    rewrite a ^ 2 - a ^ 2 * sin(u) ^ 2 to a ^ 2 * (1 - sin(u) ^ 2)
    rewrite 1 - sin(u) ^ 2 to cos(u) ^ 2
    rewrite sqrt(a ^ 2 * cos(u) ^ 2) to abs(a) * abs(cos(u))
    simplify
    rewrite 1 / sin(u) ^ 2 to csc(u) ^ 2
    apply integral identity
    simplify
    replace substitution
    simplify
    rewrite 1 - x ^ 2 / a ^ 2 to (a ^ 2 - x ^ 2) / a ^ 2
    rewrite sqrt((a ^ 2 - x ^ 2) / a ^ 2) to sqrt(a ^ 2 - x ^ 2) / a
    simplify
    rewrite SKOLEM_CONST(C) - sqrt(a ^ 2 - x ^ 2) / (a ^ 2 * x) to -sqrt(a ^ 2 - x ^ 2) / (a ^ 2 * x) + SKOLEM_CONST(C)
done

# General formula
prove (INT x. 1 / (x ^ n * sqrt(a ^ 2 - x ^ 2))) = -sqrt(a ^ 2 - x ^ 2) / ((n - 1) * a ^ 2 * x ^ (n - 1)) - (n - 2) / ((n - 1) * a ^ 2) * (INT x. 1 / (x ^ (n - 2) * sqrt(a ^ 2 - x ^ 2))) for n >= 2, x != 0, a != 0, a ^ 2 - x ^ 2 > 0
sorry
