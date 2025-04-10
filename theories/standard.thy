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