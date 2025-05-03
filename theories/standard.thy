# Standard integrals

imports base

prove (INT x. 1 / (x + a)) = log(abs(x + a)) + SKOLEM_CONST(C) for x != -a
lhs:
    substitute u for x + a
    apply integral identity
    replace substitution
done

prove (INT x. 1 / (x - a)) = log(abs(x - a)) + SKOLEM_CONST(C) for x != a
lhs:
    substitute u for x - a
    apply integral identity
    replace substitution
done

prove (INT x. 1 / (-x + a)) = -log(abs(-x + a)) + SKOLEM_CONST(C) for x != a
lhs:
    substitute u for -x + a
    apply integral identity
    replace substitution
done

prove (INT x. 1 / (a * x + b)) = log(abs(a * x + b)) / a + SKOLEM_CONST(C) for a != 0, a * x + b != 0
lhs:
    substitute u for a * x + b
    apply integral identity
    replace substitution
done

prove (INT x. 1 / (a * x)) = log(abs(a * x)) / a + SKOLEM_CONST(C) for a != 0, x != 0
lhs:
    substitute u for a * x
    apply integral identity
    replace substitution
done

prove (INT x. 1 / (a * x - b)) = log(abs(a * x - b)) / a + SKOLEM_CONST(C) for a != 0, a * x - b != 0
lhs:
    substitute u for a * x - b
    apply integral identity
    replace substitution
done

prove (INT x. 1 / (-(a * x) + b)) = log(abs(-(a * x) + b)) / (-a) + SKOLEM_CONST(C) for a != 0, -(a * x) + b != 0
lhs:
    substitute u for -(a * x) + b
    apply integral identity
    replace substitution
done

prove (INT x. 1 / (-(a * x) - b)) = log(abs(-(a * x) - b)) / (-a) + SKOLEM_CONST(C) for a != 0, -(a * x) - b != 0
lhs:
    substitute u for -(a * x) - b
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

prove (INT x. exp(a * x + b)) = exp(a * x + b) / a + SKOLEM_CONST(C) for a != 0
lhs:
    substitute u for a * x + b
    apply integral identity
    replace substitution
    rewrite 1/a * exp(a * x + b) to exp(a * x + b) / a
done

prove (INT x. exp(a * x - b)) = exp(a * x - b) / a + SKOLEM_CONST(C) for a != 0
lhs:
    rewrite exp(a * x - b) to exp(-b) * exp(a * x)
    substitute u for a * x
    simplify
    rewrite exp(-b + u) to exp(-b) * exp(u)
    apply integral identity
    replace substitution
done

prove (INT x. exp(-(a * x))) = - exp(-(a * x)) / a + SKOLEM_CONST(C) for a != 0
lhs:
    substitute u for -(a * x)
    apply integral identity
    replace substitution
done

prove (INT x. exp(-x)) = - exp(-x) + SKOLEM_CONST(C)
lhs:
    substitute u for -x
    apply integral identity
    replace substitution
done

prove (INT x. sin(a * x)) = -(cos(a * x) / a) + SKOLEM_CONST(C) for a != 0
lhs:
    substitute u for a * x
    apply integral identity
    replace substitution
    rewrite 1 / a * -cos(a * x) to -cos(a * x) / a
done

prove (INT x. sin(-x)) = cos(-x) + SKOLEM_CONST(C)
lhs:
    substitute u for -x
    apply integral identity
    replace substitution
    simplify
done

prove (INT x. sin(a * x + b)) = -cos(a * x + b) / a + SKOLEM_CONST(C) for a != 0
lhs:
    substitute u for a * x + b
    apply integral identity
    replace substitution
    rewrite 1 / a * -cos(a * x + b) to -cos(a * x + b) / a
done

prove (INT x. cos(a * x)) = sin(a * x) / a + SKOLEM_CONST(C) for a != 0
lhs:
    substitute u for a * x
    apply integral identity
    replace substitution
done

prove (INT x. cos(a * x + b)) = sin(a * x + b) / a + SKOLEM_CONST(C) for a != 0
lhs:
    substitute u for a * x + b
    apply integral identity
    replace substitution
done

prove (INT x. cos(-x)) = -sin(-x) + SKOLEM_CONST(C)
lhs:
    rewrite cos(-x) to cos(x)
    apply integral identity
    rewrite sin(x) to -sin(-x)
done

prove (INT x. tan(x)) = -log(abs(cos(x))) + SKOLEM_CONST(C)
sorry

prove (INT x. tan(-x)) = log(abs(cos(x))) + SKOLEM_CONST(C)
sorry

prove (INT x. cot(x)) = -log(abs(sin(-x))) + SKOLEM_CONST(C)
sorry

prove (INT x. tan(a * x)) = log(abs(sec(a*x))) / a + SKOLEM_CONST(C) for a != 0
sorry

prove (INT x. tan(a * x + b)) = log(abs(sec(a*x + b))) / a + SKOLEM_CONST(C) for a != 0
sorry

prove (INT x. tan(a * x - b)) = log(abs(sec(a*x - b))) / a + SKOLEM_CONST(C) for a != 0
sorry

prove (INT x. 1 / (x - a)) = log(abs(x - a)) + SKOLEM_CONST(C) for x - a != 0
lhs:
    substitute u for x - a
    apply integral identity
    replace substitution
done

prove (INT x. 1 / (x + a)) = log(abs(x + a)) + SKOLEM_CONST(C) for x + a != 0
lhs:
    substitute u for x + a
    apply integral identity
    replace substitution
done

prove (INT x. 1 / (b - a * x)) = -log(abs(b - a * x)) / a + SKOLEM_CONST(C) for a != 0, b - a * x != 0
lhs:
    substitute u for b - a * x
    apply integral identity
    replace substitution
done

prove (INT x. 1 / (b + a * x)) = log(abs(b + a * x)) / a + SKOLEM_CONST(C) for a != 0, b + a * x != 0
lhs:
    substitute u for b + a * x
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
done

prove (INT x. 1 / (a + x ^ 2)) = (1 / sqrt(a)) * arctan(x / sqrt(a)) + SKOLEM_CONST(C) for a > 0
lhs:
    rewrite 1/(a + x^2) to (1/a)/(1 + (x/sqrt(a))^2)
    substitute u for x/sqrt(a)
    apply integral identity
    replace substitution
done

prove (INT x. 1 / (a + b * x ^ 2)) = (1 / sqrt(a*b)) * arctan(sqrt(a/b)*x) + SKOLEM_CONST(C) for a > 0, b > 0
sorry

prove (INT x. 1 / (b * x ^ 2 + a)) = (1 / sqrt(a*b)) * arctan(sqrt(a/b)*x) + SKOLEM_CONST(C) for a > 0, b > 0
sorry

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

prove (INT x. sec(x)) = log(abs(sec(x)+tan(x))) + SKOLEM_CONST(C)
sorry

prove (INT x. sec(-x)) = -log(abs(sec(-x)+tan(-x))) + SKOLEM_CONST(C)
sorry

prove (INT x. sec(a*x)) = log(abs(sec(a*x)+tan(a*x))) / a + SKOLEM_CONST(C) for a != 0
sorry

prove (INT x. sec(a*x + b)) = log(abs(sec(a*x+b)+tan(a*x+b))) / a + SKOLEM_CONST(C) for a != 0
sorry

prove (INT x. sec(a*x - b)) = log(abs(sec(a*x-b)+tan(a*x-b))) / a + SKOLEM_CONST(C) for a != 0
sorry

prove (INT x. 1/sqrt(1-x^2)) = arcsin(x) + SKOLEM_CONST(C)
sorry

prove (INT x. 1/sqrt(-(x^2)+1)) = arcsin(x) + SKOLEM_CONST(C)
sorry

prove (INT x. 1/sqrt(a-x^2)) = arcsin(x/sqrt(a)) + SKOLEM_CONST(C) for a > 0
sorry

prove (INT x. 1/sqrt(-(x^2)+a)) = arcsin(x/sqrt(a))+ SKOLEM_CONST(C) for a > 0
sorry

prove (INT x. 1/sqrt(x^2-1)) = arccos(x) + SKOLEM_CONST(C)
sorry

prove (INT x. 1/sqrt(-1+x^2)) = arccos(x) + SKOLEM_CONST(C)
sorry

prove (INT x. 1/sqrt(x^2-a)) = arccos(x/sqrt(a)) + SKOLEM_CONST(C) for a > 0
sorry

prove (INT x. 1/sqrt(-a+x^2)) = arccos(x/sqrt(a)) + SKOLEM_CONST(C) for a > 0
sorry
