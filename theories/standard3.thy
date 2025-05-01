## 2.7 Logarithms and Inverse-Hyperbolic Functions

### 2.7.2 - 2.7.3 Combinations of logarithms and algebraic functions
// Page 238

### 2.7.2.5

//1
prove (INT x. (a + b * x) ^ m * log(x)) = 1/((m + 1) * b) * ((a + b * x) ^ (m + 1) * log(x) - (INT x. (a + b * x) ^ (m + 1)/x)) + SKOLEM_CONST(C)
sorry

// 2
prove (INT x. (a + b * x) ^ m * log(x)) = 1/((m + 1) * b) * ((a + b * x) ^ (m + 1) - a ^ (m + 1)) * log(x) - SUM(k, 0, m, binom(m, k) * a ^ (m - k) * b ^ k * x ^ (k + 1)/(k + 1) ^ 2) + SKOLEM_CONST(C)
sorry

### 2.7.2.6
// 1
prove (INT x. (a + b * x) * log(x)) = ((a + b * x) ^ 2/(2 * b) - a ^ 2/(2 * b)) * log(x) - (a * x + 1/4 * b * x ^ 2) + SKOLEM_CONST(C) for x > 0, b != 0
lhs:
    integrate by parts with u = log(x), v = a*x + (b*x^2)/2
    simplify
    rewrite (1/x)*(b*x^2/2 + a*x) to (b*x)/2 + a
    apply integral identity
    simplify
    rewrite log(x) * (b * x ^ 2 / 2 + a * x) to (b * x ^ 2 / 2 + a * x) * log(x)
    rewrite b * x ^ 2 / 2 + a * x to (b * x ^ 2 + 2 * a * x)/2
    rewrite (b * x ^ 2 + 2 * a * x)/2 to (x * (b * x + 2 * a))/2
    rewrite (x * (b * x + 2 * a)) / 2 to (a + b * x) ^ 2 / (2 * b) - a ^ 2 / (2 * b)
    simplify
done

// 2
prove (INT x. (a + b * x) ^ 2 * log(x)) = 1/(3 * b) * ((a + b * x) ^ 3 - a ^ 3) * log(x) - (a ^ 2 * x + a * b * x ^ 2/2 + b ^ 2 * x ^ 3/9 ) + SKOLEM_CONST(C) for x > 0, b != 0
lhs:
    integrate by parts with u = log(x), v = (a + b * x)^3 / (3 * b)
    simplify
    rewrite (b * x + a)^3 / x to (b^3 * x^3 + 3 * a * b^2 * x^2 + 3 * a^2 * b * x + a^3) / x
    rewrite (b^3 * x^3 + 3 * a * b^2 * x^2 + 3 * a^2 * b * x + a^3) / x to b^3 * x^2 + 3 * a * b^2 * x + 3 * a^2 * b + a^3 / x
    apply integral identity
    simplify
     rewrite (b * x + a)^3 to (a + b * x)^3
    rewrite -(1 / (3 * b) * (3 * a * b ^ 2 * x ^ 2 / 2 + b ^ 3 * x ^ 3 / 3 + 3 * a ^ 2 * b * x + a ^ 3 * log(x))) to - (a^2 * x + a * b * x^2 / 2 + b^2 * x^3 / 9 + a^3 * log(x) / (3 * b))
    rewrite log(x) * (a + b * x) ^ 3 / (3 * b) to 1/(3 * b) * ((a + b * x) ^ 3 - a ^ 3) * log(x) + a^3 * log(x) / (3 * b)
    simplify
done

// 3
prove (INT x. (a + b * x) ^ 3 * log(x)) = 1/(4 * b) * ((a + b * x) ^ 4 - a ^ 4) * log(x) - (a ^ 3 * x + 3/4 * a ^ 2 * b * x ^ 2 + 1/3 * a * b ^ 2 * x ^ 3 + 1/16 * b^3 * x^4) + SKOLEM_CONST(C) for x > 0, b != 0
lhs:
    integrate by parts with u = log(x), v = (a + b*x)^4 / (4*b)
    simplify
    expand polynomial
    apply integral identity
    simplify
    rewrite -(1 / (4 * b) * (4 * a * b ^ 3 * x ^ 3 / 3 + 3 * a ^ 2 * b ^ 2 * x ^ 2 + b ^ 4 * x ^ 4 / 4 + 4 * a ^ 3 * b * x + a ^ 4 * log(x))) to -(1 / (4 * b)) * (a ^ 4 * log(x) + 4 * a ^ 3 * b * x + 3 * a ^ 2 * b ^ 2 * x ^ 2 + (4/3) * a * b ^ 3 * x ^ 3 + (1/4) * b ^ 4 * x ^ 4)
    rewrite -(1 / (4 * b)) * (a ^ 4 * log(x) + 4 * a ^ 3 * b * x + 3 * a ^ 2 * b ^ 2 * x ^ 2 + (4/3) * a * b ^ 3 * x ^ 3 + (1/4) * b ^ 4 * x ^ 4) to -a^3*x - (3/4)*a^2*b*x^2 - (1/3)*a*b^2*x^3 - (1/16)*b^3*x^4 - (a^4*log(x))/(4*b)
    rewrite -(a^3)*x - 3/4*a^2*b*x^2 - 1/3*a*b^2*x^3 - 1/16*b^3*x^4 - a^4*log(x)/(4*b) + 3*a^2*b*x^2*log(x)/2 + a*b^2*x^3*log(x) + b^3*x^4*log(x)/4 + a^4*log(x)/(4*b) + a^3*x*log(x) + SKOLEM_CONST(C) to (a^3*x + (3/2)*a^2*b*x^2 + a*b^2*x^3 + (1/4)*b^3*x^4)*log(x) - (a^3*x + (3/4)*a^2*b*x^2 + (1/3)*a*b^2*x^3 + (1/16)*b^3*x^4) + SKOLEM_CONST(C)
    rewrite (a^3*x + (3/2)*a^2*b*x^2 + a*b^2*x^3 + (1/4)*b^3*x^4)*log(x) to (1/(4*b))*(4*a^3*b*x + 6*a^2*b^2*x^2 + 4*a*b^3*x^3 + b^4*x^4)*log(x)
    rewrite (4*a^3*b*x + 6*a^2*b^2*x^2 + 4*a*b^3*x^3 + b^4*x^4) to (a + b*x)^4 - a^4
done

### 2.7.2.7

// 1
prove (INT x. log(x)/(a + b * x) ^ m) = 1/(b * (m - 1)) * (-log(x)/(a + b * x) ^ (m - 1) + prove (INT x. 1/x * (a + b * x) ^ (m - 1))) + SKOLEM_CONST(C) for a + b*x > 0, b != 0
sorry

// 2
prove (INT x. log(x)/(a + b * x)) = 1/b * log(x) * log(a + b * x) - 1/b * (INT x. log(a + b * x)/x) + SKOLEM_CONST(C) for a + b*x > 0, x > 0, b != 0
sorry

// 3
prove (INT x. log(x)/(a + b * x) ^ 2) = -log(x)/(b * (a + b * x)) + 1/(a * b) * log(x/(a + b * x)) + SKOLEM_CONST(C) for a + b*x > 0, x > 0, b != 0
sorry

// 4
prove (INT x. log(x)/(a + b * x) ^ 3) = -log(x)/(2 * b * (a + b * x) ^ 2) + 1/(2 * a * b * (a + b * x)) + 1/(2 * a ^ 2 * b) * log(x/(a + b * x)) + SKOLEM_CONST(C) for a + b*x > 0, x > 0, b != 0, a != 0
sorry

// 5a
prove (INT x. log(x)/sqrt(a + b * x)) = 2/b * ((log(x) - 2) * sqrt(a + b * x) - 2 * sqrt(a) * log((sqrt(a + b * x) - sqrt(a))/sqrt(x))) + SKOLEM_CONST(C) for a > 0, x > 0, a + b * x > 0, b > 0
sorry

# 5b
prove (INT x. log(x)/sqrt(a + b * x)) = 2/b * ((log(x) - 2) * sqrt(a + b * x) - 2 * sqrt(-a) * arctan(sqrt((a + b * x)/(-a)))) + SKOLEM_CONST(C) for a < 0, x > 0, a + b * x > 0, b > 0
sorry

### 2.7.2.9

// 1
prove (INT x. x ^ m * log(a + b * x)) = 1 / (m + 1) * (x ^ (m + 1) - (-a) ^ (m + 1)/b^ (m + 1)) * log(a + b * x) + 1 / (m + 1) * SUM(k, 1, m + 1, (-1)^k * x ^ (m - k + 2) * a ^ (k - 1) / (m - k + 2) * b^(k - 1)) + SKOLEM_CONST(C) for a + b*x > 0, x^2 + a^2 > 0, m >= 0, isInt(m)
sorry

// 2
prove (INT x. x * log(a + b * x)) = 1/2 * (x^2 - a^2 / b^2) * log(a + b * x) - 1/2 * (x^2/2 - a*x/b) + SKOLEM_CONST(C) for a + b*x > 0
sorry

// 3
prove (INT x. x^2 * log(a + b * x)) = 1/3 * (x^3 + a^3 / b^3) * log(a + b*x) - 1/3 * (x^3 / 3 - a * x^2 / 2 * b + a^2 * x / b^2) + SKOLEM_CONST(C) for a + b*x > 0, b != 0
sorry

// 4
prove (INT x. x^3 * log(a + b*x)) = 1/4 * (x^4 + a^4 / b^4) * log(a + b*x) - 1/4 * (x^4 / 4 - a * x^4 / 4 * b + a^2 * x^2 / 2 * b^2 - a^3 * x / b^3) + SKOLEM_CONST(C) for a + b*x > 0, b != 0
sorry


#### 2.7.3.3

// 1
prove (INT x. log(x ^ 2 + a ^ 2)) = x * log(x ^ 2 + a ^ 2) - 2 * x + 2 * a * arctan(x/a) + SKOLEM_CONST(C) for a != 0, x^2 + a^2 > 0
lhs:
    rewrite log(x ^ 2 + a ^ 2) to 1 * log(x ^ 2 + a ^ 2)
    integrate by parts with u = log(x^2 + a^2), v = x
    simplify
    rewrite -(2 * (INT x. x ^ 2 / (a ^ 2 + x ^ 2))) to -(2 * (INT x. (x ^ 2 + a ^ 2 - a ^ 2) / (a ^ 2 + x ^ 2)))
    rewrite (x ^ 2 + a ^ 2 - a ^ 2) / (a ^ 2 + x ^ 2) to 1 - a ^ 2 / (a ^ 2 + x ^ 2)
    simplify
    apply integral identity
    rewrite 2 * a ^ 2 * (INT x. 1 / (a ^ 2 + x ^ 2)) to 2 * a * (INT x. a / (a ^ 2 + x ^ 2))
    substitute u for x/a
    simplify
    apply integral identity
    simplify
    rewrite 1 / (a^2 * u^2 + a^2) to (1/a^2) * (1 / (u^2 + 1))
    simplify
    apply integral identity
    replace substitution
done

// 2
prove (INT x. x * log(x ^ 2 + a ^ 2)) = 1/2 * ((x ^ 2 + a ^ 2) * log(x ^ 2 + a ^ 2) - x ^ 2) + SKOLEM_CONST(C) for x^2 + a^2 > 0
lhs:
    integrate by parts with u = log(x^2 + a^2), v = x^2/2
    simplify
    rewrite x^3/(a^2 + x^2) to x - (a^2*x)/(x^2 + a^2)
    simplify
    apply integral identity
    simplify
    rewrite x^2*log(a^2 + x^2)/2 to (1/2)*(x^2 + a^2)*log(x^2 + a^2) - (a^2/2)*log(x^2 + a^2)
    simplify
    substitute u for x^2 + a^2
    simplify
    apply integral identity
    simplify
    rewrite sqrt(-(a^2) + u)/(u * sqrt(-(a^2) + u)) to 1/u
    simplify
    apply integral identity
    simplify
    replace substitution
    simplify
    rewrite log(a^2 + x^2) * (a^2/2 + x^2/2) to (x^2 + a^2)/2 * log(x^2 + a^2)
    rewrite (x^2 + a^2)/2 * log(x^2 + a^2) - x^2/2 to 1/2 * ((x^2 + a^2) * log(x^2 + a^2) - x^2)
done

// 3
prove (INT x. x ^ 2 * log(x ^ 2 + a ^ 2)) = 1/3 * (x ^ 3 * log(x ^ 2 + a ^ 2) - 2/3 * x ^ 3 + 2 * a ^ 2 * x - 2 * a ^ 3 * arctan(x/a)) + SKOLEM_CONST(C) for a != 0, x^2 + a^2 > 0
lhs:
    integrate by parts with u = log(x^2 + a^2), v = x^3 / 3
    simplify
    rewrite x^4 / (3*a^2 + 3*x^2) to (x^2/3 - a^2/3 + a^4/(3*(x^2 + a^2)))
    simplify
    apply integral identity
    simplify
    substitute u for x/a
    simplify
    apply integral identity
    simplify
    rewrite 1/(a^2*u^2 + a^2) to 1/(a^2*(u^2 + 1))
    simplify
    apply integral identity
    simplify
    replace substitution
done

// 4
prove (INT x. x ^ 3 * log(x ^ 2 + a ^ 2)) = 1/4 * ((x ^ 4 - a ^ 4) * log(x ^ 2 + a ^ 2) - x ^ 4/2 + a ^ 2 * x ^ 2) + SKOLEM_CONST(C) for x^2 + a^2 > 0
sorry

// 5
prove (INT x. x ^ 4 * log(x ^ 2 + a ^ 2)) = 1/5 * (x ^ 5 * log(x ^ 2 + a ^ 2) - 2/5 * x ^ 5 + 2/3 * a ^ 2 * x ^ 3 - 2 * a ^ 4 * x + 2 * a ^ 5 * arctan(x/a)) + SKOLEM_CONST(C) for a != 0, x^2 + a^2 > 0
lhs:
    integrate by parts with u = log(x^2 + a^2), v = x^5/5
    simplify
    rewrite x^6 / (5*a^2 + 5*x^2) to (x^4/5 - a^2*x^2/5 + a^4/5 - a^6/(5*(a^2 + x^2)))
    simplify
    apply integral identity
    simplify
    simplify
    substitute u for x/a
    simplify
    apply integral identity
    rewrite 1/(a^2*u^2 + a^2) to (1/a^2)/(u^2 + 1)
    simplify
    apply integral identity
    simplify
    replace substitution
    simplify
done
