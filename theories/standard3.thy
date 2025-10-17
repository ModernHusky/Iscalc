imports standard2

## 2.7 Logarithms and Inverse-Hyperbolic Functions

### 2.7.2 - 2.7.3 Combinations of logarithms and algebraic functions
// Page 238

### 2.7.2.5

// 1
prove (INT x. (a + b * x) ^ m * log(x)) = 1/((m + 1) * b) * ((a + b * x) ^ (m + 1) * log(x) - (INT x. (a + b * x) ^ (m + 1)/x)) for b != 0, m > 0, x > 0, x:real, b:real, m:int, a:real
lhs:
    integrate by parts with u=log(x), v=(a+b*x)^(m+1) / (b*(m+1))
    simplify
    expand polynomial
rhs:
    expand polynomial
done

// 2
prove (INT x. (a + b * x) ^ m * log(x)) = 1/((m + 1) * b) * ((a + b * x) ^ (m + 1) - a ^ (m + 1)) * log(x) - SUM(k, 0, m, binom(m, k) * a ^ (m - k) * b ^ k * x ^ (k + 1)/(k + 1) ^ 2) + SKOLEM_CONST(C) for b != 0, m > 0, x > 0, x:real, b:real, m:int, a:real
lhs:
    integrate by parts with u=log(x), v=(a+b*x)^(m+1) / (b*(m+1))
    apply series expansion on (b*x+a)^(m+1) index k (at 2)
    split sum region at 0
    simplify
    change sum lower to 0
    simplify
    expand polynomial
    rewrite 1 / x * SUM(k, 0, m, a ^ (m - k) * (b * x) ^ (k + 1) * binom(m + 1,k + 1)) to SUM(k, 0, m, 1 / x * a ^ (m - k) * (b * x) ^ (k + 1) * binom(m + 1,k + 1))
    rewrite 1 / x * a^(m-k) * (b*x)^(k+1) to a^(m-k) * (b*x)^(k+1) / x
    apply integral identity
    rewrite (b*x)^(k+1) to b^(k+1)*x^(k+1)
    rewrite a ^ (m - k) * (b ^ (k + 1) * x ^ (k + 1)) / x to a ^ (m - k) * (b ^ (k + 1) * x ^k)
    exchange integral and sum
    apply integral identity
    rewrite binom(m+1, k+1) to (m+1)/(k+1) * binom(m, k)
    rewrite a ^ (m - k) * b ^ (k + 1) * ((m + 1) / (k + 1) * binom(m,k)) * (x ^ (k + 1) / (k + 1)) to b * (m+1) * a ^ (m - k) * b^k * binom(m,k) * x ^ (k + 1) / (k + 1)^2
    simplify
    expand polynomial
    simplify
    rewrite -(b * m / (b * m + b) * SUM(k, 0, m, a ^ (m - k) * b ^ k * x ^ (k + 1) / (k + 1) ^ 2 * binom(m,k))) - b / (b * m + b) * SUM(k, 0, m, a ^ (m - k) * b ^ k * x ^ (k + 1) / (k + 1) ^ 2 * binom(m,k)) to -((b * m + b) / (b * m + b)) * SUM(k, 0, m, a ^ (m - k) * b ^ k * x ^ (k + 1) / (k + 1) ^ 2 * binom(m,k))
    rewrite ((b * m + b) / (b * m + b)) to 1
rhs:
    expand polynomial
done


### 2.7.2.6
// 1
prove (INT x. (a + b * x) * log(x)) = ((a + b * x) ^ 2/(2 * b) - a ^ 2/(2 * b)) * log(x) - (a * x + 1/4 * b * x ^ 2) + SKOLEM_CONST(C) for x > 0, b: real, b != 0
lhs:
    integrate by parts with u = log(x), v = a*x + (b*x^2)/2
    rewrite (1/x)*(b*x^2/2 + a*x) to (b*x)/2 + a
    apply integral identity
    rewrite log(x) * (b * x ^ 2 / 2 + a * x) to (b * x ^ 2 / 2 + a * x) * log(x)
    rewrite b * x ^ 2 / 2 + a * x to (b * x ^ 2 + 2 * a * x)/2
    rewrite (b * x ^ 2 + 2 * a * x)/2 to (x * (b * x + 2 * a))/2
    rewrite (x * (b * x + 2 * a)) / 2 to (a + b * x) ^ 2 / (2 * b) - a ^ 2 / (2 * b)
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
    rewrite log(x) * (a + b * x) ^ 3 / (3 * b) - 1 / (3 * b) * (3 * a * b ^ 2 * x ^ 2 / 2 + b ^ 3 * x ^ 3 / 3 + 3 * a ^ 2 * b * x + a ^ 3 * log(x)) to log(x) * (a + b * x) ^ 3 / (3 * b) - (a^2 * x + a * b * x^2 / 2 + b^2 * x^3 / 9 + a^3 * log(x) / (3 * b))
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
    rewrite 3 * a ^ 2 * b * x ^ 2 * log(x) / 2 - 1 / (4 * b) * (4 * a * b ^ 3 * x ^ 3 / 3 + 3 * a ^ 2 * b ^ 2 * x ^ 2 + b ^ 4 * x ^ 4 / 4 + 4 * a ^ 3 * b * x + a ^ 4 * log(x)) + a * b ^ 2 * x ^ 3 * log(x) + b ^ 3 * x ^ 4 * log(x) / 4 + a ^ 4 * log(x) / (4 * b) + a ^ 3 * x * log(x) to -a^3*x - (3/4)*a^2*b*x^2 - (1/3)*a*b^2*x^3 - (1/16)*b^3*x^4 - (a^4*log(x))/(4*b) + 3*a^2*b*x^2*log(x)/2 + a*b^2*x^3*log(x) + b^3*x^4*log(x)/4 + a^4*log(x)/(4*b) + a^3*x*log(x)
    rewrite -(a^3)*x - 3/4*a^2*b*x^2 - 1/3*a*b^2*x^3 - 1/16*b^3*x^4 - a^4*log(x)/(4*b) + 3*a^2*b*x^2*log(x)/2 + a*b^2*x^3*log(x) + b^3*x^4*log(x)/4 + a^4*log(x)/(4*b) + a^3*x*log(x) + SKOLEM_CONST(C) to (a^3*x + (3/2)*a^2*b*x^2 + a*b^2*x^3 + (1/4)*b^3*x^4)*log(x) - (a^3*x + (3/4)*a^2*b*x^2 + (1/3)*a*b^2*x^3 + (1/16)*b^3*x^4) + SKOLEM_CONST(C)
    rewrite (a^3*x + (3/2)*a^2*b*x^2 + a*b^2*x^3 + (1/4)*b^3*x^4)*log(x) to (1/(4*b))*(4*a^3*b*x + 6*a^2*b^2*x^2 + 4*a*b^3*x^3 + b^4*x^4)*log(x)
    rewrite (4*a^3*b*x + 6*a^2*b^2*x^2 + 4*a*b^3*x^3 + b^4*x^4) to (a + b*x)^4 - a^4
done

### 2.7.2.7

// 1
prove (INT x. log(x)/(a + b * x) ^ m) = 1/(b * (m - 1)) * (-log(x)/(a + b * x) ^ (m - 1) + (INT x. 1/(x * (a + b * x) ^ (m - 1)))) for x > 0, m > 1, b!=0, m:int, x:real, a:real, b:real, a+b*x != 0
lhs:
    integrate by parts with u=log(x), v=(a+b*x)^(1-m)/(b*(1-m))
    expand polynomial
    rewrite -(b * m * x) + b * x to (-b*m + b) * x
    rewrite (b * x + a) ^ (-m + 1) / ((-b*m + b) * x) to 1 / (-b*m + b) * (b * x + a) ^ (-m + 1) / x
    simplify
rhs:
    expand polynomial
    rewrite log(x) * (b * x + a) ^ (-m + 1) / (b * m - b) to -log(x) * (b * x + a) ^ (-m + 1) / (-b * m + b)
    rewrite 1 / (b * m - b) to -1 / (-b * m + b)
    simplify
done

// 2
prove (INT x. log(x)/(a + b * x)) = 1/b * log(x) * log(a + b * x) - 1/b * (INT x. log(a + b * x)/x) for x > 0, b!=0, x:real, a:real, b:real, a+b*x > 0
lhs:
    integrate by parts with u=log(x), v=log(a+b*x)/b
    rewrite log(b * x + a) / (b * x) to 1/b * log(b*x+a)/x
    simplify
rhs:
    simplify
done

// 3
prove (INT x. log(x)/(a + b * x) ^ 2) = -log(x)/(b * (a + b * x)) + 1/(a * b) * log(x/(a + b * x)) + SKOLEM_CONST(C) for a + b*x > 0, x > 0, a != 0, b != 0, a:real, b:real, x:real
lhs:
    integrate by parts with u=log(x), v=(a+b*x)^(-1) / (-b)
    rewrite 1 / (b * x * (b * x + a)) to 1/a * (1/(b*x) - 1/(b*x+a))
    expand polynomial
    apply integral identity
    simplify
rhs:
    rewrite log(x / (a+b*x)) to log(x) - log(a + b*x)
    expand polynomial
done

// 4
prove (INT x. log(x)/(a + b * x) ^ 3) = -log(x)/(2 * b * (a + b * x) ^ 2) + 1/(2 * a * b * (a + b * x)) + 1/(2 * a ^ 2 * b) * log(x/(a + b * x)) + SKOLEM_CONST(C) for  a + b*x > 0, x > 0, a != 0, b != 0, a:real, b:real, x:real
lhs:
    integrate by parts with u=log(x), v=(a+b*x)^(-2)/(-2*b)
    simplify
    rewrite 1 / (x * (b * x + a) ^ 2) to 1/a^2*1/x - b/a^2*1/(a+b*x) - b/a * 1/(a+b*x)^2
    apply integral identity
    simplify
    substitute u for b*x+a
    simplify
    apply integral identity
    replace substitution
    simplify
    rewrite 1 / (2 * b) * (1 / (a * (b * x + a)) - log(b * x + a) / a ^ 2 + log(x) / a ^ 2)  to -log(b * x + a) / (2 * a ^ 2 * b) + 1 / (2 * a * (b * x + a) * b) + log(x) / (2 * a ^ 2 * b)
    simplify
rhs:
    rewrite log(x/(a+b*x)) to log(x) - log(a+b*x)
    rewrite 1 / (2 * a ^ 2 * b) * (log(x) - log(a + b * x)) to -log(b * x + a) / (2 * a ^ 2 * b) + log(x) / (2 * a ^ 2 * b)
    simplify
done

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
prove (INT x. x * log(a + b * x)) = 1/2 * (x^2 - a^2 / b^2) * log(a + b * x) - 1/2 * (x^2/2 - a*x/b) + SKOLEM_CONST(C) for a + b*x > 0, b != 0
lhs:
    integrate by parts with u = log(a + b * x), v = x^2/2
    simplify
    rewrite x^2 / (2*b*x + 2*a) to x / (2*b) - a / (2*b^2) + a^2 / (2*b^2*(b*x + a))
    apply integral identity
    simplify
    expand polynomial
    simplify
    rewrite x ^ 2 * log(b * x + a) / 2 - a ^ 2 * log(b * x + a) / (2 * b ^ 2) to (x^2 - a^2/b^2) / 2 * log(a + b*x)
rhs:
    expand polynomial
    simplify
    rewrite x ^ 2 * log(b * x + a) / 2 - a ^ 2 * log(b * x + a) / (2 * b ^ 2) to (x^2 - a^2/b^2) / 2 * log(a + b*x)
done

// 3

prove (INT x. x^2 * log(a + b * x)) = 1/3 * (x^3 + a^3 / b^3) * log(a + b*x) - 1/3 * (x^3 / 3 - a * x^2 / (2 * b) + a^2 * x / b^2) + SKOLEM_CONST(C) for a + b*x > 0, b != 0
lhs:
    integrate by parts with u = log(a + b * x), v = x^3/3
    simplify
    rewrite x^3 / (3*b*x + 3*a) to x^2 / (3*b) - a*x / (3*b^2) + a^2 / (3*b^3) - a^3 / (3*b^3*(b*x+a))
    apply integral identity
    simplify
    expand polynomial
    simplify
rhs:
    expand polynomial
    simplify
done

// 4
prove (INT x. x^3 * log(a + b*x)) = 1/4 * (x^4 - a^4 / b^4) * log(a + b*x) - 1/4 * (x^4 / 4 - a * x^3 / (3 * b) + a^2 * x^2 / (2 * b^2) - a^3 * x / b^3) + SKOLEM_CONST(C) for a + b*x > 0, b != 0
lhs:
    integrate by parts with u = log(a + b * x), v = x^4/4
    simplify
    rewrite x^4 / (4*b*x + 4*a) to x^3 / (4*b) - a*x^2 / (4*b^2) + a^2*x / (4*b^3) - a^3 / (4*b^4) + a^4 / (4*b^4*(b*x+a))
    apply integral identity
    simplify
    expand polynomial
    simplify
rhs:
    expand polynomial
    simplify
done


#### 2.7.3.3

// 1
prove (INT x. log(x ^ 2 + a ^ 2)) = x * log(x ^ 2 + a ^ 2) - 2 * x + 2 * a * arctan(x/a) + SKOLEM_CONST(C) for a != 0, x^2 + a^2 > 0
lhs:
    rewrite log(x ^ 2 + a ^ 2) to 1 * log(x ^ 2 + a ^ 2)
    integrate by parts with u = log(x^2 + a^2), v = x
    simplify
    rewrite x * log(a ^ 2 + x ^ 2) - 2 * (INT x. x ^ 2 / (a ^ 2 + x ^ 2)) to x * log(a ^ 2 + x ^ 2) - 2 * (INT x. (x ^ 2 + a ^ 2 - a ^ 2) / (a ^ 2 + x ^ 2))
    rewrite (x ^ 2 + a ^ 2 - a ^ 2) / (a ^ 2 + x ^ 2) to 1 - a ^ 2 / (a ^ 2 + x ^ 2)
    simplify
    apply integral identity
    simplify
done

// 2
prove (INT x. x * log(x ^ 2 + a ^ 2)) = 1/2 * ((x ^ 2 + a ^ 2) * log(x ^ 2 + a ^ 2) - x ^ 2) + SKOLEM_CONST(C) for a: real, x^2 + a^2 > 0
lhs:
    integrate by parts with u = log(x^2 + a^2), v = x^2/2
    rewrite x^3/(a^2 + x^2) to x - (a^2*x)/(x^2 + a^2)
    simplify
    apply integral identity
    rewrite x^2*log(a^2 + x^2)/2 to (1/2)*(x^2 + a^2)*log(x^2 + a^2) - (a^2/2)*log(x^2 + a^2)
    simplify
    substitute u for x^2 + a^2
    simplify
    rewrite sqrt(u - a ^ 2) / (u * sqrt(u - a ^ 2)) to 1/u
    apply integral identity
    simplify
    replace substitution
    simplify
    rewrite log(a ^ 2 + x ^ 2) * (a ^ 2 / 2 + x ^ 2 / 2) - x ^ 2 / 2 to 1/2 * ((x^2 + a^2) * log(x^2 + a^2) - x^2)
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
done

// 4
prove (INT x. x ^ 3 * log(x ^ 2 + a ^ 2)) = 1/4 * ((x ^ 4 - a ^ 4) * log(x ^ 2 + a ^ 2) - x ^ 4/2 + a ^ 2 * x ^ 2) + SKOLEM_CONST(C) for x^2 + a^2 > 0
lhs:
    integrate by parts with u = log(x^2 + a^2), v = x^4/4
    rewrite x^5 / (2*a^2 + 2*x^2) to (x^3/2 - a^2*x/2 + a^4*x/(2*(x^2 + a^2)))
    simplify
    apply integral identity
    simplify
    substitute u for x^2 + a^2
    simplify
    rewrite sqrt(u - a ^ 2) / (u * sqrt(u - a ^ 2)) to 1/u
    apply integral identity
    simplify
    replace substitution
    simplify
    rewrite x ^ 4 * log(a ^ 2 + x ^ 2) / 4 - a ^ 4 * log(a ^ 2 + x ^ 2) / 4 to (x ^ 4 - a ^ 4) * log(a ^ 2 + x ^ 2) / 4
rhs:
    expand polynomial
    rewrite -(a ^ 4 * log(a ^ 2 + x ^ 2) / 4) + x ^ 4 * log(a ^ 2 + x ^ 2) / 4 to (x ^ 4 - a ^ 4) * log(a ^ 2 + x ^ 2) / 4
done

// 5
prove (INT x. x ^ 4 * log(x ^ 2 + a ^ 2)) = 1/5 * (x ^ 5 * log(x ^ 2 + a ^ 2) - 2/5 * x ^ 5 + 2/3 * a ^ 2 * x ^ 3 - 2 * a ^ 4 * x + 2 * a ^ 5 * arctan(x/a)) + SKOLEM_CONST(C) for a != 0, x^2 + a^2 > 0
lhs:
    integrate by parts with u = log(x^2 + a^2), v = x^5/5
    simplify
    rewrite x^6 / (5*a^2 + 5*x^2) to (x^4/5 - a^2*x^2/5 + a^4/5 - a^6/(5*(a^2 + x^2)))
    simplify
    apply integral identity
    simplify
done


# Sqrt integral lemmas needed for arcsin/arccos proofs in standard4
# Note: These could be proven but are used as lemmas here

## INT 1/sqrt(a^2-x^2) = arcsin(x/a)
prove (INT x. 1 / sqrt(a ^ 2 - x ^ 2)) = arcsin(x / a) + SKOLEM_CONST(C) for a > 0, x > -a, x < a, a ^ 2 - x ^ 2 > 0
lhs:
    substitute y for x / a
    rewrite a / sqrt(a ^ 2 - a ^ 2 * y ^ 2) to 1 / sqrt(1 - y ^ 2)
    apply integral identity
    replace substitution
done

# Version with relaxed conditions - used by standard4
prove (INT x. 1 / sqrt(a ^ 2 - x ^ 2)) = arcsin(x / a) + SKOLEM_CONST(C) for a != 0, a ^ 2 - x ^ 2 >= 0
sorry

## INT x^2/sqrt(a^2-x^2)
# Proven using integration by parts twice and algebraic manipulation
prove (INT x. x ^ 2 / sqrt(a ^ 2 - x ^ 2)) = a ^ 2 / 2 * arcsin(x / a) - x / 2 * sqrt(a ^ 2 - x ^ 2) + SKOLEM_CONST(C) for a > 0, x > -a, x < a, a ^ 2 - x ^ 2 > 0
lhs:
    # Use x^2 = a^2 - (a^2 - x^2) to split the integral
    rewrite x ^ 2 to a ^ 2 - (a ^ 2 - x ^ 2)
    rewrite (a ^ 2 - (a ^ 2 - x ^ 2)) / sqrt(a ^ 2 - x ^ 2) to a ^ 2 / sqrt(a ^ 2 - x ^ 2) - sqrt(a ^ 2 - x ^ 2)
    simplify
    # Apply the known lemma for INT 1/sqrt(a^2-x^2)
    apply integral identity
    # Now need to compute INT sqrt(a^2-x^2) using integration by parts
    integrate by parts with u = sqrt(a ^ 2 - x ^ 2), v = x
    # Simplify and rewrite the resulting integral
    rewrite -(x ^ 2 / sqrt(a ^ 2 - x ^ 2)) to -(a ^ 2 - (a ^ 2 - x ^ 2)) / sqrt(a ^ 2 - x ^ 2)
    simplify
    # Solve for INT x^2/sqrt(a^2-x^2) which appears on both sides
    solve integral INT x. x ^ 2 / sqrt(a ^ 2 - x ^ 2)
done

# Version with relaxed conditions - used by standard4
prove (INT x. x ^ 2 / sqrt(a ^ 2 - x ^ 2)) = a ^ 2 / 2 * arcsin(x / a) - x / 2 * sqrt(a ^ 2 - x ^ 2) + SKOLEM_CONST(C) for a != 0, a ^ 2 - x ^ 2 >= 0
sorry