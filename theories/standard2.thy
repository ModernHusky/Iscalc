## 9.2 Integrands involving powers of x and powers of sin x or cos x

### 9.2.1 Integrands involving x^n sin^m x

#### 9.2.1.1

// 3a
prove (INT x. sin(x) ^ 3) = 1/12 * cos(3 * x) - 3/4 * cos(x) + SKOLEM_CONST(C)
sorry

// 3b
prove (INT x. sin(x) ^ 3)  = 1/3 * cos(x) ^ 3 - cos(x) + SKOLEM_CONST(C)
lhs:
    rewrite sin(x)^3 to sin(x)*sin(x)^2
    rewrite sin(x)^2 to 1 - cos(x)^2
    simplify
    integrate by parts with u=1 - cos(x)^2, v=-cos(x)
    simplify
    substitute u for cos(x)
    apply integral identity
    replace substitution
    simplify
    rewrite -(cos(x) * (-(cos(x) ^ 2) + 1)) to cos(x)^3 - cos(x)
    rewrite cos(x)^3 - cos(x) - 2*cos(x)^3/3 to (cos(x)^3 - 2*cos(x)^3/3) - cos(x)
    rewrite cos(x)^3 - 2*cos(x)^3/3 to (1/3)*cos(x)^3
    rewrite (1/3)*cos(x)^3 - cos(x) to (1/3)*cos(x)^3 - cos(x)
done

// 4a
prove (INT x. sin(x) ^ 4) = 1/32 * sin(4 * x) - 1/4 * sin(2 * x) + 3/8 * x + SKOLEM_CONST(C)
sorry

// 4b
prove (INT x. sin(x) ^ 4) = -1/4 * sin(x) ^ 3  * cos(x) - 3/8 * sin(x) * cos(x) + 3/8 * x + SKOLEM_CONST(C)
sorry

// 5a
prove (INT x. sin(x) ^ 5) = -1/80 * cos(5 * x) + 5/48 * cos(3 * x) - 5/8 * cos(x) + SKOLEM_CONST(C)
sorry

// TODO: "apply integral identity" has a bug
// 5b
prove (INT x. sin(x) ^ 5) = -1/5 * sin(x) ^ 4 * cos(x) + 4/15 * cos(x) ^ 3 - 4/5 * cos(x) + SKOLEM_CONST(C)
lhs:
    rewrite sin(x)^5 to sin(x)^4 * sin(x)
    integrate by parts with u = sin(x)^4, v = -cos(x)
    simplify
    rewrite cos(x)^2 to 1 - sin(x)^2
    expand polynomial
    simplify
    integrate by parts with u = sin(x)^2, v = -cos(x)
    simplify
    rewrite sin(x)^2 to 1 - cos(x)^2
    expand polynomial
    simplify
    apply integral identity
    simplify
    substitute u for cos(x)
    simplify
    apply integral identity
    replace substitution
    solve integral INT x. sin(x)^5
done

// 6
prove (INT x. sin(x) ^ (2 * n)) = 1/(2 ^ (2 * n)) * binom(2 * n, n) * x + (-1) ^ n/(2 ^ (2*n - 1)) * SUM(k, 0, n-1, (-1) ^ k * binom(2 * n, k) * sin(2 * n-2 * k) * x / (2 * n-2 * k)) + SKOLEM_CONST(C)
sorry

// 7
prove (INT x. sin(x) ^ (2 * n + 1)) = 1/(2 ^ (2 * n)) * (-1) ^ (n + 1) * SUM(k, 0, n, (-1)^k * binom(2 * n + 1, k) * cos(2 * n + 1 -2 * k) * x / (2 * n + 1 -2 * k)) + SKOLEM_CONST(C)
sorry

// 8
prove (INT x. x * sin(x)) = sin(x) - x * cos(x) + SKOLEM_CONST(C)
lhs:
    integrate by parts with u = x, v = -cos(x)
    simplify
    apply integral identity
done

// 9
prove (INT x. x ^ 2 * sin(x)) = 2 * x * sin(x) - (x ^ 2 - 2) * cos(x) + SKOLEM_CONST(C)
sorry

// 10
prove (INT x. x ^ 3 * sin(x)) = (3 * x ^ 2 - 6) * sin(x) - (x ^ 3 - 6 * x) * cos(x) + SKOLEM_CONST(C)
lhs:
    integrate by parts with u = x^3, v = -cos(x)
    simplify
    integrate by parts with u = x^2, v = sin(x)
    simplify
    integrate by parts with u = x, v = -cos(x)
    simplify
    apply integral identity
    simplify
    rewrite 3 * x ^ 2 * sin(x) - x ^ 3 * cos(x) + 6 * x * cos(x) - 6 * sin(x) to (3 * x ^ 2 - 6) * sin(x) - (x ^ 3 - 6 * x) * cos(x)
done

// 11
prove (INT x. x ^ 4 * sin(x)) = (4 * x ^ 3 - 24 * x) * sin(x) - (x ^ 4 - 12 * x ^ 2 + 24) * cos(x) + SKOLEM_CONST(C)
sorry

// 12
prove (INT x. x ^ (2 * n) * sin(x)) = factorial(2 * n) * SUM(k, 0, n, (-1) ^ (k + 1) * (x ^ (2 * n - 2 * k))/factorial(2 * n -2 * k) *  cos(x)) +  SUM(k, 0, n-1, (-1)^k * (x ^ (2 * n - 2 * k - 1))/factorial(2 * n - 2 * k -1) *  sin(x)) + SKOLEM_CONST(C) for isInt(n)
sorry

// 13
prove (INT x. x ^ (2 * (n + 1)) * sin(x)) = factorial(2 * n + 1) * SUM(k, 0, n, (-1) ^ (k + 1) * (x ^ (2 * n - 2 * k + 1))/factorial(2 * n -2 * k + 1) *  cos(x)) + SUM(k, 0, n-1, (-1)^k * (x ^ (2 * n - 2 * k))/factorial(2 * n - 2 * k) *  sin(x)) + SKOLEM_CONST(C) for isInt(n)
sorry

// 14
prove (INT x. sin(x) ^ 2) = 1/2 * x - 1/4 * sin(2 * x) + SKOLEM_CONST(C)
lhs:
    rewrite sin(x)^2 to (1 - cos(2*x))/2
    simplify
    apply integral identity
    simplify
    substitute u for 2*x
    apply integral identity
    replace substitution
done

// 15
prove (INT x. x * sin(x) ^ 2) = 1/4 * x ^ 2 -  1/4 * x * sin(2 * x) -1/8 * cos(2 * x) + SKOLEM_CONST(C)
sorry

// 16
prove (INT x. x ^ 2 * sin(x) ^ 2) = 1/6 * x ^ 3 -1/4 * x * cos(2 * x) -1/4 * (x ^ 2 - 1/2) * sin(2 * x) + SKOLEM_CONST(C)
sorry

// 17
prove (INT x. x ^ m * sin(x) ^ n) = x ^ (m - 1) * sin(x) ^ (n - 1) * x / n ^ 2 * (m * sin(x) - n * x * cos(x)) + (n - 1) / n * (INT x. x ^ m * sin(x) ^ (n - 2)) - m * (m -1)/n ^ 2 * (INT x. x ^ (m - 2) * sin(x) ^ n)
sorry

### 9.2.2 Integrands involving x^-n sin^m x

#### 9.2.2.1

// 1
prove (INT x. sin(x) / x) = SUM(k, 0, oo, (-1) ^ k * x ^ (2 * k+1) /((2 * k+1) * factorial(2 * k+1))) + SKOLEM_CONST(C)
sorry

// 2
prove (INT x. sin(x) / x ^ 2) = -sin(x)/x + (INT x. cos(x)/x)
sorry

// 3
prove (INT x. sin(x) / x ^ 3) = -sin(x)/(2 * x ^ 2) - cos(x)/(2 * x) - 1/2 * (INT x. sin(x)/x)
sorry

// 4
prove (INT x. sin(x) / x ^ n) = -sin(x)/((n - 1) * x ^ (n - 1)) - cos(x)/((n - 1) * (n - 2) * x ^ (2 * n) ) - 1/((n - 1) * (n - 2)) * (INT x. sin(x)/x ^ (n - 2)) for n > 2
sorry

// 5
prove (INT x. sin(x) ^ m / x ^ n) = -(sin(x) ^ (m-1) * x * ((n - 2)  * sin(x) + m * x * cos(x)))/((n - 1) * (n - 2) * x ^ (n - 1)) - m ^ 2/((n - 1) * (n - 2)) * (INT x. sin(x) ^ m/x ^ (n - 2)) + m * (m -1)/((n - 1) * (n - 2)) * (INT x. sin(x) ^ (m - 2)/x ^ (n - 2)) for n != 1, n != 2
sorry

### 9.2.3 Integrands involving x^n sin^-m x

#### 9.2.3.1

// 1a
prove (INT x. 1 / sin(x)) = log(abs(tan(x/2))) + SKOLEM_CONST(C) for abs(x) < pi
sorry

// 1b
prove (INT x. 1 / sin(x)) = -1/2 * log((1 + cos(x))/(1 - cos(x))) + SKOLEM_CONST(C) for abs(x) < pi
sorry

// 2
prove (INT x. 1 / sin(x) ^ 2) = -cot(x) + SKOLEM_CONST(C) for abs(x)<pi
sorry

// 3
prove (INT x. 1 / sin(x) ^ 3) = -cos(x)/(2 * sin(x) ^ 2) + 1/2 * log(abs(tan(x/2))) + SKOLEM_CONST(C) for abs(x)<pi
sorry

// 4
prove (INT x. 1 / sin(x) ^ n) = -cos(x)/((n - 1) * sin(x) ^ (n - 1)) + (n - 2) / (n - 1) * (INT x. 1/sin(x) ^ (n - 2)) + SKOLEM_CONST(C) for abs(x) < pi, n > 1
sorry

// The rest involves Bernoulli numbers

### 9.2.4 Integrands involving x^n cos^m x

#### 9.2.4.1

// 3a
prove (INT x. cos(x) ^ 3)  = 1/12 * sin(3 * x) + 3/4 * sin(x) + SKOLEM_CONST(C)
sorry

// 3b
prove (INT x. cos(x) ^ 3)  = sin(x) - 1/3 * sin(x) ^ 3 + SKOLEM_CONST(C)
lhs:
    rewrite cos(x)^3 to cos(x)^2 * cos(x)
    rewrite cos(x)^2 to 1 - sin(x)^2
    simplify
    substitute u for sin(x)
    apply integral identity
    replace substitution
done

// 4a
prove (INT x. cos(x) ^ 4)  = 1/32 * sin(4 * x) + 1/4 * sin(2 * x) + 3/8 * x + SKOLEM_CONST(C)
lhs:
    rewrite cos(x)^4 to (cos(x)^2)^2
    rewrite cos(x)^2 to (cos(2*x) + 1)/2
    expand polynomial
    rewrite cos(2*x)^2 to (cos(4*x) + 1)/2
    simplify
    apply integral identity
    simplify
    substitute u for 2*x
    apply integral identity
    replace substitution
    substitute v for 4*x
    apply integral identity
    replace substitution
    simplify
done

// 4b
prove (INT x. cos(x) ^ 4)  = 1/4 * sin(x) * cos(x) ^ 3 + 3/8 * sin(x) * cos(x) + 3/8 * x + SKOLEM_CONST(C)
sorry

// 5a
prove (INT x. cos(x) ^ 5)  = 1/80 * sin(5 * x) + 5/48 * sin(3 * x) + 5/8 * sin(x) + SKOLEM_CONST(C)
sorry

// 5b
prove (INT x. cos(x) ^ 5)  = 1/5 * cos(x) ^ 4 * sin(x) - 4/15 * sin(x) ^ 3 + 4/5 * sin(x) + SKOLEM_CONST(C)
sorry

// 6
prove (INT x. cos(x) ^ (2 * n))  = 1/2^(2 * n) * binom(2 * n, n) * x + 1/2^(2 * n - 1) * SUM(k, 0, n - 1, binom(2 * n, k) * sin(2 * n -2 * k) * x/(2 * n -2 * k)) + SKOLEM_CONST(C)
sorry

// 7
prove (INT x. cos(x) ^ (2 * n + 1))  = 1/2 ^(2 * n) * SUM(k, 0, n, binom(2 * n + 1, k) * sin(2 * n -2 * k + 1) * x/(2 * n -2 * k + 1)) + SKOLEM_CONST(C)
sorry

// 8
prove (INT x. x * cos(x))  = cos(x) + x * sin(x) + SKOLEM_CONST(C)
lhs:
    integrate by parts with u = x, v = sin(x)
    simplify
    apply integral identity
    simplify
done

// 9
prove (INT x. x ^ 2 * cos(x))  = 2 * x * cos(x) + (x ^ 2 - 2) * sin(x) + SKOLEM_CONST(C)
lhs:
    integrate by parts with u = x^2, v = sin(x)
    simplify
    integrate by parts with u = x, v = -cos(x)
    simplify
    apply integral identity
    simplify
    rewrite x^2 * sin(x) + 2 * x * cos(x) - 2 * sin(x) to (x^2 - 2) * sin(x) + 2 * x * cos(x)
done

// 10
prove (INT x. x ^ 3 * cos(x))  = (3 * x ^ 2 - 6) * cos(x) + (x ^ 3 - 6 * x) * sin(x) + SKOLEM_CONST(C)
lhs:
    integrate by parts with u = x^3, v = sin(x)
    simplify
    integrate by parts with u = x^2, v = -cos(x)
    simplify
    integrate by parts with u = x, v = sin(x)
    simplify
    apply integral identity
    simplify
    rewrite 3 * x^2 * cos(x) + x^3 * sin(x) - 6 * x * sin(x) - 6 * cos(x) + SKOLEM_CONST(C) to (3 * x^2 * cos(x) - 6 * cos(x)) + (x^3 * sin(x) - 6 * x * sin(x)) + SKOLEM_CONST(C)
    rewrite (3 * x^2 * cos(x) - 6 * cos(x)) to (3 * x^2 - 6) * cos(x)
    rewrite (x^3 * sin(x) - 6 * x * sin(x)) to (x^3 - 6 * x) * sin(x)
done

// 11
prove (INT x. x ^ 4 *cos(x))  = (4 * x ^ 3 - 24 * x) * cos(x) + (x ^ 4 - 12 * x ^ 2 + 24) * sin(x) + SKOLEM_CONST(C)
lhs:
    integrate by parts with u = x^4, v = sin(x)
    simplify
    integrate by parts with u = x^3, v = -cos(x)
    simplify
    integrate by parts with u = x^2, v = sin(x)
    simplify
    integrate by parts with u = x, v = -cos(x)
    simplify
    apply integral identity
    simplify
    rewrite -(12 * x ^ 2 * sin(x)) + 4 * x ^ 3 * cos(x) + x ^ 4 * sin(x) - 24 * x * cos(x) + 24 * sin(x) + SKOLEM_CONST(C) to (4 * x ^ 3 - 24 * x) * cos(x) + (x ^ 4 - 12 * x ^ 2 + 24) * sin(x) + SKOLEM_CONST(C)
done

// 12
prove (INT x. x ^ (2 * n) * cos(x)) = factorial(2 * n) * SUM(k, 0, n, (-1) ^ k * (x ^ (2 * n - 2 * k))/factorial(2 * n -2 * k) *  sin(x)) +  SUM(k, 0, n-1, (-1)^k * (x ^ (2 * n - 2 * k - 1))/factorial(2 * n - 2 * k -1) *  cos(x)) + SKOLEM_CONST(C) for isInt(n)
sorry

// 13
prove (INT x. x ^ (2 * (n + 1)) * cos(x)) = factorial(2 * n + 1) * SUM(k, 0, n, (-1) ^ k * (x ^ (2 * n - 2 * k + 1))/factorial(2 * n -2 * k+ 1) *  sin(x)) + SUM(k, 0, n, (-1)^k * (x ^ (2 * n - 2 * k))/factorial(2 * n - 2 * k) *  cos(x)) + SKOLEM_CONST(C) for isInt(n)
sorry
