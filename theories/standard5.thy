imports base

# Handbook of mathematical formulas and integrals
## Chapter 4 Indefinite Integrals of Algebraic Functions

// Page 154

### 4.2.2 Integrands Involving a + bx
#### 4.2.2.1

// 1
prove (INT x. (a + b * x) ^ n) = (a + b * x) ^ (n + 1) / (b * (n + 1)) + SKOLEM_CONST(C) for n: int, n != -1, n >= 0, b != 0
lhs:
    substitute u for a + b * x
    simplify
    apply integral identity
    replace substitution
    simplify
done

// 2
prove (INT x. 1 / (a + b * x)) = 1 / b * log(a + b * x) + SKOLEM_CONST(C) for a + b * x > 0, b != 0
lhs:
    substitute u for a + b * x
    simplify
    apply integral identity
    replace substitution
    simplify
done

// 3
prove (INT x. 1 / (a + b * x) ^ 2) = -1 / (b * (a + b * x)) + SKOLEM_CONST(C) for a + b * x != 0, b != 0
lhs:
    substitute u for a + b*x
    simplify
    apply integral identity
    replace substitution
done

// 4
prove (INT x. 1 / (a + b * x) ^ 3) = -1 / (2 * b * (a + b * x) ^ 2) + SKOLEM_CONST(C) for a + b * x != 0, b != 0
lhs:
    substitute u for a + b * x
    simplify
    apply integral identity
    replace substitution
done

// 5
prove (INT x. 1 / (a + b * x) ^ n) = -1 / (b * (n - 1) * (a + b * x) ^ (n - 1)) + SKOLEM_CONST(C) for n: int, a + b * x != 0, b != 0, n != 1, n >= 2
lhs:
    substitute u for a + b * x
    simplify
    apply integral identity
    replace substitution
    simplify
    rewrite (b * x + a) ^ (1 - n) / (b * (1 - n)) to -1 / (b * (n - 1) * (a + b * x) ^ (n - 1))
done

#### 4.2.2.2
// 1
prove (INT x. x / (a + b * x)) = x / b - a / b ^ 2 * log(a + b * x) + SKOLEM_CONST(C) for a + b * x > 0, b != 0
lhs:
    rewrite x / (a + b * x) to 1 / b - a / (b * (a + b * x))
    simplify
    apply integral identity
    substitute u for a + b * x
    apply integral identity
    replace substitution
    simplify
done

// 2
prove (INT x. x / (a + b * x) ^ 2) = -x / (b * (a + b * x)) + 1 / (b ^ 2) * log(a + b * x) + SKOLEM_CONST(C) for a + b * x > 0, b != 0
lhs:
    integrate by parts with u = x, v = -1 / (b * (a + b * x))
    simplify
    substitute w for a + b * x
    simplify
    apply integral identity
    replace substitution
    simplify
done

// 3
prove (INT x. x / (a + b * x) ^ 3) = -(x / b + a / (2 * b ^ 2)) / (a + b * x) ^ 2 + SKOLEM_CONST(C) for a + b * x != 0, b != 0
lhs:
    substitute u for a + b * x
    simplify
    expand polynomial
    apply integral identity
    simplify
    replace substitution
    simplify
    rewrite 1 / b ^ 2 * (a / (2 * (b * x + a) ^ 2) - 1 / (b * x + a)) to (a / (2 * b ^ 2) - (b * x + a) / b ^ 2) / (b * x + a) ^ 2
    simplify
    rewrite 1 / (b * x + a) ^ 2 * (a / (2 * b ^ 2) - (b * x + a) / b ^ 2) to -(x / b + a / (2 * b ^ 2)) / (b * x + a) ^ 2
done

// 4
// This reduction formula requires complex algebraic manipulation that is beyond// the current automation capabilities of the system.
prove (INT x. x / (a + b * x) ^ n) = x / (b * (2 - n) * (a + b * x) ^ (n - 1)) - a / (b * (2 - n)) * INT x. 1 / (a + b * x) ^ n for a + b * x != 0, n != 2, b != 0
sorry

#### 4.2.2.3
// 1
prove (INT x. x ^ 2 / (a + b * x)) = x ^ 2 / (2 * b) - a * x / b ^ 2 + a ^ 2 / b ^ 3 * log(abs(a + b * x)) + SKOLEM_CONST(C) for a + b * x != 0, b != 0
lhs:
    rewrite x ^ 2 / (a + b * x) to (x / b - a / b ^ 2) + (a ^ 2 / b ^ 2) / (a + b * x)
    simplify
    apply integral identity
    simplify
    substitute u for a + b * x
    apply integral identity
    replace substitution
    simplify
done

// 2
// Requires complex algebraic rewriting to match the expected form
prove (INT x. x ^ 2 / (a + b * x) ^ 2) = x / b ^ 2 - a ^ 2 / (b ^ 3 * (a + b * x)) - 2 * a / b ^ 3 * log(abs(a + b * x)) + SKOLEM_CONST(C) for a + b * x != 0, b != 0
sorry

// 3
prove (INT x. x ^ 2 / (a + b * x) ^ 3) = (2 * a * x / b ^ 2 + 3 * a ^ 2 / (2 * b ^ 3)) / (a + b * x) + log(abs(a + b * x)) / b ^ 3 + SKOLEM_CONST(C) for a + b * x != 0, b != 0
sorry

// 4
prove (INT x. x ^ 2 / (a + b * x) ^ n) = x ^ 2 / (b * (3 - n) * (a + b * x) ^ (n - 1)) - 2 * a / (b * (3 - n)) * INT x. x / (a + b * x) ^ n for n = 3, a + b * x != 0, b != 0
sorry

#### 4.2.2.4
// 1
prove (INT x. x ^ 3 / (a + b * x)) = x ^ 3 / (3 * b) - a * x ^ 2 / (2 * b ^ 2) + a ^ 2 * x / b ^ 3 - a ^ 3 / b ^ 4 * log(abs(a + b * x)) + SKOLEM_CONST(C) for a + b * x != 0, b != 0
lhs:
    rewrite x^3 / (a + b*x) to (x^2 / b - a*x / b^2 + a^2 / b^3 - a^3 / (b^3 * (a + b*x)))
    simplify
    apply integral identity
    simplify
    substitute u for a + b*x
    simplify
    apply integral identity
    simplify
    replace substitution
done

// 2
prove (INT x. x ^ 3 / (a + b * x) ^ 2) = x ^ 2 / (2 * b ^ 2) - 2 * a * x / b ^ 3 + a ^ 3 / (b ^ 4 * (a + b * x)) + 3 * a ^ 2 / b ^ 4 * log(abs(a + b * x)) + SKOLEM_CONST(C) for a + b * x > 0, b != 0
sorry

// 3
prove (INT x. x ^ 3 / (a + b * x) ^ 3) = (x ^ 3 / b + 2 * a * x ^ 2 / b ^ 2 - 2 * a ^ 2 * x / b ^ 3 - 5 * a ^ 3 / (2 * b ^ 4)) / (a + b * x) ^ 2 - 3 * a / b ^ 4 * log(abs(a + b * x)) + SKOLEM_CONST(C) for a + b * x != 0, b != 0
sorry

// 4
prove (INT x. x ^ 3 / (a + b * x) ^ n) = x ^ 3 / (b * (4 - n) * (a + b * x) ^ (n - 1)) - 3 * a / (b * (4 - n)) * INT x. x ^ 2 / (a + b * x) ^ n for n != 4, b != 0, a + b * x != 0
sorry

#### 4.2.2.5

// 1
prove (INT x. x ^ m / (a + b * x) ^ n) = -(x ^ m) / (b * (m + 1 - n) * (a + b * x) ^ (n - 1)) - m * a / (b * (m + 1 - n)) * INT x. x ^ (m - 1) / (a + b * x) ^ n for isInt(m), isInt(n), m > 0, n > 0, m != n - 1, b != 0, a + b * x != 0
sorry

// 2
prove (INT x. x ^ (n - 1) / (a + b * x) ^ n) = x ^ (n - 1) / (b * (n - 1) * (a + b * x) ^ (n - 1)) + 1 / b * INT x. x ^ (n - 2) / (a + b * x) ^ (n - 1) for isInt(n), m > 0, n > 0, n != 1, b != 0, a + b * x != 0
sorry

#### 4.2.2.6
// 1
prove (INT x. x ^ n / (a + b * x)) = SUM(k, 1, n, x ^ k * a ^ (n - k) * (-1) ^ (n - k) / (k * b ^ (n - k + 1))) + (-1) ^ n * a ^ n / b ^ (n + 1) * log(abs(a + b * x)) + SKOLEM_CONST(C) for a + b * x != 0, b != 0
sorry

// 2
prove (INT x. x ^ n / (a + b * x) ^ 2) = SUM(k, 1, n-1, (-1) ^ (k - 1) * k * a ^ (k - 1) * x ^ (n - k) / ((n - k) * b ^ (k + 1))) + (-1) ^ (n - 1) * a ^ n / (b ^ (n + 1) * (a + b * x)) + (-1) ^ (n + 1) * n * a ^ (n - 1) / b ^ (n + 1) * log(abs(a + b * x)) + SKOLEM_CONST(C) for b != 0, a + b * x != 0
sorry

#### 4.2.2.7
// 1
prove (INT x. 1 / (x * (a + b * x))) = -1 / a * log(abs((a + b * x) / x)) + SKOLEM_CONST(C) for x != 0 ,a + b * x != 0, a != 0
lhs:
    rewrite 1 / (x * (a + b * x)) to 1 / (a * x) - b / (a * (a + b * x))
    apply integral identity
    simplify
    substitute u for a + b * x
    apply integral identity
    replace substitution
    simplify
    rewrite log(abs(x)) / a - b * log(abs(b * x + a)) / (a * b) to (-log(abs(b * x + a)) + log(abs(x))) / a
    simplify
    rewrite log(abs(x) / abs(b * x + a)) to -log(abs(b * x + a) / abs(x))
    simplify
done

// 2
prove (INT x. 1 / (x * (a + b * x) ^ 2)) = 1 / (a * (a + b * x)) - 1 / a ^ 2 * log(abs((a + b * x) / x)) + SKOLEM_CONST(C) for a + b * x != 0, x != 0, a != 0
lhs:
    partial fraction decomposition
    apply integral identity
    simplify
    substitute u for a + b * x
    apply integral identity
    replace substitution
    simplify
    substitute v for a + b * x
    apply integral identity
    replace substitution
    simplify
    rewrite b / (a * b * (b * x + a)) - b * log(abs(b * x + a)) / (a ^ 2 * b) + log(abs(x)) / a ^ 2 to b / (a * b * (b * x + a)) - log(abs(a + b * x)) / a ^ 2 + log(abs(x)) / a ^ 2
    rewrite b / (a * b * (b * x + a)) to 1 / (a * (a + b * x))
    simplify
    rewrite 1 / (a * (b * x + a)) - log(abs(b * x + a)) / a ^ 2 + log(abs(x)) / a ^ 2 to 1 / (a * (a + b * x)) + (log(abs(x)) - log(abs(a + b * x))) / a ^ 2
    simplify
    rewrite log(abs(x) / abs(b * x + a)) to -log(abs((a + b * x) / x))
    rewrite 1 / (a * (b * x + a)) to 1 / (a * (a + b * x))
    simplify
done

// 3
prove (INT x. 1 / (x * (a + b * x) ^ 3)) = (3 / (2 * a) + b * x / a ^ 2) - 1 / a ^ 3 * log(abs((a + b * x) / x)) + SKOLEM_CONST(C) for a + b * x != 0, x != 0, a != 0
sorry

// 4
prove (INT x. 1 / (x * (a + b * x) ^ n)) = 1 / (a * (n - 1) * (a + b * x) ^ (n - 1)) + 1 / a * INT x. 1 / (x * (a + b * x) ^ (n - 1)) for n != 1, a + b * x != 0, x != 0, a != 0
sorry

#### 4.2.2.8
// 1
prove (INT x. 1 / (x ^ 2 * (a + b * x))) = -1 / (a * x) + b / a ^ 2 * log(abs((a + b * x) / x)) + SKOLEM_CONST(C) for x != 0, a + b * x != 0, a != 0
lhs:
    rewrite 1 / (x ^ 2 * (a + b * x)) to 1 / (a * x ^ 2) - b / (a * x * (a + b * x))
    apply integral identity
    simplify
done

// 2
prove (INT x. 1 / (x ^ 2 * (a + b * x) ^ 2)) = -(1 / (a * x) + 2 * b / a ^ 2) / (a + b * x) + 2 * b / a ^ 3 * log(abs((a + b * x) / x)) + SKOLEM_CONST(C) for x != 0, a + b * x != 0, a != 0
sorry

// 3
prove (INT x. 1 / (x ^ 2 * (a + b * x) ^ 3)) = -(1 / (a * x) + 9 * b / (2 * a ^ 2) + 3 * b ^ 2 * x / a ^ 3) / (a + b * x) ^ 2 + 3 * b / a ^ 4 * log(abs((a + b * x) / x)) + SKOLEM_CONST(C) for x != 0, a + b * x != 0, a != 0
sorry

// 4
prove (INT x. 1 / (x ^ 2 * (a + b * x) ^ n)) = -1 / (a * x * (a + b * x) ^ (n - 1)) - n * b / a * INT x. 1 / (x * (a + b * x) ^ n) for x != 0, a + b * x != 0, a != 0
sorry

### 4.2.3 Integrands Involving Linear Factors
// Page 157
#### 4.2.3.2
// 1
prove (INT x. (a + b * x) / (c + d * x)) = b * x / d + ((a * d - b * c) / d ^ 2) * log(abs(c + d * x)) + SKOLEM_CONST(C) for d != 0, c + d * x != 0
lhs:
    rewrite (a + b*x)/(c + d*x) to b/d + (a*d - b*c)/d^2 * d/(c + d*x)
    simplify
    apply integral identity
    substitute u for c + d*x
    simplify
    apply integral identity
    replace substitution
done

// 2
prove (INT x. 1 / ((x - a) * (x - b))) = 1 / (a - b) * log(abs((x - a) / (x - b))) + SKOLEM_CONST(C) for x != a, x != b, a != b
lhs:
    rewrite 1 / ((x - a) * (x - b)) to 1 / ((a - b) * (x - a)) - 1 / ((a - b) * (x - b))
    apply integral identity
    substitute u for x - a
    apply integral identity
    replace substitution
    substitute v for x - b
    apply integral identity
    replace substitution
    simplify
    rewrite log(abs(x - a)) / (a - b) - log(abs(x - b)) / (a - b) to (log(abs(x - a)) - log(abs(x - b))) / (a - b)
    simplify
done

// 3
prove (INT x. x / ((x - a) * (x - b))) = a / (a - b) * log(abs(x - a)) - b / (a - b) * log(abs(x - b)) + SKOLEM_CONST(C) for x != a, x != b, a != b
lhs:
    rewrite x / ((x - a) * (x - b)) to a / ((a - b) * (x - a)) - b / ((a - b) * (x - b))
    apply integral identity
    substitute u for x - a
    apply integral identity
    replace substitution
    substitute v for x - b
    apply integral identity
    replace substitution
    simplify
done

// 4
prove (INT x. 1 / ((x - a) ^ 2 * (x - b))) = -1 / ((x - a) * (a - b)) - 1 / (a - b) ^ 2 * log(abs((x - a) / (x - b))) + SKOLEM_CONST(C) for x != a, a != b, x != b
lhs:
    rewrite 1 / ((x - a) ^ 2 * (x - b)) to 1 / ((a - b) ^ 2 * (x - b)) - 1 / ((a - b) ^ 2 * (x - a)) + 1 / ((a - b) * (x - a) ^ 2)
    apply integral identity
    substitute u for x - b
    apply integral identity
    replace substitution
    substitute v for x - a
    apply integral identity
    replace substitution
    substitute w for x - a
    apply integral identity
    replace substitution
    simplify
    rewrite log(abs(x - b)) / (a - b) ^ 2 - log(abs(x - a)) / (a - b) ^ 2 to (log(abs(x - b)) - log(abs(x - a))) / (a - b) ^ 2
    rewrite (log(abs(x - b)) - log(abs(x - a))) / (a - b) ^ 2 to log(abs((x - b) / (x - a))) / (a - b) ^ 2
    rewrite log(abs((x - b) / (x - a))) to -log(abs((x - a) / (x - b)))
    simplify
done

// 5
prove (INT x. 1 / ((x - a) ^ 2 * (x - b) ^ 2)) = (a + b - 2 * x) / ((x - a) * (x - b) * (a - b) ^ 2) - 2 / (a - b) ^ 3 * log(abs((x - a) / (x - b))) + SKOLEM_CONST(C) for x != a, a != b, x != b
sorry

// 6
prove (INT x. x / ((x - a) ^ 2 * (x - b))) = -a / ((x - a) * (a - b)) + b / (a - b) ^ 2 * log(abs((x - b) / (x - a))) + SKOLEM_CONST(C) for x != a, x != b, a != b
sorry

// 7
prove (INT x. x / ((x - a) ^ 2 * (x - b) ^ 2)) = (2 * a * b - (a + b) * x) / ((x - a) * (x - b) * (a - b) ^ 2) + (a + b) / (a - b) ^ 3 * log(abs((x - b) / (x - a))) + SKOLEM_CONST(C) for x != a, x != b, a != b
sorry

// 8
prove (INT x. x ^ 2 / ((x - a) ^ 2 * (x - b) ^ 2)) = (a * b * (a + b) - (a ^ 2 + b ^ 2) * x) / ((x - a) * (x - b) * (a - b) ^ 2) + 2 * a * b / (a - b) ^ 3 * log(abs((x - b) / (x - a))) + SKOLEM_CONST(C) for x != a, x != b, a != b
sorry

### 4.2.4 Integrands Involving a ^ 2 +- b ^ 2 * x ^ 2
#### 4.2.4.1