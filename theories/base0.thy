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

# Common series expansion

axiom exp(x) = SUM(n, 0, oo, x ^ n / factorial(n))

axiom sin(x) = SUM(n, 0, oo, (-1)^n * x^(2*n+1) / factorial(2*n+1))

axiom cos(x) = SUM(n, 0, oo, (-1)^n * x^(2*n) / factorial(2*n))

axiom arctan(x) = SUM(n, 0, oo, (-1)^n * x^(2*n+1) / (2*n+1)) for x >= -1, x <= 1

axiom (1 + x) ^ -1 = SUM(n, 0, oo, (-1)^n * x ^ n) for x > -1, x < 1

axiom (1 - x) ^ -1 = SUM(n, 0, oo, x ^ n) for x > -1, x < 1

axiom log(1 + x) = SUM(n, 0, oo, (-1)^n * x ^ (n+1) / (n + 1)) for x > -1, x <= 1

axiom log(1 - x) = SUM(n, 0, oo, (-1)^n * (-x)^(n+1) / (n + 1)) for x >= -1, x < 1

# Common series evaluations

axiom SUM(n, 0, oo, 1 / (n+1)^2) = (pi^2) / 6

axiom SUM(n, a, oo, x ^ n) = (x^a)/(1-x) for x > -1, x < 1

axiom SUM(n, 0, oo, (-1)^n / (n+1)^2) = (pi^2) / 12

axiom SUM(n, 0, oo, (-1)^n * x^(2*n+1) / factorial(2*n+1)) = sin(x)

axiom SUM(n, 0, oo, (-1)^n * x^(2*n) / factorial(2*n)) = cos(x)

axiom SUM(n, 0, oo, x^n/factorial(n)) = exp(x)

# Common identities

## Absolute value

axiom [simp] abs(x) = x for x: real, x >= 0

axiom [simp] abs(x) = -x for x: real, x <= 0

axiom [simp] abs(x * y) = abs(x) * abs(y)

axiom [simp] abs(x / y) = abs(x) / abs(y) for y != 0

axiom [simp] abs(x ^ n) = abs(x) ^ n for n: int

## Complex numbers

axiom [simp] abs(b * i + a) = sqrt(a^2 + b^2) for a b: real

axiom [simp] abs(-(b * i) + a) = sqrt(a^2 + b^2) for a b: real

axiom [simp] abs(b * i - a) = sqrt(a^2 + b^2) for a b: real

axiom [simp] abs(-(b * i) - a) = sqrt(a^2 + b^2) for a b: real

axiom [simp] abs(-a + i) = sqrt(a^2 + 1) for a: real

axiom [simp] abs(a + i) = sqrt(a^2 + 1) for a: real

axiom [simp] abs(-a - i) = sqrt(a^2 + 1) for a: real

axiom [simp] abs(a - i) = sqrt(a^2 + 1) for a: real

axiom [simp] abs(-(a*i)) = a for a: real

axiom [simp] abs(a*i) = a for a: real

axiom [simp] abs(i) = 1

axiom [simp] conj(-i + a) = a + i for a: real

axiom [simp] conj(i + a) = a - i for a: real

axiom [simp] conj(-(b * i) + a) = b * i + a for a b: real

axiom [simp] conj(-(b * i) - a) = b * i - a for a b: real

axiom [simp] conj((b * i) - a) = -b * i - a for a b: real

axiom [simp] conj((b * i) + a) = -b * i + a for a b: real

axiom [simp] conj(b * i) = -b * i for b: real

axiom [simp] conj(-(b * i)) = b * i for b: real

axiom [simp] conj(i) = -i

axiom [simp] conj(-i) = i

axiom [simp] conj(a) = a for a: real

## Euler's formula

axiom exp(i*x) = cos(x) + i * sin(x) for x: real

axiom sin(x) = (exp(i*x) - exp(-i * x)) / (2*i)

## Power

axiom [simp] 0 ^ x = 0 for x: real, x > 0

axiom [bidirectional] (a * b) ^ k = a ^ k * b ^ k for k: int, a != 0, b != 0

axiom [bidirectional] (a * b) ^ k = a ^ k * b ^ k for k: real, a > 0, b > 0

axiom [bidirectional] (a / b) ^ k = a ^ k / b ^ k for k: int, a != 0, b != 0

axiom [bidirectional] (a / b) ^ k = a ^ k / b ^ k for k: real, a > 0, b > 0

axiom [bidirectional] (-a) ^ k = (-1) ^ k * a ^ k for k: int, a != 0

axiom [bidirectional] a ^ x ^ y = a ^ y ^ x for a > 0, x y: real

axiom [bidirectional] a ^ x ^ y = a ^ y ^ x for a != 0, x y: int

axiom [simp, bidirectional] a ^ x ^ y = a ^ (x * y) for a > 0, x y: real

axiom [simp, bidirectional] a ^ x ^ y = a ^ (x * y) for a != 0, x: int, y: real

axiom [bidirectional] x ^ a * x ^ b = x ^ (a + b) for x > 0, a b: real

axiom [bidirectional] x ^ a * x ^ b = x ^ (a + b) for x != 0, a b: int

axiom [bidirectional] x ^ a / x ^ b = x ^ (a - b) for x > 0, a b: real

axiom [bidirectional] x ^ a / x ^ b = x ^ (a - b) for x != 0, a b: int

axiom [simp] (-1) ^ (2 * n) = 1 for n: int

axiom [simp] (-x) ^ (2 * n) = x ^ (2 * n) for n: int

axiom [bidirectional] a ^ (-x) = (1 / a) ^ x for a != 0
