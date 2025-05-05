# Common integrals

axiom (INT x. c) = c * x + SKOLEM_CONST(C)

axiom (INT x. x) = x ^ 2 / 2 + SKOLEM_CONST(C)

axiom (INT x. 1 / x) = log(abs(x)) + SKOLEM_CONST(C) for x != 0

axiom (INT x. x ^ n) = x ^ (n + 1) / (n + 1) + SKOLEM_CONST(C) for n != -1

axiom (INT x. 1 / x ^ n) = 1 / (-((n - 1) * (x ^ (n - 1)))) + SKOLEM_CONST(C) for x != 0, n != 1

axiom (INT x. sqrt(x)) = 2/3 * x ^ (3/2) + SKOLEM_CONST(C) for x >= 0

axiom (INT x. 1 / sqrt(x)) = 2 * sqrt(x) + SKOLEM_CONST(C) for x > 0

axiom (INT x. exp(x)) = exp(x) + SKOLEM_CONST(C)

axiom (INT x. sin(x)) = -cos(x) + SKOLEM_CONST(C)

axiom (INT x. cos(x)) = sin(x) + SKOLEM_CONST(C)

axiom (INT x. 1 / (x^2 + 1)) = arctan(x) + SKOLEM_CONST(C)

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

## Exponential and Logarithm

axiom [bidirectional] exp(a) ^ b = exp(a * b)

axiom exp(a + b) = exp(a) * exp(b)

axiom exp(a - b) = exp(a) * exp(-b)

axiom [simp] exp(0) = 1

axiom [simp] log(1) = 0

axiom [simp] exp(log(b)) = b for b > 0

axiom [simp] log(exp(b)) = b

axiom [bidirectional] log(a * b) = log(a) + log(b) for a > 0, b > 0

axiom [bidirectional] log(a / b) = log(a) - log(b) for a > 0, b > 0

axiom [simp] log(1 / x) = -log(x) for x > 0

axiom [simp] log(x ^ a) = a * log(x) for x > 0, a: real

## Trigonometric identities

### Function tables

axiom [simp] sin(0) = 0
axiom [simp] sin(pi/6) = 1/2
axiom [simp] sin(pi/4) = sqrt(2)/2
axiom [simp] sin(pi/3) = sqrt(3)/2
axiom [simp] sin(pi/2) = 1
axiom [simp] sin(2*pi/3) = sqrt(3)/2
axiom [simp] sin(3*pi/4) = sqrt(2)/2
axiom [simp] sin(5*pi/6) = 1/2
axiom [simp] sin(pi) = 0

axiom [simp] cos(0) = 1
axiom [simp] cos(pi/6) = sqrt(3)/2
axiom [simp] cos(pi/4) = sqrt(2)/2
axiom [simp] cos(pi/3) = 1/2
axiom [simp] cos(pi/2) = 0
axiom [simp] cos(2*pi/3) = -1/2
axiom [simp] cos(3*pi/4) = -(sqrt(2)/2)
axiom [simp] cos(5*pi/6) = -(sqrt(3)/2)
axiom [simp] cos(pi) = -1

axiom [simp] tan(0) = 0
axiom [simp] tan(pi/6) = sqrt(3)/3
axiom [simp] tan(pi/4) = 1
axiom [simp] tan(pi/3) = sqrt(3)
axiom [simp] tan(2*pi/3) = -sqrt(3)
axiom [simp] tan(3*pi/4) = -1
axiom [simp] tan(5*pi/6) = -(sqrt(3)/3)
axiom [simp] tan(pi) = 0

axiom [simp] cot(pi/6) = sqrt(3)
axiom [simp] cot(pi/4) = 1
axiom [simp] cot(pi/3) = sqrt(3)/3
axiom [simp] cot(pi/2) = 0
axiom [simp] cot(2*pi/3) = -(sqrt(3)/3)
axiom [simp] cot(3*pi/4) = -1
axiom [simp] cot(5*pi/6) = -sqrt(3)

axiom [simp] csc(pi/6) = 2
axiom [simp] csc(pi/4) = sqrt(2)
axiom [simp] csc(pi/3) = 2*sqrt(3)/3
axiom [simp] csc(pi/2) = 1
axiom [simp] csc(2*pi/3) = 2*sqrt(3)/3
axiom [simp] csc(3*pi/4) = sqrt(2)
axiom [simp] csc(5*pi/6) = 2

axiom [simp] sec(0) = 1
axiom [simp] sec(pi/6) = 2*sqrt(3)/3
axiom [simp] sec(pi/4) = sqrt(2)
axiom [simp] sec(pi/3) = 2
axiom [simp] sec(2*pi/3) = -2
axiom [simp] sec(3*pi/4) = -sqrt(2)
axiom [simp] sec(5*pi/6) = -(2*sqrt(3)/3)
axiom [simp] sec(pi) = -1

axiom [simp] arcsin(-(sqrt(3)/2)) = -(pi/3)
axiom [simp] arcsin(-(sqrt(2)/2)) = -(pi/4)
axiom [simp] arcsin(-1) = -(pi/2)
axiom [simp] arcsin(-1/2) = -(pi/6)
axiom [simp] arcsin(0) = 0
axiom [simp] arcsin(1/2) = pi/6
axiom [simp] arcsin(1) = pi/2
axiom [simp] arcsin(sqrt(2)/2) = pi/4
axiom [simp] arcsin(sqrt(3)/2) = pi/3

axiom [simp] arccos(-(sqrt(3)/2)) = 5*pi/6
axiom [simp] arccos(-(sqrt(2)/2)) = 3*pi/4
axiom [simp] arccos(-1) = pi
axiom [simp] arccos(-1/2) = 2*pi/3
axiom [simp] arccos(0) = pi/2
axiom [simp] arccos(1/2) = pi/3
axiom [simp] arccos(1) = 0
axiom [simp] arccos(sqrt(2)/2) = pi/4
axiom [simp] arccos(sqrt(3)/2) = pi/6

axiom [simp] arctan(-sqrt(3)) = -(pi/3)
axiom [simp] arctan(-(sqrt(3)/3)) = -(pi/6)
axiom [simp] arctan(-1) = -(pi/4)
axiom [simp] arctan(0) = 0
axiom [simp] arctan(1) = pi/4
axiom [simp] arctan(sqrt(3)/3) = pi/6
axiom [simp] arctan(sqrt(3)) = pi/3

axiom [simp] arccot(-sqrt(3)) = 5*pi/6
axiom [simp] arccot(-1) = 3*pi/4
axiom [simp] arccot(-(sqrt(3)/3)) = 2*pi/3
axiom [simp] arccot(0) = pi/2
axiom [simp] arccot(sqrt(3)/3) = pi/3
axiom [simp] arccot(1) = pi/4
axiom [simp] arccot(sqrt(3)) = pi/6

axiom [simp] arccsc(-2) = -(pi/6)
axiom [simp] arccsc(-sqrt(2)) = -(pi/4)
axiom [simp] arccsc(-(2*sqrt(3)/3)) = -(pi/3)
axiom [simp] arccsc(-1) = -(pi/2)
axiom [simp] arccsc(1) = pi/2
axiom [simp] arccsc(2*sqrt(3)/3) = pi/3
axiom [simp] arccsc(sqrt(2)) = pi/4
axiom [simp] arccsc(2) = pi/6

axiom [simp] arcsec(-2) = 2*pi/3
axiom [simp] arcsec(-sqrt(2)) = 3*pi/4
axiom [simp] arcsec(-(2*sqrt(3)/3)) = 5*pi/6
axiom [simp] arcsec(-1) = pi
axiom [simp] arcsec(1) = 0
axiom [simp] arcsec(2*sqrt(3)/3) = pi/6
axiom [simp] arcsec(sqrt(2)) = pi/4
axiom [simp] arcsec(2) = pi/3

### Simple relations between trigonometric functions

axiom [simp] sin(-u) = -sin(u)

axiom [simp] cos(-u) = cos(u)

axiom [simp] tan(-u) = -tan(u)

axiom [simp] cot(-u) = -cot(u)

axiom [simp] sec(-u) = sec(u)

axiom [simp] csc(-u) = -csc(u)

axiom [simp] sin(pi / 2 - u) = cos(u)

axiom [simp] cos(pi / 2 - u) = sin(u)

axiom [simp] sin(-u + pi) = sin(u)

axiom [simp] cos(-u + pi) = -cos(u)

### Identity for sin(x)^2 + cos(x)^2 and variations

axiom [simp] sin(x)^2 + cos(x)^2 = 1

axiom [bidirectional] cot(x)^2 = csc(x)^2 - 1

axiom [bidirectional] cot(x)^2 = -1 + csc(x)^2

axiom [bidirectional] a * cot(x)^2 = a * csc(x)^2 - a

axiom [bidirectional] a * cot(x)^2 = -a + a * csc(x)^2

axiom [bidirectional] tan(x)^2 = sec(x)^2 - 1

axiom [bidirectional] tan(x)^2 = -1 + sec(x)^2

axiom [bidirectional] a * tan(x)^2 = a * sec(x)^2 - a

axiom [bidirectional] a * tan(x)^2 = -a + a * sec(x)^2

axiom [bidirectional] sin(x)^2 = 1 - cos(x)^2

axiom [bidirectional] sin(x)^2 = -cos(x)^2 + 1

axiom [bidirectional] a * sin(x)^2 = a - a * cos(x)^2

axiom [bidirectional] a * sin(x)^2 = -(a * cos(x)^2) + a

axiom [bidirectional] cos(x)^2 = 1 - sin(x)^2

axiom [bidirectional] cos(x)^2 = -sin(x)^2 + 1

axiom [bidirectional] a * cos(x)^2 = a - a * sin(x)^2

axiom [bidirectional] a * cos(x)^2 = -(a * sin(x)^2) + a

axiom [bidirectional] sec(x)^2 = tan(x)^2 + 1

axiom [bidirectional] sec(x)^2 = 1 + tan(x)^2

axiom [bidirectional] a * sec(x)^2 = a * tan(x)^2 + a

axiom [bidirectional] a * sec(x)^2 = a + a * tan(x)^2

### Double-angle formulas and variations

axiom sin(x)^2 = 1/2 * (1 - cos(2*x))

axiom cos(x)^2 = 1/2 * (1 + cos(2*x))

axiom sin(x)^2 = (1 - cos(2*x)) / 2

axiom cos(x)^2 = (1 + cos(2*x)) / 2

axiom [bidirectional] sin(2*x) = 2 * sin(x) * cos(x)

axiom [bidirectional] sin(x) = 2 * sin(x/2) * cos(x/2)

axiom cos(2*x) = 2 * cos(x)^2 - 1

axiom cos(x) = cos(x/2)^2 - sin(x/2)^2

axiom cos(2*x) = 1 - 2 * sin(x) ^ 2

axiom cos(x) = 1 - 2 * sin(x/2) ^ 2

axiom cos(2*x) = cos(x)^2 - sin(x)^2

### Sum to product formulas and variations

axiom sin(a) + sin(b) = 2 * sin((a + b) / 2) * cos((a - b) / 2)

axiom sin(a) - sin(b) = 2 * cos((a + b) / 2) * sin((a - b) / 2)

axiom cos(a) + cos(b) = 2 * cos((a + b) / 2) * cos((a - b) / 2)

axiom cos(a) - cos(b) = -2 * sin((a + b) / 2) * sin((a - b) / 2)

axiom cos(a) * sin(b) = (1 / 2) * (sin(a + b) - sin(a - b))

axiom sin(a) * cos(b) = (1 / 2) * (sin(a + b) + sin(a - b))

axiom cos(a) * cos(b) = (1 / 2) * (cos(a - b) + cos(a + b))

axiom sin(a) * sin(b) = -(1 / 2) * (cos(a - b) - cos(a + b))

axiom sin(a + b) = sin(a) * cos(b) + cos(a) * sin(b)

axiom sin(a - b) = sin(a) * cos(b) - cos(a) * sin(b)

axiom cos(a + b) = cos(a) * cos(b) - sin(a) * sin(b)

axiom cos(a - b) = cos(a) * cos(b) + sin(a) * sin(b)

### Relation with inverse trigonometric functions

axiom [simp] arcsin(sin(x)) = x for x >= -pi/2, x <= pi/2

axiom [simp] sin(arccos(x)) = sqrt(1-x^2)

axiom [simp] cos(arcsin(x)) = sqrt(1-x^2)

axiom [simp] tan(arcsec(x)) = sqrt(x ^ 2 - 1)

axiom [simp] tan(arcsin(x)) = x / sqrt(1 - x ^ 2)

axiom [simp] cos(arctan(x)) = 1 / sqrt(x ^ 2 + 1)

axiom [simp] sin(arctan(x)) = x / sqrt(x ^ 2 + 1)

axiom [simp] sec(arctan(x)) = sqrt(x ^ 2 + 1)

axiom [simp] cot(arctan(x)) = 1 / x

axiom [simp] csc(arctan(x)) = sqrt(x^2 + 1) / x

axiom [simp] cot(arcsin(x)) = sqrt(1 - x^2) / x

axiom [simp] csc(arcsin(x)) = 1 / x

axiom [simp] sin(arcsec(x)) = sqrt(x^2 - 1) / x

### Other identities

axiom arctan(a) - arctan(b) = arctan((a - b) / (1 + a * b))

axiom arctan(sqrt((1 - x) / (1 + x))) = arccos(x) / 2

axiom tan(a - b) = (tan(a) - tan(b)) / (1 + tan(a) * tan(b))

axiom [simp] cos(2 * arctan(z)) = (1 - z^2) / (1 + z^2)

axiom [simp] sin(2 * arctan(z)) = (2 * z) / (1 + z^2)

axiom arctan(x ^ -1) = pi/2 - arctan(x) for x != 0

axiom [simp] sin(2 * arcsin(x)) = 2 * x * sqrt(1 - x^2)

axiom arctan(-x) = -arctan(x)

axiom [bidirectional] 1 + sin(x) = (sin(x/2) + cos(x/2)) ^ 2

axiom [bidirectional] sin(x) + 1 = (sin(x/2) + cos(x/2)) ^ 2

axiom [bidirectional] sin(x) + cos(x) = sqrt(2) * sin(x+pi/4)

## Euler's Formula and variations

axiom [bidirectional] sin(x) = (exp(i*x) - exp(-i*x)) / (2*i) for x: real

axiom [bidirectional] exp(i*x) = cos(x) + i * sin(x) for x: real

axiom sin(x)^(2*n-1) = 1/(2^(2*n-2)) * SUM(k, 0, n-1, (-1)^(n+k-1) * binom(2*n-1, k) * sin((2*n - 2*k - 1) * x))

axiom sin(x)^(2*n) = 1/(2^(2*n)) * binom(2*n, n) + 1/(2^(2*n)) * SUM(k, 0, n-1, (-1)^(n-k) * 2 * binom(2*n, k) * cos(2*(n-k)*x))

## Factorial and binomial coefficient

axiom (m + 1) * factorial(m) = factorial(m + 1)

axiom [bidirectional] m * factorial(m - 1) = factorial(m)

define binom(n, m) = factorial(n) / (factorial(m) * factorial(n-m))

axiom binom(2*k+2, k+1) = 2 * binom(2*k, k) * ((2*k+1) / (k+1))

axiom (x + y) ^ n = SUM(k, 0, n, binom(n,k) * x^k * y^(n-k))

## Discrete functions

axiom [simp] sgn(x) = 1 for x > 0

axiom [simp] sgn(x) = -1 for x < 0

axiom [simp] sgn(0) = 0

axiom [simp] sgn(a) = 0 for a = 0

## Hyperbolic functions

define cosh(x) = (exp(x) + exp(-x)) / 2

define sinh(x) = (exp(x) - exp(-x)) / 2

## Results from contour integration

// Inside interesting Integrals, Section 8.10, C8.2

axiom (INT x:[0, oo]. sin(m * x) / (x * (a ^ 2 + x ^ 2))) = (pi * (1 - exp(-a * m))) / (2 * a ^ 2) for a > 0, m > 0

// Inside interesting Integrals, Section 3.1.7

axiom (INT x:[0, oo]. cos(a * x) / (x ^ 2 + b ^ 2)) = (pi / (2 * b)) * exp(-a * b) for a > 0, b > 0 

### Splitting rules for summation and product

axiom [split] SUM(n, l, u, f(n)) = f(l) + SUM(n, l+1, u, f(n)) for l < u

axiom [split] SUM(n, l, u, f(n)) = SUM(n, l, u-1, f(n)) + f(u) for l < u

axiom [split] SUM(n, l, u, f(n)) = SUM(n, l, a-1, f(n)) + SUM(n, a, u, f(n)) for a > l, a <= u

axiom [split] SUM(n, l, u, f(n)) = SUM(n, l, a, f(n)) + SUM(n, a+1, u, f(n)) for a >= l, a < u

axiom [split] SUM(n, l, u, f(n)) = SUM(n, l, a-1, f(n)) + f(a) + SUM(n, a+1, u, f(n)) for a > l, a < u

axiom [split] SUM(n, 0, oo, f(n)) = SUM(n, 0, oo, f(2*n+1)) + SUM(n, 0, oo, f(2*n))

axiom [split] MUL(n, l, u, f(n)) = MUL(n, l, u-1, f(n)) * f(u) for l < u

axiom [split] MUL(n, l, u, f(n)) = f(l) * MUL(n, l+1, u, f(n)) for l < u
