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

axiom (INT x. 1 / (x - i)) = log(x - i) + SKOLEM_CONST(C)

axiom (INT x:[a,b]. 1 / (x - i)) = [log(x - i)]_x=a,b for a b: real

axiom (INT x:[a,b]. 1 / (c * x - i)) = 1 / c * [log(c * x - i)]_x=a,b for a b c: real, c != 0

axiom (INT x:[a,b]. 1 / (x - d * i)) = [log(x - d * i)]_x=a,b for a b d: real

axiom (INT x:[a,b]. 1 / (c * x - d * i)) = 1 / c * [log(c * x - d * i)]_x=a,b for a b c d: real

axiom (INT x. 1 / (a * (x - i))) = 1 / a * log(x - i) + SKOLEM_CONST(C) for a: real, a != 0

axiom (INT x:[a,b]. 1 / (c * (x - i))) = 1 / c * [log(x - i)]_x=a,b for a b c: real, c != 0

axiom (INT x:[a,b]. 1 / (d * (c * x - i))) = 1 / (d * c) * [log(c * x - i)]_x=a,b for a b c d: real, c != 0, d != 0

axiom (INT x:[a,b]. 1 / (c * (x - d * i))) = 1 / c * [log(x - d * i)]_x=a,b for a b c d: real, c != 0

axiom (INT x:[a,b]. 1 / (e * (c * x - d * i))) = 1 / c * [log(c * x - d * i)]_x=a,b for a b c d e: real, e != 0

axiom (INT x. 1 / (x + i)) = log(x + i) + SKOLEM_CONST(C)

axiom (INT x:[a,b]. 1 / (x + i)) = [log(x + i)]_x=a,b for a b: real

axiom (INT x:[a,b]. 1 / (c * x + i)) = 1 / c * [log(c * x + i)]_x=a,b for a b c: real, c != 0

axiom (INT x:[a,b]. 1 / (x + d * i)) = [log(x + d * i)]_x=a,b for a b d: real

axiom (INT x:[a,b]. 1 / (c * x + d * i)) = 1 / c * [log(c * x + d * i)]_x=a,b for a b c d: real

axiom (INT x. 1 / (a * (x + i))) = 1 / a * log(x + i) + SKOLEM_CONST(C) for a: real, a != 0

axiom (INT x:[a,b]. 1 / (c * (x + i))) = 1 / c * [log(x + i)]_x=a,b for a b c: real, c != 0

axiom (INT x:[a,b]. 1 / (d * (c * x + i))) = 1 / (d * c) * [log(c * x + i)]_x=a,b for a b c d: real, c != 0, d != 0

axiom (INT x:[a,b]. 1 / (c * (x + d * i))) = 1 / c * [log(x + d * i)]_x=a,b for a b c d: real, c != 0

axiom (INT x:[a,b]. 1 / (e * (c * x + d * i))) = 1 / c * [log(c * x + d * i)]_x=a,b for a b c d e: real, e != 0

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

### Simple relations between trigonometric functions

axiom sin(-u) = -sin(u)

axiom cos(-u) = cos(u)

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

