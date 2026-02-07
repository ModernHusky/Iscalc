imports standard3
imports sqrt_integrals

# Indefinite Integrals of Inverse Trigonometric Functions
# Handbook of mathematical formulas and integrals
# Page 225

## 10.1 INTEGRANDS INVOLVING POWERS OF x AND POWERS OF INVERSE TRIGONOMETRIC FUNCTIONS
### 10.1.1 Integrands Involving x ^ n * arcsin(x/a) ^ m

// find an error in the book, the original answer is `arcsin(x/a) + sqrt(a ^ 2 - x ^ 2) + SKOLEM_CONST(C)`
// 1
prove (INT x. arcsin(x / a)) = x * arcsin(x / a) + sqrt(a ^ 2 - x ^ 2) + SKOLEM_CONST(C) for abs(x / a) < 1, a != 0, a ^ 2 - x ^ 2 >= 0
lhs:
    integrate by parts with u = arcsin(x / a), v = x
    simplify
    rewrite x / sqrt(1 - x ^ 2 / a ^ 2) to x / sqrt((a ^ 2 - x ^ 2) / a ^ 2)
    rewrite x / sqrt((a ^ 2 - x ^ 2) / a ^ 2) to (a * x) / sqrt(a ^ 2 - x ^ 2)
    rewrite 1 / a * (INT x. (a * x) / sqrt(a ^ 2 - x ^ 2)) to (INT x. x / sqrt(a ^ 2 - x ^ 2))
    substitute u for a ^ 2 - x ^ 2
    simplify
    rewrite 1 / sqrt(u) to u ^ (-1/2)
    apply integral identity
    replace substitution
    simplify
done

// 2
prove (INT x. arcsin(x / a) ^ 2) = x * arcsin(x / a) ^ 2 + 2 * sqrt(a ^ 2 - x ^ 2) * arcsin(x / a) - 2 * x + SKOLEM_CONST(C) for abs(x / a) < 1, a > 0, a ^ 2 - x ^ 2 >= 0
lhs:
    integrate by parts with u = arcsin(x / a) ^ 2, v = x
    simplify
    rewrite 1 - x ^ 2 / a ^ 2 to (a ^ 2 - x ^ 2) / a ^ 2
    rewrite x / sqrt((a ^ 2 - x ^ 2) / a ^ 2) to a * x / sqrt(a ^ 2 - x ^ 2)
    rewrite 2 / a * (INT x. a * x / sqrt(a ^ 2 - x ^ 2) * arcsin(x / a)) to 2 * (INT x. x / sqrt(a ^ 2 - x ^ 2) * arcsin(x / a))
    integrate by parts with u = arcsin(x / a), v = -sqrt(a ^ 2 - x ^ 2)
    simplify
    rewrite sqrt(a ^ 2 - x ^ 2) / sqrt(1 - x ^ 2 / a ^ 2) to sqrt(a ^ 2 - x ^ 2) / sqrt((a ^ 2 - x ^ 2) / a ^ 2)
    rewrite sqrt(a ^ 2 - x ^ 2) / sqrt((a ^ 2 - x ^ 2) / a ^ 2) to sqrt(a ^ 2 - x ^ 2) * sqrt(a ^ 2 / (a ^ 2 - x ^ 2))
    rewrite sqrt(a ^ 2 - x ^ 2) * sqrt(a ^ 2 / (a ^ 2 - x ^ 2)) to sqrt(a ^ 2)
    rewrite sqrt(a ^ 2) to abs(a)
    rewrite abs(a) to a
    apply integral identity
    simplify
done

// 3
prove (INT x. arcsin(x / a) ^ 3) = x * arcsin(x / a) ^ 3 + 3 * sqrt(a ^ 2 - x ^ 2) * arcsin(x / a) ^ 2 - 6 * x * arcsin(x / a) - 6 * sqrt(a ^ 2 - x ^ 2) + SKOLEM_CONST(C) for abs(x / a) < 1, a > 0, a ^ 2 - x ^ 2 >= 0
lhs:
    integrate by parts with u = arcsin(x / a) ^ 3, v = x
    simplify
    rewrite 1 - x ^ 2 / a ^ 2 to (a ^ 2 - x ^ 2) / a ^ 2
    rewrite x / sqrt((a ^ 2 - x ^ 2) / a ^ 2) to a * x / sqrt(a ^ 2 - x ^ 2)
    rewrite 3 / a * (INT x. a * x / sqrt(a ^ 2 - x ^ 2) * arcsin(x / a) ^ 2) to 3 * (INT x. x / sqrt(a ^ 2 - x ^ 2) * arcsin(x / a) ^ 2)
    integrate by parts with u = arcsin(x / a) ^ 2, v = -sqrt(a ^ 2 - x ^ 2)
    simplify
    rewrite sqrt(a ^ 2 - x ^ 2) / sqrt(1 - x ^ 2 / a ^ 2) to sqrt(a ^ 2 - x ^ 2) / sqrt((a ^ 2 - x ^ 2) / a ^ 2)
    rewrite sqrt(a ^ 2 - x ^ 2) / sqrt((a ^ 2 - x ^ 2) / a ^ 2) to sqrt(a ^ 2 - x ^ 2) * sqrt(a ^ 2 / (a ^ 2 - x ^ 2))
    rewrite sqrt(a ^ 2 - x ^ 2) * sqrt(a ^ 2 / (a ^ 2 - x ^ 2)) to sqrt(a ^ 2)
    rewrite sqrt(a ^ 2) to abs(a)
    rewrite abs(a) to a
    rewrite 6 / a * (INT x. a * arcsin(x / a)) to 6 * (INT x. arcsin(x / a))
    apply integral identity
    simplify
done

// 4
prove (INT x. arcsin(x / a) ^ n) = x * arcsin(x / a) ^ n + n * sqrt(a ^ 2 - x ^ 2) * arcsin(x / a) ^ (n - 1) - n * (n - 1) * (INT x. arcsin(x / a) ^ (n - 2)) for abs(x / a) <= 1, n != 1, n >= 2, isInt(n), a > 0, a ^ 2 - x ^ 2 >= 0
lhs:
    integrate by parts with u = arcsin(x / a) ^ n, v = x
    simplify
    rewrite 1 - x ^ 2 / a ^ 2 to (a ^ 2 - x ^ 2) / a ^ 2
    rewrite x / sqrt((a ^ 2 - x ^ 2) / a ^ 2) to a * x / sqrt(a ^ 2 - x ^ 2)
    rewrite n / a * (INT x. a * x / sqrt(a ^ 2 - x ^ 2) * arcsin(x / a) ^ (n - 1)) to n * (INT x. x / sqrt(a ^ 2 - x ^ 2) * arcsin(x / a) ^ (n - 1))
    integrate by parts with u = arcsin(x / a) ^ (n - 1), v = -sqrt(a ^ 2 - x ^ 2)
    simplify
    rewrite sqrt(a ^ 2 - x ^ 2) / sqrt(1 - x ^ 2 / a ^ 2) to sqrt(a ^ 2 - x ^ 2) / sqrt((a ^ 2 - x ^ 2) / a ^ 2)
    rewrite sqrt(a ^ 2 - x ^ 2) / sqrt((a ^ 2 - x ^ 2) / a ^ 2) to sqrt(a ^ 2 - x ^ 2) * sqrt(a ^ 2 / (a ^ 2 - x ^ 2))
    rewrite sqrt(a ^ 2 - x ^ 2) * sqrt(a ^ 2 / (a ^ 2 - x ^ 2)) to sqrt(a ^ 2)
    rewrite sqrt(a ^ 2) to abs(a)
    rewrite abs(a) to a
    rewrite (n - 1) / a * (INT x. a * arcsin(x / a) ^ (n - 2)) to (n - 1) * (INT x. arcsin(x / a) ^ (n - 2))
    expand polynomial
    simplify
    rewrite n * sqrt(a ^ 2 - x ^ 2) * arcsin(x / a) ^ (n - 1) - n ^ 2 * (INT x. arcsin(x / a) ^ (n - 2)) + n * (INT x. arcsin(x / a) ^ (n - 2)) + x * arcsin(x / a) ^ n to x * arcsin(x / a) ^ n + n * sqrt(a ^ 2 - x ^ 2) * arcsin(x / a) ^ (n - 1) + (n - n ^ 2) * (INT x. arcsin(x / a) ^ (n - 2))
    rewrite n - n ^ 2 to -n * (n - 1)
    rewrite x * arcsin(x / a) ^ n + n * sqrt(a ^ 2 - x ^ 2) * arcsin(x / a) ^ (n - 1) + -n * (n - 1) * (INT x. arcsin(x / a) ^ (n - 2)) to x * arcsin(x / a) ^ n + n * sqrt(a ^ 2 - x ^ 2) * arcsin(x / a) ^ (n - 1) - n * (n - 1) * (INT x. arcsin(x / a) ^ (n - 2))
done

// 5
prove (INT x. x * arcsin(x / a)) = (x ^ 2 / 2 -  a ^ 2 / 4) * arcsin(x / a) + x / 4 * sqrt(a ^ 2 - x ^ 2) + SKOLEM_CONST(C) for abs(x / a) <= 1, a != 0, a ^ 2 - x ^ 2 >= 0, sqrt(a ^ 2 - x ^ 2) != 0
lhs:
    integrate by parts with u = arcsin(x / a), v = x ^ 2 / 2
    simplify
    rewrite 1 - x ^ 2 / a ^ 2 to (a ^ 2 - x ^ 2) / a ^ 2
    rewrite x ^ 2 / sqrt((a ^ 2 - x ^ 2) / a ^ 2) to x ^ 2 * a / sqrt(a ^ 2 - x ^ 2)
    simplify
    apply integral identity
    simplify
    rewrite x * sqrt(a ^ 2 - x ^ 2) / 4 - a ^ 2 / 4 * arcsin(x / a) + x ^ 2 / 2 * arcsin(x / a) to (x ^ 2 / 2 - a ^ 2 / 4) * arcsin(x / a) + x / 4 * sqrt(a ^ 2 - x ^ 2)
done

// 6
prove (INT x. x ^ 2 * arcsin(x / a)) = x ^ 3 / 3 * arcsin(x / a) + 1 / 9 * (x ^ 2 + 2 * a ^ 2) * sqrt(a ^ 2 - x ^ 2) + SKOLEM_CONST(C) for abs(x / a) <= 1, a != 0, a ^ 2 - x ^ 2 >= 0, sqrt(a ^ 2 - x ^ 2) != 0
lhs:
    integrate by parts with u = arcsin(x / a), v = x ^ 3 / 3
    simplify
    rewrite 1 - x ^ 2 / a ^ 2 to (a ^ 2 - x ^ 2) / a ^ 2
    rewrite x ^ 3 / sqrt((a ^ 2 - x ^ 2) / a ^ 2) to x ^ 3 * a / sqrt(a ^ 2 - x ^ 2)
    simplify
    apply integral identity
    simplify
    rewrite sqrt(a ^ 2 - x ^ 2) * (-(2 * a ^ 2) - x ^ 2) / 9 to sqrt(a ^ 2 - x ^ 2) * (x ^ 2 + 2 * a ^ 2) / 9
    rewrite x ^ 3 / 3 * arcsin(x / a) - sqrt(a ^ 2 - x ^ 2) * (x ^ 2 + 2 * a ^ 2) / 9 + SKOLEM_CONST(C) to x ^ 3 / 3 * arcsin(x / a) + 1/9 * (x ^ 2 + 2 * a ^ 2) * sqrt(a ^ 2 - x ^ 2) + SKOLEM_CONST(C)
done

// 7
prove (INT x. x ^ 3 * arcsin(x / a)) = (x ^ 4 / 4 -  3 * a ^ 4 / 32) * arcsin(x / a) + 1 / 32 * (2 * x ^ 3 + 3 * a ^ 2 * x) * sqrt(a ^ 2 - x ^ 2) + SKOLEM_CONST(C) for abs(x / a) <= 1, a > 0, a ^ 2 - x ^ 2 >= 0, sqrt(a ^ 2 - x ^ 2) != 0
lhs:
    integrate by parts with u = arcsin(x / a), v = x ^ 4 / 4
    simplify
    rewrite 1 - x ^ 2 / a ^ 2 to (a ^ 2 - x ^ 2) / a ^ 2
    rewrite x ^ 4 / sqrt((a ^ 2 - x ^ 2) / a ^ 2) to x ^ 4 * a / sqrt(a ^ 2 - x ^ 2)
    simplify
    apply integral identity
    simplify
    rewrite 3 * a ^ 2 * x * sqrt(a ^ 2 - x ^ 2) / 32 + x ^ 3 * sqrt(a ^ 2 - x ^ 2) / 16 - 3 * a ^ 4 / 32 * arcsin(x / a) + x ^ 4 / 4 * arcsin(x / a) to (x ^ 4 / 4 - 3 * a ^ 4 / 32) * arcsin(x / a) + 1 / 32 * (2 * x ^ 3 + 3 * a ^ 2 * x) * sqrt(a ^ 2 - x ^ 2)
done

// 8
prove (INT x. x ^ n * arcsin(x / a)) = x ^ (n + 1) / (n + 1) * arcsin(x / a) - 1 / (n + 1) * (INT x. x ^ (n + 1) / sqrt(a ^ 2 - x ^ 2)) for abs(x / a) < 1, a != 0, n >= 0, isInt(n), a ^ 2 - x ^ 2 > 0
lhs:
    integrate by parts with u = arcsin(x / a), v = x ^ (n + 1) / (n + 1)
    simplify
    rewrite 1 - x ^ 2 / a ^ 2 to (a ^ 2 - x ^ 2) / a ^ 2
    rewrite x ^ (n + 1) / sqrt((a ^ 2 - x ^ 2) / a ^ 2) to x ^ (n + 1) * a / sqrt(a ^ 2 - x ^ 2)
    simplify
done

## 10.1.2 Integrands Involving x ^ (-n) * arcsin(x/a)
#### 10.1.2.1

// 2
prove (INT x. 1 / x ^ 2 * arcsin(x / a)) = -1 / x * arcsin(x / a) - 1 / a * log(abs((a + sqrt(a ^ 2 - x ^ 2))/x)) + SKOLEM_CONST(C) for abs(x / a) < 1, a > 0, x > 0, a ^ 2 - x ^ 2 > 0
lhs:
    integrate by parts with u = arcsin(x / a), v = -1 / x
    simplify
    rewrite 1 - x ^ 2 / a ^ 2 to (a ^ 2 - x ^ 2) / a ^ 2
    rewrite 1 / (x * sqrt((a ^ 2 - x ^ 2) / a ^ 2)) to a / (x * sqrt(a ^ 2 - x ^ 2))
    simplify
    apply integral identity
done

// 3
prove (INT x. 1 / x ^ 3 * arcsin(x / a)) = -1 / (2 * x ^ 2 ) * arcsin(x / a) - sqrt(a ^ 2 - x ^ 2)/(2 * a ^ 2 * x) + SKOLEM_CONST(C) for abs(x / a) < 1, a > 0, x != 0, a ^ 2 - x ^ 2 > 0, sqrt(a ^ 2 - x ^ 2) != 0
lhs:
    integrate by parts with u = arcsin(x / a), v = -1 / (2 * x ^ 2)
    simplify
    rewrite 1 - x ^ 2 / a ^ 2 to (a ^ 2 - x ^ 2) / a ^ 2
    rewrite 1 / (x ^ 2 * sqrt((a ^ 2 - x ^ 2) / a ^ 2)) to a / (x ^ 2 * sqrt(a ^ 2 - x ^ 2))
    simplify
    apply integral identity
    simplify
done

// 4
prove (INT x. 1 / x ^ n * arcsin(x / a)) = -1 / ((n - 1) * x ^ (n - 1)) * arcsin(x / a) + 1 / (n - 1) * (INT x. 1 / (x ^ (n - 1) * sqrt(a ^ 2 - x ^ 2))) for abs(x / a) < 1, n >= 2, isInt(n), a != 0, x > 0, a ^ 2 - x ^ 2 > 0
lhs:
    rewrite 1 / x ^ n to x ^ (-n)
    integrate by parts with u = arcsin(x / a), v = x ^ (1 - n) / (1 - n)
    simplify
    rewrite 1 - x ^ 2 / a ^ 2 to (a ^ 2 - x ^ 2) / a ^ 2
    rewrite x ^ (1 - n) / sqrt((a ^ 2 - x ^ 2) / a ^ 2) to x ^ (1 - n) * a / sqrt(a ^ 2 - x ^ 2)
    simplify
    rewrite x ^ (1 - n) / (1 - n) * arcsin(x / a) - 1 / (1 - n) * (INT x. x ^ (1 - n) / sqrt(a ^ 2 - x ^ 2)) to x ^ (1 - n) / (1 - n) * arcsin(x / a) + 1 / (n - 1) * (INT x. x ^ (1 - n) / sqrt(a ^ 2 - x ^ 2))
    rewrite x ^ (1 - n) to 1 / x ^ (n - 1)
    rewrite 1 / x ^ (n - 1) / (1 - n) to -1 / (x ^ (n - 1) * (n - 1))
done

## 10.1.3 Integrands Involving x ^ n * arccos(x/a) ^ m
### 10.1.3.1

// 1
prove (INT x. arccos(x / a)) = x * arccos(x / a) - sqrt(a ^ 2 - x ^ 2) + SKOLEM_CONST(C) for abs(x / a) <= 1, a != 0, a ^ 2 - x ^ 2 >= 0
lhs:
    integrate by parts with u = arccos(x / a), v = x
    simplify
    rewrite x / sqrt(1 - x ^ 2 / a ^ 2) to x / sqrt((a ^ 2 - x ^ 2) / a ^ 2)
    rewrite x / sqrt((a ^ 2 - x ^ 2) / a ^ 2) to (a * x) / sqrt(a ^ 2 - x ^ 2)
    rewrite 1 / a * (INT x. (a * x) / sqrt(a ^ 2 - x ^ 2)) to (INT x. x / sqrt(a ^ 2 - x ^ 2))
    substitute u for a ^ 2 - x ^ 2
    simplify
    rewrite 1 / sqrt(u) to u ^ (-1/2)
    apply integral identity
    replace substitution
    simplify
done

// 2
prove (INT x. arccos(x / a) ^ 2) = x * arccos(x / a) ^ 2 - 2 * sqrt(a ^ 2 - x ^ 2) * arccos(x / a) - 2 * x + SKOLEM_CONST(C) for abs(x / a) <= 1, a > 0, a ^ 2 - x ^ 2 >= 0
lhs:
    integrate by parts with u = arccos(x / a) ^ 2, v = x
    simplify
    rewrite 1 - x ^ 2 / a ^ 2 to (a ^ 2 - x ^ 2) / a ^ 2
    rewrite x / sqrt((a ^ 2 - x ^ 2) / a ^ 2) to a * x / sqrt(a ^ 2 - x ^ 2)
    rewrite 2 / a * (INT x. a * x / sqrt(a ^ 2 - x ^ 2) * arccos(x / a)) to 2 * (INT x. x / sqrt(a ^ 2 - x ^ 2) * arccos(x / a))
    integrate by parts with u = arccos(x / a), v = -sqrt(a ^ 2 - x ^ 2)
    simplify
    rewrite sqrt(a ^ 2 - x ^ 2) / sqrt(1 - x ^ 2 / a ^ 2) to sqrt(a ^ 2 - x ^ 2) / sqrt((a ^ 2 - x ^ 2) / a ^ 2)
    rewrite sqrt(a ^ 2 - x ^ 2) / sqrt((a ^ 2 - x ^ 2) / a ^ 2) to sqrt(a ^ 2 - x ^ 2) * sqrt(a ^ 2 / (a ^ 2 - x ^ 2))
    rewrite sqrt(a ^ 2 - x ^ 2) * sqrt(a ^ 2 / (a ^ 2 - x ^ 2)) to sqrt(a ^ 2)
    rewrite sqrt(a ^ 2) to abs(a)
    rewrite abs(a) to a
    apply integral identity
    simplify
done

// 3
prove (INT x. arccos(x / a) ^ 3) = x * arccos(x / a) ^ 3 - 3 * sqrt(a ^ 2 - x ^ 2) * arccos(x / a) ^ 2 - 6 * x * arccos(x / a) + 6 * sqrt(a ^ 2 - x ^ 2) + SKOLEM_CONST(C) for abs(x / a) <= 1, a > 0, a ^ 2 - x ^ 2 >= 0
lhs:
    integrate by parts with u = arccos(x / a) ^ 3, v = x
    simplify
    rewrite 1 - x ^ 2 / a ^ 2 to (a ^ 2 - x ^ 2) / a ^ 2
    rewrite x / sqrt((a ^ 2 - x ^ 2) / a ^ 2) to a * x / sqrt(a ^ 2 - x ^ 2)
    rewrite 3 / a * (INT x. a * x / sqrt(a ^ 2 - x ^ 2) * arccos(x / a) ^ 2) to 3 * (INT x. x / sqrt(a ^ 2 - x ^ 2) * arccos(x / a) ^ 2)
    integrate by parts with u = arccos(x / a) ^ 2, v = -sqrt(a ^ 2 - x ^ 2)
    simplify
    rewrite sqrt(a ^ 2 - x ^ 2) / sqrt(1 - x ^ 2 / a ^ 2) to sqrt(a ^ 2 - x ^ 2) / sqrt((a ^ 2 - x ^ 2) / a ^ 2)
    rewrite sqrt(a ^ 2 - x ^ 2) / sqrt((a ^ 2 - x ^ 2) / a ^ 2) to sqrt(a ^ 2 - x ^ 2) * sqrt(a ^ 2 / (a ^ 2 - x ^ 2))
    rewrite sqrt(a ^ 2 - x ^ 2) * sqrt(a ^ 2 / (a ^ 2 - x ^ 2)) to sqrt(a ^ 2)
    rewrite sqrt(a ^ 2) to abs(a)
    rewrite abs(a) to a
    rewrite 6 / a * (INT x. a * arccos(x / a)) to 6 * (INT x. arccos(x / a))
    apply integral identity
    simplify
done

// 4
prove (INT x. arccos(x / a) ^ n) = x * arccos(x / a) ^ n -  n * sqrt(a ^ 2 - x ^ 2) * arccos(x / a) ^ (n - 1) - n * (n - 1) * (INT x. arccos(x / a) ^ (n - 2)) for abs(x / a) <= 1, n != 1, n >= 2, isInt(n), a > 0, a ^ 2 - x ^ 2 >= 0
lhs:
    integrate by parts with u = arccos(x / a) ^ n, v = x
    simplify
    rewrite 1 - x ^ 2 / a ^ 2 to (a ^ 2 - x ^ 2) / a ^ 2
    rewrite x / sqrt((a ^ 2 - x ^ 2) / a ^ 2) to a * x / sqrt(a ^ 2 - x ^ 2)
    rewrite n / a * (INT x. a * x / sqrt(a ^ 2 - x ^ 2) * arccos(x / a) ^ (n - 1)) to n * (INT x. x / sqrt(a ^ 2 - x ^ 2) * arccos(x / a) ^ (n - 1))
    integrate by parts with u = arccos(x / a) ^ (n - 1), v = -sqrt(a ^ 2 - x ^ 2)
    simplify
    rewrite sqrt(a ^ 2 - x ^ 2) / sqrt(1 - x ^ 2 / a ^ 2) to sqrt(a ^ 2 - x ^ 2) / sqrt((a ^ 2 - x ^ 2) / a ^ 2)
    rewrite sqrt(a ^ 2 - x ^ 2) / sqrt((a ^ 2 - x ^ 2) / a ^ 2) to sqrt(a ^ 2 - x ^ 2) * sqrt(a ^ 2 / (a ^ 2 - x ^ 2))
    rewrite sqrt(a ^ 2 - x ^ 2) * sqrt(a ^ 2 / (a ^ 2 - x ^ 2)) to sqrt(a ^ 2)
    rewrite sqrt(a ^ 2) to abs(a)
    rewrite abs(a) to a
    rewrite (n - 1) / a * (INT x. a * arccos(x / a) ^ (n - 2)) to (n - 1) * (INT x. arccos(x / a) ^ (n - 2))
    expand polynomial
    rewrite -(n * sqrt(a ^ 2 - x ^ 2) * arccos(x / a) ^ (n - 1)) - n ^ 2 * (INT x. arccos(x / a) ^ (n - 2)) + n * (INT x. arccos(x / a) ^ (n - 2)) + x * arccos(x / a) ^ n to x * arccos(x / a) ^ n - n * sqrt(a ^ 2 - x ^ 2) * arccos(x / a) ^ (n - 1) + (n - n ^ 2) * (INT x. arccos(x / a) ^ (n - 2))
    rewrite n - n ^ 2 to -n * (n - 1)
    rewrite x * arccos(x / a) ^ n - n * sqrt(a ^ 2 - x ^ 2) * arccos(x / a) ^ (n - 1) + -n * (n - 1) * (INT x. arccos(x / a) ^ (n - 2)) to x * arccos(x / a) ^ n - n * sqrt(a ^ 2 - x ^ 2) * arccos(x / a) ^ (n - 1) - n * (n - 1) * (INT x. arccos(x / a) ^ (n - 2))
done

// 5
prove (INT x. x * arccos(x / a)) = (x ^ 2 / 2 -  a ^ 2 / 4) * arccos(x / a) + x / 4 * sqrt(a ^ 2 - x ^ 2) + SKOLEM_CONST(C) for abs(x / a) < 1, a != 0, a ^ 2 - x ^ 2 >= 0, sqrt(a ^ 2 - x ^ 2) != 0
lhs:
    integrate by parts with u = arccos(x / a), v = x ^ 2 / 2
    simplify
    rewrite 1 - x ^ 2 / a ^ 2 to (a ^ 2 - x ^ 2) / a ^ 2
    rewrite x ^ 2 / sqrt((a ^ 2 - x ^ 2) / a ^ 2) to x ^ 2 * a / sqrt(a ^ 2 - x ^ 2)
    simplify
    apply integral identity
    simplify
    rewrite a ^ 2 / 4 * arcsin(x / a) - x * sqrt(a ^ 2 - x ^ 2) / 4 + x ^ 2 / 2 * arccos(x / a) + SKOLEM_CONST(C) to a ^ 2 / 4 * (pi / 2 - arccos(x / a)) - x * sqrt(a ^ 2 - x ^ 2) / 4 + x ^ 2 / 2 * arccos(x / a) + SKOLEM_CONST(C)
    expand polynomial
    simplify
    rewrite -(x * sqrt(a ^ 2 - x ^ 2) / 4) - a ^ 2 / 4 * arccos(x / a) + x ^ 2 / 2 * arccos(x / a) + a ^ 2 * pi / 8 + SKOLEM_CONST(C) to (x ^ 2 / 2 - a ^ 2 / 4) * arccos(x / a) + x / 4 * sqrt(a ^ 2 - x ^ 2) + (a ^ 2 * pi / 8 + SKOLEM_CONST(C))
    rewrite (a ^ 2 * pi / 8 + SKOLEM_CONST(C)) to SKOLEM_CONST(C)
done

// 6
prove (INT x. x ^ 2 * arccos(x / a)) = x ^ 3 / 3 * arccos(x / a) - 1 / 9 * (x ^ 2 + 2 * a ^ 2) * sqrt(a ^ 2 - x ^ 2) + SKOLEM_CONST(C) for abs(x / a) < 1, a != 0, a ^ 2 - x ^ 2 >= 0, sqrt(a ^ 2 - x ^ 2) != 0
lhs:
    integrate by parts with u = arccos(x / a), v = x ^ 3 / 3
    simplify
    rewrite 1 - x ^ 2 / a ^ 2 to (a ^ 2 - x ^ 2) / a ^ 2
    rewrite x ^ 3 / sqrt((a ^ 2 - x ^ 2) / a ^ 2) to x ^ 3 * a / sqrt(a ^ 2 - x ^ 2)
    simplify
    apply integral identity
    simplify
    rewrite sqrt(a ^ 2 - x ^ 2) * (-(2 * a ^ 2) - x ^ 2) / 9 + x ^ 3 / 3 * arccos(x / a) + SKOLEM_CONST(C) to x ^ 3 / 3 * arccos(x / a) - 1/9 * (x ^ 2 + 2 * a ^ 2) * sqrt(a ^ 2 - x ^ 2) + SKOLEM_CONST(C)
done

// 7
prove (INT x. x ^ 3 * arccos(x / a)) = (x ^ 4 / 4 -  3 * a ^ 4 / 32) * arccos(x / a) - 1 / 32 * (2 * x ^ 3 + 3 * a ^ 2 * x) * sqrt(a ^ 2 - x ^ 2) + SKOLEM_CONST(C) for abs(x / a) <= 1, a > 0, a ^ 2 - x ^ 2 >= 0, sqrt(a ^ 2 - x ^ 2) != 0
lhs:
    integrate by parts with u = arccos(x / a), v = x ^ 4 / 4
    simplify
    rewrite 1 - x ^ 2 / a ^ 2 to (a ^ 2 - x ^ 2) / a ^ 2
    rewrite x ^ 4 / sqrt((a ^ 2 - x ^ 2) / a ^ 2) to x ^ 4 * a / sqrt(a ^ 2 - x ^ 2)
    simplify
    apply integral identity
    simplify
    rewrite -(3 * a ^ 2 * x * sqrt(a ^ 2 - x ^ 2) / 32) - x ^ 3 * sqrt(a ^ 2 - x ^ 2) / 16 + 3 * a ^ 4 / 32 * arcsin(x / a) + x ^ 4 / 4 * arccos(x / a) + SKOLEM_CONST(C) to (x ^ 4 / 4 - 3 * a ^ 4 / 32) * arccos(x / a) - 1 / 32 * (2 * x ^ 3 + 3 * a ^ 2 * x) * sqrt(a ^ 2 - x ^ 2) + 3 * a ^ 4 / 32 * arcsin(x / a) + SKOLEM_CONST(C)
    rewrite (x ^ 4 / 4 - 3 * a ^ 4 / 32) * arccos(x / a) - 1 / 32 * (2 * x ^ 3 + 3 * a ^ 2 * x) * sqrt(a ^ 2 - x ^ 2) + 3 * a ^ 4 / 32 * arcsin(x / a) + SKOLEM_CONST(C) to (x ^ 4 / 4 - 3 * a ^ 4 / 32) * arccos(x / a) - 1 / 32 * (2 * x ^ 3 + 3 * a ^ 2 * x) * sqrt(a ^ 2 - x ^ 2) + 3 * a ^ 4 / 32 * (pi / 2 - arccos(x / a)) + SKOLEM_CONST(C)
    expand polynomial
    simplify
    rewrite -(3 * a ^ 2 * x * sqrt(a ^ 2 - x ^ 2) / 32) - x ^ 3 * sqrt(a ^ 2 - x ^ 2) / 16 - 3 * a ^ 4 / 16 * arccos(x / a) + x ^ 4 / 4 * arccos(x / a) + 3 * a ^ 4 * pi / 64 + SKOLEM_CONST(C) to (x ^ 4 / 4 - 3 * a ^ 4 / 16) * arccos(x / a) - 1 / 32 * (2 * x ^ 3 + 3 * a ^ 2 * x) * sqrt(a ^ 2 - x ^ 2) + (3 * a ^ 4 * pi / 64 + SKOLEM_CONST(C))
    rewrite (3 * a ^ 4 * pi / 64 + SKOLEM_CONST(C)) to SKOLEM_CONST(C)
    rewrite (x ^ 4 / 4 - 3 * a ^ 4 / 16) to (x ^ 4 / 4 - 3 * a ^ 4 / 32 - 3 * a ^ 4 / 32)
    simplify
    rewrite arccos(x / a) * (x ^ 4 / 4 - 3 * a ^ 4 / 16) - sqrt(a ^ 2 - x ^ 2) * (3 * a ^ 2 * x / 32 + x ^ 3 / 16) + SKOLEM_CONST(C) to (x ^ 4 / 4 - 3 * a ^ 4 / 32) * arccos(x / a) - 1 / 32 * (2 * x ^ 3 + 3 * a ^ 2 * x) * sqrt(a ^ 2 - x ^ 2) + SKOLEM_CONST(C)
done

// 8
prove (INT x. x ^ n * arccos(x / a)) = x ^ (n + 1) / (n + 1) * arccos(x / a) + 1 / (n + 1) * (INT x. x ^ (n + 1) / sqrt(a ^ 2 - x ^ 2)) for abs(x / a) < 1, a != 0, n >= 0, isInt(n), a ^ 2 - x ^ 2 > 0
lhs:
    integrate by parts with u = arccos(x / a), v = x ^ (n + 1) / (n + 1)
    simplify
    rewrite 1 - x ^ 2 / a ^ 2 to (a ^ 2 - x ^ 2) / a ^ 2
    rewrite x ^ (n + 1) / sqrt((a ^ 2 - x ^ 2) / a ^ 2) to x ^ (n + 1) * a / sqrt(a ^ 2 - x ^ 2)
    simplify
done

## 10.1.4 Integrands Involving x ^ (-n) * arccos(x/a)
### 10.1.4.1

// 2
prove (INT x. 1 / x ^ 2 * arccos(x / a)) = -1 / x * arccos(x / a) + 1 / a * log(abs((a + sqrt(a ^ 2 - x ^ 2)) / x)) + SKOLEM_CONST(C) for abs(x / a) < 1, a > 0, x > 0, a ^ 2 - x ^ 2 > 0, sqrt(a ^ 2 - x ^ 2) != 0
lhs:
    integrate by parts with u = arccos(x / a), v = -1 / x
    simplify
    rewrite 1 - x ^ 2 / a ^ 2 to (a ^ 2 - x ^ 2) / a ^ 2
    rewrite 1 / (x * sqrt((a ^ 2 - x ^ 2) / a ^ 2)) to a / (x * sqrt(a ^ 2 - x ^ 2))
    simplify
    apply integral identity
done

// 3
prove (INT x. 1 / x ^ 3 * arccos(x / a)) = -1 / (2 * x ^ 2) * arccos(x / a) + sqrt(a ^ 2 - x ^ 2) / (2 * a ^ 2 * x) + SKOLEM_CONST(C) for abs(x / a) < 1, a > 0, x != 0, a ^ 2 - x ^ 2 > 0, sqrt(a ^ 2 - x ^ 2) != 0
lhs:
    integrate by parts with u = arccos(x / a), v = -1 / (2 * x ^ 2)
    simplify
    rewrite 1 - x ^ 2 / a ^ 2 to (a ^ 2 - x ^ 2) / a ^ 2
    rewrite 1 / (x ^ 2 * sqrt((a ^ 2 - x ^ 2) / a ^ 2)) to a / (x ^ 2 * sqrt(a ^ 2 - x ^ 2))
    simplify
    apply integral identity
    simplify
done

// 4
prove (INT x. 1 / x ^ n * arccos(x / a)) = -1 / ((n - 1) * x ^ (n - 1)) * arccos(x / a) - 1 / (n - 1) * (INT x. 1 / (x ^ (n - 1) * sqrt(a ^ 2 - x ^ 2))) for abs(x / a) < 1, n >= 2, isInt(n), a != 0, x > 0, a ^ 2 - x ^ 2 > 0
lhs:
    rewrite 1 / x ^ n to x ^ (-n)
    integrate by parts with u = arccos(x / a), v = x ^ (1 - n) / (1 - n)
    simplify
    rewrite 1 - x ^ 2 / a ^ 2 to (a ^ 2 - x ^ 2) / a ^ 2
    rewrite x ^ (1 - n) / sqrt((a ^ 2 - x ^ 2) / a ^ 2) to x ^ (1 - n) * a / sqrt(a ^ 2 - x ^ 2)
    simplify
    rewrite 1 / (1 - n) * (INT x. x ^ (1 - n) / sqrt(a ^ 2 - x ^ 2)) + x ^ (1 - n) / (1 - n) * arccos(x / a) to x ^ (1 - n) / (1 - n) * arccos(x / a) - 1 / (n - 1) * (INT x. x ^ (1 - n) / sqrt(a ^ 2 - x ^ 2))
    rewrite x ^ (1 - n) to 1 / x ^ (n - 1)
    rewrite 1 / x ^ (n - 1) / (1 - n) to -1 / (x ^ (n - 1) * (n - 1))
done

## 10.1.5 Integrands Involving x ^ n * arctan(x/a)
### 10.1.5.1

// 1
prove (INT x. arctan(x / a)) = x * arctan(x / a) - a / 2 * log(x ^ 2 + a ^ 2) + SKOLEM_CONST(C) for a: real, a != 0, x != 0
lhs:
    integrate by parts with u = arctan(x / a), v = x
    simplify
    rewrite x / (x^2 / a^2 + 1) to (x * a^2) / (x^2 + a^2)
    rewrite (x * a^2) / (x^2 + a^2) to (a / 2) * (2 * x * a) / (x^2 + a^2)
    substitute u for x^2 + a^2
    apply integral identity
    replace substitution
    rewrite log(abs(x^2 + a^2)) to log(x^2 + a^2)
done

// 2
prove (INT x. x * arctan(x / a)) = 1 / 2 * (x ^ 2 + a ^ 2) * arctan(x / a) - 1 / 2 * a * x + SKOLEM_CONST(C) for a != 0, x != 0
lhs:
    integrate by parts with u = arctan(x / a), v = x ^ 2 / 2
    simplify
    rewrite x ^ 2 / a ^ 2 + 1 to (x ^ 2 + a ^ 2) / a ^ 2
    rewrite x ^ 2 / ((x ^ 2 + a ^ 2) / a ^ 2) to a ^ 2 * x ^ 2 / (x ^ 2 + a ^ 2)
    rewrite a ^ 2 * x ^ 2 / (x ^ 2 + a ^ 2) to a ^ 2 * (1 - a ^ 2 / (x ^ 2 + a ^ 2))
    simplify
    apply integral identity
    simplify
    rewrite x ^ 2 / 2 * arctan(x / a) - a / 2 * (x - a * arctan(x / a)) to (a^2 / 2 + x^2 / 2) * arctan(x / a) - a * x / 2
    simplify
done

// 3
prove (INT x. x ^ 2 * arctan(x / a)) = 1 / 3 * x ^3 * arctan(x / a) - 1 / 6 * a * x ^ 2 + 1 / 6 * a ^ 3 * log(x ^ 2 + a ^ 2) + SKOLEM_CONST(C) for a: real, a != 0, x != 0
lhs:
    integrate by parts with u = arctan(x / a), v = x^3 / 3
    simplify
    rewrite x ^ 3 / (x ^ 2 / a ^ 2 + 1) to (a ^ 2 * x ^ 3) / (x ^ 2 + a ^ 2)
    simplify
    rewrite x ^ 3 / (a ^ 2 + x ^ 2) to x - (a ^ 2 * x) / (a ^ 2 + x ^ 2)
    simplify
    apply integral identity
    substitute u for a^2 + x^2
    simplify
    apply integral identity
    replace substitution
    simplify
    rewrite x ^ 3 / 3 * arctan(x / a) - a / 3 * (x ^ 2 / 2 - a ^ 2 * log(a ^ 2 + x ^ 2) / 2) to x ^ 3 / 3 * arctan(x / a) + a ^ 3 * log(a ^ 2 + x ^ 2) / 6 - a * x ^ 2 / 6
    simplify
done

// 4
prove (INT x. x ^ 3 * arctan(x / a)) = 1 / 4 * (x ^ 4 - a ^ 4) * arctan(x / a) - 1 / 12 * a * x ^ 3 + 1 / 4 * a ^ 3 * x + SKOLEM_CONST(C) for a != 0, x != 0
lhs:
    integrate by parts with u = arctan(x / a), v = x^4 / 4
    simplify
    rewrite x ^ 4 / (x ^ 2 / a ^ 2 + 1) to (a ^ 2 * x ^ 4) / (x ^ 2 + a ^ 2)
    simplify
    rewrite x ^ 4 / (a ^ 2 + x ^ 2) to x ^ 2 - a ^ 2 + a ^ 4 / (a ^ 2 + x ^ 2)
    simplify
    apply integral identity
    simplify
    rewrite x ^ 4 / 4 * arctan(x / a) - a / 4 * (a ^ 3 * arctan(x / a) - a ^ 2 * x + x ^ 3 / 3) to (x^4 / 4 - a^4 / 4) * arctan(x / a) + a^3 * x / 4 - a * x^3 / 12
    simplify
done

// 5
prove (INT x. x ^ n * arctan(x / a)) = x ^ (n + 1) / (n + 1) * arctan(x / a) - a / (n + 1) * (INT x. x ^ (n + 1) / (x ^ 2 + a ^ 2)) for a != 0, n >= 0, isInt(n)
lhs:
    integrate by parts with u = arctan(x / a), v = x ^ (n + 1) / (n + 1)
    simplify
    rewrite x ^ 2 / a ^ 2 + 1 to (x ^ 2 + a ^ 2) / a ^ 2
    rewrite x ^ (n + 1) / ((x ^ 2 + a ^ 2) / a ^ 2) to x ^ (n + 1) * a ^ 2 / (x ^ 2 + a ^ 2)
    simplify
done

## 10.1.6 Integrands Involving x ^ (-n) * arctan(x/a)
### 10.1.6.1

// 1a
prove (INT x. 1 / x * arctan(x / a)) = SUM(k, 0, oo, (-1) ^ k / (2 * k + 1) ^ 2 * (x / a) ^ (2 * k + 1)) + SKOLEM_CONST(C) for abs(x / a) < 1, a != 0
sorry

// 1b
prove (INT x. 1 / x * arctan(x / a)) = pi / 2 * log(abs(x)) + SUM(k, 0, oo, (-1) ^ k / (2 * k + 1) ^ 2 * (x / a) ^ (2 * k + 1)) + SKOLEM_CONST(C) for x / a > 1, a != 0
sorry

// 1c
prove (INT x. 1 / x * arctan(x / a)) = -pi / 2 * log(abs(x)) + SUM(k, 0, oo, (-1) ^ k / (2 * k + 1) ^ 2 * (x / a) ^ (2 * k + 1)) + SKOLEM_CONST(C) for x / a < -1, a != 0
sorry

// 2
prove (INT x. 1 / x ^ 2 * arctan(x / a)) = -1 / x * arctan(x / a) + 1 / (2 * a) * log(x ^ 2 / (x ^ 2 + a ^ 2)) + SKOLEM_CONST(C) for a: real, a != 0, x != 0
lhs:
    integrate by parts with u = arctan(x / a), v = -1 / x
    simplify
    rewrite x^2 / a^2 + 1 to (x^2 + a^2) / a^2
    rewrite 1 / (x * ((x^2 + a^2) / a^2)) to a^2 / (x * (x^2 + a^2))
    rewrite a^2 / (x * (x^2 + a^2)) to a^2 * (1/x - x/(x^2 + a^2)) / a^2
    simplify
    apply integral identity
    rewrite x/(a^2 + x^2) to x/(x^2 + a^2)
    substitute u for x^2 + a^2
    simplify
    apply integral identity
    replace substitution
    simplify
    rewrite log(abs(x)) - log(a ^ 2 + x ^ 2) / 2 to log(abs(x)) - log(abs(x^2 + a^2)^(1/2))
    rewrite log(abs(x)) - log(abs(x^2 + a^2)^(1/2)) to log(abs(x)/abs(x^2 + a^2)^(1/2))
    rewrite abs(x^2 + a^2)^(1/2) to sqrt(x^2 + a^2)
    rewrite abs(x)/sqrt(x^2 + a^2) to sqrt(x^2)/sqrt(x^2 + a^2)
    rewrite sqrt(x^2)/sqrt(x^2 + a^2) to sqrt(x^2/(x^2 + a^2))
    rewrite sqrt(x^2/(x^2 + a^2)) to (x^2/(x^2 + a^2))^(1/2)
    rewrite log((x^2/(x^2 + a^2))^(1/2)) to (1/2)*log(x^2/(x^2 + a^2))
    simplify
done

// 3
prove (INT x. 1 / x ^ 3 * arctan(x / a)) = -1 / 2 * (1 / x ^ 2 + 1 / a ^ 2) * arctan(x / a) - 1 / (2 * a * x) + SKOLEM_CONST(C) for a != 0, x != 0
lhs:
    integrate by parts with u = arctan(x / a), v = -1 / (2 * x ^ 2)
    simplify
    rewrite x^2 / a^2 + 1 to (x^2 + a^2) / a^2
    rewrite 1 / (x ^ 2 * ((x^2 + a^2) / a^2)) to a^2 / (x ^ 2 * (x^2 + a^2))
    partial fraction decomposition
    apply integral identity
    simplify
    rewrite 1 / (2 * a) * (-(1 / a * arctan(x / a)) - 1 / x) - 1 / (2 * x ^ 2) * arctan(x / a) to -1/2 * (1 / x ^ 2 + 1 / a ^ 2) * arctan(x / a) - 1 / (2 * a * x)
done

// 4
prove (INT x. 1 / x ^ n * arctan(x / a)) = -1 / ((n - 1) * x ^ (n - 1)) * arctan(x / a) + a / (n -1) * (INT x. 1 / (x ^ (n - 1) * (x ^ 2 + a ^ 2))) for n != 1, a != 0
sorry

## 10.1.7 Integrands Involving x ^ n * arccot(x/a)
### 10.1.7.1

// 1
prove (INT x. arccot(x / a)) = x * arccot(x / a) + 1 / 2 * a * log(x ^ 2 + a ^ 2) + SKOLEM_CONST(C) for a != 0
lhs:
    integrate by parts with u = arccot(x / a), v = x
    simplify
    rewrite x / (x^2 / a^2 + 1) to (x * a^2) / (x^2 + a^2)
    rewrite (x * a^2) / (x^2 + a^2) to (a / 2) * (2 * x * a) / (x^2 + a^2)
    substitute u for x^2 + a^2
    apply integral identity
    replace substitution
    rewrite log(abs(x^2 + a^2)) to log(x^2 + a^2)
done

// 2
prove (INT x. x * arccot(x / a)) = 1 / 2 * (x ^ 2 + a ^ 2) * arccot(x / a) + 1 / 2 * a * x + SKOLEM_CONST(C) for a != 0
lhs:
    integrate by parts with u = arccot(x / a), v = x ^ 2 / 2
    simplify
    rewrite x ^ 2 / a ^ 2 + 1 to (x ^ 2 + a ^ 2) / a ^ 2
    rewrite x ^ 2 / ((x ^ 2 + a ^ 2) / a ^ 2) to a ^ 2 * x ^ 2 / (x ^ 2 + a ^ 2)
    rewrite a ^ 2 * x ^ 2 / (x ^ 2 + a ^ 2) to a ^ 2 * (1 - a ^ 2 / (x ^ 2 + a ^ 2))
    simplify
    apply integral identity
    simplify
    rewrite arctan(x / a) to pi / 2 - arccot(x / a)
    rewrite a / 2 * (x - a * (pi / 2 - arccot(x / a))) to a / 2 * x - a ^ 2 / 2 * (pi / 2 - arccot(x / a))
    rewrite a ^ 2 / 2 * (pi / 2 - arccot(x / a)) to a ^ 2 * pi / 4 - a ^ 2 / 2 * arccot(x / a)
    simplify
    rewrite a ^ 2 / 2 * arccot(x / a) + x ^ 2 / 2 * arccot(x / a) to (a ^ 2 / 2 + x ^ 2 / 2) * arccot(x / a)
    simplify
    rewrite (a ^ 2 / 2 + x ^ 2 / 2) to 1 / 2 * (a ^ 2 + x ^ 2)
    rewrite a * x / 2 to 1 / 2 * a * x
    rewrite arccot(x / a) * (1 / 2 * (a ^ 2 + x ^ 2)) to 1 / 2 * (a ^ 2 + x ^ 2) * arccot(x / a)
    rewrite 1 / 2 * (a ^ 2 + x ^ 2) * arccot(x / a) to 1 / 2 * (x ^ 2 + a ^ 2) * arccot(x / a)
    rewrite to 1 / 2 * (x ^ 2 + a ^ 2) * arccot(x / a) + 1 / 2 * a * x + SKOLEM_CONST(C)
done

// 3
prove (INT x. x ^ 2 * arccot(x / a)) = 1 / 3 * x ^ 3 * arccot(x / a) + 1 / 6 * a * x ^ 2 - 1 / 6 * a ^ 3 * log(x ^ 2 + a ^ 2) + SKOLEM_CONST(C) for a != 0
lhs:
    integrate by parts with u = arccot(x / a), v = x^3 / 3
    simplify
    rewrite x ^ 3 / (x ^ 2 / a ^ 2 + 1) to (a ^ 2 * x ^ 3) / (x ^ 2 + a ^ 2)
    simplify
    rewrite x ^ 3 / (a ^ 2 + x ^ 2) to x - (a ^ 2 * x) / (a ^ 2 + x ^ 2)
    simplify
    apply integral identity
    substitute u for a^2 + x^2
    simplify
    apply integral identity
    replace substitution
    simplify
    rewrite a / 3 * (x ^ 2 / 2 - a ^ 2 * log(a ^ 2 + x ^ 2) / 2) + x ^ 3 / 3 * arccot(x / a) to x ^ 3 / 3 * arccot(x / a) + a * x ^ 2 / 6 - a ^ 3 * log(a ^ 2 + x ^ 2) / 6
done

// 4
prove (INT x. x ^ 3 * arccot(x / a)) = 1 / 4 * (x ^ 4 - a ^ 4) * arccot(x / a) + 1 / 12 * a * x ^ 3 - 1 / 4 * a ^ 3 * x + SKOLEM_CONST(C) for a != 0
lhs:
    integrate by parts with u = arccot(x / a), v = x^4 / 4
    simplify
    rewrite x ^ 4 / (x ^ 2 / a ^ 2 + 1) to (a ^ 2 * x ^ 4) / (x ^ 2 + a ^ 2)
    simplify
    rewrite x ^ 4 / (a ^ 2 + x ^ 2) to x ^ 2 - a ^ 2 + a ^ 4 / (a ^ 2 + x ^ 2)
    simplify
    apply integral identity
    simplify
    rewrite arctan(x / a) to pi / 2 - arccot(x / a)
    rewrite a / 4 * (a ^ 3 * (pi / 2 - arccot(x / a)) - a ^ 2 * x + x ^ 3 / 3) + x ^ 4 / 4 * arccot(x / a) to (x^4 / 4 - a^4 / 4) * arccot(x / a) + a * x^3 / 12 - a^3 * x / 4 + (a^4 * pi / 8)
    rewrite (x^4 / 4 - a^4 / 4) * arccot(x / a) + a * x^3 / 12 - a^3 * x / 4 + (a^4 * pi / 8) + SKOLEM_CONST(C) to (x^4 / 4 - a^4 / 4) * arccot(x / a) + a * x^3 / 12 - a^3 * x / 4 + SKOLEM_CONST(C)
    simplify
done

// 5
prove (INT x. x ^ n * arccot(x / a)) = x ^ (n + 1) / (n + 1) * arccot(x / a) + a / (n + 1) * (INT x. x ^ (n + 1) / (x ^ 2 + a ^ 2)) for a != 0, n >= 0, isInt(n)
lhs:
    integrate by parts with u = arccot(x / a), v = x ^ (n + 1) / (n + 1)
    simplify
    rewrite x ^ 2 / a ^ 2 + 1 to (x ^ 2 + a ^ 2) / a ^ 2
    rewrite x ^ (n + 1) / ((x ^ 2 + a ^ 2) / a ^ 2) to x ^ (n + 1) * a ^ 2 / (x ^ 2 + a ^ 2)
    simplify
done

## 10.1.8 Integrands Involving x ^ (-n) * arccot(x/a)
### 10.1.8.1

// 1a
prove (INT x. 1 / x * arccot(x / a)) = pi / 2 * log(abs(x)) - SUM(k, 0, oo, (-1) ^ k/(2 * k + 1) ^ 2 * (x / a) ^ (2 * k + 1)) + SKOLEM_CONST(C) for abs(x/a) < 1, a != 0
sorry

// 1b
prove (INT x. 1 / x * arccot(x / a)) = -SUM(k, 0, oo, (-1) ^ k / (2 * k + 1) ^ 2 * (x / a) ^ (2 * k + 1)) + SKOLEM_CONST(C) for x / a > 1, a != 0
sorry

// 1c
prove (INT x. 1 / x * arccot(x / a)) = pi * log(abs(x)) - SUM(k, 0, oo, (-1) ^ k / (2 * k + 1) ^ 2 * (x / a) ^ (2 * k + 1)) + SKOLEM_CONST(C) for x / a < -1, a != 0
sorry

// 2
prove (INT x. 1 / x ^ 2 * arccot(x / a)) = -1 / x * arccot(x / a) - 1 / (2 * a) * log(x ^2 / (x ^ 2 + a ^ 2)) + SKOLEM_CONST(C) for x != 0, a != 0
lhs:
    integrate by parts with u = arccot(x / a), v = -1 / x
    simplify
    rewrite x^2 / a^2 + 1 to (x^2 + a^2) / a^2
    rewrite 1 / (x * ((x^2 + a^2) / a^2)) to a^2 / (x * (x^2 + a^2))
    rewrite a^2 / (x * (x^2 + a^2)) to a^2 * (1/x - x/(x^2 + a^2)) / a^2
    simplify
    apply integral identity
    rewrite x/(a^2 + x^2) to x/(x^2 + a^2)
    substitute u for x^2 + a^2
    simplify
    apply integral identity
    replace substitution
    simplify
    rewrite log(abs(x)) - log(a ^ 2 + x ^ 2) / 2 to log(abs(x)) - log(abs(x^2 + a^2)^(1/2))
    rewrite log(abs(x)) - log(abs(x^2 + a^2)^(1/2)) to log(abs(x)/abs(x^2 + a^2)^(1/2))
    rewrite abs(x^2 + a^2)^(1/2) to sqrt(x^2 + a^2)
    rewrite abs(x)/sqrt(x^2 + a^2) to sqrt(x^2)/sqrt(x^2 + a^2)
    rewrite sqrt(x^2)/sqrt(x^2 + a^2) to sqrt(x^2/(x^2 + a^2))
    rewrite sqrt(x^2/(x^2 + a^2)) to (x^2/(x^2 + a^2))^(1/2)
    rewrite log((x^2/(x^2 + a^2))^(1/2)) to (1/2)*log(x^2/(x^2 + a^2))
    simplify
done

// 3
prove (INT x. 1 / x ^ 3 * arccot(x / a)) = -1 / 2 * (1 / x ^ 2 + 1 / a ^ 2) * arccot(x / a) + 1 / (2 * a * x) + SKOLEM_CONST(C) for x != 0, a != 0
lhs:
    integrate by parts with u = arccot(x / a), v = -1 / (2 * x ^ 2)
    simplify
    rewrite x^2 / a^2 + 1 to (x^2 + a^2) / a^2
    rewrite 1 / (x ^ 2 * ((x^2 + a^2) / a^2)) to a^2 / (x ^ 2 * (x^2 + a^2))
    partial fraction decomposition
    apply integral identity
    simplify
    rewrite arctan(x / a) to pi / 2 - arccot(x / a)
    rewrite to -1/2 * (1 / x ^ 2 + 1 / a ^ 2) * arccot(x / a) + 1 / (2 * a * x) + SKOLEM_CONST(C)
done

// 4
prove (INT x. 1 / x ^ n * arccot(x / a)) = -1 / ((n - 1) * x ^ (n - 1)) * arccot(x / a) - a / (n -1) * (INT x. 1 / (x ^ (n - 1) * (x ^ 2 + a ^ 2))) for n != 1, a != 0, x != 0
sorry

## 10.1.9 Integrands Involving Products of Rational Functions and arccot(x/a)
## 10.1.9.1

// 1
prove (INT x. 1 / (x ^ 2 + a ^ 2) * arccot(x / a)) = -1 / (2 * a) * arccot(x / a) ^ 2 + SKOLEM_CONST(C) for a != 0, x != 0
lhs:
    rewrite x ^ 2 + a ^ 2 to a ^ 2 * (x ^ 2 / a ^ 2 + 1)
    rewrite 1 / (a ^ 2 * (x ^ 2 / a ^ 2 + 1)) * arccot(x / a) to (1/a) * arccot(x / a) * (1/(a * (x ^ 2 / a ^ 2 + 1)))
    rewrite (1/a) * arccot(x / a) * (1/(a * (x ^ 2 / a ^ 2 + 1))) to (-1/a) * arccot(x / a) * (-(1/(a * (x ^ 2 / a ^ 2 + 1))))
    substitute u for arccot(x / a)
    apply integral identity
    replace substitution
    simplify
done

// 2
prove (INT x. x ^ 2 / (x ^ 2 + a ^ 2) * arccot(x / a)) = x * arccot(x / a) +1 / 2 * a * log(x ^ 2 + a ^ 2) + 1 / 2 * a * arccot(x / a) ^ 2 + SKOLEM_CONST(C) for a != 0, x != 0
sorry

// 3
prove (INT x. 1 / (x ^ 2 + a ^ 2) ^ 2 * arccot(x / a)) = x / (2 * a ^ 2 * (x ^ 2 + a ^ 2)) * arccot(x / a) - 1 / (4 * a ^ 3) * arccot(x / a) ^ 2 -1 / (4 * a * (x ^ 2 + a ^ 2)) + SKOLEM_CONST(C) for x != 0, a != 0
sorry

// 4
prove (INT x. 1 / (x ^ 2 + a ^ 2) * arccot(x / a) ^ n) = -1 / ((n + 1) * a) * arccot(x / a) ^ (n + 1) + SKOLEM_CONST(C) for x != 0, a != 0, n >= 0, isInt(n)
lhs:
    rewrite x ^ 2 + a ^ 2 to a ^ 2 * (x ^ 2 / a ^ 2 + 1)
    rewrite 1 / (a ^ 2 * (x ^ 2 / a ^ 2 + 1)) * arccot(x / a) ^ n to (1/a) * arccot(x / a) ^ n * (1/(a * (x ^ 2 / a ^ 2 + 1)))
    rewrite (1/a) * arccot(x / a) ^ n * (1/(a * (x ^ 2 / a ^ 2 + 1))) to (-1/a) * arccot(x / a) ^ n * (-(1/(a * (x ^ 2 / a ^ 2 + 1))))
    substitute u for arccot(x / a)
    apply integral identity
    replace substitution
    simplify
done