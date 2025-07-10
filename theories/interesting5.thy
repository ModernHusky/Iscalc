imports interesting4

# Inside Interesting Integrals, Chapter 5

## Chapter 5, Section 1, Catalan's Constant

define G = SUM(n, 0, oo, (-1)^n / (2*n+1)^2)

// (5.1.1)

prove (INT x:[0,1]. arctan(x) / x) = G
subgoal 1: converges(SUM(n, 0, oo, INT x:[0,1]. x ^ (2 * n) / (2 * n + 1)))
arg:
    simplify
    apply integral identity
    simplify
done
lhs:
    apply series expansion on arctan(x) index n
    rewrite x ^ (2 * n + 1) to x ^ (2 * n) * x
    simplify
    exchange integral and sum
    apply integral identity
    simplify
rhs:
    expand definition for G
done

prove (INT x:[1,oo]. log(x) / x ^ k) = 1 / (k - 1) ^ 2 for k > 1
lhs:
    improper integral to limit creating t
    integrate by parts with u = log(x), v = x ^ (1 - k) / (1 - k)
    simplify
    apply integral identity
    simplify
rhs:
    rewrite 1 / (k - 1) ^ 2 to 1 / (-k + 1) ^ 2
done

// (5.1.2)

prove (INT x:[1,oo]. log(x) / (x ^ 2 + 1)) = G
subgoal 1: converges(SUM(n, 0, oo, INT x:[1,oo]. x ^ (-(2 * n) - 2) * log(x)))
arg:
    rewrite x ^ (-(2 * n) - 2) * log(x) to log(x) / x ^ (2 * n + 2)
    apply integral identity
    simplify
done
lhs:
    rewrite log(x) / (x ^ 2 + 1) to log(x) * x ^ (-2) * (1 + 1 / x ^ 2) ^ (-1)
    apply series expansion on (1 + 1 / x ^ 2) ^ (-1) index n
    rewrite log(x) * x ^ (-2) * SUM(n, 0, oo, (-1) ^ n * (1 / x ^ 2) ^ n) to SUM(n, 0, oo, (-1) ^ n * (1 / x ^ 2) ^ n * log(x) * x ^ (-2))
    rewrite (1 / x ^ 2) ^ n to x ^ (-(2 * n))
    rewrite (-1) ^ n * x ^ (-(2 * n)) * log(x) * x ^ (-2) to (-1) ^ n * x^(-(2*n)-2) * log(x)
    exchange integral and sum
    simplify
    rewrite x ^ (-(2 * n) - 2) * log(x) to log(x) / x ^ (2 * n + 2)
    apply integral identity
    simplify
rhs:
    expand definition for G
done

// (5.1.3)

prove (INT x:[0,oo]. log(x + 1) / (x ^ 2 + 1)) = pi / 4 * log(2) + G
lhs:
    split region at 1
    apply integral identity
    rewrite x + 1 to x * (1 + 1 / x)
    rewrite log(x * (1 + 1 / x)) to log(x) + log(1 + 1 / x)
    rewrite (log(x) + log(1 + 1 / x)) / (x ^ 2 + 1) to log(x) / (x ^ 2 + 1) + log(1 + 1 / x) / (x ^ 2 + 1)
    simplify
    apply integral identity
    substitute u for 1 / x
    rewrite u ^ 2 * (1 / u ^ 2 + 1) to u ^ 2 + 1
    apply integral identity
    simplify
    rewrite pi * log(2) / 4 to pi / 4 * log(2)
done

// (5.1.4)

prove (INT x:[0,pi]. x * sin(x) / (a + b * cos(x) ^ 2)) = pi / sqrt(a * b) * arctan(sqrt(b / a)) for a b: real, a > 0, b > 0
let I(a,b) = (INT x:[0,pi]. x * sin(x) / (a + b * cos(x) ^ 2))
subgoal 1: I(a,b) = (INT x:[0,pi]. (pi - x) * sin(x) / (a + b * cos(x) ^ 2))
lhs:
    expand definition for I
    substitute x for pi - x
    rewrite sin(x) * (-x + pi) / (b * cos(x) ^ 2 + a) to (pi - x) * sin(x) / (a + b * cos(x) ^ 2)
done
lhs:
    fold definition for I (all)
    rewrite I(a,b) to 1/2 * (I(a,b) + I(a,b))
    expand definition for I (at 1)
    apply 1 on I(a,b)
    rewrite (INT x:[0,pi]. x * sin(x) / (b * cos(x) ^ 2 + a)) + (INT x:[0,pi]. (pi - x) * sin(x) / (a + b * cos(x) ^ 2)) to INT x:[0,pi]. x * sin(x) / (a + b * cos(x) ^ 2) + (pi - x) * sin(x) / (a + b * cos(x) ^ 2)
    rewrite x * sin(x) / (a + b * cos(x) ^ 2) + (pi - x) * sin(x) / (a + b * cos(x) ^ 2) to pi * sin(x) / (a + b * cos(x) ^ 2)
    substitute u for cos(x)
    substitute x for sqrt(b / a) * u
    apply integral identity
    rewrite 1 / (-(a * x ^ 2) - a) to -1/a * (1/(x^2+1))
    apply integral identity
    simplify
    rewrite arctan(-(sqrt(b) / sqrt(a))) to -arctan(sqrt(b) / sqrt(a))
    simplify
done

## Chapter 5, Section 2, Power Series for the Log Functions

// (5.2.1)

prove (INT x:[0,1]. log(1 + x) / x) = pi ^ 2 / 12
subgoal 1: converges(SUM(n, 0, oo, INT x:[0,1]. x ^ n / (n + 1)))
arg:
    simplify
    apply integral identity
    simplify
done
lhs:
    apply series expansion on log(1 + x) index n
    rewrite SUM(n, 0, oo, (-1) ^ n * x ^ (n + 1) / (n + 1)) / x to SUM(n, 0, oo, (-1) ^ n * x ^ (n + 1) / (n + 1) * (1 / x))
    exchange integral and sum
    simplify
    apply integral identity
    simplify
    apply series evaluation
done

// (5.2.2)

prove (INT x:[0, 1]. log(1 - x) / x) = -(pi ^ 2 / 6)
subgoal 1:(INT x:[0,1]. (x ^ n / (n + 1))) = (1/(n+1)^2) for n>=0
lhs:
    apply integral identity
    simplify
done
subgoal 2:SUM(n, 0, oo, INT x:[0, 1]. x^n / (n+1)) = pi ^ 2 / 6
lhs:
    apply 1 on INT x:[0, 1]. x^n / (n+1)
    apply series evaluation
done
subgoal 3: SUM(n, 0, oo, 1 / (n+1)^2) = pi ^ 2 / 6
lhs:
    apply series evaluation
done
lhs:
    apply series expansion on log(1-x) index n
    rewrite SUM(n, 0, oo, (-1) ^ n * (-x) ^ (n + 1) / (n + 1)) / x to SUM(n, 0, oo, (-1) ^ n * (-x) ^ (n + 1) / (n + 1))*(1/x)
    rewrite (-x) ^ (n + 1) to (-1) ^ (n + 1) * x ^ (n + 1)
    rewrite SUM(n, 0, oo, (-1) ^ n * ((-1) ^ (n + 1) * x ^ (n + 1)) / (n + 1)) * (1 / x) to SUM(n, 0, oo, (-1) ^ n * ((-1) ^ (n + 1) * x ^ (n + 1)) / (n + 1) * (1 / x))
    rewrite (-1) ^ n * ((-1) ^ (n + 1) * x ^ (n + 1)) / (n + 1) * (1 / x) to (-1) ^ n * ((-1) ^ (n + 1) * x ^ n) / (n + 1)
    rewrite SUM(n, 0, oo, (-1) ^ n * ((-1) ^ (n + 1) * x ^ n) / (n + 1)) to SUM(n, 0, oo, (-1) ^ (2*n+1) * (x ^ n) / (n + 1))
    simplify
    exchange integral and sum
    apply 1 on (INT x:[0,1]. x ^ n / (n + 1))
    apply 3 on SUM(n, 0, oo, 1 / (n + 1) ^ 2)
done

// (5.2.4)

prove (INT x:[0, pi/2]. cos(x)/sin(x) * log(1/cos(x))) = pi^2/24
subgoal 1: (-log(1-x) - log(1+x)) = -SUM(k,0,oo,(-1)^k*(-x)^(k+1) / (k+1))-SUM(k,0,oo,(-1)^k*x^(k+1)/(k+1)) for x != 0, abs(x) < 1
lhs:
    apply series expansion on log(1-x) index k
    apply series expansion on log(1+x) index k
done
subgoal 2:x / (-(x ^ 2) + 1) = 1/2 * SUM(k, 0, oo, x ^ k) - 1/2 * SUM(k, 0, oo, x ^ k * (-1) ^ k) for x != 0, abs(x) < 1
from 1:
    differentiate both sides at x
    simplify
    rewrite 1 / (-x + 1) - 1 / (x + 1) to 2 * (x / (1-x^2))
    solve equation for x / (1-x^2)
    rewrite (-1) ^ k * (-x) ^ k to x ^ k
    rewrite (SUM(k, 0, oo, x ^ k) - SUM(k, 0, oo, x ^ k * (-1) ^ k)) / 2 to 1/2 * SUM(k, 0, oo, x ^ k) - 1/2 * SUM(k, 0, oo, x ^ k * (-1) ^ k)
done
subgoal 7: converges(SUM(k, 0, oo, INT y:[0,1]. log(y) * y ^ k))
arg:
    integrate by parts with u=log(y),v=y^(k+1)/(k+1)
    simplify
    apply integral identity
    simplify
done
subgoal 3:(INT y:[0,1]. (SUM(k, 0, oo, log(y) * y ^ k * (-1) ^ k))) = -SUM(k, 0, oo, (-1) ^ k / (k + 1) ^ 2)
lhs:
    exchange integral and sum
    apply integral identity
    integrate by parts with u=log(y),v=y^(k+1)/(k+1)
    simplify
    apply integral identity
    simplify
done
subgoal 4:(INT y:[0,1]. SUM(k, 0, oo, log(y) * y ^ k)) = -SUM(k, 0, oo, 1 / (k + 1) ^ 2)
lhs:
    exchange integral and sum
    integrate by parts with u=log(y),v=y^(k+1)/(k+1)
    simplify
    apply integral identity
    simplify
done
subgoal 5: SUM(k, 0, oo, (-1) ^ k / (k + 1) ^ 2) = pi^2/12
lhs:
    apply series evaluation
done
subgoal 6: SUM(k, 0, oo, 1 / (k + 1) ^ 2) = pi^2/6
lhs:
    apply series evaluation
done
lhs:
    substitute t for cos(x)
    simplify
    substitute y for t
    rewrite y * log(y) / (-(y ^ 2) + 1) to log(y) * (y / (-(y ^ 2) + 1))
    apply 2 on y / (-(y ^ 2) + 1)
    rewrite log(y) * (1/2 * SUM(k, 0, oo, y ^ k) - 1/2 * SUM(k, 0, oo, y ^ k * (-1) ^ k)) to 1/2 * log(y) * SUM(k, 0, oo, y ^ k) - 1/2 * log(y) * SUM(k, 0, oo, y ^ k * (-1) ^ k)
    expand polynomial 
    simplify
    rewrite (log(y) * SUM(k, 0, oo, y ^ k)) to SUM(k, 0, oo, log(y) * y ^ k)
    rewrite log(y) * SUM(k, 0, oo, y ^ k * (-1) ^ k) to SUM(k, 0, oo, log(y) * y ^ k * (-1) ^ k)
    apply 3 on (INT y:[0,1]. (SUM(k, 0, oo, log(y) * y ^ k * (-1) ^ k)))
    apply 4 on (INT y:[0,1]. SUM(k, 0, oo, log(y) * y ^ k))
    apply 5 on SUM(k, 0, oo, (-1) ^ k / (k + 1) ^ 2)
    apply 6 on SUM(k, 0, oo, 1 / (k + 1) ^ 2)
    simplify
done

## Chapter 5, Section 4, Euler's Constant and Related Integrals

// Use (5.4.3) as definition

define EulerConstant = -(INT x:[0, oo]. exp(-x) * log(x))

// 5.4.1

prove (INT x:[0,1]. (-exp(-x) + 1) / x) - (INT x:[1,oo]. exp(-x) / x) = EulerConstant
lhs:
    integrate by parts with u = exp(-x), v = log(x) (at 2)
    integrate by parts with u = 1 - exp(-x), v = log(x)
    simplify
rhs:
    expand definition for EulerConstant
    split region at 1
    simplify
done
