imports interesting5

# Inside Interesting Integrals, Chapter 6

## Chapter 6, Section 1, Bernoulli's Integral

// (6.1.2)

prove (INT x:[0,1]. x ^ (c * x ^ a)) = SUM(k, 0, oo, (-c) ^ k / (k * a + 1) ^ (k + 1)) for a c: real, a > 0, c != 0
subgoal 1: converges(SUM(k, 0, oo, abs(INT x:[0,1]. (c * x ^ a * log(x)) ^ k / factorial(k))))
arg:
    simplify
    rewrite (c * x ^ a * log(x)) ^ k to (c * x ^ a) ^ k * log(x) ^ k
    rewrite (c * x ^ a) ^ k to c ^ k * x ^ a ^ k
    simplify
    apply integral identity
    simplify
done
lhs:
    rewrite x ^ (c * x ^ a) to exp(log(x ^ (c * x ^ a)))
    apply series expansion on exp(log(x ^ (c * x ^ a))) index k
    exchange integral and sum
    rewrite log(x ^ (c * x ^ a)) to c * x ^ a * log(x)
    rewrite (c * x ^ a * log(x)) ^ k to (c * x ^ a) ^ k * log(x) ^ k
    rewrite (c * x ^ a) ^ k to c ^ k * x ^ a ^ k
    simplify
    apply integral identity
    simplify
    rewrite c ^ k * (-1) ^ k to (-c) ^ k
done

// (6.1.3)

prove (INT x:[0,1]. x ^ x) = SUM(k, 0, oo, (-1) ^ k * (k + 1) ^ (-k - 1))
lhs:
    rewrite x ^ x to x ^ (1 * x ^ 1)
    apply integral identity
    simplify
done

// (6.1.4)

prove (INT x:[0,1]. x ^ -x) = SUM(k, 0, oo, (k + 1) ^ (-k - 1))
lhs:
    rewrite x ^ -x to x ^ (-1 * x ^ 1)
    apply integral identity
    simplify
done

// (6.1.5)

prove (INT x:[0,1]. x ^ (x ^ 2)) = SUM(k, 0, oo, (-1) ^ k * (2 * k + 1) ^ (-k - 1))
lhs:
    rewrite x ^ (x ^ 2) to x ^ (1 * x ^ 2)
    apply integral identity
    simplify
done

// (6.1.6)

prove (INT x:[0,1]. x ^ sqrt(x)) = SUM(k, 0, oo, (-1) ^ k * (2 / (k + 2)) ^ (k + 1))
lhs:
    rewrite x ^ sqrt(x) to x ^ (1 * x ^ (1/2))
    apply integral identity
    simplify
    rewrite k / 2 + 1 to (2 / (k + 2)) ^ (-1)
    rewrite (2 / (k + 2)) ^ (-1) ^ (-k - 1) to (2 / (k + 2)) ^ (k + 1)
done

## Chapter 6, Section 2, Ahmed's Integral

prove (INT x:[0,1]. arctan(sqrt(2 + x ^ 2)) / ((1 + x ^ 2) * sqrt(2 + x ^ 2))) = 5 * pi ^ 2 / 96
define I(u) = (INT x:[0,1]. arctan(u * sqrt(2 + x ^ 2)) / ((1 + x ^ 2) * sqrt(2 + x ^ 2))) for u: real, u > 0
subgoal 1: I(1) = (INT x:[0,1]. arctan(sqrt(x ^ 2 + 2)) / ((x ^ 2 + 1) * sqrt(x ^ 2 + 2)))
lhs:
    expand definition for I
done
subgoal 2: (D u. I(u)) = 1 / (1 + u ^ 2) * (pi / 4 - u / sqrt(1 + 2 * u ^ 2) * arctan(u / sqrt(1 + 2 * u ^ 2))) for u: real, u > 0
lhs:
    expand definition for I (all)
    exchange derivative and integral
    simplify
    rewrite 1 / ((x ^ 2 + 1) * (u ^ 2 * (x ^ 2 + 2) + 1)) to 1 / (u ^ 2 + 1) * (1 / (1 + x ^ 2) - u ^ 2 / (1 + 2 * u ^ 2 + u ^ 2 * x ^ 2))
    simplify
    rewrite 1 / (u ^ 2 * x ^ 2 + 2 * u ^ 2 + 1) to u ^ (-2) * (x ^ 2 + (2 * u ^ 2 + 1) / u ^ 2) ^ (-1)
    simplify
    substitute y * sqrt(u ^ (-2) * (2 * u ^ 2 + 1)) for x
    simplify
    rewrite 1 / (y ^ 2 * (2 * u ^ 2 + 1) / u ^ 2 + (2 * u ^ 2 + 1) / u ^ 2) to 1 / (y ^ 2 + 1) * (u ^ 2 / (2 * u ^ 2 + 1))
    apply integral identity
    simplify
done
subgoal 3: (INT u:[1,oo]. D u. I(u)) = pi ^ 2 / 12 - I(1)
lhs:
    simplify
    expand definition for I (at 1)
    simplify
    integrate by parts with u = 1, v = arctan(x / sqrt(2 + x ^ 2)) / 2
    simplify
done
subgoal 4: (INT u:[1,oo]. D u. I(u)) = -(pi ^ 2 / 48) + I(1)
lhs:
    apply 2 on D u. I(u)
    expand polynomial
    simplify
    substitute 1 / x for u
    simplify
    rewrite x ^ 3 * (1 / x ^ 2 + 1) * sqrt(2 / x ^ 2 + 1) to sqrt((1 + x ^ 2) ^ 2 * (2 + x ^ 2))
    rewrite x * sqrt(2 / x ^ 2 + 1) to sqrt(x ^ 2 + 2)
    simplify
    rewrite 1 / sqrt(x ^ 2 + 2) to sqrt(x ^ 2 + 2) ^ (-1)
    rewrite arctan(sqrt(x ^ 2 + 2) ^ (-1)) to pi / 2 - arctan(sqrt(x ^ 2 + 2))
    expand polynomial
    simplify
    rewrite arctan(sqrt(x ^ 2 + 2)) / (x ^ 2 * sqrt(x ^ 2 + 2) + sqrt(x ^ 2 + 2)) to arctan(sqrt(x ^ 2 + 2)) / ((x ^ 2 + 1) * sqrt(x ^ 2 + 2))
    apply 1 on INT x:[0,1]. arctan(sqrt(x ^ 2 + 2)) / ((x ^ 2 + 1) * sqrt(x ^ 2 + 2))
    integrate by parts with u = 1, v = arctan(x / sqrt(2 + x ^ 2))
    apply integral identity
    simplify
done
from 3:
    apply 4 on INT u:[1,oo]. D u. I(u)
    solve equation for I(1)
    expand definition for I (all)
done
