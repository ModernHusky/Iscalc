
// Chapter 1, Section 5

prove (INT x:[0,oo]. log(x) / (x ^ 2 + 1)) = 0
lhs:
    split region at 1
    substitute 1 / u for x
    simplify
    rewrite u ^ 2 * (1 / u ^ 2 + 1) to u ^ 2 + 1
    simplify
done

// Chapter 1, Section 7

prove (INT x:[0,1]. (x^4*(1-x)^4)/(1+x^2)) = 22/7 - pi
lhs:
    rewrite (x^4*(1-x)^4)/(1+x^2) to (x^6-4*x^5+5*x^4-4*x^2+4)-4/(1+x^2)
    simplify
    apply integral identity
    simplify
done

// Chapter 2, Section 1

prove (INT x:[1,oo]. 1 / ((x+a)*sqrt(x-1))) = pi / sqrt(a+1) for a: real, a > -1
lhs:
    substitute t for sqrt(x - 1)
    simplify
    substitute y for t / sqrt(a + 1)
    rewrite y ^ 2 * (a + 1) + a + 1 to (a + 1) * (y^2 + 1)
    apply integral identity
    simplify
done

prove (INT x:[0, oo]. log(1 + a^2 / x^2)) = a * pi for a: real, a > 0
lhs:
    integrate by parts with u = log(1+a^2/x^2), v = x
    simplify
    rewrite x^2 * (a^2 / x^2 + 1) to a^2 + x^2
    apply integral identity
    simplify
done

prove (INT x:[0, oo]. log(x) / (x^2+b^2)) = pi * log(b) / (2*b) for b: real, b > 0
lhs:
    substitute 1/t for x
    rewrite log(1/t) to -log(t)
    rewrite -log(t) / ((1/t)^2 + b^2) * -(1/t^2) to log(t) / (1 + b^2*t^2)
    substitute s/b for t
    rewrite log(s/b) to log(s) - log(b)
    expand polynomial
    apply integral identity
    simplify
done

prove (INT x:[1,oo]. log(x) / (x+1)^2) = log(2)
subgoal 1: (INT x:[0,oo]. 1 / (1 + exp(a*x))) = log(2) / a for a: real, a > 0
lhs:
    substitute u for exp(a * x)
    simplify
    rewrite 1 / (u * (u+1)) to 1/u - 1/(u+1)
    simplify
    substitute y for u + 1
    apply integral identity
    simplify
done
subgoal 2: (INT x:[1,oo]. log(x) / (a ^ 2 * (x + 1) ^ 2)) = log(2) / a^2 for a: real, a > 0
from 1:
    differentiate both sides at a
    simplify
    substitute y for exp(a * x)
    solve equation for INT y:[1,oo]. log(y) / (a ^ 2 * (y + 1) ^ 2)
done
lhs:
    rewrite (x+1) ^ 2 to 1^2 * (x+1)^2
    apply 2 on INT x:[1,oo]. log(x) / (1 ^ 2 * (x + 1) ^ 2)
    simplify
done

prove (INT x:[sqrt(2),oo]. 1 / (x + x ^ sqrt(2))) = (1 + sqrt(2)) * log(1 + 2 ^ (1/2 * (1 - sqrt(2))))
lhs:
    rewrite 1 / (x + x ^ sqrt(2)) to x ^ -sqrt(2) / (x ^ (-sqrt(2) + 1) + 1)
    substitute u for log(x ^ (1 - sqrt(2)) + 1)
    apply integral identity
    simplify
    rewrite -sqrt(2) + 1 to -1 / (1 + sqrt(2)) (at 2)
    simplify
    rewrite sqrt(2) to 2 ^ (1/2) (at 2)
    rewrite 2 ^ (1/2) ^ (-sqrt(2) + 1) to 2 ^ (1/2 * (-sqrt(2) + 1))
    rewrite (sqrt(2) + 1) * log(2 ^ (1/2 * (-sqrt(2) + 1)) + 1) to (1 + sqrt(2)) * log(1 + 2 ^ (1/2 * (1 - sqrt(2))))
done

prove (INT x:[-oo,oo]. 1 / cosh(x)) = pi
lhs:
    expand definition for cosh (all)
    substitute t for exp(x)
    rewrite t * (1 / t + t) to 1 + t ^ 2
    apply integral identity
    simplify
done

// Chapter 2, Section 2

calculate INT x:[0,pi / 2]. sqrt(sin(x)) / (sqrt(sin(x)) + sqrt(cos(x)))
    substitute y for pi / 2 - x
    rewrite sqrt(cos(y)) / (sqrt(cos(y)) + sqrt(sin(y))) to 1 - sqrt(sin(y)) / (sqrt(cos(y)) + sqrt(sin(y)))
    apply integral identity
    solve integral INT x:[0,pi / 2]. sqrt(sin(x)) / (sqrt(sin(x)) + sqrt(cos(x)))
done

calculate INT x:[0,pi]. x * sin(x) / (1 + cos(x) ^ 2)
    substitute y for pi - x
    expand polynomial
    simplify
    solve integral INT x:[0,pi]. x * sin(x) / (1 + cos(x) ^ 2)
    substitute u for cos(y)
    rewrite -(1 / (-(u ^ 2) - 1)) to 1/(u^2+1)
    apply integral identity
    simplify
done

prove (INT x:[0,pi / 2]. sin(x) ^ 2 / (sin(x) + cos(x))) = sqrt(2) / 4 * log(3 + 2 * sqrt(2))
subgoal 1: (INT x:[0,pi / 2]. sin(x) ^ 2 / (sin(x) + cos(x))) = (INT x:[0,pi / 2]. cos(x) ^ 2 / (sin(x) + cos(x)))
lhs:
    substitute y for pi / 2 - x
done
subgoal 2: (INT x:[0,pi / 2]. sin(x) ^ 2 / (sin(x) + cos(x))) = 1/2 * (INT x:[0,pi / 2]. 1 / (sin(x) + cos(x)))
rhs:
    simplify
    rewrite 1 to sin(x) ^ 2 + cos(x) ^ 2
    expand polynomial
    simplify
    rewrite cos(x) + sin(x) to sin(x) + cos(x) (at 1)
    apply 1 on INT x:[0,pi / 2]. cos(x) ^ 2 / (sin(x) + cos(x))
    simplify
done
lhs:
    apply 2 on INT x:[0,pi / 2]. sin(x) ^ 2 / (sin(x) + cos(x))
    substitute z for tan(x / 2)
    simplify
    rewrite (-(z ^ 2) + 1) / (z ^ 2 + 1) + 2 * z / (z ^ 2 + 1) to (2 - (z - 1) ^ 2) / (z ^ 2 + 1)
    rewrite (z ^ 2 + 1) * ((2 - (z - 1) ^ 2) / (z ^ 2 + 1)) to 2 - (z - 1) ^ 2
    rewrite 2 - (z - 1) ^ 2 to (sqrt(2) + (z - 1)) * (sqrt(2) - (z - 1))
    rewrite 1 / ((sqrt(2) + (z - 1)) * (sqrt(2) - (z - 1))) to sqrt(2) / 4 * (1 / (sqrt(2) + (z - 1)) + 1 / (sqrt(2) - (z - 1)))
    simplify
    substitute u for sqrt(2) + 1 - z (at 1)
    substitute u for sqrt(2) - 1 + z (at 2)
    apply integral identity
    simplify
    rewrite sqrt(2) * (log(sqrt(2) + 1) - log(sqrt(2) - 1)) / 4 to 1/4 * sqrt(2) * (log(sqrt(2) + 1) - log(sqrt(2) - 1))
    rewrite log(sqrt(2) + 1) - log(sqrt(2) - 1) to log((sqrt(2) + 1) / (sqrt(2) - 1))
    rewrite (sqrt(2) + 1) / (sqrt(2) - 1) to 3 + 2 * sqrt(2)
    simplify
done

prove (INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)) = pi / 8 * log(2)
subgoal 1: (INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)) = (INT x:[0,pi / 4]. log(tan(x) + 1))
lhs:
    substitute tan(u) for x
    rewrite sec(u) ^ 2 to tan(u) ^ 2 + 1
    simplify
done
subgoal 2: (INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)) = pi / 4 * log(2) - (INT x:[0,1]. log(x + 1) / (x ^ 2 + 1))
lhs:
    apply 1 on INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)
    substitute pi / 4 - y for x
    simplify
    rewrite tan(pi / 4 - y) to (tan(pi / 4) - tan(y)) / (1 + tan(pi / 4) * tan(y))
    simplify
    rewrite (-tan(y) + 1) / (tan(y) + 1) + 1 to 2 / (1 + tan(y))
    rewrite log(2 / (1 + tan(y))) to log(2) - log(1 + tan(y))
    apply integral identity
    simplify
    apply 1 on INT x:[0,pi / 4]. log(tan(x) + 1)
done
from 2:
    solve equation for INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)
done

prove (INT t:[0,a]. log(t + a) / (t ^ 2 + a ^ 2)) = pi / (8 * a) * log(2 * a ^ 2) for a: real, a > 0
subgoal 1: (INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)) = a * (INT t:[0,a]. log(t + a) / (t ^ 2 + a ^ 2)) - pi / 4 * log(a)
lhs:
    substitute t / a for x
    simplify
    rewrite 1 / (t ^ 2 / a ^ 2 + 1) * log(t / a + 1) to log(t / a + 1) * a ^ 2 / (t ^ 2 + a ^ 2)
    rewrite t / a + 1 to (t + a) / a
    simplify
    rewrite log((a + t) / a) to log(a + t) - log(a)
    rewrite 1 / (a ^ 2 + t ^ 2) * (log(a + t) - log(a)) to log(a + t) / (a ^ 2 + t ^ 2) - log(a) / (a ^ 2 + t ^ 2)
    simplify
    apply integral identity
    simplify
    expand polynomial
done
subgoal 2: a * (INT t:[0,a]. log(t + a) / (t ^ 2 + a ^ 2)) - pi / 4 * log(a) = pi * log(2) / 8
lhs:
    apply 1 on a * (INT t:[0,a]. log(t + a) / (t ^ 2 + a ^ 2)) - pi / 4 * log(a)
    apply integral identity
done
from 2:
    solve equation for INT t:[0,a]. log(t + a) / (t ^ 2 + a ^ 2)
    rewrite pi * log(a) / 4 to 1/8 * pi * (2 * log(a))
    rewrite 2 * log(a) to log(a ^ 2)
    rewrite 1/8 * pi * log(a ^ 2) + pi * log(2) / 8 to 1/8 * pi * (log(2) + log(a ^ 2))
    rewrite log(2) + log(a ^ 2) to log(2 * a ^ 2)
done
