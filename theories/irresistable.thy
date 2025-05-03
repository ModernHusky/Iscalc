# Irresistable Integrals, Section 2.3

prove (INT x:[0,oo]. 1 / (x ^ 2 + b) ^ (m + 1)) = pi / 2 ^ (2 * m + 1) * binom(2 * m,m) * (1 / b ^ ((2 * m + 1) / 2)) for m: int, b: real, b > 0, m >= 0
let I(m,b) = (INT x:[0,oo]. 1 / (x ^ 2 + b) ^ (m + 1)) for b > 0, m >= 0
subgoal 1: (D b. I(m,b)) = -(m + 1) * I(m + 1,b)
lhs:
    expand definition for I (all)
    exchange derivative and integral
    simplify
rhs:
    expand definition for I (all)
    simplify
done
subgoal 2: I(m,b) = pi / 2 ^ (2 * m + 1) * binom(2 * m,m) * (1 / b ^ ((2 * m + 1) / 2)) for m: int
induction on m
    base:
        lhs:
            expand definition for I
            substitute sqrt(b) * u for x
            simplify
            rewrite 1 / (b * u ^ 2 + b) to 1 / b * (1 / (1 ^ 2 + u ^ 2))
            apply integral identity
            simplify
    done
    induct:
        lhs:
            apply 1 on I(m + 1,b)
            apply induction hypothesis (all)
            simplify
            rewrite -((2 * m + 1) / 2) - 1 to -m - 3/2
            rewrite to b ^ (-m - 3/2) * 2 ^ -(2 * m) * pi * (2 * m + 1) / (4 * m + 4) * binom(2 * m,m)
        rhs:
            rewrite binom(2 * m + 2,m + 1) to 2 * binom(2 * m,m) * ((2 * m + 1) / (m + 1))
            rewrite -((2 * m + 3) / 2) to -m - 3/2
            simplify
    done
done
lhs:
    fold definition for I (all)
    apply 2 on I(m,b)
done