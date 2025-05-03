imports interesting3

# Inside Interesting Integrals, Chapter 4

## Chapter 4, Section 1, Euler's Gamma functions

define Gamma(n) = (INT x:[0,oo]. exp(-x) * x^(n-1)) for n: real, n > 0

prove [bidirectional] Gamma(n) = (n - 1) * Gamma(n - 1) for n: real, n > 1
lhs:
    expand definition for Gamma
    integrate by parts with u = x ^ (n - 1), v = -exp(-x)
    simplify
rhs:
    expand definition for Gamma (all)
done

prove Gamma(n) = factorial(n - 1) for n: int, n >= 1
induction on n starting from 1
    base:
    lhs:
        expand definition for Gamma
        apply integral identity
        simplify
    done
    induct:
    lhs:
        rewrite Gamma(n + 1) to n * Gamma(n)
        apply induction hypothesis (all)
        rewrite n * factorial(n - 1) to factorial(n)
    done
done

calculate INT x:[0,oo]. exp(-(x ^ 3))
    substitute y for x ^ 3
    simplify
    rewrite exp(-y) / y ^ (2/3) to exp(-y) * y ^ (1/3 - 1)
    fold definition for Gamma (all)
    rewrite to (4/3 - 1) * Gamma(4/3 - 1)
    rewrite (4/3 - 1) * Gamma(4/3 - 1) to Gamma(4/3)
done

## Chapter 4, Section 2, Wallis' Integral and the Beta functions

define Beta(m,n) = Gamma(m) * Gamma(n) / Gamma(m+n)

