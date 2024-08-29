## Examples on exponential functions

// Source:
// The Calculus Page Problems List by D. A. Kouba
// Problems on exponential functions
// URL: https://www.math.ucdavis.edu/~kouba/CalcTwoDIRECTORY/expondirectory/Exponentials.html

calculate INT x. 5 * exp(x)
    apply integral identity
done

calculate INT x. 2 - 3 * exp(x)
    apply integral identity
    simplify
done

calculate INT x:[0,log(2) / 7]. 14 * exp(7 * x)
    apply integral identity
    simplify
done

calculate INT x. 7 ^ (2 * x + 3)
    substitute u for 2 * x
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. exp(5 * x) * (exp(2 * x) / 7 + 3 / exp(3 * x))
    expand polynomial
    apply integral identity
    simplify
done

calculate INT x. exp(x) * (1 + 2 * exp(x)) ^ 4
    substitute u for 1 + 2 * exp(x)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. (exp(4 * x) - exp(-4 * x)) ^ 2
    expand polynomial
    apply integral identity
    simplify
done

calculate INT x. exp(x) * (1 - exp(x)) * (1 + exp(x)) ^ 10
    substitute u for exp(x) + 1
    simplify
    expand polynomial
    apply integral identity
    replace substitution
    simplify
done

calculate INT x:[0,1]. (3 ^ x + 4 ^ x) / 5 ^ x
    rewrite (3 ^ x + 4 ^ x) / 5 ^ x to 3 ^ x / 5 ^ x + 4 ^ x / 5 ^ x
    rewrite 3 ^ x / 5 ^ x to (3/5) ^ x
    rewrite 4 ^ x / 5 ^ x to (4/5) ^ x
    simplify
    apply integral identity
    simplify
done

calculate INT x. 30 * exp(-3 * x) * (1 + 3 * exp(-x)) ^ 5
    substitute u for 3 * exp(-x) + 1
    simplify
    rewrite u ^ 5 * (u - 1) ^ 2 to u ^ 7 - 2 * u ^ 6 + u ^ 5
    simplify
    apply integral identity
    replace substitution
    simplify
done

prove (INT x. (27 * exp(9 * x) + exp(12 * x)) ^ (1/3)) = (exp(3 * x) + 27) ^ (4/3) / 4 + SKOLEM_CONST(C)
lhs:
    rewrite 27 * exp(9 * x) + exp(12 * x) to exp(9 * x) * (27 + exp(3 * x))
    rewrite (exp(9 * x) * (27 + exp(3 * x))) ^ (1/3) to exp(3 * x) * (27 + exp(3 * x)) ^ (1/3)
    substitute u for 27 + exp(3 * x)
    simplify
    apply integral identity
    simplify
    replace substitution
done
