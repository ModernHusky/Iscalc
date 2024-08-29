## Examples on power substitution

// Source:
// The Calculus Page Problems List by D. A. Kouba
// Problems on power substitution
// URL: https://www.math.ucdavis.edu/~kouba/CalcTwoDIRECTORY/powersubdirectory/PowerSub.html

calculate INT x. 1 / (1 + sqrt(x))
    substitute u for sqrt(x)
    partial fraction decomposition
    simplify
    substitute v for u + 1
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. (2 + sqrt(x)) / (3 - sqrt(x))
    substitute u for sqrt(x)
    partial fraction decomposition
    simplify
    substitute v for u - 3
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 3 / (4 + x ^ (1/3))
    substitute u for x ^ (1/3)
    simplify
    partial fraction decomposition
    simplify
    substitute v for u + 4
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. x * (x - 1) ^ (1/6) for x > 1
    substitute u for (x - 1) ^ (1/6)
    simplify
    expand polynomial
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. (3 * x + 2) / sqrt(x - 9) for x > 9
    substitute u for sqrt(x - 9)
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1 / (x ^ (2/3) - x ^ (1/3))
    substitute u for x ^ (1/3)
    partial fraction decomposition
    substitute v for u - 1
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. (sqrt(x) + 1) / (sqrt(x) * (x ^ (1/3) + 1)) for x > 0
    substitute u for x ^ (1/6)
    simplify
    partial fraction decomposition
    simplify
    expand polynomial
    simplify
    substitute v for u ^ 2 + 1
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. sqrt(5 + sqrt(x)) for x >= 0
    substitute u for sqrt(5 + sqrt(x))
    simplify
    expand polynomial
    simplify
    apply integral identity
    simplify
    replace substitution
    simplify
done
    
calculate INT x. (1 + sqrt(x - 3)) ^ (1/3) for x > 3
    substitute u for (1 + sqrt(x - 3)) ^ (1/3)
    simplify
    expand polynomial
    simplify
    apply integral identity
    simplify
    replace substitution
    simplify
done

calculate INT x. sqrt(2 + sqrt(4 + sqrt(x))) for x > 0
    substitute u for sqrt(2 + sqrt(4 + sqrt(x)))
    simplify
    expand polynomial
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1 / sqrt(2 + sqrt(1 + sqrt(x))) for x > 0
    substitute u for sqrt(2 + sqrt(1 + sqrt(x)))
    simplify
    expand polynomial
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. sqrt(x) / (x - 1)
    substitute u for sqrt(x)
    simplify
    partial fraction decomposition
    simplify
    substitute v for 2 * u + 2
    substitute w for 2 * u - 2 (at 2)
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. sqrt(4 - x) / x ^ 2 for x < 4
    substitute u for sqrt(4 - x)
    simplify
    partial fraction decomposition
    simplify
    substitute v1 for 8 * u + 16
    substitute v2 for 8 * u - 16 (at 2)
    substitute v3 for u + 2 (at 3)
    substitute v4 for u - 2 (at 4)
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. (x ^ (1/4) + 5) / (x - 16) for x > 0
    substitute u for x ^ (1/4)
    simplify
    partial fraction decomposition
    simplify
    expand polynomial
    simplify
    substitute v1 for 2 * u ^ 2 + 8
    substitute v2 for u / 2 (at 2)
    rewrite 8 * v2 ^ 2 + 8 to 8 * (v2 ^ 2 + 1)
    substitute v3 for 4 * u + 8 (at 3)
    substitute v4 for 4 * u - 8 (at 4)
    simplify
    apply integral identity
    replace substitution
    simplify
done
