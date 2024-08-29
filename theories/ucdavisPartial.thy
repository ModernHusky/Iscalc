## Examples on partial fractions

// Source:
// The Calculus Page Problems List by D. A. Kouba
// Problems on partial fractions 
// URL: https://www.math.ucdavis.edu/~kouba/CalcTwoDIRECTORY/partialfracdirectory/PartialFrac.html

calculate INT x:[3, 4]. 1 / (x ^ 2 - 4)
    partial fraction decomposition
    simplify
    substitute u for 4 * x + 8
    substitute u for 4 * x - 8 (at 2)
    apply integral identity
    simplify
done

calculate INT x:[1, 2]. (2 * x + 3) / (x ^ 2 - 9)
    partial fraction decomposition
    simplify
    substitute u for 2 * x - 6
    substitute u for 2 * x + 6 (at 2)
    apply integral identity
    simplify
done

calculate INT x:[1, 2]. (2 - x) / (x ^ 2 + 5 * x)
    partial fraction decomposition
    apply integral identity
    simplify
done

calculate INT x:[1, 2]. (x^2 - 1) / (x^2 - 16)
    partial fraction decomposition
    simplify
    apply integral identity
    simplify
done

calculate INT x:[3, 4]. (x^4 + x^3 + x^2 + 1)/(x^2 + x - 2)
    partial fraction decomposition
    simplify
    apply integral identity
    simplify
done

calculate INT x:[2, 4]. (x^2 + x - 1) / (x * (x^2 - 1))
    partial fraction decomposition
    simplify
    apply integral identity
    simplify
done

calculate INT x:[2, 4]. (x + 7) / (x ^ 2 * (x + 2))
    partial fraction decomposition
    apply integral identity
    simplify
done

calculate INT x:[2, 4]. (x^5 + 1) / (x ^ 3 * (x + 2))
    partial fraction decomposition
    apply integral identity
    simplify
done

calculate INT x:[0, 1]. (x ^ 2 - x + 1) / (x + 1)^3
    partial fraction decomposition
    substitute u for x + 1
    apply integral identity
    simplify
done

calculate INT x:[2, 3]. (x^3 + 4) / ((x ^ 2 - 1) * (x ^ 2 + 3 * x + 2))
    partial fraction decomposition
    simplify
    substitute u for 3 * x + 6
    substitute u for 4 * x + 4 (at 2)
    substitute u for 12 * x - 12 (at 3)
    substitute u for x + 1 (at 4)
    apply integral identity
    simplify
done

calculate INT x:[3, 4]. (x ^ 3 + 2 * x - 1) / (x ^ 2 - x - 2) ^ 2
    partial fraction decomposition
    simplify
    substitute u for 27 * x + 27
    substitute u for 27 * x - 54 (at 2)
    substitute u for x + 1 (at 3)
    substitute u for x - 2 (at 4)
    apply integral identity
    simplify
done

calculate INT x:[pi/12, pi/6]. sec(x) ^ 2 / (tan(x) ^ 3 - tan(x) ^ 2)
    substitute u for tan(x)
    partial fraction decomposition
    apply integral identity
    simplify
done

calculate INT x:[3, 4]. (x ^ 3 + 8) / ((x ^ 2 - 1) * (x - 2))
    partial fraction decomposition
    simplify
    substitute u for 6 * x + 6
    substitute u for 2 * x - 2 (at 2)
    substitute u for 3 * x - 6 (at 3)
    apply integral identity
    simplify
done

calculate INT x:[1, 2]. exp(x) / ((exp(x) - 1) * (exp(x) + 3))
    substitute u for exp(x)
    partial fraction decomposition
    simplify
    substitute v for 4 * u + 12
    substitute v for 4 * u - 4 (at 2)
    apply integral identity
    simplify
done

calculate INT x:[0, 1]. 1 / (exp(x) + 1)
    rewrite 1 / (exp(x) + 1) to exp(x) / (exp(x) * (exp(x) + 1))
    substitute u for exp(x)
    partial fraction decomposition
    apply integral identity
    simplify
done

calculate INT x:[1, 2]. (3 - x) / (x * (x ^2 + 1))
    partial fraction decomposition
    apply integral identity
    simplify
    rewrite (3 * x + 1) / (x ^ 2 + 1) to 3 * x / (x ^ 2 + 1) + 1 / (x ^ 2 + 1)
    apply integral identity
    simplify
    substitute u for x ^ 2 + 1
    apply integral identity
    simplify
done

calculate INT x:[1, 2]. (3 * x + 1) / (x ^ 2 * (x ^ 2 + 25))
    partial fraction decomposition
    apply integral identity
    simplify
    expand polynomial
    simplify
    substitute u for 25 * x ^ 2 + 625
    apply integral identity
    simplify
done

calculate INT x:[0, 1]. 1 / (x ^ 4 - 16)
    partial fraction decomposition
    simplify
    substitute u for x / 2
    rewrite 2 / (32 * u ^ 2 + 32) to 1/16 * (1 / (u ^ 2 + 1))
    apply integral identity
    simplify
done

calculate INT x:[pi/6, pi/3]. cos(x) / (sin(x) ^ 3 + sin(x))
    substitute u for sin(x)
    partial fraction decomposition
    simplify
    substitute v for u ^ 2 + 1
    apply integral identity
    simplify
done

calculate INT x:[0, 1]. 1 / (x ^ 4 + 4)
    partial fraction decomposition
    simplify
    rewrite 8 * x ^ 2 + 16 * x + 16 to 8 * ((x + 1) ^ 2 + 1)
    substitute u for x + 1
    rewrite 8 * x ^ 2 - 16 * x + 16 to 8 * ((x - 1) ^ 2 + 1)
    substitute u for x - 1 (at 2)
    expand polynomial
    simplify
    substitute v for u ^ 2 + 1
    apply integral identity
    simplify
    substitute w for u^2
    apply integral identity
    simplify
done
