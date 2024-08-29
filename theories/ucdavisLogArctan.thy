## Examples on rational functions, resulting in log or arctangent functions

// Source:
// The Calculus Page Problems List by D. A. Kouba
// Problems on rational functions, resulting in log or actangent functions
// URL: https://www.math.ucdavis.edu/~kouba/CalcTwoDIRECTORY/logarctandirectory/LogArctan.html

calculate INT x:[1, 2]. (3 + x^2)/x
    partial fraction decomposition
    apply integral identity
    simplify
done

calculate INT x:[0, 1]. 7/(x + 5)
    substitute u for x + 5
    apply integral identity
    simplify
done

calculate INT x:[-1, 1]. 3/(x^2 + 4)
    substitute u for x / 2
    rewrite 4 * u ^ 2 + 4 to 4 * (u ^ 2 + 1)
    apply integral identity
    simplify
done

calculate INT x:[-1, 0]. x/(x^2 + 4)
    substitute u for x ^ 2 + 4
    apply integral identity
    simplify
done

calculate INT x:[-1, 1].  x^3/(x + 2)
    partial fraction decomposition
    apply integral identity
    simplify
done

calculate INT x:[-1, 1]. (x^4 + x^3)/(x^2 + 1)
    partial fraction decomposition
    apply integral identity
    simplify
    expand polynomial
    simplify
    split region at 0
    substitute u for x ^ 2 + 1
    substitute u for x ^ 2 + 1 (at 2)
    apply integral identity
    simplify
done

calculate INT x:[1, exp(1)]. 1 / (x * (3 + log(x)))
    substitute u for 3 + log(x)
    apply integral identity
    simplify
done

calculate INT x:[1, exp(1)]. (7 - log(x)) / (x * (3 + log(x)))
    substitute u for log(x)
    substitute v for u + 3
    expand polynomial
    apply integral identity
    simplify
done

calculate INT x:[0, 1]. 1 / (x^2 + 4 * x + 5)
    rewrite x ^ 2 + 4 * x + 5 to (x + 2) ^ 2 + 1
    substitute u for x + 2
    apply integral identity
    simplify
done

calculate INT x:[1, 2]. 1 / (2*x^2 + 12*x + 68)
    rewrite 2*x^2 + 12*x + 68 to 2 * (x+3)^2 + 50
    substitute u for x + 3
    substitute v for u / 5
    rewrite 50 * v ^ 2 + 50 to 50 * (v ^ 2 + 1)
    apply integral identity
    simplify
done

calculate INT x:[9, 16]. 1 / (sqrt(x) * (x + 9))
    substitute u for sqrt(x)
    substitute v for u / 3
    rewrite 9 * v^2 + 9 to 9 * (v^2 + 1)
    apply integral identity
    simplify
done

calculate INT x:[0, pi/4]. sin(2*x)/(1 + cos(2*x))
    substitute u for 1 + cos(2*x)
    apply integral identity
    simplify
done

calculate INT x:[0, pi/2]. cos(x) / (1 + sin(x) ^ 2)
    substitute u for sin(x)
    apply integral identity
    simplify
done

calculate INT x:[0, 1]. exp(2*x)/(1 + exp(2*x))
    substitute u for exp(2*x)
    substitute v for 2 * u + 2
    apply integral identity
    simplify
done

calculate INT x:[0, 1]. exp(x)/(1 + exp(2*x))
    rewrite exp(2*x) to exp(x) * exp(x)
    substitute u for exp(x)
    apply integral identity
    simplify
done

calculate INT x:[0, 1]. exp(4*x)/(1 + exp(2*x))
    substitute u for exp(2*x)
    partial fraction decomposition
    simplify
    substitute v for 2 * u + 2
    apply integral identity
    simplify
done

calculate INT x:[0, 1]. x / (x^2 + 2 * x + 1)
    rewrite x^2 + 2*x + 1 to (x+1)^2
    substitute u for x + 1
    expand polynomial
    apply integral identity
    simplify
done

calculate INT x:[0, 1]. x / (x ^ 2 + 2 * x + 5)
    rewrite x^2 + 2*x + 5 to (x+1)^2 + 4
    substitute u for x + 1
    expand polynomial
    simplify
    substitute v for u ^ 2 + 4
    apply integral identity
    substitute v for u / 2
    rewrite 4 * v ^ 2 + 4 to 4 * (v ^ 2 + 1)
    apply integral identity
    simplify
done

calculate INT x:[-1, 1]. (7 * x - 2) / (2*x^2 - 16 * x + 42)
    rewrite 2*x^2 - 16*x + 42 to 2*(x-4)^2 + 10
    substitute u for x - 4
    expand polynomial
    simplify
    substitute v for 2 * u ^ 2 + 10
    apply integral identity
    simplify
done

calculate INT x:[0, 1]. 1 / (1 + exp(2*x))
    rewrite 1 / (1 + exp(2*x)) to exp(-2*x) / (exp(-2*x) + 1)
    substitute u for exp(-2*x) + 1
    apply integral identity
    simplify
done

calculate INT x:[0, 1]. exp(3*x)/(1 + exp(2*x))
    rewrite exp(3*x) / (1 + exp(2*x)) to exp(x) * exp(x) * exp(x) / (1 + exp(x) * exp(x))
    substitute u for exp(x)
    partial fraction decomposition
    apply integral identity
    simplify
done

calculate INT x:[1, 4]. (9 + 6*sqrt(x) + x)/(4*sqrt(x) + x)
    rewrite 4*sqrt(x) + x to sqrt(x) * (4 + sqrt(x))
    substitute (u - 4)^2 for x
    substitute t for u - 4
    simplify
    partial fraction decomposition
    apply integral identity
    simplify
done
