## Indefinite integrals, rational functions

// Source:
// Tongji 7'th edition
// Chapter 4, Section 4

calculate INT x. x ^ 3 / (x + 3)
    partial fraction decomposition
    apply integral identity
    simplify
done

calculate INT x. (2*x + 3) / (x^2 + 3*x - 10)
    partial fraction decomposition
    apply integral identity
    simplify
done

calculate INT x. (x + 1) / (x^2 - 2*x + 5)
    rewrite x^2 - 2*x + 5 to (x - 1)^2 + 4
    substitute u for (x - 1) / 2
    rewrite 4 * u ^ 2 + 4 to 4 * (u^2 + 1)
    simplify
    expand polynomial
    apply integral identity
    substitute v for u^2 + 1
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1 / (x * (x^2 + 1))
    partial fraction decomposition
    apply integral identity
    substitute u for x^2 + 1
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 3 / (x^3 + 1)
    partial fraction decomposition
    apply integral identity
    rewrite x^2 - x + 1 to (x - 1/2)^2 + 3/4
    substitute u for (x - 1/2)
    expand polynomial
    simplify
    substitute v for 2*u/sqrt(3)
    rewrite 3*v^2/2 + 3/2 to 3/2*(v^2+1)
    simplify
    apply integral identity
    substitute w for u^2 + 3/4
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. (x^2 + 1) / ((x+1)^2 * (x-1))
    partial fraction decomposition
    apply integral identity
    substitute u for x + 1
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. x / ((x + 1) * (x + 2) * (x + 3))
    partial fraction decomposition
    apply integral identity
    simplify
done

calculate INT x. (x^5 + x^4 - 8) / (x^3 - x)
    partial fraction decomposition
    apply integral identity
    simplify
done

calculate INT x. 1 / ((x^2 + 1) * (x^2 + x))
    partial fraction decomposition
    apply integral identity
    rewrite 2 * x^2 + 2 to 2 * (x^2 + 1)
    simplify
    expand polynomial
    apply integral identity
    substitute u for x^2 + 1
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1 / (x^4 - 1)
    partial fraction decomposition
    rewrite 2 * x^2 + 2 to 2 * (x^2 + 1)
    apply integral identity
    simplify
done

calculate INT x. 1 / ((x^2 + 1) * (x^2 + x + 1))
    partial fraction decomposition
    simplify
    rewrite x^2 + x + 1 to (x+1/2) ^ 2 + 3/4
    substitute u for x + 1/2
    expand polynomial
    simplify
    substitute v for 2*u/sqrt(3)
    rewrite 3*v^2/2 + 3/2 to 3/2 * (v^2 + 1)
    simplify
    apply integral identity
    substitute w for u^2 + 3/4
    apply integral identity
    substitute y for x^2 + 1
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. (x + 1)^2 / (x^2 + 1)^2
    partial fraction decomposition
    apply integral identity
    substitute u for x^2 + 1
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. (-x^2 - 2) / (x^2 + x + 1)^2
    partial fraction decomposition
    sorry

calculate INT x. 1 / (3 + sin(x)^2)
    rewrite 1 / (3 + sin(x)^2) to 1 / (3*sec(x)^2 + tan(x)^2) * (1 / cos(x)^2)
    rewrite sec(x)^2 to tan(x)^2 + 1
    rewrite 3 * (tan(x)^2 + 1) + tan(x)^2 to 3 + 4*tan(x)^2
    substitute u for tan(x)
    substitute v for 2*u/sqrt(3)
    rewrite 6 * v^2 + 6 to 6 * (v^2 + 1)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1 / (3 + cos(x))
    substitute u for tan(x/2)
    rewrite 2 / ((u ^ 2 + 1) * ((-(u ^ 2) + 1) / (u ^ 2 + 1) + 3)) to 1 / (2 + u^2)
    substitute v for u/sqrt(2)
    rewrite 2*v^2 + 2 to 2*(v^2 + 1)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1 / (2 + sin(x))
    substitute u for tan(x/2)
    rewrite 2 / ((u ^ 2 + 1) * (2 * u / (u ^ 2 + 1) + 2)) to 1 / (u^2 + u + 1)
    rewrite u^2 + u + 1 to (u + 1/2)^2 + 3/4
    substitute v for 2*(u + 1/2)/sqrt(3)
    rewrite 3 * v^2/2 + 3/2 to 3/2*(v^2 + 1)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1 / (1 + sin(x) + cos(x))
    substitute u for tan(x/2)
    rewrite 2 / ((u ^ 2 + 1) * ((-(u ^ 2) + 1) / (u ^ 2 + 1) + 2 * u / (u ^ 2 + 1) + 1)) to 1 / (1 + u)
    apply integral identity
    simplify
    apply integral identity
    replace substitution
done

calculate INT x. 1 / (2*sin(x) - cos(x) + 5)
    substitute u for tan(x/2)
    rewrite 2 / ((u ^ 2 + 1) * (-((-(u ^ 2) + 1) / (u ^ 2 + 1)) + 4 * u / (u ^ 2 + 1) + 5)) to 1 / (3*u^2 + 2*u + 2)
    rewrite 3 * u^2 + 2*u + 2 to 3 * ((u+1/3)^2 + 5/9)
    substitute v for 3 * (u + 1/3) / sqrt(5)
    rewrite 5 * v^2 + 5 to 5 * (v^2 + 1)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1 / (1 + (x + 1)^(1/3))
    substitute u for (x + 1)^(1/3)
    partial fraction decomposition
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. (x^(2/3) - 1) / (sqrt(x) + 1) for x > 0
    sorry

calculate INT x. (sqrt(x + 1) - 1) / (sqrt(x + 1) + 1)
    substitute t for sqrt(x + 1)
    partial fraction decomposition
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1 / (sqrt(x) + x^(1/4)) for x > 0
    substitute u for x^(1/4)
    simplify
    partial fraction decomposition
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. sqrt((1-x)/(1+x)) / x
    sorry

calculate INT x. 1 / ((x+1)^2*(x-1)^4) ^ (1/3)
    sorry
