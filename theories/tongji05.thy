## Definite integrals

// Source:
// Tongji 7'th edition
// Chapter 5

calculate INT x:[2,3]. 2 * x + x ^ 2
    apply integral identity
    simplify
done

calculate INT x:[0,1]. (3 * x + 1) ^ (-2)
    substitute u for 3 * x + 1
    apply integral identity
    simplify
done

calculate INT x:[0,1]. exp(6 * x)
    apply integral identity
    simplify
done

calculate INT x:[-1,2]. x * exp(x)
    integrate by parts with u = x, v = exp(x)
    apply integral identity
    simplify
done

calculate INT x:[0,pi/4]. sin(x)
    apply integral identity
    simplify
done

calculate INT x:[0,1]. 3*x^2 - x + 1
    apply integral identity
    simplify
done

calculate INT x:[1,2]. x^2 + 1/x^4
    apply integral identity
    simplify
done

calculate INT x:[pi/3, pi]. sin(2*x + pi/3)
    substitute u for 2*x + pi/3
    apply integral identity
    simplify
done

calculate INT x:[4, 9]. x ^ (1 / 3) * (x ^ (1 / 2) + 1)
    expand polynomial
    apply integral identity
    simplify
done

calculate INT x:[-1, 0]. (3 * x ^ 4 + 3 * x ^ 2 + 1) / (x ^ 2 + 1)
    partial fraction decomposition
    apply integral identity
    simplify
done

calculate INT x:[4, exp(1) + 3]. (x ^ 3 - 12 * x ^ 2 - 42) / (x - 3)
    partial fraction decomposition
    apply integral identity
    simplify
done

calculate INT x:[0, pi / 2]. sin(x) * cos(x) ^ 3
    substitute u for cos(x)
    apply integral identity
    simplify
done

calculate INT x:[0, pi]. 1 - sin(x) ^ 3
    simplify
    rewrite sin(x) ^ 3 to sin(x) * sin(x) ^ 2
    rewrite sin(x) ^ 2 to 1 - cos(x) ^ 2
    substitute u for cos(x)
    apply integral identity
    simplify
done

calculate INT x:[pi/6, pi/2]. cos(x) ^ 2
    apply integral identity
    simplify
done

calculate INT x:[0, 1]. (1 - x^2) ^ (1/2)
    substitute sin(u) for x
    rewrite 1 - sin(u) ^ 2 to cos(u) ^ 2
    simplify
    apply integral identity
    simplify
done

calculate INT x:[0, sqrt(2)]. sqrt(2 - x^2)
    substitute sqrt(2) * sin(u) for x
    simplify
    rewrite sin(u) ^ 2 to 1 - cos(u) ^ 2
    rewrite -(2 * (1 - cos(u) ^ 2)) + 2 to 2 * cos(u)^2
    simplify
    apply integral identity
    simplify
done

calculate INT y:[-sqrt(2), sqrt(2)]. sqrt(8 - 2*y^2)
    substitute 2 * sin(u) for y
    simplify
    rewrite sin(u) ^ 2 to 1 - cos(u) ^ 2
    rewrite -(8 * (1 - cos(u) ^ 2)) + 8 to 8*cos(u)^2
    simplify
    apply integral identity
    expand polynomial
done

calculate INT x:[1/sqrt(2), 1]. sqrt(1 - x^2) / x ^ 2
    substitute sin(u) for x
    simplify
    rewrite sin(u) ^ 2 to 1 - cos(u) ^ 2
    simplify
    rewrite cos(u) ^ 2 to 1 - sin(u) ^ 2
    expand polynomial
    simplify
    rewrite 1 / sin(u) ^ 2 to csc(u) ^ 2
    apply integral identity
    simplify
done

calculate INT x:[-1, 1]. x / sqrt(5 - 4 * x)
    substitute u for 5 - 4 * x
    expand polynomial
    simplify
    apply integral identity
    simplify
done

calculate INT x:[1,4]. 1 / (1 + sqrt(x))
    substitute u for sqrt(x)
    substitute v for u + 1
    expand polynomial
    simplify
    apply integral identity
    simplify
done

calculate INT x:[3/4, 1]. 1 / (sqrt(1-x) - 1)
    substitute u for sqrt(1 - x)
    substitute v for u - 1
    expand polynomial
    simplify
    apply integral identity
    simplify
done

calculate INT t:[0, 1]. t * exp(-t ^ 2 / 2)
    substitute u for t ^ 2 / 2
    apply integral identity
    simplify
done

calculate INT x:[1, exp(2)]. 1 / (x * sqrt(1 + log(x)))
    substitute u for 1 + log(x)
    apply integral identity
    simplify
done

calculate INT x:[-2, 0]. (x + 2) / (x^2 + 2*x + 2)
    rewrite x^2 + 2*x + 2 to (x + 1) ^ 2 + 1
    substitute u for x + 1
    expand polynomial
    simplify
    split region at 0
    substitute v for u ^ 2 + 1
    apply integral identity
    substitute v for u ^ 2 + 1
    apply integral identity
    simplify
done

calculate INT x:[-pi/2, pi/2]. cos(x) ^ 4
    rewrite cos(x) ^ 4 to (cos(x) ^ 2) ^ 2
    rewrite cos(x) ^ 2 to (1 + cos(2*x)) / 2
    expand polynomial
    substitute u for 2 * x
    simplify
    apply integral identity
    simplify
done

calculate INT x:[-pi/2, pi/2]. sqrt(cos(x) - cos(x)^3)
    rewrite cos(x) - cos(x)^3 to cos(x) * (1 - cos(x)^2)
    rewrite 1 - cos(x)^2 to sin(x)^2
    simplify
    split region at 0
    simplify
    substitute u for cos(x)
    apply integral identity
    substitute u for cos(x)
    apply integral identity
    simplify
done

calculate INT x:[0, pi]. sqrt(1 + cos(2*x))
    rewrite cos(2 * x) to 2 * cos(x) ^ 2 - 1
    simplify
    split region at pi / 2
    simplify
    apply integral identity
    simplify
done

calculate INT x:[0, 1]. x * exp(-x)
    integrate by parts with u = x, v = -exp(-x)
    simplify
    apply integral identity
    simplify
done

calculate INT x:[1, exp(1)]. x * log(x)
    integrate by parts with u = log(x) / 2, v = x ^ 2
    apply integral identity
    simplify
done

calculate INT x:[pi/4, pi/3]. x / sin(x)^2
    rewrite sin(x) ^ 2 to csc(x) ^ -2
    integrate by parts with u = x, v = -cot(x)
    simplify
    rewrite cot(x) to cos(x) / sin(x)
    substitute u for sin(x)
    apply integral identity
    simplify
done

calculate INT x:[1, 4]. log(x) / sqrt(x)
    integrate by parts with u = 2 * log(x), v = sqrt(x)
    apply integral identity
    simplify
done

calculate INT x:[0, 1]. x * arctan(x)
    integrate by parts with u = arctan(x) / 2, v = x ^ 2
    simplify
    rewrite x^2 / (2 * x^2 + 2) to (1 - 1 / (x^2 + 1)) / 2
    apply integral identity
    simplify
done

calculate INT x:[0, pi/2]. exp(2*x) * cos(x)
    integrate by parts with u = exp(2*x), v = sin(x)
    simplify
    integrate by parts with u = exp(2*x), v = -cos(x)
    simplify
    solve integral INT x:[0, pi/2]. exp(2*x)*cos(x)
done

calculate INT x:[0,pi]. (x * sin(x)) ^ 2
    simplify
    rewrite sin(x) ^ 2 to (1 - cos(2*x)) / 2
    expand polynomial
    simplify
    integrate by parts with u = x^2 / 2, v = sin(2*x)
    simplify
    integrate by parts with u = x / 2, v = -cos(2*x)
    simplify
    apply integral identity
    simplify
done

calculate INT x:[1, exp(1)]. sin(log(x))
    substitute u for log(x)
    integrate by parts with u = -exp(u), v = cos(u)
    simplify
    integrate by parts with u = exp(u), v = sin(u)
    simplify
    solve integral INT u:[0,1]. exp(u) * sin(u)
done

calculate INT x:[1/exp(1), exp(1)]. abs(log(x))
    split region at 1
    simplify
    integrate by parts with u = log(x), v = x
    simplify
    integrate by parts with u = log(x), v = x
    apply integral identity
    simplify
done
