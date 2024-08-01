## Indefinite integrals, integration by parts

// Source:
// Tongji 7'th edition
// Chapter 4, Section 3

calculate INT x. x * sin(x)
    integrate by parts with u = x, v = -cos(x)
    apply integral identity
    simplify
done

calculate INT x. log(x) for x > 0
    integrate by parts with u = log(x), v = x
    simplify
    apply integral identity
    simplify
done

calculate INT x. asin(x)
    integrate by parts with u = asin(x), v = x
    substitute u for -(x^2) + 1
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. x * exp(-x)
    integrate by parts with u = x, v = -exp(-x)
    apply integral identity
    simplify
done

calculate INT x. x^2 * log(x)
    integrate by parts with u = log(x), v = 1/3*x^3
    apply integral identity
    simplify
done

calculate INT x. exp(-x) * cos(x)
    integrate by parts with u = exp(-x), v = sin(x)
    simplify
    integrate by parts with u = exp(-x), v = -cos(x)
    simplify
    solve integral INT x. exp(-x) * cos(x)
    simplify
done

calculate INT x. exp(-2*x) * sin(x/2)
    integrate by parts with u = exp(-2*x), v = -2*cos(x/2)
    simplify
    integrate by parts with u = exp(-2*x), v = 2*sin(x/2)
    simplify
    solve integral INT x. exp(-2*x) * sin(x/2)
    simplify
done

calculate INT x. x * cos(x/2)
    integrate by parts with u = x, v = 2*sin(x/2)
    substitute u for x / 2
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. x^2 * atan(x)
    integrate by parts with u = atan(x), v = 1/3*x^3
    partial fraction decomposition
    apply integral identity
    simplify
    substitute u for 3*x^2 + 3
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. x * tan(x)^2
    rewrite tan(x)^2 to sec(x)^2 - 1
    expand polynomial
    apply integral identity
    integrate by parts with u = x, v = tan(x)
    apply integral identity
    simplify
done

calculate INT x. x^2 * cos(x)
    integrate by parts with u = x^2, v = sin(x)
    simplify
    integrate by parts with u = x, v = -cos(x)
    apply integral identity
    simplify
done

calculate INT t. t * exp(-t)
    integrate by parts with u = t, v = -exp(-t)
    apply integral identity
    simplify
done

calculate INT x. log(x)^2 for x > 0
    integrate by parts with u = log(x)^2, v = x
    simplify
    integrate by parts with u = log(x), v = x
    apply integral identity
    simplify
done

calculate INT x. x * sin(x) * cos(x)
    rewrite x * sin(x) * cos(x) to x * (2 * sin(x) * cos(x)) / 2
    rewrite 2 * sin(x) * cos(x) to sin(2*x)
    integrate by parts with u = x/4, v = -cos(2*x)
    simplify
    apply integral identity
    simplify
done

calculate INT x. x^2 * cos(x/2)^2
    rewrite cos(x/2)^2 to (1+cos(x)) / 2
    expand polynomial
    apply integral identity
    integrate by parts with u = x^2, v = sin(x)
    simplify
    integrate by parts with u = x, v = -cos(x)
    simplify
    apply integral identity
    simplify
done

calculate INT x. x * log(x - 1) for x > 1
    integrate by parts with u = log(x - 1), v = x^2/2
    simplify
    substitute u for 2*x-2
    expand polynomial
    simplify
    apply integral identity
    replace substitution
    expand polynomial
    simplify
done

calculate INT x. (x^2-1) * sin(2*x)
    expand polynomial
    apply integral identity
    integrate by parts with u = x^2, v = -cos(2*x)/2
    simplify
    integrate by parts with u = x, v = sin(2*x)/2
    simplify
    apply integral identity
    simplify
done

calculate INT x. log(x)^3 / x^2
    integrate by parts with u = log(x)^3, v = -1/x
    simplify
    integrate by parts with u = log(x)^2, v = -1/x
    simplify
    integrate by parts with u = log(x), v = -1/x
    simplify
    apply integral identity
    simplify
done

calculate INT x. exp(x^(1/3))
    substitute u^3 for x
    simplify
    integrate by parts with u = u^2, v = exp(u)
    simplify
    integrate by parts with u = u, v = exp(u)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. cos(log(x)) for x > 0
    integrate by parts with u = cos(log(x)), v = x
    simplify
    integrate by parts with u = sin(log(x)), v = x
    simplify
    solve integral INT x. cos(log(x))
done

calculate INT x. asin(x)^2 for x > -1, x < 1
    integrate by parts with u = asin(x)^2, v = x
    simplify
    integrate by parts with u = asin(x), v = -sqrt(-(x^2) + 1)
    simplify
    apply integral identity
    simplify
done

calculate INT x. exp(x) * sin(x)^2
    rewrite sin(x)^2 to (1 - cos(2*x)) / 2
    expand polynomial
    apply integral identity
    integrate by parts with u = exp(x), v = sin(2*x)/2
    simplify
    integrate by parts with u = exp(x), v = -cos(2*x)/2
    simplify
    sorry

calculate INT x. x * log(x)^2
    integrate by parts with u = log(x)^2, v = x^2/2
    integrate by parts with u = log(x), v = x^2/2
    apply integral identity
    simplify
done

calculate INT x. exp(sqrt(3*x+9))
    substitute u for sqrt(3*x+9)
    simplify
    integrate by parts with u = u, v = exp(u)
    simplify
    apply integral identity
    replace substitution
    simplify
done
