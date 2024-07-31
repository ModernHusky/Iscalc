## Indefinite integrals, method of substitution

// Source:
// Tongji 7'th edition
// Chapter 4, Section 2

calculate INT t. exp(5*t)
    apply integral identity
done

calculate INT x. (3 - 2*x)^3
    substitute u for 3 - 2*x
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1 / (1 - 2*x)
    simplify
    apply integral identity
    simplify
done

calculate INT x. 1 / (2-3*x)^(1/3)
    substitute u for 2-3*x
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. sin(a*x) - exp(x/b)
    apply integral identity
    substitute u for x / b
    apply integral identity
    replace substitution
done

calculate INT t. sin(sqrt(t)) / sqrt(t) for t > 0
    substitute u for sqrt(t)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. x * exp(-x^2)
    substitute u for x^2
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. x * cos(x^2)
    substitute u for x^2
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. x / sqrt(2 - 3*x^2)
    substitute u for 2 - 3*x^2
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 3*x^3 / (1 - x^4)
    substitute u for 1 - x^4
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. (x + 1) / (x^2 + 2*x + 5)
    rewrite x + 1 to (2*x + 2) / 2
    substitute u for x^2 + 2*x + 5
    apply integral identity
    replace substitution
    simplify
done

calculate INT t. cos(w * t + phi) ^ 2 * sin(w * t + phi)
    simplify
    substitute u for cos(t * w + phi)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. sin(x) / cos(x)^3
    substitute u for cos(x)
    apply integral identity
    simplify
    replace substitution
    simplify
done

calculate INT x. (sin(x) + cos(x)) / (sin(x) - cos(x))^(1/3)
    rewrite sin(x) + cos(x) to cos(x) + sin(x)
    substitute u for sin(x) - cos(x)
    apply integral identity
    simplify
    replace substitution
    simplify
done

calculate INT x. tan(x)^10 * sec(x)^2
    substitute u for tan(x)
    apply integral identity
    simplify
    replace substitution
done

calculate INT x. 1 / (x * log(x) * log(log(x)))
    substitute u for log(log(x))
    apply integral identity
    replace substitution
done

calculate INT x. 1 / (asin(x)^2 * sqrt(1 - x^2)) for x > -1, x < 1
    substitute u for asin(x)
    rewrite -(sin(u)^2) + 1 to 1 - sin(u)^2
    rewrite 1 - sin(u)^2 to cos(u)^2
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 10^(2*acos(x)) / sqrt(1 - x^2) for x > -1, x < 1
    substitute u for acos(x)
    rewrite -(cos(u)^2) + 1 to 1 - cos(u)^2
    rewrite 1 - cos(u)^2 to sin(u)^2
    simplify
    substitute v for 2 * u
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. tan(sqrt(1 + x^2)) * x / sqrt(1 + x^2)
    simplify
    substitute u for sqrt(x^2 + 1)
    apply integral identity
    replace substitution
done

calculate INT x. atan(sqrt(x)) / (sqrt(x) * (1 + x)) for x > 0
    substitute u for atan(sqrt(x))
    simplify
    rewrite tan(u)^2 + 1 to sec(u)^2
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. (1 + log(x)) / (x * log(x)) ^ 2
    rewrite 1 + log(x) to log(x) + 1
    substitute u for x * log(x)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1 / (sin(x) * cos(x))
    rewrite 1 / (sin(x) * cos(x)) to 2 / (2 * sin(x) * cos(x))
    rewrite 2 * sin(x) * cos(x) to sin(2*x)
    rewrite 2 / sin(2*x) to 2 * csc(2*x)
    substitute u for 2 * x
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. log(tan(x)) / (cos(x) * sin(x))
    substitute u for log(tan(x))
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. cos(x)^3
    rewrite cos(x)^3 to cos(x) * cos(x)^2
    rewrite cos(x)^2 to 1 - sin(x)^2
    expand polynomial
    apply integral identity
    substitute u for sin(x)
    apply integral identity
    replace substitution
    simplify
done

calculate INT t. cos(w*t + phi)^2
    substitute u for w*t + phi
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. sin(2*x) * cos(3*x)
    rewrite sin(2*x) * cos(3*x) to 1/2 * (sin(5*x) - sin(1*x))
    apply integral identity
    simplify
done

calculate INT x. cos(x) * cos(x/2)
    rewrite cos(x) * cos(x/2) to 1/2 * (cos(3/2*x) + cos(1/2*x))
    apply integral identity
    simplify
done

calculate INT x. sin(5*x) * sin(7*x)
    rewrite sin(5*x) * sin(7*x) to -1/2 * (cos(-2*x) - cos(12*x))
    apply integral identity
    simplify
done

calculate INT x. tan(x)^3 * sec(x) for x > -pi/2, x < pi/2
    rewrite tan(x)^3 * sec(x) to tan(x)^2 * (tan(x) * sec(x))
    rewrite tan(x)^2 to sec(x)^2 - 1
    substitute u for sec(x)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1 / (exp(x) + exp(-x))
    rewrite 1 / (exp(x) + exp(-x)) to exp(x) / (exp(x)^2 + 1)
    substitute u for exp(x)
    apply integral identity
    replace substitution
done

calculate INT x. (1 - x) / sqrt(9 - 4 * x^2) for x > -3/2, x < 3/2
    substitute u for 2*x/3
    simplify
    rewrite sqrt(-(9 * u^2) + 9) to 3 * sqrt(-u^2 + 1)
    substitute sin(t) for u
    rewrite -(sin(t)^2) + 1 to 1 - sin(t)^2
    rewrite 1 - sin(t)^2 to cos(t)^2
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. x^3 / (9 + x^2) for x != 0
    partial fraction decomposition
    substitute u for x ^ 2 + 9
    expand polynomial
    simplify
    apply integral identity
    simplify
    replace substitution
    simplify
done

calculate INT x. 1 / (2*x^2 - 1)
    substitute u for sqrt(2) * x
    rewrite 2*u^2 - 2 to 2 * (u^2-1)
    simplify
    partial fraction decomposition
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1 / ((x + 1) * (x - 2))
    partial fraction decomposition
    apply integral identity
    simplify
done

calculate INT x. x / (x^2 - x - 2)
    partial fraction decomposition
    apply integral identity
    simplify
done

calculate INT x. x^2 / sqrt(a^2 - x^2) for a > 0, x > -a, x < a
    substitute a*sin(t) for x
    rewrite a^2 - (a*sin(t))^2 to a^2 * (1 - sin(t)^2)
    rewrite 1-sin(t)^2 to cos(t)^2
    simplify
    rewrite sin(t)^2 to (1 - cos(2*t)) / 2
    apply integral identity
    replace substitution
    expand polynomial
    simplify
done

calculate INT x. 1 / (x * sqrt(x^2 - 1)) for x > 1
    substitute t for sqrt(x^2-1)
    simplify
    apply integral identity
    replace substitution
done

calculate INT x. 1 / sqrt((x^2 + 1)^3)
    substitute tan(t) for x
    rewrite tan(t)^2 + 1 to sec(t)^2
    simplify
    rewrite 1 / sec(t) to cos(t)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. sqrt(x^2 - 9) / x for x > 3
    substitute 3*sec(u) for x
    simplify
    rewrite 9 * sec(u)^2 - 9 to 9 * (sec(u)^2 - 1)
    rewrite sec(u)^2 - 1 to tan(u)^2
    simplify
    rewrite tan(u)^2 to sec(u)^2 - 1
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1 / (1 + sqrt(2*x)) for x > 0
    substitute t for sqrt(2*x)
    simplify
    partial fraction decomposition
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1 / (1 + sqrt(1 - x^2)) for x > -1, x < 1
    substitute sin(t) for x
    rewrite 1 - sin(t)^2 to cos(t)^2
    simplify
    rewrite cos(t) to cos(2*(t/2))
    rewrite cos(2*(t/2)) to 2*cos(t/2)^2 - 1
    rewrite cos(t) to cos(2*(t/2))
    rewrite cos(2*(t/2)) to 2*cos(t/2)^2 - 1
    expand polynomial
    simplify
    rewrite 1 / cos(t/2)^2 to sec(t/2)^2
    substitute u for t / 2
    apply integral identity
    replace substitution
    simplify
    sorry

calculate INT x. 1 / (x + sqrt(1 - x^2)) for x > -1, x < 1
    substitute sin(t) for x
    rewrite 1 - sin(t)^2 to cos(t)^2
    simplify
    rewrite cos(t) / (cos(t) + sin(t)) to 1/2 * (1 + (cos(t) - sin(t)) / (cos(t) + sin(t)))
    simplify
    substitute u for cos(t) + sin(t)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. (x - 1) / (x^2 + 2*x + 3)
    rewrite x^2 + 2*x + 3 to (x + 1)^2 + 2
    substitute y for x + 1
    substitute z for y / sqrt(2)
    rewrite 2 * z^2 + 2 to 2 * (z^2 + 1)
    simplify
    expand polynomial
    apply integral identity
    substitute u for z^2 + 1
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. (x^3 + 1) / (x^2 + 1)^2
    substitute tan(t) for x
    rewrite tan(t)^2 + 1 to sec(t)^2
    simplify
    expand polynomial
    rewrite tan(t)^3 / sec(t)^2 to sin(t)^3 / cos(t)
    rewrite 1 / sec(t)^2 to cos(t)^2
    rewrite sin(t)^3 to sin(t) * sin(t)^2
    rewrite sin(t)^2 to 1 - cos(t)^2
    expand polynomial
    apply integral identity
    substitute u for cos(t)
    apply integral identity
    substitute u for cos(t)
    apply integral identity
    replace substitution
    simplify
done
