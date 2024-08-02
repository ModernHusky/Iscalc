## Examples on trigonometric functions

// Source:
// The Calculus Page Problems List by D. A. Kouba
// Problems on trigonometric functions
// URL: https://www.math.ucdavis.edu/~kouba/CalcTwoDIRECTORY/trigintdirectory/TrigInt.html

calculate INT x. tan(5 * x) for x > 0, x < pi/10
    rewrite tan(5 * x) to sin(5 * x) / cos(5 * x)
    substitute u for cos(5 * x)
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 5 * sec(4 * x) * tan(4 * x)
    rewrite sec(4 * x) to 1 / cos(4 * x)
    rewrite tan(4 * x) to sin(4 * x) / cos(4 * x)
    simplify
    substitute u for cos(4 * x)
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. (sin(x) + cos(x)) ^ 2
    expand polynomial
    rewrite 2 * cos(x) * sin(x) + cos(x) ^ 2 + sin(x) ^ 2 to 2 * cos(x) * sin(x) + (sin(x) ^ 2 + cos(x) ^ 2)
    rewrite sin(x) ^ 2 + cos(x) ^ 2 to 1
    simplify
    substitute u for sin(x)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 3 * cos(5 * x) ^ 2
    substitute u for 5 * x
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. (2 + tan(x)) ^ 2
    expand polynomial
    simplify
    apply integral identity
    rewrite tan(x) ^ 2 to sec(x) ^ 2 - 1
    simplify
    apply integral identity
    simplify
done

calculate INT x. sin(x) ^ 3
    rewrite sin(x) ^ 3 to sin(x) * sin(x) ^ 2
    rewrite sin(x) ^ 2 to 1 - cos(x) ^ 2
    expand polynomial
    apply integral identity
    substitute u for cos(x)
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. cos(5 * x) / (3 + sin(5 * x)) for x > 0, x < pi/10
    substitute u for sin(5 * x)
    substitute v for u + 3
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. cos(x) ^ 2 / (1 + sin(x))
    rewrite cos(x) ^ 2 to 1 - sin(x) ^ 2
    rewrite 1 - sin(x) ^ 2 to (1 + sin(x)) * (1 - sin(x))
    simplify
    apply integral identity
    simplify
done

calculate INT x. sin(x) / (1 + sin(x))
    rewrite sin(x) / (1 + sin(x)) to sin(x) * (1 - sin(x)) / (1 - sin(x) ^ 2)
    rewrite 1 - sin(x) ^ 2 to cos(x) ^ 2
    rewrite sin(x) * (1 - sin(x)) / cos(x) ^ 2 to 1 / cos(x) * (sin(x) / cos(x)) - (sin(x) / cos(x)) ^ 2
    rewrite 1 / cos(x) to sec(x)
    rewrite sin(x) / cos(x) to tan(x)
    rewrite sin(x) / cos(x) to tan(x)
    rewrite tan(x) ^ 2 to sec(x) ^ 2 - 1
    simplify
    apply integral identity
    simplify
done

calculate INT x. (csc(3 * x) + cot(3 * x)) ^ 2
    expand polynomial
    rewrite cot(3 * x) ^ 2 to csc(3 * x) ^ 2 - 1
    simplify
    substitute u for 3 * x
    simplify
    apply integral identity
    substitute v for 3 * x
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. sec(x) ^ 2 * sqrt(5 + tan(x))
    substitute u for 5 + tan(x)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. tan(x) ^ 5
    rewrite tan(x) ^ 5 to tan(x) ^ 3 * tan(x) ^ 2
    rewrite tan(x) ^ 2 to sec(x) ^ 2 - 1
    expand polynomial
    simplify
    substitute u for tan(x)
    apply integral identity
    simplify
    rewrite tan(x) ^ 3 to tan(x) * tan(x) ^ 2
    rewrite tan(x) ^ 2 to sec(x) ^ 2 - 1
    rewrite tan(x) * (sec(x) ^ 2 - 1) to tan(x) * sec(x) ^ 2 - tan(x)
    simplify
    substitute u for tan(x)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. (sin(x) - cos(x)) * (sin(x) + cos(x)) ^ 5
    rewrite (sin(x) - cos(x)) * (sin(x) + cos(x)) ^ 5 to -((sin(x) + cos(x)) ^ 5) * (cos(x) - sin(x))
    substitute u for sin(x) + cos(x)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. cos(x) * exp(4 + sin(x))
    rewrite cos(x) * exp(4 + sin(x)) to exp(4 + sin(x)) * cos(x)
    substitute u for 4 + sin(x)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. sin(3 * x) * sin(cos(3 * x))
    rewrite sin(3 * x) * sin(cos(3 * x)) to -sin(cos(3 * x)) / 3 * -(3 * sin(3 * x))
    substitute u for cos(3 * x)
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. cos(x) * log(sin(x)) / sin(x)
    rewrite cos(x) * log(sin(x)) / sin(x) to log(sin(x)) * (cos(x) / sin(x))
    substitute u for log(sin(x))
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. sec(x) * tan(x) * sqrt(4 + 3 * sec(x))
    rewrite sec(x) * tan(x) * sqrt(4 + 3 * sec(x)) to sqrt(4 + 3 * sec(x)) / 3 * (3 * sec(x) * tan(x))
    substitute u for 4 + 3 * sec(x)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. exp(x) * cos(exp(x))
    substitute u for exp(x)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. x * sin(3 * x)
    integrate by parts with u = x, v = -1/3 * cos(3 * x)
    simplify
    apply integral identity
    simplify
done

calculate INT x. x ^ 2 * cos(x)
    integrate by parts with u = x ^ 2, v = sin(x)
    simplify
    integrate by parts with u = x, v = -cos(x)
    simplify
    apply integral identity
    simplify
done

calculate INT x. sin(x) * cos(x) * exp(sin(x))
    rewrite sin(x) * cos(x) * exp(sin(x)) to sin(x) * exp(sin(x)) * cos(x)
    substitute u for sin(x)
    integrate by parts with u = u, v = exp(u)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. exp(x) * sin(x)
    integrate by parts with u = exp(x), v = -cos(x)
    simplify
    integrate by parts with u = exp(x), v = sin(x)
    solve integral INT x. exp(x) * sin(x)
done

calculate INT x. cos(4 * x) * sin(3 * x)
    integrate by parts with u = sin(3 * x), v = sin(4 * x) / 4
    simplify
    integrate by parts with u = cos(3 * x), v = -cos(4 * x) / 4
    simplify
    solve integral INT x. cos(4 * x) * sin(3 * x)
done

calculate INT x. sec(x) * sqrt(sec(x) + tan(x))
    rewrite sec(x) * sqrt(sec(x) + tan(x)) to 2 * ((sec(x) * tan(x) + sec(x) ^ 2) / (2 * sqrt(sec(x) + tan(x))))
    substitute u for sqrt(sec(x) + tan(x))
    apply integral identity
    replace substitution
done

calculate INT x. (sin(2 * x) - cos(2 * x)) / (sin(2 * x) + cos(2 * x))
    rewrite (sin(2 * x) - cos(2 * x)) / (sin(2 * x) + cos(2 * x)) to -1 / (2 * (sin(2 * x) + cos(2 * x))) * (2 * cos(2 * x) - 2 * sin(2 * x))
    substitute u for sin(2 * x) + cos(2 * x)
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. (sin(x) + cos(x)) / (exp(-x) + sin(x))
    rewrite (sin(x) + cos(x)) / (exp(-x) + sin(x)) to 1 / (1 + exp(x) * sin(x)) * (cos(x) * exp(x) + exp(x) * sin(x))
    substitute u for 1 + exp(x) * sin(x)
    apply integral identity
    replace substitution
    simplify
done
