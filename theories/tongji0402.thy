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
