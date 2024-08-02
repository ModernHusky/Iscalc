# 4.2 Substitution
# Section B
// page 175

calculate INT x. exp(exp(x)*sin(x))*(sin(x)+cos(x))*exp(x) for cos(x) + sin(x) != 0
    substitute u for exp(x)*sin(x)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. exp(asin(sqrt(x)))/sqrt(x-x^2) for x < 1, x > 0
    rewrite sqrt(x-x^2) to sqrt(x) * sqrt(1-x)
    substitute u for asin(sqrt(x))
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. exp(tan(1/x))/x^2 * sec(1/x)^2
    substitute u for tan(1/x)
    apply integral identity
    replace substitution
    simplify
done
