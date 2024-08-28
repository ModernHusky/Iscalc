# 4.3 Inverse Substitution
# Section A
// page 180

calculate INT x. x^2 / sqrt(a^2-x^2) for a > 0, x > -a, x < a
    substitute a*sin(u) for x
    rewrite sqrt(a ^ 2 - (a * sin(u)) ^ 2) to a * sqrt(1-sin(u)^2)
    rewrite 1-sin(u)^2 to cos(u)^2
    simplify
    integrate by parts with u=-sin(u), v=cos(u)
    rewrite cos(u)^2 to 1-sin(u)^2
    simplify
    rewrite (INT u. sin(u) ^ 2) to 1/a^2 * (a ^ 2 * (INT u. sin(u) ^ 2))
    solve integral a ^ 2 * (INT u. sin(u) ^ 2)
    expand polynomial
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1 / sqrt(x^2 + 1)
    substitute tan(u) for  x
    rewrite tan(u)^2 + 1 to sec(u)^2
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. sqrt(x^2 - 9) / x for x > 3
    substitute 3 * sec(u) for x
    simplify
    rewrite sqrt(9 * sec(u) ^ 2 - 9) to 3 * sqrt(sec(u)^2 - 1)
    rewrite sec(u)^2 - 1 to tan(u)^2
    simplify
    rewrite tan(u)^2 to sec(u)^2 - 1
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1 / (x * sqrt(1-x^2)) for x > 0, x < 1
    substitute sin(u) for x
    rewrite 1 - sin(u)^2 to cos(u)^2
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1/(1+sqrt(2*x)) for x >= 0
    substitute t for sqrt(2*x)
    partial fraction decomposition
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1/(x^2*sqrt(1-x^2)) for x > 0, x < 1
    substitute sin(u) for x
    rewrite 1-sin(u)^2 to cos(u)^2
    simplify
    apply integral identity
    replace substitution
    simplify
done

// page 181

calculate INT x. 1/((x^2+1)*sqrt(1-x^2)) for x>0, x<1
    substitute sin(u) for x
    rewrite 1-sin(u)^2 to cos(u)^2
    simplify
    rewrite 1/(sin(u)^2+1) to (1/cos(u))^2 / ((sin(u)/cos(u))^2 + (1/cos(u))^2)
    rewrite 1/cos(u) to sec(u)
    rewrite 1/cos(u) to sec(u)
    rewrite sin(u)/cos(u) to tan(u)
    rewrite sec(u)^2 to 1+tan(u)^2 (at 2)
    simplify
    substitute v for tan(u)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1/((x^2+1)*sqrt(x^2-1)) for x > 1
    substitute sec(u) for x
    rewrite sec(u)^2-1 to tan(u)^2
    simplify
    rewrite sec(u) to 1/cos(u)
    rewrite sec(u) to 1/cos(u)
    rewrite 1 / cos(u) / ((1 / cos(u)) ^ 2 + 1) to cos(u) / (1+cos(u)^2)
    rewrite cos(u)^2 to 1-sin(u)^2
    substitute v for sin(u)
    substitute sqrt(2)*w for v
    simplify
    partial fraction decomposition
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1/(x+sqrt(1-x^2)) for x > 0, x < 1
    substitute sin(u) for x
    rewrite 1-sin(u)^2 to cos(u)^2
    simplify
    rewrite cos(u) / (cos(u) + sin(u)) to 1/2 * ((cos(u) + sin(u) + cos(u)-sin(u))/ (cos(u) + sin(u)))
    rewrite (cos(u) + sin(u) + cos(u)-sin(u))/ (cos(u) + sin(u)) to 1 + (cos(u)-sin(u)) / (cos(u)+sin(u))
    simplify
    apply integral identity
    substitute v for cos(u)+sin(u)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. (x^3+1) / (x^2+1)^2
    substitute tan(u) for x
    rewrite tan(u)^2+1 to  sec(u)^2
    simplify
    rewrite tan(u) to sin(u)/cos(u)
    rewrite sec(u) to 1/cos(u)
    expand polynomial
    simplify
    rewrite sin(u)^3 to sin(u)^2 * sin(u)
    rewrite sin(u)^2 to 1 - cos(u)^2
    substitute v for cos(u)
    simplify
    expand polynomial
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1 / (sqrt(x)+x^(1/4)) for x > 0
    substitute u for x^(1/4)
    simplify
    rewrite u^2+u to u*(u+1)
    simplify
    partial fraction decomposition
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. sqrt((1-x) / (x^2*(1+x))) for x > 0, x < 1
    rewrite (1-x) / (x^2*(1+x)) to ((1-x)*(1+x)) / (x^2 * (1+x)^2)
    rewrite (1-x)*(1+x) to 1-x^2
    rewrite sqrt((1-x^2) / (x ^ 2 * (1 + x) ^ 2)) to sqrt(1-x^2) / (x*(1+x))
    substitute sin(u) for x
    rewrite 1-sin(u)^2 to cos(u)^2
    simplify
    rewrite cos(u)^2 to 1-sin(u)^2
    rewrite 1-sin(u)^2 to (1+sin(u))*(1-sin(u))
    simplify
    expand polynomial
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1 / (exp(x/2) + exp(x))
    substitute u for exp(x/2)
    simplify
    partial fraction decomposition
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. exp(arctan(x))/(1+x^2)^(3/2)
    substitute u for arctan(x)
    rewrite tan(u)^2+1 to sec(u)^2
    simplify
    rewrite sec(u) to 1/cos(u)
    simplify
    integrate by parts with u=exp(u), v=sin(u)
    integrate by parts with u=-exp(u), v=cos(u)
    simplify
    solve integral INT u. cos(u)*exp(u)
    replace substitution
    simplify
done

// calculate INT x. 2/(3-x) * sqrt((5-x)/(x-1))
