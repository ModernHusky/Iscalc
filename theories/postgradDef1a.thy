# 5.1 Definite integrals
# Section A

# page 209

calculate INT phi:[0, pi/2].sin(phi)*cos(phi)^3
    substitute u for cos(phi)
    apply integral identity
    simplify
done

calculate INT theta:[0, pi]. 1-sin(theta)^3
    apply integral identity
    rewrite sin(theta)^3 to sin(theta)^2*sin(theta)
    rewrite sin(theta)^2 to 1-cos(theta)^2
    substitute u for cos(theta)
    apply integral identity
    simplify
done

calculate INT u:[pi/6, pi/2]. cos(u)^2
    rewrite cos(u)^2 to 1/2*(1+cos(2*u))
    apply integral identity
    simplify
done

calculate INT x:[0, sqrt(2)]. sqrt(2-x^2)
    substitute sqrt(2)*sin(u) for x
    rewrite 2 - (sqrt(2) * sin(u)) ^ 2 to 2 * (1-sin(u)^2)
    rewrite 1-sin(u)^2 to cos(u)^2
    simplify
    rewrite cos(u)^2 to 1/2*(1+cos(2*u))
    apply integral identity
    simplify
done

calculate INT x:[sqrt(2)/2, 1]. sqrt(1-x^2)/x^2
    substitute sin(u) for x
    rewrite 1-sin(u)^2 to cos(u)^2
    simplify
    rewrite cos(u)^2 / sin(u)^2 to (cos(u)/sin(u))^2
    rewrite cos(u)/sin(u) to cot(u)
    rewrite cot(u)^2 to csc(u)^2-1
    apply integral identity
    simplify
done

calculate INT x:[0, a]. x^2*sqrt(a^2-x^2) for a>0
    substitute a*sin(u) for x
    rewrite a^2-(a*sin(u))^2 to a^2*(1-sin(u)^2)
    rewrite 1-sin(u)^2 to cos(u)^2
    simplify
    rewrite cos(u)^2 to 1/2*(1+cos(2*u))
    rewrite sin(u)^2 to 1/2*(1-cos(2*u))
    expand polynomial
    rewrite cos(2*u)^2 to 1/2*(1+cos(4*u))
    apply integral identity
    simplify
done

calculate INT x:[1, sqrt(3)]. 1/(x^2*sqrt(1+x^2))
    substitute tan(u) for x
    rewrite 1+tan(u)^2 to sec(u)^2
    simplify
    rewrite sec(u) to 1/cos(u)
    rewrite tan(u) to sin(u)/cos(u)
    simplify
    substitute v for sin(u)
    apply integral identity
    simplify
done

calculate INT x:[1, 4]. 1/(1+sqrt(x))
    substitute u for sqrt(x)
    partial fraction decomposition
    apply integral identity
    simplify
done

calculate INT x:[1, exp(1)^2]. 1/(x*sqrt(1+log(x)))
    substitute u for sqrt(1+log(x))
    apply integral identity
    simplify
done

calculate INT x:[-2, 0]. (x+2)/(x^2+2*x+2)
    substitute u for x+1
    expand polynomial
    apply integral identity
    substitute v for u^2
    simplify
done

calculate INT theta:[-pi/2, pi/2]. 4*cos(theta)^4
    split region at 0
    substitute u for -theta
    simplify
    rewrite cos(u)^4 to cos(u)^2^2
    rewrite cos(u)^2 to 1/2*(1+cos(2*u))
    expand polynomial
    rewrite cos(2*u)^2 to 1/2*(1+cos(4*u))
    apply integral identity
    simplify
done
