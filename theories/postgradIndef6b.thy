# 4.6 Partial fraction
# Section B

# page 200

calculate INT x. (x^2+1)/((x-1)*(x+1)^2)
    partial fraction decomposition
    apply integral identity
    substitute u for x+1
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. x/((x-1)*(x^2+1))
    partial fraction decomposition
    expand polynomial
    apply integral identity
    substitute u for x^2
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. x^3/(x+3)
    partial fraction decomposition
    apply integral identity
    simplify
done

calculate INT x. 1/(x*(x^6+3))
    partial fraction decomposition
    apply integral identity
    substitute u for x^6
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. x^2/(1+x^2)^2
    partial fraction decomposition
    apply integral identity
    substitute tan(u) for x
    rewrite tan(u)^2+1 to sec(u)^2
    simplify
    rewrite sec(u) to 1/cos(u)
    simplify
    rewrite cos(u)^2 to 1/2*(1+cos(2*u))
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1/(x*(1+x^4))
    partial fraction decomposition
    simplify
    apply integral identity
    substitute u for x^4
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1/(x^4*(1+x^2))
    partial fraction decomposition
    apply integral identity
    simplify
done

calculate INT x. (x^2+1)/(1+x^4) for x != 0
    rewrite (x^2+1)/(1+x^4) to (1+1/x^2)/(1/x^2+x^2)
    rewrite (1/x^2+x^2) to (x-1/x)^2+2
    substitute u for x-1/x
    substitute v for u/sqrt(2)
    simplify
    apply integral identity
    replace substitution
    simplify
done