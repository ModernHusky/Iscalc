# 4.4 Integration by parts
# Section A

calculate INT x. x*sin(x)
    integrate by parts with u=-x, v=cos(x)
    apply integral identity
    simplify
done

calculate INT x. log(x) for x > 0
    integrate by parts with u=log(x), v=x
    apply integral identity
    simplify
done

calculate INT x. arcsin(x) for x > -1, x < 1
    integrate by parts with u=arcsin(x), v=x
    substitute sin(u) for x
    rewrite -(sin(u)^2)+1 to cos(u)^2
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. x*exp(-x)
    integrate by parts with u=x, v=-exp(-x)
    apply integral identity
    simplify
done

calculate INT x. x*log(x)^2 for x > 0
    integrate by parts with u=log(x)^2, v=1/2*x^2
    integrate by parts with u=log(x), v=1/2*x^2
    apply integral identity
    simplify
done

calculate INT x. exp(-x) * cos(x)
    integrate by parts with u = cos(x),  v=-exp(-x)
    simplify
    integrate by parts with u = sin(x), v = -exp(-x)
    simplify
    solve integral INT x. exp(-x) * cos(x)
    simplify
done

calculate INT x. x*cos(x/2)
    integrate by parts with u=x, v=sin(x/2)*2
    substitute u for x/2
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. x^2*arctan(x)
    integrate by parts with u=arctan(x), v=x^3/3
    simplify
    rewrite 3*x^2+3 to 3*(x^2+1)
    substitute tan(u) for x
    rewrite tan(u)^2 + 1 to sec(u)^2
    simplify
    rewrite tan(u)^3 to tan(u)^2 * tan(u)
    rewrite tan(u)^2 to sec(u)^2 - 1
    expand polynomial
    apply integral identity
    substitute v for tan(u)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. x*tan(x)^2
    rewrite tan(x)^2 to sec(x)^2-1
    expand polynomial
    apply integral identity
    integrate by parts with u=x,v=tan(x)
    apply integral identity
    simplify
done

calculate INT x. x^2*cos(x)
    integrate by parts with u=x^2, v=sin(x)
    integrate by parts with u=2*x, v=-cos(x)
    apply integral identity
    simplify
done

calculate INT x. log(x)^2 for x>0
    integrate by parts with u=log(x)^2, v=x
    integrate by parts with u=2*log(x), v=x
    apply integral identity
    simplify
done

calculate INT x. x^2*cos(x/2)^2
    integrate by parts with u=cos(x/2)^2, v=x^3/3
    rewrite cos(x/2)*sin(x/2) to sin(x)/2
    simplify
    integrate by parts with u=x^3, v=-cos(x)
    simplify
    integrate by parts with u=x^2, v=sin(x)
    simplify
    integrate by parts with u=x, v=-cos(x)
    apply integral identity
    simplify
done

calculate INT x. x*log(x-1)
    integrate by parts with u=log(x-1),v=x^2/2
    simplify
    partial fraction decomposition
    apply integral identity
    simplify
done

calculate INT x. cos(sqrt(x))^2
    substitute u for sqrt(x)
    integrate by parts with u=cos(u)^2, v=u^2
    rewrite 2 * u ^ 2 * cos(u) * sin(u) to 2 * u^2 * (cos(u)*sin(u))
    rewrite cos(u)*sin(u) to sin(2*u)/2
    simplify
    integrate by parts with u=u^2,v=-cos(2*u)/2
    integrate by parts with u=-u, v=sin(2*u)/2
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. log(x+2)/(x+1)^2
    integrate by parts with u=log(x+2),v=-1/(x+1)
    partial fraction decomposition
    apply integral identity
    simplify
done