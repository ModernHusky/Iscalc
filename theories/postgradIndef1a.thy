# 4.1 Application of basic formulas
# Section A
// page 166

calculate INT x. 1 / sqrt(x-x^2)
    rewrite x - x^2 to 1/4 - (x-1/2)^2
    substitute u for arcsin(x-1/2)
    substitute v for sin(u)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1/x^2
    apply integral identity
    simplify
done

calculate INT x. arcsin(2*x-1) for x > 0, x < 1
    substitute u for 2*x-1
    simplify
    integrate by parts with u=arcsin(u), v=u
    substitute v for arcsin(u)
    rewrite -sin(v)^2 + 1 to 1-sin(v)^2
    rewrite 1-sin(v)^2 to cos(v)^2
    simplify
    expand polynomial
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. arccos(1-2*x) for x > 0, x < 1
    substitute u for 1-2*x
    simplify
    integrate by parts with u=arccos(u), v=u
    substitute v for arcsin(u)
    rewrite -sin(v)^2 + 1 to 1-sin(v)^2
    rewrite 1-sin(v)^2 to cos(v)^2
    simplify
    expand polynomial
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. x * sqrt(x)
    rewrite x * sqrt(x) to x^(3/2)
    apply integral identity
    simplify
done

calculate INT x. 1 / sqrt(x)
    rewrite 1 / sqrt(x) to x^(-1/2)
    apply integral identity
    simplify
done

calculate INT x. x^2 * x^(1/3)
    rewrite x^2 * x^(1/3) to x^(7/3)
    apply integral identity
    simplify
done

calculate INT x. 1/(x^2*sqrt(x))
    rewrite 1/(x^2*sqrt(x)) to x^(-5/2)
    apply integral identity
    simplify
done

calculate INT x. (x^2)^(1/3)
    rewrite (x^2)^(1/3) to x^(2/3)
    apply integral identity
    simplify
done

calculate INT x. x^2-3*x+2
    apply integral identity
    simplify
done

calculate INT x. (x^2+1)^2
    expand polynomial
    apply integral identity
    simplify
done

// page 167

calculate INT x. (1-x)^2 / sqrt(x)
    expand polynomial
    simplify
    apply integral identity
    simplify
done

calculate INT x. 2*exp(x)+3/x
    apply integral identity
    simplify
done

calculate INT x. (3/(1+x^2)-2/sqrt(1-x^2))
    apply integral identity
    simplify
done

calculate INT x. exp(x)*(1-exp(-x)/sqrt(x))
    expand polynomial
    simplify
    apply integral identity
done

calculate INT x. 3^x * exp(x)
    rewrite exp(x) to exp(1)^x
    rewrite 3^x * exp(1)^x to (3*exp(1)) ^ x
    apply integral identity
    simplify
done

calculate INT x. (2*3^x - 5*2^x) / 3^x
    expand polynomial
    simplify
    rewrite 3^-x to (1/3)^x
    rewrite 2^x * (1/3)^x to (2/3)^x
    apply integral identity
    simplify
done

calculate INT x. sin(x/2)^2
    integrate by parts with u=-2*sin(x/2),v=cos(x/2)
    simplify
    rewrite cos(x/2)^2 to 1 - sin(x/2)^2
    simplify
    solve integral INT x. sin(x/2)^2
    apply integral identity
    simplify
done

calculate INT x. cos(x/2)^2
    integrate by parts with u=2*cos(x/2),v=sin(x/2)
    simplify
    rewrite sin(x/2)^2 to 1 - cos(x/2)^2
    simplify
    apply integral identity
    solve integral INT x. cos(x/2)^2
    simplify
done

calculate INT x. sec(x)*(sec(x)-tan(x))
    expand polynomial
    apply integral identity
    simplify
done

calculate INT x. tan(x)^2
    rewrite tan(x)^2 to sec(x)^2 - 1
    apply integral identity
    simplify
done

calculate INT x. cot(x)^2
    rewrite cot(x)^2 to csc(x)^2 - 1
    apply integral identity
    simplify
done

calculate INT x. cos(2*x)/(cos(x)-sin(x)) for x > 0, x < pi/4
    rewrite cos(2*x) to cos(x)^2-sin(x)^2
    rewrite cos(x)^2-sin(x)^2 to (cos(x)+sin(x))*(cos(x)-sin(x))
    simplify
    apply integral identity
    simplify
done

calculate INT x. cos(2*x)/(cos(x)^2*sin(x)^2) for x > 0, x < pi/2
    rewrite cos(2*x) to cos(x)^2-sin(x)^2
    expand polynomial
    simplify
    apply integral identity
    simplify
done
