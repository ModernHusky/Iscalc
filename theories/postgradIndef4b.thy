# 4.4 Integration by parts
# Section B

# page 186 

calculate INT x. arctan(sqrt(x))
    substitute u for sqrt(x)
    integrate by parts with u=arctan(u), v=u^2
    partial fraction decomposition
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. arctan(x)/x^2
    integrate by parts with u=arctan(x), v=-1/x
    partial fraction decomposition
    apply integral identity
    substitute u for x^2
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. arctan(exp(x))/exp(x)
    substitute u for exp(x)
    integrate by parts with u=arctan(u), v=-1/u
    partial fraction decomposition
    apply integral identity
    substitute v for u^2
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. arctan(x) / (x^2*(1+x^2)) for x != 0
    substitute tan(u) for x
    simplify
    rewrite tan(u)^2 + 1 to sec(u)^2
    simplify
    rewrite tan(u) to sin(u)/cos(u)
    simplify
    rewrite cos(u)^2 to 1-sin(u)^2
    expand polynomial
    simplify
    apply integral identity
    integrate by parts with u=-u, v=cot(u)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. arctan(sqrt(x)) / (sqrt(x)+sqrt(x^3)) for x > 0
    substitute u for sqrt(x)
    simplify
    rewrite u^3+u to u*(u^2+1)
    simplify
    substitute v for arctan(u)
    apply integral identity
    replace substitution
    simplify
done 
    
// page 187

calculate INT x. x*exp(x) / sqrt(exp(x)-1) for x > 0
    integrate by parts with u=2*x, v=sqrt(exp(x)-1)
    substitute u for sqrt(exp(x)-1)
    simplify
    partial fraction decomposition
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. x*exp(x) / (1+x)^2 for x != -1
    integrate by parts with u=-x*exp(x), v=1/(1+x)
    rewrite -(x * exp(x)) - exp(x) to -exp(x) * (x+1)
    simplify
    apply integral identity
    simplify
done

calculate INT x. x*exp(-x) / (1-x)^2 for x != 1
    integrate by parts with u=x*exp(-x),v=1/(1-x)
    rewrite -(x * exp(-x)) + exp(-x) to exp(-x) * (-x + 1)
    simplify
    apply integral identity
    simplify
done

calculate INT x. x^2*exp(x) / (x+2)^2 for x != -2
    integrate by parts with u = -x^2*exp(x), v=1/(x+2)
    rewrite -(x ^ 2 * exp(x)) - 2 * x * exp(x) to -x*exp(x) * (x+2)
    simplify
    integrate by parts with u=x, v=exp(x)
    apply integral identity
    simplify
done

calculate INT x. x * exp(x) / (1+exp(x))^2
    integrate by parts with u=-x, v=1/(exp(x) + 1)
    simplify
    substitute u for exp(x)
    partial fraction decomposition
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. x*exp(x) / sqrt(exp(x)-2) for exp(x)-2>0
    integrate by parts with u=2*x, v=sqrt(exp(x)-2)
    substitute u for sqrt(exp(x)-2)
    simplify
    partial fraction decomposition
    apply integral identity
    substitute v for u/sqrt(2)
    apply integral identity 
    replace substitution
    simplify
done

// page 188

calculate INT x. sqrt(a^2+x^2) for a > 0
    integrate by parts with u=sqrt(a^2+x^2),v=x
    rewrite x^2 to x^2 + a^2 - a^2 (at 2)
    rewrite x^2 + a^2 to sqrt(x^2+a^2)^2
    rewrite (sqrt(x ^ 2 + a ^ 2) ^ 2 - a ^ 2) / sqrt(a ^ 2 + x ^ 2) to sqrt(x^2+a^2) - a^2 / sqrt(a^2+x^2)
    apply integral identity
    substitute a*tan(u) for x (at 2)
    rewrite a ^ 2 + (a * tan(u)) ^ 2 to a^2 * (tan(u)^2+1)
    rewrite tan(u)^2+1 to sec(u)^2
    simplify
    apply integral identity
    solve integral INT x. sqrt(a^2+x^2)
    replace substitution
    simplify
done

calculate INT x. sqrt(a^2-x^2) for a > 0
    integrate by parts with u=sqrt(a^2-x^2),v=x
    simplify
    rewrite x^2 to -(sqrt(a^2-x^2)^2 - a^2)
    rewrite -(sqrt(a ^ 2 - x ^ 2) ^ 2 - a ^ 2) / sqrt(a ^ 2 - x ^ 2) to -sqrt(a^2-x^2)+a^2/sqrt(a^2-x^2)
    apply integral identity
    solve integral INT x. sqrt(a^2-x^2)
    simplify
done

calculate INT x. sec(x)^3
    integrate by parts with u=sec(x), v=tan(x)
    rewrite tan(x)^2 to sec(x)^2-1
    expand polynomial
    apply integral identity
    solve integral INT x. sec(x)^3
    simplify
done

calculate INT x. csc(x)^3
    integrate by parts with u=-csc(x), v=cot(x)
    rewrite cot(x)^2 to csc(x)^2-1
    expand polynomial
    apply integral identity
    solve integral INT x. csc(x)^3
    simplify
done

calculate INT x. exp(a*x)*sin(b*x) for a > 0, b > 0
    integrate by parts with u=sin(b*x), v=exp(a*x)/a
    simplify
    integrate by parts with u=cos(b*x), v=exp(a*x)/a
    simplify
    solve integral INT x. exp(a*x)*sin(b*x)
    expand polynomial
    expand polynomial
    simplify
done

calculate INT x. exp(a*x)*cos(b*x) for a > 0, b > 0
    integrate by parts with u=cos(b*x), v=exp(a*x)/a
    simplify
    integrate by parts with u=sin(b*x), v=exp(a*x)/a
    simplify
    solve integral INT x. exp(a*x)*cos(b*x)
    expand polynomial
    expand polynomial
    simplify
done

calculate INT x. exp(x) * ((1-x)/(1+x^2))^2
    rewrite ((1-x)/(1+x^2))^2 to 1/(1+x^2) - 2*x/(1+x^2)^2
    expand polynomial
    simplify
    integrate by parts with u=-exp(x)/2, v=1/(1+x^2)
    simplify
    rewrite exp(x) / (2 * x ^ 2 + 2) to exp(x)/(x^2+1) * (1/2)
    simplify
done

// missing skolem variable
calculate INT x. exp(-x) * (1+sin(x)) / (1-cos(x)) for x > 0, x < pi/2
    rewrite sin(x) to 2*sin(x/2)*cos(x/2)
    rewrite cos(x) to 1-2*sin(x/2)^2
    rewrite 1 - (1 - 2 * sin(x / 2) ^ 2) to 2*sin(x/2)^2
    rewrite exp(-x)*(1+2*sin(x/2)*cos(x/2))/(2*sin(x/2)^2) to exp(-x)/(2*sin(x/2)^2)+exp(-x)*cos(x/2)/sin(x/2)
    simplify
    rewrite exp(-x) / sin(x / 2) * cos(x / 2) to exp(-x) * (cos(x/2)/sin(x/2))
    rewrite cos(x/2)/sin(x/2) to cot(x/2)
    rewrite exp(-x)/sin(x/2)^2 to exp(-x) * (1/sin(x/2))^2
    rewrite 1/sin(x/2) to csc(x/2)
    integrate by parts with u=exp(-x), v=-cot(x/2)*2 (at 2)
    simplify
done

// page 189

// missing skolem variable
calculate INT x. exp(sin(x))*(x*cos(x)^3-sin(x))/(cos(x)^2) for x > -pi/2, x < pi / 2
    expand polynomial
    simplify
    integrate by parts with u=x, v=exp(sin(x)) (at 2)
    simplify
    integrate by parts with u=exp(sin(x)), v=1/cos(x)
    simplify
done

calculate INT x. sqrt(1-x^2)*arcsin(x) for x > -1, x < 1
    substitute sin(u) for x
    simplify
    rewrite -(sin(u)^2)+1 to cos(u)^2
    simplify
    rewrite cos(u)^2 to (cos(2*u)+1)/2
    expand polynomial
    apply integral identity
    integrate by parts with u=u,v=sin(2*u)/2
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. arcsin(x)^2 for x > -1, x < 1
    substitute sin(u) for x
    simplify
    integrate by parts with u=u^2, v=sin(u)
    integrate by parts with u=2*u, v=-cos(u)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. arcsin(x)*x for x > -1, x < 1
    substitute sin(u) for x
    simplify
    rewrite u*cos(u)*sin(u) to u*(sin(u)*cos(u))
    rewrite sin(u)*cos(u) to sin(2*u)/2
    simplify
    integrate by parts with u=-u, v=cos(2*u)/2
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. arcsin(sqrt(x))/sqrt(x) for x > 0, x < 1
    substitute u for sqrt(x)
    substitute sin(v) for u
    simplify
    integrate by parts with u=v,v=sin(v)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. arcsin(exp(x)) / exp(x)
    substitute u for exp(x)
    substitute sin(v) for u
    simplify
    integrate by parts with u=v, v=-1/sin(v)
    rewrite 1/sin(v) to csc(v)
    apply integral identity
    replace substitution
    simplify
done

// page 190

calculate INT x. arccos(x) / sqrt((1-x^2)^3) for x > -1, x < 1
    substitute cos(u) for x
    rewrite 1-cos(u)^2 to sin(u)^2
    simplify
    rewrite u/sin(u)^2 to u * (1/sin(u))^2
    rewrite 1/sin(u) to csc(u)
    integrate by parts with u=-u, v=cot(u)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. x*arccos(x)/sqrt(1-x^2) for x > -1, x < 1
    substitute cos(u) for x
    rewrite 1-cos(u)^2 to sin(u)^2
    simplify
    integrate by parts with u=u, v=sin(u)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. arcsin(x) / x^2 * (1+x^2) / sqrt(1-x^2) for x > 0, x < 1
    substitute sin(u) for x
    rewrite 1-sin(u)^2 to cos(u)^2
    simplify
    expand polynomial
    rewrite u/sin(u)^2 to u * (1/sin(u))^2
    rewrite 1/sin(u) to csc(u)
    apply integral identity
    integrate by parts with u=u, v=-cot(u)
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. arcsin(x)-arccos(x) for x > -1, x < 1
    simplify
    integrate by parts with u=arccos(x), v=x
    substitute sin(u) for x
    rewrite -(sin(u)^2)+1 to cos(u)^2
    simplify
    apply integral identity
    integrate by parts with u=arcsin(x), v=x
    substitute sin(v) for x
    rewrite -(sin(v)^2)+1 to cos(v)^2
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. arcsin(x)*arccos(x) for x>-1, x<1
    integrate by parts with u=arcsin(x)*arccos(x), v=x
    expand polynomial
    simplify
    substitute cos(u) for x
    simplify
    rewrite -(cos(u)^2)+1 to sin(u)^2
    simplify
    integrate by parts with u=u, v=sin(u) (at 2)
    apply integral identity
    substitute sin(v) for x
    rewrite -(sin(v)^2)+1 to cos(v)^2
    simplify
    integrate by parts with u=-v,v=cos(v)
    apply integral identity
    replace substitution
    simplify
done

// page 191

calculate INT x. log(1+x)/x^2
    integrate by parts with u=-log(1+x), v=1/x
    partial fraction decomposition
    apply integral identity
    simplify
done

calculate INT x. log(1+x^2)/x^2 for x!=0
    integrate by parts with u=log(1+x^2), v=-1/x
    apply integral identity
    simplify
done

calculate INT x. log(x)^2/x^3 for x > 0
    integrate by parts with u=log(x)^2, v=1/x^2/-2
    integrate by parts with u=log(x), v=1/x^2/2
    apply integral identity
    simplify
done

calculate INT x. x*log(x)/(1+x^2)^2 for x>0
    integrate by parts with u=log(x), v=1/(1+x^2)/-2
    partial fraction decomposition
    apply integral identity
    substitute u for x^2
    apply integral identity
    replace substitution
    simplify
done

// missing skolem variable
calculate INT x. log(log(x))+1/log(x)
    substitute u for log(x)
    expand polynomial
    simplify
    integrate by parts with u=log(u), v=exp(u)
    simplify
    replace substitution
    simplify
done

calculate INT x. log(x+sqrt(x^2+1)) for x > -1
    integrate by parts with u=log(x+sqrt(x^2+1)), v=x
    rewrite x/sqrt(x^2+1)+1 to (x+sqrt(x^2+1))/sqrt(x^2+1)
    simplify
    substitute tan(u) for x
    rewrite tan(u)^2+1 to sec(u)^2
    simplify
    apply integral identity
    replace substitution
    simplify
done

// page 192

calculate INT x. tan(x)^4
    rewrite tan(x)^4 to tan(x)^2^2
    rewrite tan(x)^2 to sec(x)^2-1
    expand polynomial
    apply integral identity
    integrate by parts with u=sec(x)^2,v=tan(x)
    rewrite sec(x)^2 to tan(x)^2+1 (at 2)
    expand polynomial
    apply integral identity
    rewrite tan(x)^2 to sec(x)^2-1
    apply integral identity
    solve integral INT x. tan(x)^4
    expand polynomial
    simplify
done

calculate INT x. 1/sin(x)^3 for x>0, x<pi
    rewrite 1/sin(x)^3 to csc(x)^3
    integrate by parts with u=csc(x), v=-cot(x)
    rewrite cot(x)^2 to csc(x)^2-1
    expand polynomial
    apply integral identity
    solve integral INT x. csc(x)^3
    expand polynomial
    simplify
done

calculate INT x. arcsin(x)^3 for x>-1, x<1
    substitute sin(u) for x
    simplify
    integrate by parts with u=u^3,v=sin(u)
    integrate by parts with u=3*u^2,v=-cos(u)
    simplify
    integrate by parts with u=u,v=sin(u)
    apply integral identity
    replace substitution
    simplify
done