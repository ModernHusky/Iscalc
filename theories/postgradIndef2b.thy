# 4.2 Substitution
# Section B
// page 175

# 4.2 Substitution
# Section B
// page 175

calculate INT x. exp(exp(x)*sin(x))*(sin(x)+cos(x))*exp(x) for cos(x) + sin(x) != 0
    substitute u for exp(x)*sin(x)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. exp(arcsin(sqrt(x)))/sqrt(x-x^2) for x < 1, x > 0
    rewrite sqrt(x-x^2) to sqrt(x) * sqrt(1-x)
    substitute u for arcsin(sqrt(x))
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. exp(tan(1/x))/x^2 * sec(1/x)^2 for sec(1/x) != 0, x != 0
    substitute u for tan(1/x)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. exp(sqrt(1+sin(x)))*cos(x) / sqrt(1+sin(x))
    substitute u for sqrt(1+sin(x))
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. (1+2*x^2)*exp(x^2) / (2-3*x*exp(x^2))
    substitute u for 2-3*x*exp(x^2)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. sin(x)*cos(x)^3 / (1+cos(x)^2) for sin(x) != 0, cos(x) != 0
    substitute u for cos(x)
    substitute v for u^2
    simplify
    partial fraction decomposition
    apply integral identity
    replace substitution
    simplify
done

// page 176

calculate INT x. 1 / (arcsin(sqrt(x))*sqrt(x-x^2)) for x > 0, x < 1
    rewrite sqrt(x-x^2) to sqrt(x)*sqrt(1-x)
    rewrite 1 / (arcsin(sqrt(x)) * (sqrt(x) * sqrt(1 - x))) to 1 / arcsin(sqrt(x)) / (sqrt(x) * sqrt(1 - x))
    substitute u for arcsin(sqrt(x))
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1 / (sqrt(1+x)+(1+x)^(3/2)) for x > -1
    rewrite (1+x)^(3/2) to (1+x) * sqrt(1+x)
    rewrite sqrt(1 + x) + (1 + x) * sqrt(1 + x) to sqrt(1+x)*(1+1+x)
    substitute u for sqrt(1+x)
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1 / (exp(x)*(1+exp(2*x)))
    rewrite exp(2*x) to exp(x)^2
    substitute u for exp(x)
    substitute tan(v) for u
    rewrite tan(v)^2 + 1 to sec(v)^2
    simplify
    rewrite tan(v) to sin(v)/cos(v)
    simplify
    rewrite cos(v)^2 to 1 - sin(v)^2
    expand polynomial
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1 / ((2-x)*sqrt(1-x)) for x < 1
    rewrite 2-x to 1 + (1-x)
    substitute u for sqrt(1-x)
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1/(a^2*sin(x)^2+b^2*cos(x)^2) for a != 0, b != 0
    rewrite 1/(a^2*sin(x)^2+b^2*cos(x)^2) to (1/cos(x)^2) / (a^2*(sin(x)/cos(x))^2 + b^2)
    rewrite 1/cos(x)^2 to sec(x)^2
    rewrite sin(x)/cos(x) to tan(x) 
    substitute u for tan(x)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1/(2+cos(x))
    rewrite cos(x) to cos(2*(x/2))
    rewrite cos(2*(x/2)) to cos(x/2)^2 - sin(x/2)^2
    rewrite 2 + (cos(x / 2) ^ 2 - sin(x / 2) ^ 2) to  1 - sin(x/2)^2 + 1 + cos(x/2)^2
    rewrite  1 - sin(x/2)^2 to cos(x/2)^2
    simplify
    rewrite 1 to sin(x/2)^2 + cos(x/2)^2 (at 2)
    simplify
    rewrite 1 / (3 * cos(x / 2) ^ 2 + sin(x / 2) ^ 2) to 1/cos(x/2)^2 / (3 + (sin(x/2)/cos(x/2))^2)
    rewrite 1/cos(x/2)^2 to sec(x/2)^2
    rewrite sin(x/2)/cos(x/2) to tan(x/2)
    substitute u for tan(x/2)
    rewrite u ^ 2 / 2 to (1/2) * u^2
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 2^x * 3^x / (9^x+4^x)
    rewrite 2^x * 3^x to 6^x
    rewrite 6 ^ x / (9 ^ x + 4 ^ x) to (6^x/4^x)/(9^x/4^x + 1)
    rewrite 6^x/4^x to (3/2)^x
    rewrite 9^x/4^x to (3/2)^2^x
    rewrite (3/2)^2^x to (3/2)^x^2
    substitute u for (3/2)^x
    apply integral identity
    replace substitution
    simplify
done

// page 177

calculate INT x. sqrt(x/(1-x^3)) for  x != 0
    rewrite x^3 to x^(3/2)^2
    substitute u for x^(3/2)
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. exp(x/2) / (16-exp(x))     
    rewrite exp(x) to exp(x/2)^2
    substitute u for exp(x/2)
    substitute v for u/4
    rewrite 4 / (-(8 * v ^ 2) + 8) to 1/2 * (1 / (1-v^2))
    apply integral identity
    partial fraction decomposition
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. cos(x) / sqrt(2+cos(2*x))
    rewrite cos(2*x) to cos(x)^2 - sin(x)^2
    rewrite cos(x)^2 to 1-sin(x)^2
    simplify
    substitute u for sin(x)
    apply integral identity
    substitute v for sqrt(2)*u/sqrt(3)
    rewrite sqrt(-(3 * v ^ 2) + 3) to sqrt(3) * sqrt(1-v^2)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1/sqrt((x-a)*(x-b)) for a < b
    rewrite (x-a)*(x-b) to (x-(a+b)/2)^2 - ((a-b)/2)^2
    substitute u for (x-(a+b)/2)/((a-b)/2)
    simplify
    rewrite sqrt(u ^ 2 * (a - b) ^ 2 / 4 - (a - b) ^ 2 / 4) to (a-b)/2 * sqrt(u^2-1)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. (1+x) / (x*(1+x*exp(x))) for x != -1, x != 0
    rewrite (1+x) / (x*(1+x*exp(x))) to ((1+x)*exp(x))/((x*exp(x)) * (1+x*exp(x)))
    substitute u for x*exp(x)
    partial fraction decomposition
    apply integral identity
    replace substitution
    simplify
done

// page 178

calculate INT x. exp(x)*(1+x)/(1-x*exp(x)) for x != -1
    substitute u for x*exp(x)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. exp(x)*(x-1)/(x-exp(x))^2 for x != 0
    rewrite exp(x)*(x-1)/(x-exp(x))^2 to (exp(x)*(x-1)/x^2) / (1 - exp(x)/x)^2
    substitute u for exp(x)/x
    substitute v for 1-u
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. (x + sin(x)*cos(x)) / (cos(x) - x*sin(x))^2
    rewrite (x + sin(x)*cos(x)) / (cos(x) - x*sin(x))^2 to (x * (1/cos(x)^2) + sin(x)/cos(x)) / (1 - x * (sin(x)/cos(x)))^2
    rewrite 1/cos(x)^2 to sec(x)^2
    rewrite sin(x)/cos(x) to tan(x)
    rewrite sin(x)/cos(x) to tan(x)
    substitute u for x*tan(x)
    substitute v for -u+1
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. (1-log(x))/(x-log(x))^2
    substitute u for log(x)
    rewrite exp(u) * (-u + 1) / (exp(u) - u) ^ 2 to -exp(u) * (u-1)/u^2 / (exp(u)/u - 1) ^ 2
    substitute v for exp(u) / u
    substitute w for v-1
    apply integral identity
    replace substitution
    simplify
done
