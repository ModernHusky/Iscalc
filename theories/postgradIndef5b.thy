# 4.5 Trigonometric functions
# Section B

# page 196

// found an error in the book's answer
// 1/2*sin(x)^2+log(abs(csc(2*x)-cot(2*x))) + C
calculate INT x. 1/(cos(x)*sin(x)^3) for x > 0, x < pi/2
    rewrite 1 to sin(x)^2+cos(x)^2
    expand polynomial
    simplify
    rewrite cos(x)*sin(x) to 1/2*sin(2*x)
    simplify
    rewrite 1/sin(2*x) to csc(2*x)
    substitute u for 2*x
    apply integral identity
    substitute v for sin(x)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. (3*cos(x)-sin(x)) / (cos(x)+sin(x))
    rewrite 3*cos(x)-sin(x) to (cos(x)+sin(x))+2*(cos(x)-sin(x))
    rewrite (cos(x) + sin(x) + 2 * (cos(x) - sin(x))) / (cos(x) + sin(x)) to 1 + 2 * (cos(x) - sin(x)) / (cos(x) + sin(x))
    simplify
    substitute u for cos(x)+sin(x)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. tan(x)*cos(x)^6/sin(x)^4 for x > 0 , x < pi/2
    rewrite tan(x) to sin(x)/cos(x)
    simplify
    rewrite cos(x)^5 to cos(x)^2^2*cos(x)
    rewrite cos(x)^2 to 1-sin(x)^2
    substitute u for sin(x)
    expand polynomial
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. cos(x)*sin(x)^3 / (1+cos(x)^2) for x > 0, x < pi / 2
    rewrite sin(x)^3 to sin(x)^2*sin(x)
    rewrite sin(x)^2 to 1-cos(x)^2
    substitute u for cos(x)
    substitute v for u^2
    partial fraction decomposition
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. sin(x)*cos(x)^3 / (1+sin(x)^2) for x > 0, x < pi / 2
    rewrite cos(x)^3 to cos(x)^2*cos(x)
    rewrite cos(x)^2 to 1-sin(x)^2
    substitute u for sin(x)
    substitute v for u^2
    partial fraction decomposition
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1/(1+4*cos(x))
    substitute t for tan(x/2)
    rewrite (-(4 * t ^ 2) + 4) / (t ^ 2 + 1) + 1 to (-4*t^2+4)/(t^2+1)+(t^2+1)/(t^2+1) 
    rewrite (-4*t^2+4)/(t^2+1)+(t^2+1)/(t^2+1) to (5-3*t^2)/(t^2+1)
    simplify
    substitute u for sqrt(3/5)*t
    simplify
    partial fraction decomposition
    apply integral identity
    replace substitution
    expand polynomial
done

calculate INT x. 1/(sqrt(sin(x))*cos(x)) for x > 0, x < pi
    substitute u for sqrt(sin(x))
    simplify
    partial fraction decomposition
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1/sqrt(tan(x)^2+2) for x > -pi/2, x < pi/2
    rewrite tan(x) to sin(x)/cos(x)
    rewrite 1 / sqrt((sin(x) / cos(x)) ^ 2 + 2) to cos(x)/(cos(x)*sqrt((sin(x) / cos(x)) ^ 2 + 2))
    rewrite cos(x)*sqrt((sin(x) / cos(x)) ^ 2 + 2) to sqrt(sin(x) ^ 2 + 2*cos(x)^2)
    substitute u for sin(x)
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1/(sin(x)^4+cos(x)^4)
    rewrite 1/(sin(x)^4+cos(x)^4) to (1/cos(x))^4 / ((sin(x)/cos(x))^4+1)
    rewrite 1/cos(x) to sec(x)
    rewrite sin(x)/cos(x) to tan(x)
    rewrite sec(x)^4 to sec(x)^2*sec(x)^2
    rewrite sec(x)^2 to tan(x)^2+1
    substitute u for tan(x)
    rewrite (u ^ 2 + 1) / (u ^ 4 + 1) to (1 / u ^ 2 + 1)/(u^2+1/u^2)
    rewrite u^2+1/u^2 to (u-1/u)^2+2
    substitute v for u-1/u
    substitute w for v/sqrt(2)
    apply integral identity
    replace substitution
    simplify
done