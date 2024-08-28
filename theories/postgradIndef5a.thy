# 4.5 Trigonometric functions
# Section A

# page 195

calculate INT x. cos(x)^3 for x > 0, x < pi
    integrate by parts with u=cos(x)^2,v=sin(x)
    rewrite sin(x)^2 to 1-cos(x)^2
    expand polynomial
    simplify
    solve integral INT x. cos(x)^3
    apply integral identity
    expand polynomial
    simplify
done

calculate INT x. sin(x)^2
    integrate by parts with u=sin(x), v=-cos(x)
    rewrite cos(x)^2 to 1-sin(x)^2
    simplify
    apply integral identity
    solve integral INT x. sin(x)^2
    expand polynomial
    simplify
done

calculate INT x. cos(3*x)^2
    integrate by parts with u=cos(3*x), v=sin(3*x)/3
    rewrite sin(3*x)^2 to 1-cos(3*x)^2
    simplify
    apply integral identity
    solve integral INT x. cos(3*x)^2
    simplify
done

calculate INT x. sin(x)^2*cos(x)^2
    rewrite sin(x)^2*cos(x)^2 to (sin(x)*cos(x))^2
    rewrite sin(x)*cos(x) to 1/2*sin(2*x)
    simplify
    integrate by parts with u=sin(2*x),v=-cos(2*x)/2
    rewrite cos(2*x)^2 to 1-sin(2*x)^2
    simplify
    apply integral identity
    solve integral (INT x. sin(2*x)^2)/4
    simplify
done

calculate INT x. sin(x)^2*cos(x)^3
    rewrite cos(x)^3 to cos(x)^2 * cos(x)
    rewrite cos(x)^2 to 1-sin(x)^2
    substitute u for sin(x)
    expand polynomial
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. sin(2*x)*cos(3*x)
    rewrite sin(2*x)*cos(3*x) to (1/2)*(sin(5*x)-sin(x))
    apply integral identity
    simplify
done

calculate INT x. cos(4*x)*cos(2*x)
    rewrite cos(4*x)*cos(2*x) to 1/2 * (cos(2*x)+cos(6*x))
    apply integral identity
    simplify
done

calculate INT x. sin(x)^4
    rewrite sin(x)^4 to sin(x)^2^2
    rewrite sin(x)^2 to 1/2*(1-cos(2*x))
    expand polynomial
    rewrite cos(2*x)^2 to 1/2 * (1+cos(4*x))
    apply integral identity
    simplify
done

calculate INT x. tan(x)^4
    rewrite tan(x)^4 to tan(x)^2^2
    rewrite tan(x)^2 to sec(x)^2-1
    expand polynomial
    apply integral identity
    integrate by parts with u=sec(x)^2, v=tan(x)
    substitute u for tan(x)
    apply integral identity
    replace substitution
    simplify
done
