## Examples on substitution

// Source:
// The Calculus Page Problems List by D. A. Kouba
// Problems on u-substitution 
// URL: https://www.math.ucdavis.edu/~kouba/CalcTwoDIRECTORY/usubdirectory/USubstitution.html 

calculate INT x:[0, 1]. (2*x + 5)*(x^2+5*x)^7
    substitute u for x^2 + 5*x
    apply integral identity
    simplify
done

calculate INT x:[0, 1]. (3 - x)^10
    substitute u for 3 - x
    apply integral identity
    simplify
done

calculate INT x:[0, 1]. sqrt(7*x+9)
    substitute u for 7 * x + 9
    apply integral identity
    simplify
done

calculate INT x:[0, 1]. x^3/(1+x^4)^(1/4)
    substitute u for x ^ 4
    simplify
    substitute v for 1 + u
    apply integral identity
    simplify
done

calculate INT x:[0, 1]. exp(5*x+2)
    substitute u for 5*x+2
    apply integral identity
    simplify
done

calculate INT x:[0, pi/6]. 4*cos(3*x)
    substitute u for 3*x
    apply integral identity
    simplify
done

calculate INT x:[1, exp(2*pi)]. sin(log(x)) / x
    substitute u for log(x)
    apply integral identity
    simplify
done

calculate INT x:[1, 2]. (3*x+6)/(x^2 + 4 * x - 3)
    rewrite x^2 + 4*x - 3 to (x+2)^2 - 7
    substitute u for x + 2
    substitute v for u^2 - 7
    apply integral identity
    simplify
done

calculate INT x:[0, 1]. x*3^(x^2 + 1)
    substitute u for x^2+1
    apply integral identity
    simplify
done

calculate INT x:[exp(1), exp(2)]. 3/(x*log(x))
    substitute u for log(x)
    apply integral identity
    simplify
done

calculate INT x:[0, pi/10]. cos(5*x)/exp(sin(5*x))
    substitute u for sin(5*x)
    apply integral identity
    simplify
done

calculate INT x:[0, sqrt(pi)]. x * sin(x^2)
    substitute u for x^2
    apply integral identity
    simplify
done

calculate INT x:[0, 1]. (x + 3) * (x - 1) ^ 5
    substitute u for x - 1
    expand polynomial
    apply integral identity
    simplify
done

calculate INT x:[0, 4].x*sqrt(4-x)
    substitute u for 4 - x
    expand polynomial
    apply integral identity
    simplify
done

calculate INT x:[0, 1]. (x+5)/(2*x + 3)
    rewrite (x+5)/(2*x + 3) to 1/2 + 7 / (2*(2*x+3))
    simplify
    substitute u for 2*x+3
    apply integral identity
    simplify
done

calculate INT x:[0, 1]. (x^2+4)/(x+2)
    rewrite (x^2+4)/(x+2) to x-2 + 8/(x+2)
    simplify
    substitute u for x + 2
    apply integral identity
    simplify
done

calculate INT x:[1, exp(2)]. ((3+log(x))^2 * (2 - log(x)))/(4 * x)
    substitute u for log(x)
    expand polynomial
    apply integral identity
    simplify
done

calculate INT x:[0, 9]. sqrt(4 - sqrt(x))
    substitute u for sqrt(x)
    substitute v for -u + 4
    expand polynomial
    apply integral identity
    simplify
done
