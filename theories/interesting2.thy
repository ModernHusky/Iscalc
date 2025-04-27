imports standard

# Chapter 2

## Chapter 2, Section 1, Six 'Easy' Warm-Ups

prove (INT x:[1,oo]. 1 / ((x+a)*sqrt(x-1))) = pi / sqrt(a+1) for a: real, a > -1
lhs:
    substitute t for sqrt(x - 1)
    simplify
    substitute y for t / sqrt(a + 1)
    rewrite y ^ 2 * (a + 1) + a + 1 to (a + 1) * (y^2 + 1)
    apply integral identity
    simplify
done

prove (INT x:[0, oo]. log(1 + a^2 / x^2)) = a * pi for a: real, a > 0
lhs:
    integrate by parts with u = log(1+a^2/x^2), v = x
    simplify
    rewrite x^2 * (a^2 / x^2 + 1) to a^2 + x^2
    apply integral identity
    simplify
done

prove (INT x:[0, oo]. log(x) / (x^2+b^2)) = pi * log(b) / (2*b) for b: real, b > 0
lhs:
    substitute 1/t for x
    rewrite log(1/t) to -log(t)
    rewrite -log(t) / ((1/t)^2 + b^2) * -(1/t^2) to log(t) / (1 + b^2*t^2)
    substitute s/b for t
    rewrite log(s/b) to log(s) - log(b)
    expand polynomial
    apply integral identity
    simplify
done

prove (INT x:[1,oo]. log(x) / (x+1)^2) = log(2)
subgoal 1: (INT x:[0,oo]. 1 / (1 + exp(a*x))) = log(2) / a for a: real, a > 0
lhs:
    substitute u for exp(a * x)
    simplify
    rewrite 1 / (u * (u+1)) to 1/u - 1/(u+1)
    simplify
    substitute y for u + 1
    apply integral identity
    simplify
done
subgoal 2: (INT x:[1,oo]. log(x) / (a ^ 2 * (x + 1) ^ 2)) = log(2) / a^2 for a: real, a > 0
from 1:
    differentiate both sides at a
    simplify
    substitute y for exp(a * x)
    solve equation for INT y:[1,oo]. log(y) / (a ^ 2 * (y + 1) ^ 2)
done
lhs:
    rewrite (x+1) ^ 2 to 1^2 * (x+1)^2
    apply 2 on INT x:[1,oo]. log(x) / (1 ^ 2 * (x + 1) ^ 2)
    simplify
done

prove (INT x:[sqrt(2),oo]. 1 / (x + x ^ sqrt(2))) = (1 + sqrt(2)) * log(1 + 2 ^ (1/2 * (1 - sqrt(2))))
lhs:
    rewrite 1 / (x + x ^ sqrt(2)) to x ^ -sqrt(2) / (x ^ (-sqrt(2) + 1) + 1)
    substitute u for log(x ^ (1 - sqrt(2)) + 1)
    apply integral identity
    simplify
    rewrite -sqrt(2) + 1 to -1 / (1 + sqrt(2)) (at 2)
    simplify
    rewrite sqrt(2) to 2 ^ (1/2) (at 2)
    rewrite 2 ^ (1/2) ^ (-sqrt(2) + 1) to 2 ^ (1/2 * (-sqrt(2) + 1))
    rewrite (sqrt(2) + 1) * log(2 ^ (1/2 * (-sqrt(2) + 1)) + 1) to (1 + sqrt(2)) * log(1 + 2 ^ (1/2 * (1 - sqrt(2))))
done

prove (INT x:[-oo,oo]. 1 / cosh(x)) = pi
lhs:
    expand definition for cosh (all)
    substitute t for exp(x)
    rewrite t * (1 / t + t) to 1 + t ^ 2
    apply integral identity
    simplify
done

## Chapter 2, Section 2, A New Trick

calculate INT x:[0,pi / 2]. sqrt(sin(x)) / (sqrt(sin(x)) + sqrt(cos(x)))
    substitute y for pi / 2 - x
    rewrite sqrt(cos(y)) / (sqrt(cos(y)) + sqrt(sin(y))) to 1 - sqrt(sin(y)) / (sqrt(cos(y)) + sqrt(sin(y)))
    apply integral identity
    solve integral INT x:[0,pi / 2]. sqrt(sin(x)) / (sqrt(sin(x)) + sqrt(cos(x)))
done

calculate INT x:[0,pi]. x * sin(x) / (1 + cos(x) ^ 2)
    substitute y for pi - x
    expand polynomial
    simplify
    solve integral INT x:[0,pi]. x * sin(x) / (1 + cos(x) ^ 2)
    substitute u for cos(y)
    rewrite -(1 / (-(u ^ 2) - 1)) to 1/(u^2+1)
    apply integral identity
    simplify
done

prove (INT x:[0,pi / 2]. sin(x) ^ 2 / (sin(x) + cos(x))) = sqrt(2) / 4 * log(3 + 2 * sqrt(2))
subgoal 1: (INT x:[0,pi / 2]. sin(x) ^ 2 / (sin(x) + cos(x))) = (INT x:[0,pi / 2]. cos(x) ^ 2 / (sin(x) + cos(x)))
lhs:
    substitute y for pi / 2 - x
done
subgoal 2: (INT x:[0,pi / 2]. sin(x) ^ 2 / (sin(x) + cos(x))) = 1/2 * (INT x:[0,pi / 2]. 1 / (sin(x) + cos(x)))
rhs:
    simplify
    rewrite 1 to sin(x) ^ 2 + cos(x) ^ 2
    expand polynomial
    simplify
    rewrite cos(x) + sin(x) to sin(x) + cos(x) (at 1)
    apply 1 on INT x:[0,pi / 2]. cos(x) ^ 2 / (sin(x) + cos(x))
    simplify
done
lhs:
    apply 2 on INT x:[0,pi / 2]. sin(x) ^ 2 / (sin(x) + cos(x))
    substitute z for tan(x / 2)
    simplify
    rewrite (-(z ^ 2) + 1) / (z ^ 2 + 1) + 2 * z / (z ^ 2 + 1) to (2 - (z - 1) ^ 2) / (z ^ 2 + 1)
    rewrite (z ^ 2 + 1) * ((2 - (z - 1) ^ 2) / (z ^ 2 + 1)) to 2 - (z - 1) ^ 2
    rewrite 2 - (z - 1) ^ 2 to (sqrt(2) + (z - 1)) * (sqrt(2) - (z - 1))
    rewrite 1 / ((sqrt(2) + (z - 1)) * (sqrt(2) - (z - 1))) to sqrt(2) / 4 * (1 / (sqrt(2) + (z - 1)) + 1 / (sqrt(2) - (z - 1)))
    simplify
    substitute u for sqrt(2) + 1 - z (at 1)
    substitute u for sqrt(2) - 1 + z (at 2)
    apply integral identity
    simplify
    rewrite sqrt(2) * (log(sqrt(2) + 1) - log(sqrt(2) - 1)) / 4 to 1/4 * sqrt(2) * (log(sqrt(2) + 1) - log(sqrt(2) - 1))
    rewrite log(sqrt(2) + 1) - log(sqrt(2) - 1) to log((sqrt(2) + 1) / (sqrt(2) - 1))
    rewrite (sqrt(2) + 1) / (sqrt(2) - 1) to 3 + 2 * sqrt(2)
    simplify
done

prove (INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)) = pi / 8 * log(2)
subgoal 1: (INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)) = (INT x:[0,pi / 4]. log(tan(x) + 1))
lhs:
    substitute tan(u) for x
    rewrite sec(u) ^ 2 to tan(u) ^ 2 + 1
    simplify
done
subgoal 2: (INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)) = pi / 4 * log(2) - (INT x:[0,1]. log(x + 1) / (x ^ 2 + 1))
lhs:
    apply 1 on INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)
    substitute pi / 4 - y for x
    simplify
    rewrite tan(pi / 4 - y) to (tan(pi / 4) - tan(y)) / (1 + tan(pi / 4) * tan(y))
    simplify
    rewrite (-tan(y) + 1) / (tan(y) + 1) + 1 to 2 / (1 + tan(y))
    rewrite log(2 / (1 + tan(y))) to log(2) - log(1 + tan(y))
    apply integral identity
    simplify
    apply 1 on INT x:[0,pi / 4]. log(tan(x) + 1)
done
from 2:
    solve equation for INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)
done

prove (INT t:[0,a]. log(t + a) / (t ^ 2 + a ^ 2)) = pi / (8 * a) * log(2 * a ^ 2) for a: real, a > 0
subgoal 1: (INT x:[0,1]. log(x + 1) / (x ^ 2 + 1)) = a * (INT t:[0,a]. log(t + a) / (t ^ 2 + a ^ 2)) - pi / 4 * log(a)
lhs:
    substitute t / a for x
    simplify
    rewrite 1 / (t ^ 2 / a ^ 2 + 1) * log(t / a + 1) to log(t / a + 1) * a ^ 2 / (t ^ 2 + a ^ 2)
    rewrite t / a + 1 to (t + a) / a
    simplify
    rewrite log((a + t) / a) to log(a + t) - log(a)
    rewrite 1 / (a ^ 2 + t ^ 2) * (log(a + t) - log(a)) to log(a + t) / (a ^ 2 + t ^ 2) - log(a) / (a ^ 2 + t ^ 2)
    simplify
    apply integral identity
    simplify
    expand polynomial
done
subgoal 2: a * (INT t:[0,a]. log(t + a) / (t ^ 2 + a ^ 2)) - pi / 4 * log(a) = pi * log(2) / 8
lhs:
    apply 1 on a * (INT t:[0,a]. log(t + a) / (t ^ 2 + a ^ 2)) - pi / 4 * log(a)
    apply integral identity
done
from 2:
    solve equation for INT t:[0,a]. log(t + a) / (t ^ 2 + a ^ 2)
    rewrite pi * log(a) / 4 to 1/8 * pi * (2 * log(a))
    rewrite 2 * log(a) to log(a ^ 2)
    rewrite 1/8 * pi * log(a ^ 2) + pi * log(2) / 8 to 1/8 * pi * (log(2) + log(a ^ 2))
    rewrite log(2) + log(a ^ 2) to log(2 * a ^ 2)
done

## Chapter 2, Section 3, Two Old Tricks, Plus a New One

prove (INT x:[0,oo]. 1 / (x ^ 4 + 2 * x ^ 2 * cosh(2 * a) + 1)) = pi / (4 * cosh(a))
lhs:
    expand definition for cosh (all)
    rewrite x ^ 4 + 2 * x ^ 2 * ((exp(-(2 * a)) + exp(2 * a)) / 2) + 1 to (x ^ 2 + exp(2 * a)) * (x ^ 2 + exp(-(2 * a)))
    rewrite 1 / ((x ^ 2 + exp(2 * a)) * (x ^ 2 + exp(-(2 * a)))) to 1 / (exp(2 * a) - exp(-(2 * a))) * (1 / (x ^ 2 + exp(-(2 * a))) - 1 / (x ^ 2 + exp(2 * a)))
    simplify
    rewrite exp(-(2 * a)) to exp(-a) ^ 2
    rewrite exp(-(2 * a)) to exp(-a) ^ 2
    rewrite exp(2 * a) to exp(a) ^ 2
    rewrite exp(2 * a) to exp(a) ^ 2
    apply integral identity
    simplify
    rewrite to pi / (4 * ((exp(a) + exp(-a)) / 2))
    fold definition for cosh (all)
done

prove (INT x:[0,oo]. 1/(x^4+2*x^2*cos(2*a)+1)) = pi/abs((4*cos(a))) for a: real, cos(a) != 0
subgoal c1: x^4 + 2*x^2*cos(2*a) + 1 != 0 for cos(a) != 0
case analysis on x != 0
    case true:
    lhs:
        rewrite to (x^2 - 1)^2 + 2*x^2*(1 + cos(2*a))
        rewrite cos(2*a) to 2*cos(a)^2 - 1
        simplify
    done
    case false:
    lhs:
        simplify                    
    done
done
subgoal c2: (x^2 - 2*x*sin(a) + 1) * (x^2 + 2*x*sin(a) + 1) != 0 for cos(a) != 0
case analysis on x != 0
    case true:
    lhs:
        expand polynomial
        rewrite sin(a)^2 to 1 - cos(a)^2
        simplify
        rewrite to (x^2 - 1) ^ 2 + 4*x^2*(cos(a)^2)
    done
    case false:
    lhs:
        simplify
    done
done
subgoal 1: (INT x:[0,oo]. x^2 / (x ^ 4 + 2 * x^2* cos(2 * a) + 1)) = (INT x:[0,oo]. 1 / (x^4+2*x^2*cos(2*a)+1))
rhs:
    substitute y for 1/x
    rewrite 1 / (y ^ 2 * (2 * cos(2 * a) / y ^ 2 + 1 / y ^ 4 + 1)) to (1/y ^ 2)/ (2 * cos(2 * a) / y ^ 2 + 1 / y ^ 4 + 1)
    rewrite 1 / y ^ 2 / (2 * cos(2 * a) / y ^ 2 + 1 / y ^ 4 + 1) to (y^4*(1 / y ^ 2)) / (y^4*(2 * cos(2 * a) / y ^ 2 + 1 / y ^ 4 + 1))
    rewrite y ^ 4 * (1 / y ^ 2) / (y ^ 4 * (2 * cos(2 * a) / y ^ 2 + 1 / y ^ 4 + 1)) to y^2/(y^4+2*y^2*cos(2*a)+1)
    substitute x for y
done
subgoal 2: (INT x:[0,oo]. 1/(x^4+2*x^2*cos(2*a)+1)) = 1/2*(INT x:[0,oo]. (1 + x^2)/(x^4+2*x^2*cos(2*a)+1))
rhs:
    rewrite (1 + x^2)/(x^4+2*x^2*cos(2*a)+1) to (1/(x^4+2*x^2*cos(2*a)+1) + x^2/(x^4+2*x^2*cos(2*a)+1))
    simplify
    rewrite (2 * x ^ 2 * cos(2 * a) + x ^ 4 + 1) to (x ^ 4 + 2 * x^2* cos(2 * a) + 1)
    apply 1 on (INT x:[0,oo]. x ^ 2 / (x ^ 4 + 2 * x^2* cos(2 * a) + 1))
    simplify
done
subgoal 3: (INT x:[0,oo]. 1/(x^4+2*x^2*cos(2*a)+1)) = 1/4*(INT x:[-oo,oo]. (1 + x^2)/(x^4+2*x^2*cos(2*a)+1))
rhs:
    split region at 0
    substitute u for -x
    substitute x for u
    simplify
    rewrite (INT x:[0,oo]. (x ^ 2 + 1) / (2 * x ^ 2 * cos(2 * a) + x ^ 4 + 1)) to (INT x:[0,oo]. (1 + x^2)/(x^4+2*x^2*cos(2*a)+1))
    apply 2 on (INT x:[0,oo]. (1 + x ^ 2) / (x ^ 4 + 2 * x ^ 2 * cos(2 * a) + 1))
    simplify
    rewrite to (INT x:[0,oo]. 1/(x^4+2*x^2*cos(2*a)+1))
done
subgoal 4: (INT x:[-oo,oo]. 2*x*sin(a) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1))) = -(INT x:[-oo,oo]. 2*x*sin(a) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1)))
lhs:
    substitute u for -x
    substitute x for u
    rewrite (INT x:[-oo,oo]. -(2 * x * sin(a) / ((2 * x * sin(a) + x ^ 2 + 1) * (-(2 * x * sin(a)) + x ^ 2 + 1)))) to -(INT x:[-oo,oo]. 2*x*sin(a) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1)))
done
subgoal 5: (INT x:[-oo,oo]. 2*x*sin(a) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1))) = 0
lhs:
    rewrite to 1/2*(INT x:[-oo,oo]. 2*x*sin(a) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1)))+1/2*(INT x:[-oo,oo]. 2*x*sin(a) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1)))
    apply 4 on (INT x:[-oo,oo]. 2*x*sin(a) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1)))
    rewrite to 1/2*((INT x:[-oo,oo]. 2*x*sin(a) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1))) - (INT x:[-oo,oo]. 2*x*sin(a) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1))))
    simplify
done
subgoal 6: (INT x:[-oo,oo]. (1 + x ^ 2) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1))) = (INT x:[-oo,oo]. (1 + 2*x*sin(a) + x ^ 2) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1)))
rhs:
    expand polynomial
    rewrite (-(4 * x ^ 2 * sin(a) ^ 2) + 2 * x ^ 2 + x ^ 4 + 1) to ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1))
    simplify
    rewrite 2 * sin(a) * (INT x:[-oo,oo]. x / ((2 * x * sin(a) + x ^ 2 + 1) * (-(2 * x * sin(a)) + x ^ 2 + 1))) to (INT x:[-oo,oo]. (2*x*sin(a)) / ((2 * x * sin(a) + x ^ 2 + 1) * (-(2 * x * sin(a)) + x ^ 2 + 1)))
    rewrite ((2 * x * sin(a) + x ^ 2 + 1) * (-(2 * x * sin(a)) + x ^ 2 + 1)) to ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1))
    apply 5 on (INT x:[-oo,oo]. 2 * x * sin(a) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1)))
    rewrite to (INT x:[-oo,oo]. (x ^ 2 / (-(4 * x ^ 2 * sin(a) ^ 2) + 2 * x ^ 2 + x ^ 4 + 1) + 1 / (-(4 * x ^ 2 * sin(a) ^ 2) + 2 * x ^ 2 + x ^ 4 + 1)))
    rewrite x ^ 2 / (-(4 * x ^ 2 * sin(a) ^ 2) + 2 * x ^ 2 + x ^ 4 + 1) + 1 / (-(4 * x ^ 2 * sin(a) ^ 2) + 2 * x ^ 2 + x ^ 4 + 1) to (1 + x ^ 2) / (-(4 * x ^ 2 * sin(a) ^ 2) + 2 * x ^ 2 + x ^ 4 + 1)
    rewrite (-(4 * x ^ 2 * sin(a) ^ 2) + 2 * x ^ 2 + x ^ 4 + 1) to ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1))
done
subgoal 7: (INT x:[0,oo]. 1/(x^4+2*x^2*cos(2*a)+1)) = pi/(4*cos(a)) for cos(a)>0
lhs:
    apply 3 on (INT x:[0,oo]. 1/(x^4+2*x^2*cos(2*a)+1))
    rewrite cos(2*a) to 1 - 2*(sin(a))^2
    rewrite 2 * x ^ 2 * (1 - 2 * sin(a) ^ 2) to 2*x^2 - 4*x^2*sin(a)^2
    rewrite (x ^ 4 + (2 * x ^ 2 - 4 * x ^ 2 * sin(a) ^ 2) + 1) to (x^2 - 2*x*sin(a)+1)*(x^2+2*x*sin(a)+1)
    apply 6 on (INT x:[-oo,oo]. (1 + x ^ 2) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1)))
    rewrite (1 + 2 * x * sin(a) + x ^ 2) to (x ^ 2 + 2 * x * sin(a) + 1)
    rewrite (x ^ 2 + 2 * x * sin(a) + 1) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1)) to 1 / (x ^ 2 - 2 * x * sin(a) + 1)
    rewrite 1 to sin(a)^2 + cos(a)^2
    rewrite 1 to sin(a)^2 + cos(a)^2
    rewrite (x ^ 2 - 2 * x * sin(a) + (sin(a) ^ 2 + cos(a) ^ 2)) to (x ^ 2 - 2 * x * sin(a) + sin(a) ^ 2 + cos(a) ^ 2)
    rewrite x ^ 2 - 2 * x * sin(a) + sin(a) ^ 2 to (x-sin(a))^2
    rewrite sin(a)^2 + cos(a)^2 to 1
    substitute u for (x - sin(a))
    apply integral identity
    simplify
    rewrite to 1 / (4 * cos(a))*((LIM {u -> oo}. arctan(u / cos(a)))-(LIM {u -> oo}. arctan(-(u / cos(a)))))
    rewrite arctan(-(u / cos(a))) to -arctan((u / cos(a)))
    simplify
done
subgoal 8: (INT x:[0,oo]. 1/(x^4+2*x^2*cos(2*a)+1)) = -pi/(4*cos(a)) for cos(a)<0
lhs:
    apply 3 on (INT x:[0,oo]. 1/(x^4+2*x^2*cos(2*a)+1))
    rewrite cos(2*a) to 1 - 2*(sin(a))^2
    rewrite 2 * x ^ 2 * (1 - 2 * sin(a) ^ 2) to 2*x^2 - 4*x^2*sin(a)^2
    rewrite (x ^ 4 + (2 * x ^ 2 - 4 * x ^ 2 * sin(a) ^ 2) + 1) to (x^2 - 2*x*sin(a)+1)*(x^2+2*x*sin(a)+1)
    apply 6 on (INT x:[-oo,oo]. (1 + x ^ 2) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1)))
    rewrite (1 + 2 * x * sin(a) + x ^ 2) to (x ^ 2 + 2 * x * sin(a) + 1)
    rewrite (x ^ 2 + 2 * x * sin(a) + 1) / ((x ^ 2 - 2 * x * sin(a) + 1) * (x ^ 2 + 2 * x * sin(a) + 1)) to 1 / (x ^ 2 - 2 * x * sin(a) + 1)
    rewrite 1 to sin(a)^2 + cos(a)^2
    rewrite 1 to sin(a)^2 + cos(a)^2
    rewrite (x ^ 2 - 2 * x * sin(a) + (sin(a) ^ 2 + cos(a) ^ 2)) to (x ^ 2 - 2 * x * sin(a) + sin(a) ^ 2 + cos(a) ^ 2)
    rewrite x ^ 2 - 2 * x * sin(a) + sin(a) ^ 2 to (x-sin(a))^2
    rewrite sin(a)^2 + cos(a)^2 to 1
    substitute u for (x - sin(a))
    apply integral identity
    simplify
    rewrite to 1 / (4 * cos(a))*((LIM {u -> oo}. arctan(u / cos(a)))-(LIM {u -> oo}. arctan(-(u / cos(a)))))
    rewrite arctan(-(u / cos(a))) to -arctan((u / cos(a)))
    simplify
done
case analysis on cos(a)
    case negative:
    lhs:
        apply 8 on (INT x:[0,oo]. 1/(x^4+2*x^2*cos(2*a)+1))
    rhs:
        simplify
    done
    case positive:
    lhs:
        apply 7 on (INT x:[0,oo]. 1/(x^4+2*x^2*cos(2*a)+1))
    rhs:
        simplify
    done
done

## Chapter 2, Section 4, Euler's Log-Sine integral

prove (INT x:[0,pi/2]. log(a * sin(x))) = pi/2 * log(a/2) for a: real, a>0
subgoal 1: (INT x:[0,pi/2]. log(a * sin(x))) = (INT x:[0,pi/2]. log(a * cos(x)))
lhs:
    substitute y for pi/2-x
done
subgoal 2: (INT x:[0,pi/2]. log(a * sin(2*x))) = (INT x:[0,pi/2]. log(a * sin(x)))
lhs:
    substitute t for 2*x
    simplify
    split region at pi/2
    simplify
    substitute u for pi-t
    simplify
    substitute x for pi-u
done
subgoal 3: 2*cos(x)*sin(x) = sin(2*x)
rhs:
    rewrite to 2*cos(x)*sin(x)
done
subgoal 4: (INT x:[0,pi/2]. log(a * sin(x)))=1/2 * (INT x:[0,pi / 2]. log(a * sin(x))) + pi * log(a) / 4 - pi * log(2) / 4
lhs:
    rewrite to 1/2*((INT x:[0,pi/2]. log(a * sin(x)))+(INT x:[0,pi/2]. log(a * sin(x))))
    apply 1 on (INT x:[0,pi/2]. log(a * sin(x)))
    rewrite to 1/2*(INT x:[0,pi/2]. (log(a * sin(x))+log(a*cos(x))))
    rewrite log(a*cos(x)) to log(a )+ log(cos(x))
    rewrite to 1/2 * (INT x:[0,pi / 2]. log(a * sin(x)) + log(a) + log(cos(x)))
    rewrite log(a * sin(x)) + log(a) to log(a * sin(x)*a)
    rewrite log(a * sin(x) * a) + log(cos(x)) to log(a * sin(x) * a*cos(x))
    rewrite to 1/2 * (INT x:[0,pi / 2]. log(a ^ 2 *1/2*(2 * cos(x) * sin(x))))
    apply 3 on (2 * cos(x) * sin(x))
    rewrite log(a ^ 2 * 1 / 2 * sin(2 * x)) to log(a*1/2*a*sin(2*x))
    rewrite log(a*1/2*a*sin(2*x)) to log(a*sin(2*x)*a*1/2)
    rewrite log(a*sin(2*x)*a*1/2) to log(a*sin(2*x)*a)+log(1/2)
    rewrite log(a * sin(2 * x) * a) to log(a * sin(2 * x))+log(a)
    apply integral identity
    simplify
    apply 2 on (INT x:[0,pi / 2]. log(a * sin(2 * x)))
done
subgoal 5: (INT x:[0,pi / 2]. log(a * sin(x))) = pi * log(a) / 2 - pi * log(2) / 2 
from 4:
    solve equation for INT x:[0,pi / 2]. log(a * sin(x))
done
lhs:
    apply 5 on (INT x:[0,pi / 2]. log(a * sin(x)))
    rewrite to pi/2*(log(a)-log(2))
    rewrite log(a) - log(2) to log(a/2)
done

prove (INT x:[0,pi / 2]. log(sin(x) / x)) = pi / 2 * (1 - log(pi))
lhs:
    rewrite log(sin(x) / x) to log(sin(x)) - log(x)
    simplify
    rewrite log(sin(x)) to log(1 * sin(x))
    apply integral identity
    integrate by parts with u = log(x), v = x
    apply integral identity
    simplify
    expand polynomial
    simplify
rhs:
    expand polynomial
done

prove (INT x:[0,1]. log(x + 1 / x) / (x ^ 2 + 1)) = pi / 2 * log(2)
subgoal 1: (INT x:[0,oo]. log(x ^ 2 + 1) / (x ^ 2 + 1)) = pi * log(2)
lhs:
    substitute tan(u) for x
    rewrite sec(u) ^ 2 to tan(u) ^ 2 + 1
    simplify
    rewrite tan(u) ^ 2 + 1 to sec(u) ^ 2
    rewrite sec(u) to cos(u) ^ (-1)
    simplify
    substitute x for pi / 2 - u
    rewrite sin(x) to 1 * sin(x)
    apply integral identity
    simplify
done
from 1:
    split region at 1
    substitute y for 1 / x (at 2)
    rewrite y ^ 2 * (1 / y ^ 2 + 1) to y ^ 2 + 1
    rewrite 1 / (y ^ 2 + 1) * log(1 / y ^ 2 + 1) to log(1 / y ^ 2 + 1) / (y ^ 2 + 1)
    rewrite (INT y:[0,1]. log(y ^ 2 + 1) / (y ^ 2 + 1)) + (INT y:[0,1]. log(1 / y ^ 2 + 1) / (y ^ 2 + 1)) to INT y:[0,1]. log(y ^ 2 + 1) / (y ^ 2 + 1) + log(1 / y ^ 2 + 1) / (y ^ 2 + 1)
    rewrite log(y ^ 2 + 1) / (y ^ 2 + 1) + log(1 / y ^ 2 + 1) / (y ^ 2 + 1) to (log(y ^ 2 + 1) + log(1 / y ^ 2 + 1)) / (y ^ 2 + 1)
    rewrite log(y ^ 2 + 1) + log(1 / y ^ 2 + 1) to log((y ^ 2 + 1) * (1 / y ^ 2 + 1))
    rewrite (y ^ 2 + 1) * (1 / y ^ 2 + 1) to (y + 1 / y) ^ 2
    rewrite log((y + 1 / y) ^ 2) to 2 * log(y + 1 / y)
    simplify
    rewrite 1 / (y ^ 2 + 1) * log(1 / y + y) to log(y + 1 / y) / (y ^ 2 + 1)
    solve equation for INT y:[0,1]. log(y + 1 / y) / (y ^ 2 + 1)
done

prove (INT x:[0,oo]. log(x) / (x ^ 2 - b * x + 1)) = 0 for b: real, b > -2, b < 2
subgoal 1: x ^ 2 - b * x + 1 != 0
lhs:
    rewrite x ^ 2 - b * x + 1 to (x - 1/2 * b) ^ 2 + 1 - 1/4 * b ^ 2
done
subgoal 2: (INT x:[0,oo]. log(x ^ a + 1) / (x ^ 2 - b * x + 1)) = (INT x:[0,oo]. log(x ^ a + 1) / (x ^ 2 - b * x + 1)) - a * (INT x:[0,oo]. log(x) / (x ^ 2 - b * x + 1)) for a > 0
lhs:
    substitute 1 / u for x
    simplify
    expand polynomial
    rewrite (1 / u) ^ a to 1 ^ a / u ^ a
    rewrite 1 ^ a / u ^ a + 1 to (1 + u ^ a) / u ^ a
    rewrite log((1 + u ^ a) / u ^ a) to log(1 + u ^ a) - log(u ^ a)
    expand polynomial
    simplify
done
from 2:
    solve equation for INT x:[0,oo]. log(x) / (x ^ 2 - b * x + 1)
done

prove (INT x:[0,1]. (1 - x) / (1 + x + x ^ 2)) = sqrt(3) * pi / 6 - log(3) / 2
lhs:
    rewrite 1 + x + x ^ 2 to (x + 1/2) ^ 2 + 3/4
    substitute u for 2 * (x + 1/2) / sqrt(3)
    rewrite 3 * u ^ 2 / 2 + 3/2 to 3/2 * (u ^ 2 + 1)
    simplify
    rewrite 1 / (u ^ 2 + 1) * (-(u * sqrt(3) / 2) + 3/2) to -sqrt(3) / 2 * (u / (u ^ 2 + 1)) + 3/2 * (1 / (u ^ 2 + 1))
    apply integral identity
    simplify
    substitute t for u ^ 2 + 1
    apply integral identity
    simplify
    expand polynomial
    simplify
done

## Chapter 2, Section 5, Challenge Problems

// Problem C2.1

prove (INT x:[0,4]. log(x) / sqrt(4 * x - x ^ 2)) = 0
subgoal 1: (INT y:[0,1]. 1 / (sqrt(y) * sqrt(1 - y))) = pi
lhs:
    substitute sin(x) ^ 2 for y
    rewrite sin(x) ^ 2 to 1 - cos(x) ^ 2 (at 2)
    simplify
    apply integral identity
    simplify
done
subgoal 2: (INT y:[0,1]. log(y) / (sqrt(y) * sqrt(1 - y))) = -(2 * pi * log(2))
lhs:
    substitute sin(x) ^ 2 for y
    rewrite log(sin(x) ^ 2) to 2 * log(sin(x))
    rewrite sin(x) ^ 2 to 1 - cos(x) ^ 2 (at 2)
    simplify
    rewrite sin(x) to 1 * sin(x)
    apply integral identity
    simplify
done
subgoal 3: 4 * x - x ^ 2 >= 0 for x > 0, x < 4
lhs:
    rewrite 4 * x - x ^ 2 to x * (4 - x)
done
subgoal 4: sqrt(4 * x - x ^ 2) != 0 for x > 0, x < 4
lhs:
    rewrite 4 * x - x ^ 2 to x * (4 - x)
done
lhs:
    substitute y for x / 4
    rewrite log(4 * y) to log(4) + log(y)
    rewrite sqrt(-(16 * y ^ 2) + 16 * y) to 4 * sqrt(-(y ^ 2) + y)
    rewrite sqrt(-(y ^ 2) + y) to sqrt(y) * sqrt(1 - y)
    expand polynomial
    simplify
    rewrite -y + 1 to 1 - y
    rewrite -y + 1 to 1 - y
    apply 1 on INT y:[0,1]. 1 / (sqrt(y) * sqrt(1 - y))
    apply 2 on INT y:[0,1]. log(y) / (sqrt(y) * sqrt(1 - y))
    simplify
done

// Problem C2.2

prove (INT x:[0,1]. (x - 2) / (x ^ 2 - x + 1)) = -pi/sqrt(3)
subgoal 1:(INT u:[-1/2,0]. u / (u ^ 2 + 3/4)) = -(INT u:[0,1/2]. u / (u ^ 2 + 3/4))
lhs:
    substitute t for -u
    simplify
done
subgoal 2:(INT u:[-1/2,1/2]. u/(u^2+3/4)) = 0
lhs:
    split region at 0
    apply 1 on INT u:[-1/2,0]. u / (u ^ 2 + 3/4)
    simplify
done
subgoal 3:3/2*(INT u:[-1/2,1/2]. 1/(u^2+3/4)) = pi/sqrt(3)
lhs:
    simplify
    rewrite 1 / (u ^ 2 + 3/4) to (4/3)/((4/3)*u^2+(4/3)*(3/4))
    simplify
    substitute t for (2*u)/sqrt(3)
    rewrite sqrt(3) / (2 * t ^ 2 + 2) to (sqrt(3))/2*(1/(t^2+1))
    simplify
    apply integral identity
    simplify
    rewrite to pi/sqrt(3)
done
lhs:
    substitute u for x-1/2
    rewrite (u + 1/2) ^ 2 - u + 1/2 to u^2+3/4
    expand polynomial
    simplify
    apply 2 on INT u:[-1/2,1/2]. u / (u ^ 2 + 3/4)
    apply 3 on 3/2 * (INT u:[-1/2,1/2]. 1 / (u ^ 2 + 3/4))
    rewrite to -pi/sqrt(3)
done

// Problem C2.3

prove (INT x:[0,oo]. 1 / (x ^ 4 + 1) ^ (m + 1)) = (4 * m - 1) / (4 * m) * (INT x:[0,oo]. 1 / (x ^ 4 + 1) ^ m) for m: int, m >= 1
subgoal 1: (INT x:[0,oo]. (x ^ 4 + 1) ^ -m) = 4 * m * ((INT x:[0,oo]. 1 / (x ^ 4 + 1) ^ m) - (INT x:[0,oo]. 1 / (x ^ 4 + 1) ^ (m + 1)))
lhs:
    integrate by parts with u = 1 / (x ^ 4 + 1) ^ m, v = x
    simplify
    rewrite x ^ 4 * (x ^ 4 + 1) ^ (-m - 1) to (x ^ 4 + 1) / (x ^ 4 + 1) ^ (m + 1) - 1 / (x ^ 4 + 1) ^ (m + 1)
    rewrite INT x:[0,oo]. (x ^ 4 + 1) / (x ^ 4 + 1) ^ (m + 1) - 1 / (x ^ 4 + 1) ^ (m + 1) to (INT x:[0,oo]. 1 / (x ^ 4 + 1) ^ m) - (INT x:[0,oo]. 1 / (x ^ 4 + 1) ^ (m + 1))
done
from 1:
    solve equation for INT x:[0,oo]. 1 / (x ^ 4 + 1) ^ (m + 1)
    rewrite -(1 / (4 * m) * (INT x:[0,oo]. (x ^ 4 + 1) ^ -m)) + (INT x:[0,oo]. (x ^ 4 + 1) ^ -m) to (4 * m - 1) / (4 * m) * (INT x:[0,oo]. 1 / (x ^ 4 + 1) ^ m)
done

// Problem C2.5

prove (INT x:[0,oo]. log(x + 1) / x ^ (3/2)) = 2 * pi
lhs:
    integrate by parts with u = log(1 + x), v = -2 / sqrt(x)
    simplify
    substitute t for sqrt(x)
    simplify
    apply integral identity
    simplify
done
