imports interesting2

# Inside Interesting Integrals, Chapter 3

## Chapter 3, Section 1, Leibniz's formula

prove (INT x:[0,oo]. 1 / (x ^ 2 + a ^ 2)) = pi / (2 * a) for a: real, a > 0
lhs:
    substitute a * u for x
    simplify
    rewrite 1 / (a ^ 2 * u ^ 2 + a ^ 2) to 1 / (a ^ 2 * (u ^ 2 + 1))
    simplify
    apply integral identity
    simplify
done

prove (INT x:[0,oo]. 1 / (x ^ 2 + a ^ 2) ^ 2) = pi / (4 * a ^ 3) for a: real, a > 0
lhs:
    substitute u for x / a
    simplify
    substitute theta for arctan(u)
    rewrite (a ^ 2 * tan(theta) ^ 2 + a ^ 2) ^ 2 to (a ^ 2 * (tan(theta) ^ 2 + 1)) ^ 2
    rewrite tan(theta) ^ 2 + 1 to sec(theta) ^ 2
    simplify
    rewrite 1 / sec(theta) ^ 2 to cos(theta) ^ 2
    rewrite cos(theta) ^ 2 to (1 + cos(2 * theta)) / 2
    simplify
    apply integral identity
    simplify
done

prove (INT x:[0,oo]. 1 / (x ^ 2 + a ^ 2) ^ 2) = pi / (4 * a ^ 3) for a: real, a > 0
subgoal 1: (INT x:[0,oo]. 1 / (x ^ 2 + a ^ 2)) = pi / (2 * a)
lhs:
    rewrite x^2 + a^2 to a^2 + x^2
    apply integral identity
done
from 1:
    differentiate both sides at a
    simplify
    solve equation for INT x:[0,oo]. 1 / (a ^ 2 + x ^ 2) ^ 2
done

prove (INT x:[0,oo]. 1 / (x ^ 2 + a ^ 2) ^ 3) = 3 * pi / (16 * a ^ 5) for a: real, a > 0
subgoal 1: (INT x:[0,oo]. 1 / (x ^ 2 + a ^ 2) ^ 2) = pi / (4 * a ^ 3)
lhs:
    apply integral identity
done
from 1:
    differentiate both sides at a
    simplify
    solve equation for INT x:[0,oo]. 1 / (a ^ 2 + x ^ 2) ^ 3
done

// Probability integral, (3.1.4)

prove (INT x:[-oo,oo]. exp(-(x ^ 2) / 2)) = sqrt(2 * pi)
let g(t) = (INT x:[0,t]. exp(-(x ^ 2) / 2)) ^ 2
subgoal 1: (INT x:[-oo,oo]. exp(-(x ^ 2) / 2)) = 2 * (LIM {t -> oo}. sqrt(g(t)))
lhs:
    split region at 0
    substitute y for -x
    substitute x for y
    simplify
rhs:
    expand definition for g (all)
    simplify
done
subgoal 2: (D t. g(t) + 2 * (INT y:[0,1]. exp(-(1 + y ^ 2) * t ^ 2 / 2) / (1 + y ^ 2))) = 0 for t > 0
lhs:
    expand definition for g (all)
    simplify
    substitute y for x / t (at 2)
    rewrite exp(t ^ 2 * (-(y ^ 2) - 1) / 2) to exp(1/2 * t ^ 2 * (-(y ^ 2) - 1))
    rewrite 1/2 * t ^ 2 * (-(y ^ 2) - 1) to -1/2 * t ^ 2 * y ^ 2 + 1/2 * t ^ 2 * -1
    simplify
    rewrite exp(-(t ^ 2 * y ^ 2 / 2) - t ^ 2 / 2) to exp(-1/2 * t ^ 2 * y ^ 2) * exp(-1/2 * t ^ 2)
    simplify
    rewrite (-(y ^ 2) - 1) / (y ^ 2 + 1) to -1
    simplify
done
subgoal 3: 2 * (INT y:[0,1]. exp(1/2 * t ^ 2 * (-(y ^ 2) - 1)) * (y ^ 2 + 1) ^ (-1)) + g(t) = SKOLEM_CONST(C) for t > 0
from 2:
    integrate both sides
    apply integral identity
    simplify
done
subgoal 4: pi / 2 = SKOLEM_CONST(C)
from 3:
    apply limit t -> 0 both sides
    simplify
    expand definition for g (all)
    apply integral identity
    simplify
done
subgoal 5: g(t) = -(2 * (INT y:[0,1]. 1 / (y ^ 2 + 1) * exp(t ^ 2 * (-(y ^ 2) - 1) / 2))) + pi / 2 for t > 0
lhs:
    apply 3 on g(t)
    apply 4 on SKOLEM_CONST(C)
    simplify
done
lhs:
    apply 1 on INT x:[-oo,oo]. exp(-(x ^ 2) / 2)
    apply 5 on g(t)
    simplify
done

prove (INT x:[0,oo]. exp(-(x ^ 2) / 2)) = sqrt(2) * sqrt(pi) / 2
subgoal 1: (INT x:[-oo,oo]. exp(-(x ^ 2) / 2)) = sqrt(2 * pi)
lhs:
    apply integral identity
done
from 1:
    split region at 0
    substitute y for -x
    substitute x for y
    simplify
    solve equation for INT x:[0,oo]. exp(-(x ^ 2 / 2))
done

prove (INT x:[-oo,oo]. exp(-(a * x ^ 2))) = sqrt(pi / a) for a > 0
lhs:
    substitute u for sqrt(2 * a) * x
    simplify
    substitute x for u
    rewrite -(x ^ 2 / 2) to -(x ^ 2) / 2
    apply integral identity
    simplify
done

prove (INT x:[0,oo]. exp(-(x ^ 2))) = sqrt(pi) / 2
subgoal 1: (INT x:[0,oo]. exp(-(x ^ 2) / 2)) = sqrt(2) * sqrt(pi) / 2
lhs:
    apply integral identity
done
from 1:
    substitute x for x / sqrt(2)
    simplify
    solve equation for INT x:[0,oo]. exp(-(x ^ 2))
done

// Application, (3.1.8)

prove (INT x:[0,1]. 1 / sqrt(-log(x))) = sqrt(pi)
subgoal 1: (INT x:[0,oo]. exp(-(x ^ 2))) = sqrt(pi) / 2
lhs:
    apply integral identity
done
from 1:
    substitute t for exp(-(x ^ 2))
    simplify
    solve equation for INT x:[0,1]. 1 / sqrt(-log(x))
done

## Chapter 3, Section 2, An Amazing Integral

prove (INT x:[0,oo]. sin(a*x)/x) = pi/2 * sgn(a)
let g(y,a) = INT x:[0,oo]. exp(-x * y) * sin(a * x) / x for y >= 0
subgoal 1: (D y. g(y, a)) = - a / (a ^ 2 + y ^ 2) for y > 0, a != 0
lhs:
    expand definition for g(all)
    exchange derivative and integral
    simplify
    apply integral identity
rhs:
    simplify
done
subgoal 2: g(y, a) = -arctan(y / a) + SKOLEM_FUNC(C(a)) for y >= 0, a > 0
from 1:
    integrate both sides
    apply integral identity
    simplify
done
subgoal 3: g(y, a) = -arctan(y / a) + SKOLEM_FUNC(C(a)) for y >= 0, a < 0
from 1:
    integrate both sides
    apply integral identity
    simplify
done
subgoal 4: (LIM {y -> oo}. g(y, a)) = 0 for y >= 0
lhs:
    expand definition for g(all)
    simplify
done
subgoal 5: SKOLEM_FUNC(C(a)) = pi / 2 for a > 0
from 2:
    apply limit y -> oo both sides
    apply 4 on LIM {y -> oo}. g(y,a)
    simplify
    solve equation for SKOLEM_FUNC(C(a))
done
subgoal 6: SKOLEM_FUNC(C(a)) = -pi / 2 for a < 0
from 3:
    apply limit y -> oo both sides
    apply 4 on LIM {y -> oo}. g(y,a)
    simplify
    solve equation for SKOLEM_FUNC(C(a))
done
subgoal 7: g(0,a) = pi / 2 for a > 0
from 2:
    apply limit y -> 0 both sides
    simplify
    apply 5 on SKOLEM_FUNC(C(a))
done
subgoal 8: g(0,a) = -pi / 2 for a < 0
from 3:
    apply limit y -> 0 both sides
    simplify
    apply 6 on SKOLEM_FUNC(C(a)) 
done
case analysis on a
case positive:
    lhs:
        rewrite sin(a*x)/x to exp(-x*0)*sin(a*x)/x
        fold definition for g
        apply 7 on g(0,a)
    done
case zero:
    lhs:
        simplify
    done
case negative:
    lhs:
        rewrite sin(a*x)/x to exp(-x*0)*sin(a*x)/x
        fold definition for g
        apply 8 on g(0,a)
    done
done

## Chapter 3, Section 3, Frullani's Integral

prove (INT x:[0,oo]. (arctan(a * x) - arctan(b * x)) / x) = pi * log(a) / 2 - pi * log(b) / 2 for a b: real, a > 0, b > 0
let I(a,b) = (INT x:[0,oo]. (arctan(a * x) - arctan(b * x)) / x)
subgoal 1: (D a. I(a,b)) = pi / (2 * a)
lhs:
    expand definition for I (all)
    exchange derivative and integral
    simplify
    substitute u for a * x
    apply integral identity
    simplify
done
subgoal 2: I(a,b) = pi * log(a) / 2 + SKOLEM_FUNC(C(b))
from 1:
    integrate both sides
    simplify
    apply integral identity
    simplify
done
subgoal 3: SKOLEM_FUNC(C(a)) = -(pi * log(a) / 2)
from 2:
    substitute b for a in equation
    solve equation for SKOLEM_FUNC(C(a))
    expand definition for I (all)
    simplify
done
lhs:
    fold definition for I (all)
    apply 2 on I(a,b)
    apply 3 on SKOLEM_FUNC(C(b))
    simplify
done

## Chapter 3, Section 4, The Flip-Side of Feynman's Trick

// (3.4.2)

prove (INT t:[0,oo]. (exp(-p*t^2)-exp(-q*t^2))/t^2) = sqrt(pi)*(sqrt(q) - sqrt(p)) for p q: real, p > 0, q > 0
subgoal 1: (INT t:[0,oo]. (exp(-p*t^2)-exp(-q*t^2))/t^2) = (INT t:[0,oo]. (INT a:[p,q]. exp(-a*t^2)))
rhs:
    substitute x for -a*t (at 2)
    apply integral identity
    simplify
    rewrite 1 / t * (exp(-(q * t ^ 2)) / t - exp(-(p * t ^ 2)) / t) to exp(-q * t ^ 2) / t^2 - exp(-p * t ^ 2) / t^2
    rewrite to (INT t:[0,oo]. exp(-p * t ^ 2) / t ^ 2 - exp(-q * t ^ 2) / t ^ 2)
    rewrite exp(-p * t ^ 2) / t ^ 2 - exp(-q * t ^ 2) / t ^ 2 to (exp(-p*t^2)-exp(-q*t^2))/t^2
done
lhs:
    apply 1 on (INT t:[0,oo]. (exp(-p*t^2)-exp(-q*t^2))/t^2)
    exchange integral and integral
    substitute x for sqrt(2*a)*t (at 2)
    apply integral identity
    simplify
    rewrite -(x^2/2) to -(x^2)/2
    apply integral identity
    rewrite to sqrt(pi)*(sqrt(q) - sqrt(p))
done

// (3.4.3)

prove (INT x:[0,1]. (x ^ a - 1) / log(x)) = log(a + 1) for a: real, a > -1
let I(a) = INT x:[0, 1]. (x ^ a - 1) / log(x)
subgoal 1: (D a. I(a)) = 1 / (a + 1)
lhs:
    expand definition for I (all)
    exchange derivative and integral
    simplify
    apply integral identity
    simplify
done
subgoal 2: I(a) = log(a + 1) + SKOLEM_CONST(C)
from 1:
    integrate both sides
    apply integral identity
    simplify
done
subgoal 3: SKOLEM_CONST(C) = 0
from 2:
    substitute a for 0 in equation
    expand definition for I (all)
    simplify
    solve equation for SKOLEM_CONST(C)
done
lhs:
    fold definition for I (all)
    apply 2 on I(a)
    apply 3 on SKOLEM_CONST(C)
    simplify
done

// (3.4.4)

prove (INT x:[0,1]. (x ^ a - x ^ b) / log(x)) = log((a + 1) / (b + 1)) for a b: real, a > -1, b > -1
lhs:
    rewrite x ^ a - x ^ b to x ^ a - 1 - (x ^ b - 1)
    rewrite (x ^ a - 1 - (x ^ b - 1)) / log(x) to (x ^ a - 1) / log(x) - (x ^ b - 1) / log(x)
    simplify
    apply integral identity
    rewrite log(a + 1) - log(b + 1) to log((a + 1) / (b + 1))
done

// (3.4.5)

prove (INT x:[0, oo]. exp(-(t*x)) * (cos(a*x) - cos(b*x)) / x) = log(sqrt((t^2+b^2)/(t^2+a^2))) for a b t: real, a>0, b>0, t>0
subgoal 1: (INT s:[a,b]. sin(x*s)) = (cos(a*x)-cos(b*x))/x for x>0
lhs:
    apply integral identity
    simplify
    rewrite cos(a * x) / x - cos(b * x) / x to (cos(a*x)-cos(b*x))/x
done
subgoal 2: log(sqrt((b ^ 2 + t ^ 2) / (a ^ 2 + t ^ 2))) = 1/2*log((b ^ 2 + t ^ 2) / (a ^ 2 + t ^ 2))
lhs:
    rewrite sqrt((b ^ 2 + t ^ 2) / (a ^ 2 + t ^ 2)) to ((b ^ 2 + t ^ 2) / (a ^ 2 + t ^ 2))^(1/2)
    simplify
done
lhs:
    rewrite exp(-(t*x)) * (cos(a*x) - cos(b*x)) / x to exp(-(t*x)) * ((cos(a*x) - cos(b*x)) / x)
    apply 1 on (cos(a*x)-cos(b*x))/x
    rewrite INT x:[0,oo]. exp(-(t * x)) * (INT s:[a,b]. sin(x * s)) to INT x:[0,oo]. (INT s:[a,b]. exp(-(t * x)) * sin(s * x))
    exchange integral and integral
    rewrite -(t*x) to -(x*t)
    apply integral identity
    substitute u for s^2 + t^2
    apply integral identity
    simplify
    rewrite log(b ^ 2 + t ^ 2) / 2 - log(a ^ 2 + t ^ 2) / 2 to 1/2*log(b ^ 2 + t ^ 2) - 1/2*log(a ^ 2 + t ^ 2)
    rewrite to 1/2*(log(b ^ 2 + t ^ 2) - log(a ^ 2 + t ^ 2))
    rewrite log(b ^ 2 + t ^ 2) - log(a ^ 2 + t ^ 2) to log((b ^ 2 + t ^ 2)/(a ^ 2 + t ^ 2))
    apply 2 on log((b ^ 2 + t ^ 2) / (a ^ 2 + t ^ 2))
    simplify
    rewrite to log(sqrt((t^2+b^2)/(t^2+a^2)))
done

// (3.4.7)

prove (INT x:[0,1]. x^a * (log(x))^2) = 2/(a+1)^3 for a: real, a > -1
subgoal 1: (D a. (D a. (INT x:[0,1]. x^a))) = 2/(a+1)^3
lhs:
    apply integral identity
    simplify
done
from 1:
    simplify
done

## Chapter 3, Section 10, Challenge Problems

// C3.1

prove (INT x:[0,oo]. log(1 + a ^ 2 * x ^ 2) / (b ^ 2 + x ^ 2)) = pi / b * log(1 + a * b) for a b: real, a > 0, b > 0
let I(a,b) = (INT x:[0,oo]. log(1 + a ^ 2 * x ^ 2) / (b ^ 2 + x ^ 2))
subgoal 1: (D a. I(a,b)) = pi / (1 + a * b)
lhs:
    expand definition for I (all)
    exchange derivative and integral
    simplify
    rewrite x ^ 2 / ((b ^ 2 + x ^ 2) * (a ^ 2 * x ^ 2 + 1)) to 1 / (1 - a ^ 2 * b ^ 2) * (1 / (1 + a ^ 2 * x ^ 2) - b ^ 2 / (b ^ 2 + x ^ 2))
    simplify
    apply integral identity
    simplify
    rewrite to pi / (1 + a * b)
done
subgoal 2: I(a,b) = pi / b * log(1 + a * b) + SKOLEM_FUNC(C(b))
from 1:
    integrate both sides
    substitute u for 1 + a * b
    simplify
    apply integral identity
    replace substitution
    rewrite abs(1 + a * b) to 1 + a * b
done
subgoal 3: I(0,b) = 0
lhs:
    expand definition for I
    simplify
done
subgoal 4: SKOLEM_FUNC(C(b)) = 0
from 2:
    apply limit a -> 0 both sides
    simplify
    apply 3 on I(0,b)
    solve equation for SKOLEM_FUNC(C(b))
done
from 2:
    apply 4 on SKOLEM_FUNC(C(b))
    expand definition for I (all)
    simplify
done

// C3.5

prove (INT x:[0, oo]. cos(a * x) * sin(b * x) / x) = pi/4 + pi/4 * sgn(b-a) for a b: real, a > 0, b > 0
let I(a, b) = (INT x:[0, oo]. cos(a * x) * sin(b * x) / x)
subgoal 1: I(a, b) = 1/2 * (INT x:[0, oo]. sin((b + a) * x) / x) + 1/2 * (INT x:[0, oo]. sin((-a + b) * x) / x)
lhs:
    expand definition for I
    rewrite cos(a * x) * sin(b * x) to 1/2 * (sin(b * x + a * x) - sin(a * x - b * x))
    rewrite 1/2 * (sin(b * x + a * x) - sin(a * x - b * x)) / x to 1/2 * sin((b + a) * x) / x - 1/2 * sin(-((b - a) * x)) / x
    simplify
done
case analysis on -a + b
    case positive:
    lhs:
        fold definition for I
        apply 1 on I(a,b)
        apply integral identity
    done
    case zero:
    lhs:
        fold definition for I
        apply 1 on I(a,b)
        apply integral identity
    done
    case negative:
    lhs:
        fold definition for I
        apply 1 on I(a,b)
        apply integral identity
    done
done

// C3.6

prove (INT x:[-1,1]. sqrt((1 + x) / (1 - x))) = pi
lhs:
    substitute cos(2 * u) for x
    rewrite cos(2 * u) to 2 * cos(u) ^ 2 - 1 (at 1)
    rewrite cos(2 * u) to 1 - 2 * sin(u) ^ 2
    simplify
    rewrite sin(2 * u) to 2 * sin(u) * cos(u)
    simplify
    rewrite cos(u) ^ 2 to 1/2 * (1 + cos(2*u))
    apply integral identity
    simplify
done

// C3.7a

prove (INT x:[-oo,oo]. x * exp(-(x ^ 2) - x)) = -1/2 * sqrt(pi * sqrt(exp(1))) for x:real
let I(a,b) = (INT x:[-oo,oo]. exp(-a * x ^ 2 + b * x)) for a b: real, a > 0
subgoal 1: I(a,b) = exp(b ^ 2 / (4 * a)) * sqrt(pi / a) for a b: real, a > 0
lhs:
    expand definition for I
    rewrite b * x - a * x ^ 2 to b ^ 2 / (4 * a) - a * (x - b / (2 * a)) ^ 2
    rewrite exp(b ^ 2 / (4 * a) - a * (x - b / (2 * a)) ^ 2) to exp(b ^ 2 / (4 * a)) * exp(-a * (x - b / (2 * a)) ^ 2)
    simplify
    substitute y for x - b / (2 * a)
    apply integral identity
    simplify
done
subgoal 2: (D b. I(a,b)) = b / (2 * a) * exp(b ^ 2 / (4 * a)) * sqrt(pi / a) for a b: real, a > 0
lhs:
    apply 1 on I(a,b)
    simplify
    rewrite a ^ (3/2) to a*sqrt(a)
    rewrite (2 * (a * sqrt(a))) to (2 * a * sqrt(a))
    rewrite b * sqrt(pi) / (2 * a * sqrt(a)) to b * 1/sqrt(a) * sqrt(pi) / (2 * a)
    rewrite b * 1 / sqrt(a) * sqrt(pi) to 1 / sqrt(a) * sqrt(pi) * b
    rewrite 1 / sqrt(a) * sqrt(pi) to sqrt(pi/a)
    rewrite sqrt(pi / a) * b / (2 * a) * exp(b ^ 2 / (4 * a)) to b / (2 * a) * exp(b ^ 2 / (4 * a)) * sqrt(pi / a)
done
subgoal 3: (INT x:[-oo,oo]. x * exp(-(a * x ^ 2) + b * x)) = b / (2 * a) * exp(b ^ 2 / (4 * a)) * sqrt(pi / a) for a b: real, a > 0
from 2:
    expand definition for I (all)
    simplify
done
lhs:
    rewrite x * exp(-(x ^ 2) - x) to x * exp(-(1 * x ^ 2) + -1 * x)
    apply 3 on INT x:[-oo,oo]. x * exp(-(1 * x ^ 2) + -1 * x)
    rewrite to -1/2 * sqrt(pi * sqrt(exp(1)))
done

// C3.7b

prove (INT x:[-oo,oo]. x ^ 2 * exp(-(x ^ 2) - x)) = 3/4 * sqrt(pi * sqrt(exp(1)))
let I(a,b) = (INT x:[-oo,oo]. exp(-a * x ^ 2 + b * x)) for a > 0
subgoal 1: I(a,b) = exp(b ^ 2 / (4 * a)) * sqrt(pi / a) for a > 0
lhs:
    expand definition for I
    rewrite b * x - a * x ^ 2 to b ^ 2 / (4 * a) - a * (x - b / (2 * a)) ^ 2
    rewrite exp(b ^ 2 / (4 * a) - a * (x - b / (2 * a)) ^ 2) to exp(b ^ 2 / (4 * a)) * exp(-a * (x - b / (2 * a)) ^ 2)
    simplify
    substitute y for x - b / (2 * a)
    apply integral identity
    simplify
done
subgoal 2: (D a. I(a,b)) = -(b ^ 2 / (4 * a ^ 2)) * exp(b ^ 2 / (4 * a)) * sqrt(pi / a) - 1 / (2 * a) * exp(b ^ 2 / (4 * a)) * sqrt(pi / a) for a > 0
lhs:
    apply 1 on I(a,b)
    simplify
rhs:
    simplify
done
subgoal 3: (INT x:[-oo,oo]. x ^ 2 * exp(-(a * x ^ 2) + b * x)) = b ^ 2 / (4 * a ^ 2) * exp(b ^ 2 / (4 * a)) * sqrt(pi / a) + 1 / (2 * a) * exp(b ^ 2 / (4 * a)) * sqrt(pi / a) for a > 0
from 2:
    expand definition for I (all)
    simplify
    rewrite INT x:[-oo,oo]. x ^ 2 * exp(b * x - a * x ^ 2) to INT x:[-oo,oo]. x ^ 2 * exp(-(a * x ^ 2) + b * x)
    solve equation for INT x:[-oo,oo]. x ^ 2 * exp(-(a * x ^ 2) + b * x)
done
lhs:
    rewrite x ^ 2 * exp(-(x ^ 2) - x) to x ^ 2 * exp(-(1 * x ^ 2) + -1 * x)
    apply 3 on INT x:[-oo,oo]. x ^ 2 * exp(-(1 * x ^ 2) + -1 * x)
    rewrite to 3/4 * sqrt(pi * sqrt(exp(1)))
done

// C3.8

// Note: rely on integral of sin(m*x) / (x * (x^2 + a^2))

prove (INT x:[0,oo]. sin(m * x) / (x * (a ^ 2 + x ^ 2) ^ 2)) = pi / (2 * a ^ 4) * (1 - (2 + m * a) / 2 * exp(-a * m)) for a > 0, m > 0
subgoal 1: (INT x:[0,oo]. sin(m * x) / (x * (a ^ 2 + x ^ 2))) = pi * (1 - exp(-a * m)) / (2 * a ^ 2)
lhs:
    apply integral identity
done
from 1:
    differentiate both sides at a
    exchange derivative and integral (all)
    simplify
    solve equation for INT x:[0,oo]. sin(m * x) / (x * (a ^ 2 + x ^ 2) ^ 2)
    rewrite -((2 * a ^ 2 * m * pi * exp(-(a * m)) - 4 * a * pi * (1 - exp(-(a * m)))) / (8 * a ^ 5)) to pi / (2 * a ^ 4) * (1 - (2 + m * a) / 2 * exp(-a * m))
done

// C3.9

prove (INT x:[0,1]. x / (a * x + b * (1 - x)) ^ 3) = 1 / (2 * a ^ 2 * b) for a b: real, a > 0, b > 0, a > b
subgoal 1: (INT x:[0,1]. 1 / (a * x + b * (1 - x)) ^ 2) = 1 / (a * b)
lhs:
    substitute u for (a - b) * x + b
    rewrite 1 / ((a - b) * (b * (1 - (u - b) / (a - b)) + a * (u - b) / (a - b)) ^ 2) to 1 / (u ^ 2 * (a - b))
    apply integral identity
    simplify
    rewrite 1 / (a - b) * (1 / b - 1 / a) to 1 / (a * b)
done
from 1:
    differentiate both sides at a
    exchange derivative and integral (all)
    simplify
    rewrite (b * (1 - x) + a * x) ^ 3 to (a * x + b * (1 - x)) ^ 3
    solve equation for INT x:[0,1]. x / (a * x + b * (1 - x)) ^ 3
done
