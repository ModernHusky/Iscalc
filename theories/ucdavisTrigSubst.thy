## Examples on trigonometric substitution

// Source:
// The Calculus Page Problems List by D. A. Kouba
// Problems on trigonometric substitution
// URL: https://www.math.ucdavis.edu/~kouba/CalcTwoDIRECTORY/trigsubdirectory/TrigSub.html

calculate INT x. sqrt(1 - x^2) for x > -1, x < 1
    substitute u for asin(x)
    rewrite -(sin(u) ^ 2) + 1 to 1 - sin(u)^2
    rewrite 1 - sin(u)^2 to cos(u)^2
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. (x^2 - 1)^(3/2) / x for x > 1
    substitute u for asec(x)
    rewrite sec(u)^2 - 1 to tan(u)^2
    simplify
    rewrite tan(u)^4 to tan(u)^2 * tan(u)^2
    rewrite tan(u)^2 to sec(u)^2 - 1 (at 2)
    expand polynomial
    rewrite tan(u)^2 to sec(u)^2 - 1 (at 2)
    simplify
    apply integral identity
    substitute v for tan(u)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1 / (1 - x^2) ^ (3/2) for x > -1, x < 1
    substitute u for asin(x)
    rewrite -(sin(u)^2) + 1 to 1 - sin(u)^2
    rewrite 1 - sin(u)^2 to cos(u)^2
    simplify
    rewrite 1 / cos(u)^2 to (1/cos(u))^2
    rewrite 1/cos(u) to sec(u)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. sqrt(x^2 + 1) / x
    substitute u for atan(x)
    rewrite tan(u)^2 + 1 to sec(u)^2
    simplify
    rewrite sec(u)^3 to sec(u) * sec(u)^2
    rewrite sec(u)^2 to 1 + tan(u)^2
    rewrite sec(u) * (1 + tan(u)^2) / tan(u) to (1 + tan(u)^2) / sin(u)
    expand polynomial
    rewrite 1 / sin(u) to csc(u)
    rewrite tan(u) ^ 2 / sin(u) to sec(u) * tan(u)
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. x^3 * sqrt(4 - 9*x^2) for x > -2/3, x < 2/3
    substitute u for asin(3*x/2)
    simplify
    rewrite -(4 * sin(u)^2) + 4 to 4 * (1 - sin(u)^2)
    rewrite 1 - sin(u)^2 to cos(u)^2
    simplify
    rewrite cos(u)^2 * sin(u)^3 to sin(u) * sin(u)^2 * cos(u)^2
    rewrite sin(u)^2 to 1 - cos(u)^2
    expand polynomial
    simplify
    substitute v for cos(u)
    substitute v2 for cos(u) (at 2)
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. sqrt(1 - x^2) / x for x > 0, x < 1
    substitute u for asin(x)
    rewrite -(sin(u)^2) + 1 to 1 - sin(u)^2
    rewrite 1 - sin(u)^2 to cos(u)^2
    simplify
    rewrite cos(u) ^ 2 to 1 - sin(u)^2
    expand polynomial
    rewrite 1 / sin(u) to csc(u)
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. sqrt(x^2 - 9) / x^2 for x > 3
    substitute u for asec(x/3)
    rewrite 9 * sec(u)^2 - 9 to 9 * (sec(u)^2 - 1)
    rewrite sec(u)^2 - 1 to tan(u)^2
    simplify
    rewrite tan(u)^2 to sec(u)^2 - 1
    expand polynomial
    rewrite 1 / sec(u) to cos(u)
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. sqrt(x^2 + 1) / x^2 for x > 0
    substitute u for atan(x)
    rewrite tan(u)^2 + 1 to sec(u)^2
    simplify
    rewrite sec(u)^3 to sec(u) * sec(u)^2
    rewrite sec(u)^2 to 1 + tan(u)^2
    expand polynomial
    simplify
    rewrite sec(u) / tan(u)^2 to cot(u) * csc(u)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. sqrt(x^2+25)
    substitute u for atan(x/5)
    simplify
    rewrite sqrt(25 * tan(u) ^ 2 + 25) to 5 * sqrt(tan(u)^2+1)
    rewrite tan(u)^2+1 to sec(u)^2
    simplify
    integrate by parts with u=sec(u), v=tan(u)
    rewrite tan(u)^2 to sec(u)^2 - 1
    expand polynomial
    apply integral identity
    solve integral 25 * (INT u. sec(u)^3)
    expand polynomial
    simplify
    replace substitution
    simplify
done

calculate INT x. sqrt(x^2-4) for x > 2
    substitute u for asec(x/2)
    rewrite sqrt(4 * sec(u) ^ 2 - 4) to 2 * sqrt(sec(u) ^ 2 - 1)
    rewrite sec(u) ^ 2 - 1 to tan(u)^2
    simplify
    rewrite tan(u)^2 to sec(u)^2 - 1
    expand polynomial
    simplify
    integrate by parts with u=sec(u), v=tan(u)
    solve integral 4*(INT u. sec(u)*tan(u)^2)
    expand polynomial
    apply integral identity
    simplify
    replace substitution
    simplify
done

calculate INT x. x / sqrt(x^4-16) for x > 2
    substitute u for x^2
    simplify
    substitute v for asec(u/4)
    rewrite sqrt(16 * sec(v) ^ 2 - 16) to 4*sqrt(sec(v)^2-1)
    rewrite sec(v)^2-1 to tan(v)^2
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1 / sqrt(x^2-4*x) for x > 4
    rewrite x^2-4*x to (x-2)^2 - 4
    substitute u for asec((x-2)/2)
    rewrite sqrt(4 * sec(u) ^ 2 - 4) to 2 * sqrt(sec(u) ^ 2 - 1)
    rewrite sec(u) ^ 2 - 1 to tan(u)^2
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. x/sqrt(x^2 + 4*x + 5)
    rewrite x^2+4*x+5 to (x+2)^2 + 1
    substitute u for (x+2)
    substitute v for atan(u)
    rewrite tan(v)^2+1 to sec(v)^2
    simplify
    expand polynomial
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. x * sqrt(10*x - x^2) for x > 0, x < 10
    rewrite 10*x - x^2 to 25 - (x-5)^2
    substitute u for (x-5) / 5
    substitute v for asin(u)
    rewrite sqrt(-(25*sin(v)^2)+25) to 5*sqrt(1-sin(v)^2)
    rewrite 1-sin(v)^2 to cos(v)^2
    simplify
    expand polynomial
    simplify
    substitute w for cos(v)
    apply integral identity
    simplify
    replace substitution
    simplify
done

calculate INT x. sqrt((x-1) / x) for x > 1
    substitute u for sqrt(x)
    simplify
    substitute v for asec(u)
    rewrite sec(v)^2-1 to tan(v)^2
    simplify
    integrate by parts with u=tan(v),v=sec(v)
    rewrite sec(v)^3 to sec(v)*sec(v)^2
    rewrite sec(v)^2 to tan(v)^2 + 1
    expand polynomial
    expand polynomial
    simplify
    solve integral 2 * INT v. sec(v)*tan(v)^2
    expand polynomial
    apply integral identity
    simplify
    replace substitution
    simplify
done

calculate INT x. sqrt(1-x)*sqrt(x+3) for x < 1, x > -3
    rewrite sqrt(1-x)*sqrt(x+3) to sqrt((1-x)*(x+3))
    rewrite (1-x)*(x+3) to 4 - (x+1)^2
    substitute u for(x+1)/2
    rewrite sqrt(-(4 * u ^ 2) + 4) to 2 * sqrt(1 - u ^ 2)
    substitute v for asin(u)
    rewrite -(sin(v) ^ 2) + 1 to 1-sin(v)^2
    rewrite 1-sin(v) ^ 2 to cos(v)^2
    simplify
    integrate by parts with u=cos(v),v=sin(v)
    simplify
    rewrite sin(v)^2 to 1-cos(v)^2
    simplify
    solve integral 4 * (INT v. cos(v) ^ 2)
    expand polynomial
    apply integral identity
    replace substitution
    simplify
done
