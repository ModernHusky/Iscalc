## Examples on integration by parts

// Source:
// The Calculus Page Problems List by D. A. Kouba
// Problems on integration by parts
// URL: https://www.math.ucdavis.edu/~kouba/CalcTwoDIRECTORY/intbypartsdirectory/IntByParts.html

calculate INT x:[0, 1]. x*exp(x)
    integrate by parts with u = x, v = exp(x)
    apply integral identity
    simplify
done

calculate INT x:[0, pi/2]. x * sin(x)
    integrate by parts with u = x, v = -cos(x)
    apply integral identity
    simplify
done

calculate INT x:[1, exp(1)]. x * log(x)
    integrate by parts with u = log(x) / 2, v = x ^ 2
    apply integral identity
    simplify
done

calculate INT x:[0, pi/4]. x * cos(3*x)
    substitute u for 3 * x
    simplify
    integrate by parts with u = u, v = sin(u)
    apply integral identity
    simplify
done

calculate INT x:[1, exp(1)]. log(x)/x^5
    integrate by parts with u = log(x), v = -1/(4 * x^4)
    apply integral identity
    simplify
done

calculate INT x:[0, 1/3]. arcsin(3*x)
    integrate by parts with u = arcsin(3*x), v = x
    simplify
    substitute u for -(9 * x ^ 2) + 1
    apply integral identity
    simplify
done

calculate INT x:[1, exp(1)]. log(x)
    integrate by parts with u = log(x), v = x
    apply integral identity
    simplify
done

calculate INT x:[0, 1]. 2*x*arctan(x)
    integrate by parts with u = arctan(x), v = x^2
    simplify
    rewrite x ^ 2 / (x ^ 2 + 1) to 1 - 1 / (x ^ 2 + 1)
    apply integral identity
    simplify
done

calculate INT x:[1, 2]. x^2*exp(3*x)
    integrate by parts with u = x^2/3, v = exp(3*x)
    simplify
    integrate by parts with u = x/3, v = exp(3*x)
    apply integral identity
    simplify
done

calculate INT x:[1, exp(1)]. x^3*log(5*x)
    integrate by parts with u = log(5*x) / 4, v = x^4
    apply integral identity
    simplify
done

calculate INT x:[1, exp(1)]. log(x) ^ 2
    integrate by parts with u = log(x) ^ 2, v = x
    simplify
    integrate by parts with u = log(x), v = x
    apply integral identity
    simplify
done

calculate INT x:[1, exp(1)]. log(x) ^ 3
    integrate by parts with u = log(x) ^ 3, v = x
    simplify
    integrate by parts with u = log(x) ^ 2, v = x
    simplify
    integrate by parts with u = log(x), v = x
    apply integral identity
    simplify
done

calculate INT x:[-2, 1]. x * sqrt(x + 3)
    integrate by parts with u = x, v = 2/3 * (x + 3) ^ (3/2)
    simplify
    substitute u for x + 3
    apply integral identity
    simplify
done

calculate INT x:[0, pi/4]. x * sin(x) * cos(x)
    integrate by parts with u = x, v = 1/2 * sin(x) ^ 2
    simplify
    rewrite sin(x) ^ 2 to (1 - cos(2*x)) / 2
    apply integral identity
    simplify
done

calculate INT x:[1, exp(1)]. (log(x)/x)^2
    integrate by parts with u = log(x) ^ 2, v = -1/x
    simplify
    integrate by parts with u = log(x), v = -1/x
    apply integral identity
    simplify
done

calculate INT x:[0, 1]. x ^ 5 * exp(x ^ 3)
    integrate by parts with u = x ^ 3, v = exp(x^3) / 3
    simplify
    substitute u for x ^ 3
    apply integral identity
    simplify
done

calculate INT x:[0, pi/4]. x ^ 3 * cos(x ^ 2)
    integrate by parts with u = x ^ 2, v = sin(x^2) / 2
    simplify
    substitute u for x ^ 2
    apply integral identity
    simplify
done

calculate INT x:[0, 1]. x ^ 7 * sqrt(5 + 3 * x ^ 4)
    integrate by parts with u = x ^ 4, v = 1/18 * (5 + 3*x^4) ^ (3/2)
    substitute u for 3 * x ^ 4 + 5
    apply integral identity
    simplify
done

calculate INT x:[1, 2]. x ^ 3 / (x ^ 2 + 5) ^ 2
    integrate by parts with u = x^2 / 2, v = -1/(x^2 + 5)
    simplify
    substitute u for x ^ 2 + 5
    apply integral identity
    simplify
done

calculate INT x:[0, 1]. exp(6 * x) * sin(exp(3 * x))
    integrate by parts with u = exp(3*x), v = -cos(exp(3*x)) / 3
    simplify
    substitute u for exp(3 * x)
    apply integral identity
    simplify
done

calculate INT x:[0, 1]. (x ^ 3 * exp(x^2))/(x^2 + 1) ^ 2
    integrate by parts with u = x^2 * exp(x^2) / 2, v = -1/(x^2 + 1)
    simplify
    rewrite (2 * x ^ 3 * exp(x ^ 2) + 2 * x * exp(x ^ 2)) / (2 * x ^ 2 + 2) to x * exp(x^2)
    substitute u for x ^ 2
    apply integral identity
    simplify
done

calculate INT x:[0, 1]. exp(x) * cos(x)
    integrate by parts with u = exp(x), v = sin(x)
    simplify
    integrate by parts with u = exp(x), v = -cos(x)
    simplify
    solve integral INT x:[0, 1]. exp(x) * cos(x)
done

calculate INT x:[0, pi/4]. sin(3*x) * cos(5*x)
    integrate by parts with u = sin(3*x), v = 1/5*sin(5*x)
    simplify
    integrate by parts with u = cos(3*x), v = -1/5*cos(5*x)
    simplify
    solve integral INT x:[0, pi/4]. sin(3*x) * cos(5*x)
done
