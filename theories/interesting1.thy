imports standard

# Inside Interesting Integrals, Chapter 1

## Chapter 1, Section 5, Some Examples of Tricks

prove (INT x:[0,oo]. log(x) / (x ^ 2 + 1)) = 0
lhs:
    split region at 1
    substitute 1 / u for x
    simplify
    rewrite u ^ 2 * (1 / u ^ 2 + 1) to u ^ 2 + 1
    simplify
done

## Chapter 1, Section 7, Dalzell's Integral

prove (INT x:[0,1]. (x^4*(1-x)^4)/(1+x^2)) = 22/7 - pi
lhs:
    rewrite (x^4*(1-x)^4)/(1+x^2) to (x^6-4*x^5+5*x^4-4*x^2+4)-4/(1+x^2)
    simplify
    apply integral identity
    simplify
done

## Chapter 1, Section 10, Challenge Problems

// C1.5

prove (INT x:[0,pi / 3]. 1 / cos(x)) = log(2 + sqrt(3))
lhs:
    rewrite 1 / cos(x) to cos(x) / cos(x) ^ 2
    rewrite cos(x) ^ 2 to 1 - sin(x) ^ 2
    substitute u for sin(x)
    rewrite 1 / (1 - u ^ 2) to 1/2 * (1 / (1 - u) + 1 / (1 + u))
    simplify
    apply integral identity
    simplify
    rewrite 1 / (1 - u) to -1 / (u - 1)
    apply integral identity
    simplify
    rewrite 1/2 * log(sqrt(3) / 2 + 1) - 1/2 * log(1 - sqrt(3) / 2) to 1/2 * (log(sqrt(3) / 2 + 1) - log(1 - sqrt(3) / 2))
    rewrite log(sqrt(3) / 2 + 1) - log(1 - sqrt(3) / 2) to log((sqrt(3) / 2 + 1) / (1 - sqrt(3) / 2))
    rewrite (sqrt(3) / 2 + 1) / (1 - sqrt(3) / 2) to (2 + sqrt(3)) ^ 2
    simplify
    rewrite sqrt(3) + 2 to 2 + sqrt(3)
done
