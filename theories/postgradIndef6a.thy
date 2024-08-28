# 4.6 Partial fraction
# Section A

# page 199

calculate INT x. 1/(x*(x^2+1))
    partial fraction decomposition
    apply integral identity
    substitute u for x^2
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1/(x^4-1)
    partial fraction decomposition
    apply integral identity
    simplify
done

calculate INT x. (x^2+1)/((x+1)^2*(x-1))
    partial fraction decomposition
    apply integral identity
    substitute u for x+1
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1/((x^2+1)*(x^2+x+1))
    partial fraction decomposition
    simplify
    rewrite x^2+x+1 to (x+1)^2-(x+1)+1
    substitute u for x+1
    substitute v for u-1/2
    expand polynomial
    expand polynomial
    simplify
    substitute w for 2*v/sqrt(3)
    rewrite sqrt(3) / (3 * w ^ 2 + 3) to 1/sqrt(3) * (1/(w^2+1))
    simplify
    apply integral identity
    substitute t for v^2
    apply integral identity
    substitute o for x^2
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. (2*x+3)/(x^2+3*x-10)
    partial fraction decomposition
    apply integral identity
    simplify
done

calculate INT x. (x+1)/(x^2-2*x+5)
    substitute u for x-1
    expand polynomial
    simplify
    substitute v for u/2
    rewrite 2 / (4 * v ^ 2 + 4) to 1/2 * (1/(v^2+1))
    apply integral identity
    substitute w for u^2
    apply integral identity
    replace substitution
    simplify
done