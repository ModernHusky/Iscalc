# 4.1 Application of basic formulas
# Section B
// page 168

calculate INT x. x^2 / (1 + x^2)
    rewrite x^2/(1+x^2) to 1 - 1/(1+x^2)
    apply integral identity
    simplify
done

calculate INT x. x^4 / (1 + x^2)
    rewrite x^4/(1+x^2) to (x^4-1)/(1+x^2) + 1/(1+x^2)
    rewrite (x^4-1) to (x^2-1)*(x^2+1)
    simplify
    apply integral identity
    simplify
done

calculate INT x. 1 / (2-3*x^2) for x != sqrt(2/3), x != -sqrt(2/3)
    substitute u for (sqrt(3/2)*x)
    substitute sin(v) for u
    rewrite -(6 * sin(v) ^ 2) + 6 to 6*cos(v)^2
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1 / (2+3*x^2) for x != sqrt(2/3)
    substitute u for (sqrt(3/2)*x)
    substitute tan(v) for u
    rewrite 6 * tan(v) ^ 2 + 6 to 6 * sec(v)^2
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. (sqrt(x^2+1)-sqrt(x^2-1)) / sqrt(x^4-1) for x > 1
    rewrite sqrt(x^4-1) to sqrt(x^2+1)*sqrt(x^2-1)
    expand polynomial
    simplify
    apply integral identity
    substitute tan(u) for x
    rewrite tan(u)^2 + 1 to sec(u)^2
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. (exp(3*x) + 1) / (exp(x) + 1)
    rewrite exp(3*x) + 1 to (exp(x) + 1) * (exp(2*x) - exp(x) + 1)
    simplify
    apply integral identity
    simplify
done

calculate INT x. x / (4 + x^2) for x>0
    integrate by parts with u= 1/(2*(4+x^2)), v=4+x^2
    substitute u for 4+x^2
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1 / (1 + sin(x))
    rewrite 1 + sin(x) to (sin(x/2) + cos(x/2))^2
    rewrite sin(x/2)+cos(x/2) to sqrt(2)*sin(x/2+pi/4)
    simplify
    substitute u for x/2+pi/4
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. sqrt((1-x)/(1+x)) for x > -1, x < 1
    rewrite (1-x)/(1+x) to (1-x)^2/(1-x^2)
    simplify
    substitute sin(u) for x
    rewrite -(sin(u) ^ 2) + 1 to cos(u)^2
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1 / (16-x^4)
    partial fraction decomposition
    apply integral identity
    simplify
done

calculate INT x. 1 / (exp(x)-exp(-x)) for x > 0
    rewrite 1 / (exp(x)-exp(-x)) to exp(x)/(exp(2*x)-1)
    substitute u for exp(x)
    simplify
    substitute sec(v) for u
    rewrite sec(v)^2-1 to tan(v)^2
    simplify
    rewrite sec(v) to 1/cos(v)
    rewrite tan(v) to sin(v)/cos(v)
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1/sqrt(exp(2*x)+1)
    substitute u for exp(x)
    simplify
    substitute tan(v) for u
    rewrite tan(v)^2 + 1 to sec(v)^2
    simplify
    rewrite sec(v) to 1/cos(v)
    rewrite tan(v) to sin(v)/cos(v)
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1/sqrt(exp(2*x)-1) for x > 0
    substitute u for exp(x)
    simplify
    substitute sec(v) for u
    rewrite sec(v)^2 - 1 to tan(v)^2
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. sqrt((exp(x)-1)/(exp(x)+1)) for x > 0
    rewrite (exp(x)-1)/(exp(x)+1) to (exp(x)-1)^2/(exp(2*x)-1)
    simplify
    substitute u for exp(x)
    simplify
    substitute sec(v) for u
    rewrite sec(v)^2 - 1 to tan(v)^2
    simplify
    apply integral identity
    replace substitution
    simplify
done
