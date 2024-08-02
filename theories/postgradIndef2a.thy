# 4.2 Substitution
# Section A
// page 171

calculate INT x. exp(-5*x)
    apply integral identity
    simplify
done

calculate INT x. exp(sin(x))*cos(x)
    substitute u for sin(x)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. exp(x^2)*x
    substitute u for x^2
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. exp(tan(x))*sec(x)^2
    substitute u for tan(x)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. exp(sqrt(x)) / sqrt(x)
    substitute u for sqrt(x)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. sin(sin(x))*cos(x)
    substitute u for sin(x)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. sin(log(x))/x
    substitute u for log(x)
    apply integral identity
    replace substitution
    simplify
done

// page 172

calculate INT x. sin(sqrt(x)) / sqrt(x)
    substitute u for sqrt(x)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. sin(exp(x))*exp(x)
    substitute u for exp(x)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. sin(sqrt(1+x^2))*x/sqrt(1+x^2) for x != 0
    substitute u for sqrt(1+x^2)
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1/(2-3*x)^(1/3)
    substitute u for 2-3*x
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. sin(x) / cos(x)^3
    substitute u for cos(x)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. x / sqrt(2-3*x^2)
    substitute u for 2-3*x^2
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. (1+log(x))/(x*log(x))^2 for x > exp(-1)
    substitute u for x*log(x)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. tan(x)^9*sec(x)^2
    substitute u for tan(x)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1 / (asin(x)^2*sqrt(1-x^2)) for x < 1, x > -1
    substitute u for asin(x)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. (sin(x)+cos(x))/(sin(x)-cos(x))^(1/3) for x>0, x < pi/2
    substitute u for sin(x)-cos(x)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1/(1-3*x)
    apply integral identity
    simplify
done

calculate INT x. (x+1)/(x^2+2*x+5) for x != -1
    rewrite x^2+2*x+5 to (x+1)^2+4
    substitute u for (x+1)^2
    simplify
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 3*x^2/(1-x^4)
    partial fraction decomposition
    apply integral identity
    simplify
done

calculate INT x. 3*x^3/(1-x^4)
    substitute u for x^4
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. tan(x)
    rewrite tan(x) to sin(x)/cos(x)
    substitute u for cos(x)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. cot(x)
    rewrite cot(x) to cos(x)/sin(x)
    substitute u for sin(x)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1 / (x*log(x)*log(log(x))) for x > 1
    substitute u for log(x)
    substitute v for log(u)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1 / sqrt(4-9*x^2)
    rewrite 4-9*x^2 to 4*(1 - (3/2*x)^2)
    rewrite sqrt(4*(1 - (3/2*x)^2)) to 2 * sqrt(1 - (3/2*x)^2)
    substitute u for 3/2*x
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1/sqrt(5-2*x-x^2)
    rewrite 5-2*x-x^2 to 6*(1-((x+1)/sqrt(6))^2)
    rewrite sqrt(6*(1-((x+1)/sqrt(6))^2)) to sqrt(6)*sqrt(1-((x+1)/sqrt(6))^2)
    substitute u for (x+1)/sqrt(6)
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1/sqrt(x*(4-x))
    rewrite x*(4-x) to 4*(1-((x-2)/2)^2)
    rewrite sqrt(4*(1-((x-2)/2)^2)) to 2 * sqrt(1-((x-2)/2)^2)
    substitute u for (x-2)/2
    apply integral identity
    replace substitution
    simplify
done

calculate INT x. 1/(x*sqrt(1-log(x)^2))
    substitute u for log(x)
    apply integral identity
    replace substitution
    simplify
done
