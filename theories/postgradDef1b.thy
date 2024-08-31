// page 210

calculate INT x:[0, 2]. x^2*sqrt(4-x^2)
    substitute 2*sin(u) for x
    rewrite 4-(2*sin(u))^2 to 4*(1-sin(u)^2)
    rewrite 1-sin(u)^2 to cos(u)^2
    simplify
    rewrite cos(u) ^ 2 * sin(u) ^ 2 to (sin(u)*cos(u))^2
    rewrite sin(u)*cos(u) to 1/2*sin(2*u)
    simplify
    rewrite sin(2*u)^2 to 1/2*(1-cos(4*u))
    apply integral identity
    simplify
done

calculate INT x:[0, 1]. x*arctan(x)
    integrate by parts with u=arctan(x), v=1/2*x^2
    simplify
    partial fraction decomposition
    apply integral identity
    simplify
done

calculate INT x:[1, exp(1)]. sin(log(x))
    substitute u for log(x)
    integrate by parts with u=exp(u), v=-cos(u)
    simplify
    integrate by parts with u=exp(u), v=sin(u)
    simplify
    solve integral INT u:[0,1]. exp(u)*sin(u)
    simplify
done

// 2*arcsin(sqrt(2/3))-pi/2 is the answer in the book
calculate INT x:[1/2, 2/3]. 1/sqrt(x*(1-x))
    rewrite x*(1-x) to 1/4 - (x-1/2)^2
    substitute u for x-1/2
    apply integral identity
    simplify
done

// could further optimize the simplification
// exp(-(2 * log(2))) = exp(-2*log(2)) = exp(log(2))^(-2) = 2 ^ (-2) = 1/4
calculate INT x:[0, -log(2)]. sqrt(1-exp(2*x))
    substitute 1/2*log(1-t^2) for x
    simplify
    partial fraction decomposition
    apply integral identity
    simplify
done