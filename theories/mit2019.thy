imports standard

calculate INT x:[0,pi / 100]. (sin(20 * x) + sin(19 * x)) / (cos(20 * x) + cos(19 * x))
    rewrite sin(20 * x) + sin(19 * x) to 2 * cos(1/2 * x) * sin(39/2 * x)
    rewrite cos(20 * x) + cos(19 * x) to 2 * cos(1/2 * x) * cos(39/2 * x)
    simplify
    substitute u for cos(39/2 * x)
    apply integral identity
    simplify
done
