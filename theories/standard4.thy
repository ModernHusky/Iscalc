# Indefinite Integrals of Inverse Trigonometric Functions
# Handbook of mathematical formulas and integrals
# Page 225

prove (INT x. arcsin(x/a)) = arcsin(x/a) + sqrt(a ^ 2 - x ^ 2) + SKOLEM_CONST(C) for abs(x/a)<=1
sorry

prove (INT x. arcsin(x/a) ^ 2) = x * arcsin(x/a) ^ 2 + 2 * sqrt(a ^ 2 - x ^ 2) * arcsin(x/a) - 2 * x + SKOLEM_CONST(C) for abs(x/a)<=1
sorry

prove (INT x. arcsin(x/a) ^ 3) = x * arcsin(x/a) ^ 3 + 3 * sqrt(a ^ 2 - x ^ 2) * arcsin(x/a) ^ 2 - 6 * x * arcsin(x/a) - 6 * sqrt(a ^ 2 - x ^ 2) + SKOLEM_CONST(C) for abs(x/a)<=1
sorry

prove (INT x. arcsin(x/a) ^ n) = x * arcsin(x/a) ^ n + n * sqrt(a ^ 2 - x ^ 2) * arcsin(x/a) ^ (n -1) - n * (n -1) * (INT x. arcsin(x/a) ^ (n -2)) + SKOLEM_CONST(C) for abs(x/a)<=1, n !=1
sorry

prove (INT x. x * arcsin(x/a)) = (x ^ 2/2 -  a ^ 2/4) * arcsin(x/a) + x/4 * sqrt(a ^ 2 - x ^ 2) + SKOLEM_CONST(C) for abs(x/a)<=1
sorry

prove (INT x. x ^ 2 * arcsin(x/a)) = x ^ 3/3 * arcsin(x/a) + 1/9 * (x ^ 2 + 2 * a ^ 2) * sqrt(a ^ 2 - x ^ 2) + SKOLEM_CONST(C) for abs(x/a)<=1
sorry

prove (INT x. x ^ 3 * arcsin(x/a)) = (x ^ 4/4 -  3 * a ^ 4/32) * arcsin(x/a) + 1/32 * (2 * x ^ 3 + 3 * a ^ 2 * x) * sqrt(a ^ 2 - x ^ 2) + SKOLEM_CONST(C) for abs(x/a)<=1
sorry

prove (INT x. x ^ n * arcsin(x/a)) = x ^ (n + 1)/(n + 1) * arcsin(x/a) - 1/(n + 1) * (INT x. x ^ (n + 1)/sqrt(a ^ 2 - x ^ 2)) + SKOLEM_CONST(C) for abs(x/a)<=1
sorry

prove (INT x. 1/x ^ 2 * arcsin(x/a)) = -1/x * arcsin(x/a) - 1/a * log(abs((a + sqrt(a ^ 2 - x ^ 2))/x)) + SKOLEM_CONST(C) for abs(x/a)<1
sorry

prove (INT x. 1/x ^ 3 * arcsin(x/a)) = -1/(2 * x ^ 2 ) * arcsin(x/a) - sqrt(a ^ 2 - x ^ 2)/(2 * a ^ 2 * x) + SKOLEM_CONST(C) for abs(x/a)<1
sorry

prove (INT x. 1/x ^ n * arcsin(x/a)) = -1/((n - 1) * x ^ (n - 1)) * arcsin(x/a) + 1/(n - 1) * (INT x. 1/(x ^ (n - 1) * sqrt(a ^ 2 - x ^ 2))) + SKOLEM_CONST(C) for abs(x/a)<1, n !=1
sorry

prove (INT x. arccos(x/a)) = x * arccos(x/a) - sqrt(a ^ 2 - x ^ 2) + SKOLEM_CONST(C) for abs(x/a)<=1
sorry

prove (INT x. arccos(x/a) ^ 2) = x * arccos(x/a) ^ 2 - 2 * sqrt(a ^ 2 - x ^ 2) * arccos(x/a) - 2 * x + SKOLEM_CONST(C) for abs(x/a)<=1
sorry

prove (INT x. arccos(x/a) ^ 3) = x * arccos(x/a) ^ 3 - 3 * sqrt(a ^ 2 - x ^ 2) * arccos(x/a) ^ 2 - 6 * x * arccos(x/a) + 6 * sqrt(a ^ 2 - x ^ 2) + SKOLEM_CONST(C) for abs(x/a)<=1
sorry

prove (INT x. arccos(x/a) ^ n) = x * arccos(x/a) ^ n -  n * sqrt(a ^ 2 - x ^ 2) * arccos(x/a) ^ (n - 1) - n * (n - 1) * (INT x. arccos(x/a) ^ (n - 2)) + SKOLEM_CONST(C) for abs(x/a)<=1, n !=1
sorry

prove (INT x. x * arccos(x/a)) = (x ^ 2/2 -  a ^ 2/4) * arccos(x/a) + x/4 * sqrt(a ^ 2 - x ^ 2) + SKOLEM_CONST(C) for abs(x/a)<1
sorry

prove (INT x. x ^ 2 * arccos(x/a)) = x ^ 3/3 * arccos(x/a) - 1/9 * (x ^ 2 + 2 * a ^ 2) * sqrt(a ^ 2 - x ^ 2) + SKOLEM_CONST(C) for abs(x/a)<1
sorry

prove (INT x. x ^ 3 * arccos(x/a)) = (x ^ 4/4 -  3 * a ^ 4/32) * arccos(x/a) - 1/32 * (2 * x ^ 3 + 3 * a ^ 2 * x) * sqrt(a ^ 2 - x ^ 2) + SKOLEM_CONST(C) for abs(x/a)<=1
sorry

prove (INT x. x ^ n * arccos(x/a)) = x ^ (n + 1)/(n + 1) * arccos(x/a) +1/(n + 1) * INT x. x ^ ( n + 1/sqrt(a ^ 2 - x ^ 2)) + SKOLEM_CONST(C) for abs(x/a)<=1, n !=1
sorry

prove (INT x. 1/x ^ 2 * arccos(x/a)) = -1/x * arccos(x/a) + a * log(abs((a + sqrt(a ^ 2 - x ^ 2))/x)) + SKOLEM_CONST(C) for abs(x/a)<=1
sorry

prove (INT x. 1/x ^ 3 * arccos(x/a)) = -1/(2 * x ^ 2) * arccos(x/a) + sqrt(a ^ 2 - x ^ 2)/(2 * a ^ 2 * x) + SKOLEM_CONST(C) for abs(x/a)<=1
sorry

prove (INT x. 1/x ^ n * arccos(x/a)) = -1/((n - 1) * x ^ (n - 1)) * arccos(x/a) - 1/(n - 1) * (INT x. 1/(x ^ (n - 1) * sqrt(a ^ 2 - x ^ 2))) + SKOLEM_CONST(C) for abs(x/a)<=1, n !=1
sorry

prove (INT x. arctan(x/a)) = x * arctan(x/a) - a/2 * log(x ^ 2 + a ^ 2) + SKOLEM_CONST(C)
sorry

prove (INT x. x * arctan(x/a)) = 1/2 * (x ^ 2 + a ^ 2) * arctan(x/a) - 1/2 * a * x + SKOLEM_CONST(C)
sorry

prove (INT x. x ^ 2 * arctan(x/a)) = 1/3 * x ^3 * arctan(x/a) - 1/6 * a * x ^ 2 + 1/6 * a ^ 3 * log(x ^ 2 + a ^ 2) + SKOLEM_CONST(C)
sorry

prove (INT x. x ^ 3 * arctan(x/a)) = 1/4 * (x ^ 4 - a ^ 4) * arctan(x/a) - 1/12 * a * x ^ 3 + 1/4 * a ^ 3 * x + SKOLEM_CONST(C)
sorry

prove (INT x. x ^ n * arctan(x/a)) = x ^ (n + 1)/(n + 1) * arctan(x/a) - a/(n + 1) * (INT x. x ^ (n + 1)/(x ^ 2 + a ^ 2)) + SKOLEM_CONST(C)
sorry

prove (INT x. 1/x * arctan(x/a)) = SUM(k, 0, oo, (-1) ^ k/(2 * k+1) ^ 2 * (x/a) ^ (2 * k+1)) + SKOLEM_CONST(C) for abs(x/a)<1
sorry

prove (INT x. 1/x * arctan(x/a)) = pi/2 * log(abs(x)) + SUM(k, 0, oo, (-1) ^ k/(2 * k+1) ^ 2 * (x/a) ^ (2 * k+1)) + SKOLEM_CONST(C) for x/a>1
sorry

prove (INT x. 1/x * arctan(x/a)) = -pi/2 * log(abs(x)) + SUM(k, 0, oo, (-1) ^ k/(2 * k+1) ^ 2 * (x/a) ^ (2 * k+1)) + SKOLEM_CONST(C) for x/a<-1
sorry

prove (INT x. 1/x ^ 2 * arctan(x/a)) = -1/x * arctan(x/a) + 1/(2 * a) * log(x ^ 2/(x ^ 2 + a ^ 2)) + SKOLEM_CONST(C)
sorry

prove (INT x. 1/x ^ 3 * arctan(x/a)) = -1/2 * (1/x ^2 +1/a ^2) * arctan(x/a) - 1/(2 * a * x) + SKOLEM_CONST(C)
sorry

prove (INT x. 1/x ^ n * arctan(x/a)) = -1/((n - 1) * x ^ (n - 1)) * arctan(x/a) + a/(n -1) * prove (INT x. 1/(x ^ (n - 1) * (x ^ 2 + a ^ 2))) + SKOLEM_CONST(C) for n !=1
sorry

prove (INT x. arccot(x/a)) = x * arccot(x/a) + 1/2 * a * log(x ^ 2 + a ^ 2) + SKOLEM_CONST(C)
sorry

prove (INT x. x * arccot(x/a)) = 1/2 * (x ^ 2 + a ^ 2) * arccot(x/a) + 1/2 * a * x + SKOLEM_CONST(C)
sorry

prove (INT x. x ^ 2 * arccot(x/a)) = 1/3 * x ^ 3 * arccot(x/a) + 1/6 * a * x ^ 2 - 1/6 * a ^ 3 * log(x ^ 2 + a ^ 2) + SKOLEM_CONST(C)
sorry

prove (INT x. x ^ 3 * arccot(x/a)) = 1/4 * (x ^ 4 - a ^ 4) * arccot(x/a) + 1/12 * a * x ^ 3 - 1/4 * a ^ 3 * x + SKOLEM_CONST(C)
sorry

prove (INT x. x ^ n * arccot(x/a)) = x ^ (n + 1)/(n + 1) * arccot(x/a) + a/(n + 1) * (INT x. x ^ (n - 1)/(x ^ 2 + a ^ 2)) + SKOLEM_CONST(C)
sorry

prove (INT x. 1/x * arccot(x/a)) = pi/2 * log(abs(x)) - SUM(k, 0, oo, (-1) ^ k/(2 * k+1) ^ 2 * (x/a) ^ (2 * k+1)) + SKOLEM_CONST(C) for abs(x/a) < 1
sorry

prove (INT x. 1/x * arccot(x/a)) = -SUM(k, 0, oo, (-1) ^ k/(2 * k+1) ^ 2 * (x/a) ^ (2 * k+1)) + SKOLEM_CONST(C) for x/a > 1
sorry

prove (INT x. 1/x * arccot(x/a)) = pi * log(abs(x))-SUM(k, 0, oo, (-1) ^ k/(2 * k+1) ^ 2 * (x/a) ^ (2 * k+1)) + SKOLEM_CONST(C) for x/a<-1
sorry

prove (INT x. 1/x ^ 2 * arccot(x/a)) = -1/x * arccot(x/a) - 1/(2 * a) * log(x ^2/(x ^ 2 + a ^ 2)) + SKOLEM_CONST(C)
sorry

prove (INT x. 1/x ^ 3 * arccot(x/a)) = -1/2 * (1/x ^ 2 + 1/a ^ 2) * arccot(x/a) + 1/(2 * a * x) + SKOLEM_CONST(C)
sorry

prove (INT x. 1/x ^ n * arccot(x/a)) = -1/((n - 1) * x ^ (n - 1)) * arccot(x/a) - a/(n -1) * prove (INT x. 1/(x ^ (n - 1) * (x ^ 2 + a ^ 2))) + SKOLEM_CONST(C) for n !=1
sorry

prove (INT x. 1/(x ^ 2 + a ^ 2) * arccot(x/a)) = -1/(2 * a) * arccot(x/a) ^ 2 + SKOLEM_CONST(C)
sorry

prove (INT x. x ^ 2/(x ^ 2 + a ^ 2) * arccot(x/a)) = x * arccot(x/a) +1/2 * a * log(x ^ 2 + a ^ 2) + 1/2 * a * arccot(x/a) ^ 2 + SKOLEM_CONST(C)
sorry

prove (INT x. 1/(x ^ 2 + a ^ 2) ^ 2 * arccot(x/a)) = x / (2 * a ^ 2 * (x ^ 2 + a ^ 2)) * arccot(x/a) - 1/(4 * a ^ 3) * arccot(x/a) ^ 2 -1/(4 * a * (x ^ 2 + a ^ 2)) + SKOLEM_CONST(C)
sorry

prove (INT x. 1/(x ^ 2 + a ^ 2) * arccot(x/a) ^ n) = -1/((n + 1) * a) * arccot(x/a) ^ (n + 1) + SKOLEM_CONST(C)
sorry