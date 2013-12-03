import math
import random

rules = [ "Rectangle", "Midpoint", "Trapezoid", "Simpson's",
          "Monte Carlo" ]

def regular(f, a, b, rule):
    x_mid = (a + b) / 2.0
    if rule == "Rectangle":
        return f(a) * (b - a)
    elif rule == "Midpoint":
        return f(x_mid) * (b - a)
    elif rule == "Trapezoid":
        return (f(a) + f(b)) * (b - a) / 2.0
    elif rule == "Simpson's":
        return (f(a) + 4*f(x_mid) + f(b)) * (b - a) / 6.0
    else:
        raise ValueError(rule + ' not implemented')

def integral(f, a, b, rule, n_max=10**7, max_error = 1e-6):
    n = 1
    while True:
        if n > n_max:
            print " N exceeds", n_max, "aborting ..."
            break
        d = (b - a) / float(n)
        if rule == "Monte Carlo":
            approx = sum(f(random.uniform(a, b)) for i in range(n))
            approx *= (b - a) / float(n)
        else:
            approx = sum(regular(f, a + i*d, a + (i+1)*d, rule)
                         for i in range(n))
        error = approx - exact
        print " " + repr(n) + "\t" + repr(error)
        if abs(error) < max_error:
            break
        n *= 2
    return approx

def f(x):
    return 1 / (1.0 + x**2)

a = 0.0
b = 5.0
exact = math.atan(b) - math.atan(a)

print " Quadrature of 1 / (1 + x**2) from", a, "to", b
print " Exact Integral =", exact
print " ------------------------------"

for rule in rules:
    print " Quadrature rule:", rule
    print " N    \tError"
    answer =  integral(f, a, b, rule)
    print " Integral = " + repr(answer)
    print " ------------------------------"
