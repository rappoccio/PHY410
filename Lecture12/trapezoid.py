#!/usr/bin/env python

def trapezoid(f, a, b, n):
    """Approximates the definite integral of f from a to b by
    the composite trapezoidal rule, using n subintervals.
    From http://en.wikipedia.org/wiki/Trapezoidal_rule
    """
    h = (b - a) / n
    s = f(a) + f(b)
    for i in xrange(1, n):
        s += 2 * f(a + i * h)
    return s * h / 2


def adaptive_trapezoid(f, a, b, acc, output=False):
    """
    Uses the adaptive trapezoidal method to compute the definite integral
    of f from a to b to desired accuracy acc. 
    """
    
    old_s = -1e-30
    h = b - a
    n = 1
    s = (f(a) + f(b)) / 2
    if output == True : 
        print "N = " + str(n+1) + ",  Integral = " + str( h*s )
    while abs(h * (old_s - s/2)) > acc :
        old_s = s
        for i in xrange(n) :
            s += f(a + (i + 0.5) * h)
        n *= 2
        h /= 2
        if output == True :
            print "N = " + str(n+1) + ",  Integral = " + str( h*s )
    return h * s

