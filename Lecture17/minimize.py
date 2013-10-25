from cpt import *
from math import *
import numpy as np


def f (p) :
    x = p[0]
    y = p[1]
    return x*x + 2 * y*y + 0.3 * cos(3 * pi * x) + 0.4 * cos(4 * pi * y)

def df(p) :
    x = p[0]
    y = p[1]
    x = 2 * x - 0.9 * pi * sin(3 * pi * x)
    y = 4 * y - 1.6 * pi * sin(4 * pi * y)
    return np.array( [x,y] )


print " Minimization using Broyden-Fletcher-Goldfarb-Shanno Algorithm"
print " Find minimum of f(x,y) given an initial guess for x, y"
p = input(" Enter starting point coordinates x y: ")
gtol = input( " Enter desired accuracy: ")
f_min = 0.0
iterations = 0
res = scipy.optimize.fmin_bfgs(f=f, fprime=df,x0=p, gtol=gtol)
print res
