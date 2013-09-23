from root_finding import *
from math import *

def f ( x ) :
    return exp(x) * log(x) - x * x


print(" Algorithms for root of exp(x)*log(x) - x*x")
print(" ------------------------------------------------")

print(" 1. Simple search")
x0 = float ( input(" Enter initial guess x_0 : ") )
dx = float ( input(" Enter step dx : ") )
acc = float ( input(" Enter accuracy : ") )
answer = root_simple(f, x0, dx, acc,1000,True)
print  str ( answer ) + "\n\n"

print(" 2. Bisection search")
x0 = float ( input(" Enter bracketing guess x_0 : ") )
x1 = float ( input(" Enter bracketing guess x_1 : ") )
acc = float ( input(" Enter accuracy : ") )
answer = root_bisection(f, x0, x1, acc,1000,True)
print  str ( answer ) + "\n\n"


