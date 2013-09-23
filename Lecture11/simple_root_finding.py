#!/usr/bin/env python
from math import exp, log

def f (x) :
    return exp(x) * log(x) - x * x;


print ' Simple search for root of exp(x)*log(x) - x*x'
print ' ---------------------------------------------'
x = float (input(' Enter guess x: '))
dx = float (input(' Enter step dx: ') )
acc = float (input(' Enter accuracy: ' ) )


print " Step            x                    dx"
print " ----    ------------------    ------------------"
step = 0
print  ' {0:4.0f}    {1:20.15f} {2:20.15f}'.format(step, x, dx)
f_old = f(x)

while abs(dx) > abs(acc) :
    x += dx
    if f_old * f(x) < 0 :
        x -= dx
        dx /= 2
    
    step += 1
    print  ' {0:4.0f}    {1:20.15f} {2:20.15f}'.format(step, x, dx)
