import sys
import math
import numpy as np
import random

n = 15
t = 0.
m = 1.

print n
print t

for i in xrange(n):
    print ' %6d %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e' % ( m, random.uniform(-1.0,1.0), random.uniform(-1.0,1.0), 0, random.uniform(-1.0,1.0), random.uniform(-1.0,1.0), 0 )
