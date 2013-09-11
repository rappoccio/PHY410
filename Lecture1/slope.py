from math import *


# "Dumb" slope calculator : No error handling!
def getSlope0 ( x0,y0, x1, y1) : 
    slope = (y1-y0) / (x1-x0)
    return slope

# Better slope calculator : Error handling when
# denominator is zero
def getSlope ( x0,y0, x1, y1) :
    num = y1-y0
    den = x1-x0
    if abs(den) > 0.00001 :
        slope = num / den
    else :
        print '----> invalid input, setting to None'
        slope = None
    return slope


# Unit tests!

print 'Unit test 1'
x0 = 0.0
y0 = 0.0
x1 = 1.0
y1 = 1.0
slope = getSlope( x0,y0,x1,y1 )
print slope


print 'Unit test 2'
x0 = 0.0
y0 = 0.0
x1 = 0.0
y1 = 1.0
slope = getSlope( x0,y0,x1,y1 )
print slope


