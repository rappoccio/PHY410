from math import sqrt, sin, cos, exp
from numpy import array, nan_to_num, inf
 
def ridder(f,x,h):
    ''' f     : name of function to be differentiated
        x     : the point at which df/dx is required
        h     : suggestion for an initial step size
        error : estimate of error by algorithm
    '''

    if h == 0.0 :
        print "diff_Ridders: h must be non-zero"
        exit

    n = 10           # dimension of extrapolation table
    a = array(  [[0.0] * n] *n )             # extrapolation table

    a[0][0] = (f(x + h) - f(x - h)) / (2 * h)
    answer = 0.0
    error = nan_to_num( inf ) / 2.0  # get a large value for the error
    for i in xrange(n) : 
        h /= 1.4
        a[0][i] = (f(x + h) - f(x - h)) / (2 * h)
        fac = 1.4 * 1.4
        for j in range(1,i+1) :
            a[j][i]=(a[j-1][i] * fac - a[j-1][i-1]) / (fac - 1)
            fac *= 1.4 * 1.4
            err = max(abs(a[j][i] - a[j-1][i]),
                      abs(a[j][i] - a[j-1][i-1]))
            if err <= error :
                error = err
                answer = a[j][i]

        if abs(a[i][i] - a[i-1][i-1]) >= 2 * error :
            break
    return answer,error


def test_ridders_x2 ( x ) :
    return x*x

print "Testing Ridder's algorithm for computing numeric derivatives"

error = 0.0
h = 0.01

for x in [2.0,3.0,4.0,5.0] :
    xprime,error =  ridder( exp, x, h)
    print 'exp({0:1.0f}) : True value = {1:38.32f}, Ridder = {2:38.32f}, exp error = {3:38.32f}, obs error = {4:38.32f}'.format( x, exp(x), xprime, h**4, error)
    xprime,error =  ridder( sin, x, h)
    print 'sin({0:1.0f}) : True value = {1:38.32f}, Ridder = {2:38.32f}, exp error = {3:38.32f}, obs error = {4:38.32f}'.format( x, cos(x), xprime, h**4, error)
    xprime,error =  ridder( test_ridders_x2, x, h)
    print 'x*2({0:1.0f}) : True value = {1:38.32f}, Ridder = {2:38.32f}, exp error = {3:38.32f}, obs error = {4:38.32f}'.format( x, 2*x, xprime, h**4, error)
    
