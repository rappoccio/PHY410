import math
from cpt import *

# the restricted circular planar 3-body problem has one parameter


# represent a point in the extended phase space by a 5-component vector
# trv = [ t, r, v ] = [ t, x, y, vx, vy ]


class Rcp3Body :
    # alpha = m2/(m1+m2) in the webnotes
    a = 0.0
    
    # switch to zero in on Poincare section point
    # use y instead of t as independent variable
    step_using_y = False
    def __init__(self, a, step_using_y ) :
        self.a = a
        self.step_using_y=step_using_y

    def set_step_using_y( self, step_using_y ) :
        self.step_using_y = step_using_y
    def __call__(self, trv):         # equations in co-rotating frame

        t = trv[0]
        x = trv[1] ; y = trv[2] ; vx = trv[3] ; vy = trv[4]

        d1 = math.pow( (x - self.a)**2 + y**2, 1.5 )
        d2 = math.pow( (x + 1 - self.a)**2 + y**2, 1.5 )

        ax = - (1 - self.a) * (x - self.a) / d1 - self.a * (x + 1 -self.a) / d2 + x + 2 * vy
        ay = - (1 - self.a) * y / d1 - self.a * y / d2 + y - 2 * vx

        flow = [ 1.0, vx, vy, ax, ay ]

        if self.step_using_y:        # change integration variable from t to y
            for i in range(len(flow)):
                flow[i] /= vy

        return flow

    def Jacobi(self, trv):    # Jacobi Integral

        t = trv[0]
        x = trv[1] ; y = trv[2] ; vx = trv[3] ; vy = trv[4]

        r1 = math.sqrt( (x - self.a)**2 + y**2 )
        r2 = math.sqrt( (x + 1 - self.a)**2 + y**2 )

        return x**2 + y**2 + 2 * (1 - self.a) / r1 + 2 * self.a / r2 - vx**2 - vy**2

    def f_x(self, x):     # effective x component of force on the x-axis
        return ( x - (1 - self.a) * (x - self.a) / abs(x - self.a) / (x - self.a)**2
                   - self.a * (x + 1 - self.a) / abs(x + 1 - self.a) / (x + 1 - self.a)**2 )

    def zero(self, f, x_lower, x_upper, accuracy=1.0e-6, max_steps=1000):
        # use bisection search to solve f(x) = 0 in interval [x_lower, x_upper]
        assert x_lower < x_upper, " zero(f, x_lower, x_upper) not bracketed"
        x_mid = (x_upper + x_lower) / 2
        dx = x_upper - x_lower
        f_lower = f(x_lower)
        step = 0
        while abs(dx) > accuracy:
            f_mid = f(x_mid)
            if f_mid == 0:
                dx = 0
            else:
                if f_lower * f_mid > 0:
                    x_lower = x_mid
                    f_lower = f_mid
                else:
                    x_upper = x_mid
                x_mid = (x_upper + x_lower) / 2
                dx = x_upper - x_lower
            step += 1
            assert step < max_steps, " zero(f) too many steps " + str(max_steps)
        return x_mid

# get parameters from user
print " Restricted circular planar 3-body problem"

while True:
    a = input(" Enter alpha = m2/(m1+m2) > 0 and < 0.5: ")
    if a <= 0 or a > 0.5:
        print " Bad alpha, please try again"
        continue
    break

rcp3Body = Rcp3Body (a, False)

# find Lagrangian points for this alpha
eps = 1e-6
print " Lagrangian points:"
print " L1 : x = ", rcp3Body.zero(rcp3Body.f_x, rcp3Body.a - 1 + eps, rcp3Body.a - eps), " y = ", 0.0
print " L2 : x = ", rcp3Body.zero(rcp3Body.f_x, -1.5,rcp3Body. a - 1 - eps), " y = ", 0.0
print " L3 : x = ", rcp3Body.zero(rcp3Body.f_x, rcp3Body.a + eps, 1.5), " y = ", 0.0
print " L4 : x = ", rcp3Body.a - 0.5, " y = ", math.sqrt(3.0) / 2
print " L5 : x = ", rcp3Body.a - 0.5, " y = ", math.sqrt(3.0) / 2

# get initial values
one = input(" Enter 0 to specify [x,y,vx,vy] or 1 to specify C and [x,y,vx]: ")
if one:
    while True:
        C = input(" Enter value of the Jacobi integral C: ")
        x, y, vx = input(" Enter x, y, vx: ")
        r1 = math.sqrt((x - a)**2 + y**2)
        r2 = math.sqrt((x + 1 - rcp3Body.a)**2 + y**2)
        vy_sqd = - C + x**2 + y**2 + 2 * (1 - rcp3Body.a) / r1 + 2 * rcp3Body.a / r2 - vx**2
        if vy_sqd < 0:
            print " Sorry C too large, cannot solve for vy"
        else:
            vy = math.sqrt(vy_sqd)
            print " vy = ", vy
            break
else:
    x, y, vx, vy = input(" Enter x, y, vx, vy: ")
    print " Jacobi Integral C = ", rcp3Body.Jacobi([0, x, y, vx, vy])

t_max = input(" Enter maximum integration time t_max: ")


crossing = 0
dt = 0.01
t = 0
trv = [ t, x, y, vx, vy ]

x_plot = []
y_plot = []

while t < t_max:


    # write trajectory point
    t = trv[0]
    x = trv[1] ; y = trv[2] ; vx = trv[3] ; vy = trv[4]
    y_save = y      # remember y to check section crossing
    x_plot.append( x )
    y_plot.append( y )

    # use adaptive Runge-Kutta with default accuracy
    dt = RK4_adaptive_step(trv, dt, rcp3Body)

    # Poincare section at y = 0 and vy positive
    x = trv[1] ; y = trv[2] ; vx = trv[3] ; vy = trv[4]
    if y_save < 0 and y >= 0 and vy >= 0:
        rcp3Body.set_step_using_y(True)
        dy = -y
        RK4_step(trv, dy, rcp3Body)
        t = trv[0] ; x = trv[1] ; y = trv[2] ; vx = trv[3] ; vy = trv[4]
        crossing += 1
        print " Crossing No.", crossing, " at t = ", t, " C = ", rcp3Body.Jacobi(trv)
        rcp3Body.set_step_using_y(False)



import matplotlib
matplotlib.rcParams['legend.fancybox'] = True
import matplotlib.pyplot as plt

plt.scatter( x_plot, y_plot, c='blue' )

plt.show()
