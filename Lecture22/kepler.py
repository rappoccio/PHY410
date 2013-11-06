# pendul - Program to compute the motion of a simple pendulum
# using the Euler or Verlet method
from cpt import *
from math import *
from numpy import array


class Kepler :
    G_m1_plus_m2 = 0.
    step_using_y = False
    def __init__( self, G_m1_plus_m2, step_using_y ):
        self.G_m1_plus_m2 = G_m1_plus_m2
        self.step_using_y = step_using_y

    def __call__(self, p ) :

        t = p[0]
        x = p[1]
        y = p[2]
        vx = p[3]
        vy = p[4]
        r = math.sqrt(x**2 + y**2)
        ax = - self.G_m1_plus_m2 * x / r**3
        ay = - self.G_m1_plus_m2 * y / r**3
        flow = [ 1, vx, vy, ax, ay ]
        if self.step_using_y:            # change independent variable from t to y
            for i in range(5):
                flow[i] /= vy
        return flow


def main () : 


    kepler = Kepler( G_m1_plus_m2 = 4 * math.pi**2, step_using_y=False  )

    #* Select the numerical method to use: Euler or Verlet
    method = input( "Choose a numerical method : RK4 (0), or Adaptive RK4 (1): ")
    print " Kepler orbit using fixed and then adaptive Runge-Kutta"
    r_aphelion = float(input(" Enter aphelion distance in AU: "))
    eccentricity = float(input(" Enter eccentricity: "))
    a = r_aphelion / (1 + eccentricity)
    T = a**1.5
    vy0 = math.sqrt(kepler.G_m1_plus_m2 * (2 / r_aphelion - 1 / a))
    print " Semimajor axis a = ", a, " AU"
    print " Period T = ", T, " yr"
    print " v_y(0) = ", vy0, " AU/yr"
    tau = float(input(" Enter step size dt: "))
    accuracy = float(input(" Enter desired accuracy for adaptive integration: "))



    # Initialize vector
    xv = [0.0] * 5
    xv[0] = 0.0
    xv[1] = r_aphelion
    xv[2] = 0.0
    xv[3] = 0.0
    xv[4] = vy0

    #* Loop over desired number of steps with given time step
    #    and numerical method
    nStep = input( "Enter number of time steps: " )

    t_plot = []
    x_plot = []
    y_plot = []
    dt_min = tau
    dt_max = tau
    
    for  iStep in xrange(nStep) :

        #if ( xv[0] > T ) :
        #    break

        #* Record angle and time for plotting
        t_plot.append( xv[0] )
        x_plot.append( xv[1] )
        y_plot.append( xv[2] )

        if  method == 0  :
            RK4_step( xv, tau, kepler )
        elif method == 1 : 
            tau = RK4_adaptive_step(xv, tau, kepler, accuracy)
            dt_min = min(tau, dt_min)
            dt_max = max(tau, dt_max)


    import matplotlib
    matplotlib.rcParams['legend.fancybox'] = True
    import matplotlib.pyplot as plt

    plt.scatter( x_plot, y_plot )

    plt.show()



if __name__ == "__main__" :
    main()
