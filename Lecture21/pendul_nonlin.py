# pendul - Program to compute the motion of a simple pendulum
# using the Euler or Verlet method
from cpt import *
from math import *

class PendulumFunction :
    L=0.                           # length in meters
    omega_0=0.                     # natural frequency
    gamma=0.                       # damping constant
    F_D=0.                         # driving force amplitude
    omega_D=0.                     # driving force frequency  
    g = 9.8                        # gravitational constant
    def __init__( self, L=g, gamma=0, F_D=0, omega_D=0 ):
        self.L=L
        self.omega_0 = sqrt(self.L / self.g)
        self.gamma=gamma
        self.F_D=F_D
        self.omega_D=omega_D

class PendulumAcceleration( PendulumFunction ):
    def __init__( self, L=PendulumFunction.g, gamma=0, F_D=0, omega_D=0 ):
        PendulumFunction.__init__( self, L, gamma, F_D, omega_D)

    def __call__( self, p ) :
        t=p[0]
        theta=p[1]
        dtheta_dt=p[2]
        a = - self.omega_0 * self.omega_0 * sin(theta)     # due to gravity
        if self.gamma > 0.0 :
            a += - self.gamma * dtheta_dt                        # damping
        if self.F_D > 0.0 :
            a += self.F_D * cos(self.omega_D * t)                     # driving force
        return a


class PendulumDiffEq( PendulumFunction ):
    def __init__( self, L=PendulumFunction.g, gamma=0, F_D=0, omega_D=0 ):
        PendulumFunction.__init__(self, L, gamma, F_D, omega_D)

    def __call__(self, p ) :
        t = float(p[0])
        x = float(p[1])
        y = float(p[2])
        out =[0.]* 3
        out[0] = 1.0
        out[1] = y
        out[2] = -self.omega_0*self.omega_0*sin(x) - self.gamma*y + self.F_D*cos(self.omega_D*t)
        #print 'p is '
        #print p
        #print ' out is '
        #print out
        #junk = input("enter junk : ")
        return out


def main () : 

    #* Select the numerical method to use: Euler or Verlet
    method = input( "Choose a numerical method : Euler (0), RK4 (1), or Adaptive RK4 (2): ")
					   
    #* Set initial position and velocity of pendulum
    theta0 = input( "Enter initial angle (in degrees): ")
    theta = theta0*pi/180   # Convert angle to radians
    omega = input("Enter initial angular frequency (degrees/s): ")
    omega = omega*pi/180.0
    gamma = input( "Enter damping coefficient : " )
    F_D = input( "Enter driving amplitude : ")
    omega_D = input( "Enter driving frequency : " )

    #* Set the physical constants and other variables
    L = PendulumDiffEq.g  # The constant g/L
    time = 0.0            # Initial time
    time_old = 0.0        # Time of previous reversal
    irev = 0               # Used to count number of reversals
    tau = input( "Enter time step: ")
    
    pendulumAcc = PendulumAcceleration ( L, gamma, F_D, omega_D)

    pendulumDiffEq = PendulumDiffEq( L, gamma, F_D, omega_D)



    T_p = 0.0               # Period for Poincare section

    if omega_D > 0.0 :
        T_D = 2*pi/omega_D      # Period of driving force
        T_p = T_D
    else :
        T_p = 2*pi/sqrt(PendulumDiffEq.g / L)


    accuracy = 1e-6

    #* Take one backward step to start Verlet

    xv = [0.0] * 3
    xv[0] = time
    xv[1] = theta
    xv[2] = omega
    accel = pendulumAcc( xv )
    theta_old = theta - omega*tau + 0.5*tau*tau*accel    


    xv[0] = time
    xv[1] = theta_old
    xv[2] = omega

    iperiod = 0
    iperiod_old = -1

    #* Loop over desired number of steps with given time step
    #    and numerical method
    nStep = input( "Enter number of time steps: " )


    plotTrajectory = True # plot Poincare and trajectory if true.
                          # plot only Poincare if false.

    t_plot = []
    th_plot = []
    om_plot = []
    period = []
    poincare_map = []
    dt_min = tau
    dt_max = tau
    for  iStep in xrange(nStep) :
        theta_old = xv[1]

        #* Record angle and time for plotting
        if plotTrajectory  : 
            t_plot.append( xv[0] )

            thplot = xv[1]
            omplot = xv[2]

            while thplot > pi :
                thplot -= 2*pi
            while thplot <= -pi:
                thplot += 2*pi

            th_plot.append( thplot*180/pi )   # Convert angle to degrees
            om_plot.append( omplot*180/pi )

        iperiod = int(  float(xv[0]) / float(T_p) )
        # Calculate a Poincare map of the dynamics
        if  iperiod != iperiod_old :
            iperiod_old = iperiod

            thplot = xv[1]
            omplot = xv[2]

            while thplot > pi :
                thplot -= 2*pi
            while thplot <= -pi:
                thplot += 2*pi


            poincare_map.append( [thplot*180/pi, omplot*180/pi] )
            print "Period reached. theta = "  + str(  thplot*180/pi ) + ", omega = " + str(  omplot*180/pi )



        if  method == 0 :
            Euler_step( xv, tau, pendulumAcc )
        elif  method == 1  :
            RK4_step( xv, tau, pendulumDiffEq )
        elif method == 2 : 
            tau = RK4_adaptive_step(xv, tau, pendulumDiffEq, accuracy)
            dt_min = min(tau, dt_min)
            dt_max = max(tau, dt_max)


        #* Test if the pendulum has passed through theta = 0
        #    if yes, use time to estimate period
        if xv[1]*theta_old < 0 : # Test position for sign change
            if irev == 0:           # If this is the first change,
                time_old = xv[0]       # just record the time
            else  :
                period.append( 2*(xv[0] - time_old) )
                time_old = xv[0]
            irev += 1       # Increment the number of reversals
        nPeriod = irev-1    # Number of times period is measured

    #* Estimate period of oscillation, including error bar
    AvePeriod = 0.0
    ErrorBar = 0.0
    for i in xrange( nPeriod ) : 
        AvePeriod += period[i]


    AvePeriod /= nPeriod
    for i in xrange( nPeriod ) :
        ErrorBar += (period[i] - AvePeriod)*(period[i] - AvePeriod)
    ErrorBar = sqrt(ErrorBar/(nPeriod*(nPeriod-1)))
    print "Average period = " + str( AvePeriod ) + " +/- " + str( ErrorBar )

    #* Print out the plotting variables: t_plot, th_plot
    if plotTrajectory :
        plotOut = file("pend_plot.txt", 'w')
        for i in xrange( len(t_plot) ) :
            s = str( t_plot[i] ) + " " + str ( th_plot[i] ) +  " "  + str( om_plot[i]) + '\n'
            plotOut.write( s ) 

        plotOut.close(  )

    poincareOut = file("poincare_plot.txt", 'w')
    xi = [ poincare_map[i][0] for i in xrange(len(poincare_map)) ]
    yi = [ poincare_map[i][1] for i in xrange(len(poincare_map)) ]
    for i in xrange( len(xi) ) : 
        s = str( xi[i] ) + ' ' + str( yi[i] ) + '\n'
        poincareOut.write( s ) 
    poincareOut.close (  )

    import matplotlib
    matplotlib.rcParams['legend.fancybox'] = True
    import matplotlib.pyplot as plt
    from matplotlib import legend


    plt.scatter( xi, yi )
    plt.ylim( [-180,180] )
    plt.xlim( [-180,180] )

    plt.show()



if __name__ == "__main__" :
    main()
