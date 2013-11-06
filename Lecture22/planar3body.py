# pendul - Program to compute the motion of a simple pendulum
# using the Euler or Verlet method
from cpt import *
from math import *
from numpy import array


class Planar3Body :
    m1 = 0.01           # lightest body
    m2 = 0.10           # next heavier body
    m3 = 0.89           # heaviest body
    # use units such that G*(m1+m2+m3) = 4*pi**2
    G = 4 * math.pi**2 / (m1 + m2 + m3)

    def __init__( self, m1, m2, m3 ):
        self.m1 = m1
        self.m2 = m2
        self.m3 = m3
        self.G = 4 * math.pi**2 / (m1 + m2 + m3)

    def __call__(self, trv ) :

        t = trv[0]
        x1 = trv[1] ; y1 = trv[2]  ;  vx1 = trv[3]  ;  vy1 = trv[4]
        x2 = trv[5] ; y2 = trv[6]  ;  vx2 = trv[7]  ;  vy2 = trv[8]
        x3 = trv[9] ; y3 = trv[10] ;  vx3 = trv[11] ;  vy3 = trv[12]

        r12 = math.sqrt( (x1 - x2)**2 + (y1 - y2)**2 )
        r13 = math.sqrt( (x1 - x3)**2 + (y1 - y3)**2 )
        r23 = math.sqrt( (x2 - x3)**2 + (y2 - y3)**2 )

        ax1 = - self.G * self.m2 * (x1 - x2) / r12**3 - self.G * self.m3 * (x1 - x3) / r13**3
        ay1 = - self.G * self.m2 * (y1 - y2) / r12**3 - self.G * self.m3 * (y1 - y3) / r13**3

        ax2 = - self.G * self.m1 * (x2 - x1) / r12**3 - self.G * self.m3 * (x2 - x3) / r23**3
        ay2 = - self.G * self.m1 * (y2 - y1) / r12**3 - self.G * self.m3 * (y2 - y3) / r23**3

        ax3 = - self.G * self.m1 * (x3 - x1) / r13**3 - self.G * self.m2 * (x3 - x2) / r23**3
        ay3 = - self.G * self.m1 * (y3 - y1) / r13**3 - self.G * self.m2 * (y3 - y2) / r23**3

        #print [1.0, vx1, vy1, ax1, ay1, vx2, vy2, ax2, ay2, vx3, vy3, ax3, ay3]

        return [1.0, vx1, vy1, ax1, ay1, vx2, vy2, ax2, ay2, vx3, vy3, ax3, ay3]


def main () : 

    method = input( "Choose a numerical method : RK4 (0), or Adaptive RK4 (1): ")
    m1, m2, m3 = input (" Enter m1, m2, m3: ")
    x1, y1, vx1, vy1 = input(" Enter x1, y1, vx1, vy1: ")
    x2, y2, vx2, vy2 = input(" Enter x2, y2, vx2, vy2: ")

    tau = input(" Enter time step dt = ")
    t_max = input(" Enter total time : ")


    planar3Body = Planar3Body( m1, m2, m3  )

    # compute position and velocity of m3 assuming center of mass at origin
    x3 = - (planar3Body.m1 * x1 + planar3Body.m2 * x2) / planar3Body.m3
    y3 = - (planar3Body.m1 * y1 + planar3Body.m2 * y2) / planar3Body.m3
    vx3 = - (planar3Body.m1 * vx1 + planar3Body.m2 * vx2) / planar3Body.m3
    vy3 = - (planar3Body.m1 * vy1 + planar3Body.m2 * vy2) / planar3Body.m3
    print " x3 = ", x3, " y3 = ", y3, " vx3 = ", vx3, " vy3 = ", vy3

    # compute net angular velocity
    I = planar3Body.m1 * (x1**2 + y1**2) + planar3Body.m2 * (x2**2 + y2**2) + planar3Body.m3 * (x3**2 + y3**3)
    L = planar3Body.m1 * (x1*vy1 - y1*vx1) + planar3Body.m2 * (x2*vy2 - y2*vx2) + planar3Body.m3 * (x3*vy3 - y3*vx3)
    omega = L / I
    print " Net angular velocity = ", omega


    accuracy = 1e-6
    # Initialize vector
    t = 0.0
    xv = [ t, x1, y1, vx1, vy1, x2, y2, vx2, vy2, x3, y3, vx3, vy3 ]


    t_plot = []
    x1_plot = []
    y1_plot = []
    x2_plot = []
    y2_plot = []
    x3_plot = []
    y3_plot = []
    dt_min = tau
    dt_max = tau

    print ''
    print ''
    while xv[0] < t_max :

        #* Record angle and time for plotting
        t_plot.append( xv[0] )
        x1_plot.append( xv[1] )
        y1_plot.append( xv[2] )
        x2_plot.append( xv[5] )
        y2_plot.append( xv[6] )
        x3_plot.append( xv[9] )
        y3_plot.append( xv[10] )
        for ival in xv :
            print '{0:12.6f}'.format( ival),
        print ''
        

        if  method == 0  :
            RK4_step( xv, tau, planar3Body )
        elif method == 1 : 
            tau = RK4_adaptive_step(xv, tau, planar3Body, accuracy)
            dt_min = min(tau, dt_min)
            dt_max = max(tau, dt_max)


    import matplotlib
    matplotlib.rcParams['legend.fancybox'] = True
    import matplotlib.pyplot as plt

    plt.scatter( x1_plot, y1_plot, c='blue', s=1 )
    plt.scatter( x2_plot, y2_plot, c='green', s=11*1 )
    plt.scatter( x3_plot, y3_plot, c='red', s=109*1 )

    plt.show()



if __name__ == "__main__" :
    main()
