import math
import cmath
import time
import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from cpt import RK4_step, RK4_adaptive_step
import copy
class RKHelper:
    def __init__(self, flow):
        self.flow = flow

    def __call__(self, i) :
        return self.flow

    
    
## Modified from https://stackoverflow.com/questions/9401658/how-to-animate-a-scatter-plot
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

class AnimatedScatter(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""
    def __init__(self, numpoints=50, dt=0.1, G=1.):
        self.numpoints = numpoints
        self.dt = dt
        self.G=G
        self.stream = self.data_stream()

        # Setup the figure and axes...
        self.fig, self.ax = plt.subplots()
        # Then setup FuncAnimation.
        self.ani = animation.FuncAnimation(self.fig, self.update, interval=5, 
                                           init_func=self.setup_plot, blit=True)

    def setup_plot(self):
        """Initial drawing of the scatter plot."""
        x, y, vx, vy, ax, ay = next(self.stream)
        self.scat = self.ax.scatter(x, y, animated=True)
        self.ax.axis([-10, 10, -10, 10])

        # For FuncAnimation's sake, we need to return the artist we'll be using
        # Note that it expects a sequence of artists, thus the trailing comma.
        return self.scat,

    def data_stream(self):
        """Generate a random walk (brownian motion). Data is scaled to produce
        a soft "flickering" effect."""
        data = np.random.random((6, self.numpoints))
        xy = data[:2, :]
        vxy = data[2:4, :]
        axy = data[4:, :]
        xy -= 0.5
        xy *= 10
        vxy -= 0.5
        axy -= 0.5
        #axy *= 0.001
        t = 0.0
        while True:
            for ipoint in xrange( self.numpoints ) :
                xi  = xy[0][ipoint]
                yi  = xy[1][ipoint]
                vxi = vxy[0][ipoint]
                vyi = vxy[1][ipoint]
                axi = axy[0][ipoint] = 0.0
                ayi = axy[1][ipoint] = 0.0
                for jpoint in xrange ( self.numpoints ) :
                    if ipoint == jpoint :
                        continue
                    xj  = xy[0][jpoint]
                    yj  = xy[1][jpoint]
                    vxj = vxy[0][jpoint]
                    vyj = vxy[1][jpoint]
                    axj = axy[0][jpoint]
                    ayj = axy[1][jpoint]

                    r = math.sqrt(  (xi-xj)**2 + (yi-yj)**2  )

                    axi += -self.G * (xi-xj) / r**3
                    ayi += -self.G * (yi-yj) / r**3
                flow = np.array([1, vxi, vyi, axi, ayi])
                vec = np.array([t, xi, yi, vxi, vyi])
                helper = RKHelper(flow)
                RK4_step( vec, self.dt, helper )
                xy[0][ipoint] = vec[1]
                xy[1][ipoint] = vec[2]
                vxy[0][ipoint] = vec[3]
                vxy[1][ipoint] = vec[4]
                axy[0][ipoint] = axi
                axy[1][ipoint] = axj                            
            for ival in xrange(len(xy)):
                for jval in xrange(len(xy[ival])) : 
                    if xy[ival][jval] > 10. :
                        xy[ival][jval] = -10.
                    if xy[ival][jval] < -10. :
                        xy[ival][jval] = 10.
            
            yield data

    def update(self, i):
        """Update the scatter plot."""
        data = next(self.stream)

        # Set x and y data...
        self.scat.set_offsets(data[:2, :])

        # We need to return the updated artist for FuncAnimation to draw..
        # Note that it expects a sequence of artists, thus the trailing comma.
        return self.scat,

    def show(self):
        plt.show()

if __name__ == '__main__':
    a = AnimatedScatter(dt=0.01, numpoints = 10, G=0.1)
    a.show()
