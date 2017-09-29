import matplotlib
matplotlib.use('TkAgg')

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from read_orbit import read_orbit

class AnimatedScatter(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""
    def __init__(self, numpoints=50,filename="figure8.out"):
        self.numpoints = numpoints
        self.stream = self.data_stream()
        self.nobj, self.snapshots = read_orbit(filename)

        # Setup the figure and axes...
        self.fig, self.ax = plt.subplots()
        # Then setup FuncAnimation.
        self.ani = animation.FuncAnimation(self.fig, self.update, interval=50, 
                                           init_func=self.setup_plot, blit=True)

    def setup_plot(self):
        """Initial drawing of the scatter plot."""
        x, y = next(self.stream)
        self.scat = self.ax.scatter(x, y, animated=True)
        self.ax.axis([-10, 10, -10, 10])

        # For FuncAnimation's sake, we need to return the artist we'll be using
        # Note that it expects a sequence of artists, thus the trailing comma.
        return self.scat,

    def data_stream(self):
        i = 0
        x = [0.] * self.nobj
        y = [0.] * self.nobj
        while True:
            if i >= len(self.snapshots):
                i = 0
            snapshot = self.snapshots[i]
            for iparticle,particle in enumerate(snapshot):
                x[iparticle] = particle.pos[0]
                y[iparticle] = particle.pos[1]
            i += 1
            yield [x,y]

    def update(self, i):
        """Update the scatter plot."""
        data = next(self.stream)

        # Set x and y data...
        self.scat.set_offsets(data)

        # We need to return the updated artist for FuncAnimation to draw..
        # Note that it expects a sequence of artists, thus the trailing comma.
        return self.scat,

    def show(self):
        plt.show()

import sys

if __name__ == '__main__':
    if len(sys.argv) < 2 : 
        a = AnimatedScatter()        
    else :
        a = AnimatedScatter(filename=sys.argv[1])
    a.show()
