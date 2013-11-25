import math
import time
import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys

class Animator :


    def __init__(self, ninit=500):
        self.fig, self.ax = plt.subplots()
        self.ax.set_ylim(-0.2,0.5)
        self.lastvals = [0.] * ninit
        self.line, = self.ax.plot(self.lastvals)
        

    def update(self, data) :
        self.line.set_ydata(data)
        return self.line,
        
    def time_step(self):
        try : 
            y = input()
            if y != '\n' :
                lastvals = y
                yield [float(iy) for iy in y]
        except EOFError :
            yield [float(iy) for iy in self.lastvals]

    def create_widgets(self):
        self.QUIT = Button(self, text="QUIT", command=self.quit)
        self.QUIT.pack(side=BOTTOM)

        self.draw = Canvas(self, width="600", height="400")
        self.draw.pack(side=TOP)

    def animate(self) :
        self.ani = animation.FuncAnimation( self.fig,        # Animate our figure
                                            self.update,     # Update function draws our data
                                            self.time_step,  # "frames" function does the time step, each iteration
                                            interval=50,     # 50 ms between iterations
                                            blit=False       # don't blit anything
                                            )


def main() : 
    animator = Animator(ninit = int(sys.argv[1]) )
    animator.animate()
    plt.show()

if __name__ == "__main__" :
    main()
