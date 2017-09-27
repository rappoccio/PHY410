import math
import time
import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


class Advection : 


    def __init__(self, method="FTCS()", N=500, L=5.0, c=1.0, step_wave=False ):


        self.L = float(L)           # system size
        self.N = 500                # number of cells in x
        self.dx = float(L) / float(N)  # grid spacing
        self.c = 1.0                # wave speed
        self.t = 0.0                # time
        self.dt = self.dx / self.c  # time step
        self.step_number = 0        # integration step number
        self.method = method        # integration algorithm function
        self.x = []                 # grid points
        self.u = []                 # wave amplitude
        self.u_new = []             # amplitude at next step

        self.step_wave = step_wave  # True: step waveform  False: Gaussian cosine


        self.dx = L / float(N)
        self.x = [ i * self.dx for i in range(N+1) ]
        self.u = [ self.f_0(self.x[i]) for i in range(N+1) ]
        self.u_new = [ self.u[i] for i in range(N+1) ]

    def f_0(self, x):           # initial waveform
        self.x_0 = self.L / 2.0           # position
        self.sigma = 0.1 * self.L         # width
        if self.step_wave:
            if abs(x - self.x_0) < self.sigma:
                return 1.0
            else:
                return 0.0
        else:                   # Gaussian modulated cosine waveform
            k = math.pi / self.sigma
            gaussian = math.exp(-(x - self.x_0)**2 / (2 * self.sigma**2))
            return math.cos(k * (x - self.x_0)) * gaussian


    def FTCS(self):
        for i in range(self.N+1):
            i_plus_1 = i + 1
            if i == self.N:
                i_plus_1 = 0
            i_minus_1 = i - 1
            if i == 0:
                i_minus_1 = self.N
            self.u_new[i] = self.u[i] - self.c * self.dt / (2 * self.dx) * (self.u[i_plus_1] - self.u[i_minus_1])

    def Lax(self):
        for i in range(self.N+1):
            i_plus_1 = i + 1
            if i == self.N:
                i_plus_1 = 0
            i_minus_1 = i - 1
            if i == 0:
                i_minus_1 = self.N
            self.u_new[i] = (self.u[i_plus_1] + self.u[i_minus_1]) / 2.0
            self.u_new[i] -= self.c * self.dt / (2 * self.dx) * (self.u[i_plus_1] - self.u[i_minus_1])

    def Lax_Wendroff(self):
        D = (self.c * self.dt / self.dx)**2 / 2.0
        for i in range(self.N+1):
            i_plus_1 = i + 1
            if i == self.N:
                i_plus_1 = 0
            i_minus_1 = i - 1
            if i == 0:
                i_minus_1 = self.N
            self.u_new[i] = self.u[i] - self.c * self.dt / (2 * self.dx) * (self.u[i_plus_1] - self.u[i_minus_1])
            self.u_new[i] += D * (self.u[i_plus_1] + self.u[i_minus_1] - 2 * self.u[i])

    def take_step(self):
        eval('self.' + self.method )
        swap = self.u
        self.u = self.u_new
        self.u_new = swap
        self.t += self.dt
        self.step_number += 1

    def save_u(self, plot_number):
        file_name = "u_" + repr(plot_number) + ".data"
        file = open(file_name, "w")
        for i in range(self.N):
            file.write(repr(self.x[i]) + "\t" + repr(self.u[i]) + "\n")
        file.close()
        print " saved u(x,t) at t =", self.t, "in", file_name

class Animator :


    def __init__(self, periodic=True,advection=None):
        self.avg_times = []
        self.advection = advection
        self.t = 0.        
        self.fig, self.ax = plt.subplots()
        self.ax.set_ylim(-2.,2.)
        initvals = [ ix for ix in self.advection.u]
        self.line, = self.ax.plot(initvals)
        

    def update(self, data) :
        self.line.set_ydata(data)
        return self.line,
        
    def time_step(self):
        advection.take_step()
        yield [ix for ix in self.advection.u]

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



m = input(" Choose a numerical method: 1) FTCS, 2) Lax, 3) Lax-Wendroff : ")
if m == 1:
    method = "FTCS()"
elif m == 2:
    method = "Lax()"
else:
    method = "Lax_Wendroff()"
N = input(" Enter number of grid cells N = ")
L = input(" Enter L = " )
c = input(" Enter c = " )

#print " Time for wave to move one grid spacing = " + repr(L / N / c)
#dt = input(" Enter time step dt: ")

advection = Advection( method=method, N=N, L=L, c=c )
animator = Animator(advection=advection)
animator.animate()
plt.show()
