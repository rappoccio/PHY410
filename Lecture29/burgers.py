import math
import time
import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


class Burgers :

    def __init__(self, N=500, L=1.0, CFL_ratio=1.0) : 

        self.L = L                         # size of periodic region
        self.N = N                         # number of grid points
        self.dx = L / float(N)             # lattice spacing
        self.t = 0.0                       # time
        self.u_max = 1.0                   # maximum wave amplitude
        self.CFL_ratio = CFL_ratio         # Courant-Friedrichs-Lewy ratio
        self.dt = CFL_ratio * self.dx      # time step
        self.initial_waveform = "Sine"

        self.nu = 1.0e-6                   # kinematic viscosity
        self.u = [] ;
        self.u_new = []                    # the solution and its update
        self.next = [] ;
        self.prev = []                     # for periodic boundary conditions

        self.method = "Godunov"            # integration method
        self.step = 0                      # integration step number


        # create arrays for lattice vectors
        self.u = [ 0.0 for j in range(N) ]
        self.u_new = [ 0.0 for j in range(N) ]

        # initialize arrays for periodic boundary conditions
        self.next = [ j+1 for j in range(N) ]
        self.prev = [ j-1 for j in range(N) ]
        self.next[N-1] = 0
        self.prev[0] = N-1

        # reset lattice spacing and initialize waveform
        self.dx = L / float(N)
        self.u_max = 0.0
        for j in range(N):
            x = j * self.dx
            if self.initial_waveform == "Sine":
                self.u[j] = math.sin(2 * math.pi * x) + 0.5 * math.sin(math.pi * x)
            elif self.initial_waveform == "Step":
                if x > L/4.0 and x < 3*L/4.0:
                    self.u[j] = 1.0
                else:
                    self.u[j] = 0.0
            else:
                self.u[j] = 1.0
            self.u_max = max(abs(self.u[j]), self.u_max)

        # set time, step and step number
        self.t = 0.0
        self.dt = self.CFL_ratio * self.dx / self.u_max
        self.step = 0

        self.T = 1.0                 # time to travel length L
        self.frames_per_sec = 25     # animation rate for screen redraws
        

    def FTCS(self):
        for j in range(self.N):
            self.u_new[j]  = self.u[j] * (1.0 - self.dt / (2.0 * self.dx) * (self.u[self.next[j]] - self.u[self.prev[j]]))
            self.u_new[j] += self.nu * self.dt / self.dx**2 * (self.u[self.next[j]] + self.u[self.prev[j]] - 2 * self.u[j])

    def Lax(self):
        for j in range(self.N):
            self.u_new[j]  = (self.u[self.next[j]] + self.u[self.prev[j]]) / 2.0
            self.u_new[j] -= self.u[j] * self.dt / (2.0 * self.dx) * (self.u[self.next[j]] - self.u[self.prev[j]])
            self.u_new[j] += self.nu * self.dt / self.dx**2 * (self.u[self.next[j]] + self.u[self.prev[j]] - 2 * self.u[j])

    def Lax_Wendroff(self):
        F = [ self.u[j]**2 / 2.0 for j in range(self.N) ]     # flux vector
        for j in range(self.N):
            self.u_new[j]  = (self.u[j] + self.u[self.next[j]]) / 2.0
            self.u_new[j] -= self.dt / (2.0 * self.dx) * (F[self.next[j]] - F[j])
            self.u_new[j] += self.nu * self.dt / (2.0 * self.dx**2) * (
                        (self.u[self.next[j]] + self.u[self.prev[j]] - 2 * self.u[j])/ 2.0 +
                        (self.u[self.next[self.next[j]]] + self.u[j] - 2 * self.u[self.next[j]])/ 2.0 )
        for j in range(self.N):
            F[j] = self.u_new[j]**2 / 2.0
        for j in range(self.N):
            self.u_new[j]  = self.u[j] - self.dt / self.dx * (F[j] - F[self.prev[j]])
            self.u_new[j] += self.nu * self.dt / self.dx**2 * (self.u[self.next[j]] + self.u[self.prev[j]] - 2 * self.u[j])

    def Godunov(self):
        self.u_plus = [ self.u[j] for j in range(self.N) ]
        self.u_minus = [ self.u[j] for j in range(self.N) ]
        for j in range(self.N):
            if self.u_plus[j] < 0.0:
                self.u_plus[j] = 0.0
            if self.u_minus[j] > 0.0:
                self.u_minus[j] = 0.0
        F = [ 0.0 ] *self.N
        for j in range(self.N):
            f1 = self.u_plus[self.prev[j]]**2 / 2.0
            f2 = self.u_minus[j]**2 / 2.0
            if f1 > f2:
                F[self.prev[j]] = f1
            else:
                F[self.prev[j]] = f2
            f1 = self.u_plus[j]**2 / 2.0
            f2 = self.u_minus[self.next[j]]**2 / 2.0
            if f1 > f2:
                F[j] = f1
            else:
                F[j] = f2
            self.u_new[j]  = self.u[j]
            self.u_new[j] += self.nu * self.dt / self.dx**2 * (self.u[self.next[j]] + self.u[self.prev[j]] - 2 * self.u[j])
            self.u_new[j] -= self.dt / self.dx * (F[j] - F[self.prev[j]])


    def time_step(self):
        eval( 'self.' + self.method + "()")
        swap = self.u
        self.u = self.u_new
        self.u_new = swap
        self.t += self.dt;

class Animator :


    def __init__(self, periodic=True,burgers=None):
        self.avg_times = []
        self.burgers = burgers       
        self.t = 0.        
        self.fig, self.ax = plt.subplots()

        self.ax.set_ylim(-2.,2.)
        initvals = [ ix for ix in self.burgers.u]
        self.line, = self.ax.plot(initvals)
        

    def update(self, data) :
        self.line.set_ydata(data)
        return self.line,
        
    def time_step(self):
        while True :
            start_time = time.clock()
            self.burgers.time_step()
            end_time = time.clock()
            #print 'Tridiagnonal step in ' + str(end_time - start_time) 
            yield [ix for ix in self.burgers.u]

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


burgers = Burgers()
animator = Animator(burgers=burgers)
animator.animate()
plt.show()
