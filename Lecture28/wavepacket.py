import math
import cmath
import time
import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


class Wavepacket : 

    def __init__(self, N=600, L=100., dt=0.1, periodic=True):

        self.h_bar = 1.0             # Planck's constant / 2pi in natural units
        self.mass = 1.0              # particle mass in natural units

        # The spatial grid
        self.N = N                   # number of interior grid points
        self.L = L                   # system extends from x = 0 to x = L
        self.dx = L / float(N + 1)   # grid spacing
        self.dt = dt                 # time step
        self.x = []                  # vector of grid points
        self.periodic = periodic     # True = periodic, False = Dirichlet boundary conditions

        # The potential V(x)
        self.V_0 = -0.5               # height of potential barrier
        self.V_width = 10.0          # width of potential barrier
        self.V_center = 0.75 * L     # center of potential barrier
        self.gaussian = True         # True = Gaussian potential, False = step potential


        # Initial wave packet
        self.x_0 = L / 4.0           # location of center
        self.E = 1.0                 # average energy
        self.sigma_0 = L / 10.0      # initial width of wave packet
        self.psi_norm = 1.0          # norm of psi
        self.k_0 = 0.0               # average wavenumber
        self.velocity = 0.0          # average velocity

        self.t = 0.0                 # time
        self.psi = []                # complex wavefunction
        self.chi = []                # wavefunction for simplified Crank-Nicholson
        self.a = []
        self.b = []                  # to represent tridiagonal elements of matrix Q
        self.c = []
        self.alpha = 0.0
        self.beta = 0.0              # corner elements of matrix Q

        # reset global vectors
        self.psi = [0 + 1j*0 for j in range(N)]
        self.chi = [0 + 1j*0 for j in range(N)]

        # reset time and the lattice
        self.t = 0.0
        self.dx = L / float(self.N + 1)
        self.x = [ float(j * self.dx) for j in range(self.N) ]

        # initialize the packet
        self.k_0 = math.sqrt(2*self.mass*self.E - self.h_bar**2 / 2 / self.sigma_0**2) / self.h_bar
        self.velocity = self.k_0 / self.mass
        self.psi_norm = 1 / math.sqrt(self.sigma_0 * math.sqrt(math.pi))
        for j in range(self.N):
            exp_factor = math.exp( - (self.x[j] - self.x_0)**2 / (2 * self.sigma_0**2))
            self.psi[j] = (math.cos(self.k_0 * self.x[j]) + 1j * math.sin(self.k_0 * self.x[j]))
            self.psi[j] *= exp_factor * self.psi_norm

        # elements of tridiagonal matrix Q = (1/2)(1 + i dt H / (2 hbar))
        for j in range(self.N):
            self.a.append( - 1j * self.dt * self.h_bar / (8 * self.mass * self.dx**2) )
            self.b.append( 0.5 + 1j * self.dt / (4 * self.h_bar) *
                      (self.V(self.x[j]) + self.h_bar**2 / (self.mass * self.dx**2)) )
            self.c.append( - 1j * self.dt * self.h_bar / (8 * self.mass * self.dx**2) )
        self.alpha = self.c[N-1]
        self.beta = self.a[0]

        


    def V(self, x):
        half_width = abs(0.5 * self.V_width)
        if self.gaussian:
            return self.V_0 * math.exp(-(x - self.V_center)**2 / (2 * half_width**2))
        else:
            if abs(x - self.V_center) <= half_width:
                return self.V_0
            else:
                return 0.0


    def solve_tridiagonal(self, a, b, c, r, u):
        n = len(r)
        gamma = [ 0 + 1j*0 for j in range(n) ]
        beta = b[0]
        u[0] = r[0] / beta
        for j in range(1, n):
            gamma[j] = c[j-1] / beta
            beta = b[j] - a[j] * gamma[j]
            u[j] = (r[j] - a[j] * u[j-1]) / beta
        for j in range(n-2, -1, -1):
            u[j] -= gamma[j+1] * u[j+1]

    def solve_tridiagonal_cyclic(self, a, b, c, alpha, beta, r, x):
        n = len(r)
        bb = [0 + 1j*0 for j in range(self.N)]
        u = [0 + 1j*0 for j in range(self.N)]
        z = [0 + 1j*0 for j in range(self.N)]
        gamma = -b[0]
        bb[0] = b[0] - gamma
        bb[n-1] = b[n-1] - alpha * beta / gamma
        for i in range(1, n-1):
            bb[i] = b[i]
        self.solve_tridiagonal(a, bb, c, r, x)
        u[0] = gamma
        u[n-1] = alpha
        for i in range(1, n-1):
            u[i] = 0
        self.solve_tridiagonal(a, bb, c, u, z)
        fact = x[0] + beta * x[n-1] / gamma
        fact /= 1.0 + z[0] + beta * z[n-1] / gamma
        for i in range(n):
            x[i] -= fact * z[i]

#    T = 5.0                 # time to travel length L


class Animator :


    def __init__(self, periodic=True,wavepacket=None):
        self.avg_times = []
        self.periodic = periodic
        self.wavepacket = wavepacket       
        self.t = 0.        
        self.fig, self.ax = plt.subplots()

        self.myline = plt.axvline( x=(self.wavepacket.V_center - 0.5 * self.wavepacket.V_width)/ self.wavepacket.dx,
                                   color='r'
            )
        self.myline = plt.axvline( x=(self.wavepacket.V_center + 0.5 * self.wavepacket.V_width)/ self.wavepacket.dx,
                                   color='r'
            )
        self.ax.set_ylim(0,0.5)
        initvals = [ abs(ix) for ix in self.wavepacket.psi]
        self.line, = self.ax.plot(initvals)
        

    def update(self, data) :
        self.line.set_ydata(data)
        return self.line,
        
    def time_step(self):
        while True :
            start_time = time.clock()
            if self.periodic:
                self.wavepacket.solve_tridiagonal_cyclic(self.wavepacket.a, self.wavepacket.b,
                                                         self.wavepacket.c, self.wavepacket.alpha, self.wavepacket.beta,
                                                         self.wavepacket.psi, self.wavepacket.chi)
            else:
                self.wavepacket.solve_tridiagonal(self.wavepacket.a, self.wavepacket.b,
                                                  self.wavepacket.c, self.wavepacket.psi, self.wavepacket.chi)
            for j in range(self.wavepacket.N):
                self.wavepacket.psi[j] = self.wavepacket.chi[j] - self.wavepacket.psi[j]
            self.t += self.wavepacket.dt;
            end_time = time.clock()
            print 'Tridiagnonal step in ' + str(end_time - start_time) 
            yield [abs(ix) for ix in self.wavepacket.psi]

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


wavepacket = Wavepacket(N=128)
animator = Animator(periodic=True,wavepacket=wavepacket)
animator.animate()
plt.show()
