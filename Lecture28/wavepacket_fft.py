import math
import cmath
import cpt
import time
import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


class WavepacketFFT : 

    def __init__( self, N=128, L=100., dt=0.1 ):


        self.h_bar = 1.0             # Planck's constant / 2pi in natural units
        self.mass = 1.0              # particle mass in natural units

        # The spatial grid
        self.N = N                   # number of interior grid points
        self.L = L                   # system extends from x = 0 to x = L
        self.dx = L / float(N)       # grid spacing
        self.dt = dt                 # time step
        self.x = []                  # vector of grid points

        # The potential V(x)
        self.V_0 = 0.5               # height of potential barrier
        self.V_width = 10.0          # width of potential barrier
        self.V_center = 0.75 * L     # center of potential barrier



        # Initial wave packet
        self.x_0 = L / 4.0           # location of center
        self.E = 1.0                 # average energy
        self.sigma_0 = L / 10.0      # initial width of wave packet
        self.psi_norm = 1.0          # norm of psi
        self.k_0 = 0                 # average wavenumber
        self.velocity = 0            # average velocity

        self.t = 0.0                 # time
        self.psi = []                # complex wavefunction
        self.T_exp_factor = []       # precomputed phase rotation from kinetic energy
        self.V_exp_factor = []       # precomputed phase rotation from potential energy


        # reset global vectors
        czero = 0.0 + 0.0j
        self.psi = [czero for j in range(N)]
        self.T_exp_factor = [czero for j in range(N)]
        self.V_exp_factor = [czero for j in range(N)]

        # reset the lattice
        self.dx = L / float(N)
        self.x = [ j * self.dx for j in range(N) ]

        # initialize the packet
        self.k_0 = math.sqrt(2*self.mass*self.E - self.h_bar**2 / 2.0 / self.sigma_0**2) / self.h_bar
        self.velocity = self.k_0 / self.mass
        self.psi_norm = 1 / math.sqrt(self.sigma_0 * math.sqrt(math.pi))
        for j in range(N):
            exp_factor = math.exp( - (self.x[j] - self.x_0)**2 / (2.0 * self.sigma_0**2))
            self.psi[j] = (math.cos(self.k_0 * self.x[j]) + 1j * math.sin(self.k_0 * self.x[j]))
            self.psi[j] *= exp_factor * self.psi_norm

        # initialize the phase rotation factors
        for j in range (N):

            # kinetic factor exp(-iT/h_bar dt)
            if j < int(N / 2):
                p = j
            else:
                p = j - N
            p *= self.h_bar * 2 * math.pi / L
            theta = - p**2 / (2 * self.mass) / self.h_bar * self.dt
            self.T_exp_factor[j] = math.cos(theta) + 1j * math.sin(theta)

            # potential factor exp[-iV(x)/(2h_bar) dt]
            theta = - self.V(self.x[j]) / 2.0 / self.h_bar * self.dt
            self.V_exp_factor[j] = math.cos(theta) + 1j * math.sin(theta)


    def V(self,x):
        half_width = abs(0.5 * self.V_width)
        if abs(x - self.V_center) <= half_width:
            return self.V_0
        else:
            return 0

    def take_step(self):

        # first half of potential phase rotation
        for j in range(self.N):
            self.psi[j] *= self.V_exp_factor[j]

        # FFT to momentum space
        self.psi = cpt.fft(self.psi)

        # kinetic phase rotation
        for j in range(self.N):
            self.psi[j] *= self.T_exp_factor[j]

        # FFT back to position space
        do_inverse = True
        self.psi = cpt.fft(self.psi, do_inverse)

        # second half of potential phase rotation
        for j in range(self.N):
            self.psi[j] *= self.V_exp_factor[j]

        self.t += self.dt





class Animator :


    def __init__(self, periodic=True,wavepacketFFT=None):
        self.periodic = periodic
        self.wavepacketFFT = wavepacketFFT       
        self.t = 0.        
        self.fig, self.ax = plt.subplots()

        self.myline = plt.axvline( x=(self.wavepacketFFT.V_center - 0.5 * self.wavepacketFFT.V_width)/ self.wavepacketFFT.dx,
                                   color='r'
            )
        self.myline = plt.axvline( x=(self.wavepacketFFT.V_center + 0.5 * self.wavepacketFFT.V_width)/ self.wavepacketFFT.dx,
                                   color='r'
            )
        self.ax.set_ylim(0,0.5)
        initvals = [ abs(ix) for ix in self.wavepacketFFT.psi]
        self.line, = self.ax.plot(initvals)
        

    def update(self, data) :
        self.line.set_ydata(data)
        return self.line,
        
    def time_step(self):
        while True :
            start_time = time.clock()
            self.wavepacketFFT.take_step()
            end_time = time.clock()
            print 'FFT step in ' + str(end_time - start_time)             
            yield [abs(ix) for ix in self.wavepacketFFT.psi]

    def create_widgets(self):
        self.QUIT = Button(self, text="QUIT", command=self.quit)
        self.QUIT.pack(side=BOTTOM)

        self.draw = Canvas(self, width="600", height="400")
        self.draw.pack(side=TOP)

    def animate(self) :
        self.ani = animation.FuncAnimation( self.fig, self.update,   self.time_step, interval=50, blit=False )


wavepacketFFT = WavepacketFFT(N=600)
animator = Animator(periodic=True,wavepacketFFT=wavepacketFFT)
animator.animate()
plt.show()
