import math
import time
import cpt

class Poisson :
    def __init__ (self, L=50):
        self.L = L                   # number of interior points in x and y
        self.omega = 1.88177         # over-relaxation parameter for L = 50
        self.N = L + 2               # interior plus two boundary points
        N=self.N
        self.V = cpt.Matrix(N, N)    # potential to be found
        self.rho = cpt.Matrix(N, N)  # given charge density
        self.VNew = cpt.Matrix(N, N) # new potential after each step
        self.h = 1.0 / (L + 1)       # lattice spacing assuming size in x and y = 1
        self.q = 10.0                # point charge
        i = N / 2                    # center of lattice
        self.rho[i][i] = self.q / self.h**2    # charge density


    def Jacobi(self) :               # Jacobi algorithm for a single iterative step
        VNew = self.VNew
        V    = self.V         #avoid lots of typing
        rho  = self.rho
        h    = self.h
        for i in range(1, self.L+1):
            for j in range(1, self.L+1):
                VNew[i][j] = 0.25 * (V[i-1][j] + V[i+1][j] +
                                     V[i][j-1] + V[i][j+1] +
                                     h**2 * rho[i][j])

    def GaussSeidel(self):          # Gauss-Seidel algorithm for one iterative step
        L = self.L
        VNew = self.VNew
        V    = self.V         #avoid lots of typing
        rho  = self.rho
        h    = self.h        

        # copy V to VNew
        for i in range(1, self.L+1):
            for j in range(1, self.L+1):
                VNew[i][j] =  V[i][j]

        # perform Gauss-Seidel update
        for i in range(1, self.L):
            for j in range(1, self.L+1):
                VNew[i][j] = 0.25 * (VNew[i-1][j] + VNew[i+1][j] +
                                     VNew[i][j-1] + VNew[i][j+1] +
                                     h**2 * rho[i][j])


    def SuccessiveOverRelaxation(self):
        L = self.L
        VNew = self.VNew
        V    = self.V         #avoid lots of typing
        rho  = self.rho
        h    = self.h
        omega= self.omega

        # update even sites in red-black scheme
        for i in range(1, self.L+1):
            for j in range(1, self.L+1):
                if (i + j) % 2 == 0:
                    VNew[i][j] = (1 - omega) * V[i][j] + omega / 4 * (
                                 V[i-1][j] + V[i+1][j] + V[i][j-1] +
                                 V[i][j+1] + h**2 * rho[i][j] )

        # update odd sites in red-black scheme
        for i in range(1, self.L+1):
            for j in range(1, self.L+1):
                if (i + j) % 2 != 0:
                    VNew[i][j] = (1 - omega) * V[i][j] + omega / 4 * (
                                 VNew[i-1][j] + VNew[i+1][j] + VNew[i][j-1] +
                                 VNew[i][j+1] + h**2 * rho[i][j] )

    def relativeError(self):
        L = self.L
        VNew = self.VNew
        V    = self.V         #avoid lots of typing
        rho  = self.rho
        h    = self.h
        omega= self.omega
        
        error = 0
        n = 0
        for i in range(1, self.L+1):
            for j in range(1, self.L+1):
                if VNew[i][j] != 0 and VNew[i][j] != V[i][j]:
                    error += abs(1 - V[i][j] / VNew[i][j])
                    n += 1
        if n != 0:
            error /= n
        return error





print " Iterative solution of Poisson's equation"
print " ----------------------------------------"
L = int(input(" Enter number of interior points in x or y: "))
poisson = Poisson(L=L)
accuracy = float(input(" Enter desired accuracy in solution: "))
choice = int( input(" Enter choice of algorithm, Jacobi (0), Gauss-Seidel (1) or Successive overrelaxation (2): ") )

start_time = time.clock()

steps = 0

while True:
    if choice == 0 : 
        poisson.Jacobi()
    elif choice == 1 : 
        poisson.GaussSeidel()
    else : 
        poisson.SuccessiveOverRelaxation()
    if poisson.relativeError() < accuracy:
        break
    for i in range(1, poisson.L+1):
        for j in range(1, poisson.L+1):
            poisson.V[i][j] = poisson.VNew[i][j]
    steps += 1

print " Number of steps =", steps
print " CPU time =", time.clock() - start_time, "sec"


# Convert x, y, V(x,y) to a surface plot
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib.mlab import griddata
import matplotlib.pyplot as plt
import numpy as np
# Define the axes
x = np.arange(0, poisson.h*(poisson.L+2), poisson.h)
y = np.arange(0, poisson.h*(poisson.L+2), poisson.h)
# Get the grid
X, Y = np.meshgrid(x, y)
# Set Z to the poisson V[i][j]
Z = np.array( poisson.V )

fig = plt.figure()
ax = fig.gca(projection='3d')
scat = ax.plot_surface( X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                        linewidth=0, antialiased=False )

plt.show()

