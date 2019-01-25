import math
import time
from cpt import *


class PoissonMG : 

    def __init__(self, accuracy = 0.001, L = 64, n_smooth = 5):
        self.accuracy = accuracy    # desired relative accuracy in solution
        self.L = L                  # number of interior points in each dimension
        self.n_smooth = n_smooth    # number of pre and post smoothing iterations

        self.h = 1.0 / (L + 1)           # assuming physical size of region in x and y = 1
        self.steps = 0                   # global variable for number of iteration steps
        
        # check that L is a power of 2 as required by multigrid
        power_of_2 = 1
        while (power_of_2 < L):
            power_of_2 *= 2
        if power_of_2 != L:
            self.L = power_of_2
            print " Setting L = " + repr(L) + " (must be a power of 2)"

        # recreate (L+2)x(L+2) matrices with elements set to zero
        self.psi = Matrix(L+2, L+2)
        self.psi_new = Matrix(L+2, L+2)
        self.rho = Matrix(L+2, L+2)

        # assume that the domain has unit size in x and y
        self.h = 1.0 / (L + 1)       # assume physical size in x and y = 1
        self.q = 10.0                # point charge
        i = L / 2                    # approximately at center of lattice
        self.rho[i][i] = self.q / self.h**2    # charge density


    def Gauss_Seidel(self, h, u, f):
        # find L for this grid
        self.L = len(u) - 2
        L = self.L

        # use checkerboard updating
        for color in range(2):
            for i in range(1, L+1):
                for j in range(1, L+1):
                    if (i + j) % 2 == color:
                        u[i][j] = 0.25 * ( u[i - 1][j] + u[i + 1][j] +
                                           u[i][j - 1] + u[i][j + 1] +
                                           h**2 * f[i][j] )

    def two_grid(self, h, u, f):

        # find L for this grid
        self.L = len(u) - 2
        L = self.L

        # solve exactly if there is only one interior point
        if L == 1:
            u[1][1] = 0.25 * ( u[0][1] + u[2][1] +
                               u[1][0] + u[1][2] + h**2 * f[1][1] )
            return

        # do a few pre-smoothing Gauss-Seidel steps
        for i in range(n_smooth):
            self.Gauss_Seidel(h, u, f)

        # find the residual
        r = Matrix(L+2, L+2)
        for i in range(1, L+1):
            for j in range(1, L+1):
                r[i][j] = f[i][j] + ( u[i + 1][j] + u[i - 1][j] +
                                      u[i][j + 1] + u[i][j - 1]- 4 * u[i][j] ) / h**2

        # restrict residual to coarser grid
        L2 = L / 2
        R = Matrix(L2+2, L2+2)
        for I in range(1, L2+1):
            i = 2 * I - 1
            for J in range(1, L2+1):
                j = 2 * J - 1
                R[I][J] = 0.25 * ( r[i][j] + r[i + 1][j] + r[i][j + 1] +
                                   r[i + 1][j + 1] )

        # initialize correction V on coarse grid to zero
        V = Matrix(L2+2, L2+2, 0.0)

        # call two_grid recursively
        H = 2 * h
        self.two_grid(H, V, R)

        # prolongate V to fine grid using simple injection
        v = Matrix(L+2, L+2)
        for I in range(1, L2+1):
            i = 2 * I - 1
            for J in range(1, L2+1):
                j = 2 * J - 1
                v[i][j] = v[i+1][j] = v[i][j+1] = v[i+1][j+1] = V[I][J]

        # correct u
        for i in range(1, L+1):
            for j in range(1, L+1):
                u[i][j] += v[i][j]

        # do a few post-smoothing Gauss-Seidel steps
        for i in range(n_smooth):
            self.Gauss_Seidel(h, u, f)

    def relative_error(self):
        error = 0.0             # average relative error per lattice point
        n = 0                   # number of non-zero differences
        for i in range(1, L+1):
            for j in range(1, L+1):
                if self.psi_new[i][j] != 0.0:
                    if self.psi_new[i][j] != self.psi[i][j]:
                        error += abs(1 - self.psi[i][j] / self.psi_new[i][j])
                        n += 1
        if n > 0:
            error /= n
        return error

print " Multigrid solution of Poisson's equation"
print " ----------------------------------------"
accuracy = input(" Enter desired accuracy in the solution: ")
L = input(" Enter number of interior points in x or y: ")
n_smooth = input(" Enter number of smooting iterations: ")

poissonmg = PoissonMG(accuracy,L,n_smooth)
start_time = time.clock()
while True:
    for i in range(L+2):
        for j in range(L+2):
            poissonmg.psi_new[i][j] = poissonmg.psi[i][j]
    poissonmg.two_grid(poissonmg.h, poissonmg.psi, poissonmg.rho)
    poissonmg.steps += 1
    error = poissonmg.relative_error()
    print " Step No. " + repr(poissonmg.steps) + "\tError = " + repr(error)
    if poissonmg.steps > 1 and error < accuracy:
        break

print " CPU time =", time.clock() - start_time, "sec"


# Convert x, y, V(x,y) to a surface plot
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib.mlab import griddata
import matplotlib.pyplot as plt
import numpy as np
# Define the axes
x = np.arange(0, poissonmg.h*(poissonmg.L+2), poissonmg.h)
y = np.arange(0, poissonmg.h*(poissonmg.L+2), poissonmg.h)
# Get the grid
X, Y = np.meshgrid(x, y)
# Set Z to the poissonmg V[i][j]
Z = np.array( poissonmg.psi )

fig = plt.figure()
ax = fig.gca(projection='3d')
scat = ax.plot_surface( X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                        linewidth=0, antialiased=False )

plt.show()

