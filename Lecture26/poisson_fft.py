import math
import cmath
import time
import cpt

print " FFT solution of Poisson's equation"
print " ----------------------------------"
N = input(" Enter number points in x or y: ")

h = 1.0 / (N - 1)               # assume physical size in x and y = 1

start_time = time.clock()

q = 10.0                        # point charge
jq = kq = int(N / 2)            # at center of lattice
rho = cpt.Matrix(N, N)
for j in range(N):
    for k in range(N):
        if j == jq and k == kq:
            rho[j][k] = q / h**2
        else:
            rho[j][k] = 0.0

# FFT rows of rho
f = [ 0.0 ] * N                 # to store rows and columns
for j in range(N):
    for k in range(N):
        f[k] = rho[j][k]
    f = cpt.fft(f)
    for k in range(N):
        rho[j][k] = f[k]

# FFT columns of rho
for k in range(N):
    for j in range(N):
        f[j] = rho[j][k]
    f = cpt.fft(f)
    for j in range(N):
        rho[j][k] = f[j]

# Solve equation in Fourier space
V = cpt.Matrix(N, N)
W = cmath.exp(1.0j * 2 * math.pi / N)
Wm = Wn = 1.0 + 0.0j
for m in range(N):
    for n in range(N):
        denom = 4.0 - Wm - 1 / Wm - Wn - 1 / Wn
        if abs(denom) != 0.0:
            V[m][n] = rho[m][n] * h**2 / denom
        Wn *= W
    Wm *= W

# Inverse FFT rows of V
need_inverse = True             # to store rows and columns
for j in range(N):
    for k in range(N):
        f[k] = V[j][k]
    f = cpt.fft(f, need_inverse)
    for k in range(N):
        V[j][k] = f[k]

# Inverse FFT columns of V
for k in range(N):
    for j in range(N):
        f[j] = V[j][k]
    f = cpt.fft(f, need_inverse)
    for j in range(N):
        V[j][k] = f[j]

print " CPU time =", time.clock() - start_time, "sec"



# Convert x, y, V(x,y) to a surface plot
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib.mlab import griddata
import matplotlib.pyplot as plt
import numpy as np
# Define the axes
x = np.arange(0, h*(N), h)
y = np.arange(0, h*(N), h)
# Get the grid
X, Y = np.meshgrid(x, y)
# Set Z to the poisson V[i][j]
Z = np.array( V ).real

print X
print Y
print Z

fig = plt.figure()
ax = fig.gca(projection='3d')
scat = ax.plot_surface( X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                        linewidth=0, antialiased=False )

plt.show()

