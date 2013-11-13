from read_2d_plot import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib.mlab import griddata
import matplotlib.pyplot as plt
import numpy as np

s = "poisson_fft.data"
ix, iy, iz = read_plot(s)
L = input(" Enter number of grid points: " )
h = 1.0 / float(L)
N = L 

V = np.array( [[0.0] * N] * N)
for i in xrange(N) :
    for j in xrange(N) :
        print 'i,j = ' + str(i) + ', ' + str(j)
        V[i][j] = iz[i*N + j]

print V


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

