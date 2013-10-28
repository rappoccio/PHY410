import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from read_plot_3d import read_plot_3d

x,y,z = read_plot_3d( 'Na2Cl2.data' )

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(x, y, z)

plt.show()

