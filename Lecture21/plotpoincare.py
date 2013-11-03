from read_plot import *
import matplotlib
matplotlib.rcParams['legend.fancybox'] = True
import matplotlib.pyplot as plt
from matplotlib import legend

x,y = read_plot("poincare_plot.txt")

plt.scatter(x,y)
plt.ylim( [-180,180] )
plt.xlim( [-180,180] )

plt.show()
