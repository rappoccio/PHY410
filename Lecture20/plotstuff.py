from read_plot import *
import matplotlib
matplotlib.rcParams['legend.fancybox'] = True
import matplotlib.pyplot as plt
from matplotlib import legend

x,y = read_plot("plot.txt")
xt,yt = read_plot("noAirPlot.txt")

ax1 = plt.subplot(1,1,1)

p1, = ax1.plot(x,y)
p2, = ax1.plot(xt,yt)

ax1.legend([p1, p2], ["Numerical", "Exact (no air)"])

plt.show()
