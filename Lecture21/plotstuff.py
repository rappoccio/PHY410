from read_multi_plot import *
import matplotlib
matplotlib.rcParams['legend.fancybox'] = True
import matplotlib.pyplot as plt
from matplotlib import legend

x,y = read_plot("pend_plot.txt")
y0=[]
y1=[]

for iy in y :
    y0.append( iy[0] )
    y1.append( iy[1] )

ax1 = plt.subplot(2,1,1)

p1, = ax1.plot(x,y0)


ax2 = plt.subplot(2,1,2)

p2, = ax2.plot(y0,y1)


plt.show()
