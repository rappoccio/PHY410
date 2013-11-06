import matplotlib
matplotlib.rcParams['legend.fancybox'] = True
import matplotlib.pyplot as plt
from read_multi_plot import *


t, v = read_plot( "kepler_fixed.data" )

x1_plot = []
y1_plot = []

for iv in v :
    x1_plot.append( iv[0] )
    y1_plot.append( iv[1] )



plt.scatter( x1_plot, y1_plot )

plt.show()

