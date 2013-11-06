import matplotlib
matplotlib.rcParams['legend.fancybox'] = True
import matplotlib.pyplot as plt
from read_multi_plot import *

t, v = read_plot( "earthmoonsun_cpp.txt" )

x1_plot = []
y1_plot = []
x2_plot = []
y2_plot = []
x3_plot = []
y3_plot = []

for iv in v :
    x1_plot.append( iv[0] )
    y1_plot.append( iv[1] )
    x2_plot.append( iv[4] )
    y2_plot.append( iv[5] )
    x3_plot.append( iv[8] )
    y3_plot.append( iv[9] )



plt.scatter( x1_plot, y1_plot, c='blue', s=1 )
plt.scatter( x2_plot, y2_plot, c='green', s=11*1 )
plt.scatter( x3_plot, y3_plot, c='red', s=109*1 )

plt.show()

