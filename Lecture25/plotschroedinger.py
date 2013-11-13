import matplotlib.pyplot as plt

from read_plot import *
from read_plot_tokened import *

E_data, F_data = read_plot("F.data")

x_data, phi_data = read_plot_tokened("phi.data")

s1 = plt.subplot(2,1,1)
for i in range(len(x_data)) :
	plt.plot( x_data[i], phi_data[i] )

s2 = plt.subplot(2,1,2)
plt.plot( E_data, F_data )
plt.ylim( [-10, 10] )

plt.show()
