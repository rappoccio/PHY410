from read_plot import *
import matplotlib.pyplot as plt

m,e = read_plot("ising.data")
x = [i for i in xrange(len(m)) ]

plt.plot(x,m)
plt.show()
