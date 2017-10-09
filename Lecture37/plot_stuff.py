#!/usr/bin/env python

from read_plot import *
import matplotlib.pyplot as plt

[x,y] = read_plot("psiSqd.data")
plt.plot(x,y)
plt.show()
