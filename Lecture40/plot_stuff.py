#!/usr/bin/env python

from read_plot import *
import matplotlib.pyplot as plt
import sys

if len(sys.argv) >= 2 :
    [x,y] = read_plot(sys.argv[1])
    plt.plot(x,y)
    plt.show()
