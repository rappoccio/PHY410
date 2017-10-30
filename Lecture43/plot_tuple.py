#!/usr/bin/env python

#from read_tuple import *
import matplotlib.pyplot as plt
import sys
import numpy
from matplotlib import legend
import matplotlib.pyplot as plt

if len(sys.argv) >= 2 :
    t,V,n,m,h = numpy.loadtxt(sys.argv[1], unpack=True)

    ax1 = plt.subplot(2,1,1)
    p1, = plt.plot(t,V)
    ax1.legend( [p1], ["Voltage"] )

    ax2 = plt.subplot(2,1,2)
    p2, = plt.plot(t,n)
    p3, = plt.plot(t,m)
    p4, = plt.plot(t,h)
    ax2.legend( [p2,p3,p4], ["n", "m", "h"] )
    
    plt.show()
