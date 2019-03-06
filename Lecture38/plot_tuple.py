#!/usr/bin/env python

import matplotlib.pyplot as plt
import sys
import numpy
from matplotlib import legend
import matplotlib.pyplot as plt

if len(sys.argv) >= 2 :
    x,y1,y2 = numpy.loadtxt(sys.argv[1], unpack=True)

    ax = plt.subplot(1,1,1)
    p1, = plt.plot(x,y1, '.')
    p2, = plt.plot(x,y2)
    ax.legend( [p1, p2], ["Numerical", "Exact"] )
    plt.show()
