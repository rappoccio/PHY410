#!/usr/bin/env python

#from read_tuple import *
import matplotlib.pyplot as plt
import sys
import numpy
from matplotlib import legend
import matplotlib.pyplot as plt

if len(sys.argv) >= 2 :

    labels = ['u', 'd', 'ubar', 'dbar', 's', 'sbar', 'g']
    x,uv,dv,usea,dsea,s,sbar,glu = numpy.loadtxt(sys.argv[1], unpack=True)

    ax = plt.subplot(1,1,1)
    plt.title('Proton PDFs')
    plt.xlabel('x')
    plt.ylabel('x f(x)')
    
    ii = 0
    for iy in [uv,dv,usea,dsea,s,sbar,glu]: 
        p1, = plt.plot(x,iy, label=labels[ii]) 
        ii += 1
    ax.set_yscale("log", nonposy='clip')
    ax.legend(loc=3) ## lower left
    plt.show()
