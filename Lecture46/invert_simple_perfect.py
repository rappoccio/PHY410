#!/usr/bin/env python
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import copy

def smear(x, mu, s ):
    return x
    #return (x*mu) + s * np.random.rand()

mu, sigma = 1.1, 0.1
N_pe = 100000
xtrue = np.random.rand(N_pe)
xreco = copy.copy(xtrue)
xvals = np.arange(0,1.0,0.1)
xbins = np.linspace(0,1.0, num=11)
w = np.ones_like( xtrue )
w /= N_pe * 0.1

for i,x in enumerate(xreco):
    xreco[i] = smear(x,mu,sigma)

plt.figure(1)
G, xbins2d, ybins2d, patches2d = plt.hist2d(xreco, xtrue, bins=[xbins,xbins], weights=w, cmap='Blues')
plt.xlabel('True Length (cm)')
plt.ylabel('Reconstructed Length (cm)')
plt.colorbar()

#data = [0.35, 0.32, 0.33, 0.33, 0.30]
data = [0.2, 0.3, 0.4, 0.5 ]

plt.figure(2)
d,xbinsout,patches = plt.hist(data, bins=xbins-0.05, histtype='stepfilled')
plt.xlabel('Length (cm)')
plt.ylabel('Number')

plt.figure(3)
Ginv = np.linalg.inv( np.matrix(G) )
m = np.matmul(Ginv, d).flatten().tolist()[0]
meas, = plt.plot(xvals, d)
inv, = plt.plot(xvals, m, '--')
plt.legend( [meas,inv], ['Measured', 'Reconstructed'])
plt.xlabel('Length (cm)')
plt.ylabel('Number')

plt.show()
