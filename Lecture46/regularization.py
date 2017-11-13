#!/usr/bin/env python
import numpy.linalg as la
import scipy.sparse.linalg as sla

def rlsq_solution_V1(X, y, l):
    n, m = X.shape
    I = np.identity(m)
    w = np.dot(np.dot(la.inv(np.dot(X.T, X) + l*I), X.T), y)
    return w
def rlsq_solution_V2(X, y, l):
    w = sla.lsqr(X, y, damp=l)[0]
    return w
def rlsq_solution_V3(X, y, l):
    w = sla.lsmr(X, y, damp=l)[0]
    return w



import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import copy

def smear(x, mu, s ):
    return s * np.random.randn() + mu + x

mu, sigma = 0.0, 0.1
N_pe = 100000
xtrue = np.random.rand(N_pe)
xreco = copy.copy(xtrue)
xvals = np.arange(0,1.0,0.1)
xbins = np.linspace(0,1.0, num=11)
w = np.ones_like( xtrue )
w /= N_pe * 0.1

for i,x in enumerate(xreco):
    xreco[i] = smear(x,mu=mu,s=sigma)

    
plt.figure(1)
G, xbins2d, ybins2d, patches2d = plt.hist2d(xtrue, xreco, bins=[xbins,xbins], weights=w, cmap='Blues')
plt.xlabel('True Length (cm)')
plt.ylabel('Reconstructed Length (cm)')
plt.colorbar()

#data = [0.35, 0.32, 0.33, 0.33, 0.30]
data = np.random.normal( loc=0.5, scale=0.2, size=1000 )

plt.figure(2)
d,xbinsout,patches = plt.hist(data, bins=xbins-0.05, histtype='stepfilled')
plt.xlabel('Length (cm)')
plt.ylabel('Number')

plt.figure(3)
l = 0.1
m = rlsq_solution_V1(G, d, l)
meas, = plt.plot(xvals, d)
inv, = plt.plot(xvals, m, '--')
plt.legend( [meas,inv], ['Measured', 'True'])
plt.xlabel('Length (cm)')
plt.ylabel('Number')

plt.show()
