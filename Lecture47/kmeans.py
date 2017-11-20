#!/usr/bin/env python
# Inspired by https://github.com/mubaris/friendly-fortnight/blob/master/kmeans-from-scratch-and-sklearn.py
import numpy as np
import math
import copy

import matplotlib.pyplot as plt

k = 8
n = 1000
eps = 1e-3



points = np.random.rand( n, 2 )
centroids = np.random.rand( k, 2 )
clusters = np.zeros( len(points) )

dcentroids = np.ones( k ) * 9999.
dcentroid =  math.sqrt( sum ( val**2  for i,val in enumerate(dcentroids) ) )

print 'points : '
print points

while dcentroid > eps :


    oldcentroids = copy.deepcopy( centroids)
    for i in xrange(n):
        distances = np.linalg.norm( points[i] - centroids, axis=1 )
        cluster = np.argmin( distances )
        clusters[i] = cluster
    for i in xrange(k):
        newvals = [ points[j] for j in range(len(points)) if clusters[j] == i ]        
        centroids[i] = np.mean( newvals, axis=0 )
    dcentroid = np.linalg.norm( centroids - oldcentroids )
    print 'dcentroid = ', dcentroid, ' centroids = '
    print centroids    


colors = ['r', 'g', 'b', 'y', 'c', 'm', 'darkviolet', 'brown']
fig, ax = plt.subplots()
for i in range(k):
    p = np.array([points[j] for j in range(len(points)) if clusters[j] == i])
    ax.scatter(p[:, 0], p[:, 1], s=20, c=colors[i])
ax.scatter(centroids[:, 0], centroids[:, 1], marker='*', s=500, c=colors)

plt.show()
