#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

mu, sigma = 1.1, 0.1
x2 = mu + sigma * np.random.randn(10000)
n, bins, patches = plt.hist(x2,50)
print n
print bins
print patches
plt.ylabel('Number')
plt.xlabel('Reconstructed / Truth')
plt.show()
