import matplotlib.pyplot as plt

from read_plot import *


x1, y1 = read_plot("poor.data")
x2, y2 = read_plot("park_miller.data")
x3, y3 = read_plot("random.random.data")

s1 = plt.subplot(3,1,1)
plt.scatter(x1,y1 )

s2 = plt.subplot(3,1,2)
plt.scatter(x2,y2 )

s3 = plt.subplot(3,1,3)
plt.scatter(x3,y3 )

plt.show()
