from math import *
from bvp import *

print " Relaxing solution of u'' = - pi^2 (u+1) / 4"
print " with boundary conditions u(0) = 0, u(1) = 1";
acc = input(" Enter desired accuracy: ")

print ""
print " Iteration          Change              Error"
print " -----------------------------------------------------"

relaxer = BVPRelax( x0=0., x1=1., u0=0., u1=1., N=100 )

relaxer.guess()
relaxer.relax()
iteration = 1
while relaxer.error() > acc:
    relaxer.relax();
    iteration += 1
    print " %6d      %18.16g  %18.16g" % (iteration, relaxer.change(), relaxer.error())

file = open("relaxing.data", "w")
for i in range(relaxer.N+1):
    x = relaxer.x0 + i * relaxer.dx
    file.write(repr(x) + "\t" + repr(relaxer.u[i]) + "\n")
file.close()
print " N =", relaxer.N, "step solution in file relaxing.data"
