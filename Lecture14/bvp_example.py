from cpt import *
from math import *

# Solves a boundary value problem : 
#
#  d^2u / dx^2 = -pi^2/4 ( u+1)
#
# with Dirichlet boundary conditions u(0)=0 and u(1)=1.
# This is discretized as : 
# (2 u_i - u_(i+1) - u_(i-1)) / h^2 = pi^2/4 * (u_i + 1)
# 

n = 2                  # Number of masses
h = 0.01               # Discretization width
k = h*h * pi*pi * 0.25 # constant in the tridiagonal matrix

u = [0.0] * n
a = [0.0] * n
b = [0.0] * n
c = [0.0] * n
r = [0.0] * n

# To check : fill the whole matrix
A = Matrix (n,n)

# Dirichlet boundary conditions
u0 = 0.0
un = 1.0

# Set up our tridiagonal matrix
for i in xrange(n) :

  b[i] = 2.0 - k

  if  i != 0 :
    a[i] = -1.0
  else :
    a[i] = 0.0


  if  i != n-1 :
    c[i] = -1.0
  else :
    c[i] = 0.0

  if  i == 0 :
    r[i] = u0 + k
  elif i == n-1 :
    r[i] = un + k
  else  :
    r[i] = k


for i in xrange( n ) : 
  A[i][i] = 2-k
  if i == 0 :
    A[i][i+1] = -1.0
  elif i == n-1 :
    A[i][i-1] = -1.0
  else :
    A[i][i+1] = -1.0
    A[i][i-1] = -1.0



print "a = "
print a
print "b = "
print b
print "c = "
print c

print "r = "
print r

u = solve_tridiag( a, b, c, r, n)

print "u = "
print u

print "Another way : "

print "A = "
print A

[det, Ainv] = inverse(A)

print "Ainv = "
print Ainv

