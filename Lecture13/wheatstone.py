from cpt import *

print " Unbalanced Wheatstone bridge equations"
print " --------------------------------------"

v0 = 1.5
r1 = 100.
r2 = r1
r3 = 150.
rx = 120.
ra = 1000.
rv = 10.

v = Matrix(3, 1)       # column vector with 3 rows
v[0][0] = v0
print 'v = '
print v 

R = Matrix(3, 3)       # 3x3 resistance matrix
R[0][0] = r1 + rv      # set components using slicing notation
R[0][1] = r2
R[0][2] = rv

R[1][0] = r1 + ra
R[1][1] = -ra
R[1][2] = -r3

R[2][0] = rx + ra       # or use subscripting notation
R[2][1] = -r2 - rx - ra
R[2][2] = rx

print 'R = '
print R

# the solve_Gauss_Jordan replaces R by R^-1 and v by i
# so save the original R and copy v into a vector i
R_save = Matrix_copy(R)
i = Matrix_copy(v)

solve_Gauss_Jordan(R, i)
print " Solution using Gauss-Jordan elimination"
print " i = "
print i

# find the other currents
print " i_a = i_1 - i_2 = " + str(i[0][0] - i[1][0])
print " i_v = i_1 + i_3 = " + str(i[0][0] + i[2][0])
print " i_x = i_1 - i_2 + i_3 = " + str(i[0][0] - i[1][0] + i[2][0])


# see whether LU decomposition gives the same answer
i = Vector(v)
solve_LU_decompose(R_save, i)
print " Solution using LU Decompositon"
print " i = "
print i

