from cpt import *


m1 = 1.0
m2 = 2.0
m3 = 1.0

# Matrix with masses as the diagonal elements
M = Matrix(3,3)
M[0][0] = m1
M[1][1] = m2
M[2][2] = m3
print "M = "
Matrix_print(M)

[detM, Minv] = inverse(M)
print "Minv = "
Matrix_print(Minv)


# "Spring" constants affecting each mass
k12 = 1.0
k23 = 1.0

K = [
    [   k12,    -k12,       0.0     ],
    [   -k12,   k12 + k23,  -k23    ],
    [   0.0,    -k23,       k23     ]
]

print "K = "
Matrix_print(K)




# Solve with generalized eigenvector solution
[ eigenvalues, eigenvectors ] = solve_eigen_generalized(K, M)

print "Eigenvalues ="
print eigenvalues
print "Eigenvectors ="
Matrix_print(eigenvectors)
