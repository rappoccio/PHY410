from math import exp, pi, sqrt
from cpt import minimize_BFGS

# physical constants
hbar = 1.0                  # Planck's constant / 2pi
m = 1.0                     # electron mass
e = 1.0                     # proton charge

# LCAO variational wave function
# psi = sum( d_i g(alpha_i, r) ) for i = 0, 1, 2, ...
# assume d_0 = 1 and vary alpha_0, d_1, alpha_1, d_2, alpha_2, ...
# vector of variational parameters
p = [ 1.0, 1.0, 0.5 ]       # initial guess for [ alpha_0, d_1, alpha_1 ]
N = int( (len(p) + 1) / 2 ) # number of Gaussians

accuracy = 1.0e-6           # desired accuracy for numerical operations

def g(alpha, r):            # normalized s-wave Gaussian orbital
    return (2.0 * alpha / pi)**(3.0/4.0) * exp(-alpha * r**2)

def Sij(alpha_i, alpha_j):  # matrix elements of S
    return (pi / (alpha_i + alpha_j))**(3.0/2.0)

def Tij(alpha_i, alpha_j):  # matrix elements of T
    return (3.0 * hbar**2 / m * alpha_i * alpha_j *
            pi**(3.0/2.0) / (alpha_i + alpha_j)**(5.0/2.0))

def Vij(alpha_i, alpha_j):  # matrix elements of V
    return - 2.0 * e**2 * pi / (alpha_i + alpha_j)

def E(alpha, d):            # energy as function of N alpha_i and d_i
    S = H = 0.0
    for i in range(len(alpha)):
        for j in range(len(alpha)):
            fac = (alpha[i] * alpha[j])**(3.0/4.0)* d[i] * d[j]
            H += fac * (Tij(alpha[i], alpha[j]) + Vij(alpha[i], alpha[j]))
            S += fac * Sij(alpha[i], alpha[j])
    return H / S

def func(p):                # function for BFGS minimization
    # assume p = [ alpha_0, d_1, alpha_1, d_2, alpha_2, ... ]
    alpha = [ max(p[2 * i], accuracy) for i in range(N) ]
    d = [ 1.0 ]
    d.extend(p[2 * i + 1] for i in range(N - 1))
    return E(alpha, d)

def dfunc(p, g):            # gradient of func for BFGS minimization
    # use symmetric finite difference f'(x) = (f(x+eps) - f(x-eps)) / (2 eps)
    eps = 0.5 * accuracy    # finite difference
    for i in range(len(p)):
        p1 = list(p)
        p1[i] += eps
        p2 = list(p)
        p2[i] -= eps
        g[i] = (func(p1) - func(p2)) / (2 * eps)
    return

def norm(p):                # norm of LCAO
    alpha = [ p[2 * i] for i in range(N) ]
    d = [ 1.0 ]
    d.extend(p[2 * i + 1] for i in range(N - 1))
    norm = 0.0
    for i in range(N):
        for j in range(N):
            norm += Sij(alpha[i], alpha[j]) * d[i] * d[j]
    return sqrt(norm)

print(" Variational method for Hydrogen using Gaussian LCAO")
print(" Minimize <psi|H|psi>/<psi|psi> using BFGS algorithm")
gtol = accuracy
iterations, e = minimize_BFGS(p, gtol, func, dfunc)
print(" number of Gaussians N = " + repr(N))
print(" number of iterations = " + repr(iterations))
print(" energy E = " + repr(e))
print(" i\talpha_i\t\t\td_i")
exit
for i in range(N):
    alpha_i = p[2 * i]
    if i == 0:
        d_i = 1.0 / norm(p)
    else:
        d_i = p[2 * i - 1] / norm(p)
    print(" " + repr(i) + "\t" + repr(p[2*i]) + "\t" + repr(d_i))
