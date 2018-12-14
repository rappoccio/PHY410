#!/usr/bin/env python

# Computational Physics Toolkit cpt.py

import math
import cmath
import sys
import copy
import scipy.optimize

from cmath import exp, pi
from math import sin, cos
import matplotlib.pyplot as plt

from numpy import array, conj, divide


def dft(data, inverse=False):
    """Return Discrete Fourier Transform (DFT) of a complex data vector"""
    N = len(data)
    transform = [ 0 ] * N
    for k in range(N):
        for j in range(N):
            angle = 2 * math.pi * k * j / float(N)
            if inverse:
                angle = -angle
            transform[k] += data[j] * cmath.exp(1j * angle)
    if inverse:
        for k in range(N):
            transform[k] /= float(N)
    return transform

def fft(data, inverse=False):
    """Return Fast Fourier Transform (FFT) using Danielson-Lanczos Lemma"""
    N = len(data)
    if N == 1:               # transform is trivial
        return [data[0]]
    elif N % 2 == 1:         # N is odd, lemma does not apply
        return dft(data, inverse)
    # perform even-odd decomposition and transform recursively
    even = fft([data[2*j] for j in range(N//2)], inverse)
    odd  = fft([data[2*j+1] for j in range(N//2)], inverse)
    W = cmath.exp(1j * 2 * math.pi / N)
    if inverse:
        W = 1.0 / W
    Wk = 1.0
    transform = [ 0 ] * N
    for k in range(N):
        transform[k] = even[k % (N//2)] + Wk * odd[k % (N//2)]
        Wk *= W
    if inverse:
        for k in range(N):
            transform[k] /= 2.0
    return transform

def sine_fft(data, inverse=False):
    """Return Fast Sine Transform of N data values starting with zero."""
    N = len(data)
    if data[0] != 0.0:
        raise Exception("data[0] != 0.0")
    extend_data = [ 0.0 ] * (2*N)
    for j in range(1, N):
        extend_data[j] = data[j]
        extend_data[2*N-j] = -data[j]
    transform = fft(extend_data)
    sineft = [ 0.0 ] * N
    for k in range(N):
        sineft[k] = transform[k].imag / 2.0
        if inverse:
            sineft[k] *= 2.0 / N
    return sineft

def fft_power(x) :
    N = len(x)
    if N <=1 : return x
    power = [ 0.0 ] * (N/2+1)
    power[0] = abs(x[0])**2
    for i in range(1,N/2) :
        power[i] = abs(x[i])**2 + abs(x[N-i])**2
    power[N/2] = abs(x[N/2])**2
    for i in range(0,N/2+1) :
        power[i] /= float(N)
    return power


def cosine_fft(data, inverse=False):
    """Return Fast Cosine Transform of (N+1) data values
       including two boundary points with index 0, N.
    """
    N = len(data)-1
    extend_data = [ 0.0 ] * (2*N)
    extend_data[0] = data[0]
    for j in range(1, N):
        extend_data[j] = data[j]
        extend_data[2*N-j] = data[j]
    extend_data[N] = data[N]
    transform = fft(extend_data)
    cosineft = [ 0.0 ] * (N+1)
    for k in range(N+1):
        cosineft[k] = transform[k].real / 2.0
        if inverse:
            cosineft[k] *= 2.0 / N
    return cosineft





def ifft(x) :
    # conjugate the complex numbers
    x = conj(x)
 
    # forward fft
    X = fft( x );
 
    # conjugate the complex numbers again
    X = conj(X)
 
    # scale the numbers
    X = divide(X, len(X))

    return X
    



# Define an alias for compatibility with C++

Vector = list

def Matrix(rows, cols, initial_value=0.0):
    """ Constructs new matrix with specified number of rows and col(umn)s.
        The matrix elements are initialized to initial_value
    """
    m = []
    for row in range(rows):
        m.append([initial_value for col in range(cols)])
    return m

def Matrix_copy(m):
    """ Return a copy.
    """
    c = Matrix(len(m), len(m[0]))
    for i in range(len(m)):
        for j in range(len(m[0])):
            c[i][j] = m[i][j]
    return c

def Matrix_multiply(A, B):
    """ Returns the matrix product of matrices A and B.
        Number of columns of A must equal number of rows of B.
    """
    if len(A[0]) != len(B):
        raise Exception("Number of columns of A not = number of rows of B")
    C = Matrix(len(A), len(B[0]))
    for i in range(len(A)):
        for j in range(len(B[0])):
            for k in range(len(B)):
                C[i][j] += A[i][k] * B[k][j]
    return C

def Matrix_print(m):
    """ Prints rows and columns in table form.
    """
    for row in range(len(m)):
        values = []
        line = ''
        for col in range(len(m[row])):
            line += '{' + str( len(values) ) + ':6.2f} '
            values.append( m[row][col] )
        s = line.format( *values )        
        print s

def Matrix_transpose(m):
    """ Returns a new matrix with [i][j] element = m[j][i].
    """
    transpose = Matrix(len(m[0]), len(m))
    for i in range(len(m[0])):
        for j in range(len(m)):
            transpose[i][j] = m[j][i]
    return transpose

def solve_Gauss_Jordan(A, B):
    """ Solve A x = B using Gauss-Jordan elimination.
        A must be a square matrix.  B has one or more columns.
        The solution x replaces B and A is destroyed.
    """
    n = len(A)
    if n != len(A[0]):
        raise Exception("A must be square")
    if n != len(B):
        raise Exception("Number of rows of B not equal to colums of A")
    m = len(B[0])
    indxc = [ 0 for i in range(n) ]
    indxr = [ 0 for i in range(n) ]
    ipiv = [ 0 for i in range(n) ]
    for i in range(n):
        big = 0.0
        for j in range(n):
            if ipiv[j] != 1:
                for k in range(n):
                    if ipiv[k] == 0:
                        Ajk = A[j][k];
                        if abs(Ajk) >= big:
                            big = abs(Ajk)
                            irow = j
                            icol = k
        ipiv[icol] += 1
        if irow != icol:
            for j in range(n):
                swap= A[irow][j]
                A[irow][j] = A[icol][j]
                A[icol][j] = swap
            for j in range(m):
                swap = B[irow][j]
                B[irow][j] = B[icol][j]
                B[icol][j] = swap
        indxr[i] = irow
        indxc[i] = icol
        if A[icol][icol] == 0.0:
            raise Exception("solve_Gauss_Jordan singular matrix")
        pivinv = 1.0 / A[icol][icol]
        A[icol][icol] = 1.0
        for j in range(n):
            A[icol][j] *= pivinv
        for j in range(m):
            B[icol][j] *= pivinv
        for j in range(n):
            if j != icol:
                mult = A[j][icol]
                for k in range(n):
                    A[j][k] -= A[icol][k] * mult
                for k in range(m):
                    B[j][k] -= B[icol][k] * mult
    for i in range(n - 1, -1, -1):
        if indxr[i] != indxc[i]:
            for j in range(n):
                swap = A[j][indxr[i]]
                A[j][indxr[i]] = A[j][indxc[i]]
                A[j][indxc[i]] = swap

def ludcmp(A, indx):
    """ LU Decomposes Matrix A.
    """
    TINY = 1.0e-20
    n = len(A)
    vv = [ 0.0 for i in range(n) ]
    d = 1.0
    for i in range(n):
        big = 0.0
        for j in range(n):
            big = max(abs(A[i][j]), big)
        if big == 0.0:
            raise Exception("Singular matrix in routine ludcmp")
        vv[i] = 1.0 / big
    for j in range(n):
        for i in range(j):
            sm = A[i][j]
            for k in range(i):
                sm -= A[i][k] * A[k][j]
            A[i][j] = sm
        big = 0.0
        for i in range(j, n):
            sm = A[i][j]
            for k in range(j):
                sm -= A[i][k] * A[k][j]
            A[i][j] = sm
            dum = vv[i] * abs(sm)
            if dum >= big:
                big = dum
                imax = i
        if j != imax:
            for k in range(n):
                dum = A[imax][k]
                A[imax][k] = A[j][k]
                A[j][k] = dum
            d = -d
            vv[imax] = vv[j]
        indx[j] = imax
        if A[j][j] == 0.0:
            A[j][j]=TINY;
        if j != n-1:
            dum = 1.0 / A[j][j]
            for i in range(j+1, n):
                A[i][j] *= dum
    return d

def lubksb(A, indx, B):
    """ LU Decomposition.
    """
    n = len(A)
    ii = 0
    for k in range(len(B[0])):
        for i in range(n):
            ip = indx[i]
            sm = B[ip][k]
            B[ip][k] = B[i][k]
            if ii != 0:
                for j in range(ii-1, i):
                    sm -= A[i][j] *B [j][k]
            elif sm != 0.0:
                ii = i+1
            B[i][k] = sm
        for i in range(n-1, -1, -1):
            sm = B[i][k]
            for j in range(i+1, n):
                sm -= A[i][j] * B[j][k]
            B[i][k] = sm / A[i][i]

def solve_LU_decompose(A, B):
    """ Solve A x = B using LU Decomposition.
        A must be a square matrix.  B has one or more columns.
        The solution x replaces B and A is destroyed.
    """
    n = len(A)
    if n != len(A[0]):
        raise Exception("A must be square")
    if n != len(B):
        raise Exception("Number of rows of B not equal to colums of A")
    m = len(B[0])
    indx = [ 0 for i in range(n) ]
    d = ludcmp(A, indx)
    lubksb(A, indx, B)


# Adapted from Numerical Recipes in C, 1992, 
# "tridag" method, Section 2.4. 
def solve_tridiag( a, b, c, r, n) :
    bet = 0.0
    gam = [0.0] * n
    u = [0.0] * n
    
    if b[0] == 0.0  :
        return u

    bet = b[0]
    u[0] = r[0] / bet
    for j in range(1,n) :
        gam[j] = c[j-1] / bet
        bet = b[j] - a[j] * gam[j]
        if bet == 0.0 :
            return u
        u[j] = (r[j]-a[j]*u[j-1]) / bet

    for j in range(n-2,j,-1) : 
        u[j] -= gam[j+1]*u[j+1]

    u[0] -= gam[1]*u[1]

    return u



def reduce_Householder(A, d, e):
    """ A is a matrix and d and e are vectors
    """
    n = len(d)
    for i in range(n-1, 0, -1):
        l = i - 1
        h = 0.0
        scale = 0.0
        if l > 0:
            for k in range(l+1):
                scale += abs(A[i][k])
            if scale == 0.0:
                e[i] = A[i][l]
            else:
                for k in range(l+1):
                    A[i][k] /= scale
                    h += A[i][k] * A[i][k]
                f = A[i][l]
                if f < 0.0:
                    g = math.sqrt(h)
                else:
                    g = -math.sqrt(h)
                e[i] = scale * g
                h -= f * g
                A[i][l] = f - g
                f = 0.0
                for j in range(l+1):
                    A[j][i] = A[i][j] / h
                    g = 0.0
                    for k in range(j+1):
                       g += A[j][k] * A[i][k]
                    for k in range(j+1, l+1):
                       g += A[k][j] * A[i][k]
                    e[j] = g / h
                    f += e[j] * A[i][j]
                hh = f / (h + h)
                for j in range(l+1):
                    f = A[i][j]
                    e[j] = g = e[j] - hh * f
                    for k in range(j+1):
                        A[j][k] -= f * e[k] + g * A[i][k]
        else:
            e[i] = A[i][l]
        d[i] = h
    d[0] = 0.0
    e[0] = 0.0
    for i in range(n):
        if d[i] != 0.0:
            for j in range(i):
                g = 0.0
                for k in range(i):
                    g += A[i][k] * A[k][j]
                for k in range(i):
                    A[k][j] -= g * A[k][i]
        d[i] = A[i][i]
        A[i][i] = 1.0
        for j in range(i):
            A[j][i] = A[i][j] = 0.0

def pythag(a, b):
    absa = abs(a)
    absb = abs(b)
    if absa > absb:
        ratio = absb / absa
        return absa * math.sqrt(1 + ratio**2)
    else:
        if absb == 0.0:
            return 0.0
        else:
            ratio = absa / absb
            return absb * math.sqrt(1 + ratio**2)

def sign(a, b):
    if b >= 0.0:
        if a >= 0.0:
            return a
        else:
            return -a
    else:
        if a >= 0.0:
            return -a
        else:
            return a

def solve_TQLI(d, e, z):
    """
    """
    n = len(d)
    for i in range(1, n):
        e[i-1] = e[i]
    e[n-1] = 0.0
    for l in range(n):
        iteration = 0
        m = l
        while True:
            for m_loop in range(l, n-1):
                m = m_loop
                dd = abs(d[m]) + abs(d[m+1])
                if (abs(e[m]) + dd) == dd:
                    break
                else:
                    m += 1
            if m != l:
                if iteration == 30:
                    raise Exception("Too many iterations in solve_TQLI")
                iteration += 1
                g = (d[l+1] - d[l]) / (2.0 * e[l])
                r = pythag(g, 1.0)
                g = d[m] - d[l] + e[l] / (g + sign(r, g))
                s = 1.0
                c = 1.0
                p = 0.0
                for i in range(m-1, l-1, -1):
                    f = s * e[i]
                    b = c * e[i]
                    e[i+1] = r = pythag(f, g)
                    if r == 0.0:
                        d[i+1] -= p
                        e[m] = 0.0
                        break
                    s = f / r
                    c = g / r
                    g = d[i+1] - p
                    r = (d[i] - g) * s + 2.0 * c * b
                    p = s * r
                    d[i+1] = g + p
                    g = c * r - b
                    for k in range(n):
                        f = z[k][i+1]
                        z[k][i+1] = s * z[k][i] + c * f
                        z[k][i] = c * z[k][i] - s * f
                if r == 0.0 and i >= l:
                    continue
                d[l] -= p
                e[l] = g
                e[m] = 0.0
            if m == l:
                break

def sort_eigen(d, z, ascending=False):
    """ Sorts eigenvalues and eigenvectors in descending order.
    """
    n = len(d)
    if len(z) != n or len(z[0]) != n:
        raise Exception("Bad vector, matrix dimensions in sort_eigen")
    for i in range(n-1):
        k = i
        p = d[k]
        for j in range(i, n):
            if d[j] >= p:
                k = j
                p = d[k]
        if k != i:
            d[k] = d[i]
            d[i] = p
            for j in range(n):
                p = z[j][i]
                z[j][i] = z[j][k]
                z[j][k] = p

def solve_eigen_symmetric(A):
    """ Returns tuple d, C where e is the vector of eigenvalues and
        C is the matrix of eigenvectors.
        Eigenvalues are sorted in decreasing order.
    """
    n = len(A)
    d = [ 0.0 ] * n
    e = [ 0.0 ] * n
    C = Matrix_copy(A)
    reduce_Householder(C, d, e)
    solve_TQLI(d, e, C)
    #sort_eigen(d, C)
    return d, C

def solve_eigen_generalized(A, S):
    """ Returns tuple e, C where e is the vector of eigenvalues and
        C is the matrix of eigenvectors.
        Eigenvalues are sorted in decreasing order.
        Uses Householder and solve_TQLI
    """
    n = len(A)
    V = Matrix_copy(S)
    s, V = solve_eigen_symmetric(V)
    for i in range(n):
        for j in range(n):
            V[i][j] /= math.sqrt(s[j])
    VT = Matrix_transpose(V)
    C = Matrix_multiply( A, V)
    C = Matrix_multiply( VT, C )
    e, C = solve_eigen_symmetric(C)
    C = Matrix_multiply(V, C)
    return e, C


def find_maximum( a, b, f, tol ) :
    negf = lambda x : -1.0*f(x)
    xmin = scipy.optimize.minimize( fun=negf, x0=a, tol=tol )
    fxmin = negf(xmin.x)
    return [xmin.x,-1*fxmin]



def find_minimum( a, b, f, tol ) :
    xmin = scipy.optimize.minimize( fun=f, x0=a, tol=tol )
    fxmin = f(xmin.x)
    return [xmin.x, fxmin]



def line_search(xold, fold, g, p, stpmax, func):
    """ Returns x, f, check
        Input : xold is a list of n floats representing a point in n-dimensions
                fold is the value of the function at xold
                g is a list of n floats representing the gradient at xold
                p is a list of n floats representing a direction
                stpmax limits length of steps
                func is name of function
        Output: x is a list of n floats representing the new point along p
                f is the value of func(x)
                check is False on normal exit, True if x is too close to xold
    """
    ALF = 1.0e-4 ; TOLX = 2.22e-16
    alam2 = 0.0 ; f2 = 0.0
    n = len(xold)
    x = list(xold)
    check = False
    s = sum(p[i]**2 for i in range(n))
    s = math.sqrt(s)
    if s > stpmax:
        for i in range(n):
            p[i] *= stpmax / s
    slope = sum(g[i] * p[i] for i in range(n))
    if slope >= 0.0:
        raise Exception("Roundoff problem in line_search")
    test = 0.0
    for i in range(n):
        temp = abs(p[i]) / max(abs(xold[i]), 1.0)
        if temp > test:
            test=temp
    alamin = TOLX / test
    alam = 1.0
    while True:
        for i in range(n):
            x[i] = xold[i] + alam * p[i]
        f = func(x)
        if alam < alamin:
            for i in range(n):
                x[i] = xold[i]
            check = True
            break
        elif f <= fold + ALF * alam * slope:
            break
        else:
            if alam == 1.0:
                tmplam = -slope / (2.0 * (f - fold - slope))
            else:
                rhs1 = f - fold - alam * slope
                rhs2 = f2 - fold - alam2 * slope
                a = (rhs1 / alam**2 - rhs2 / alam2**2) / (alam - alam2)
                b = (- alam2 * rhs1 / alam**2 +
                       alam * rhs2 / alam2**2) / (alam - alam2)
                if a == 0.0:
                    tmplam = -slope / (2.0 * b)
                else:
                    disc = b**2 - 3.0 * a * slope
                    if disc < 0.0:
                        tmplam = 0.5 * alam
                    elif b <= 0.0:
                        tmplam = ( - b + math.sqrt(disc) ) / (3.0 * a)
                    else:
                        tmplam = - slope / (b + math.sqrt(disc))
                if tmplam > 0.5 * alam:
                     tmplam = 0.5 * alam;
        alam2 = alam
        f2 = f
        alam = max(tmplam, 0.1 * alam)
    return x, f, check

def minimize_BFGS(p, gtol, func, dfunc):
    """ Returns iterations, fret
        Inputs: p is list of n floats representing a point in n-dimensions
                gtol is float convergence required for zeroing dfunc
                func is function of p
                dfunc is computes the gradient
        Output: p replaced by minimum point
                iterations is number of iterations used
                fret is value of func at p
        Uses the Broyden-Fletcher-Goldfarb-Shanno algorithm to minimize func
    """
    ITMAX = 200 ; EPS = 2.22e-16 ; TOLX = 4 * EPS ; STPMX = 100.0
    n = len(p)
    g = [ 0.0 ] * n
    dg = list(g) ; hdg = list(g) ; pnew = list(g) ; xi = list(g)
    hessin = Matrix(n, n)

    fp = func(p)
    dfunc(p, g)
    svm = 0.0
    for i in range(n):
        for j in range(n):
            hessin[i][j] = 0.0
        hessin[i][i] = 1.0
        xi[i] = -g[i]
        svm += p[i]**2
    stpmax = STPMX * max(math.sqrt(svm), float(n))
    for its in range(ITMAX):
        iterations = its
        pnew, fret, check = line_search(p, fp, g, xi, stpmax, func)
        fp = fret
        for i in range(n):
            xi[i] = pnew[i] - p[i]
            p[i] = pnew[i]
        test = 0.0
        for i in range(n):
            temp = abs(xi[i] ) / max(abs(p[i]), 1.0)
            if temp > test:
                test=temp
        if test < TOLX:
            return iterations, fret
        for i in range(n):
            dg[i] = g[i]
        dfunc(p, g)
        test = 0.0;
        den = max(fret, 1.0)
        for i in range(n):
            temp = abs(g[i]) * max(abs(p[i]), 1.0) / den
            if temp > test:
                test = temp
        if test < gtol:
            return iterations, fret
        for i in range(n):
            dg[i] = g[i] - dg[i]
        for i in range(n):
            hdg[i] = 0.0
            for j in range(n):
                hdg[i] += hessin[i][j] * dg[j]
        fac = fae = sumdg = sumxi = 0.0
        for i in range(n):
            fac += dg[i] * xi[i]
            fae += dg[i] * hdg[i]
            sumdg += dg[i] * dg[i]
            sumxi += xi[i]**2
        if fac > math.sqrt(EPS * sumdg * sumxi):
            fac = 1.0 / fac
            fad = 1.0 / fae
            for i in range(n):
                dg[i] = fac * xi[i] - fad * hdg[i]
            for i in range(n):
                for j in range(i, n):
                    hessin[i][j] += (fac * xi[i] * xi[j] -
                        fad * hdg[i] * hdg[j] + fae * dg[i] * dg[j] )
                    hessin[j][i] = hessin[i][j]
        for i in range(n):
            xi[i] = 0.0
            for j in range(n):
                xi[i] -= hessin[i][j] * g[j]
    raise Exception("Too many iterations in minimize_BFGS")

def least_squares_fit(x, y):
    """ Fits y[i] = a + b * x[i] using least squares algorithm
        Returns a, b, sigma_a, sigma_b, sigma
        a = intercept, b = slope, sigma_a,b = estimated error in a,b
        sigma = estimated error bar on y
    """
    n = len(x)
    s_x  = sum(x)
    s_y  = sum(y)
    s_xx = sum(x_i**2 for x_i in x)
    s_xy = sum(x[i]*y[i] for i in range(n))
    denom = n * s_xx - s_x**2
    a = (s_xx * s_y - s_x * s_xy) / denom
    b = (n * s_xy - s_x * s_y) / denom
    variance = sum((y[i] - (a + b*x[i]))**2 for i in range(n))
    sigma = math.sqrt(variance/(n-2))
    sigma_a = math.sqrt(sigma**2 * s_xx / denom)
    sigma_b = math.sqrt(sigma**2 * n / denom)
    return [a, b, sigma_a, sigma_b, sigma]

def chi_square_fit(x, y, err):
    """ Fits y[i] = a + b * x[i] using chi square algorithm
        err[i] is the error in y[i]
        Returns a, b, sigma_a, sigma_b, chi_square
        a = intercept, b = slope, sigma_a,b = estimated error in a,b
    """
    n = len(x)
    S = sum(1/err[i]**2 for i in range(n))
    S_x = sum(x[i]/err[i]**2 for i in range(n))
    S_y = sum(y[i]/err[i]**2 for i in range(n))
    t = [(x[i] - S_x/S) / err[i] for i in range(n)]
    S_tt = sum(t_i**2 for t_i in t)
    b = sum(t[i]*y[i]/err[i] for i in range(n)) / S_tt
    a = (S_y - S_x * b) / S
    sigma_a = math.sqrt((1 + S_x**2/S/S_tt) / S)
    sigma_b = math.sqrt(1/S_tt)
    chi_square = sum(((y[i] - a - b*x[i]) / err[i])**2 for i in range(n))
    return(a, b, sigma_a, sigma_b, chi_square)


def Euler_step( p, dt, acceleration ) :
    t = p[0]
    theta = p[1]
    dtheta_dt = p[2]
    dtheta_dt_old = dtheta_dt

    #apply Euler's algorithm
    dtheta_dt += acceleration(p) * dt
    theta += dtheta_dt_old * dt
    t += dt
    p[0] = t
    p[1] = theta
    p[2] = dtheta_dt
    

def RK4_step(x, dt, flow):
    """replaces x(t) by x(t + dt) using fourth order Runge-Kutta
    with derivative vector flow
    """
    n = len(x)
    k1 = [ dt * k for k in flow(x) ]
    x_temp = [ x[i] + k1[i] / 2.0 for i in range(n) ]
    k2 = [ dt * k for k in flow(x_temp) ]
    x_temp = [ x[i] + k2[i] / 2.0 for i in range(n) ]
    k3 = [ dt * k for k in flow(x_temp) ]
    x_temp = [ x[i] + k3[i] for i in range(n) ]
    k4 = [ dt * k for k in flow(x_temp) ]
    for i in range(n):
        x[i] += (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0

def RK4_adaptive_step(x, dt, flow, accuracy=1e-6):
    """returns the adaptive step size that ensures the requested accuracy
    and replaces x(t) by x(t + dt_adapted)
    uses the algorithm in Numerical Recipes to adapt the step size
    returns the adapted step that ensures the request accuracy
    """
    SAFETY = 0.9; PGROW = -0.2; PSHRINK = -0.25;
    ERRCON = 1.89E-4; TINY = 1.0E-30
    n = len(x)
    scale = flow(x)
    scale = [ abs(x[i]) + abs(scale[i] * dt) + TINY for i in range(n) ]
    while True:
        dt /= 2.0
        x_half = [ x[i] for i in range(n) ]
        RK4_step(x_half, dt, flow)
        RK4_step(x_half, dt, flow)
        dt *= 2.0
        x_full = [ x[i] for i in range(n) ]
        RK4_step(x_full, dt, flow)
        Delta = [ x_half[i] - x_full[i] for i in range(n) ]
        error = max( abs(Delta[i] / scale[i]) for i in range(n) ) / accuracy
        if error <= 1.0:
            break;
        dt_temp = SAFETY * dt * error**PSHRINK
        if dt >= 0.0:
            dt = max(dt_temp, 0.1 * dt)
        else:
            dt = min(dt_temp, 0.1 * dt)
        if abs(dt) == 0.0:
            raise OverflowError("step size underflow")
    if error > ERRCON:
        dt *= SAFETY * error**PGROW
    else:
        dt *= 5.0
    for i in range(n):
        x[i] = x_half[i] + Delta[i] / 15.0
    return dt

def RK4_integrate(x, dt, flow, Delta_t, accuracy=1e-6):
    """replaces x(t) by x(t+Delta_t) using adaptive Runge-Kutta with
    suggested initial step size dt and returns the adapted step size
    assumes that t = x[0]
    """
    if dt * Delta_t <= 0.0:
        raise Exception("dt * Delta_t <= 0.0")
    if abs(dt) > abs(Delta_t):
        dt = Delta_t / 5.0
    t0 = x[0]
    while abs(x[0] - t0) < abs(Delta_t):
        dt = RK4_adaptive_step(x, dt, flow, accuracy)
    RK4_step(x, t0 + Delta_t - x[0], flow)
    return dt



def root_print_header(algorithm, accuracy):
    sys.stdout.write("\n ROOT FINDING using " + algorithm +
                     "\n Requested accuracy = " +repr(accuracy) +
                     "\n Step     Guess For Root          Step Size      " +
                     "     Function Value" +
                     "\n ----  --------------------  --------------------" +
                     "  --------------------" + "\n")

def root_print_step(step, x, dx, f_of_x):
    sys.stdout.write(repr(step).rjust(5))
    for val in [x, dx, f_of_x]:
        sys.stdout.write("  " + repr(val).ljust(20))
    sys.stdout.write("\n")

def root_max_steps(algorithm, max_steps):
    raise Exception(" " + algorithm + ": maximum number of steps " +
                    repr(max_steps) + " exceeded\n")

def root_simple(f, x, dx, accuracy=1.0e-6, max_steps=1000, root_debug=False):
    """Return root of f(x) given guess x and step dx with specified accuracy.
    Step must be in direction of root: dx must have same sign as (root - x).
    """
    f0 = f(x)
    fx = f0
    step = 0
    if root_debug:
        root_print_header("Simple Search with Step Halving", accuracy)
        root_print_step(step, x, dx, f0)
    while abs(dx) > abs(accuracy) and f0 != 0.0:
        x += dx
        fx = f(x)
        if f0 * fx < 0.0:   # stepped past root
            x -= dx         # step back
            dx /= 2.0       # use smaller step
        step += 1
        if step > max_steps:
            root_max_steps("root_simple", max_steps)
        if root_debug:
            root_print_step(step, x, dx, fx)
    return x

def root_bisection(f, x1, x2, accuracy=1.0e-6, max_steps=1000, root_debug=False):
    """Return root of f(x) in bracketed by x1, x2 with specified accuracy.
    Assumes that f(x) changes sign once in the bracketed interval.
    Uses bisection root-finding algorithm.
    """
    f1 = f(x1)
    f2 = f(x2)
    if f1 * f2 > 0.0:
        raise Exception("f(x1) * f(x2) > 0.0")
    x_mid = (x1 + x2) / 2.0
    f_mid = f(x_mid)
    dx = x2 - x1
    step = 0
    if root_debug:
        root_print_header("Bisection Search", accuracy)
        root_print_step(step, x_mid, dx, f_mid)
    while abs(dx) > accuracy:
        if f_mid == 0.0:
            dx = 0.0
        else:
            if f1 * f_mid > 0:
                x1 = x_mid
                f1 = f_mid
            else:
                x2 = x_mid
                f2 = f_mid
            x_mid = (x1 + x2) / 2.0
            f_mid = f(x_mid)
            dx = x2 - x1
        step += 1
        if step > max_steps:
            warning = "Too many steps (" + repr(step) + ") in root_bisection"
            raise Exception(warning)
        if root_debug:
            root_print_step(step, x_mid, dx, f_mid)
    return x_mid

def root_secant(f, x0, x1, accuracy=1.0e-6, max_steps=20,root_debug=False):
    """Return root of f(x) given guesses x0 and x1 with specified accuracy.
    Uses secant root-finding algorithm.
    """
    f0 = f(x0)
    dx = x1 - x0
    step = 0
    if root_debug:
        root_print_header("Secant Search", accuracy)
        root_print_step(step, x0, dx, f0)
    if f0 == 0:
        return x0
    while abs(dx) > abs(accuracy):
        f1 = f(x1)
        if f1 == 0:
            return x1
        if f1 == f0:
            raise Exception("Secant horizontal f(x0) = f(x1) algorithm fails")
        dx *= - f1 / (f1 - f0)
        x0 = x1
        f0 = f1
        x1 += dx
        step += 1
        if step > max_steps:
            root_max_steps("root_secant", max_steps)
        if root_debug:
            root_print_step(step, x1, dx, f1)
    return x1

def root_tangent(f, fp, x0, accuracy=1.0e-6, max_steps=20, root_debug=False):
    """Return root of f(x) with derivative fp = df(x)/dx
    given initial guess x0, with specified accuracy.
    Uses Newton-Raphson (tangent) root-finding algorithm.
    """
    f0 = f(x0)
    fp0 = fp(x0)
    if fp0 == 0.0:
        raise Exception(" root_tangent df/dx = 0 algorithm fails")
    dx = - f0 / fp0
    step = 0
    if root_debug:
        root_print_header("Tangent Search", accuracy)
        root_print_step(step, x0, dx, f0)
    if f0 == 0.0:
        return x0
    while True:
        fp0 = fp(x0)
        if fp0 == 0.0:
            raise Exception(" root_tangent df/dx = 0 algorithm fails")
        dx = - f0 / fp0
        x0 += dx
        f0 = f(x0)
        if abs(dx) <= accuracy or f0 == 0.0:
            return x0
        step += 1
        if step > max_steps:
            root_max_steps("root_tangent", max_steps)
        if root_debug:
            root_print_step(step, x0, dx, f0)
    return x0




def polyfit( x, y, sigma, M ) :
    '''
    Adapted from example code in 
    Garcia, "Numerical Methods for Physics", Second Edition.
    Download available at : 
    http://www.mathworks.com/matlabcentral/fileexchange/2270-numerical-methods-for-physics-2e/all_files
    This adapts the usage in two ways : 


    inputs :
      * x : dependent variable (data x values 0... n_data-1)
      * y : independent variable (data y values 0...n_data-1)
      * sigma : uncertainties on y (0...n_data-1)
      * M : order of polynomial to fit

    returns :
      * a_fit : coefficients of the polynomial fit (0...M-1)
      * sig_a : uncertainties on a(0...M-1)
      * yy    : fitted y values (0...n_data-1)
      * chisqr: chisquared value
    '''

    # Form the vector b and design matrix A
    N = len(x)
    b = [0.0] * N
    A = Matrix(N,M)
    a_fit = [0.0] * M
    sig_a = [0.0] * M
    yy = [0.0] * N
    chisq = 0.0
    
    for i in xrange(N) :
        b[i] = y[i] / sigma[i]
        for j in xrange(M) :
            if sigma[i] > 0.0 :
                A[i][j] = pow(x[i],float(j))/sigma[i]
            else :
                A[i][j] = 0.0

    # Compute the correlation matrix C
    C = Matrix(M,M)
    Cinv = Matrix(M,M)
    for i in xrange(M) :
        for j in xrange(M) :
            Cinv[i][j] = 0.0
            for k in xrange(N) :
                Cinv[i][j] += A[k][i] * A[k][j]


    [determinant, C] = inverse( Cinv )
    print ' Cinv is : '
    print Cinv
    print ' C is : '
    print C
 
    # Compute the least squares polynomial coefficients a_fit
    for k in xrange(M) :
        a_fit[k] = 0.0
        for j in xrange(M) :
            for i in xrange(N) :
                a_fit[k] += C[k][j] * A[i][j] * b[i]

    # Compute the estimated error bars for the coefficients
    for j in xrange(M) :
        sig_a[j] = math.sqrt( C[j][j] )
        chisqr = 0.0
        for i in xrange(N) :
            yy[i] = 0.0
            for j in xrange(M) :
                yy[i] += a_fit[j] * pow( x[i], float(j) )
            delta = (y[i] - yy[i]) / sigma[i]
            chisqr += delta**2

    return [ a_fit, sig_a, yy, chisq]


def inverse(A) :
    '''
    # Input
    #    A    -    Matrix A (N by N)
    # Outputs
    #   Ainv  -    Inverse of matrix A (N by N)
    #  determ -    Determinant of matrix A (return value)
    '''


    N = len(A)
    M = len(A[0])
    if N != M :
        return None
  
    Ainv = copy.deepcopy(A)
    
    scale = [0.0] * N
    b = Matrix(N,N)
    index = [0] * N


    # Matrix b is initialized to the identity matrix
    for i in xrange(N) :
        b[i][i] = 1.0

    # Set scale factor, scale(i) = max( |a(i,j)| ), for each row
    for i in xrange(N) :
        index[i] = i                         # Initialize row index list
        scalemax = 0.
        for j in xrange(N) :
            imax = abs(A[i][j])
            if imax > scalemax :
                scalemax = imax
        scale[i] = scalemax

    # Loop over rows k = 0, ..., (N-2)
    signDet = 1
    for k in xrange(N-1) :
        #* Select pivot row from max( |a(j,k)/s(j)| )
        ratiomax = 0.0
        jPivot = k
        for i in range(k,N) : 
            ratio = abs(A[index[i]][k])/scale[index[i]]
            if ratio > ratiomax :
                jPivot = i
                ratiomax = ratio

        #* Perform pivoting using row index list
        indexJ = index[k]
        if  jPivot != k :               # Pivot
            indexJ = index[jPivot]
            index[jPivot] = index[k]   # Swap index jPivot and k
            index[k] = indexJ
            signDet *= -1                          # Flip sign of determinant

        #* Perform forward elimination
        for i in range(k+1,N) : 
            coeff = A[index[i]][k]/A[indexJ][k]
            for j in range(k+1,N) : 
                A[index[i]][j] -= coeff*A[indexJ][j]
            A[index[i]][k] = coeff
            for j in xrange(N) :
                b[index[i]][j] -= A[index[i]][k]*b[indexJ][j]

    # Compute determinant as product of diagonal elements
    determ = signDet         # Sign of determinant
    for i in xrange(N) : 
        determ *= A[index[i]][i]

    # Perform backsubstitution
    for k in xrange(N) : 
        Ainv[N-1][k] = b[index[N-1]][k]/A[index[N-1]][N-1]
        for i in range(N-1,-1,-1):
            sumi = b[index[i]][k]
            for j in range(i+1,N) :
                sumi -= A[index[i]][j]*Ainv[j][k]
            Ainv[i][k] = sumi / A[index[i]][i]
           

    return [determ, Ainv]
