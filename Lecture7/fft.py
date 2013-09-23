from cmath import exp, pi
from math import sin, cos
import matplotlib.pyplot as plt

import numpy
from numpy.fft import fft
from numpy import array, conj, divide



def discrete_transform(data):
    """Return Discrete Fourier Transform (DFT) of a complex data vector"""
    N = len(data)
    transform = [ 0 ] * N
    for k in range(N):
        for j in range(N):
            angle = 2 * pi * k * j / N
            transform[k] += data[j] * exp(1j * angle)
    return transform

def fft(x):
    N = len(x)
    if N <= 1: return x
    elif N % 2 == 1:         # N is odd, lemma does not apply
        print 'N is ' + str(N) + ', fall back to discrete transform'
        return discrete_transform(x)
    even = fft(x[0::2])
    odd =  fft(x[1::2])
    return [even[k] + exp(-2j*pi*k/N)*odd[k] for k in xrange(N/2)] + \
           [even[k] - exp(-2j*pi*k/N)*odd[k] for k in xrange(N/2)]


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
    
