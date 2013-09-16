import matplotlib.pyplot as plt

from fft import fft, fft_power, ifft
from numpy import array, real
import math
import time


# data downloaded from http://www.sidc.be/DATA/yearssn.dat
print ' Sunspot data from Solar Influences Data Analysis Center'
data_file_name = 'monthssn.dat'
file = open(data_file_name, 'r')
lines = file.readlines()
file.close()
print ' read', len(lines), 'lines from', data_file_name

window = False

yinput = []
xinput = []

for line in lines :
    if line[0] != '#' :
        try:
            words = line.split()
            xval = float(words[1])
            yval = float( words[2] )
            yinput.append( yval )
            xinput.append( xval )
        except ValueError :
            print 'bad data:',line


N = len(yinput)
log2N = math.log(N, 2)
if log2N - int(log2N) > 0.0 :
    print 'Padding with zeros!'
    pads = [0.0] * (pow(2, int(log2N)+1) - N)
    yinput = yinput + pads
    N = len(yinput)
    print 'Padded : '
    print len(yinput)
    # Apply a window to reduce ringing from the 2^n cutoff
    if window : 
        for iy in xrange(len(yinput)) :
            yinput[iy] = yinput[iy] * (0.5 - 0.5 * math.cos(2*math.pi*iy/float(N-1)))

y = array( yinput ) 
x = array([ float(i) for i in xrange(len(y)) ] )
Y = fft(y)

maxfreq = 50
# Now smooth the data
for iY in range(maxfreq, len(Y)-maxfreq ) :
    Y[iY] = complex(0,0)
    #Y[iY] = Y[iY] * (0.5 - 0.5 * math.cos(2*math.pi*iY/float(N-1))) 

    #for iY in range(0,N) : 
    #    Y[iY] = Y[iY] * math.exp(-1.0*iY / 50.0)

powery = fft_power(Y)
powerx = array([ float(i) for i in xrange(len(powery)) ] )

Yre = [math.sqrt(Y[i].real**2+Y[i].imag**2) for i in xrange(len(Y))]

ysmoothed = ifft(Y)
ysmoothedreal = real(ysmoothed)

ax1 = plt.subplot(2, 1, 1)
p1, = plt.plot( x, y )
p2, = plt.plot( x, ysmoothedreal )
ax1.legend( [p1,p2], ['Original', 'Smoothed'] )

ax2 = plt.subplot(2, 1, 2)
p3, = plt.plot( powerx, powery )
p4, = plt.plot( x, Yre )
ax2.legend( [p3, p4], ["Power", "Magnitude"] )
plt.yscale('log')


plt.show()
