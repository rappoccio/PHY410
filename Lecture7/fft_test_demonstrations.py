import matplotlib.pyplot as plt

from fft import fft
from fft import fft_power
from numpy import array
import math
import time

plotfirst = True

if plotfirst == True : 
    # make some fake data as a single-frequency sinusoid

    N = 1024      # number of samples
    f = 10.0      # frequency / N
    m = 0.00      # linear slope, if desired

    # For demonstrations : 
    clip = 0      # "clip" so we can pad with zeros
    window = False # window or not? 
    
    
    x = array([ float(i) for i in xrange(N) ] )
    if window : 
        y = array([ math.sin(-2*math.pi*f* xi / float(N)) *(0.5 - 0.5 * math.cos(2*math.pi*xi/float(N-1)))  + m*xi  for xi in x ])
    else :
        y = array([ math.sin(-2*math.pi*f* xi / float(N)) + m*xi  for xi in x ])

    for i in range(N-clip, N) :
        y[i] = 0.0

    start_time = time.time()
    Y = fft(y)
    print time.time() - start_time, " seconds" 

    powery = fft_power(Y)
    powerx = array([ float(i) for i in xrange(len(powery)) ] )


    Yre = [math.sqrt(Y[i].real**2 + Y[i].imag**2) for i in xrange(N)]

    s1 = plt.subplot(2, 1, 1)
    plt.plot( x, y )

    s2 = plt.subplot(2, 1, 2)
    s2.set_autoscalex_on(False)
    p1, = plt.plot( powerx, powery )
    p2, = plt.plot( x, Yre )
    s2.legend( [p1, p2], ["Power", "Magnitude"] )
    plt.xlim([0,N/4])
    plt.yscale('log')
    plt.show()


else : 
    # data downloaded from ftp://ftp.cmdl.noaa.gov/ccg/co2/trends/co2_mm_mlo.txt
    print ' C02 Data from Mauna Loa'
    data_file_name = 'co2_mm_mlo.txt'
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
                xval = float(words[2])
                yval = float( words[4] )
                yinput.append( yval )
                xinput.append( xval )
            except ValueError :
                print 'bad data:',line


    N = len(yinput)
    log2N = math.log(N, 2)
    if log2N - int(log2N) > 0.0 :
        print 'Padding with zeros!'
        pads = [300.0] * (pow(2, int(log2N)+1) - N)
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

    powery = fft_power(Y)
    powerx = array([ float(i) for i in xrange(len(powery)) ] )

    Yre = [math.sqrt(Y[i].real**2+Y[i].imag**2) for i in xrange(len(Y))]


    plt.subplot(2, 1, 1)
    plt.plot( x, y )

    ax = plt.subplot(2, 1, 2)
    p1, = plt.plot( powerx, powery )
    p2, = plt.plot( x, Yre )
    ax.legend( [p1, p2], ["Power", "Magnitude"] )
    plt.yscale('log')


    plt.show()
