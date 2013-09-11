import math
import matplotlib.pyplot as plt

# define a function to linear fit x-y data without error bars
def least_squares_fit(x, y):
    """Perform a least-squares fit to data (x,y)

    Args :
       x : x values
       y : y values

    Returns :
       a : intercept
       b : slope
       sigma : total uncertainty (sqrt(variance/(n-2)))
       sigma_a : uncertainty on a
       sigma_b : uncertainty on b

    """
    n = len(x)
    s_x  = sum(x)
    s_y  = sum(y)
    s_xx = sum(x_i**2 for x_i in x)
    s_xy = sum(x[i]*y[i] for i in range(n))
    denom = n * s_xx - s_x**2
    if abs(denom) > 0.00001 : 
        a = (s_xx * s_y - s_x * s_xy) / denom
        b = (n * s_xy - s_x * s_y) / denom
        variance = sum((y[i] - (a + b*x[i]))**2 for i in range(n))
        sigma = math.sqrt(variance/(n-2))
        sigma_a = math.sqrt(sigma**2 * s_xx / denom)
        sigma_b = math.sqrt(sigma**2 * n / denom)
        return [a, b, sigma, sigma_a, sigma_b]
    else :
        print 'error : divided by zero!'
        return None



"""Plot data for the Gutenberg-Richter Model.

Here, we plot the curve of the number of earthquakes
greater than magnitude M, for each M value.

So, we loop over the earthquakes, and store the
frequency of each magnitude. At the end of the loop,
we compute the cumulative distribution such that the
value at magnitude M will be the integral of the frequency
distribution for >= M. This is what the Gutenberg-Richter
Model predicts. 

"""


# data downloaded from http://earthquake.usgs.gov/earthquakes/search/
print ' Earthquake data: Gutenberg-Richter Model'
data_file_name = 'california_earthquakes_2010_to_2013.csv'
file = open(data_file_name, 'r')
lines = file.readlines()
file.close()
print ' read', len(lines), 'lines from', data_file_name

# store event data in two ways for demonstration
# 1. in a python dictionary object with
#    key = magnitude starting in column 50
#    value = number of events with this magnitude
# 2. keeping a "tuple" of the results for later usage with matplotlib
histogram = dict()
magvalues = []
for line in lines:
    if line[0] != 't' :
        try:
            words = line.split(',')
            [latitude,longitude,depth,mag] = [float(s) for s in words[1:5] ]
            magvalues.append( mag )
            # For debugging : 
            #print 'read : {0:s} , ({1:6.3f},{2:6.3f}), d = {3:4.1f}, m = {4:6.3f}'.format(
            #    words[0], latitude, longitude, depth, mag
            #    )
            histogram[mag] += 1
        except KeyError : 
            histogram.setdefault(mag, 1)
        except ValueError:
            print 'bad data:', line

num_events = sum(histogram[M] for M in histogram.keys())
num_bins = len(histogram)
print ' stored', num_events, 'events in', num_bins, 'bins'

# x data = M values sorted in increasing order
# y data = log_10(N) where N = number of events with magnitude >= M
M_values = sorted(histogram.keys())
dN_values = [histogram[M] for M in M_values]
log10N_values = [ math.log10(sum(dN_values[i:]))
                  for i in range(len(M_values)) ]

print 'log10N_values is '
for i in xrange(len(log10N_values)) :
    print str(M_values[i]) + '  ' + str(log10N_values[i])
# First plot our "home-grown" values with our least-squares
# fit in place. 
plt.subplot( 2, 1, 1)
plt.plot( M_values, log10N_values, 'v')
plt.xlabel( 'Magnitude (M)' )
plt.ylabel( 'log(N)' )

# Next also plot matplotlib's version of the same thing.
plt.subplot( 2, 1, 2)
plt.hist( magvalues, bins=90, range=[1.0,10.0], log=True, bottom=0.1,cumulative=-1)
plt.xlabel( 'Magnitude (M)' )
plt.ylabel( 'N' )

# perform a least square fit
fit = least_squares_fit(M_values, log10N_values)
print ' least_squares fit to data:'
print ' slope b = {0:6.3f} +- {1:6.3f}'.format( fit[1], fit[4])
print ' intercept a = {0:6.3f} +- {1:6.3f}'.format( fit[0], fit[3])
print ' log_10(N) error bar = {0:6.3f}'.format( fit[2] )

plt.show()
