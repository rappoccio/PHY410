# pendul - Program to compute the motion of a simple pendulum
# using the Euler or Verlet method

from math import *

#* Select the numerical method to use: Euler or Verlet
method = input( "Choose a numerical method : Euler (0) or Verlet (1): " )

					   
#* Set initial position and velocity of pendulum
theta0 = input( "Enter initial angle (in degrees): " )
theta = theta0*pi/180   # Convert angle to radians
omega = 0.0             # Set the initial velocity

#* Set the physical constants and other variables
g_over_L = 1.0          # The constant g/L
time = 0.0              # Initial time
time_old = 0.0          # Time of previous reversal
irev = 0                   # Used to count number of reversals
tau = input("Enter time step: ")


#* Take one backward step to start Verlet
accel = -g_over_L*sin(theta)    # Gravitational acceleration
theta_old = theta - omega*tau + 0.5*tau*tau*accel    

#* Loop over desired number of steps with given time step
#    and numerical method
nStep = input( "Enter number of time steps: " )

t_plot = []
th_plot = []
period = []

for iStep in xrange(nStep) :

    #* Record angle and time for plotting
    t_plot.append( time )            
    th_plot.append( theta*180/pi )   # Convert angle to degrees
    time += tau
  
    #* Compute new position and velocity using 
    #    Euler or Verlet method
    accel = -g_over_L*sin(theta)    # Gravitational acceleration
    if  method == 0 :
        theta_old = theta        # Save previous angle
        theta += tau*omega       # Euler method
        omega += tau*accel 
    else :
        theta_new = 2*theta - theta_old + tau*tau*accel
        theta_old = theta	    # Verlet method
        theta = theta_new  
    
  
    #* Test if the pendulum has passed through theta = 0
    #    if yes, use time to estimate period
    if theta*theta_old < 0  : # Test position for sign change
        print "Turning point at time t = " + str(time)
        if irev == 0 :          # If this is the first change,
            time_old = time       # just record the time
        else : 
            period.append( 2*(time - time_old) )
            time_old = time
      
        irev += 1       # Increment the number of reversals
    
  
    nPeriod = irev-1    # Number of times period is measured

#* Estimate period of oscillation, including error bar
AvePeriod = 0.0
ErrorBar = 0.0

for i in xrange(nPeriod) :
    AvePeriod += period[i]
  
AvePeriod /= nPeriod
for i in xrange(nPeriod) :
    ErrorBar += (period[i] - AvePeriod)*(period[i] - AvePeriod)
  
ErrorBar = sqrt(ErrorBar/(nPeriod*(nPeriod-1)))
print "Average period = " + str(AvePeriod) + " +/- " + str( ErrorBar )

import matplotlib.pyplot as plt
plt.plot( t_plot, th_plot )
plt.show()
