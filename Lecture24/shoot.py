import math
import cpt
import bvp


    
print " Steping solution of u'' = - pi^2 (u+1) / 4"
x0, x1, u0, u1 = input( "Input x0, x1, u0, u1 : ")
delta = input( "Input initial guess : ")

stepper = bvp.BVPStep( x0=x0, x1=x1, u0=u0, u1=u1, delta=delta, N= 100)
shooter = bvp.BVPShoot( step=stepper, x0=x0, x1=x1, u0=u0, u1=u1, delta=delta, N= 100)

slope = shooter.solve()
print ''
print " Shooting method found slope = ", slope
print ""
print ""
x_plot = []
u_plot = []
#print " Trajectory : "
stepper.xy = [ x0, u0, slope ]
#print '{0:>8s} {1:>8s} {2:>8s}'.format( 'x', 'u', "u'" )
#print '{0:8.4f} {1:8.4f} {2:8.4f}'.format( stepper.xy[0], stepper.xy[1], stepper.xy[2] )
x_plot.append( stepper.xy[0] )
u_plot.append( stepper.xy[1] )
for i in range(stepper.N):
    stepper.step()
    #print '{0:8.4f} {1:8.4f} {2:8.4f}'.format( stepper.xy[0], stepper.xy[1], stepper.xy[2] )
    x_plot.append( stepper.xy[0] )
    u_plot.append( stepper.xy[1] )


import matplotlib
matplotlib.rcParams['legend.fancybox'] = True
import matplotlib.pyplot as plt

plt.plot( x_plot, u_plot )

plt.show()
