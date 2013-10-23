from cpt import *

def V(x) : 
    return -pow(x, 2.0)/ 2 + pow(x, 4.0)/4



acc = 1e-6
guess1 = 0.1
guess2 = -0.1
[x_min, V_min] = find_minimum(guess1, guess2, V, acc)
print "V(x) = " + str(V_min)+  " at x = " + str(x_min)

[x_max, V_max] = find_maximum (guess1, guess2, V, acc)
print "V(x) = " + str(V_max)+  " at x = " + str(x_max)
