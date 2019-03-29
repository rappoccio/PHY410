import math
import sys
from cpt import *

# Membrane constants from Table 3

C_M = 1.0           # membrane capacitance per unit area
V_Na = -115         # sodium Nernst potential
V_K = +12           # potassium Nernst potential
V_l = -10.613       # leakage potential
g_Na = 120          # sodium conductance
g_K = 36            # potassium conductance
g_l = 0.3           # leakage conductance

# Voltage-dependent rate constants (constant in time)

def alpha_n(V):
    return 0.01 * (V + 10) / (math.exp((V + 10) / 10) - 1)

def beta_n(V):
    return 0.125 * math.exp(V / 80)

def alpha_m(V):
    return 0.1 * (V + 25)  / (math.exp((V + 25) / 10) - 1)

def beta_m(V):
    return 4 * math.exp(V / 18)

def alpha_h(V):
    return 0.07 * math.exp(V / 20)

def beta_h(V):
    return 1 / (math.exp((V + 30) / 10) + 1)

# Membrane current as function of time

def I(t):
    # In a voltage clamp experiment I = 0, see page 522 of H-H article
    return 0
    # For propagated action potential see Eqs. (30,31) in the article

# Hodgkin-Huxley equations

def HH_equations(Vnmht):
    # returns flow vector given extended solution vector [V, n, m, h, t]
    V = Vnmht[0]
    n = Vnmht[1]
    m = Vnmht[2]
    h = Vnmht[3]
    t = Vnmht[4]
    flow = [0] * 5
    flow[0] = ( I(t) - g_K * n**4 * (V - V_K) - g_Na * m**3 * h * (V - V_Na) -
                g_l * (V - V_l) ) / C_M
    flow[1] = alpha_n(V) * (1 - n) - beta_n(V) * n
    flow[2] = alpha_m(V) * (1 - m) - beta_m(V) * m
    flow[3] = alpha_h(V) * (1 - h) - beta_h(V) * h
    flow[4] = 1
    return flow

print(" Hodgkin-Huxley Fig. 12")
# resting state is defined by V = 0, dn/dt = dm/dt = dh/dt = 0
# calculate the resting conductances
n_0 = alpha_n(0) / (alpha_n(0) + beta_n(0))
m_0 = alpha_m(0) / (alpha_m(0) + beta_m(0))
h_0 = alpha_h(0) / (alpha_h(0) + beta_h(0))
V_0 = -90
print(" Initial depolarization V(0) =", V_0, "mV")
t = 0
Vnmht = [ V_0, n_0, m_0, h_0, t ]
t_max = 6
dt = 0.01
dt_min, dt_max = [dt, dt]
print(" Integrating using RK4 with adaptive step size dt =", dt)
print("  t          V(t)          n(t)          m(t)          h(t)")
print(" ----------------------------------------------------------------")
skip_steps = 10
step = 0
file = open("hodgkin-huxley.data", "w")
while t < t_max + dt:
    if step % skip_steps == 0:
        print(" ", t, Vnmht[0], Vnmht[1], Vnmht[2], Vnmht[3])
    data = str(t)
    for i in range(4):
        data += '   ' + str(Vnmht[i])
    file.write(data + '\n')
    dt = RK4_adaptive_step(Vnmht, dt, HH_equations)
    dt_min = min(dt_min, dt)
    dt_max = max(dt_max, dt)
    t = Vnmht[4]
    step += 1
file.close()
print(" min, max adaptive dt =", dt_min, dt_max)
print(" Data in file hodgkin-huxley.data")
