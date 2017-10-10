# Variational Monte Carlo for the harmonic oscillator

import math
import random

def zero_accumulators():
    global psi_sqd, e_sum, e_sqd_sum
    e_sum = e_sqd_sum = 0.0
    for i in range(n_psi_sqd):
        psi_sqd[i] = 0.0

def initialize():
    global x, n_psi_sqd, psi_sqd
    x = []
    for i in range(N):
        x.append(random.random() - 0.5)
    n_psi_sqd = int((x_max - x_min) / dx)
    psi_sqd = [ 0.0 ] * n_psi_sqd
    zero_accumulators()

# the trial function is exp(-alpha x^2)

def p(x_trial, x_i):
    # compute ratio of rho(x_trial) / rho(x_i)
    return math.exp( - 2 * alpha * (x_trial**2 - x_i**2))

def e_local(x_i):
    # compute the local energy
    return alpha + x_i**2 * (0.5 - alpha**2)

def Metropolis_step():
    global x, n_accept, e_sum, e_sqd_sum, psi_sqd
    # choose a walker at random
    i = random.randrange(N)
    # make a trial move
    x_trial = x[i] + delta * random.normalvariate(0.0, 1.0)
    # Metropolis test
    if p(x_trial, x[i]) > random.random():
        x[i] = x_trial
        n_accept += 1
    # accumulate energy and wave function
    e = e_local(x[i])
    e_sum += e
    e_sqd_sum += e**2
    bin = int((x[i] - x_min) / dx)
    if bin >= 0 and bin < n_psi_sqd:
        psi_sqd[bin] += 1.0

def one_Monte_Carlo_step():
    # perform N Metropolis steps
    for i in range(N):
        Metropolis_step()

print(" Variational Monte Carlo for the Harmonic Oscillator")
print(" ---------------------------------------------------")
N = int(input(" Enter number of walkers: "))
alpha = float(input(" Enter parameter alpha: "))
MC_steps = int(input(" Enter number of Monte Carlo steps: "))

x_max = 10.0        # maximum x for histogram
x_min = -x_max
dx = 0.1            # histogram bin width
delta = 1.0         # Metropolis step size

initialize()

# perform 20% of MC_steps as thermalization steps
# and adjust step size so acceptance ratio ~50%
therm_steps = int(0.2 * MC_steps)
adjust_interval = int(0.1 * therm_steps)
n_accept = 0
print(" Performing", therm_steps, "thermalization steps ...")
for i in range(therm_steps):
    one_Monte_Carlo_step()
    if (i + 1) % adjust_interval == 0:
        delta *= n_accept / (0.5 * N * adjust_interval)
        n_accept = 0
print(" Adjusted Gaussian step size =", delta)

# production steps
zero_accumulators()
n_accept = 0
print(" Performing", MC_steps, "production steps ...")
for i in range(MC_steps):
    one_Monte_Carlo_step()

# compute and print energy
e_ave = e_sum / float(N) / MC_steps
e_var = e_sqd_sum / float(N) / MC_steps - e_ave**2
error = math.sqrt(e_var) / math.sqrt(float(N) * MC_steps)
print(" <Energy> =", e_ave, "+/-", error, "Variance =", e_var)

# write wave function squared in file
file = open("psi_sqd.data", "w")
psi_norm = 0.0
for i in range(n_psi_sqd):
    psi_norm += psi_sqd[i] * dx
for i in range(n_psi_sqd):
    x_i = x_min + i * dx
    file.write(" " + str(x_i) + "\t" + str(psi_sqd[i] / psi_norm) + "\n")
file.close()
print(" Probability density written in file psi_sqd.data")
