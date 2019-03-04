import math
import random
import cpt
import matplotlib.pyplot as plt

print " Random walk in 1 dimension"
print " --------------------------"
n_walkers = int(input(" Enter number of walkers: "))
n_steps = int(input(" Enter number of steps: "))

# walker positions initialized at x = 0
x = [ 0.0 ] * n_walkers

steps = [ 0.0 ] * n_steps       # to save step number i (time)
x2ave = [ 0.0 ] * n_steps       # to accumulate x^2 values
sigma = [ 0.0 ] * n_steps       # to accumulate fluctuations in x^2

# loop over walkers
for walker in range(n_walkers):

    # loop over number of steps
    for step in range(n_steps):

        # take a random step
        x[walker] += random.choice( (-1, 1) )

        # accumulate data
        steps[step] = step + 1.0
        x2ave[step] += x[walker]**2
        sigma[step] += x[walker]**4

# average the squared displacements and their variances
for step in range(n_steps):
    x2ave[step] /= float(n_walkers)
    sigma[step] /= float(n_walkers)
for step in range(n_steps):
    sigma[step] = math.sqrt(sigma[step] - x2ave[step]**2)
    if sigma[step] == 0.0:      # happens for first step!
        sigma[step] = 1.0
    if n_walkers > 1:
        sigma[step] /= math.sqrt(n_walkers - 1.0)

# fit data to a straight line
a, b, sigma_a, sigma_b, chisqr = cpt.chi_square_fit(steps, x2ave, sigma)
print " Fit to straight line <x^2> = a + b n"
print " Intercept a = " + repr(a) + " +- " + repr(sigma_a)
print " Slope     b = " + repr(b) + " +- " + repr(sigma_b)
print " Chisqr/dof  = " + repr(chisqr / (n_steps - 2.0))

# store in file for plotting
file = open("rwalk.data", "w")
for step in range(n_steps):
    file.write(repr(steps[step]) + "\t" + repr(x2ave[step]) +
                    "\t" + repr(sigma[step]) + "\n")
file.close()
print " t, <x^2>, sigma in file rwalk.data"

plt.scatter( steps, x2ave )

plt.show()
