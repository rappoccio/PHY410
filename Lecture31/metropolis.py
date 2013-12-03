import math
import random

class Metropolis :
    def __init__(self, delta=1.0) : 
        self.x = 0.0             # initial position of walker
        self.delta = delta       # step size
        self.accepts = 0         # number of steps accepted

    def step(self):
        x_trial = self.x + self.delta * random.uniform(-1, 1)
        ratio = self.P(x_trial) / self.P(self.x)
        if (ratio > random.random()):
            self.x = x_trial
            self.accepts += 1
                
    def P(self, x):
        # normalized Gaussian function
        return math.exp(-x**2 / 2.0) / math.sqrt(2 * math.pi)

    def f_over_w(self):
        # integrand divided by weight function
        return float(self.x**2)


print " Monte Carlo Quadrature using Metropolis et al. Algorithm"
print " ========================================================"
delta = float(input(" Enter step size delta: "))
M = int(input(" Enter number of trials M: "))
N = int(input(" Enter number of Metropolis steps per trial N: "))

f_sum = 0.0         # accumulator for f(x) values
f2_sum = 0.0        # [f(x)]**2 values
err_sum = 0.0       # error estimates

metropolis = Metropolis( delta=delta )

for i in range(M):
    avg = 0.0
    var = 0.0
    for j in range(N):
        metropolis.step()
        fx = metropolis.f_over_w()
        avg += fx
        var += fx**2
    avg /= float(N)
    var /= float(N)
    var = var - avg**2
    err = math.sqrt(var / N)
    f_sum += avg
    f2_sum += avg**2
    err_sum += err
ans = f_sum / float(M)
std_dev = math.sqrt(f2_sum / M - ans * ans)
std_dev /= math.sqrt(M-1.0)
err = err_sum / float(M)
err /= math.sqrt(M)

print ""
print " Exact answer =", 1.0
print "     Integral =", ans, "+-", err
print "    Std. Dev. =", std_dev
print " Accept ratio =", metropolis.accepts / float(N*M)
