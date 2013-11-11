import math
import cmath
import matplotlib.pyplot as plt

class KronigPenney : 
	def __init__(self) : 
		self.a = 1.0                             # size of unit cell - lattice spacing
		self.V_0 = -5.0                          # height of potential barrier
		self.Delta = 0.2                         # width of potential barrier

	def solve_for_E(self, E, k):              # to solve 2x2 eigenvalue problem
	    # E is the desired enegy (input)
	    # k is a list of the two solutions
	    q = math.sqrt(2 * E)
	    kappa = math.sqrt(2 * (E - self.V_0))
	    i = 1.0j
	    T11 = ( cmath.exp(i * q * (self.a - self.Delta)) / (4 * q * kappa) *
		    ( cmath.exp(i * kappa * self.Delta)  * (q + kappa)**2 -
		      cmath.exp(-i * kappa * self.Delta) * (q - kappa)**2   ) )
	    T22 = T11.conjugate()
	    T12 = ( -i * cmath.exp(i * q * (self.a - self.Delta)) / (2 * q * kappa) *
		    (q**2 - kappa**2) * math.sin(kappa * self.Delta)  )
	    T21 = T12.conjugate()

	    # solve quadratic determinantal equation
	    b = - (T11 + T22)
	    c = (T11 * T22 - T12 * T21)
	    k[0] = (- b + cmath.sqrt(b**2 - 4*c)) / 2.0
	    k[1] = (- b - cmath.sqrt(b**2 - 4*c)) / 2.0
	    for j in range(2):
		k[j] = cmath.log(k[j]) / (i * self.a)

	def compute_bands(self, dE, steps, band_file_name):
	    # dE = step size in E for search
	    # steps = number of steps
	    E_values = []
	    rq_values = []
	    file = open(band_file_name, "w")
	    E = dE
	    for step in range(steps):
		q = [ 0.0 + 0.0j, 0.0 + 0.0j ]
		self.solve_for_E(E, q)
		for j in range(2):
		    rq = q[j].real
		    if rq > 0.0 and rq < math.pi / self.a:
			rq_values.append( rq )
			rq_values.append( -rq )
			E_values.append( E )
			E_values.append( E )
		E += dE
	    return rq_values, E_values



kp = KronigPenney()
print " Kronig-Penney Model"
dE = 0.01
steps = 3000

print " V_0 =", kp.V_0, " Delta =", kp.Delta
rq_values1, E_values1 = kp.compute_bands(dE, steps, "band.data")
kp.V_0 = 0.0
print " V_0 =", kp.V_0, " Delta =", kp.Delta
rq_values2, E_values2 = kp.compute_bands(dE, steps, "band0.data")

ax1 = plt.subplot(2,1,1)
plt.scatter( rq_values1, E_values1 )

ax2 = plt.subplot(2,1,2)
plt.scatter( rq_values2, E_values2 )


plt.show()
