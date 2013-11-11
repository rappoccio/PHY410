import math
import cpt
import matplotlib.pyplot as plt

class Schroedinger : 

    def __init__(self) :
	self.hbar = 1.0                  # Planck's constant / 2pi
	self.m = 1.0                     # particle mass
	self.omega = 1.0                 # oscillator frequency
	self.E = 0.0                     # current energy in search
	self.N = 500                     # number of lattice points = N+1
	self.x_left = -5.0               # left boundary
	self.x_right = 5.0               # right boundary
	self.h = (self.x_right - self.x_left) / self.N  # grid spacing

	self.phi_left = [0.0]*(self.N+1)      # wave function integrating from left
	self.phi_right = [0.0]*(self.N+1)     # wave function integrating from right
	self.phi = [0.0]*(self.N+1)           # whole wave function

	self.sign = 1                    # current sign used to make F(E) continuous
	self.nodes = 0                   # current number of nodes in wavefunction


    def V(self, x):                   # harmonic oscillator potential
        return 0.5 * self.m * self.omega**2 * x**2


    def q(self, x):                   # Sturm-Liouville q function
        return 2 * self.m / self.hbar**2 * (self.E - self.V(x))


    def F(self, energy):              # eigenvalue at F(E) = 0

        # set energy needed by the q(x) function
        self.E = energy

        # find the right turning point
        i_match = self.N
        x = self.x_right             # start at right boundary
        while self.V(x) > self.E:         # in forbidden region
            i_match -= 1
            x -= self.h
            if i_match < 0:
                raise Exception("can't find right turning point")

        # integrate self.phi_left using Numerov algorithm
        self.phi_left[0] = 0.0
        self.phi_left[1] = 1.0e-10
        c = self.h**2 / 12.0         # constant in Numerov formula
        for i in range(1, i_match+1):
            x = self.x_left + i * self.h
            self.phi_left[i+1]  = 2 * (1 - 5 * c * self.q(x)) * self.phi_left[i]
            self.phi_left[i+1] -= (1 + c * self.q(x - self.h)) * self.phi_left[i-1]
            self.phi_left[i+1] /= 1 + c * self.q(x + self.h)

        # integrate self.phi_right
        self.phi[self.N] = self.phi_right[self.N] = 0.0
        self.phi[self.N-1] = self.phi_right[self.N-1] = 1.0e-10
        for i in range(self.N - 1, i_match - 1, -1):
            x = self.x_right - i * self.h
            self.phi_right[i-1]  = 2 * (1 - 5 * c * self.q(x)) * self.phi_right[i]
            self.phi_right[i-1] -= (1 + c * self.q(x + self.h)) * self.phi_right[i+1]
            self.phi_right[i-1] /= 1 + c * self.q(x - self.h)
            self.phi[i-1] = self.phi_right[i-1]

        # rescale self.phi_left
        scale = self.phi_right[i_match] / self.phi_left[i_match]
        for i in range(i_match + 2):
            self.phi_left[i] *= scale
            self.phi[i] = self.phi_left[i]

        # make F(E) continuous
        # count number of nodes in self.phi_left
        n = 0
        for i in range(1, i_match+1):
            if self.phi_left[i-1] * self.phi_left[i] < 0.0:
                n += 1

        # flip its sign when a new node develops

        if n != self.nodes:
            self.nodes = n
            self.sign = -self.sign

        return ( self.sign *
		 ( self.phi_right[i_match-1] - self.phi_right[i_match+1] - 
		   self.phi_left [i_match-1] + self.phi_left[i_match+1] ) /
		(2 * self.h * self.phi_right[i_match]) )

    def normalize(self):
        norm = 0.0
        for i in range(self.N):
            norm += self.phi[i]**2
        norm /= self.N
        norm = math.sqrt(norm)
        for i in range(self.N):
            self.phi[i] /= norm

print " Eigenvalues of the Schroedinger equation"
print " for the harmonic oscillator V(x) = 0.5 x^2"
print " ------------------------------------------"
E_max = input(" Enter maximum energy E: ")

phi_file = open("phi.data", "w")


level = 0                   # level number
E_old = 0.0                 # previous energy eigenvalue

schroedinger = Schroedinger()
schroedinger.E = 0.1 # guess and E below the ground state

# draw the potential
for i in range(schroedinger.N+1):
	x = schroedinger.x_left + i * schroedinger.h
swrite = '{0:12.6f} {1:12.6f}'.format( x, schroedinger.V(x) )
print swrite
print ''

# find the energy levels
print ""
print " Level       Energy           Simple Steps   Secant Steps"
print " -----   ------------------   ------------   ------------"

x_data = []
phi_data = []

while True:                 # loop over levels

	# estimate next E and dE
	dE = 0.5 * (schroedinger.E - E_old)
	E_old = schroedinger.E
	schroedinger.E += dE

	# use simple search to locate root with relatively low accuracy
	accuracy = 0.01
	schroedinger.E = cpt.root_simple(schroedinger.F, schroedinger.E, dE, accuracy)
	#simple_steps = cpt.root_steps()

	# use secant search with relatively high accuracy
	accuracy = 1.0e-6
	E1 = schroedinger.E + 100 * accuracy # guess second point required
	schroedinger.E = cpt.root_secant(schroedinger.F, schroedinger.E, E1, accuracy)
	#secant_steps = cpt.root_steps()

	#ans = " " + repr(level).rjust(3) + " "*5 + repr(E).ljust(20) + " "*6
	#ans += repr(simple_steps).rjust(3) + " "*11 + repr(secant_steps).rjust(3)
	#print ans
	level += 1

	accuracy = 0.001
	x = cpt.root_simple(schroedinger.q, schroedinger.x_left, schroedinger.h, accuracy)
	swrite = '{0:12.6f} {1:12.6f}'.format( x, schroedinger.E )
	print swrite
	
	x = cpt.root_simple(schroedinger.q, schroedinger.x_right, -schroedinger.h, accuracy)
	swrite = '{0:12.6f} {1:12.6f}'.format( x, schroedinger.E )
	print swrite

	print ''
	
	schroedinger.normalize()
	x_data.append( [] )
	phi_data.append( [] )
	iphi = len(x_data) - 1
	for i in range(schroedinger.N+1):
		x = schroedinger.x_left + i * schroedinger.h
		x_data[iphi].append(x)
		phi_data[iphi].append( schroedinger.phi[i] )


	if schroedinger.E >= E_max:          # we are done
		break




# print the search function
schroedinger.E = 0.1
dE = 0.01
E_data = []
F_data = []
while schroedinger.E < E_max:
	E_data.append( schroedinger.E )
	F_data.append(  schroedinger.F(schroedinger.E) )
	schroedinger.E += dE


s1 = plt.subplot(2,1,1)
for i in range(len(x_data)) :
	plt.plot( x_data[i], phi_data[i] )

s2 = plt.subplot(2,1,2)
plt.plot( E_data, F_data )
plt.ylim( [-10, 10] )

plt.show()
