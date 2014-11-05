import math
import cpt


class ParticleInBox :
    E = 0.0             # current value of E for root finding
    accuracy = 1.0e-6   # accuracy for adaptive RK4 integration
    def __init__( self, E, accuracy ) :
        self.E = E
        self.accuracy = accuracy


class ParticleInBoxEnergy( ParticleInBox ) :
    def __init__( self, E, accuracy ):
        ParticleInBox.__init__( self, E, accuracy )
    def __call__(self, xy):   # time independent Schroedinger equation
        x = xy[0] ; psi = xy[1] ; psi_prime = xy[2]
        dx_dx = 1.0
        dpsi_dx = psi_prime
        dpsi_prime_dx = ( self.V(x) - self.E ) * psi
        return [ dx_dx, dpsi_dx, dpsi_prime_dx ]
    def V(self, x):           # potential function
        if x >= 0.0 and x <= 1.0:
            return 0.0
        else:
            return 1.0e30

class ParticleInBoxWavefunction( ParticleInBox ) :
    def __init__( self, E, accuracy ):
        ParticleInBox.__init__( self, E, accuracy )
        self.particleE = ParticleInBoxEnergy( E, accuracy)
    def __call__(self, energy):      # search function for root finding
        self.E = energy
        self.particleE.E = energy
        x = 0.0         # left wall
        psi = 0.0       # psi(0)
        psi_prime = 1.0 # arbitrary value of slope
        xy = [ x, psi, psi_prime ]
        Delta_x = 1.0   # width of well
        dx = 0.001      # suggested step size of Runge-Kutta
        cpt.RK4_integrate(xy, dx, self.particleE, Delta_x, self.accuracy)
        return xy[1]    # searching for psi(x=1) = 0



print " Bound state energies of infinitely deep potential well"
print " ------------------------------------------------------"
E0, E1, accuracy = input(" Enter guess E_0, E_1, and desired accuracy: ")
particlePsi = ParticleInBoxWavefunction( E0, accuracy )
E = cpt.root_bisection(particlePsi, E0, E1, particlePsi.accuracy, 1000, True)
print " Eigenvalue =", E
