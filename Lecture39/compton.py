from cmath import *
import vegas

class Compton:
    def __init__(self, omega):
        self.omega = omega             # incident photon angular frequency
        self.ep_10 = complex(0,0)      # x-component of incident photon polarization
        self.ep_20 = complex(0,0)      # y-component
        self.ep_1 = complex(0,0)       # in-plane scattered photon polarization
        self.ep_2 = complex(0,0)       # perpendicular component

    def __call__(self,x):
        theta = x[0]
        phi = x[1]
        omega_prime = self.omega / (1 + self.omega * (1 - cos(theta)))

        # compute polarization dot products
        ep_star_dot_ep0 = ((self.ep_1 * cos(theta) * cos(phi) - self.ep_2 * sin(phi)) * self.ep_10).conjugate() + ((self.ep_1 * cos(theta) * sin(phi) + self.ep_2 * cos(phi)) * self.ep_20).conjugate()

        ep_dot_ep0 = (self.ep_1 * cos(theta) * cos(phi) - self.ep_2 * sin(phi)) * self.ep_10 + (self.ep_1 * cos(theta) * sin(phi) + self.ep_2 * cos(phi)) * self.ep_20

        return abs( pow(omega_prime / self.omega, 2.0) * ( abs(ep_star_dot_ep0) + (self.omega - omega_prime)**2 / (4 * self.omega * omega_prime) * (1 + abs(ep_star_dot_ep0) - abs(ep_dot_ep0))) )


def main():

    print " Relativistic Compton Scattering Cross Sections using VEGAS\n"


    regn = [ [0.0,  pi],       # 0 < theta < pi
             [0.0, 2 * pi]]    # 0 < phi < 2 pi
    compton = Compton(omega=0.0024)    
    # define and normalize polarization vectors
    compton.ep_10 = 1.0
    compton.ep_20 = 0.0
    compton.ep_1 = 1.0
    compton.ep_2 = 0.0
    n0 = sqrt( abs(compton.ep_10) + abs(compton.ep_20))
    n = sqrt( abs(compton.ep_1) + abs(compton.ep_2))
    compton.ep_10 /= n0
    compton.ep_20 /= n0
    compton.ep_1 /= n
    compton.ep_2 /= n

    init = 0                       # initialize grid
    ncall = 50000                  # number of function evaluations
    itmx = 10                      # number of iterations
    nprn = 0                       # prsome info after each iteration

    # assign integration volume to integrator
    integ = vegas.Integrator(regn)

    # adapt to the integrand; discard results
    integ(compton, nitn=5, neval=1000)

    # do the final integral
    result = integ(compton, nitn=10, neval=10000)

    print result

if __name__ == "__main__":
    main()
