import math
import cpt


class Bloch :

    def __init__(self) : 
        self.a = 1.0                         # size of unit cell

        # potential functions

        self.V_0 = -5.0                      # height of potential barrier
        self.Delta = 0.2                     # width of potential barrier

        self.cutoff_coulomb = 0              # coulomb with cutoff at ion core radius
        self.sine_step = 1                   # sine in first half, step in second
        self.kronig_penney = 2               # as in kronig-penney.py
        self.free_particle = 3

        self.E = 0.0                         # current value of energy needed by flow


        self.potential_type = 0              # defaults to 0 = cutoff_coulomb

        self.potential_name = [
            "Cutoff Coulomb", "Sine Step", "Kronig-Penney", "Free Particle" ]

    def V(self, x):
        # Just to avoid writing "self." all the time
        a = self.a
        potential_type = self.potential_type
        V_0 = self.V_0
        Delta = self.Delta
        cutoff_coulomb = self.cutoff_coulomb
        sine_step = self.sine_step
        kronig_penney = self.kronig_penney
        free_particle = self.free_particle
        dx = x - self.a / 2.0            # displacement from center of cell
        if potential_type == cutoff_coulomb:
            dx = abs(dx)            # symmetric about center of cell
            if dx < Delta/2:        # inside ion core
                return V_0          # potential is constant
            else:
                return (V_0 *       # value at the core radius
                (Delta / 2) / dx)   # Coulomb 1/r outside ion core
        elif potential_type == sine_step:
            if dx < a/2:            # left half of cell
                return math.sin(    # inverted (dx < 0) half sine wave
                2 * math.pi *dx / a)
            else:
                d = abs(dx - a/4)   # distance from center of barrier
                if d < Delta/2:     # less than 1/2 barrier width
                    return V_0      # inside barrier
                return 0.0          # outside barrier
        elif potential_type == kronig_penney:
            dx = abs(dx)            # symmetric about center of cell
            if dx < Delta/2:        # inside step
                return V_0          # barrier height
            return 0.0              # outside barrier
        return 0.0                  # default to free particle


    # flow function for Runge-Kutta integration across unit cell
    def flow(self, y):
        x = y[0] ; phi = y[1]; d_phi_dx = y[2]
        f = [ 0.0 ] * 3
        f[0] = 1.0
        f[1] = d_phi_dx
        f[2] = 2 * (self.V(x) - self.E) * phi
        return f

    def find_bloch_momentum(self):
        # returns a tuple (in_band, k)
        # in_band is set True if in band, False if in gap
        # k = value of Bloch momentum if in band

        # first fundamental solution
        y_1 = [ 0.0 ] * 3
        x = 0.0
        psi = 1.0
        psi_prime = 0.0
        y_1[0] = x
        y_1[1] = psi
        y_1[2] = psi_prime
        dx = 0.01                   # suggested step size for Runge-Kutta
        cpt.RK4_integrate(y_1, dx, self.flow, self.a)

        # second fundamental solution
        y_2 = [ 0.0 ] * 3
        x = 0.0
        psi = 0.0
        psi_prime = 1.0 / self.a
        y_2[0] = x
        y_2[1] = psi
        y_2[2] = psi_prime
        cpt.RK4_integrate(y_2, dx, self.flow, self.a)

        # find k from dispersion relation
        f = (y_1[1] + self.a * y_2[2]) / (2 * self.a)
        if abs(f) <= 1.0:
            in_band = True
            k = math.acos(f)
        else:
            in_band = False
            k = 0.0
        return in_band, k




bloch = Bloch()
print " Band Structure in 1-D For General Potential"
print " -------------------------------------------"
potential_type = input(" Enter potential type 0, 1, 2, or 3: ")
print " Potential type =", bloch.potential_name[potential_type]
E_min, E_max, n_levels = input(" Enter E_min, E_max, and number of levels: ")


file_name = "bloch.data"
file = open(file_name, "w")
new_band = False
for n in range(n_levels):
    bloch.E = E_min + n * (E_max - E_min) / (n_levels - 1.0)
    in_band, k = bloch.find_bloch_momentum()
    if in_band and not new_band:
        file.write("\n")        # separate band values with blank line
        new_band = True
    if in_band:
        file.write(repr(k) + "\t" + repr(bloch.E) + "\n")
    else:
        new_band = False
file.close()
print " Levels in file", file_name
