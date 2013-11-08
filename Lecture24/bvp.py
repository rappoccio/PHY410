import math
import cpt


class BVP :
    x0 = 0.0        # left boundary point
    x1 = 1.0        # right boundary point
    u0 = 0.0        # boundary condition at x0
    u1 = 1.0        # boundary condition at x1
    delta = 1.0     # guess for initial slope u'(0)
    N = 100         # number of integration steps
    dx = 0.0        # step size
    xy = []         # current point in extended space (x, u, delta)
    def __init__( self, x0, x1, u0, u1, delta, N ) : 
        self.x0 = x0
        self.x1 = x1
        self.u0 = u0
        self.u1 = u1
        self.delta = delta
        self.xy = [x0, u0, delta]
        self.N = N
        self.dx = (x1-x0)/N
        

class BVPStep( BVP ) :
    def __init__(self, x0, x1, u0, u1, delta, N ) :
        BVP.__init__( self, x0, x1, u0, u1, delta, N )
    def __call__(self, xy):   # derivative vector
        x = xy[0] ; u = xy[1] ; uprime = xy[2]
        dx_dx = 1.0
        du_dx = uprime
        duprime_dx = - math.pi**2 * (u + 1) / 4
        return [ dx_dx, du_dx, duprime_dx ]
    def step( self ) :
        cpt.RK4_step(self.xy, self.dx, self)

class BVPShoot( BVP ) :
    step = None
    def __init__(self, step, x0, x1, u0, u1, delta, N ) :
        BVP.__init__( self, x0, x1, u0, u1, delta, N )
        self.step = step
    def __call__( self, x ):               # function whose root is to be found
        self.delta = x
        self.xy = [ self.x0, self.u0, self.delta ]
        for i in range(self.N):
            cpt.RK4_step(self.xy, self.dx, self.step )
        return self.xy[1] - self.u1
    def solve( self ) :
        return cpt.root_secant( self, self.delta, self.delta + self.dx )

class BVPRelax( BVP ) :
    u = None
    uOld = None
    def __init__( self, x0, x1, u0, u1, N) :
        BVP.__init__( self, x0, x1, u0, u1, 0.0, N )
        self.u = [0.0] * (N+1)
        self.uOld = [0.0] * (N+1)



    def g(self, u, x):        # RHS function
        return math.pi**2 / 4 * (u + 1)

    def uExact(self, x):      # exact solution
        return math.cos(math.pi * x / 2) + 2 * math.sin(math.pi * x / 2) - 1


    def guess(self):        # linear fit to boundary conditions
        for i in range(self.N+1):
            x = self.x0 + i * self.dx
            self.u[i] = self.u0 + (x - self.x0) / (self.x1 - self.x0) * (self.u1 - self.u0)

    def relax(self):        # Jacobi algorithm
        for i in range(self.N+1):
            self.uOld[i] = self.u[i]  # save old values
        for i in range(1, self.N):
            x = self.x0 + i * self.dx
            self.u[i] = 0.5 * ( self.uOld[i+1] + self.uOld[i-1] + self.dx**2 * self.g(self.uOld[i], x) )

    def change(self):       # compute maximum change
        maxDiff = 0.
        for i in range(1, self.N):
            diff = abs(self.u[i] - self.uOld[i])
            if diff > maxDiff:
                maxDiff = diff
        return maxDiff

    def error(self):        # compute maximum change
        maxDiff = 0.
        for i in range(1, self.N):
            x = self.x0 + i * self.dx
            diff = abs(self.u[i] - self.uExact(x))
            if diff > maxDiff:
                maxDiff = diff
        return maxDiff

