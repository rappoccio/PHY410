from cpt import *
from math import *
import scipy
import numpy as np

class Ion :                     # defines an ion object
                                # all properties are
    """
    Defines the base class for an ion object.
    The parameters are :
      - charge ( in e )
      - mass (in AMU)
      - radius (in nm)
      - position (x,y,z in nm)
    This is a class that will be used by Na and Cl
    below to store the common data members.
    The three-vectors are stored as numpy array objects
    to facilitate vector-like operations. 
    """
    name=""
    q=1.0                       # charge = +1 for sodium, -1 for chloride
    m=0.0                       # mass in atomic units
    a=0.0                       # radius in nm
    r=None                      # position vector (x,y,z)
    def __init__( self ) :
        r = np.array(ndmin=3)
        m = 1.0
        a = 0.0
        q = 0.0

class Na (Ion) :                # Na inherits members of Ion
    """
    Na sets the parameters of an Ion with
    parameters appropriate to sodium    
    """
    def set_constants(self) :
        self.q = +1.0
        self.m = 22.990
        self.a = 0.116
        self.name = "Na+"

    def __init__( self, r = None ) : 
        if r != None :
            self.r = np.array(r)
        self.set_constants()


class Cl ( Ion ) :
    """
    Cl sets the parameters of an Ion with
    parameters appropriate to chlorine
    """
    def set_constants(self) :
        self.q = -1.0
        self.m = 35.453
        self.a = 0.169
        self.name = "Cl-"

    def __init__(self, r = None ) :
        if r != None :
            self.r = np.array(r)
        self.set_constants()


class IonPotential :
    """
    // The ion potential class represents the
    // potential energy and gradient for a set of ions that
    // obey the "Alkali-Halide" model we discussed in class.
    // The parameters are :
    //   double k ;        // Coulomb energy constant (eV nm)
    //   double alpha ;    // exponential repulsion prefactor (eV)
    //   double rho ;      // range of exponential repulsion (nm)
    //   double b ;        // 1/r^12 repulsion prefactor (eV)
    //   double c ;        // 1/r^12 repulsion range (nm)
    // The functional form is : 
    // V(rij) = -ke^2 / rij + alpha*e^(-rij/rho) + b(c/rij)^12
    """
    k = 1.44
    alpha = 1.09e3
    rho=0.0321
    b = 1.0
    c = 0.01

    def __init__( self, k=1.44, alpha=1.09e3, rho=0.0321,
                  b=1.0, c=0.01 )  :
        self.k = k
        self.alpha = alpha
        self.rho = rho
        self.b = b
        self.c = c


    def potential( self, ion1, ion2 ) :
        r12 = ion1.r - ion2.r
        r = sqrt( np.dot(r12, r12))
        V = ion1.q * ion2.q * self.k / r + self.b * pow(self.c/r, 12.0)
        if ion1.q * ion2.q == -1 :
            V += self.alpha * exp(-r/self.rho)
        return V


    def potential_derivative(self, ion1, ion2) :
        r12 = ion1.r - ion2.r
        r = sqrt(np.dot(r12, r12))
        rhat = np.array( r12 ) / r
        dV_mag = - ion1.q * ion2.q * self.k / (r * r) - 12.0 * self.b / r * pow(self.c/r, 12.0)
        if ion1.q * ion2.q == -1:
            dV_mag =  dV_mag - self.alpha / self.rho * exp(-r/self.rho)
        return rhat * dV_mag




class Cluster :
    """
    // The Cluster class represents a group of 
    // sodium and chlorine atoms in some particular
    // arrangement. 
    // This utilizes the "IonPotential" class from above
    // to compute the inter-ion potentials. 
    // The "minimize" function will use the BFGS algorithm
    // to minimize the potential energy, modifying the
    // inter-ion distances until it is reached. 
    // The data members are : 
    // ions;                 // Vector of ions
    // IonPotential pot;     // potential functor
    //
    //
    // This also defines two helper classes to use as callable objects
    // in the minimization routine, callPotential and callPotentialDerivative.
    // These are passed a copy of the self object, and they have __call__
    // methods implemented so that the fmin_bfgs method can call them.
    //
    // This uses the scipy.optimize implementation of the bfgs algorithm.
    // 
    """
    ions = []
    pot = IonPotential()    

    def __init__( self ) :
        self.ions = []
        self.pot = IonPotential()
    
    def add(self, i ) :
        self.ions.append(i)
        
    def ion(self, i=0):
        return self.ions[i]

    def size(self) :
        return len(self.ions) 

    def potential_energy(self) :
        e = 0.0
        for i in range( 0, self.size()-1 ) :
            for j in range( i+1, self.size() ) :
                e += self.pot.potential( self.ions[i], self.ions[j] )
        return e


    def max_variable_xyz(self,i) :
        return min( i, 3)

    def number_of_variables(self) :
        n = 0
        for i in xrange(self.size()) :
            n += self.max_variable_xyz(i)
        return n

    def get_variables(self) :
        n = self.number_of_variables()
        p = [0.0] * n 

        i_p = 0
        for i_ion in xrange(self.size()) :
            for i in xrange(self.max_variable_xyz(i_ion)) :
                p[i_p] = self.ions[i_ion].r[i]
                i_p += 1
        return p

    def set_variables(self, p) :
        i_p = 0
        for i_ion in xrange(self.size()) :
            for i in xrange(self.max_variable_xyz(i_ion)) :
                self.ions[i_ion].r[i] = p[i_p]
                i_p += 1



    class callPotential :
        cluster = None
        def __init__ (self, cluster) :
            self.cluster=cluster
        def __call__ ( self, p ) :
            self.cluster.set_variables(p)
            return self.cluster.potential_energy()

    class callPotentialDerivative :
        cluster = None
        def __init__ (self, cluster ) :
            self.cluster=cluster
        def __call__ (self,p ) : 
            self.cluster.set_variables(p)
            dp = np.array( [0.0] * len(p) )
            for i_p in xrange(len(p) ) : 
                i_ion = 0      # index of this ion
                i_xyz = 0      # index in (x,y,z) of this ion
                if i_p == 0 :
                    i_ion = 1
                    i_xyz = 0
                elif i_p == 1 :
                    i_ion = 2
                    i_xyz = 0
                elif i_p == 2 :
                    i_ion = 2
                    i_xyz = 1
                else  :
                    i_ion = 2 + i_p / 3
                    i_xyz = i_p % 3

                sum_dV = np.array([0.0] * 3)
                for i in range(len(self.cluster.ions)) :
                    if i != i_ion :
                        sum_dV = sum_dV + self.cluster.pot.potential_derivative(self.cluster.ions[i_ion], self.cluster.ions[i]) 
                dp[i_p] = sum_dV[i_xyz]
            return dp



    def minimize( self, accuracy ) :
        pe = 0.0
        p = self.get_variables()
        cp = self.callPotential( self )
        cpder = self.callPotentialDerivative( self )
        res = scipy.optimize.fmin_bfgs(x0=p, gtol=accuracy, f=cp, fprime=cpder, full_output=True )
        return res


    def __str__(self) :
        s = ''
        for i in xrange(self.size()) :
            for j in xrange(3) :
                s +=  '{0:8.3f}'.format( self.ion(i).r[j] )
            s += '\n'
        return s

    def convert(self) :
        x = []
        y = []
        z = []
        for i in xrange(self.size()) :
            x.append( self.ion(i).r[0] )
            y.append( self.ion(i).r[1] )
            z.append( self.ion(i).r[2] )
        return x,y,z
