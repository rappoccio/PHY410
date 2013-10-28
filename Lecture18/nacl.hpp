#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

#include "linalg.hpp"
#include "nonlin.hpp"
using namespace cpt;

class Ion {                     // defines an ion object
public:                         // all properties are
    string name;
    int q;                      // charge = +1 for sodium, -1 for chloride
    double m;                   // mass in atomic units
    double a;                   // radius in nm
    Matrix<double,1> r;         // position vector (x,y,z)
    Ion() : r(3), q(0),         // constructor initializes r = (0,0,0), q = 0,
        m(1.0), a(0.0) { }      // m = 1.0, and a = 0
};

class Na : public Ion {         // Na inherits members of Ion
    void set_constants() { q = +1;  m = 22.990;  a = 0.116;  name = "Na+"; }
public:
    Na() {                      // no argument constructor
        Ion();                  // invokes Ion constructor
        set_constants();
    }
    Na(Matrix<double,1> r) {    // argument r hides member r
        this->r = r;            // C++ keyword "this" is a pointer to object
        set_constants();
    }
};

class Cl : public Ion {
public:
    void set_constants() { q = -1;  m = 35.453;  a = 0.169;  name = "Cl-"; }
public:
    Cl() {                      // no argument constructor
        Ion();                  // invokes Ion constructor
        set_constants();
    }
    Cl(Matrix<double,1> r) {    // argument r hides member r
        this->r = r;            // C++ keyword "this" is a pointer to object
        set_constants();
    }
};

class IonPotential {

public : 
  IonPotential(double ik=1.44, double ialpha=1.09e3, double irho=0.0321,
	       double ib=1.0, double ic=0.01 )  :
    k(ik), alpha(ialpha), rho(irho), b(ib), c(ic)
  {
  }
  ~IonPotential() {}


  double potential( Ion const & ion1, Ion const & ion2 ) const {

    Matrix<double,1> r12 = ion1.r - ion2.r;
    double r = sqrt(dot(r12, r12));
    double V = ion1.q * ion2.q * k / r + b * pow(c/r, 12.0);
    if (ion1.q * ion2.q == -1)
        V += alpha * exp(-r/rho);
    return V;
    
  }


  Matrix<double,1> potential_derivative(Ion const & ion1, Ion const & ion2) const
  {
    Matrix<double,1> r12 = ion1.r - ion2.r;
    double r = sqrt(dot(r12, r12));
    Matrix<double,1> rhat = r12 / r;
    double dV_mag = - ion1.q * ion2.q * k / (r * r)
      - 12.0 * b / r * pow(c/r, 12.0);
    if (ion1.q * ion2.q == -1)
      dV_mag += - alpha / rho * exp(-r/rho);
    return rhat * dV_mag;
  }

protected : 
  double k ;        // Coulomb energy constant (eV nm)
  double alpha ;    // exponential repulsion prefactor (eV)
  double rho ;      // range of exponential repulsion (nm)
  double b ;        // 1/r^12 repulsion prefactor (eV)
  double c ;        // 1/r^12 repulsion range (nm)


};



class Cluster {
public:

    Cluster(){
    }

    void add(Ion i)
    {
        ions.push_back(i);
    }

    inline Ion const & ion( unsigned int i = 0) const { return ions[i]; }

    unsigned int size() const { return ions.size();} 

    double potential_energy()
    {
        double e = 0;
        for (int i = 0; i < ions.size() - 1; i++) {
            for (int j = i + 1; j < ions.size(); j++) {
                e += pot.potential(ions[i], ions[j]);
            }
        }
        return e;
    }

    int max_variable_xyz(int i)
    {
        if (i == 0) return 0;   // first ion is fixed
        if (i == 1) return 1;   // vary x component of second
        if (i == 2) return 2;   // vary x,y components of third
        return 3;               // vary x,y,z of others
    }

    int number_of_variables()
    {
        int n = 0;
        for (int i = 0; i < ions.size(); i++)
            n += max_variable_xyz(i);
        return n;
    }

    Matrix<double,1> get_variables()
    {
        int n = number_of_variables();
        Matrix<double,1> p(n);
        int i_p = 0;
        for (int i_ion = 0; i_ion < ions.size(); i_ion++)
            for (int i = 0; i < max_variable_xyz(i_ion); i++)
                p[i_p++] = ions[i_ion].r[i];
        return p;
    }

    void set_variables(Matrix<double,1> const & p)
    {
        int i_p = 0;
        for (int i_ion = 0; i_ion < ions.size(); i_ion++)
            for (int i = 0; i < max_variable_xyz(i_ion); i++)
                ions[i_ion].r[i] = p[i_p++];
    }

  double minimize( double accuracy, int & iterations ) {
    double pe = 0.0;
    Matrix<double,1> p = get_variables();
    minimize_BFGS(p, accuracy, iterations, pe, *this, *this);
    return pe;
  }

  double operator() ( Matrix<double,1> & p)  {
    set_variables(p);
    return potential_energy();
  }


  void operator()(Matrix<double,1> & p, Matrix<double,1> & dp)  {
    set_variables(p);
    for (int i_p = 0; i_p < p.size(); i_p++) {
      int i_ion;      // index of this ion
      int i_xyz;      // index in (x,y,z) of this ion
      switch (i_p) {
      case 0: i_ion = 1; i_xyz = 0; break;
      case 1: i_ion = 2; i_xyz = 0; break;
      case 2: i_ion = 2; i_xyz = 1; break;
      default: i_ion = 2 + i_p / 3; i_xyz = i_p % 3; break;
      }
      Matrix<double,1> sum_dV(3);
      for (int i = 0; i < ions.size(); i++)
	if (i != i_ion)
	  sum_dV = sum_dV + pot.potential_derivative(ions[i_ion], ions[i]);
      dp[i_p] = sum_dV[i_xyz];
    }
  }

protected : 
  vector<Ion> ions;
  IonPotential pot;

};

std::ostream& operator<<(std::ostream& os, const Cluster& c)
{
    for (int i = 0; i < c.size(); i++) {
        for (int j = 0; j < 3; j++)
	  os << '\t' << c.ion(i).r[j];
        os << '\n';
    }
    return os;
}
