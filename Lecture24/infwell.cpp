#include "cptstd.hpp"
#include "basalg.hpp"
#include "diffeq.hpp"
using namespace cpt;

class ParticleInBox {

public : 
  enum extended_space { _x, _psi, _dpsi_dx, dimension };


  ParticleInBox( double acc ) : accuracy(acc) { E = 0.0;}

protected : 
  double E;                       // current value of E for root finding
  double accuracy;                // for integration and root finding
};


class ParticleInBoxEnergy : public ParticleInBox {
public : 

  ParticleInBoxEnergy( double iE, double acc ) : ParticleInBox( acc), E(iE) {}

  double V(double x) const {            // potential function
    if (x >= 0 && x <= 1)
      return 0;
    else return 1e30;           // large value
  }

  // flow vector for Runge-Kutta integration
  Matrix<double,1> operator()(Matrix<double,1>& p)  const       // point in extended space
  {
    double x = p[_x], psi = p[_psi], dpsi_dx = p[_dpsi_dx];
    Matrix<double,1> flow(dimension);
    flow[_x] = 1;
    flow[_psi] = dpsi_dx;
    flow[_dpsi_dx] = (V(x) - E) * psi;
    return flow;
  }
protected : 
  double E;

};

class ParticleInBoxWavefunction : public ParticleInBox {
public : 
  ParticleInBoxWavefunction( double acc ) : ParticleInBox( acc) {}

  double operator()(double iE)   const           // search function for root finding
  {

    ParticleInBoxEnergy energy( iE, accuracy );

    Matrix<double,1> p(dimension);
    p[_x] = 0;                  // left wall
    p[_psi] = 0;                // psi(0)
    p[_dpsi_dx] = 1;            // arbitrary value of slope
    double dx = 0.001;          // suggested step size for Runge-Kutta
    double Delta_x = 1.0;       // integration interval
    RK4_integrate(p, dx, energy, Delta_x, accuracy);

    return p[_psi];             // looking for psi(1) = 0
  }

  double solve(double E0, double E1) {
    // use bisection search to find the energy eigenvalue
    BisectionSearchT<ParticleInBoxWavefunction> bs;                 // Bisection search object
    bs.set_first_root_estimate(E0);
    bs.set_second_root_estimate(E1);
    bs.set_accuracy(accuracy);
    bs.print_steps();                   // root finder will print each step
    bs.find_root(*this);
  }

};

int main() {

    cout << " Bound state energies of infinitely deep potential well\n"
         << " ------------------------------------------------------\n"
         << " Enter bracketing guesses E_0, E_1, and desired accuracy: ";
    double E0, E1;
    double accuracy;
    cin >> E0 >> E1 >> accuracy;

    ParticleInBoxWavefunction infwell(accuracy);
    double E = infwell.solve(E0, E1);
    cout << " Eigenvalue E = " << E << endl;

    // integrate to get wavefunction
    ofstream data_file("infwell.data");
    int N = 100;
    double dx = 1.0 / N;
    Matrix<double,1> p(ParticleInBox::dimension);
    p[ParticleInBox::_x] = 0, p[ParticleInBox::_psi] = 0, p[ParticleInBox::_dpsi_dx] = 1;
    ParticleInBoxEnergy energy(E,accuracy);
    data_file << p[ParticleInBox::_x] << '\t' << p[ParticleInBox::_psi] << '\t' << p[ParticleInBox::_dpsi_dx] << '\n';
    for (int i = 0; i < N; i++) {
        RK4_step(p, dx, energy);          // take one Runge-Kutta step
        data_file << p[ParticleInBox::_x] << '\t' << p[ParticleInBox::_psi] << '\t' << p[ParticleInBox::_dpsi_dx] << '\n';
    }
    data_file.close();
    cout << "\n " << N << " step solution in file infwell.data"
         << endl;
}
