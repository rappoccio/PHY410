#include "bvp.hpp"
using namespace cpt;

int main() {
  BVPStep bvpStepper;
  BVPShoot bvpShooter( bvpStepper );
    cout << " Shooting solution of u'' = - pi^2 (u+1) / 4\n"
         << " with boundary conditions u(0) = 0, u(1) = 1\n";

    cout << " using initial guess delta = " << bvpShooter.get_delta() << endl;
    double delta = bvpShooter.solve();
    cout << "\n Secant search found delta = " << delta << endl;

    bvpShooter.set_initial_values();           // to integrate the final solution
    bvpStepper.update_values( bvpShooter );    // Initialize to the solved equation
    ofstream data_file("shoot.data");

    Matrix<double,1> p = bvpStepper.get_p();

    data_file << p[BVP::_x] << '\t' << p[BVP::_u] << '\t' << p[BVP::_dudx] << '\n';

    for (int i = 0; i < bvpShooter.get_N(); i++) {
      bvpStepper.step(); 
      p = bvpStepper.get_p();
      data_file << p[BVP::_x] << '\t' << p[BVP::_u] << '\t' << p[BVP::_dudx] << '\n';
    }
    data_file.close();
    cout << "\n " << bvpStepper.get_N() << " step solution in file shoot.data"
         << endl;
}
