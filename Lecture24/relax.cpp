#include "cptstd.hpp"
#include "matrix.hpp"
using namespace cpt;
#include "bvp.hpp"

int main() {
    cout << " Relaxing solution of u'' = - pi^2 (u+1) / 4\n"
         << " with boundary conditions u(0) = 0, u(1) = 1\n";
    cout << " Enter desired accuracy: ";
    double acc;
    cin >> acc;

    cout << "\n Iteration          Change              Error"
         << "\n -----------------------------------------------------\n";

    BVPRelax relax;
    relax.initial_guess();
    relax.Jacobi_step();
    int iteration = 1;
    relax.print(iteration);
    while (relax.error() > acc) {
      relax.Jacobi_step();
      relax.print(++iteration);
    }

    ofstream data_file("relaxing.data");
    for (int i = 0; i <= relax.get_N(); i++)
      data_file << relax.get_u_lb() + i * relax.get_h() << '\t' << relax.get_u(i) << '\n';
    data_file.close();
    cout << "\n " << iteration << " step solution in file relax.data"
         << endl;
}
