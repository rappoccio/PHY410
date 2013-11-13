#include "cptstd.hpp"
#include "matrix.hpp"
using namespace cpt;

class Poisson {
public : 
  Poisson( int iL = 50 )  : 
    L(iL)
  {
    int N = L + 2;
    V = rho = V_new = Matrix<double,2>(N, N);


    h = 1 / double(L + 1);      // assume physical size in x and y = 1
    double q = 10;              // point charge
    int i = N / 2;              // center of lattice
    rho[i][i] = q / (h * h);    // charge density
    steps = 0;
  }

  void Jacobi() {
    // Jacobi algorithm for a single iterative step
    for (int i = 1; i <= L; i++)
      for (int j = 1; j <= L; j++)
	V_new[i][j] = 0.25 * (V[i - 1][j] + V[i + 1][j] +
			      V[i][j - 1] + V[i][j + 1] +
			      h * h * rho[i][j]);
  }
  double relative_error()
  {
    double error = 0;           // average relative error per lattice point
    int n = 0;                  // number of non-zero differences

    for (int i = 1; i <= L; i++)
      for (int j = 1; j <= L; j++) {
	if (V_new[i][j] != 0)
	  if (V_new[i][j] != V[i][j]) {
	    error += abs(1 - V[i][j] / V_new[i][j]);
	    ++n;
	  }
      }
    if (n != 0)
      error /= n;
    return error;
  }

  void Gauss_Seidel()
  {
    // copy V to V_new
    V_new = V;

    // Gauss-Seidel update in place
    for (int i = 1; i <= L; i++)
      for (int j = 1; j <= L; j++)
	V_new[i][j] = 0.25 * (V_new[i - 1][j] + V_new[i + 1][j] +
			      V_new[i][j - 1] + V_new[i][j + 1] +
			      h * h * rho[i][j]);
  }

  void successive_over_relaxation()   // using red-black checkerboard updating
  {
    // update even sites first
    for (int i = 1; i <= L; i++)
      for (int j = 1; j <= L; j++)
	if ((i + j) % 2 == 0)
	  V_new[i][j] = (1 - omega) * V[i][j] + omega / 4 *
	    (V[i - 1][j] + V[i + 1][j] +
	     V[i][j - 1] + V[i][j + 1] +
	     h * h * rho[i][j]);

    // update odd sites using updated even sites
    for (int i = 1; i <= L; i++)
      for (int j = 1; j <= L; j++)
	if ((i + j) % 2 != 0)
	  V_new[i][j] = (1 - omega) * V[i][j] + omega / 4 *
	    (V_new[i - 1][j] + V_new[i + 1][j] +
	     V_new[i][j - 1] + V_new[i][j + 1] +
	     h * h * rho[i][j]);
  }

  template< typename T> 
  void iterate( T const & method,double accuracy)
  {
    clock_t t0 = clock();

    while (true) {
      (this->*method)();
      ++steps;
      double error = relative_error();
      if (error < accuracy)
	break;
      swap(V, V_new);         // use <algorithm> std::swap
    }
    cout << " Number of steps = " << steps << endl;

    clock_t t1 = clock();
    cout << " CPU time = " << double(t1 - t0) / CLOCKS_PER_SEC
         << " sec" << endl;
  }

  void set_omega ( double i ) { omega = i; }

  double get_h() const { return h; }
  int    get_L() const { return L; }
  double get_V( int i, int j) const { return V[i][j]; }

protected : 
  int L ;                         // number of interior points in x and y
  Matrix<double,2> V,             // potential to be found
    rho,                          // given charge density
    V_new;                        // new potential after each step

  double h;                       // lattice spacing
  int steps;                      // number of iteration steps
  double accuracy;                // desired accuracy in solution
  double omega;                   // overrelaxation parameter


};

int main() {

  int L;
    cout << " Iterative solution of Poisson's equation\n"
         << " ----------------------------------------\n";
    cout << " Enter number of interior points in x or y: ";
    cin >> L;

    Poisson p(L);

    double accuracy; 
    cout << " Enter desired accuracy in solution: ";
    cin >> accuracy;
    cout << " Enter 0 for Jacobi, 1 for Gauss Seidel, 2 for SOR: ";
    int choice;
    cin >> choice;

    void (Poisson::* ptfptr) (void);
    switch (choice) {
    case 0:
      ptfptr = &Poisson::Jacobi;
      p.iterate(ptfptr, accuracy);

      break;
    case 1:
      ptfptr = &Poisson::Gauss_Seidel;
      p.iterate(ptfptr, accuracy);

      break;
    case 2: default : 
      ptfptr = &Poisson::successive_over_relaxation;
      double omega = 2 / (1 + 4 * atan(1.0) / double(L));
      p.set_omega( omega ) ;
      p.iterate(ptfptr, accuracy);

      break;
    }


    // write potential to file
    cout << " Potential in file poisson.data" << endl;
    ofstream date_file("poisson.data");
    for (int i = 0; i < p.get_L() + 2; i++) {
      double x = i * p.get_h();
      for (int j = 0; j < p.get_L() + 2; j++) {
	double y = j * p.get_h();
	char buff[1000];
	sprintf(buff, "%12.6f %12.6f %12.6f", x, y, p.get_V(i,j) );
	date_file << buff << endl;
      }
    }
    date_file.close();
}
