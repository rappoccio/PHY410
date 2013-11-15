#include "cptstd.hpp"
#include "matrix.hpp"
using namespace cpt;

class PoissonMG {
public : 
  PoissonMG( double iaccuracy = 0.001, int iL = 64, int in_smooth = 5) :
    accuracy(iaccuracy), L(iL), n_smooth(in_smooth)
  {
    // check that L is a power of 2 as required by multigrid
    int power_of_2 = 1;
    while (power_of_2 < L)
      power_of_2 *= 2;
    if (power_of_2 != L) {
      L = power_of_2;
      cout << " Setting L = " << L << " (must be a power of 2)" << endl;
    }

    // create (L+2)x(L+2) matrices and zero them
    psi = psi_new = rho = Matrix<double,2>(L+2, L+2);

    h = 1 / double(L + 1);      // assume physical size in x and y = 1
    double q = 10;              // point charge
    int i = L / 2;              // center of lattice
    rho[i][i] = q / (h * h);    // charge density

    steps = 0;
  }

  void Gauss_Seidel(double h, Matrix<double,2>& u, const Matrix<double,2>& f)
  {
    int L = u.dim1() - 2;

    // use checkerboard updating
    for (int color = 0; color < 2; color++)
      for (int i = 1; i <= L; i++)
	for (int j = 1; j <= L; j++)
	  if ((i + j) % 2 == color)
	    u[i][j] = 0.25 * (u[i - 1][j] + u[i + 1][j] +
			      u[i][j - 1] + u[i][j + 1] +
			      h * h * f[i][j]);
  }

  void two_grid(double h, Matrix<double,2>& u, Matrix<double,2>& f)
  {
    // solve exactly if there is only one interior point
    int L = u.dim1() - 2;
    if (L == 1) {
      u[1][1] = 0.25 * (u[0][1] + u[2][1] + u[1][0] + u[1][2] +
			h * h * f[1][1]);
      return;
    }

    // do a few pre-smoothing Gauss-Seidel steps
    for (int i = 0; i < n_smooth; i++)
      Gauss_Seidel(h, u, f);

    // find the residual
    Matrix<double,2> r(L+2, L+2);
    for (int i = 1; i <= L; i++)
      for (int j = 1; j <= L; j++)
	r[i][j] = f[i][j] +
	  ( u[i + 1][j] + u[i - 1][j] +
	    u[i][j + 1] + u[i][j - 1] - 4 * u[i][j]) / (h * h);

    // restrict residual to coarser grid
    int L2 = L / 2;
    Matrix<double,2> R(L2 + 2, L2 + 2);
    for (int I = 1; I <= L2; I++) {
      int i = 2 * I - 1;
      for (int J = 1; J <= L2; J++) {
	int j = 2 * J - 1;
	R[I][J] = 0.25 * ( r[i][j] + r[i + 1][j] + r[i][j + 1] +
			   r[i + 1][j + 1]);
      }
    }

    // initialize correction V on coarse grid to zero
    Matrix<double,2> V(L2 + 2, L2 + 2);

    // call twoGrid recursively
    double H = 2 * h;
    two_grid(H, V, R);

    // prolongate V to fine grid using simple injection
    Matrix<double,2> v(L + 2, L + 2);
    for (int I = 1; I <= L2; I++) {
      int i = 2 * I - 1;
      for (int J = 1; J <= L2; J++) {
	int j = 2 * J - 1;
	v[i][j] = v[i + 1][j] = v[i][j + 1] = v[i + 1][j + 1] = V[I][J];
      }
    }

    // correct u
    for (int i = 1; i <= L; i++)
      for (int j = 1; j <= L; j++)
	u[i][j] += v[i][j];

    // do a few post-smoothing Gauss-Seidel steps
    for (int i = 0; i < n_smooth; i++)
      Gauss_Seidel(h, u, f);
  }

  double relative_error()
  {
    double error = 0;           // average relative error per lattice point
    int n = 0;                  // number of non-zero differences

    for (int i = 1; i <= L; i++)
      for (int j = 1; j <= L; j++) {
	if (psi_new[i][j] != 0.0)
	  if (psi_new[i][j] != psi[i][j]) {
	    error += abs(1 - psi[i][j] / psi_new[i][j]);
	    ++n;
	  }
      }
    if (n != 0)
      error /= n;

    return error;
  }

  Matrix<double,2> const & get_psi() const { return psi; }
  Matrix<double,2> const & get_psi_new() const { return psi_new; }
  Matrix<double,2> const & get_rho() const { return rho; }


  Matrix<double,2>  & get_psi()  { return psi; }
  Matrix<double,2>  & get_psi_new()  { return psi_new; }
  Matrix<double,2>  & get_rho()  { return rho; }

  void set_psi    ( int i, int j, double f) { psi[i][j] = f; }
  void set_psi_new( int i, int j, double f) { psi_new[i][j] = f; }
  void set_rho    ( int i, int j, double f) { rho[i][j] = f; }

  double get_h () const { return h;}
  int    get_steps() const { return steps; }

  void inc_steps() { ++steps;}
protected : 

  double accuracy;        // desired relative accuracy in solution
  int L;                  // number of interior points in each dimension
  int n_smooth;           // number of pre and post smoothing iterations

  Matrix<double,2> psi,   // solution to be found
    psi_new,              // approximate solution after 1 iteration
    rho;                  // given source function

  double h;               // step size
  int steps;              // number of iteration steps


};

int main()
{
  int L, n_smooth;
  double accuracy;
  cout << " Multigrid solution of Poisson's equation\n"
       << " ----------------------------------------\n";
  cout << " Enter number of interior points in x or y: ";
  cin >> L;
  cout << " Enter desired accuracy in the solution: ";
  cin >> accuracy;
  cout << " Enter number of smoothing iterations: ";
  cin >> n_smooth;

  PoissonMG poissonmg(accuracy, L, n_smooth);
  clock_t t0 = clock();
  while (true) {
    for (int i = 0; i < L+2; i++)
      for (int j = 0; j < L+2; j++)
	poissonmg.set_psi_new( i, j, poissonmg.get_psi()[i][j]) ;
    poissonmg.two_grid(poissonmg.get_h(), poissonmg.get_psi(), poissonmg.get_rho());
    poissonmg.inc_steps();
    double error = poissonmg.relative_error();
    cout << " Step No. " << poissonmg.get_steps() << "\tError = " << error << endl;
    if (poissonmg.get_steps() > 1 && error < accuracy)
      break;
  }
  clock_t t1 = clock();
  cout << " CPU time = " << double(t1 - t0) / CLOCKS_PER_SEC
       << " sec" << endl;

  // write potential to file
  ofstream file("poisson_mg.data");
  for (int i = 0; i < L + 2; i++) {
    double x = i * poissonmg.get_h();
    for (int j = 0; j < L + 2; j++) {
      double y = j * poissonmg.get_h();
      file << x << '\t' << y << '\t' << poissonmg.get_psi()[i][j] << '\n';
    }
    file << '\n';
  }
  file.close();
  cout << " Potential in file poisson_mg.data" << endl;
}
