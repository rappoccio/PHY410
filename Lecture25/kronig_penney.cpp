#include "cptstd.hpp"

class KronigPenney {

public : 
  KronigPenney( double ia = 1.0, double iV0 = -5.0, double iDelta = 0.2) :
    a (ia), V_0 (iV0), Delta (iDelta)
  {
  }

  void solve_for_E (                  // to solve 2x2 eigenvalue problem
		    const double E,                 // desired energy E
		    complex<double> k[] )           // two possible solutions
  {

    double q = sqrt(2*E);
    double kappa = sqrt(2*(E-V_0));
    complex<double> i(0,1), T11, T12, T21, T22;

    T11 = exp(i * q * (a -  Delta)) / (4 * q * kappa)
      * (exp(i * kappa * Delta) * (q + kappa) * (q + kappa)
	 - exp(- i * kappa * Delta) * (q - kappa) * (q - kappa));
    T22 = conj(T11);
    T12 = - i * exp(i * q * (a - Delta)) / (2 * q * kappa)
      * (q * q - kappa * kappa) * sin(kappa  * Delta);
    T21 = conj(T12);

    // solve quadratic determinantal equation
    complex<double> b = - (T11 + T22);
    complex<double> c = (T11 *  T22 - T12 *  T21);
    k[0] = (- b + sqrt(b * b - 4.0 * c)) / 2.0;
    k[1] = (- b - sqrt(b * b - 4.0 * c)) / 2.0;
    for (int j = 0; j < 2; j++)
      k[j] = log(k[j])  / (i * a);
  }

  void compute_bands(
		     double dE,                      // step size in E for search
		     int steps,                      // number of steps
		     string band_file_name)          // output file name
  {
    const double pi = 4 * atan(1.0);
    ofstream file(band_file_name.c_str());
    double E = dE;
    for (int i = 0; i < steps; i++) {
      complex<double> q[2];
      solve_for_E(E, q);
      for (int j = 0; j < 2; j++) {
	double rq = real(q[j]);
	if (rq > 0 && rq < pi/a) {
	  file << rq << '\t' << E << '\n';
	  file << -rq << '\t' << E << '\n';
	}
      }
      E += dE;
    }
    file.close();
    cout << " Energy bands written in file " << band_file_name << endl;
  }

  double get_V0() const { return V_0;}
  double get_delta() const { return Delta;}
  void set_V0( double d = 0.0 ) { V_0 = d; }

protected : 


  double a;                       // size of unit cell - lattice spacing
  double V_0;                    // height of potential barrier
  double Delta;                 // width of potential barrier

};

int main()
{
    cout << " Kronig-Penney Model" << endl;
    double dE = 0.01;
    int steps = 3000;

    KronigPenney kp;

    cout << " V_0 = " << kp.get_V0() << " Delta = " << kp.get_delta() << endl;
    kp.compute_bands(dE, steps, "band.data");
    kp.set_V0(0.0);
    cout << " V_0 = " << kp.get_V0() << " Delta = " << kp.get_delta() << endl;
    kp.compute_bands(dE, steps, "band0.data");
}
