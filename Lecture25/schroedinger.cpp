#include "cptstd.hpp"
#include "basalg.hpp"
#include "matrix.hpp"
using namespace cpt;

class Schroedinger {

public : 

  Schroedinger(double iE, 
	       int iN=500, double ixl = -5.0, double ixr = 5.0,
	       double iomega = 1.0,
	       double ihbar = 1.0, double im = 1.0) :
    hbar(ihbar), m(im), omega(iomega), E(iE),
    N(iN), x_left(ixl), x_right(ixr), 
    h( (x_right-x_left) / N),
    phi_left(N+1),phi_right(N+1),phi(N+1)
  {

  }

  double V(double x)                  // harmonic oscillator potential
  {
    return 0.5 * m * omega * omega * x * x;
  }


  double q(double x)                  // Sturm-Liouville q function
  {
    return 2 * m / (hbar * hbar) * (E - V(x));
  }


  double F(double iE)                  // eigenvalues at F(E) = 0
  {

    // set global E value needed by the q(x) function
    E = iE;

    // find right turning point
    int i_match = N;
    double x = x_right;             // start at right boundary
    while (V(x) > E)  {             // in forbidden region
      --i_match;
      x -= h;
      if (i_match < 0) {
	cerr << "can't find right turning point" << endl;
	exit(EXIT_FAILURE);
      }
    }

    // integrate phi_left using Numerov algorithm
    phi_left[0] = 0;
    phi_left[1] = 1e-10;
    double c = h * h / 12;          // constant in Numerov formula
    for (int i = 1; i <= i_match; i++) {
      x = x_left + i * h;
      phi_left[i+1]  = 2 * (1 - 5 * c * q(x)) * phi_left[i];
      phi_left[i+1] -= (1 + c * q(x - h)) * phi_left[i-1];
      phi_left[i+1] /= 1 + c * q(x + h);
    }

    // integrate phi_right
    phi[N]   = phi_right[N]   = 0;
    phi[N-1] = phi_right[N-1] = 1e-10;
    for (int i = N - 1; i >= i_match; i--) {
      x = x_right - i * h;
      phi_right[i-1]  = 2 * (1 - 5 * c * q(x)) * phi_right[i];
      phi_right[i-1] -= (1 + c * q(x + h)) * phi_right[i+1];
      phi[i-1] = phi_right[i-1] /= 1 + c * q(x - h);
    }

    // rescale phi_left
    double scale = phi_right[i_match] / phi_left[i_match];
    for (int i = 0; i <= i_match + 1; i++)
      phi[i] = phi_left[i] *= scale;

    // make F(E) continuous
    static int sign = 1;            // current sign used
    static int nodes = 0;           // current number of nodes

    // count number of nodes in phi_left
    int n = 0;
    for (int i = 1; i <= i_match; i++)
      if (phi_left[i-1] * phi_left[i] < 0)
	++n;

    // flip its sign when a new node develops
    if (n != nodes) {
      nodes = n;
      sign = -sign;
    }

    return sign * ( phi_right[i_match-1] - phi_right[i_match+1]
		    - phi_left[i_match-1] + phi_left[i_match+1] )
      / (2 * h * phi_right[i_match]);
  }

  void normalize() {
    double norm = 0;
    for (int i = 0; i < N; i++)
      norm += phi[i] * phi[i];
    norm /= N;
    norm = sqrt(norm);
    for (int i = 0; i < N; i++)
      phi[i] /= norm;
  }

  void set_E(double iE) { E = iE;}

  int get_N() const { return N;}
  double get_x_left() const { return x_left;}
  double get_x_right() const { return x_right;}
  double get_h() const { return h; }

  double get_phi(int i) const { return phi[i]; }

protected : 

  double hbar ;               // Planck's constant / 2pi
  double m ;                  // particle mass
  double omega ;              // oscillator frequency

  double E;                   // current energy in search

  int N ;                     // number of lattice points = N + 1
  double x_left;              // left boundary
  double x_right;             // right boundary
  double h ;                  // grid spacing

  Matrix<double,1> phi_left;  // wave function integrating from left
  Matrix<double,1> phi_right; // wave function integrating from right
  Matrix<double,1> phi;       // whole wavefunction


};

class ScroedingerWrapFE  {
public : 
  ScroedingerWrapFE( Schroedinger & inputS) : s(inputS) {}

  double operator()( double x ) {
    return s.F(x);
  }

protected : 
  Schroedinger & s;

};


class ScroedingerWrapQ  {
public : 
  ScroedingerWrapQ( Schroedinger & inputS) : s(inputS) {}

  double operator()( double x ) {
    return s.q(x);
  }

protected : 
  Schroedinger & s;

};

int main() {
    cout << " Eigenvalues of the Schroedinger equation\n"
         << " for the harmonic oscillator V(x) = 0.5 x^2\n"
         << " ------------------------------------------\n"
         << " Enter maximum energy E: ";
    double E_max;
    cin >> E_max;


    // find the energy levels
    cout << "\n Level       Energy       Simple Steps   Secant Steps"
         << "\n -----   --------------   ------------   ------------\n";
    int level = 0;                      // level number
    double E_old = 0;                   // previous energy eigenvalue
    double E = 0.1;                            // guess an E below the ground state

    Schroedinger s(E);                 // Main Schroedinger object
    // Note : these two are holding active, non-const references to s! 
    ScroedingerWrapFE s_wrapfe( s );   // Wrapper class to call Schroedinger::F(E) from operator()
    ScroedingerWrapQ s_wrapq( s );     // Wrapper class to call Schroedinger::q(E) from operator()



    ofstream levels_file("levels.data");
    ofstream phi_file("phi.data");

    // draw the potential
    for (int i = 0; i <= s.get_N(); i++) {
      double x = s.get_x_left() + i * s.get_h();
      levels_file << x << '\t' << s.V(x) << '\n';
    }
    levels_file << '\n';

    do {


        // estimate next E and dE
        double dE = 0.5 * (E - E_old);
        E_old = E;
        E += dE;

	s.set_E( E );   // Remember, this is referenced below by s_wrapfe and s_wrapq. 


        // Simple search to locate next energy level approximately
        SimpleSearchT<ScroedingerWrapFE>  simpleFE;
        double accuracy = 0.01;         // set a relatively low accuracy
        simpleFE.set_accuracy(accuracy);
        simpleFE.set_first_root_estimate(E);
        simpleFE.set_step_estimate(dE);
        E = simpleFE.find_root(s_wrapfe);

        // Secant search to locate energy level precisely
        accuracy = 0.000001;            // can use a relatively high accuracy
        double E1 = E + 100 * accuracy; // guess second point required by secant
        SecantSearchT<ScroedingerWrapFE> secant;
        secant.set_accuracy(accuracy);
        secant.set_first_root_estimate(E);
        secant.set_step_estimate(dE);
	s.set_E(E);
        E = secant.find_root(s_wrapfe);


        // Simple search to locate zeros of Schroedinger::q
        SimpleSearchT<ScroedingerWrapQ>  simpleQ;
        simpleQ.set_accuracy(accuracy);

	
        cout.precision(10);
        cout << setw(4) << level++ << "  "
             << setw(15) << E << "  " << endl;
        simpleQ.set_first_root_estimate(s.get_x_left());
        simpleQ.set_step_estimate(s.get_h());
        levels_file << simpleQ.find_root(s_wrapq) << '\t' << E << '\n';
        simpleQ.set_first_root_estimate(s.get_x_right());
        simpleQ.set_step_estimate(-s.get_h());
        levels_file << simpleQ.find_root(s_wrapq) << '\t' << E << '\n';
        levels_file << '\n';
        s.normalize();
        for (int i = 0; i <= s.get_N(); i++) {
	  double x = s.get_x_left() + i * s.get_h();
	  phi_file << x << '\t' << s.get_phi(i) << '\n';
        }
        phi_file << '\n';
    } while (E < E_max);

    levels_file.close();
    phi_file.close();

    // output the search function to a file
    ofstream search_file("F.data");
    E = 0.1;
    double dE = 0.01;
    while (E < E_max) {
        search_file << E << '\t' << s.F(E) << '\n';
        E += dE;
    }
    search_file.close();

    cout << "\n Energy levels in file levels.data"
         << "\n Eigenfunctions in file phi.data"
         << "\n Search function in file F.data" << endl;
}
