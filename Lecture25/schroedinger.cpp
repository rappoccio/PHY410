#include "cptstd.hpp"
#include "basalg.hpp"
#include "matrix.hpp"
using namespace cpt;

const double hbar = 1;              // Planck's constant / 2pi
const double m = 1;                 // particle mass
double omega = 1;                   // oscillator frequency

double V(double x)                  // harmonic oscillator potential
{
    return 0.5 * m * omega * omega * x * x;
}

double E;                           // current energy in search

double q(double x)                  // Sturm-Liouville q function
{
    return 2 * m / (hbar * hbar) * (E - V(x));
}

int N = 500;                        // number of lattice points = N + 1
double x_left = -5;                 // left boundary
double x_right = 5;                 // right boundary
double h = (x_right - x_left) / N;  // grid spacing

Matrix<double,1> phi_left(N+1);     // wave function integrating from left
Matrix<double,1> phi_right(N+1);    // wave function integrating from right
Matrix<double,1> phi(N+1);          // whole wavefunction

double F(double E)                  // eigenvalues at F(E) = 0
{

    // set global E value needed by the q(x) function
    ::E = E;

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

int main() {
    cout << " Eigenvalues of the Schroedinger equation\n"
         << " for the harmonic oscillator V(x) = 0.5 x^2\n"
         << " ------------------------------------------\n"
         << " Enter maximum energy E: ";
    double E_max;
    cin >> E_max;

    ofstream levels_file("levels.data");
    ofstream phi_file("phi.data");

    // draw the potential
    for (int i = 0; i <= N; i++) {
        double x = x_left + i * h;
        levels_file << x << '\t' << V(x) << '\n';
    }
    levels_file << '\n';

    // find the energy levels
    cout << "\n Level       Energy       Simple Steps   Secant Steps"
         << "\n -----   --------------   ------------   ------------\n";
    int level = 0;                      // level number
    double E_old = 0;                   // previous energy eigenvalue
    E = 0.1;                            // guess an E below the ground state
    do {

        // estimate next E and dE
        double dE = 0.5 * (E - E_old);
        E_old = E;
        E += dE;

        // Simple search to locate next energy level approximately
        SimpleSearch simple;
        double accuracy = 0.01;         // set a relatively low accuracy
        simple.set_accuracy(accuracy);
        simple.set_first_root_estimate(E);
        simple.set_step_estimate(dE);
        E = simple.find_root(F);
        int simple_steps = simple.get_steps();

        // Secant search to locate energy level precisely
        accuracy = 0.000001;            // can use a relatively high accuracy
        double E1 = E + 100 * accuracy; // guess second point required by secant
        SecantSearch secant;
        secant.set_accuracy(accuracy);
        secant.set_first_root_estimate(E);
        secant.set_step_estimate(dE);
        E = secant.find_root(F);
        int secant_steps = secant.get_steps();

        cout.precision(10);
        cout << setw(4) << level++ << "  "
             << setw(15) << E << "  "
             << setw(10) << simple_steps << "     "
             << setw(10) << secant_steps << endl;
        simple.set_first_root_estimate(x_left);
        simple.set_step_estimate(h);
        levels_file << simple.find_root(q) << '\t' << E << '\n';
        simple.set_first_root_estimate(x_right);
        simple.set_step_estimate(-h);
        levels_file << simple.find_root(q) << '\t' << E << '\n';
        levels_file << '\n';
        normalize();
        for (int i = 0; i <= N; i++) {
            double x = x_left + i * h;
            phi_file << x << '\t' << phi[i] << '\n';
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
        search_file << E << '\t' << F(E) << '\n';
        E += dE;
    }
    search_file.close();

    cout << "\n Energy levels in file levels.data"
         << "\n Eigenfunctions in file phi.data"
         << "\n Search function in file F.data" << endl;
}
