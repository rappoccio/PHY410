#include "cptstd.hpp"
#include "diffeq.hpp"
using namespace cpt;

const double pi = 4 * atan(1.0);

double a = 1;                  // size of unit cell

// potential functions

double V_0 = -5;                // height of potential barrier
double Delta = 0.2;             // width of potential barrier

const int
    cutoff_coulomb = 0,         // coulomb with cutoff at ion core radius
    sine_step = 1,              // sine in first half, step in second
    kronig_penney = 2,          // as in kronig-penney.cpp
    free_particle = 3;

int potential_type;             // defaults to 0 = cutoff_coulomb;

const string potential_name[] = {
    "Cutoff Coulomb",
    "Sine Step",
    "Kronig-Penney",
    "Free Particle"
};

double V(double x)
{
    double dx = x - a / 2;      // displacement from center of cell

    if (potential_type == cutoff_coulomb) {

        dx = abs(dx);           // symmetric about center of cell
        if (dx < Delta / 2)     // inside ion core
            return V_0;         // potential is constant
        else
            return V_0 *        // value at core radius
              (Delta / 2) / dx; // Coulomb 1/r outside ion core

    } else if (potential_type == sine_step) {

        if (dx < a / 2) {       // left half of cell
            return V_0 *        // height of barrier
                sin(2 * pi * dx / a);  // inverted (dx < 0) half sine wave
        } else {                // right half of cell
            if (abs(dx - a / 4) // distance from center of barrier
                < Delta / 2)    // less than 1/2 barrier width
                return V_0;     // inside barrier
        return 0;               // outside barrier
        }
    } else if (potential_type == kronig_penney) {

        dx = abs(dx);           // symmetric about center of cell
        if (dx < Delta / 2)     // inside step
            return V_0;         // barrier height
        return 0;               // outside barrier
    }

    return 0;                   // default to free particle
}

double E;                       // current value of energy

// flow function for Runge-Kutta integration across unit cell

Matrix<double,1> flow(Matrix<double,1>& y)
{
    double x = y[0], phi = y[1], d_phi_dx = y[2];
    Matrix<double,1> f(3);
    f[0] = 1;
    f[1] = d_phi_dx;
    f[2] = 2 * (V(x) - E) * phi;
    return f;
}

void find_bloch_momentum(
    bool& in_band,              // set true if in band, false if in gap - output
    double& k)                  // value of Bloch momentum if in band - output
{
    // first fundamental solution
    Matrix<double,1> y_1(3);
    double x = 0, psi = 1, psi_prime = 0;
    y_1[0] = x;
    y_1[1] = psi;
    y_1[2] = psi_prime;
    double dx = 0.01;
    RK4_integrate(y_1, dx, flow, a);

    // second fundamental solution
    Matrix<double,1> y_2(3);
    x = 0, psi = 0, psi_prime = 1/a;
    y_2[0] = x;
    y_2[1] = psi;
    y_2[2] = psi_prime;
    RK4_integrate(y_2, dx, flow, a);

    // find k from dispersion relation
    double f = (y_1[1] + a * y_2[2]) / (2 * a);
    if (abs(f) <= 1) {
        in_band = true;
        k = acos(f);
    } else {
        in_band = false;
    }
}

int main()
{
    cout << " Band Structure in 1-D For General Potential\n"
         << " -------------------------------------------\n"
         << " Enter potential type 0, 1, 2, or 3: ";
    cin >> potential_type;
    cout << " Potential type = " << potential_name[potential_type] << endl;
    cout << " Enter E_min, E_max, and number of levels: ";
    double E_min, E_max;
    int n_levels;
    cin >> E_min >> E_max >> n_levels;

    string file_name("bloch.data");
    ofstream file("bloch.data");
    bool new_band = false;
    for (int n = 0; n < n_levels; n++) {
        E = E_min + n * (E_max - E_min) / double(n_levels - 1);
        double k;
        bool in_band;
        find_bloch_momentum(in_band, k);
        if (in_band && !new_band) {
            file << '\n';    // separate band values with blank line
            new_band = true;
        }
        if (in_band)
            file << k << '\t' << E << '\n';
        else
            new_band = false;
    }
    file.close();
    cout << " Levels in file " << file_name << endl;
}
