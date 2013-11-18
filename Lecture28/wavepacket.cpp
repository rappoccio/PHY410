#include "wavepacket.hpp"

// define names for complex numbers and vectors
typedef complex<double> cdouble;
typedef Matrix<cdouble,1> cvector;

cvector chi;                    // complex wavefunction
cvector a, b, c;                // to represent tridiagonal elements of Q
cdouble alpha, beta;            // corner elements of Q

void initialize_chi_and_Q ()
{
    // create and initialize vectors to zero
    chi = cvector(N);
    a = b = c = cvector(N);

    // elements of tridiagonal matrix Q = (1/2)(1 + i dt H / (2 hbar))
    const cdouble i(0.0, 1.0);
    for (int j = 0; j < N; j++) {
        a[j] = - i * dt * h_bar / (8 * mass * dx * dx);
        b[j] = 0.5 + 1j * dt / (4 * h_bar) *
              (V(x[j]) + h_bar * h_bar / (mass * dx * dx));
        c[j] = a[j];
    }
    alpha = c[N-1];
    beta = a[0];
}

cvector solve_tridiagonal(cvector& a, cvector& b, cvector& c, cvector& r) {
    // solve Ax = r where A is tridiagonal with diagonals (a, b, c) and return x
    int n = r.size();
    cvector x(n), gama(n);
    cdouble beta = b[0];
    x[0] = r[0] / beta;
    for (int j = 1; j < n; j++) {
        gama[j] = c[j-1] / beta;
        beta = b[j] - a[j] * gama[j];
        x[j] = (r[j] - a[j] * x[j-1]) / beta;
    }
    for (int j = n-2; j >= 0; j--)
        x[j] -= gama[j+1] * x[j+1];
    return x;
}

cvector solve_tridiagonal_cyclic(
    cvector& a, cvector& b, cvector& c, cvector& r,
    complex<double> alpha, complex<double> beta)
{
    // solve Ax = r where A is tridagonal with corner elements alpha, beta
    int n = r.size();
    cvector x(n), b_prime(n), u(n), z(n);
    cdouble gama = -b[0];
    b_prime[0] = b[0] - gama;
    b_prime[n-1] = b[n-1] - alpha * beta / gama;
    for (int j = 1; j < n; j++)
        b_prime[j] = b[j];
    x = solve_tridiagonal(a, b_prime, c, r);
    u[0] = gama;
    u[n-1] = alpha;
    for (int j = 1; j < n-1; j++)
        u[j] = 0;
    z = solve_tridiagonal(a, b_prime, c, u);
    cdouble fact = x[0] + beta * x[n-1] / gama;
    fact /= 1.0 + z[0] + beta * z[n-1] / gama;
    for (int j = 0; j < n; j++)
        x[j] -= fact * z[j];
    return x;
}

void take_step () {              // time step using sparse matrix routines
    if (periodic)
        chi = solve_tridiagonal_cyclic(a, b, c, psi, alpha, beta);
    else
        chi = solve_tridiagonal(a, b, c, psi);
    for (int j = 0; j < N; j++)
        psi[j] = chi[j] - psi[j];
    t += dt;

}

int main() {

    cout << " Quantum Wavepacket Motion" << endl;
    initialize();
    initialize_chi_and_Q();

    ofstream file("potential.data");
    for (int i = 0; i < N; i++)
        file << x[i] << '\t' << V(x[i]) << '\n';
    file.close();
    cout << " saved V(x) in file potential.data" << endl;

    save_psi(0);
    int plots = 10;
    for (int plot = 1; plot <= plots; plot++) {
        double delta_t = 0;
        while (delta_t < L / (plots * velocity)) {
            take_step();
            delta_t += dt;
        }
        save_psi(plot);
    }

    // simple gnuplot animation
    int plot = 0;
    while (true) {
        fprintf(gnupipe, "set term %s\n", terminal.c_str());
        fprintf(gnupipe, "plot \"psi_%d.data\" w l\n", plot);
        plot = (plot+1) % (plots+1);
        fflush(gnupipe);
        time_t start_time = clock();
        while (true) {
            double secs = (clock() - start_time) / double(CLOCKS_PER_SEC);
            if (secs > 0.5)
                break;
        }
    }
}
