#include "wavepacket.hpp"

Matrix<complex<double>,1>
  T_exp_factor, V_exp_factor;   // precomputed phase rotations

void initialize_phases()
{
    // initialize the phase rotation factors
    T_exp_factor = Matrix<complex<double>,1>(N);
    V_exp_factor = Matrix<complex<double>,1>(N);

    for (int j = 0; j < N; j++) {
        // kinetic factor exp[-iT/h_bar dt]
        const double pi = 4 * atan(1.0);
        double p = j < N / 2 ? j : j - N;
        p *= h_bar * 2 * pi / L;
        double theta = - p * p / (2 * mass) / h_bar * dt;
        T_exp_factor[j] = complex<double>(cos(theta), sin(theta));

        // potential factor exp[-iV(x)/(2h_bar) dt]
        theta = - V(x[j]) / 2 / h_bar * dt;
        V_exp_factor[j] = complex<double>(cos(theta), sin(theta));
    }
}

void take_step() {

    // first half potential phase rotation
    for (int j = 0; j < N; j++)
        psi[j] *= V_exp_factor[j];

    // FFT to momentum space
    fft(psi);

    // kinetic phase rotation
    for (int j = 0; j < N; j++)
        psi[j] *= T_exp_factor[j];

    // FFT back to position space
    bool do_inverse = true;
    fft_inv(psi);

    // second half potential phase rotation
    for (int j = 0; j < N; j++)
        psi[j] *= V_exp_factor[j];

    t += dt;
}

int main() {

    cout << " Quantum Wavepacket Motion using FFT" << endl;
    initialize();
    initialize_phases();

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
