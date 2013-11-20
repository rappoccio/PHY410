#include "cptstd.hpp"
#include "matrix.hpp"
using namespace cpt;

#include <cstdio>
#ifdef _WIN32
    string gnuplot = "env pgnuplot.exe ";
    string terminal = "windows";
    FILE *gnupipe = _popen(gnuplot.c_str(), "w");
#else
    string gnuplot = "/usr/bin/env gnuplot ";
    string terminal = "x11";
    FILE *gnupipe = popen(gnuplot.c_str(), "w");
#endif

const double pi = 4 * atan(1.0);

double L = 5.0;                 // system size
int N;                          // number of grid intervals
double dx;                      // grid spacing
double c = 1;                   // wave speed
double t;                       // time
double dt;                      // time step
void (*method)();               // integration algorithm
Matrix<double,1> x;             // grid points
Matrix<double,1> u, u_new;      // wave amplitude
int step_number;                // step number
bool step_wave;                 // otherwise Gaussian modulate cosine

double f_0(double x) {        // initial waveform

    double x_0 = L / 2;
    double sigma = 0.1 * L;   // width of initial waveform

    if (step_wave) {
        if (abs(x - x_0) < sigma)
            return 1;
        else
            return 0;
    } else {
        double k = pi / sigma;
        double gaussian = exp(- (x - x_0) * (x - x_0) / (2 * sigma * sigma));
        return cos(k * (x - x_0)) * gaussian;
    }
}

void initialize() {

    dx = L / double(N);         // grid spacing

    // initialize the spatial grid
    x = Matrix<double,1>(N+1);
    for (int i = 0; i <= N; i++)
        x[i] = i * dx;

    // initialize the solution vector
    u = Matrix<double,1>(N+1);
    u_new = Matrix<double,1>(N+1);
    for(int i = 0; i <= N; i++)
        u[i] = f_0(x[i]);

    t = 0;
    step_number = 0;
}

void FTCS() {
    for (int i = 0; i <= N; i++) {
        int i_plus_1 = i < N ? i + 1 : 0;
        int i_minus_1 = i > 0 ? i - 1 : N;
        u_new[i] = u[i] - c * dt / (2 * dx) * (u[i_plus_1] - u[i_minus_1]);
    }
}

void Lax() {
    for (int i = 0; i <= N; i++) {
        int i_plus_1 = i < N ? i + 1 : 0;
        int i_minus_1 = i > 0 ? i - 1 : N;
        u_new[i] = (u[i_plus_1] + u[i_minus_1]) / 2
                  - c * dt / (2 * dx) * (u[i_plus_1] - u[i_minus_1]);
    }
}

void Lax_Wendroff() {
    double D = c * c * dt * dt / (2 * dx * dx);
    for (int i = 0; i <= N; i++) {
        int i_plus_1 = i < N ? i + 1 : 0;
        int i_minus_1 = i > 0 ? i - 1 : N;
        u_new[i] = u[i] - c * dt / (2 * dx) * (u[i_plus_1] - u[i_minus_1])
                   + D * (u[i_plus_1] + u[i_minus_1] - 2 * u[i]);
    }
}

void take_step() {
    method();
    u = u_new;
    t += dt;
    ++step_number;
}

void save_u(int plot_number) {
    ostringstream os;
    os  << "u_" << plot_number << ".data";
    string file_name(os.str());
    ofstream file(file_name.c_str());
    if (file)
        cout << " writing " << file_name << endl;
    else
        cerr << " cannot open " << file_name << endl;
    for (int i = 0; i < N; i++)
        file << x[i] << '\t' << u[i] << '\n';
    file.close();
}

int main()
{
    cout << " Finite Difference solution of the Advection Equation\n"
         << " Choose a numerical method: 1) FTCS, 2) Lax, 3) Lax-Wendroff : ";
    int m;
    cin >> m;
    method = FTCS;
    if (m == 2)
        method = Lax;
    if (m == 3)
        method = Lax_Wendroff;
    cout << " Enter number of grid cells: ";
    cin >> N;
    cout << " Time for wave to move one grid spacing is " << L/(N*c) << endl;
    cout << " Enter time step dt: ";
    cin >> dt;

    initialize();
    int plots = 10;
    save_u(0);
    for (int plot = 0; plot < plots; plot++) {
        double delta_t = 0;
        while (delta_t < L / (c * plots)) {
            take_step();
            delta_t += t;
        }
        save_u(plot+1);
    }

    // simple Gnuplot pipes animation
    cout << " Enter animation time: ";
    double t_max;
    cin >> t_max;
    initialize();
    double frame_rate = 30;
    double dt_frame = 1 / frame_rate;
    int steps_per_frame = max(1, int(dt_frame / dt));
    while (t < t_max) {
        ostringstream os;
        os << "plot \'-\' with impulses\n";
        for (int i = 0; i < N; i++)
            os << x[i] << '\t' << u[i] << '\n';
        os << "e\n";
        fprintf(gnupipe, "%s", os.str().c_str());
        fflush(gnupipe);
        time_t start_time = clock();
        for (int step = 0; step < steps_per_frame; step++) {
            take_step();
            t += dt;
        }
        while (true) {
            double secs = (clock() - start_time) / double(CLOCKS_PER_SEC);
            if (secs > dt_frame)
                break;
        }
    }
    fclose(gnupipe);
}
