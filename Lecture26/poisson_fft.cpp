#include "cptstd.hpp"
#include "linalg.hpp"
#include "basalg.hpp"
using namespace cpt;

int main() {

    cout << " FFT solution of Poisson's equation\n"
         << " ----------------------------------\n";
    cout << " Enter number points in x or y: ";
    int N;
    cin >> N;
    int power_of_2 = 1;
    while (power_of_2 < N)
        power_of_2 *= 2;
    if (power_of_2 != N)
        cout << " Warning: N is not a power of 2 - FFT may be slow!" << endl;

    double h = 1 / double(N - 1);        // assume physical size in x and y = 1

    clock_t t0 = clock();

    double q = 10;                       // point charge
    Matrix<complex<double>,2> rho(N, N);
    for (int j = 0; j < N; j++) {
        for (int k = 0; k < N; k++) {
            if (j == N/2 && k == N/2)    // at center of lattice
                rho[j][k] = q / (h * h);
            else
                rho[j][k] = 0.0;
        }
    }

    // FFT rows of rho
    FFT fft;
    std::vector<std::complex<double> > f(N);      // to store rows and columns
    for (int j = 0; j < N; j++) {
        for (int k = 0; k < N; k++)
            f[k] = rho[j][k];
        fft.transform(f);
        for (int k = 0; k < N; k++)
            rho[j][k] = f[k];
    }

    // FFT columns of rho
    for (int k = 0; k < N; k++) {
        for (int j = 0; j < N; j++)
            f[j] = rho[j][k];
        fft.transform(f);
        for (int j = 0; j < N; j++)
            rho[j][k] = f[j];
    }

    // solve equation in Fourier space
    Matrix<complex<double>,2> V(N, N);
    complex<double> i(0.0, 1.0);
    double pi = 4 * atan(1.0);
    complex<double> W = exp(2.0 * pi * i / double(N));
    complex<double> W_m = 1.0, W_n = 1.0;
    for (int m = 0; m < N; m++) {
        for (int n = 0; n < N; n++) {
            complex<double> denom = 4.0;
            denom -= W_m + 1.0 / W_m + W_n + 1.0 / W_n;
            if (denom != 0.0)
                V[m][n] = rho[m][n] * h * h / denom;
            W_n *= W;
        }
        W_m *= W;
    }

    // inverse FFT rows of V
    for (int j = 0; j < N; j++) {
        for (int k = 0; k < N; k++)
            f[k] = V[j][k];
        fft.inverse_transform(f);
        for (int k = 0; k < N; k++)
            V[j][k] = f[k];
    }
    // inverse FFT columns of V
    for (int k = 0; k < N; k++) {
        for (int j = 0; j < N; j++)
            f[j] = V[j][k];
        fft.inverse_transform(f);
        for (int j = 0; j < N; j++)
            V[j][k] = f[j];
    }

    clock_t t1 = clock();
    cout << " CPU time = " << double(t1 - t0) / CLOCKS_PER_SEC
         << " sec" << endl;




    // write potential to file
    cout << " Potential in file poisson_fft.data" << endl;
    ofstream date_file("poisson_fft.data");
    for (int i = 0; i < N; i++) {
      double x = i * h;
      for (int j = 0; j < N; j++) {
	double y = j * h;
	char buff[1000];
	sprintf(buff, "%12.6f %12.6f %12.6f", x, y, real(V[i][j]) );
	date_file << buff << endl;
      }
    }
    date_file.close();

}
