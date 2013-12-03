#include "cptstd.hpp"
#include "linalg.hpp"
#include "random.hpp"
using namespace cpt;

Matrix<double,1> model(      // model function to fit a linear model
  const double step)
{
    Matrix<double,1> f(2);
    f[0] = 1;
    f[1] = step;
    return f;
}

int main() {
    cout << " Random walk in 1 dimension\n"
         << " --------------------------\n"
         << " Enter number of walkers: ";
    int n_walkers;
    cin >> n_walkers;
    cout << " Enter number of steps: ";
    int n_steps;
    cin >> n_steps;

    // Random number generator
    Random rng;
    cout << " Using " << rng.get_algorithm()
         << " with seed " << rng.get_seed() << endl;

    // vector of walker positions initialized to zero
    Matrix<double,1> x(n_walkers);

    // data vectors
    Matrix<double,1> steps(n_steps);     // to save step number i (time)
    Matrix<double,1> x2ave(n_steps);     // to accumulate x^2 values
    Matrix<double,1> sigma(n_steps);     // to accumulate fluctuations in x^2

    // Loop over walkers
    for (int walker = 0; walker < n_walkers; walker++) {

        // Loop over number of steps
        for (int step = 0; step < n_steps; step++) {

            // take a random step
            if (rng.rand() < 0.5)
                x[walker] += 1;
            else
                x[walker] -= 1;

            // accumulate data
            steps[step] = step + 1;
            double xx = x[walker] * x[walker];
            x2ave[step] += xx;
            sigma[step] += xx * xx;
        }
    }

    // average the squared displacements and their variances
    for (int step = 0; step < n_steps; step++) {
        x2ave[step] /= n_walkers;
        sigma[step] /= n_walkers;
    }
    for (int step = 0; step < n_steps; step++) {
        sigma[step] = sqrt(sigma[step] - x2ave[step] * x2ave[step]);
        if (sigma[step] == 0)   // happens for first step!
            sigma[step] = 1;
        if (n_walkers > 1)
            sigma[step] /= sqrt(double(n_walkers - 1));
    }

    // fit data to a straight line
    double a, b, sigma_a, sigma_b, chisqr;
    chi_square_fit(steps, x2ave, sigma, a, b, sigma_a, sigma_b, chisqr);

    cout << " Fit to a straight line <x^2> = a + b n\n"
         << " Intercept a = " << a << " +- " << sigma_a << '\n'
         << " Slope     b = " << b << " +- " << sigma_b << '\n'
         << " Chisqr/dof  = " << chisqr / (n_steps - 2.0) << endl;

    // store in file for plotting
    ofstream file("rwalk.data");
    for (int step = 0; step < n_steps; step++)
        file << steps[step] << '\t' << x2ave[step]
             << '\t' << sigma[step] << '\n';
    file.close();
    cout << " t, <x^2>, sigma in file rwalk.data" << endl;
}
