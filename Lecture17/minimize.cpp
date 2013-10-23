#include <cmath>
#include <istream>
#include <iomanip>
using namespace std;

#include "nonlin.hpp"
using namespace cpt;

const double pi = 4 * atan(1.0);

double f(Matrix<double,1>& p) {
    double x = p[0], y = p[1];
    return x*x + 2 * y*y + 0.3 * cos(3 * pi * x) + 0.4 * cos(4 * pi * y);
}

void df(Matrix<double,1>& p, Matrix<double,1>& dp) {
    double x = p[0], y = p[1];
    dp[0] = 2 * x - 0.9 * pi * sin(3 * pi * x);
    dp[1] = 4 * y - 1.6 * pi * sin(4 * pi * y);
    return;
}

int main()
{
    cout << " Minimization using Broyden-Fletcher-Goldfarb-Shanno Algorithm\n"
         << " Find minimum of f(x,y) given an initial guess for x, y\n"
         << " Enter starting point coordinates x y: ";
    Matrix<double,1> p(2);
    cin >> p[0] >> p[1];
    cout << " Enter desired accuracy: ";
    double gtol;
    cin >> gtol;
    double f_min;
    int iterations;
    minimize_BFGS(p, gtol, iterations, f_min, f, df);
    cout << " f(" << p[0] << ", " << p[1] << ") = " << f_min
         << " in " << iterations << " iterations" << endl;
}
