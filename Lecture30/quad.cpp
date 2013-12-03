#include <cmath>
#include <iostream>
#include <string>
using namespace std;

#include "random.hpp"

double f(double x)
{
    return 1 / (1.0 + x * x);
}


double regular(double f(double), double a, double b, string rule)
{
    double x_mid = (a + b) / 2.0;
    double dx = b - a;
    if (rule == "Rectangle")
        return f(a) * dx;
    else if (rule == "Midpoint")
        return f(x_mid) * dx;
    else if (rule == "Trapezoid")
        return (f(a) + f(b)) * dx / 2.0;
    else if (rule == "Simpson's")
        return (f(a) + 4.0 * f(x_mid) + f(b)) * dx / 6.0;
    else {
        cerr << " regular rule " << rule << " not implemented" << endl;
	return 0;
    }
    return 0;
}

double integral( cpt::Random & rng, double exact,
		 double f(double), double a, double b, string rule,
		 int n_max=10000000, double max_error = 1.0e-6)
{
    int n = 1;
    double approx = 0.0;

    while (true) {
        if (n > n_max) {
            cout << " N exceeds " << n_max << " aborting ..." << endl;
            break;
        }
        approx = 0.0;
        if (rule == "Monte Carlo") {
            for (int i = 0; i < n; i++) {
                double x = a + (b - a) * rng.rand();
                approx += f(x);
            }
            approx *= (b - a) / double(n);
        } else {
            for (int i = 0; i < n; i++) {
                double d = (b - a) / double(n);
                double a_i = a + i * d;
                double b_i = a_i + d;
                approx += regular(f, a_i, b_i, rule);
            }
        }
        double error = approx - exact;
        cout << " " << n << "\t" << error << endl;
        if (abs(error) < max_error)
            break;
        n *= 2;
    }
    return approx;
}

int main()
{


    string rules[] = {
      "Rectangle", "Midpoint", "Trapezoid", "Simpson's", "Monte Carlo" };
  
    const double pi = 4 * atan(1.0);
  
    double a = 0.0;             // lower limit
    double b = 1.0;             // upper limit
    double exact = atan(b) - atan(a);

    cpt::Random rng;
    cout << " Quadrature of 1 / (1 + x^2) from " << a << " to " << b << "\n"
         << " Exact integral = " << exact << "\n"
         << " ----------------------------------" << endl;

    for (int i = 0; i < sizeof(rules) / sizeof(rules[0]); i++) {
        cout << " Quadrature rule: " << rules[i] << "\n"
             << " N    \tError\n";
        double answer = integral(rng, exact, f, a, b, rules[i]);
        cout << " Integral = " << answer << "\n"
             << " ----------------------------------" << endl;
    }
}
