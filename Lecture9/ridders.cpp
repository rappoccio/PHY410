// This file basic.hpp provides functions for basic numerical algorithms

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

double diff_Ridders(
    double (*f)(double),        // name of function to be differentiated
    double x,                   // input: the point at which df/dx is required
    double h,                   // input: suggestion for an initial step size
    double& error)              // output: estimate of error by algorithm
{
    if (h == 0.0) {
        cerr << "diff_Ridders: h must be non-zero\n";
        exit(1);
    }
    const int n = 10;           // dimension of extrapolation table
    double a[n][n];             // extrapolation table
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            a[i][j] = 0;
    a[0][0] = (f(x + h) - f(x - h)) / (2 * h);
    double answer = 0;
    error = numeric_limits<double>::max();
    for (int i = 0; i < n; i++) {
        h /= 1.4;
        a[0][i] = (f(x + h) - f(x - h)) / (2 * h);
        double fac = 1.4 * 1.4;
        for (int j = 1; j <= i; j++) {
            a[j][i]=(a[j-1][i] * fac - a[j-1][i-1]) / (fac - 1);
            fac *= 1.4 * 1.4;
            double err = max(abs(a[j][i] - a[j-1][i]),
                             abs(a[j][i] - a[j-1][i-1]));
            if (err <= error) {
                error = err;
                answer = a[j][i];
            }
        }
        if (abs(a[i][i] - a[i-1][i-1]) >= 2 * error)
            break;
    }
    return answer;
}

double test_ridders_x2( double x ) { return x*x; }


int main(int argc, char ** argv) 
{
  cout <<  "Testing Ridder's algorithm for computing numeric derivatives" << endl;

  double error = 0.0;
  double h = 0.01;

  char buff[1000];
  for ( double x = 2.0; x <= 5.0 ; x+= 1.0) {
    double xprime = diff_Ridders( std::exp, x, h, error );
    sprintf( buff, "%2.0f : True value = %38.32f, Ridder = %38.32f, exp error = %38.32f, obs error = %38.32f", x, std::exp(x), xprime, h*h*h*h, error );
    cout << buff << endl;

    xprime = diff_Ridders( std::sin, x, h, error );
    sprintf( buff, "%2.0f : True value = %38.32f, Ridder = %38.32f, sin error = %38.32f, obs error = %38.32f", x, std::cos(x), xprime, h*h*h*h, error );
    cout << buff << endl;


    xprime = diff_Ridders( test_ridders_x2, x, h, error );
    sprintf( buff, "%2.0f : True value = %38.32f, Ridder = %38.32f, sin error = %38.32f, obs error = %38.32f", x, 2*x, xprime, h*h*h*h, error );
    cout << buff << endl;

  }

  return 0;

}
