

#include <iostream>
#include <math.h>
#include <assert.h>

// Approximates the definite integral of f from a to b
// by the composite Simpson's rule, using n subintervals
double simpson(              // performs iterative trapezoid quadrature
    double (*f)(double),     // function to be integrated
    double a,                // lower limit of integration
    double b,                // upper limit of integration
    unsigned n)              // number of subintervals
{

  assert(n > 0);
  double h = (b-a) / n;
  double s = f(a) + f(b);

  for (unsigned i = 1; i < n; i += 2 ) {
    s += 4 * f(a + i * h);
  }
  for (unsigned i = 2; i < n-1; i += 2 ) {
    s += 2 * f(a + i * h);
  }
  return s * h / 3.0;

}
