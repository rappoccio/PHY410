#ifndef CPT_BASALG_HPP
#define CPT_BASALG_HPP

#include <complex>
#include <iostream>
#include <vector>
#include <cassert>
#include <limits>
#include <fstream>
#include <iomanip>

namespace cpt {

// Class template to implemlent Ridder's algorithm for differentiation. 
template< typename T>
double diff_Ridders (
		      T const & f,                // name of function to be differentiated, satisfies double operator()(double)
		      double x,                   // input: the point at which df/dx is required
		      double h,                   // input: suggestion for an initial step size
		      double& error)              // output: estimate of error by algorithm
{
  assert( h != 0.0 );
  const int n = 10;           // dimension of extrapolation table
  double a[n][n];             // extrapolation table
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      a[i][j] = 0;
  a[0][0] = (f(x + h) - f(x - h)) / (2 * h);
  double answer = 0;
  error = std::numeric_limits<double>::max();
  for (int i = 0; i < n; i++) {
    h /= 1.4;
    a[0][i] = (f(x + h) - f(x - h)) / (2 * h);
    double fac = 1.4 * 1.4;
    for (int j = 1; j <= i; j++) {
      a[j][i]=(a[j-1][i] * fac - a[j-1][i-1]) / (fac - 1);
      fac *= 1.4 * 1.4;
      double err = std::max(std::abs(a[j][i] - a[j-1][i]),
			    std::abs(a[j][i] - a[j-1][i-1]));
      if (err <= error) {
	error = err;
	answer = a[j][i];
      }
    }
    if (std::abs(a[i][i] - a[i-1][i-1]) >= 2 * error)
      break;
  }
  return answer;
}

//Function template to approximate the definite integral of f from a to b by
//the composite trapezoidal rule, using n subintervals.
//From http://en.wikipedia.org/wiki/Trapezoidal_rule
template< typename T>
double trapezoid(            // performs trapezoid quadrature
    T const & f,             // function to be integrated
    double a,                // lower limit of integration
    double b,                // upper limit of integration
    unsigned n)              // number of subintervals
{
  assert ( n > 0 );
  double h = (b-a) / n;
  double s = f(a) + f(b);
  for ( unsigned i = 1; i < n; ++i ) {
    s += 2 * f(a + i*h);
  }
  return s * h / 2.0;
}

// Function template to implemlent the adaptive trapezoid algorithm for quadrature
template < typename T> 
double adaptive_trapezoid(                 // performs iterative trapezoid quadrature
		  T const & f,             // name of function to be differentiated, satisfies double operator()(double)
		  double a,                // lower limit of integration
		  double b,                // upper limit of integration
		  double acc,              // desired accuracy
		  bool output=false)       // output on each iteration if output == true
{
  double old_sum = -1e-30;
  double h = b - a;
  int n = 1;
  double sum = (f(a) + f(b)) / 2;
  if (output)
    std::cout << "N = " << n+1 << ",  Integral = " << h*sum << std::endl;
  while (std::abs(h * (old_sum - sum/2)) > acc) {
    old_sum = sum;
    for (int i = 0; i < n; i++)
      sum += f(a + (i + 0.5) * h);
    n *= 2;
    h /= 2;
    if (output)
      std::cout << "N = " << n+1 << ",  Integral = " << h*sum << std::endl;
  }
  return h * sum;
}


// Function template to approximate the definite integral of f from a to b
// by the composite Simpson's rule, using n subintervals
template< typename T>
double simpson(              // performs iterative trapezoid quadrature
    T const & f,             // function to be integrated
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

//      Declare some utility functions


//      Utility functions

static void printHead(const char *algorithm, double accuracy, std::ostream& os )
{
        os << "\n ROOT FINDING using " << algorithm
           << "\n Requested accuracy = " << accuracy
           << "\n Step     Guess For Root          Step Size           Function Value"
           << "\n ----  --------------------  --------------------  --------------------"
           << std::endl;
}

static void printStep(int step, double x, double dx,
                                    double f_of_x, std::ostream& os)
{
        int w = os.width();
        int p = os.precision();
	std::ios::fmtflags f = os.flags();
        os.setf(std::ios::right, std::ios::adjustfield);
        os << " " << std::setw(4) << step << "  ";
        os.setf(std::ios::left, std::ios::adjustfield);
        os << std::setprecision(14)
           << std::setw(20) << x << "  "
           << std::setw(20) << dx << "  "
           << std::setw(20) << f_of_x
           << std::endl;
        os.width(w);
        os.precision(p);
        os.setf(f);
}

static void printWarning(int max_steps)
{
        std::cerr << " Warning: maximum number of steps "
                << max_steps << " exceeded!" << std::endl;
}




  // This class hierarchy implements elementary root finding algorithms

  class RootFinder {          // base class for other algorithms
  public:
    RootFinder() {            // constructor sets default parameters
      x0 = 0;                 // initial guess
      x1 = 1;                 // second guess
      dx = 1;                 // step size
      accuracy = 1.0e-6;      // desired accuracy
      max_steps = 20;         // to check runaway algorithms
      verbose = false;        // false/true turns printing off/on
      os = &std::cout;        // pointer to print stream
    }

    void set_first_root_estimate(
      double x_guess
    );

    void set_step_estimate
      (double step_size
    );

    void set_second_root_estimate(
      double x_guess
    );

    bool bracket_root(
      double f(double x)
    );

    void set_accuracy(
      double epsilon
    );

    double get_accuracy() {
      return accuracy;
    }

    int get_steps() {
      return steps;
    }

    void set_max_steps(
      int steps
    );

    int get_max_steps() {
      return max_steps;
    }

    void print_steps(bool yn=true) {
        verbose = yn;
    }

    void toggle_printing() {
      verbose = verbose ? false : true;
    }

    void set_print_stream(
      std::ostream& ostr
    ) {
      os = &ostr;
    }

  protected:
    double x0;
    double x1;
    double dx;
    double accuracy;
    int steps;
    int max_steps;
    bool verbose;
    std::ostream *os;
  };

  template< typename T>
  class SimpleSearchT : public RootFinder {
  public:
    SimpleSearchT() { max_steps = 1000; }
    double find_root(T const & f);
    double find_root(T const & f, double x_guess, double dx_guess) {
      set_first_root_estimate(x_guess);
      set_step_estimate(dx_guess);
      return find_root(f);
    }
  };

  template< typename T>
  class BisectionSearchT : public RootFinder {
  public:
    BisectionSearchT() { max_steps = 100; }
    double find_root(T const & f);
    double find_root(
      T const & f, double x_guess, double x_second_guess
    ) {
      set_first_root_estimate(x_guess);
      set_second_root_estimate(x_second_guess);
      return find_root(f);
    }
  };

  template< typename T>
  class SecantSearchT : public RootFinder {
  public:
    SecantSearchT() { max_steps = 30; }
    double find_root(T const & f);
    double find_root(
      T const & f, double x_guess, double x_second_guess
    ) {
      set_first_root_estimate(x_guess);
      set_second_root_estimate(x_second_guess);
     return find_root(f);
    }
  };

  template< typename T>
  class TangentSearchT : public RootFinder {
  public:
    TangentSearchT() { max_steps = 20; }
    double find_root(T const & f, T const & f_prime);
    double find_root(
      T const & f, T const & f_prime, double x_guess
    ) {
      set_first_root_estimate(x_guess);
      return find_root(f, f_prime);
    }
  };

  // This class implements the FFT on a vector of N data points

  class FFT {
  public:

    FFT();                          // constructor

    void transform(                 // replaces data by its transform
      std::vector< std::complex<double> > &data
    );

    void inverse_transform(         // replaces data by its inverse transform
      std::vector< std::complex<double> > &data
    );

    std::vector<double> power(      // returns power spectrum of data
      const std::vector< std::complex<double> > &data
    );

  private:
    int N;
    std::vector< std::complex<double> > *f;
    bool inverse;

    void bit_reverse();
    void Danielson_Lanczos(int n);
  };


#include "basalg.icc"

} /* namespace cpt */


#endif /* CPT_BASALG_HPP  */
