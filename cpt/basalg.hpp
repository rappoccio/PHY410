#ifndef CPT_BASALG_HPP
#define CPT_BASALG_HPP

#include <complex>
#include <iostream>
#include <vector>

namespace cpt {

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

  class SimpleSearch : public RootFinder {
  public:
    SimpleSearch() { max_steps = 1000; }
    double find_root(double f(double x));
    double find_root(double f(double x), double x_guess, double dx_guess) {
      set_first_root_estimate(x_guess);
      set_step_estimate(dx_guess);
      return find_root(f);
    }
  };

  class BisectionSearch : public RootFinder {
  public:
    BisectionSearch() { max_steps = 100; }
    double find_root(double f(double x));
    double find_root(
      double f(double x), double x_guess, double x_second_guess
    ) {
      set_first_root_estimate(x_guess);
      set_second_root_estimate(x_second_guess);
      return find_root(f);
    }
  };

  class SecantSearch : public RootFinder {
  public:
    SecantSearch() { max_steps = 30; }
    double find_root(double f(double x));
    double find_root(
      double f(double x), double x_guess, double x_second_guess
    ) {
      set_first_root_estimate(x_guess);
      set_second_root_estimate(x_second_guess);
     return find_root(f);
    }
  };

  class TangentSearch : public RootFinder {
  public:
    TangentSearch() { max_steps = 20; }
    double find_root(double f(double x), double f_prime(double x));
    double find_root(
      double f(double x), double f_prime(double x), double x_guess
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

} /* namespace cpt */

#endif /* CPT_BASALG_HPP  */
