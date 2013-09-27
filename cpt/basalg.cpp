#include "basalg.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>
using namespace std;

namespace cpt {

//      Declare some utility functions

static void printHead(const char *algorithm, double accuracy, ostream& os);
static void printStep(int step, double x, double dx,
                                    double f_of_x, ostream& os);
static void printWarning(int max_steps);

//      class RootFinder functions

void RootFinder::set_first_root_estimate(double x_guess)
{
        x0 = x_guess;
}

void RootFinder::set_step_estimate(double step_size)
{
        dx = step_size;
        x1 = x0 + dx;
}

void RootFinder::set_second_root_estimate(double x_guess)
{
        x1 = x_guess;
        dx = x1 - x0;
}

bool RootFinder::bracket_root(double f(double x))
{
     // This function will attempt to bracket the root of a function
     // starting with the RootFinder data values x0 and x1, by
     // geometrically expanding the interval [x0, x1] until the root
     // is bracketed, or the number of steps exceeds max_steps

        //      a convenient expansion factor
        const double expansionFactor = 1.6;

        if ( x0 == x1 ) {
                cerr << "\n RootFinder::bracket_root: sorry, x0 = x1!\n"
                        << " x0 = " << x0 << endl
                        << " x1 = " << x1 << endl;
                return false;
        }

        int step = 0;
        double f0 = f(x0);
        double f1 = f(x1);

        while ( ++step <= max_steps ) {
                if ( f0 * f1 <= 0 )
                        return true;         // success .. root is bracketed
                if ( abs(f0) < abs(f1) ) { // x0 probably closer to root
                        x0 -= expansionFactor * dx;
                        f0 = f(x0);
                } else {                          // x1 probably closer to the root
                        x1 += expansionFactor * dx;
                        f1 = f(x1);
                }
                dx = x1 - x0;             // also adjust dx
        }

        printWarning(max_steps);   // max_steps has been exceeded
        return false;
}

void RootFinder::set_accuracy(double epsilon)
{
        if ( epsilon > 1e-38 )
                accuracy = epsilon;
        else {
                cerr << " RootFinder: you have requested a ridiculous"
                        << " accuracy: " << accuracy << " !!!" << endl;
        }
}

void RootFinder::set_max_steps(int steps)
{
        if ( steps > 0 )
                max_steps = steps;
        else {
                cerr << " RootFinder: you have requested a ridiculous"
                        << " number of steps: " << steps << " !!!" << endl;
        }
}

//      class SimpleSearch functions

double SimpleSearch::find_root(double f(double x))
{
        double f0 = f(x0);
        double f_of_x = f0;

        int step = 0;
        if ( verbose ) {
                printHead("Simple Search with Step-Halving", accuracy, *os);
                printStep(step, x0, dx, f_of_x, *os);
        }

        while ( abs(dx) > accuracy && f_of_x != 0 && ++step <= max_steps ) {
                x0 += dx;                         //    take a step
                f_of_x = f(x0);
                if ( f0 * f_of_x < 0 ) {        //      jumped past root
                        x0 -= dx;                       //      backup
                        dx /= 2;                        //      halve the step size
                }
                if ( verbose )
                        printStep(step, x0, dx, f_of_x, *os);
        }

        if ( step > max_steps )
                printWarning(max_steps);
        steps = step;
        return x0;
}

//      class BisectionSearch functions

double BisectionSearch::find_root(double f(double x))
{
        double f0 = f(x0);
        double f1 = f(x1);

        if ( f0 * f1 > 0 ) {
                cerr << " BisectionSearch: sorry, root not bracketed!\n"
                        << " f(" << x0 << ") = " << f0 << endl
                        << " f(" << x1 << ") = " << f1 << endl
                        << " Trying to bracket the root using bracket_root ..."
                        << flush;
                double save_x0 = x0;
                double save_x1 = x1;
                if ( bracket_root(f) ) {
                        cerr << " Bracketing succeeded !\n"
                                << " x0 = " << x0 << " x1 = " << x1
                                << " continuing ..." << endl;
                        f0 = f(x0);
                        f1 = f(x1);
                } else {
                        cerr << " Sorry, bracketing attempt failed" << endl;
                        x0 = save_x0;
                        x1 = save_x1;
                        return abs(f0) < abs(f1) ? x0 : x1;
                }
        }
        if ( f0 == 0 )
                return x0;
        if ( f1 == 0 )
                return x1;

        double xHalf, fHalf = 0.5 * (f0 + f1);
        int step = 0;
        if ( verbose ) {
                printHead("Bisection Search", accuracy, *os);
                printStep(step, x0, x1 - x0, fHalf, *os);
        }
        do {                                         //         iteration loop
                if ( ++step > max_steps )
                        break;
                xHalf = 0.5 * (x0 + x1);        //      bisection point
                fHalf = f(xHalf);
                if ( f0 * fHalf > 0 ) {  //     x0 and xHalf on same side of root
                        x0 = xHalf;             //      replace x0 by xHalf
                        f0 = fHalf;
                } else {                 //     x1 and xHalf on same side of root
                        x1 = xHalf;             //      replace x1 by xHalf
                        f1 = fHalf;
                }
                if ( verbose )
                        printStep(step, x0, x1 - x0, fHalf, *os);
        } while ( abs(x1 - x0) > accuracy && fHalf != 0);

        if ( step > max_steps )
                printWarning(max_steps);
        steps = step;
        return xHalf;
}

//      class SecantSearch functions

double SecantSearch::find_root(double f(double x))
{
        double f0 = f(x0);
        double f1 = f(x1);

        if ( f0 == 0 )
                return x0;
        if ( f1 == 0 )
                return x1;

        int step = 0;
        if ( verbose ) {
                printHead("Secant Search", accuracy, *os);
                printStep(step, x0, x1 - x0, f1, *os);
        }
        do {
                if ( ++step > max_steps )
                        break;
                if ( f0 == f1 ) {
                        cerr << " Secant Search: f(x0) = f(x1), algorithm fails!\n"
                                << " f(" << x0 << ") = " << f0 << endl
                                << " f(" << x1 << ") = " << f1 << endl;
                        break;
                }
                dx *= - f1 / ( f1 - f0 );
                x0 = x1;
                f0 = f1;
                x1 += dx;
                f1 = f(x1);
                if ( verbose )
                        printStep(step, x0, dx, f1, *os);
        } while ( abs(dx) > accuracy && f1 != 0);

        if ( step > max_steps )
                printWarning(max_steps);
        steps = step;
        return x1;
}

//      class TangentSearch functions

double TangentSearch::find_root(double f(double x),
                                                   double f_prime(double x))
{
        double f0 = f(x0);
        double f_prime0 = f_prime(x0);

        if ( f0 == 0 )
                return x0;

        if ( f_prime0 != 0 )
                dx = - f0 / f_prime0;

        int step = 0;
        if ( verbose ) {
                printHead("Tangent Search", accuracy, *os);
                printStep(step, x0, dx, f0, *os);
        }
        do {
                if ( ++step > max_steps )
                        break;
                if ( f_prime0 == 0 ) {
                        cerr << " Tangent Search: f'(x0) = 0, algorithm fails!\n"
                                << " f(" << x0 << ") = " << f0 << endl
                                << " f'(" << x0 << ") = " << f_prime0 << endl;
                        break;
                }
                dx = - f0 / f_prime0;
                x0 += dx;
                f0 = f(x0);
                f_prime0 = f_prime(x0);
                if ( verbose )
                        printStep(step, x0, dx, f0, *os);
        } while ( abs(dx) > accuracy && f0 != 0);

        if ( step > max_steps )
                printWarning(max_steps);
        steps = step;
        return x0;
}

//      Utility functions

static void printHead(const char *algorithm, double accuracy, ostream& os )
{
        os << "\n ROOT FINDING using " << algorithm
           << "\n Requested accuracy = " << accuracy
           << "\n Step     Guess For Root          Step Size           Function Value"
           << "\n ----  --------------------  --------------------  --------------------"
           << endl;
}

static void printStep(int step, double x, double dx,
                                    double f_of_x, ostream& os)
{
        int w = os.width();
        int p = os.precision();
        ios::fmtflags f = os.flags();
        os.setf(ios::right, ios::adjustfield);
        os << " " << setw(4) << step << "  ";
        os.setf(ios::left, ios::adjustfield);
        os << setprecision(14)
           << setw(20) << x << "  "
           << setw(20) << dx << "  "
           << setw(20) << f_of_x
           << endl;
        os.width(w);
        os.precision(p);
        os.setf(f);
}

static void printWarning(int max_steps)
{
        cerr << " Warning: maximum number of steps "
                << max_steps << " exceeded!" << endl;
}

// FFT class

  FFT::FFT() {
    N = 0;
    f = 0;
    inverse = false;
  }

  void FFT::transform(
    std::vector< std::complex<double> > &data
  ) {
    N = data.size();
    f = &data;
    bit_reverse();
    for (int n = 1; n < N; n *= 2)
      Danielson_Lanczos(n);
    for (int i = 0; i < N; ++i)
      (*f)[i] /= std::sqrt(double(N));
  }

  void FFT::inverse_transform(
    std::vector< std::complex<double> > &data
  ) {
    inverse = true;
    transform(data);
    inverse = false;
  }

  std::vector<double> FFT::power(
    const std::vector< std::complex<double> > &data
  ) {
    std::vector<double> P(1 + N / 2);
    P[0] = std::norm(data[0]) / double(N);
    for (int i = 1; i < N / 2; i++)
      P[i] = (std::norm(data[i]) + std::norm(data[N-i])) / double(N);
    P[N/2] = std::norm(data[N/2]) / double(N);
    return P;
  }

  void FFT::bit_reverse() {
    int j = 1;
    for (int i = 1; i < N; ++i) {
      if (i < j) {
        std::complex<double> temp = (*f)[i-1];
        (*f)[i-1] = (*f)[j-1];
        (*f)[j-1] = temp;
      }
      int k = N / 2;
      while ( k < j ) {
        j -= k;
        k /= 2;
      }
      j += k;
    }
  }

  void FFT::Danielson_Lanczos(int n) {
    const double pi = 4 * atan(1.0);
    std::complex<double> W(0, pi / n);
    W = inverse ? std::exp(W) : std::exp(-W);
    std::complex<double> W_j(1, 0);
    for (int j = 0; j < n; ++j) {
      for (int i = j; i < N; i += 2 * n) {
        std::complex<double> temp = W_j * (*f)[n+i];
        (*f)[n+i] = (*f)[i] - temp;
        (*f)[i] += temp;
      }
      W_j *= W;
    }
  }

} /* namespace cpt */
