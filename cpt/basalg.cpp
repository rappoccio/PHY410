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
