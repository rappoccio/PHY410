// Variational Monte Carlo for the harmonic oscillator

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include "random.hpp"



class QHO {
public :

  QHO(int Nin, double alphain, int MCStepsin) : 
    N(Nin), alpha(alphain), MCSteps(MCStepsin)
  {
    rng.set_seed(123456);

    x.resize(N);
    for (int i = 0; i < N; i++)
      x[i] = rng.rand() - 0.5;
    delta = 1;

    nPsiSqd = int((xMax - xMin) / dx);
    psiSqd.resize(nPsiSqd);
    nAccept = 0;

    zeroAccumulators();
  }

  void zeroAccumulators() {
    eSum = eSqdSum = 0;
    for (int i = 0; i < nPsiSqd; i++)
      psiSqd[i] = 0;
  }

  double p(double xTrial, double x) {

    // compute the ratio of rho(xTrial) / rho(x)
    return exp(- 2 * alpha * (xTrial*xTrial - x*x));
  }

  double eLocal(double x) {

    // compute the local energy
    return alpha + x * x * (0.5 - 2 * alpha * alpha);
  }

  void MetropolisStep() {

    // choose a walker at random
    int n = int(rng.rand() * N);

    // make a trial move
    double xTrial = x[n] + delta * rng.rand_gauss();

    // Metropolis test
    if (p(xTrial, x[n]) > rng.rand()) {
      x[n] = xTrial;
      ++nAccept;
    }

    // accumulate energy and wave function
    double e = eLocal(x[n]);
    eSum += e;
    eSqdSum += e * e;
    int i = int((x[n] - xMin) / dx);
    if (i >= 0 && i < nPsiSqd)
      psiSqd[i] += 1;
  }

  void oneMonteCarloStep() {

    // perform N Metropolis steps
    for (int i = 0; i < N; i++) {
      MetropolisStep();
    }
  }

  void adjustStep() {
    // perform 20% of MCSteps as thermalization steps
    // and adjust step size so acceptance ratio ~50%
    int thermSteps = int(0.2 * MCSteps);
    int adjustInterval = int(0.1 * thermSteps) + 1;

    std::cout << " Performing " << thermSteps << " thermalization steps ..."
	 << std::flush;
    for (int i = 0; i < thermSteps; i++) {
      oneMonteCarloStep();
      if ((i+1) % adjustInterval == 0) {
	delta *= nAccept / (0.5 * N * adjustInterval) ;
	nAccept = 0;
      }
    }
    std::cout << "\n Adjusted Gaussian step size = " << delta << std::endl;    
  }


  // production steps
  void doProductionSteps( ) {
    zeroAccumulators();
    nAccept = 0;
    std::cout << " Performing " << MCSteps << " production steps ..." << std::flush;
    for (int i = 0; i < MCSteps; i++)
      oneMonteCarloStep();
  }

  void print() {

    // compute and print energy
    double eAve = eSum / double(N) / MCSteps;
    double eVar = eSqdSum / double(N) / MCSteps - eAve * eAve;
    double error = sqrt(eVar) / sqrt(double(N) * MCSteps);
    std::cout << "\n <Energy> = " << eAve << " +/- " << error
	      << "\n Variance = " << eVar << std::endl;

    // write wave function squared in file
    std::ofstream file("psiSqd.data");
    double psiNorm = 0;
    for (int i = 0; i < nPsiSqd; i++)
      psiNorm += psiSqd[i] * dx;
    for (int i = 0; i < nPsiSqd; i++) {
      double x = xMin + i * dx;
      file << x << '\t' << psiSqd[i] / psiNorm << '\n';
    }
    file.close();
    std::cout << " Probability density written in file psiSqd.data" << std::endl;    
  }
  

protected :

  int N;                         // number of walkers
  std::vector<double> x;         // walker positions
  double delta;                  // step size
  double eSum;                   // accumulator to find energy
  double eSqdSum;                // accumulator to find fluctuations in E
  double xMin = -10;             // minimum x for histogramming psi^2(x)
  double xMax = +10;             // maximum x for histogramming psi^2(x)
  double dx = 0.1;               // psi^2(x) histogram step size
  std::vector<double> psiSqd;    // psi^2(x) histogram
  int nPsiSqd;                   // size of array
  double alpha;                  // trial function is exp(-alpha*x^2)
  int nAccept;                   // accumulator for number of accepted steps
  cpt::Random rng;               // random number generator
  int MCSteps;                   // number of MC steps
};
  
int main() {

  int N; double alpha; int MCSteps;

  std::cout << " Variational Monte Carlo for Harmonic Oscillator\n"
       << " -----------------------------------------------\n";
  std::cout << " Enter number of walkers: ";
  std::cin >> N;
  std::cout << " Enter parameter alpha: ";
  std::cin >> alpha;
  std::cout << " Enter number of Monte Carlo steps: ";
  std::cin >> MCSteps;

  
  QHO qho(N,alpha,MCSteps);
  
  qho.adjustStep();

  qho.doProductionSteps();
  
  qho.print();

}
