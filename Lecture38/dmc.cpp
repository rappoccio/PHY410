// Diffusion Monte Carlo program for the 3-D harmonic oscillator

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
using namespace std;

#include "random.hpp"


class DiffusionMC {
public:
  
  DiffusionMC(int nt, double idt) :
    N_T(nt), dt(idt)
  {
    
    N = N_T;                   // set N to target number specified by user
    for (int n = 0; n < N; n++) {
      ensureCapacity(n);
      for (int d = 0; d < DIM; d++)
	r[n][d] = rng.rand()- 0.5;
      alive[n] = true;
    }
    zeroAccumulators();
    E_T = 0;                   // initial guess for the ground state energy

    rng.set_seed(12345);
  }



  double V(double *r) {          // harmonic oscillator in DIM dimensions
    double rSqd = 0;
    for (int d = 0; d < DIM; d++)
      rSqd += r[d] * r[d];
    return 0.5 * rSqd;
  }

  void ensureCapacity(int index) {

    static int maxN = 0;       // remember the size of the array

    if (index < maxN)          // no need to expand array
      return;                // do nothing

    int oldMaxN = maxN;        // remember the old capacity
    if (maxN > 0)
      maxN *= 2;             // double capacity
    else
      maxN = 1;
    if (index > maxN - 1)      // if this is not sufficient
      maxN = index + 1;      // increase it so it is sufficient

    // allocate new storage
    double **rNew = new double* [maxN];
    bool *newAlive = new bool [maxN];
    for (int n = 0; n < maxN; n++) {
      rNew[n] = new double [DIM];
      if (n < oldMaxN) {     // copy old values into new arrays
	for (int d = 0; d < DIM; d++)
	  rNew[n][d] = r[n][d];
	newAlive[n] = alive[n];
	delete [] r[n];    // release old memory
      }
    }
    delete [] r;               // release old memory
    r = rNew;                  // point r to the new memory
    delete [] alive;
    alive = newAlive;
  }
  void zeroAccumulators() {
    ESum = ESqdSum = 0;
    for (int i = 0; i < NPSI; i++)
      psi[i] = 0;
  }


  void oneMonteCarloStep(int n) {

    // Diffusive step
    for (int d = 0; d < DIM; d++)
      r[n][d] += rng.rand_gauss() * sqrt(dt);

    // Branching step
    double q = exp(- dt * (V(r[n]) - E_T));
    int survivors = int(q);
    if (q - survivors > rng.rand())
      ++survivors;

    // append survivors-1 copies of the walker to the end of the array
    for (int i = 0; i < survivors - 1; i++) {
      ensureCapacity(N);
      for (int d = 0; d < DIM; d++)
	r[N][d] = r[n][d];
      alive[N] = true;
      ++N;
    }

    // if survivors is zero, then kill the walker
    if (survivors == 0)
      alive[n] = false;
  }

  void oneTimeStep() {

    // DMC step for each walker
    int N_0 = N;
    for (int n = 0; n < N_0; n++)
      oneMonteCarloStep(n);

    // remove all dead walkers from the arrays
    int newN = 0;
    for (int n = 0; n < N; n++)
      if (alive[n]) {
        if (n != newN) {
	  for (int d = 0; d < DIM; d++)
	    r[newN][d] = r[n][d];
	  alive[newN] = true;
        }
        ++newN;
      }
    N = newN;

    // adjust E_T

    if ( N > 0 )
      E_T += log(N_T / double(N)) / 10;

    // measure energy, wave function
    ESum += E_T;
    ESqdSum += E_T * E_T;
    for (int n = 0; n < N; n++) {
      double rSqd = 0;
      for (int d = 0; d < DIM; d++)
	rSqd = r[n][d] * r[n][d];
      int i = int(sqrt(rSqd) / rMax * NPSI);
      if (i < NPSI)
	psi[i] += 1;
    }

  }

  inline int getDIM() const { return DIM;}
  inline int getNPSI() const { return NPSI;}
  inline double const * const getPsi() const { return psi;}
  inline double getESum() const { return ESum;}
  inline double getESqdSum() const { return ESqdSum;}
  inline double getRMax() const { return rMax; }


  void print(std::ostream & out, int skip = -1 ) {
    if ( skip < 0 )
      skip = 1;

    double dr = getRMax() / getNPSI();
    for (int i = 0; i < getNPSI(); i+=skip) {
      double r = i * dr;
      out << r << '\t' << psi[i]  << std::endl;
    }    
       
  }
protected:
  static const int DIM = 3;             // dimensionality of space
  static const int NPSI = 100;   // number of bins for wave function

  // random walkers
  int N;                         // current number of walkers
  int N_T;                       // desired target number of walkers
  double **r;                    // x,y,z positions of walkers
  bool *alive;                   // is this walker alive?

  double dt;                     // Delta_t set by user
  double E_T;                    // target energy

  // observables
  double ESum;                   // accumulator for energy
  double ESqdSum;                // accumulator for variance
  double rMax = 4;               // max value of r to measure psi
  double psi[NPSI];              // wave function histogram

  cpt::Random rng;
};

int main() {


  double dt;
  int N_T, timeSteps;
    cout << " Diffusion Monte Carlo for the 3-D Harmonic Oscillator\n"
         << " -----------------------------------------------------\n";
    cout << " Enter desired target number of walkers: ";
    cin >> N_T;
    cout << " Enter time step dt: ";
    cin >> dt;
    cout << " Enter total number of time steps: ";
    cin >> timeSteps;

    DiffusionMC dmc(N_T,dt);


    // do 20% of timeSteps as thermalization steps
    int thermSteps = int(0.2 * timeSteps);
    for (int i = 0; i < thermSteps; i++)
        dmc.oneTimeStep();

    std::cout << "Initialization after thermalizing" << std::endl;
    dmc.print(std::cout,10);

    
    // production steps
    dmc.zeroAccumulators();
    for (int i = 0; i < timeSteps; i++) {
        dmc.oneTimeStep();

	if ( i % 100 == 0 && i > 0 ) {
	  std::cout << "i = " << i << ", Eavg = " << dmc.getESum() / i << std::endl;
	}
    }

    std::cout << "Final form" << std::endl;
    dmc.print(std::cout,10);

    std::cout << "ESum = " << dmc.getESum() << std::endl;
    std::cout << "ESqdSum = " << dmc.getESqdSum() << std::endl;
    
    // compute averages
    double EAve = dmc.getESum() / timeSteps;
    double EVar = dmc.getESqdSum() / timeSteps - EAve * EAve;
    cout << " <E> = " << EAve << " +/- " << sqrt(EVar / timeSteps) << endl;
    cout << " <E^2> - <E>^2 = " << EVar << endl;
    double psiNorm = 0, psiExactNorm = 0;
    double dr = dmc.getRMax() / dmc.getNPSI();
    double const * const psi = dmc.getPsi();
    for (int i = 0; i < dmc.getNPSI(); i++) {
        double r = i * dr;
        psiNorm += pow(r, dmc.getDIM()-1) * psi[i] * psi[i];
        psiExactNorm += pow(r, dmc.getDIM()-1) * exp(- r * r);
    }
    psiNorm = sqrt(psiNorm);
    psiExactNorm = sqrt(psiExactNorm);
    ofstream file("psi.data");
    for (int i = 0; i < dmc.getNPSI(); i++) {
        double r = i * dr;
        file << r << '\t' << pow(r, dmc.getDIM()-1) * psi[i] / psiNorm << '\t'
             << pow(r, dmc.getDIM()-1) * exp(- r * r / 2) / psiExactNorm << '\n';
    }
    file.close();
}
