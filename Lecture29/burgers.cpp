#include "cptstd.hpp"
#include "matrix.hpp"
using namespace cpt;

#include <cstdio>


class Burgers {

public : 

  enum {SINE, STEP};
  enum Method { M_FTCS, M_LAX, M_LAX_WENDROFF, M_GODUNOV };

  Burgers( Method imethod = M_FTCS, int iN=500, double iL=1.0, double iCFL=1.0, int init_form=SINE ) :
    method(imethod),
    N(iN), L(iL), CFLRatio(iCFL), u(iN), uNew(iN), F(iN), 
    uPlus(iN), uMinus(iN)
  {

    initialWaveform = init_form;         // sine function, step, etc.
    
    nu = 1e-6;                   // kinematic viscosity


    h = L / N;

    for (int i = 0; i < N; i++) {
      double x = i * h;
      switch (initialWaveform) {
      case SINE:
	u[i] = sin(2 * pi * x) + 0.5 * sin(pi * x);
        break;
      case STEP:
	u[i] = 0;
	if (x > L / 4 && x < 3 * L / 4)
	  u[i] = 1;
	break;
      default:
	u[i] = 1;
	break;
      }
      if (std::abs(u[i]) > uMax)
	uMax = std::abs(u[i]);
    }

    dt = CFLRatio * h / uMax;
    t = 0;
    step = 0;
  }

  void take_step()
  {
    double t0 = t;
    switch ( method ) {
    case M_FTCS : 
      FTCS();
      break;
    case M_LAX : 
      Lax();
      break;
    case M_LAX_WENDROFF :
      LaxWendroff();
      break;
    case M_GODUNOV : default: 
      Godunov();
      break;
    };
    u = uNew;
    t += dt;
    ++step;
  }

  void FTCS()
  {
    for (int j = 0; j < N; j++) {
      int jNext = j < N - 1 ? j + 1 : 0;
      int jPrev = j > 0 ? j - 1 : N - 1;
      uNew[j] = u[j] * (1 - dt / (2 * h) * (u[jNext] - u[jPrev])) +
	nu * dt / h / h * (u[jNext] + u[jPrev] - 2 * u[j]);
    }
  }

  void Lax()
  {
    for (int j = 0; j < N; j++) {
      int jNext = j < N - 1 ? j + 1 : 0;
      int jPrev = j > 0 ? j - 1 : N - 1;
      uNew[j] = (u[jNext] + u[jPrev]) / 2
	- u[j] * dt / (2 * h) * (u[jNext] - u[jPrev])
	+ nu * dt / h / h * (u[jNext] + u[jPrev] - 2 * u[j]);
    }
  }

  void LaxWendroff()
  {
    for (int j = 0; j < N; j++)
      F[j] = u[j] * u[j] / 2;
    for (int j = 0; j < N; j++) {
      int jMinus1 = j > 0 ? j - 1 : N - 1;
      int jPlus1 = j < N - 1 ? j + 1 : 0;
      int jPlus2 = jPlus1 < N - 1 ? jPlus1 + 1 : 0;
      uNew[j] = (u[j] + u[jPlus1]) / 2 -
	(dt / 2 / h) * (F[jPlus1] - F[j]) +
	(nu * dt / (2 * h * h)) * (
				    (u[jPlus1] + u[jMinus1] - 2 * u[j]) / 2 +
				    (u[jPlus2] + u[j] - 2 * u[jPlus1]) / 2 );
    }
    for (int j = 0; j < N; j++)
      F[j] = uNew[j] * uNew[j] / 2;
    for (int j = 0; j < N; j++) {
      int jMinus1 = j > 0 ? j - 1 : N - 1;
      int jPlus1 = j < N - 1 ? j + 1 : 0;
      uNew[j] = u[j] - (dt / h) * (F[j] - F[jMinus1]) +
	(nu * dt / (h * h)) * (u[jPlus1] + u[jMinus1] - 2 * u[j]);
    }
  }

  void Godunov()
  {
    for (int j = 0; j < N; j++) {
      uPlus[j] = u[j] > 0 ? u[j] : 0;
      uMinus[j] = u[j] < 0 ? u[j] : 0;
    }
    for (int j = 0; j < N; j++) {
      int jNext = j < N - 1 ? j + 1 : 0;
      int jPrev = j > 0 ? j - 1 : N - 1;
      double f1 = uPlus[jPrev] * uPlus[jPrev] / 2;
      double f2 = uMinus[j] * uMinus[j] / 2;
      F[jPrev] = f1 > f2 ? f1 : f2;
      f1 = uPlus[j] * uPlus[j] / 2;
      f2 = uMinus[jNext] * uMinus[jNext] / 2;
      F[j] = f1 > f2 ? f1 : f2;
      uNew[j] = u[j]  + nu * dt / h / h * (u[jNext] + u[jPrev] - 2 * u[j]);
      uNew[j] -= (dt / h) * (F[j] - F[jPrev]);
    }
  }


  double get_L() const { return L; }
  double get_t() const { return t; }
  int    get_N() const { return N; }
  double get_u(int i) const { return u[i]; }
  double get_dt() const{ return dt; }
  void inc_t( double dt) { t += dt; }

protected : 
  static const double pi;

  Method method;                      // FTCS, Lax, Lax_Wendroff
  double L;                           // size of periodic region
  int N ;                             // number of grid points
  double h;                           // lattice spacing
  double t;                           // time
  double uMax;                        // maximum wave amplitude
  double dt;                          // time step
  double CFLRatio;                    // Courant-Friedrichs-Lewy ratio tau/h
  int initialWaveform ;               // sine function, step, etc.

  double nu;                          // kinematic viscosity
  Matrix<double,1> u, uNew;           // the solution and its update
  Matrix<double,1> F;                 // the flow
  Matrix<double,1> uPlus, uMinus;     // for Godunov scheme
  int step;                           // integration step number




};

const double Burgers::pi = 4 * atan(1.0);

int main()
{
  cout << " Finite Difference solution of the Burgers Equation\n"
       << " Choose a numerical method: 0) FTCS, 1) Lax, 2) Lax-Wendroff, 3) Godunov : ";
  int method, N;
  cin >> method;
  cout << " Enter number of grid cells: ";
  cin >> N;
  double L;
  cout << " Enter L: ";
  cin >> L;

  Burgers burgers( static_cast<Burgers::Method>(method), N, L );


#ifdef _WIN32
  ostringstream pyout;
  pyout << "env python.exe animator_for_cpp.py " << N;
  FILE *pypipe = _popen(pyout.str().c_str(), "w");
#else
  ostringstream pyout;
  pyout << "/usr/bin/env python animator_for_cpp.py " << N;
  FILE *pypipe = popen(pyout.str().c_str(), "w");
#endif


  // simple Mpl pipes animation
  cout << " Enter animation time: ";
  double t_max;
  cin >> t_max;
  double frame_rate = 30;
  double dt_frame = 1 / frame_rate;
  int steps_per_frame = max(1, int(dt_frame /   burgers.get_dt()));

  while (burgers.get_t() < t_max) {
    ostringstream os;
    for (int i = 0; i < N; i++)
      os << burgers.get_u(i) << ',';
    fprintf(pypipe, "%s\n", os.str().c_str());
    fflush(pypipe);
    time_t start_time = clock();
    for (int step = 0; step < steps_per_frame; step++) {
      burgers.take_step();
      burgers.inc_t( dt_frame );
    }
    while (true) {
      double secs = (clock() - start_time) / double(CLOCKS_PER_SEC);
      if (secs > dt_frame)
	break;
    }
  }
  fclose(pypipe);
}
