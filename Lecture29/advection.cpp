#include "cptstd.hpp"
#include "matrix.hpp"
using namespace cpt;

#include <cstdio>


class Advection {

public : 


  Advection( int im = 0, double iL = 5.0, int iN = 500, double ic = 1.0, bool usestep=false ) :
    method(im), L(iL), N(iN), c(ic), x(iN+1), u(iN+1), u_new(iN+1), step_wave(usestep)
  {    
    dx = L / double(N);         // grid spacing
    dt = dx/c;
    // initialize the spatial grid
    for (int i = 0; i <= N; i++)
      x[i] = i * dx;

    double x0 = L / 2.0;
    double sigma = 0.1 * L;
    

    // initialize the solution vector
    for(int i = 0; i <= N; i++)
      u[i] = f_0(x[i]);

    t = 0;
    step_number = 0;
  }

  double f_0(double x) {        // initial waveform

    double x_0 = L / 2;
    double sigma = 0.1 * L;   // width of initial waveform

    if (step_wave) {
      if (abs(x - x_0) < sigma)
	return 1;
      else
	return 0;
    } else {
      double k = pi / sigma;
      double gaussian = exp(- (x - x_0) * (x - x_0) / (2 * sigma * sigma));
      return cos(k * (x - x_0)) * gaussian;
    }
  }


  void FTCS() {
    for (int i = 0; i <= N; i++) {
      int i_plus_1 = i < N ? i + 1 : 0;
      int i_minus_1 = i > 0 ? i - 1 : N;
      u_new[i] = u[i] - c * dt / (2 * dx) * (u[i_plus_1] - u[i_minus_1]);
    }
  }

  void Lax() {
    for (int i = 0; i <= N; i++) {
      int i_plus_1 = i < N ? i + 1 : 0;
      int i_minus_1 = i > 0 ? i - 1 : N;
      u_new[i] = (u[i_plus_1] + u[i_minus_1]) / 2
	- c * dt / (2 * dx) * (u[i_plus_1] - u[i_minus_1]);
    }
  }

  void Lax_Wendroff() {
    double D = c * c * dt * dt / (2 * dx * dx);
    for (int i = 0; i <= N; i++) {
      int i_plus_1 = i < N ? i + 1 : 0;
      int i_minus_1 = i > 0 ? i - 1 : N;
      u_new[i] = u[i] - c * dt / (2 * dx) * (u[i_plus_1] - u[i_minus_1])
	+ D * (u[i_plus_1] + u[i_minus_1] - 2 * u[i]);
    }
  }

  void take_step() {
    switch (method) {
    case 0 : 
      FTCS();
      break;
    case 1 : 
      Lax();
      break;
    case 2 : 
    default : 
      Lax_Wendroff();
      break;
    };
    u = u_new;
    t += dt;
    ++step_number;
  }

  void save_u(int plot_number) {
    ostringstream os;
    os  << "u_" << plot_number << ".data";
    string file_name(os.str());
    ofstream file(file_name.c_str());
    if (file)
      cout << " writing " << file_name << endl;
    else
      cerr << " cannot open " << file_name << endl;
    for (int i = 0; i < N; i++)
      file << x[i] << '\t' << u[i] << '\n';
    file.close();
  }

  double get_L() const { return L; }
  double get_c() const { return c; }
  double get_t() const { return t; }
  int    get_N() const { return N; }
  double get_dx() const{ return dx; }
  double get_dt() const{ return dt; }
  double get_u(int i) const { return u[i]; }
  double get_x(int i) const { return x[i]; }

  void inc_t( double dt) { t += dt; }

protected : 

  static const double pi;

  double L;                       // system size
  int N;                          // number of grid intervals
  double c;                       // wave speed
  double dx;                      // grid spacing
  double t;                       // time
  double dt;                      // time step
  int method;                     // integration algorithm
  Matrix<double,1> x;             // grid points
  Matrix<double,1> u, u_new;      // wave amplitude
  int step_number;                // step number
  bool step_wave;                 // otherwise Gaussian modulate cosine


};

const double Advection::pi = 4 * atan(1.0);

int main()
{
  cout << " Finite Difference solution of the Advection Equation\n"
       << " Choose a numerical method: 1) FTCS, 2) Lax, 3) Lax-Wendroff : ";
  int method, N;
  cin >> method;
  cout << " Enter number of grid cells: ";
  cin >> N;
  double L,c;
  cout << " Enter L: ";
  cin >> L;
  cout << " Enter c: ";
  cin >> c;

  Advection advection( method, L, N, c );

#ifdef _WIN32
  ostringstream pyout;
  pyout << "env python.exe animator_for_cpp.py " << N;
  FILE *pypipe = _popen(pyout.str().c_str(), "w");
#else
  ostringstream pyout;
  pyout << "/usr/bin/env python animator_for_cpp.py " << N;
  FILE *pypipe = popen(pyout.str().c_str(), "w");
#endif

  int plots = 10;
  for (int plot = 0; plot < plots; plot++) {
    double delta_t = 0;
    while (delta_t <   advection.get_L() / (  advection.get_c() * plots)) {
      advection.take_step();
      delta_t +=   advection.get_t();
    }
  }
  
  // simple Mpl pipes animation
  cout << " Enter animation time: ";
  double t_max;
  cin >> t_max;
  advection = Advection(method,L,N,c);
  double frame_rate = 30;
  double dt_frame = 1 / frame_rate;
  int steps_per_frame = max(1, int(dt_frame /   advection.get_dt()));

  while (advection.get_t() < t_max) {
    ostringstream os;
    for (int i = 0; i < N; i++)
      os << advection.get_u(i) << ',';
    fprintf(pypipe, "%s\n", os.str().c_str());
    fflush(pypipe);
    time_t start_time = clock();
    for (int step = 0; step < steps_per_frame; step++) {
      advection.take_step();
      advection.inc_t( advection.get_dt() );
    }
    while (true) {
      double secs = (clock() - start_time) / double(CLOCKS_PER_SEC);
      if (secs > dt_frame)
	break;
    }
  }
  fclose(pypipe);
}
