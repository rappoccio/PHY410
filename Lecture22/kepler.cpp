#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
using namespace std;

#include "diffeq.hpp"
using namespace cpt;

class Kepler {


public : 

  Kepler( bool istep_using_y = false ) :
    step_using_y(istep_using_y)
  {
  }

  void set_step_using_y( bool f = true ) { step_using_y = f; }

  //  Derivative vector for Newton's law of gravitation
  Matrix<double,1> operator()(Matrix<double,1>& trv) const {
    double t = trv[0], x = trv[1], y = trv[2], vx = trv[3], vy = trv[4];
    double r = sqrt(x*x + y*y);
    double ax = - G_m1_plus_m2 * x / (r*r*r);
    double ay = - G_m1_plus_m2 * y / (r*r*r);
    Matrix<double,1> flow(5);
    flow[0] = 1;
    flow[1] = vx;
    flow[2] = vy;
    flow[3] = ax;
    flow[4] = ay;
    if (step_using_y)           // change independent variable from t to y
      for (int i = 0; i < 5; i++)
	flow[i] /= vy;
    return flow;
  }

  void print_data(Matrix<double,1> values)
  {
    char names[][6] = { "t", "x", "y", "v_x", "v_y" },
      units[][6] = { "yr", "AU", "AU", "AU/yr", "AU/yr" };
      for (int i = 0; i < 5; i++)
        cout << " " << names[i] << "\t" << values[i]
             << " " << units[i] << "\n";
  }

  void write_data(ofstream& file, Matrix<double,1> values)
  {
    for (int i = 0; i < values.size(); i++)
      file << values[i] << "\t";
    file << "\n";
  }

  void integrate(Matrix<double,1>& trv, double dt, double t_max,
		 double accuracy=1e-6, bool adaptive=false)
  {
    string file_name("kepler_fixed.data"), f_or_a("fixed");
    if (adaptive) {
      file_name = "kepler_adapt.data";
      f_or_a  = "adaptive";
    }
    ofstream file(file_name.c_str());

    cout << "\n Initial conditions:\n";
    print_data(trv);
    cout << " Integrating with " << f_or_a << " step size ...";
    int step = 0;
    double dt_min = dt, dt_max = dt;

    while (true) {
      write_data(file, trv);
      double y_save = trv[2];
      step_using_y = false;
      if (adaptive) {
	dt = RK4_adaptive_step(trv, dt, *this, accuracy);
	dt_min = min(dt, dt_min);
	dt_max = max(dt, dt_max);
      } else
	RK4_step(trv, dt, *this);
      double t = trv[0], x = trv[1], y = trv[2], vx = trv[3], vy = trv[4];
      if (x > 0 && y * y_save < 0) {
	step_using_y = true;
	RK4_step(trv, -y, *this);
	write_data(file, trv);
	break;
      }
      if (t > 10 * t_max) {
	cerr << " t too big, quitting ..." << endl;
	break;
      }
      step += 1;
    }

    file.close();
    cout << "\n Number of " << f_or_a << " steps = " << step << endl;
    if (adaptive)
      cout << " Minimum dt = " << dt_min << "\n"
	   << " Maximum dt = " << dt_max << endl;
    print_data(trv);
    cout << " Trajectory data in " << file_name << endl;
  }

  static const double pi;
  static const double G_m1_plus_m2;


protected : 
  bool step_using_y;      //  to interpolate to y = 0


};

const double Kepler::pi = 4 * atan(1.0);
const double Kepler::G_m1_plus_m2 = 4 * pi * pi;


int main()
{
  cout << " Kepler orbit using fixed and then adaptive Runge-Kutta\n"
       << " Enter aphelion distance in AU: ";
  double r_aphelion, eccentricity, dt, accuracy;
  cin >> r_aphelion;
  cout << " Enter eccentricity: ";
  cin >> eccentricity;
  double a = r_aphelion / (1 + eccentricity);
  double T = pow(a, 1.5);
  double vy0 = sqrt(Kepler::G_m1_plus_m2 * (2 / r_aphelion - 1 / a));
  cout << " Semimajor axis a = " << a << " AU\n"
       << " Period T = " << T << " yr\n"
       << " v_y(0) = " << vy0 << " AU/yr\n"
       << " Enter step dt: ";
  cin >> dt;
  cout << " Enter desired accuracy for adaptive integration: ";
  cin >> accuracy;

  Kepler kepler;

  Matrix<double,1> trv(5);
  trv[0] = 0; trv[1] = r_aphelion; trv[2] = 0; trv[3] = 0; trv[4] = vy0;
  kepler.integrate(trv, dt, T, accuracy, false);
  trv[0] = 0; trv[1] = r_aphelion; trv[2] = 0; trv[3] = 0; trv[4] = vy0;
  kepler.integrate(trv, dt, T, accuracy, true);
}
