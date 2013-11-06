#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

#include "diffeq.hpp"
using namespace cpt;


// represent point in extended phase space by 5-component vector
// [ t, r, v ] = [ t, x, y, vx, vy ]

class Rcp3Body {
public : 

  Rcp3Body( double ia=0.1, bool istep_using_y=false ) :
    a(ia), step_using_y(istep_using_y) {}



  // equations in co-rotating frame
  Matrix<double,1> operator()(Matrix<double,1>& trv) const {

    double t = trv[0];
    double x = trv[1], y = trv[2], vx = trv[3], vy = trv[4];

    double d1 = pow((x - a)*(x - a) + y*y, 1.5),
      d2 = pow((x + 1 - a)*(x + 1 - a) + y*y, 1.5);

    double ax = - (1 - a) * (x - a) / d1 - a * (x + 1 -a) / d2 + x + 2 * vy,
      ay = - (1 - a) * y / d1 - a * y / d2 + y - 2 * vx;

    Matrix<double,1> flow(trv.size());
    flow[0] = 1;
    flow[1] = vx;
    flow[2] = vy;
    flow[3] = ax;
    flow[4] = ay;

    if (step_using_y) {     // change integration variable from t to y
      for (int i = 0; i < flow.size(); i++)
	flow[i] /= vy;
    }

    return flow;
  }

  double Jacobi(Matrix<double,1> trv) {     // Jacobi Integral

    double t = trv[0];
    double x = trv[1], y = trv[2], vx = trv[3], vy = trv[4];

    double r1 = sqrt((x - a)*(x - a) + y*y),
      r2 = sqrt((x + 1 - a)*(x + 1 - a) + y*y);

    return x*x + y*y + 2 * (1 - a) / r1 + 2 * a / r2 - vx*vx - vy*vy;
  }

  double f_x(double x) {      // effective x component of force on the x-axis
    return x - (1 - a) * (x - a) / abs(x - a) / (x - a) / (x - a)
      - a * (x + 1 - a) / abs(x + 1 - a) / (x + 1 - a) / (x + 1 - a);
  }

  // use bisection search to solve f(x) = 0 in interval [x_lower, x_upper]

  double zero( double x_lower, double x_upper,
	      double accuracy=1.0e-6, int max_steps=1000)
  {
    if (x_lower > x_upper) {
      cerr << " zero(f, x_lower, x_upper) not bracketed" << endl;
      exit(1);
    }

    double x_mid = (x_upper + x_lower) / 2;
    double dx = x_upper - x_lower;
    double f_lower = f_x(x_lower);
    int step = 0;
    while (abs(dx) > accuracy) {
      double f_mid = f_x(x_mid);
      if (f_mid == 0)
	dx = 0;
      else {
	if (f_lower * f_mid > 0) {
	  x_lower = x_mid;
	  f_lower = f_mid;
	} else
	  x_upper = x_mid;
	x_mid = (x_upper + x_lower) / 2;
	dx = x_upper - x_lower;
      }
      step += 1;
      if (step >= max_steps) {
	cerr << " zero(f) too many steps " << max_steps << endl;
	exit(1);
      }
    }
    return x_mid;
  }

  void set_step_using_y( bool f ) { step_using_y = f; }

protected : 
  // the restricted circular planar 3-body problem has one parameter
  double a;             // alpha = m2/(m1+m2) in the webnotes

  // switch to zero in on Poincare section point
  bool step_using_y;  // use y instead of t as independent variable


};

int main() {

    cout << " Restricted circular planar 3-body problem" << endl;

    // get parameters from user
    double a = 0.0;
    while (true) {
        cout << " Enter alpha = m2/(m1+m2) > 0 and < 0.5: ";
        cin >> a;
        if (a <= 0.0 || a > 0.5) {
            cout << " Bad alpha, please try again" << endl;
            continue;
        }
        break;
    }

    Rcp3Body rcp3Body( a, false );

    // find Lagrangian points for this alpha
    double eps = 1e-6;
    cout << " Lagrangian points:"
         << "\n L1 : x = " << rcp3Body.zero( a - 1 + eps, a - eps) << " y = " << 0.0
         << "\n L2 : x = " << rcp3Body.zero( -1.5, a - 1 - eps) << " y = " << 0.0
         << "\n L3 : x = " << rcp3Body.zero( a + eps, 1.5) << " y = " << 0.0
         << "\n L4 : x = " << a - 0.5 << " y = " << sqrt(3.0) / 2
         << "\n L5 : x = " << a - 0.5 << " y = " << sqrt(3.0) / 2 << endl;

    // get initial values
    Matrix<double,1> trv(5);
    double t = 0, x, y, vx, vy, C;
    cout << " Enter 0 to specify [x,y,vx,vy] or 1 to specify C and [x,y,vx]: ";
    int one;
    cin >> one;
    if (one == 1) {
        while (true) {
            cout << " Enter value of the Jacobi integral C: ";
            cin >> C;
            cout << " Enter x, y, vx: ";
            cin >> x, y, vx;
            double r1 = sqrt((x - a)*(x - a) + y*y),
                   r2 = sqrt((x + 1 - a)*(x + 1 -a ) + y*y);
            double vy_sqd = - C + x*x + y*y + 2 * (1 - a) / r1
                            + 2 * a / r2 - vx*vx;
            if (vy_sqd < 0.0)
                cout << " Sorry C too large, cannot solve for vy" << endl;
            else {
                vy = sqrt(vy_sqd);
                cout << " vy = " << vy << endl;
                break;
            }
        }
    } else {
        cout << " Enter x, y, vx, vy: ";
        cin >> x >> y >> vx >> vy;
        trv[0] = t; trv[1] = x; trv[2] = y; trv[3] = vx; trv[4] = vy;
        cout << " Jacobi Integral C = " << rcp3Body.Jacobi(trv) << endl;
    }
    double t_max;
    cout << " Enter maximum integration time t_max: ";
    cin >> t_max;

    ofstream trajectory_file("trajectory.data"), section_file("section.data");
    int crossing = 0;
    double dt = 0.01;
    t = 0;
    trv[0] = t; trv[1] = x; trv[2] = y; trv[3] = vx; trv[4] = vy;

    while (t < t_max) {

        // write trajectory point
        t = trv[0];
        x = trv[1]; y = trv[2]; vx = trv[3]; vy = trv[4];
        for (int i = 0; i < trv.size(); i++)
            trajectory_file << " " << trv[i];
        trajectory_file << "\n";
        double y_save = y;

        // use adaptive Runge-Kutta with default accuracy
        dt = RK4_adaptive_step(trv, dt, rcp3Body, 1e-6);

        // Poincare section at y = 0 and vy positive
        x = trv[1]; y = trv[2]; vx = trv[3]; vy = trv[4];
        if (y_save < 0 && y >= 0 && vy >= 0) {
	    rcp3Body.set_step_using_y ( true );
            double dy = -y;
            RK4_step(trv, dt, rcp3Body);
            t = trv[0];
            x = trv[1]; y = trv[2]; vx = trv[3]; vy = trv[4];
            for (int i = 0; i < trv.size(); i++)
                section_file << " " << trv[i];
            section_file << "\n";
            ++crossing;
            cout << " Crossing No." << crossing << " at t = " << t
                 << " C = " << rcp3Body.Jacobi(trv) << "\n";
            rcp3Body.set_step_using_y ( false );
        }
    }

    cout << " Trajectory in trajectory.data" << endl;
    trajectory_file.close();
    cout << " Poincare section points in section.data" << endl;
    section_file.close();
}
