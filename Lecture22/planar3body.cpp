#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

#include "diffeq.hpp"
using namespace cpt;


class Planar3Body {

public : 
  Planar3Body( double im1, double im2, double im3 ) :
    m1( im1), m2(im2), m3(im3), G( 4 * pi * pi / (m1 + m2 + m3) )
  {
  }


  // represent point in extended phase space by 13-component vector
  // [ t, r1, v1, r2, v2, r3, v3 ] = [ t, x1, y1, vx1, vy1, ... ]

  Matrix<double,1> operator()(Matrix<double,1>& trv) const {

    double t = trv[0];
    double x1 = trv[1], y1 = trv[2],  vx1 = trv[3],  vy1 = trv[4],
      x2 = trv[5], y2 = trv[6],  vx2 = trv[7],  vy2 = trv[8],
      x3 = trv[9], y3 = trv[10], vx3 = trv[11], vy3 = trv[12];

    double r12 = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2)),
      r13 = sqrt((x1 - x3)*(x1 - x3) + (y1 - y3)*(y1 - y3)),
      r23 = sqrt((x2 - x3)*(x2 - x3) + (y2 - y3)*(y2 - y3));

    double ax1 = - G*m2*(x1 - x2)/(r12*r12*r12) - G*m3*(x1 - x3)/(r13*r13*r13),
      ay1 = - G*m2*(y1 - y2)/(r12*r12*r12) - G*m3*(y1 - y3)/(r13*r13*r13);

    double ax2 = - G*m1*(x2 - x1)/(r12*r12*r12) - G*m3*(x2 - x3)/(r23*r23*r23),
      ay2 = - G*m1*(y2 - y1)/(r12*r12*r12) - G*m3*(y2 - y3)/(r23*r23*r23);

    double ax3 = - G*m1*(x3 - x1)/(r13*r13*r13) - G*m2*(x3 - x2)/(r23*r23*r23),
      ay3 = - G*m1*(y3 - y1)/(r13*r13*r13) - G*m2*(y3 - y2)/(r23*r23*r23);

    Matrix<double,1> f(trv.size());
    f[0] = 1.0;
    f[1] = vx1,  f[2] = vy1,  f[3] = ax1,  f[4] = ay1;
    f[5] = vx2,  f[6] = vy2,  f[7] = ax2,  f[8] = ay2;
    f[9] = vx3, f[10] = vy3, f[11] = ax3, f[12] = ay3;

    return f;
  }



  static const double pi;


protected : 

  // set physical constants

  double m1;        // lightest body
  double m2;        // next heavier body
  double m3;        // heaviest body
  // use units such that G*(m1+m2+m3) = 4*pi**2
  double G;

};

const double Planar3Body::pi = 4 * atan(1.0);

int main() {

    double m1, m2, m3;
    cout << " Gravitational planar 3-body problem" << endl;
    cout << " Enter standard (0) or adaptive (1) RK4 scheme : " << endl;
    bool stdrk4 = true;
    cin >> stdrk4;
    cout << " Enter m1 m2 m3: ";
    cin >> m1 >> m2 >> m3;
    double x1, y1, vx1, vy1;
    cout << " Enter x1 y1 vx1 vy1: ";
    cin >> x1 >> y1 >> vx1 >> vy1;
    double x2, y2, vx2, vy2;
    cout << " Enter x2 y2 vx2 vy2: ";
    cin >> x2 >> y2 >> vx2 >> vy2;

    // compute position and velocity of m3 assuming center of mass at origin
    double x3 = - (m1 * x1 + m2 * x2) / m3, y3 = - (m1 * y1 + m2 * y2) / m3,
           vx3 = - (m1 * vx1 + m2 * vx2) / m3,
           vy3 = - (m1 * vy1 + m2 * vy2) / m3;
    cout << " x3 = " << x3 << " y3 = " << y3
         << " vx3 = " << vx3 << " vy3 = " << vy3 << endl;

    // compute net angular velocity
    double
    I = m1*(x1*x1 + y1*y1) + m2*(x2*x2 + y2*y2) + m3*(x3*x3 + y3*y3),
    L = m1*(x1*vy1 - y1*vx1) +  m2*(x2*vy2 - y2*vx2) +  m3*(x3*vy3 - y3*vx3),
    omega = L / I;
    cout << " Net angular velocity = " << omega << endl;

    // We can transform to the co-rotating frame v -> v - omega x r
    // this non-inertial transformation will change the trajectories!
    cout << " Enter 1 to transform to co-rotating frame, 0 otherwise: ";
    bool yes;
    cin >> yes;
    if (yes) {
        vx1 = vx1 + omega * y1 ; vy1 = vy1 - omega * x1;
        vx2 = vx2 + omega * y2 ; vy2 = vy2 - omega * x2;
        vx3 = vx3 + omega * y3 ; vy3 = vy3 - omega * x3;
        cout << " vx1 = " << vx1 << " vy1 = " << vy1 << endl;
        cout << " vx2 = " << vx2 << " vy2 = " << vy2 << endl;
        cout << " vx3 = " << vx3 << " vy3 = " << vy3 << endl;
    }

    double t_max, dt;
    cout << " Enter time step dt: ";
    cin >> dt;
    cout << " Enter total time to integrate: ";
    cin >> t_max;
    string file_name;
    cout << " Enter data file name: ";
    cin >> file_name;
    ofstream file(file_name.c_str());


    Planar3Body planar3Body( m1, m2, m3 );

    double t = 0;
    Matrix<double,1> trv(13);
    trv[0] = t;
    trv[1] = x1, trv[2]  = y1, trv[3]  = vx1, trv[4]  = vy1;
    trv[5] = x2, trv[6]  = y2, trv[7]  = vx2, trv[8]  = vy2;
    trv[9] = x3, trv[10] = y3, trv[11] = vx3, trv[12] = vy3;
    while (t <= t_max) {
        for (int i = 0; i < trv.size(); i++)
            file << " " << trv[i];
        file << "\n";
	if ( stdrk4 ) 
	  RK4_step(trv, dt, planar3Body );
	else 
	  RK4_adaptive_step(trv, dt, planar3Body, 1e-6);
        t = trv[0];
    }
    file.close();
}
