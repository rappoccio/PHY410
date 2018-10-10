// pendul - Program to compute the motion of a simple pendulum
// using the Euler or Verlet method
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include "nonlin.hpp"
#include "diffeq.hpp"

using namespace std;
using namespace cpt;

class PendulumFunction {
public : 
  PendulumFunction( double iL=g, double igamma=0, double iF_D = 0, double iomega_D = 0 ) :
    L(iL), omega_0( sqrt(g/L) ),gamma(igamma),F_D(iF_D),omega_D(iomega_D)
  {
  }

  static constexpr double g=9.0;          // acceleration of gravity

protected : 
  double L;                           // length in meters
 
  double omega_0;                     // natural frequency
  double gamma;                       // damping constant
  double F_D;                         // driving force amplitude
  double omega_D;                     // driving force frequency  
  
};

class PendulumAcceleration : public PendulumFunction {
public : 

  PendulumAcceleration( double iL=g, double igamma=0, double iF_D = 0, double iomega_D = 0 ) : 
    PendulumFunction( iL, igamma, iF_D, iomega_D ) {}

  
  // Acceleration of the pendulum
  double operator() (Matrix<double,1> const & p ) const {
    double t=p[0];
    double theta=p[1];
    double dtheta_dt=p[2];
    double a = - omega_0 * omega_0 * sin(theta);     // due to gravity
    if ( gamma > 0.0 )
      a += - gamma * dtheta_dt;                        // damping
    if ( F_D > 0.0 )
      a += F_D * cos(omega_D * t);                     // driving force
    return a;
  }
};

class PendulumDiffEq : public PendulumFunction {
public : 

  PendulumDiffEq( double iL=g, double igamma=0, double iF_D = 0, double iomega_D = 0 ) : 
    PendulumFunction( iL, igamma, iF_D, iomega_D ) {}

  Matrix<double,1> operator()( Matrix<double,1> const & p) const {
    double t = p[0], x = p[1], y = p[2];
    
    Matrix<double,1> out(3);
    out[0] = 1.0;
    out[1] = y;
    out[2] = -omega_0*omega_0*sin(x) - gamma*y + F_D*cos(omega_D*t);

    return out;
  }

};



int main() {

  //* Select the numerical method to use: Euler or Verlet
  cout << "Choose a numerical method : Euler (0), RK4 (1), or Adaptive RK4 (2): ";
  int method; cin >> method;
					   
  //* Set initial position and velocity of pendulum
  cout << "Enter initial angle (in degrees): "; 
  double theta0; cin >> theta0;
  const double pi = 3.141592654;
  double theta = theta0*pi/180;   // Convert angle to radians
  cout << "Enter initial angular frequency (degrees/s): ";
  double omega;  cin >> omega;           // Set the initial velocity
  omega = omega*pi/180.0;
  cout << "Enter damping coefficient : ";
  double gamma; cin >> gamma;
  cout << "Enter driving amplitude : ";
  double F_D; cin >> F_D;
  cout << "Enter driving frequency (degrees/s): ";
  double omega_D; cin >> omega_D;
  omega_D = omega_D*pi/180.0;

  //* Set the physical constants and other variables
  double L = PendulumDiffEq::g;  // The constant g/L
  double time = 0.0;              // Initial time
  double time_old;                // Time of previous reversal
  int irev = 0;                   // Used to count number of reversals
  cout << "Enter time step: ";
  double tau; cin >> tau;

  PendulumAcceleration pendulumAcc( L, gamma, F_D, omega_D);

  PendulumDiffEq pendulumDiffEq( L, gamma, F_D, omega_D);

  double T_D = 2*pi/omega_D;      // Period of driving force

  double T_p = 0.0;               // Period for Poincare section

  if ( omega_D > 0.0 )
    T_p = T_D;
  else 
    T_p = 2*pi/sqrt(PendulumDiffEq::g / L);


  double accuracy = 1e-6;

  //* Take one backward step to start Verlet

  Matrix<double,1> xv(3);
  xv[0] = time;
  xv[1] = theta;
  xv[2] = omega;
  double accel = pendulumAcc( xv );
  double theta_old = theta - omega*tau + 0.5*tau*tau*accel;    


  xv[0] = time;
  xv[1] = theta_old;
  xv[2] = omega;

  int iperiod = 0;
  int iperiod_old = -1;

  //* Loop over desired number of steps with given time step
  //    and numerical method
  cout << "Enter number of time steps: ";
  int nStep;  cin >> nStep;


  bool plotTrajectory = true; // plot Poincare and trajectory if true.
                              // plot only Poincare if false.

  // if ( nStep > 10000 ) {
  //   cout << ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
  //   cout << "::::: Warning : Only plotting Poincare section, not trajectory ::::" << endl;
  //   cout << ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
  // }
  // plotTrajectory = ( nStep < 10000 );

  std::vector<double> t_plot;
  std::vector<double> th_plot;
  std::vector<double> om_plot;
  std::vector<double> period;
  std::vector< std::pair<double,double> > poincare_map;
  double dt_min = tau;
  double dt_max = tau;
  int iStep;
  for( iStep=0; iStep<nStep; iStep++ ) {  


    theta_old = xv[1];

    //* Record angle and time for plotting
    if ( plotTrajectory ) {
      t_plot.push_back( xv[0] );

      double thplot = xv[1];
      double omplot = xv[2];

      while (thplot > pi ) thplot -= 2*pi;
      while (thplot <= -pi) thplot += 2*pi;

      th_plot.push_back( thplot*180/pi );   // Convert angle to degrees
      om_plot.push_back( omplot*180/pi );
    }
    
    iperiod = static_cast<int>( xv[0] / T_p );
    // Calculate a Poincare map of the dynamics
    if ( iperiod != iperiod_old ) {
      iperiod_old = iperiod;
      
      double thplot = xv[1];
      double omplot = xv[2];

      while (thplot > pi ) thplot -= 2*pi;
      while (thplot <= -pi) thplot += 2*pi;


      poincare_map.push_back( std::pair<double,double>(thplot*180/pi, omplot*180/pi) );
      std::cout << "Period reached. theta = " <<  thplot*180/pi << ", omega = " <<  omplot*180/pi << endl;
    }
    

    if ( method == 0 ){
      cpt::Euler_step( xv, tau, pendulumAcc );
    } 
    else if( method == 1 ) {
      cpt::RK4_step( xv, tau, pendulumDiffEq );
    } else if ( method == 2) {
      
      tau = cpt::RK4_adaptive_step(xv, tau, pendulumDiffEq, accuracy);
      dt_min = min(tau, dt_min);
      dt_max = max(tau, dt_max);
    }
  
    //* Test if the pendulum has passed through theta = 0;
    //    if yes, use time to estimate period
    if( xv[1]*theta_old < 0 ) { // Test position for sign change
      //cout << "Turning point at time t = " << time << endl;
      if( irev == 0 )          // If this is the first change,
        time_old = xv[0];       // just record the time
      else {
        period.push_back( 2*(xv[0] - time_old) );
        time_old = xv[0];
      }
      irev++;       // Increment the number of reversals
    }
  }
  int nPeriod = irev-1;    // Number of times period is measured

  //* Estimate period of oscillation, including error bar
  double AvePeriod = 0.0, ErrorBar = 0.0;
  int i;
  for( i=0; i<nPeriod; i++ ) {
    AvePeriod += period[i];
  }
  AvePeriod /= nPeriod;
  for( i=0; i<nPeriod; i++ ) {
    ErrorBar += (period[i] - AvePeriod)*(period[i] - AvePeriod);
  }
  ErrorBar = sqrt(ErrorBar/(nPeriod*(nPeriod-1)));
  cout << "Average period = " << AvePeriod << " +/- " << ErrorBar << endl;

  //* Print out the plotting variables: t_plot, th_plot
  if ( plotTrajectory ) {
    ofstream plotOut("pend_plot.txt");
    for( i=0; i< t_plot.size(); i++ ) {
      plotOut << t_plot[i] << " "  << th_plot[i] << " " << om_plot[i] << endl;
    }
    plotOut.close();
  }

  ofstream poincareOut("poincare_plot.txt");
  for( i=0; i< poincare_map.size(); i++ ) {
    poincareOut << poincare_map[i].first << " "  << poincare_map[i].second << endl;
  }
  poincareOut.close();



  return 0;

}
