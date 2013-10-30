// pendul - Program to compute the motion of a simple pendulum
// using the Euler or Verlet method
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>

using namespace std;

int main() {

  //* Select the numerical method to use: Euler or Verlet
  cout << "Choose a numerical method : Euler (0) or Verlet (1): ";
  int method; cin >> method;
					   
  //* Set initial position and velocity of pendulum
  cout << "Enter initial angle (in degrees): "; 
  double theta0; cin >> theta0;
  const double pi = 3.141592654;
  double theta = theta0*pi/180;   // Convert angle to radians
  double omega = 0.0;             // Set the initial velocity

  //* Set the physical constants and other variables
  double g_over_L = 1.0;          // The constant g/L
  double time = 0.0;              // Initial time
  double time_old;                // Time of previous reversal
  int irev = 0;                   // Used to count number of reversals
  cout << "Enter time step: ";
  double tau; cin >> tau;

  //* Take one backward step to start Verlet
  double accel = -g_over_L*sin(theta);    // Gravitational acceleration
  double theta_old = theta - omega*tau + 0.5*tau*tau*accel;    

  //* Loop over desired number of steps with given time step
  //    and numerical method
  cout << "Enter number of time steps: ";
  int nStep;  cin >> nStep;
  std::vector<double> t_plot;
  std::vector<double> th_plot;
  std::vector<double> period;
  int iStep;
  for( iStep=0; iStep<nStep; iStep++ ) {  

    //* Record angle and time for plotting
    t_plot.push_back( time );            
    th_plot.push_back( theta*180/pi );   // Convert angle to degrees
    time += tau;
  
    //* Compute new position and velocity using 
    //    Euler or Verlet method
    accel = -g_over_L*sin(theta);    // Gravitational acceleration
    if( method == 0 ) {
      theta_old = theta;        // Save previous angle
      theta += tau*omega;       // Euler method
      omega += tau*accel; 
    } 
    else {  
      double theta_new = 2*theta - theta_old + tau*tau*accel;
      theta_old = theta;	    // Verlet method
      theta = theta_new;  
    }
  
    //* Test if the pendulum has passed through theta = 0;
    //    if yes, use time to estimate period
    if( theta*theta_old < 0 ) { // Test position for sign change
      cout << "Turning point at time t = " << time << endl;
      if( irev == 0 )          // If this is the first change,
        time_old = time;       // just record the time
      else {
        period.push_back( 2*(time - time_old) );
        time_old = time;
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
  ofstream plotOut("pend_plot.txt");
  for( i=0; i< t_plot.size(); i++ ) {
    plotOut << t_plot[i] << " "  << th_plot[i] << endl;
  }

  return 0;

}
