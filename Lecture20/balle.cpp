// balle - Program to compute the trajectory of a baseball
//         using the Euler method.
// Adapted from Garcia, Numerical Methods for Physics, 2nd Edition
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
using namespace std;

int main() {

  //* Set initial position and velocity of the baseball
  double y1, speed, theta;
  double r1[2], v1[2], r[2], v[2], accel[2]; 
  int euler = 0;
  cout << "Enter initial height (meters): "; cin >> y1;
  r1[0] = 0;  r1[1] = y1;     // Initial vector position
  cout << "Enter initial speed (m/s): "; cin >> speed; 
  cout << "Enter initial angle (degrees): "; cin >> theta;
  const double pi = 3.141592654; 
  v1[0] = speed*cos(theta*pi/180);   // Initial velocity (x)
  v1[1] = speed*sin(theta*pi/180);   // Initial velocity (y)
  r[0] = r1[0];  r[1] = r1[1];  // Set initial position and velocity
  v[0] = v1[0];  v[1] = v1[1];             

  //* Set physical parameters (mass, Cd, etc.)
  double Cd = 0.35;      // Drag coefficient (dimensionless)
  double area = 4.3e-3;  // Cross-sectional area of projectile (m^2)
  double grav = 9.81;    // Gravitational acceleration (m/s^2)
  double mass = 0.145;   // Mass of projectile (kg)
  double airFlag, rho;
  cout << "Air resistance? (Yes:1, No:0): "; cin >> airFlag;
  if( airFlag == 0 )
    rho = 0;      // No air resistance
  else
    rho = 1.2;    // Density of air (kg/m^3)
  double air_const = -0.5*Cd*rho*area/mass;  // Air resistance constant

  cout << "Use Euler (0), Euler-Cromer (1), or Midpoint(2) ? "; cin >> euler;

  //* Loop until ball hits ground or max steps completed
  double tau;
  cout << "Enter timestep, tau (sec): "; cin >> tau;
  int iStep, maxStep = 1000;   // Maximum number of steps
  std::vector<double> xplot (maxStep);  
  std::vector<double> yplot (maxStep);
  std::vector<double> xNoAir (maxStep);
  std::vector<double> yNoAir (maxStep);
  for( iStep=0; iStep<maxStep; iStep++ ) {

    //* Record position (computed and theoretical) for plotting
    xplot[iStep] = r[0];   // Record trajectory for plot
    yplot[iStep] = r[1];
    double t = (iStep)*tau;     // Current time
    xNoAir[iStep] = r1[0] + v1[0]*t;
    yNoAir[iStep] = r1[1] + v1[1]*t - 0.5*grav*t*t;
  
    //* Calculate the acceleration of the ball 
    double normV = sqrt( v[0]*v[0] + v[1]*v[1] );
    accel[0] = air_const*normV*v[0];   // Air resistance
    accel[1] = air_const*normV*v[1];   // Air resistance
    accel[1] -= grav;                  // Gravity
  
    //* Calculate the new position and velocity using Euler method
    if ( euler == 0 ) {       // Euler step
      r[0] += tau*v[0];                 
      r[1] += tau*v[1];                 
      v[0] += tau*accel[0];     
      v[1] += tau*accel[1];     
    } else if ( euler == 1 ) {// Euler-Cromer step
      v[0] += tau*accel[0];     
      v[1] += tau*accel[1];     
      r[0] += tau*v[0];                 
      r[1] += tau*v[1];                 
    } else {                  // Midpoint step
      double vx_last = v[0];
      double vy_last = v[1];
      v[0] += tau*accel[0];     
      v[1] += tau*accel[1];     
      r[0] += tau*0.5*(v[0] + vx_last);                 
      r[1] += tau*0.5*(v[1] + vy_last);
    }
  
    //* If ball reaches ground (y<0), break out of the loop
    if( r[1] < 0 )  {
      xplot[iStep] = r[0];  // Record last values computed
	  yplot[iStep] = r[1];
      break;                  // Break out of the for loop
    } 
  }

  //* Print maximum range and time of flight
  cout << "Maximum range is " << r[0] << " meters" << endl;
  cout << "Time of flight is " << iStep*tau << " seconds" << endl;

  //* Print out the plotting variables: 
  //    xplot, yplot, xNoAir, yNoAir
  ofstream plotOut("plot.txt"),
    noAirOut("noAirPlot.txt");
  int i;
  for( i=0; i<iStep; i++ ) {
    plotOut << xplot[i] << " ";
    plotOut << yplot[i] << endl;
  }
  for( i=0; i<iStep; i++ ) {
    noAirOut << xNoAir[i] << " ";
    noAirOut << yNoAir[i] << endl;
  }

  return 0;
}
