#include <vector>
#include <iostream>
#include <math.h>
#include "trapezoid.hpp"
#include "simpson.hpp"

using namespace std;


int main()
{
  unsigned n1=1,n2=0;
  cout << "Trapezoid quadrature of sin(x) [0, pi/2]" << endl;
  while ( n1 % 2 != 0 ) {
    cout << "Enter number of intervals desired for trapezoidal rule (must be even)" << endl;  
    cin >> n1;
  }
  double a = 0;
  double b = 2 * atan(1.0);
  double ans1 = trapezoid(sin, a, b, n1);
  cout << "Trapezoidal rule : " << ans1 << endl;


  double ans2 = adaptive_trapezoid( sin, a, b, 0.0001);
  cout << "Adaptive trapezoidal rule : " << ans2 << endl;
}
