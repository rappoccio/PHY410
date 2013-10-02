#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>

#include <matrix.hpp>
#include <linalg.hpp>

// Solves a boundary value problem : 
//
//  d^2u / dx^2 = -pi^2/4 ( u+1)
//
// with Dirichlet boundary conditions u(0)=0 and u(1)=1.
// This is discretized as : 
// (2 u_i - u_(i+1) - u_(i-1)) / h^2 = pi^2/4 * (u_i + 1)
// 

int main(int argc, char ** argv)
{
  using namespace std;

  const int n = 2;               // Number of masses
  double h = 0.01;               // Discretization width
  double pi = atan(-1.0);        // Pi
  double k = h*h * pi*pi * 0.25; // constant in the tridiagonal matrix
  
  cpt::Matrix<double,1> u(n);
  cpt::Matrix<double,1> a(n);
  cpt::Matrix<double,1> b(n);
  cpt::Matrix<double,1> c(n);
  cpt::Matrix<double,1> r(n);

  // To check : fill the whole matrix
  cpt::Matrix<double,2> A(n,n);

  // Dirichlet boundary conditions
  double u0 = 0.0;
  double un = 1.0;

  // Set up our tridiagonal matrix
  for ( unsigned int i = 0; i < n; ++i ) {

    b[i] = 2.0 - k;

    if ( i != 0 ) 
      a[i] = -1.0;
    else 
      a[i] = 0.0;


    if ( i != n-1 )
      c[i] = -1.0;
    else
      c[i] = 0.0;

    if ( i == 0 ) 
      r[i] = u0 + k;
    else if ( i == n-1 )
      r[i] = un + k;
    else 
      r[i] = k;
  }

  for ( unsigned int i = 0; i < n; ++i ) {
    A[i][i] = 2-k;
    if ( i == 0 ) {
      A[i][i+1] = -1.0;
    }
    else if ( i == n-1 ) {
      A[i][i-1] = -1.0;
    }
    else {
      A[i][i+1] = -1.0;
      A[i][i-1] = -1.0;
    }
  }


  cout << "a = " << a << endl;
  cout << "b = " << b << endl;
  cout << "c = " << c << endl;

  cout << "r = " << r << endl;

  u = cpt::solve_tridiag( a, b, c, r, n);
  
  cout << "u = " << u << endl;

  cout << "Another way : " << endl;

  cout << "A = " << A << endl;
  
  cpt::Matrix<double,2> Ainv;
  double det = inverse(A, Ainv);

  cout << "Ainv = " << Ainv << endl;

 

  // cout << "uagain = " << uagain << endl;
  
  
  return 0;
}
