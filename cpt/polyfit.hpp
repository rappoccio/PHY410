#ifndef POLYFIT_HPP
#define POLYFIT_HPP

// Adapted from example code in 
// Garcia, "Numerical Methods for Physics", Second Edition.
// Download available at : 
// http://www.mathworks.com/matlabcentral/fileexchange/2270-numerical-methods-for-physics-2e/all_files
// This adapts the usage in two ways : 
// 1. Use the normal C++ indexing of 0...n-1
// 2. Use the Matrix classes from the CPT of UB's PHY 410 class, 
//    in turn based on Bjarne Stroustrup's Matrix classes here : 
//    http://www.stroustrup.com/Programming/Matrix/


#include "matrix.hpp"
#include "linalg.hpp"

namespace cpt {

  typedef cpt::Matrix<double,1> polyfit_vector;
  typedef cpt::Matrix<double,2> polyfit_matrix;

void polyfit( polyfit_vector x, polyfit_vector y, polyfit_vector sigma, int M, 
	      polyfit_vector& a_fit, polyfit_vector& sig_a, polyfit_vector& yy, double& chisqr) {
// Function to fit a polynomial to data
// Inputs 
//   x       Independent variable
//   y       Dependent variable
//   sigma   Estimate error in y
//   M       Number of parameters used to fit data
// Outputs
//   a_fit   Fit parameters; a(0) is intercept, a(1) is slope
//   sig_a   Estimated error in the parameters a()
//   yy      Curve fit to the data
//   chisqr  Chi squared statistic

  //* Form the vector b and design matrix A
  int i, j, k, N = x.dim1();
  polyfit_vector b(N);
  polyfit_matrix A(N,M);
  for( i=0; i<N; i++ ) {
   b(i) = y(i)/sigma(i);
   for( j=0; j<M; j++ )
     if ( sigma(i) > 0.0 )
       A(i,j) = pow(x(i),(double)(j))/sigma(i);
     else 
       A(i,j) = 0.0;
  }

  //* Compute the correlation matrix C 
  polyfit_matrix C(M,M), Cinv(M,M);
  for( i=0; i<M; i++ ) {   // (C inverse) = (A transpose) * A
    for( j=0; j<M; j++ ) {   
      Cinv(i,j) = 0.0;
      for( k=0; k<N; k++ )
        Cinv(i,j) += A(k,i)*A(k,j);
    }
  }
  
  inverse( Cinv, C );  // C = ( (C inverse) inverse)
 
  //* Compute the least squares polynomial coefficients a_fit
  for( k=0; k<M; k++ ) {
    a_fit(k) = 0.0;
    for( j=0; j<M; j++ )
      for( i=0; i<N; i++ ) {
       a_fit(k) += C(k,j) * A(i,j) * b(i);
      }
  }
  //* Compute the estimated error bars for the coefficients
  for( j=0; j<M; j++ )
    sig_a(j) = sqrt(C(j,j));                          
  //* Evaluate curve fit at each data point and compute Chi^2
  chisqr = 0.0;
  for( i=0; i<N; i++ ) {
    yy(i) = 0.0;      // yy is the curve fit
    for( j=0; j<M; j++ )
      yy(i) += a_fit(j) * pow( x(i), (double)(j) );  
    double delta = (y(i)-yy(i))/sigma(i);
    chisqr += delta*delta;  // Chi square
  }
}

}

#endif
