#ifndef CPT_NONLIN_HPP
#define CPT_NONLIN_HPP

#include "matrix.hpp"

namespace cpt {

extern double find_maximum(     // returns f(x) at local maximum
    const double xa,            // guess for one point near desired maximum
    const double xb,            // guess for second point near maximum
    double (*)(                 // name of function to be maximized
       const double),           // argument
    double accuracy,            // desired accuracy in x
    double& xMax );             // value of x at local maximum

extern double find_minimum(     // returns f(x) at local minimum
    const double xa,            // guess for one point near desired minimum
    const double xb,            // guess for second point near minimum
    double (*)(                 // name of function to be minimized
      const double),            // argument
    double accuracy,            // desired accuracy in x
    double& xMin );             // value of x at local minimum

extern double golden(           // numerical recipes golden search routine
    const double,
    const double,
    const double,
    double (*)(const double),
    const double,
    double& );

extern void mnbrak(              // numerical recipes bracketing routine
    double&,
    double&,
    double&,
    double&,
    double&,
    double&,
    double (*)(const double) );

extern void minimize_BFGS(      // numerical recipes BFGS minimization
    Matrix<double,1>& p,        // input: initial guess chosen close to minimum
    const double gtol,          // input: desired accuracy
    int& iter,                  // output: number of iterations
    double& fret,               // output: value of function at minimum
    double (*func)(             // name of global function to minimize
        Matrix<double,1>&),
    void (*dfunc)(              // name of global gradient of function
        Matrix<double,1>&,
        Matrix<double,1>&));

} /* end namespace cpt */

#endif /* CPT_NONLIN_HPP */
