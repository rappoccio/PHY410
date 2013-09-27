#ifndef CPT_DIFFEQ_HPP
#define CPT_DIFFEQ_HPP

#include "matrix.hpp"

namespace cpt {

// Fourth order Runge-Kutta integration routines

extern void RK4_step(                   // replaces x(t) by x(t + dt)
    Matrix<double,1>& x,                // solution vector
    double dt,                          // fixed time step
    Matrix<double,1> flow(              // derivative vector f
        Matrix<double,1>&) );           // argument of f(x)

extern double RK4_adaptive_step(        // returns adapted time step
    Matrix<double,1>& x,                // solution vector
    double dt,                          // initial time step
    Matrix<double,1> flow(
        Matrix<double,1>&),             // derivative vector
    double accuracy=1e-6);              // desired accuracy = default value

extern double RK4_integrate(            // returns adapted time step
    Matrix<double,1>& x,                // solution vector
    double dt,                          // initial time step
    Matrix<double,1> flow(
        Matrix<double,1>&),             // derivative vector
    double Delta_t,                     // finite step in time x[0]
    double accuracy=1e-6);              // desired accuracy = default value

} /* end namespace cpt */

#endif /* CPT_DIFFEQ_HPP */
