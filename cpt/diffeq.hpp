#ifndef CPT_DIFFEQ_HPP
#define CPT_DIFFEQ_HPP

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <iostream>

#include "matrix.hpp"

namespace cpt {
// Euler ODE solver


template< typename T>
void Euler_step(            // advance trajectory point one Euler time step
    Matrix<double,1>& p,    // trajectory point - input/output variable
    const double dt,        // time step for Euler - input
    T const & acceleration) // acceleration, looks like double acceleration(Matrix<double,1>&)
{
    double t = p[0], theta = p[1], dtheta_dt = p[2];
    double dtheta_dt_old = dtheta_dt;

    // apply Euler's algorithm
    dtheta_dt += acceleration(p) * dt;
    theta += dtheta_dt_old * dt;
    t += dt;

    p[0] = t, p[1] = theta, p[2] = dtheta_dt;
}

// Fourth order Runge-Kutta integration routines
template< typename T>
void RK4_step(                          // replaces x(t) by x(t + dt)
    Matrix<double,1>& x,                // solution vector
    double dt,                          // fixed time step
    T const & flow)                     // derivative, looks like Matrix<double,1> flow(Matrix<double,1>&)
{
    int n = x.size();
    Matrix<double,1> f(n), k1(n), k2(n), k3(n), k4(n), x_temp(n);
    f = flow(x);
    for (int i = 0; i < n; i++) {
        k1[i] = dt * f[i];
        x_temp[i] = x[i] + k1[i] / 2;
    }
    f = flow(x_temp);
    for (int i = 0; i < n; i++) {
        k2[i] = dt * f[i];
        x_temp[i] = x[i] + k2[i] / 2;
    }
    f = flow(x_temp);
    for (int i = 0; i < n; i++) {
        k3[i] = dt * f[i];
        x_temp[i] = x[i] + k3[i];
    }
    f = flow(x_temp);
    for (int i = 0; i < n; i++)
        k4[i] = dt * f[i];
    for (int i = 0; i < n; i++)
        x[i] += (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
}

template< typename T>
double RK4_adaptive_step(               // returns adapted time step
    Matrix<double,1>& x,                // solution vector
    double dt,                          // initial time step
    T const & flow,                     // derivative, looks like Matrix<double,1> flow(Matrix<double,1>&)
    double accuracy)                    // desired accuracy
{
    // from Numerical Recipes
    const double SAFETY = 0.9, PGROW = -0.2, PSHRINK = -0.25,
        ERRCON = 1.89E-4, TINY = 1.0E-30;
    int n = x.size();
    Matrix<double,1> scale = flow(x), x_half(n), Delta(n);
    for (int i = 0; i < n; i++)
        scale[i] = std::abs(x[i]) + std::abs(scale[i] * dt) + TINY;
    double error = 0;
    while (true) {
        dt /= 2;
        x_half = x;
        RK4_step(x_half, dt, flow);
        RK4_step(x_half, dt, flow);
        dt *= 2;
        Matrix<double,1> x_full = x;
        RK4_step(x_full, dt, flow);
        for (int i = 0; i < n; i++)
            Delta[i] = x_half[i] - x_full[i];
        error = 0;
        for (int i = 0; i < n; i++)
	  error = std::max(std::abs(Delta[i] / scale[i]), error);
        error /= accuracy;
        if (error <= 1)
            break;
        double dt_temp = SAFETY * dt * std::pow(error, PSHRINK);
        if (dt >= 0)
	  dt = std::max(dt_temp, 0.1 * dt);
        else
	  dt = std::min(dt_temp, 0.1 * dt);
        if (std::abs(dt) == 0.0) {
	    std::cerr << " step size underflow" << std::endl;
            exit(1);
        }
    }
    if (error > ERRCON)
        dt *= SAFETY * pow(error, PGROW);
    else
        dt *= 5.0;
    for (int i = 0; i < n; i++)
        x[i] = x_half[i] + Delta[i] / 15.0;
    return dt;
}

template< typename T>
double RK4_integrate(                   // returns adapted time step
    Matrix<double,1>& x,                // solution vector
    double dt,                          // initial time step
    T const & flow,                     // derivative, looks like Matrix<double,1> flow(Matrix<double,1>&)
    double Delta_t,                     // finite step in time x[0]
    double accuracy)                    // desired accuracy
{
    if (dt * Delta_t <= 0) {
      std::cerr << "dt * Delta_t <= 0" << std::endl;
        exit(1);
    }
    if (std::abs(dt) > std::abs(Delta_t))
        dt = Delta_t / 5.0;
    double t0 = x[0];
    while (std::abs(x[0] - t0) < std::abs(Delta_t))
        dt = RK4_adaptive_step(x, dt, flow, accuracy);
    RK4_step(x, t0 + Delta_t - x[0], flow);
    return dt;
}



} /* end namespace cpt */

#endif /* CPT_DIFFEQ_HPP */
