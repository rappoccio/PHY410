#include "nonlin.hpp"
#include "diffeq.hpp"
#include "cptstd.hpp"
#include "basalg.hpp"
using namespace cpt;

#include <fstream>


class BVP {

public : 
  
  enum extended_space { _x, _u, _dudx, dimension };

  BVP(int iN = 100, 
      double ixlb=0., double ixrb = 1.0,
      double iulb=0., double iurb = 1.0,
      double idelta = 1.0
      ) : 
    N(iN),
    x_lb(ixlb), x_rb(ixrb), u_lb(iulb), u_rb(iurb),
    delta ( idelta ),
    u(N+1), u_old(N+1),
    p(dimension) 
  {
    dx = (x_rb - x_lb) / N;      // Runge-Kutta step size
    h = (x_rb - x_lb) / N;       // step size 

    set_initial_values();
  }


  void set_initial_values()
  {
    p[_x] = x_lb;
    p[_u] = u_lb;
    p[_dudx] = delta;
  }


  void update_values( BVP const & right )
  {
    p = right.p;
    delta = right.delta;
    u = right.u;
    u_old = right.u_old;
  }



  double           get_delta() const { return delta; }
  double           get_dx()    const { return dx; }
  double           get_h()     const { return h; }
  Matrix<double,1> get_p()     const { return p; }
  int              get_N()     const { return N; }
  double           get_x_lb()  const { return x_lb; }
  double           get_x_rb()  const { return x_rb; }
  double           get_u_lb()  const { return u_lb; }
  double           get_u_rb()  const { return u_rb; }
  double           get_u(int i)const { return u[i]; }

protected : 

  static const double pi;

  int N;                          // number of integration steps
  double x_lb;                    // left boundary point
  double x_rb;                    // right boundary point
  double u_lb;                    // boundary condition at x = 0
  double u_rb;                    // boundary condition at x = 1
  double delta;                   // guess for u'(0)
  double dx;                      // Runge-Kutta step size
  double h;                       // step size
  Matrix<double,1> u, u_old;      // solution at lattice points
                                  // and solution before Jacobi step
  Matrix<double,1> p;             // point in extended space
};

const double BVP::pi = 4 * atan(1.0);




// derivative vector for Runge-Kutta
class BVPStep : public BVP {
public : 
  BVPStep( int iN = 100, 
	   double ixlb=0., double ixrb = 1.0,
	   double iulb=0., double iurb = 1.0,
	   double idelta=1.0) : 
    BVP( iN, ixlb, ixrb, iulb, iurb, idelta ) {}


  Matrix<double,1> operator()(Matrix<double,1>& p) const
  {
    double x = p[_x], u = p[_u], dudx = p[_dudx];
    Matrix<double,1> f(dimension);
    f[_x] = 1;
    f[_u] = dudx;
    f[_dudx] = - pi * pi * (u + 1) / 4;
    return f;
  }

  void step() {
    RK4_step(p, dx, *this);      // take one Runge-Kutta step  
  }

  friend class BVPShoot;
};


// function whose root is to be found
class BVPShoot : public BVP {
public : 
  BVPShoot( BVPStep const & istepper) : 
  stepper(istepper),
  BVP( istepper.N, istepper.x_lb, istepper.x_rb, istepper.u_lb, istepper.u_rb, istepper.delta ) {}

  double operator()(double x)                  
  {
    delta = x;                      // set global delta for root finding
    set_initial_values();           // in global extended point p
    for (int i = 0; i < N; i++)     // integrate from x_lb to x_rb
      RK4_step(p, dx, stepper);     // using 4th order Runge-Kutta
    double u = p[_u];               // solution at right boundary
    update_values( stepper) ;       // update "this" values with what's in stepper
    return u - u_rb;                // mismatch in right boundary value
  }

  double solve() {
    cpt::SecantSearchT<BVPShoot> ss;                // Secant search object
    delta = ss.find_root(*this, delta, delta + dx);
    return delta;
  }

protected : 
  BVPStep stepper;

};



// function whose root is to be found
class BVPRelax : public BVP {
public : 
  BVPRelax( int iN = 100, 
      double ixlb=0., double ixrb = 1.0,
      double iulb=0., double iurb = 1.0,
      double idelta = 1.0
      ) : 
    BVP( iN, ixlb, ixrb, iulb, iurb, idelta ) {}

  double d2u_dx2(double u)            // RHS of second order differential equation
  {
    return pi * pi / 4 * (u + 1);
  }

  double u_exact(double x)            // exact solution
  {
    return cos(pi * x / 2) + 2 * sin(pi * x / 2) - 1;
  }


  void initial_guess()                // linear fit to boundary conditions
  {
    for (int i = 0; i <= N; i++) {
      double x = x_lb + i * h;
      u[i] = u_lb + (x - x_lb) / (x_rb - x_lb) * (u_rb - u_lb);
    }
    u_old = u;
  }

  void Jacobi_step()                  // Jacobi algorithm
  {
    u_old = u;                      // save solution before relaxing
    for (int i = 1; i < N; i++)
      u[i] = 0.5 * (u_old[i+1] + u_old[i-1] + h * h * d2u_dx2(u_old[i]));
  }

  double change()                     // compute maximum change
  {
    double max_diff = 0;
    for (int i = 1; i < N; i++) {
      double diff = abs(u[i] - u_old[i]);
      if (diff > max_diff)
	max_diff = diff;
    }
    return max_diff;
  }

  double error() {                    // error from exact solution
    double max_diff = 0;
    for (int i = 1; i < N; i++) {
      double x = x_lb + i * h;
      double diff = abs(u[i] - u_exact(x));
      if (diff > max_diff)
	max_diff = diff;
    }
    return max_diff;
  }



  void print(int iteration)
  {
    std::cout.precision(14);
    std::cout << " "
	      << std::setw(6) << iteration << "      "
	      << std::setw(20) << change() << " "
	      << std::setw(20) << error() << std::endl;
  }

};
