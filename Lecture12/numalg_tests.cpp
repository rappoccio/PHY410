#include "numalg_templates.hpp"

#include <math.h>
#include <iostream>

// Test function
double test1(double x) { return  exp(x) * log(x) - x * x;; }

// Test functor class
class test2 { 
public :
  test2() {}
  double operator() (double x) const { return  exp(x) * log(x) - x * x;; }
};

int main(void) 
{

  // Define a typedef for basic functions of the form
  //    double f(double);
  typedef double (*basicfunction)(double);
  
  // Define a functor class instance for "test2" :
  test2 mytest2;

  std::cout << "First root implementation : " << std::endl;
  double root1 = root_simple( test1, 4.0, -0.1, 1e-6, 1000, true ); 

  std::cout << std::endl << std::endl << std::endl 
	    << "Second root implementation : " << std::endl;
  double root2 = root_simple( mytest2, 4.0, -0.1, 1e-6, 1000, true ); 
  
}
