#include <iostream>
#include <cmath>
int main(){
   const double tol=1e-6;
   double x1,y1,x2,y2;
   double midpoint_x,midpoint_y,slope,y_intercept;
   std::cout << "Please input the x and y values for the first point:\n"
             << "Format (x1 y1 x2 y2)" << std::endl;
   std::cin  >> x1 >> y1 >> x2 >> y2;
   if (std::abs(x1-x2) <= tol){
      std::cout << "ERROR DIVISION BY ZERO" << std::endl
                << "Slope is infinite" << std::endl;
      return -1;
   } else {
      std::cout << "--------------------Answers--------------------------" 
                << std::endl;
      slope = (y1-y2) / (x1-x2);
      midpoint_x = (x1+x2)/2;
      midpoint_y = (y1+y2)/2;
      y_intercept = -(slope*x1)+y1;
      std::cout << "Slope = " << slope << std::endl
                << "Intercept = " << y_intercept << std::endl
                << "Midpoint = (" << midpoint_x << "," << midpoint_y << ")"
                << std::endl;
      if (y_intercept < 0.0){
         std::cout << "Equation of the line:\n" << "    y = " << slope << "x - "
                   << std::abs(y_intercept) << std::endl;
      } else {
         std::cout << "Equation of the line:\n" << "    y = " << slope << "x + "
                   << y_intercept << std::endl;
      }
      return 0;
   }
}
