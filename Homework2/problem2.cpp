#include<iostream>

int main(){
   int i1 = 2, i2 = 5, i3 = -3;
   double d1 = 2.0, d2 = 5.0, d3 = -0.5;
   std::cout << "(a)" << i1 + (i2 * i3) << std::endl
             << "(b)" << i1 * (i2 + i3) << std::endl
             << "(c)" << i1 / (i2 + i3) << std::endl
             << "(d)" << i1 / i2 + i3 << std::endl
             << "(e)" << 3 + 4 + 5 / 3 << std::endl
             << "(f)" << (3 + 4 +5) / 3 << std::endl
             << "(g)" << d1 + (d2 * d3) << std::endl
             << "(h)" << d1 + d2 * d3 << std::endl
             << "(i)" << d1 / d2 - d3 << std::endl
             << "(j)" << d1 / (d2 - d3) << std::endl
             << "(k)" << d1 + d2 + d3 / 3 << std::endl
             << "(l)" << (d1 + d2 + d3) / 3 << std::endl
             << "(m)" << d1 + d2 + (d3 / 3) << std::endl
             << "(n)" << 3 * (d1 + d2) * (d1 - d3) << std::endl;
   return 0;
}
