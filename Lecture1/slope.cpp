#include <iostream>
#include <math.h>

using std::cout;
using std::endl; 

// "Dumb" slope calculator : No error handling
float getSlope0( float x0, float y0, float x1, float y1) {
  float slope = (y1-y0) / (x1-x0);
  return slope;
}

// Better slope calculator : Error handling when
// denominator is zero
float getSlope( float x0, float y0, float x1, float y1) {
  float num = y1 - y0;
  float den = x1 - x0;
  float slope = 0.0;
  if ( fabs(den) > 0.00001 ) {
    slope = num / den;
  } else {
    cout << "----> invalid input, setting to -9999" << endl;
    slope = -9999;
  }
  return slope;
}


int main( void )
{
  // Unit tests!

  cout << "Unit test 1" << endl;
  float x0 = 0.0;
  float y0 = 0.0;
  float x1 = 1.0;
  float y1 = 1.0;
  float slope = getSlope( x0,y0,x1,y1 );
  cout << slope << endl;

  cout << "Unit test 2" << endl;
  x0 = 0.0;
  y0 = 0.0;
  x1 = 0.0;
  y1 = 1.0;
  slope = getSlope( x0,y0,x1,y1 );
  cout << slope << endl;

  return 0;
}
