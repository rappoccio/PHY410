/* slope_swig.i */
%module slope
%{
/* Put header files here or function declarations like below */
extern float getSlope( float x0, float y0, float x1, float y1);
%}

extern float getSlope( float x0, float y0, float x1, float y1);
