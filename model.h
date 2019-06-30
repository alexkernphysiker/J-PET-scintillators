#ifndef ____MODEL___
#define ____MODEL___
#include <math_h/interpolate.h>
extern const MathTemplates::LinearInterpolation<> BC420_lambda;
extern const MathTemplates::LinearInterpolation<> polyester_absorp;
extern const MathTemplates::LinearInterpolation<> Si_Photo_QE;
extern const MathTemplates::LinearInterpolation<> tube_QE;
const int N_photons=2700;
const int virtual_experiments_count=10000;
#endif