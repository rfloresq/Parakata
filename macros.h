// *****************************************************************
// Purpose: Define some usefull macros
//
// 2003-2012: Roberto Flores-Moreno
// *****************************************************************

#ifndef SMACROS_H
#define SMACROS_H

#include <math.h>

#ifndef M_PI
  #define M_PI acos(-1.0)
#endif

#define ABS(a) ((a)>=0.0?(a):-(a))

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))

#define RAD(a) ((a)*M_PI/180.0)
#define DEG(a) ((a)*180.0/M_PI)

#define MOD(a,b) (((int)(a))-((int)(b))*(((int)(a))/((int)(b)))) 

#define AngstromToBohr(a) ((a)*1.8897)
#define BohrToAngstrom(a) ((a)/1.8897)

#define TOEV(a) ((a)*27.21138)

#endif  // SMACROS_H
