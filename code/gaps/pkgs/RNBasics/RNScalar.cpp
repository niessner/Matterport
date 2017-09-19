/* Source file for GAPS scalar class  */



/* Include files */

#include "RNBasics/RNBasics.h"



/* Public scalar variables */

#if (RN_MATH_PRECISION == RN_FLOAT_PRECISION)
    RNScalar RN_SMALL_EPSILON = 1.0E-4;
    RNScalar RN_EPSILON = 1.0E-3;
    RNScalar RN_BIG_EPSILON = 1.0E-2;
    RNScalar RN_INFINITY = 2.0E4;
#else
    RNScalar RN_SMALL_EPSILON = 1.0E-9;
    RNScalar RN_EPSILON = 1.0E-6;
    RNScalar RN_BIG_EPSILON = 1.0E-3;
    RNScalar RN_INFINITY = 1.0E6;
#endif



/* Public useful variables */

const RNScalar RN_E = 2.7182818284590452354;
const RNScalar RN_PI = 3.14159265358979323846;
const RNScalar RN_PI_OVER_FOUR = 0.25 * RN_PI;
const RNScalar RN_PI_OVER_TWO = 0.5 * RN_PI;
const RNScalar RN_THREE_PI_OVER_TWO = 1.5 * RN_PI;
const RNScalar RN_TWO_PI = 2.0 * RN_PI;
const RNScalar RN_ONE_OVER_PI = 1.0 / RN_PI;
const RNScalar RN_SQRT_TWO = 1.4142136;
const RNScalar RN_SQRT_THREE = 1.7320508;



/* Private variables */

static RNBoolean random_seeded = FALSE;



int RNInitScalar() 
{
    /* Return OK status */
    return TRUE;
}



void RNStopScalar()
{
}



void
RNSetInfinity(RNScalar infinity)
{
    // Set scalar range
    RN_INFINITY = infinity;
#   if (RN_MATH_PRECISION == RN_FLOAT_PRECISION)
       RN_EPSILON = 1.0E-6 * infinity;
       RN_SMALL_EPSILON = 0.1 * RN_EPSILON;
       RN_BIG_EPSILON = 10.0 * RN_EPSILON;
#   else 
       RN_EPSILON = 1.0E-10 * infinity;
       RN_SMALL_EPSILON = 0.1 * RN_EPSILON;
       RN_BIG_EPSILON = 10.0 * RN_EPSILON;
#endif    
}



void
RNSetEpsilon(RNScalar epsilon)
{
    // Set scalar range
    RN_EPSILON = epsilon;
#   if (RN_MATH_PRECISION == RN_FLOAT_PRECISION)
       RN_INFINITY = 1.0E7 * epsilon;
       RN_SMALL_EPSILON = 0.1 * RN_EPSILON;
       RN_BIG_EPSILON = 10.0 * RN_EPSILON;
#   else 
       RN_INFINITY = 1.0E10 * epsilon;
       RN_SMALL_EPSILON = 0.1 * RN_EPSILON;
       RN_BIG_EPSILON = 10.0 * RN_EPSILON;
#endif    
}



/* Random number functions */

void 
RNSeedRandomScalar(RNScalar seed)
{
#if (RN_OS == RN_WINDOWS)
  if (seed == 0.0) srand(GetTickCount());
  else srand((int) (1.0E6 * seed));
#else
  if (seed == 0.0) { 
      struct timeval timevalue;
      gettimeofday(&timevalue, NULL);
      srand48(timevalue.tv_usec);
  }
  else {
    srand48((long) (1.0E6 * seed));
  }
#endif
  random_seeded = TRUE;
}



RNScalar
RNRandomScalar(void)
{
    if (!random_seeded) RNSeedRandomScalar();
#   if (RN_OS == RN_WINDOWS)
    int r1 = rand();
    RNScalar r2 = ((RNScalar) rand()) / ((RNScalar) (RAND_MAX + 1));
    return (r1 + r2) / ((RNScalar) (RAND_MAX + 1));
#else
    return drand48();
#   endif
}



int 
RNCompareScalars(const void *value1, const void *value2)
{
  const RNScalar *scalar1 = (const RNScalar *) value1;
  const RNScalar *scalar2 = (const RNScalar *) value2;
  if (*scalar1 < *scalar2) return -1;
  else if (*scalar1 > *scalar2) return 1;
  else return 0;
}



int 
RNCompareDoubles(const void *value1, const void *value2)
{
  const double *double1 = (const double *) value1;
  const double *double2 = (const double *) value2;
  if (*double1 < *double2) return -1;
  else if (*double1 > *double2) return 1;
  else return 0;
}



int 
RNCompareFloats(const void *value1, const void *value2)
{
  const float *float1 = (const float *) value1;
  const float *float2 = (const float *) value2;
  if (*float1 < *float2) return -1;
  else if (*float1 > *float2) return 1;
  else return 0;
}



int 
RNCompareInts(const void *value1, const void *value2)
{
  const int *int1 = (const int *) value1;
  const int *int2 = (const int *) value2;
  if (*int1 < *int2) return -1;
  else if (*int1 > *int2) return 1;
  else return 0;
}



