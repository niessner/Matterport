/* Include file for GAPS scalar class */



/* Initialization functions */

int RNInitScalar();
void RNStopScalar();



/* Scalar class definition */

#if (RN_MATH_PRECISION == RN_FLOAT_PRECISION)
    typedef float RNScalar;
#else
    typedef double RNScalar;
#endif



/* Scalar subclass definitions */

typedef RNScalar RNScalar;
typedef RNScalar RNAngle;
typedef RNScalar RNCoord;
typedef RNScalar RNMagnitude;
typedef RNMagnitude RNLength;
typedef RNMagnitude RNArea;
typedef RNMagnitude RNVolume;



/* Useful constant definitions */

extern const RNScalar RN_E;
extern const RNScalar RN_PI;
extern const RNScalar RN_PI_OVER_FOUR;
extern const RNScalar RN_PI_OVER_TWO;
extern const RNScalar RN_THREE_PI_OVER_TWO;
extern const RNScalar RN_TWO_PI;
extern const RNScalar RN_ONE_OVER_PI;
extern const RNScalar RN_SQRT_TWO;
extern const RNScalar RN_SQRT_THREE;
#define RN_DEG2RAD(a) ((a) * (RN_PI / 180.0))
#define RN_RAD2DEG(a) ((a) * (180.0 / RN_PI))



/* X/Y/Z dimension constant definitions */

typedef int RNDirection;
#define RN_LO                    0
#define RN_HI                    1
#define RN_NUM_DIRECTIONS        2

typedef int RNDimension;
#define RN_X                     0
#define RN_Y                     1
#define RN_Z                     2

typedef RNDimension RNAxis;
#define RN_XAXIS                 RN_X  
#define RN_YAXIS                 RN_Y
#define RN_ZAXIS                 RN_Z

typedef RNDimension RNAxisPlane;
#define RN_YZPLANE               RN_X  
#define RN_XZPLANE               RN_Y  
#define RN_XYPLANE               RN_Z  

typedef int RNSextant;
#define RN_NX_SEXTANT            0     
#define RN_PX_SEXTANT            1        
#define RN_NY_SEXTANT            2        
#define RN_PY_SEXTANT            3        
#define RN_NZ_SEXTANT            4        
#define RN_PZ_SEXTANT            5        

typedef int RNSide;
#define RN_LX_SIDE               RN_NX_SEXTANT
#define RN_HX_SIDE               RN_PX_SEXTANT
#define RN_LY_SIDE               RN_NY_SEXTANT
#define RN_HY_SIDE               RN_PY_SEXTANT    
#define RN_LZ_SIDE               RN_NZ_SEXTANT
#define RN_HZ_SIDE               RN_PZ_SEXTANT

// One bits per XY: 0 = LO (N), 1 = HI (P)
typedef int RNQuadrant;
#define RN_NN_QUADRANT           0        
#define RN_NP_QUADRANT           1        
#define RN_PN_QUADRANT           2        
#define RN_PP_QUADRANT           3        
#define RN_NUM_QUADRANTS         4

// One bits per XYZ: 0 = LO (N), 1 = HI (P)
typedef int RNOctant;
#define RN_NNN_OCTANT            0        
#define RN_NNP_OCTANT            1        
#define RN_NPN_OCTANT            2        
#define RN_NPP_OCTANT            3        
#define RN_PNN_OCTANT            4        
#define RN_PNP_OCTANT            5        
#define RN_PPN_OCTANT            6        
#define RN_PPP_OCTANT            7        
#define RN_NUM_OCTANTS           8

// Two bits per XYZ: 00 = Inside (Z), 01 = Below (N), 10 = Above (P)
typedef int RNBoxtant;        
#define RN_ZZZ_BOXTANT           0x00
#define RN_ZZN_BOXTANT           0x01
#define RN_ZZP_BOXTANT           0x02
#define RN_ZNZ_BOXTANT           0x04
#define RN_ZNN_BOXTANT           0x05               
#define RN_ZNP_BOXTANT           0x06
#define RN_ZPZ_BOXTANT           0x08               
#define RN_ZPN_BOXTANT           0x09               
#define RN_ZPP_BOXTANT           0x0A               
#define RN_NZZ_BOXTANT           0x10
#define RN_NZN_BOXTANT           0x11
#define RN_NZP_BOXTANT           0x12
#define RN_NNZ_BOXTANT           0x14
#define RN_NNN_BOXTANT           0x15               
#define RN_NNP_BOXTANT           0x16
#define RN_NPZ_BOXTANT           0x18               
#define RN_NPN_BOXTANT           0x19               
#define RN_NPP_BOXTANT           0x1A               
#define RN_PZZ_BOXTANT           0x20
#define RN_PZN_BOXTANT           0x21
#define RN_PZP_BOXTANT           0x22
#define RN_PNZ_BOXTANT           0x24
#define RN_PNN_BOXTANT           0x25               
#define RN_PNP_BOXTANT           0x26
#define RN_PPZ_BOXTANT           0x28               
#define RN_PPN_BOXTANT           0x29               
#define RN_PPP_BOXTANT           0x2A     



/* Relationship definitions */

#define RN_UNKNOWN               -1
#define RN_INSIDE                0
#define RN_OUTSIDE               1
#define RN_BELOW                 1
#define RN_ABOVE                 2          
#define RN_CROSSING              3



/* Scalar constant definitions */

extern RNScalar RN_INFINITY;
extern RNScalar RN_EPSILON;
extern RNScalar RN_SMALL_EPSILON;
extern RNScalar RN_BIG_EPSILON;



/* Useful constants */

#define RN_SPEED_OF_SOUND 13543.3       // (inches/second)



/* Range functions */

extern void RNSetInfinity(RNScalar infinity);
extern void RNSetEpsilon(RNScalar epsilon);



/* Random number generator */

extern void RNSeedRandomScalar(RNScalar seed = 0.0);
extern RNScalar RNRandomScalar(void);



/* Useful sorting functions (e.g., for qsort) */

extern int RNCompareScalars(const void *value1, const void *value2);
extern int RNCompareDoubles(const void *value1, const void *value2);
extern int RNCompareFloats(const void *value1, const void *value2);
extern int RNCompareInts(const void *value1, const void *value2);



/* Math functions */

#if ((RN_OS == RN_IRIX) && (RN_MATH_PRECISION == RN_FLOAT_PRECISION))
#define cos cosf
#define sin sinf
#define tan tanf
#define acos acosf
#define asin asinf
#define atan atanf
#define atan2 atan2f
#define sqrt sqrtf
#define exp expf
#define log logf
#define log10 log10f
#define pow powf
#endif



/* Tolerance functions */

inline RNScalar
RNSign(RNScalar scalar)
{
    // Return +1, -1, or 0
    if (scalar > 0.0) return 1.0;
    else if (scalar < 0.0) return -1.0;
    else return 0.0;
}

inline RNBoolean 
RNIsPositive(RNScalar scalar, RNMagnitude epsilon = RN_EPSILON)
{
    // Return whether scalar > 0
    return (scalar > epsilon);
}

inline RNBoolean 
RNIsNegative(RNScalar scalar, RNMagnitude epsilon = RN_EPSILON)
{
    // Return whether scalar < 0
    return (scalar < -epsilon);
}

inline RNBoolean 
RNIsPositiveOrZero(RNScalar scalar, RNMagnitude epsilon = RN_EPSILON)
{
    // Return whether scalar >= 0
    return (scalar >= -epsilon);
}

inline RNBoolean 
RNIsNegativeOrZero(RNScalar scalar, RNMagnitude epsilon = RN_EPSILON)
{
    // Return whether scalar <= 0
    return (scalar <= epsilon);
}

inline RNBoolean 
RNIsZero(RNScalar scalar, RNMagnitude epsilon = RN_EPSILON)
{
    // Return whether scalar == 0
    return (RNIsPositiveOrZero(scalar, epsilon) && RNIsNegativeOrZero(scalar, epsilon));
}

inline RNBoolean 
RNIsNotZero(RNScalar scalar, RNMagnitude epsilon = RN_EPSILON)
{
    // Return whether scalar != 0
    return (RNIsPositive(scalar, epsilon) || RNIsNegative(scalar, epsilon));
}

inline RNBoolean 
RNIsFinite(RNScalar scalar)
{
    // Return whether scalar is finite
    return ((-RN_INFINITY < scalar) && (scalar < RN_INFINITY));
}

inline RNBoolean 
RNIsInfinite(RNScalar scalar)
{
    // Return whether scalar is finite
    return ((scalar <= -RN_INFINITY) || (RN_INFINITY <= scalar));
}

inline RNBoolean 
RNIsEqual(RNScalar scalar1, RNScalar scalar2, RNMagnitude epsilon = RN_EPSILON)
{
    // Return whether scalar1 == scalar2
    return RNIsZero(scalar1 - scalar2, epsilon);
}

inline RNBoolean 
RNIsNotEqual(RNScalar scalar1, RNScalar scalar2, RNMagnitude epsilon = RN_EPSILON)
{
    // Return whether scalar1 != scalar2
    return RNIsNotZero(scalar1 - scalar2, epsilon);
}

inline RNBoolean 
RNIsGreater(RNScalar scalar1, RNScalar scalar2, RNMagnitude epsilon = RN_EPSILON)
{
    // Return whether scalar1 > scalar2
    return RNIsPositive(scalar1 - scalar2, epsilon);
}

inline RNBoolean 
RNIsLess(RNScalar scalar1, RNScalar scalar2, RNMagnitude epsilon = RN_EPSILON)
{
    // Return whether scalar1 < scalar2
    return RNIsNegative(scalar1 - scalar2, epsilon);
}

inline RNBoolean 
RNIsGreaterOrEqual(RNScalar scalar1, RNScalar scalar2, RNMagnitude epsilon = RN_EPSILON)
{
    // Return whether scalar1 >= scalar2
    return RNIsPositiveOrZero(scalar1 - scalar2, epsilon);
}

inline RNBoolean 
RNIsLessOrEqual(RNScalar scalar1, RNScalar scalar2, RNMagnitude epsilon = RN_EPSILON)
{
    // Return whether scalar1 <= scalar2
    return RNIsNegativeOrZero(scalar1 - scalar2, epsilon);
}



