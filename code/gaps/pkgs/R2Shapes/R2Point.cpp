/* Source file for the GAPS point class */



/* Include files */

#include "R2Shapes/R2Shapes.h"



/* Public variables */

const R2Point R2null_point(0.0, 0.0);
const R2Point R2ones_point(1.0, 1.0);
const R2Point R2posx_point(1.0, 0.0);
const R2Point R2posy_point(0.0, 1.0);
const R2Point R2negx_point(-1.0, 0.0);
const R2Point R2negy_point(0.0, -1.0);
const R2Point R2infinite_point(RN_INFINITY, RN_INFINITY);
const R2Point R2unknown_point(RN_UNKNOWN, RN_UNKNOWN);



/* Public functions */

int 
R2InitPoint()
{
    // Return success 
    return TRUE;
}



void R2StopPoint()
{
}



R2Point::
R2Point(void)
{
}



R2Point::
R2Point(RNCoord x, RNCoord y)
{
    v[0] = x; 
    v[1] = y; 
}



R2Point::
R2Point(const R2Point& point)
{
    v[0] = point.v[0]; 
    v[1] = point.v[1]; 
}



R2Point::
R2Point(const RNCoord array[2])
{
    v[0] = array[0]; 
    v[1] = array[1]; 
}



const RNBoolean R2Point::
IsFinite(void) const
{
    // Return whether point is finite
    return (RNIsFinite(v[0]) && RNIsFinite(v[1]));
}



const RNBoolean R2Point::
Collinear(const R2Point& point1, const R2Point& point2) const
{
    RNAbort("Not Implemented");
    return FALSE;
}



void R2Point::
Project(const R2Line& line)
{
    // Move point to closest point on line
    RNAbort("Not Implemented");
}



void R2Point::
Mirror(const R2Line& line)
{
    // Mirror point across line
    RNScalar d = R2SignedDistance(line, *this);
    *this += line.Normal() * (-2.0 * d);
}



void R2Point::
Rotate(const R2Point& origin, RNAngle theta)
{
    // Rotate point counterclockwise around axis through origin by radians ???
    RNAbort("Not Implemented");
}



void R2Point::
Transform(const R2Transformation& transformation)
{
    // Transform point
    transformation.Apply(*this);
}



void R2Point::
InverseTransform (const R2Transformation& transformation)
{
    // Inverse transform
    transformation.ApplyInverse(*this);
}



R2Point& R2Point::
operator=(const R2Point& point)
{
    v[0] = point[0];
    v[1] = point[1];
    return *this;
}



R2Point& R2Point::
operator+=(const R2Point& point)
{
    v[0] += point[0];
    v[1] += point[1];
    return *this;
}



R2Point& R2Point::
operator+=(const R2Vector& vector)
{
    v[0] += vector[0];
    v[1] += vector[1];
    return *this;
}



R2Point& R2Point::
operator-=(const R2Vector& vector)
{
    v[0] -= vector[0];
    v[1] -= vector[1];
    return *this;
}



R2Point& R2Point::
operator*=(const RNScalar a)
{
    v[0] *= a;
    v[1] *= a;
    return *this;
}



R2Point& R2Point::
operator*=(const R3Matrix& m)
{
    // Multiply point by matrix 
    RNCoord x = v[0], y = v[1];
    v[0] = x*m[0][0] + y*m[1][0] + m[2][0];
    v[1] = x*m[0][1] + y*m[1][1] + m[2][1];
    return *this;
}



R2Point& R2Point::
operator/=(const RNScalar a)
{
    //  assert(!zero(a)); 
    v[0] /= a;
    v[1] /= a;
    return *this;
}



R2Point 
operator+(const R2Point& point) 
{
    return point;
}



R2Point 
operator-(const R2Point& point)
{
    return R2Point(-point.X(), 
		   -point.Y());
}



R2Point 
operator+(const R2Point& point1, const R2Point& point2)
{
    return R2Point(point1.v[0] + point2.v[0], 
		   point1.v[1] + point2.v[1]);
}



R2Point 
operator+(const R2Point& point, const R2Vector& vector)
{
    return R2Point(point.X() + vector.X(), 
		   point.Y() + vector.Y());
}



R2Vector 
operator-(const R2Point& point1, const R2Point& point2)
{
    return R2Vector(point1.v[0] - point2.v[0], 
		    point1.v[1] - point2.v[1]);
}



R2Point 
operator-(const R2Point& point, const R2Vector& vector)
{
    return R2Point(point.X() - vector.X(), 
		   point.Y() - vector.Y());
}



R2Point 
operator*(const R2Point& point, const R3Matrix& m)
{
    // Multiply point by matrix 
    return R2Point(point.X()*m[0][0] + point.Y()*m[1][0] + m[2][0],
		   point.X()*m[0][1] + point.Y()*m[1][1] + m[2][1]);
}



R2Point 
operator*(const R2Point& point, const RNScalar a)
{
    return R2Point(point.X() * a, 
		   point.Y() * a);
}



R2Point 
operator/(const R2Point& point, const RNScalar a)
{
    assert(!RNIsZero(a));
    return R2Point(point.X() / a, 
		   point.Y() / a);
}



const RNBoolean R2Point::
operator==(const R2Point& point) const
{
    // Return whether point is equal
    return ((v[0] == point.v[0]) && (v[1] == point.v[1]));
}



const RNBoolean R2Point::
operator!=(const R2Point& point) const
{
    // Return whether point is not equal
    return ((v[0] != point.v[0]) || (v[1] != point.v[1]));
}




