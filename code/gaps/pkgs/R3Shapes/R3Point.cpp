/* Source file for the GAPS point class */



/* Include files */

#include "R3Shapes/R3Shapes.h"



/* Public variables */

const R3Point R3null_point(0.0, 0.0, 0.0);
const R3Point R3ones_point(1.0, 1.0, 1.0);
const R3Point R3posx_point(1.0, 0.0, 0.0);
const R3Point R3posy_point(0.0, 1.0, 0.0);
const R3Point R3posz_point(0.0, 0.0, 1.0);
const R3Point R3negx_point(-1.0, 0.0, 0.0);
const R3Point R3negy_point(0.0, -1.0, 0.0);
const R3Point R3negz_point(0.0, 0.0, -1.0);
const R3Point R3infinity_point(RN_INFINITY, RN_INFINITY, RN_INFINITY);
const R3Point R3unknown_point(RN_UNKNOWN, RN_UNKNOWN, RN_UNKNOWN);



/* Public functions */

int 
R3InitPoint()
{
    // Return success 
    return TRUE;
}



void R3StopPoint()
{
}



R3Point::
R3Point(void)
{
}



R3Point::
R3Point(RNCoord x, RNCoord y, RNCoord z)
{
    v[0] = x; 
    v[1] = y; 
    v[2] = z;
}



R3Point::
R3Point(const R3Point& point)
{
    v[0] = point.v[0]; 
    v[1] = point.v[1]; 
    v[2] = point.v[2]; 
}



R3Point::
R3Point(const RNCoord array[3])
{
    v[0] = array[0]; 
    v[1] = array[1]; 
    v[2] = array[2];
}



const RNBoolean R3Point::
IsFinite(void) const
{
    // Return whether point is finite
    return (RNIsFinite(v[0]) && RNIsFinite(v[1]) && RNIsFinite(v[2]));
}



const R3Box R3Point::
BBox (void) const
{
    // Return bounding box
    return R3Box(*this, *this);
}



const R3Sphere R3Point::
BSphere (void) const
{
    // Return bounding sphere
    return R3Sphere(*this, 0.0);
}



const RNBoolean R3Point::
operator==(const R3Point& point) const
{
    // Return whether point is equal
    return ((v[0] == point.v[0]) && (v[1] == point.v[1]) && (v[2] == point.v[2]));
}



const RNBoolean R3Point::
operator!=(const R3Point& point) const
{
    // Return whether point is not equal
    return ((v[0] != point.v[0]) || (v[1] != point.v[1]) || (v[2] != point.v[2]));
}



const RNBoolean R3Point::
Collinear(const R3Point& point1, const R3Point& point2) const
{
    // Check if two of points are same
    if ((*this == point1) || (*this == point2) || (point1 == point2)) return TRUE;

    /// Check if three points are collinear
    R3Vector v = point1 - *this;
    v.Cross(point1 - point2);
    if (RNIsZero(v.Length())) return TRUE;
    return FALSE;
}



void R3Point::
Project(const R3Line& line)
{
    // Move point to closest point on line
    const R3Point *p = &(line.Point());
    const R3Vector *v = &(line.Vector());
    RNScalar denom = v->Dot(*v);
    if (RNIsZero(denom)) return;
    RNScalar t = (v->X() * (X() - p->X()) + v->Y() *(Y() - p->Y()) + v->Z() * (Z() - p->Z())) / denom;
    *this = *p + *v * t;
}



void R3Point::
Project(const R3Plane& plane)
{
    // Move point to closest point on plane
    RNScalar d = R3SignedDistance(plane, *this);
    *this += plane.Normal() * -d;
}



void R3Point::
Mirror(const R3Plane& plane)
{
    // Mirror point across plane
    RNScalar d = R3SignedDistance(plane, *this);
    *this += plane.Normal() * (-2.0 * d);
}



void R3Point:: 
XRotate(RNAngle radians)
{
    // rotate matrix around X axis counterclockwise
    RNScalar c = cos(radians);
    RNScalar s = sin(radians);
    RNScalar y = v[1];
    RNScalar z = v[2];
    v[1] = c*y - s*z;
    v[2] = s*y + c*z;
}



void R3Point:: 
YRotate(RNAngle radians)
{
    // rotate matrix around Y axis counterclockwise
    RNScalar c = cos(radians);
    RNScalar s = sin(radians);
    RNScalar x = v[0];
    RNScalar z = v[2];
    v[0] =  c*x + s*z;
    v[2] = -s*x + c*z;
}



void R3Point:: 
ZRotate(RNAngle radians)
{
    // rotate matrix around Z axis counterclockwise 
    RNScalar c = cos(radians);
    RNScalar s = sin(radians);
    RNScalar x = v[0];
    RNScalar y = v[1];
    v[0] = c*x - s*y;
    v[1] = s*x + c*y;
}



void R3Point:: 
Rotate(RNAxis axis, RNAngle radians)
{
    // rotate point around an axis counterclockwise
    switch (axis) {
    case RN_XAXIS: 
	XRotate(radians); 
	break;

    case RN_YAXIS: 
	YRotate(radians); 
	break;

    case RN_ZAXIS: 
	ZRotate(radians); 
	break;

    default: 
	RNWarning("Matrix rotation around undefined axis");
	break;
    }
}



void R3Point:: 
Rotate(const R3Vector& radians)
{
    // Rotate first around X, then around Y, and finally around Z
    ZRotate(radians.Z());
    YRotate(radians.Y());
    XRotate(radians.X());
}



void R3Point:: 
Rotate(const R3Quaternion& quaternion)
{
    // Rotate by quaternion
    *this = quaternion * (*this);
}



void R3Point::
Rotate(const R3Vector& axis, RNAngle theta)
{
    // Rotate point counterclockwise around axis through origin by radians ???
    R3Vector v = Vector();
    v.Rotate(axis, theta);
    *this = v.Point();
}



void R3Point::
Rotate(const R3Line& axis, RNAngle theta)
{
    // Translate axis to origin
    R3Vector v = *this - axis.Point();

    // Rotate point counterclockwise around axis through origin by radians ???
    v.Rotate(axis.Vector(), theta);

    // Translate axis back from origin
    *this = axis.Point() + v;
}



void R3Point::
Transform(const R3Transformation& transformation)
{
    // Transform point
    transformation.Apply(*this);
}



void R3Point::
InverseTransform(const R3Transformation& transformation)
{
    // Transform vector
    transformation.ApplyInverse(*this);
}



R3Point& R3Point::
operator=(const R3Point& point)
{
    v[0] = point.v[0];
    v[1] = point.v[1];
    v[2] = point.v[2];
    return *this;
}



R3Point& R3Point::
operator+=(const R3Point& point)
{
    v[0] += point[0];
    v[1] += point[1];
    v[2] += point[2];
    return *this;
}



R3Point& R3Point::
operator+=(const R3Vector& vector)
{
    v[0] += vector[0];
    v[1] += vector[1];
    v[2] += vector[2];
    return *this;
}



R3Point& R3Point::
operator-=(const R3Vector& vector)
{
    v[0] -= vector[0];
    v[1] -= vector[1];
    v[2] -= vector[2];
    return *this;
}



R3Point& R3Point::
operator*=(const RNScalar a)
{
    v[0] *= a;
    v[1] *= a;
    v[2] *= a;
    return *this;
}



R3Point& R3Point::
operator/=(const RNScalar a)
{
    assert(a != 0);
    v[0] /= a;
    v[1] /= a;
    v[2] /= a;
    return *this;
}



R3Point 
operator-(const R3Point& point)
{
    return R3Point(-point.v[0], 
		   -point.v[1], 
		   -point.v[2]);
}



R3Point 
operator+(const R3Point& point1, const R3Point& point2)
{
    return R3Point(point1.v[0] + point2.v[0], 
		   point1.v[1] + point2.v[1], 
		   point1.v[2] + point2.v[2]);
}



R3Point 
operator+(const R3Point& point, const R3Vector& vector)
{
    return R3Point(point.X() + vector.X(), 
		   point.Y() + vector.Y(), 
		   point.Z() + vector.Z());
}



R3Vector 
operator-(const R3Point& point1, const R3Point& point2)
{
    return R3Vector(point1.v[0] - point2.v[0], 
		    point1.v[1] - point2.v[1], 
		    point1.v[2] - point2.v[2]);
}



R3Point 
operator-(const R3Point& point, const R3Vector& vector)
{
    return R3Point(point.X() - vector.X(), 
		   point.Y() - vector.Y(), 
		   point.Z() - vector.Z());
}



R3Point 
operator*(const R3Point& point, const RNScalar a)
{
    return R3Point(point.X() * a, 
		   point.Y() * a,
		   point.Z() * a);
}



R3Point 
operator/(const R3Point& point, const RNScalar a)
{
    assert(a != 0);
    return R3Point(point.X() / a, 
		   point.Y() / a, 
		   point.Z() / a);
}






