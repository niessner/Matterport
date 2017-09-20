/* Source file for the GAPS vector class */



/* Include files */

#include "R3Shapes/R3Shapes.h"



/* Public variables */

const R3Vector R3null_vector(0.0, 0.0, 0.0);
const R3Vector R3ones_vector(1.0, 1.0, 1.0);
const R3Vector R3posx_vector(1.0, 0.0, 0.0);
const R3Vector R3posy_vector(0.0, 1.0, 0.0);
const R3Vector R3posz_vector(0.0, 0.0, 1.0);
const R3Vector R3negx_vector(-1.0, 0.0, 0.0);
const R3Vector R3negy_vector(0.0, -1.0, 0.0);
const R3Vector R3negz_vector(0.0, 0.0, -1.0);
const R3Vector R3unknown_vector(RN_UNKNOWN, RN_UNKNOWN, RN_UNKNOWN);



/* Public functions */

int 
R3InitVector()
{
    /* Return success */
    return TRUE;
}



void 
R3StopVector()
{
}



R3Vector::
R3Vector(void)
{
}



R3Vector::
R3Vector(RNCoord x, RNCoord y, RNCoord z)
{
    v[0] = x; 
    v[1] = y; 
    v[2] = z;
}



R3Vector::
R3Vector(const R3Vector& vector)
{
    v[0] = vector.v[0]; 
    v[1] = vector.v[1]; 
    v[2] = vector.v[2]; 
}



R3Vector::
R3Vector(RNAngle pitch, RNAngle yaw)
{
    RNScalar cosine_yaw = cos(yaw);
    v[0] = cosine_yaw * cos(pitch);
    v[1] = cosine_yaw * sin(pitch); 
    v[2] = sin(yaw);
}



R3Vector::
R3Vector(const RNCoord array[3])
{
    v[0] = array[0]; 
    v[1] = array[1]; 
    v[2] = array[2];
}



const RNBoolean R3Vector::
IsFinite(void) const
{
    // Return whether vector is finite
    return (RNIsFinite(v[0]) && RNIsFinite(v[1]) && RNIsFinite(v[2]));
}



const RNLength R3Vector::
Length(void) const
{
    return (sqrt((v[0]*v[0]) + (v[1]*v[1]) + (v[2]*v[2])));
}



const R3Point R3Vector::
Point(void) const
{
    // Return point at (0,0,0) plus vector
    return R3Point(v[0], v[1], v[2]);
}



const RNDimension R3Vector::
MinDimension(void) const
{
    // Return principal dimension of vector
    if (fabs(v[0]) <= fabs(v[1])) {
        if (fabs(v[0]) <= fabs(v[2])) return RN_X;
        else return RN_Z;
    }
    else {
        if (fabs(v[1]) <= fabs(v[2])) return RN_Y;
        else return RN_Z;
    }
}



const RNDimension R3Vector::
MaxDimension(void) const
{
    // Return principal dimension of vector
    if (fabs(v[0]) >= fabs(v[1])) {
        if (fabs(v[0]) >= fabs(v[2])) return RN_X;
        else return RN_Z;
    }
    else {
        if (fabs(v[1]) >= fabs(v[2])) return RN_Y;
        else return RN_Z;
    }
}



const RNSextant R3Vector::
Sextant(void) const
{
    // Return principal direction of vector
    if (fabs(v[0]) >= fabs(v[1])) {
        if (fabs(v[0]) >= fabs(v[2])) {
            if (v[0] >= 0.0) return RN_PX_SEXTANT;
            else return RN_NX_SEXTANT;
        }
        else {
            if (v[2] >= 0.0) return RN_PZ_SEXTANT;
            else return RN_NZ_SEXTANT;
        }
    }
    else {
        if (fabs(v[1]) >= fabs(v[2])) {
            if (v[1] >= 0.0) return RN_PY_SEXTANT;
            else return RN_NY_SEXTANT;
        }
        else {
            if (v[2] >= 0.0) return RN_PZ_SEXTANT;
            else return RN_NZ_SEXTANT;
        }
    }
}



const RNOctant R3Vector::
Octant(void) const
{
    // Return octant vector points into
    if (v[0] >= 0.0) {
        if (v[1] >= 0.0) {
	    if (v[2] >= 0.0) return RN_PPP_OCTANT;
            else return RN_PPN_OCTANT;
        }
	else {
	    if (v[2] >= 0.0) return RN_PNP_OCTANT;
            else return RN_PNN_OCTANT;
        }
    }
    else {
        if (v[1] >= 0.0) {
	    if (v[2] >= 0.0) return RN_NPP_OCTANT;
            else return RN_NPN_OCTANT;
        }
	else {
	    if (v[2] >= 0.0) return RN_NNP_OCTANT;
            else return RN_NNN_OCTANT;
        }
    }
}



const RNScalar R3Vector::
Dot(const R3Vector& vector) const
{
    return((v[0]*vector.v[0]) + (v[1]*vector.v[1]) + (v[2]*vector.v[2]));
}



const RNBoolean R3Vector::
operator==(const R3Vector& vector) const
{
    // Return whether vector is equal
    return ((v[0] == vector.v[0]) && (v[1] == vector.v[1]) && (v[2] == vector.v[2]));
}



const RNBoolean R3Vector::
operator!=(const R3Vector& vector) const
{
    // Return whether vector is not equal
    return ((v[0] != vector.v[0]) || (v[1] != vector.v[1]) || (v[2] != vector.v[2]));
}



void R3Vector::
Normalize(void)
{
    RNLength length = Length();
    if (length == 0.0) return;
    v[0] /= length;
    v[1] /= length;
    v[2] /= length;
}



void R3Vector::
Flip (void) 
{
    // Flip vector direction
    v[0] = -v[0];
    v[1] = -v[1];
    v[2] = -v[2];
}



void R3Vector::
Cross(const R3Vector& vector)
{
    const RNScalar x = (v[1]*vector.v[2]) - (v[2]*vector.v[1]);
    const RNScalar y = (v[2]*vector.v[0]) - (v[0]*vector.v[2]);
    const RNScalar z = (v[0]*vector.v[1]) - (v[1]*vector.v[0]);
    v[2] = z; v[1] = y; v[0] = x; 
}



void R3Vector:: 
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



void R3Vector:: 
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



void R3Vector:: 
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



void R3Vector:: 
Rotate(RNAxis axis, RNAngle radians)
{
    // rotate vector around an axis counterclockwise
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



void R3Vector:: 
Rotate(const R3Vector& radians)
{
    // Rotate first around X, then around Y, and finally around Z
    ZRotate(radians.Z());
    YRotate(radians.Y());
    XRotate(radians.X());
}



void R3Vector:: 
Rotate(const R3Quaternion& quaternion)
{
    // Rotate by quaternion
    *this = quaternion * (*this);
}



void R3Vector::
Rotate(const R3Vector& axis, RNAngle theta)
{
    // Rotate vector counterclockwise around axis (looking at axis end-on) (rz(xaxis) = yaxis)
    // From Goldstein: v' = v cos t + a (v . a) [1 - cos t] - (v x a) sin t
    if (RNIsEqual(fabs(Dot(axis)), 1.0)) return;
    const RNScalar cos_theta = cos(theta);
    const RNScalar dot = this->Dot(axis);
    R3Vector cross = *this % axis;
    *this *= cos_theta;
    *this += axis * dot * (1.0 - cos_theta);
    *this -= cross * sin(theta); 
}



void R3Vector::
Project(const R3Vector& vector)
{
    // Project onto another vector    
    RNScalar dot = this->Dot(vector);
    RNLength length = vector.Length();
    if ((RNIsPositive(length)) && (length != 1.0)) 
	dot /= (length * length);
    *this = vector * dot;
}



void R3Vector::
Project(const R3Plane& plane) 
{
    // Project onto plane
    *this -= plane.Normal() * this->Dot(plane.Normal());
}



void R3Vector::
Mirror(const R3Plane& plane)
{
    // Mirror vector across plane
    RNScalar d = Dot(plane.Normal());
    *this += plane.Normal() * (-2.0 * d);
}



void R3Vector::
Transform(const R3Transformation& transformation)
{
    // Transform vector
    transformation.Apply(*this);
}



void R3Vector::
InverseTransform(const R3Transformation& transformation)
{
    // Transform vector
    transformation.ApplyInverse(*this);
}



R3Vector& R3Vector::
operator=(const R3Vector& vector)
{
    v[0] = vector.v[0];
    v[1] = vector.v[1];
    v[2] = vector.v[2];
    return *this;
}



R3Vector& R3Vector::
operator+=(const R3Vector& vector)
{
    v[0] += vector.v[0];
    v[1] += vector.v[1];
    v[2] += vector.v[2];
    return *this;
}



R3Vector& R3Vector::
operator-=(const R3Vector& vector)
{
    v[0] -= vector.v[0];
    v[1] -= vector.v[1];
    v[2] -= vector.v[2];
    return *this;
}



R3Vector& R3Vector::
operator*=(const RNScalar a)
{
    v[0] *= a;
    v[1] *= a;
    v[2] *= a;
    return *this;
}



R3Vector& R3Vector::
operator*=(const R3Vector& vector)
{
    // Entry by entry multiply (not dot or cross product)
    v[0] *= vector.v[0];
    v[1] *= vector.v[1];
    v[2] *= vector.v[2];
    return *this;
}



R3Vector& R3Vector::
operator/=(const RNScalar a)
{
    assert(a != 0);
    v[0] /= a;
    v[1] /= a;
    v[2] /= a;
    return *this;
}



R3Vector& R3Vector::
operator/=(const R3Vector& vector)
{
    // Entry by entry divide
    assert(vector.v[0] != 0);
    v[0] /= vector.v[0];
    assert(vector.v[1] != 0);
    v[1] /= vector.v[1];
    assert(vector.v[2] != 0);
    v[2] /= vector.v[2];
    return *this;
}



R3Vector 
operator+(const R3Vector& vector)
{
    return vector;
}



R3Vector 
operator-(const R3Vector& vector)
{
    return R3Vector(-vector.v[0], 
		    -vector.v[1], 
		    -vector.v[2]);
}



R3Vector 
operator+(const R3Vector& vector1, const R3Vector& vector2)
{
    return R3Vector(vector1.v[0] + vector2.v[0], 
		    vector1.v[1] + vector2.v[1], 
		    vector1.v[2] + vector2.v[2]);
}



R3Vector 
operator-(const R3Vector& vector1, const R3Vector& vector2)
{
    return R3Vector(vector1.v[0] - vector2.v[0], 
		    vector1.v[1] - vector2.v[1], 
		    vector1.v[2] - vector2.v[2]);
}



R3Vector 
operator*(const R3Vector& vector1, const R3Vector& vector2)
{
    // Entry by entry multiply (not dot or cross product)
    return R3Vector(vector1.v[0] * vector2.v[0], 
		    vector1.v[1] * vector2.v[1], 
		    vector1.v[2] * vector2.v[2]);
}



R3Vector 
operator*(const R3Vector& vector, const RNScalar a)
{
    return R3Vector(vector.v[0] * a, 
		    vector.v[1] * a, 
		    vector.v[2] * a);
}



R3Vector 
operator/(const R3Vector& vector1, const R3Vector& vector2)
{
    assert(vector2.v[0] != 0);
    assert(vector2.v[1] != 0);
    assert(vector2.v[2] != 0);
    return R3Vector(vector1.v[0]/vector2.v[0], 
		    vector1.v[1]/vector2.v[1], 
		    vector1.v[2]/vector2.v[2]);
}



R3Vector 
operator/(const R3Vector& vector, const RNScalar a)
{
    assert(a != 0);
    return R3Vector(vector.v[0]/a, 
		    vector.v[1]/a, 
		    vector.v[2]/a);
}



R3Vector 
operator%(const R3Vector& vector1, const R3Vector& vector2)
{
    // Cross product 
    R3Vector v = vector1;
    v.Cross(vector2);
    return v;
}



R3Vector 
R3RandomDirection(void)
{
  RNScalar phi = RN_TWO_PI * RNRandomScalar();
  RNScalar z = RNRandomScalar();
  RNScalar theta = acos(z);
  RNScalar sin_theta = sin(theta);
  RNScalar x = sin_theta * cos(phi);
  RNScalar y = sin_theta * sin(phi);
  R3Vector vector(x, y, z);
  if (RNRandomScalar() < 0.5) vector.Flip();
  return vector;
}
