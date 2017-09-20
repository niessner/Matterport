/* Source file for the R2 vector class */



/* Include files */

#include "R2Shapes/R2Shapes.h"



/* Public variables */

const R2Vector R2null_vector(0.0, 0.0);
const R2Vector R2ones_vector(1.0, 1.0);
const R2Vector R2posx_vector(1.0, 0.0);
const R2Vector R2posy_vector(0.0, 1.0);
const R2Vector R2negx_vector(-1.0, 0.0);
const R2Vector R2negy_vector(0.0, -1.0);
const R2Vector R2unknown_vector(RN_UNKNOWN, RN_UNKNOWN);



/* Public functions */

int 
R2InitVector()
{
    /* Return success */
    return TRUE;
}



void 
R2StopVector()
{
}



R2Vector::
R2Vector(void)
{
}



R2Vector::
R2Vector(RNCoord x, RNCoord y)
{
    v[0] = x; 
    v[1] = y; 
}



R2Vector::
R2Vector(const R2Vector& vector)
{
    v[0] = vector.v[0]; 
    v[1] = vector.v[1]; 
}



R2Vector::
R2Vector(RNAngle angle)
{
    v[0] = cos(angle);
    v[1] = sin(angle); 
}



R2Vector::
R2Vector(const RNCoord array[2])
{
    v[0] = array[0]; 
    v[1] = array[1]; 
}



const RNBoolean R2Vector::
IsFinite(void) const
{
    // Return whether vector is finite
    return (RNIsFinite(v[0]) && RNIsFinite(v[1]));
}



const RNLength R2Vector::
Length(void) const
{
    return sqrt((v[0]*v[0]) + (v[1]*v[1]));
}



const RNAngle R2Vector::
Angle(void) const
{
    // Return angle to positive x axis
    RNLength length = Length();
    if (length == 0) return 0;
    RNCoord x = v[0] / length;
    if (x == 0) return RN_PI_OVER_TWO;
    RNCoord y = v[1] / length;
    return atan2(y, x);
}



const R2Point R2Vector::
Point(void) const
{
    // Return point at (0,0) plus vector
    return R2Point(v[0], v[1]);
}



const RNDimension R2Vector::
MaxDimension(void) const
{
    // Return principal dimension of vector
    if (fabs(v[0]) >= fabs(v[1])) return RN_X;
    else return RN_Y;
}



const RNQuadrant R2Vector::
Quadrant(void) const
{
    // Return quadrant vector points into
    if (v[0] >= 0.0) {
        if (v[1] >= 0.0) return RN_PP_QUADRANT;
	else return RN_PN_QUADRANT;
    }
    else {
        if (v[1] >= 0.0) return RN_NP_QUADRANT;
	else return RN_NN_QUADRANT;
    }
}



const RNScalar R2Vector::
Dot(const R2Vector& vector) const
{
    return (v[0]*vector.v[0]) + (v[1]*vector.v[1]);
}



const RNScalar R2Vector::
Cross(const R2Vector& vector) const
{
    return (v[0]*vector.v[1]) - (v[1]*vector.v[0]);
}



const RNBoolean R2Vector::
operator==(const R2Vector& vector) const
{
    // Return whether vector is equal
    return ((v[0] == vector.v[0]) && (v[1] == vector.v[1]));
}



const RNBoolean R2Vector::
operator!=(const R2Vector& vector) const
{
    // Return whether vector is not equal
    return ((v[0] != vector.v[0]) || (v[1] != vector.v[1]));
}



void R2Vector::
Normalize(void)
{
    RNLength length = Length();
    if (length == 0.0) return;
    v[0] /= length;
    v[1] /= length;
}



void R2Vector::
Flip (void) 
{
    // Flip vector direction
    v[0] = -v[0];
    v[1] = -v[1];
}



void R2Vector::
Scale(RNScalar a)
{
    v[0] *= a;
    v[1] *= a;
}



void R2Vector::
Rotate(RNAngle theta)
{
    // Rotate vector counterclockwise 
    RNScalar c = cos(theta);
    RNScalar s = sin(theta);
    RNCoord x = v[0], y = v[1];
    v[0] = c*x - s*y;
    v[1] = s*x + c*y;
}



void R2Vector::
Project(const R2Vector& vector)
{
    // Project onto another vector    
    RNScalar dot = this->Dot(vector);
    RNLength length = vector.Length();
    if ((RNIsPositive(length)) && (length != 1.0)) 
	dot /= (length * length);
    *this = vector * dot;
}



void R2Vector::
Mirror(const R2Line& line)
{
    // Mirror vector across line
    RNScalar d = Dot(line.Normal());
    *this += line.Normal() * (-2.0 * d);
}



void R2Vector::
Transform(const R2Transformation& transformation)
{
    // Transform point
    transformation.Apply(*this);
}



void R2Vector::
InverseTransform (const R2Transformation& transformation)
{
    // Inverse transform
    transformation.ApplyInverse(*this);
}



R2Vector& R2Vector::
operator=(const R2Vector& vector)
{
    v[0] = vector.v[0];
    v[1] = vector.v[1];
    return *this;
}



R2Vector& R2Vector::
operator+=(const R2Vector& vector)
{
    v[0] += vector.v[0];
    v[1] += vector.v[1];
    return *this;
}



R2Vector& R2Vector::
operator-=(const R2Vector& vector)
{
    v[0] -= vector.v[0];
    v[1] -= vector.v[1];
    return *this;
}



R2Vector& R2Vector::
operator*=(const RNScalar a)
{
    v[0] *= a;
    v[1] *= a;
    return *this;
}



R2Vector& R2Vector::
operator*=(const R2Vector& vector)
{
    // Entry by entry multiply (not dot or cross product)
    v[0] *= vector.v[0];
    v[1] *= vector.v[1];
    return *this;
}



R2Vector& R2Vector::
operator*=(const R3Matrix& m)
{
    // Multiply vector by matrix 
    RNCoord x = v[0], y = v[1];
    v[0] = x*m[0][0] + y*m[1][0];
    v[1] = x*m[0][1] + y*m[1][1];
    return *this;
}



R2Vector& R2Vector::
operator/=(const RNScalar a)
{
    assert(RNIsNotZero(a));
    v[0] /= a;
    v[1] /= a;
    return *this;
}



R2Vector& R2Vector::
operator/=(const R2Vector& vector)
{
    // Entry by entry divide
    assert(RNIsNotZero(vector.v[0]));
    v[0] /= vector.v[0];
    assert(RNIsNotZero(vector.v[1]));
    v[1] /= vector.v[1];
    return *this;
}



R2Vector 
operator+(const R2Vector& vector) 
{
    return vector;
}



R2Vector 
operator-(const R2Vector& vector)
{
    return R2Vector(-vector.X(), 
		    -vector.Y());
}



R2Vector 
operator+(const R2Vector& vector1, const R2Vector& vector2)
{
    return R2Vector(vector1.v[0] + vector2.v[0], 
		    vector1.v[1] + vector2.v[1]);
}



R2Vector 
operator-(const R2Vector& vector1, const R2Vector& vector2)
{
    return R2Vector(vector1.v[0] - vector2.v[0], 
		    vector1.v[1] - vector2.v[1]);
}



R2Vector 
operator*(const R2Vector& vector1, const R2Vector& vector2)
{
    // Entry by entry multiply (not dot or cross product)
    return R2Vector(vector1.v[0] * vector2.v[0], 
		    vector1.v[1] * vector2.v[1]);
}



R2Vector 
operator*(const R2Vector& vector, const RNScalar a)
{
    return R2Vector(vector.X() * a, 
		    vector.Y() * a);
}



R2Vector 
operator*(const R2Vector& vector, const R3Matrix& m)
{
    // Multiply point by matrix 
    return R2Vector(vector.X()*m[0][0] + vector.Y()*m[1][0],
		    vector.X()*m[0][1] + vector.Y()*m[1][1]);
}



R2Vector 
operator/(const R2Vector& vector1, const R2Vector& vector2)
{
    assert(RNIsNotZero(vector2.v[0]));
    assert(RNIsNotZero(vector2.v[1]));
    return R2Vector(vector1.v[0] / vector2.v[0], 
		    vector1.v[1] / vector2.v[1]);
}



R2Vector 
operator/(const R2Vector& vector, const RNScalar a)
{
    assert(RNIsNotZero(a));
    return R2Vector(vector.X() / a, 
		    vector.Y() / a);
}



RNScalar
operator%(const R2Vector& vector1, const R2Vector& vector2)
{
    // Return cross product
    return vector1.X()*vector2.Y() - vector1.Y()*vector2.X();
}












