/* Source file for the R3 affine transformation class */



/* Include files */

#include "R3Shapes/R3Shapes.h"



/* Public variables */

const R3Affine R3null_affine(
    R4Matrix(0.0, 0.0, 0.0, 0.0,
	     0.0, 0.0, 0.0, 0.0,
	     0.0, 0.0, 0.0, 0.0,
	     0.0, 0.0, 0.0, 0.0)
);

const R3Affine R3identity_affine(
    R4Matrix(1.0, 0.0, 0.0, 0.0,
	     0.0, 1.0, 0.0, 0.0,
	     0.0, 0.0, 1.0, 0.0,
	     0.0, 0.0, 0.0, 1.0)
);



/* Public functions */

int R3InitAffine()
{
    /* Return success */
    return TRUE;
}



void R3StopAffine()
{
}



R3Affine::
R3Affine(void)
{
}



R3Affine::
R3Affine(const R4Matrix& matrix, RNBoolean mirror)
    : matrix(matrix),
      flags((mirror) ? R3_AFFINE_MIRROR_FLAG : 0)
{
}



R3Affine::
~R3Affine(void)
{
}



const R4Matrix& R3Affine::
InverseMatrix(void) const
{
    // Update inverse matrix 
    if (!flags[R3_AFFINE_INVERSE_UPTODATE_FLAG]) {
        // Update inverse matrix
        R3Affine *affine = (R3Affine *) this;
	affine->inverse = matrix.Inverse();
	affine->flags.Add(R3_AFFINE_INVERSE_UPTODATE_FLAG);
    }

    // Return inverse matrix
    return inverse;
}



const RNBoolean R3Affine::
IsIdentity(void) const
{
    // Return whether affine transform is identity
    return (matrix.IsIdentity() && !HasMirror());
}



const RNBoolean R3Affine::
IsMirrored(void) const
{
    // Return whether affine transform is mirrored
    return HasMirror();
}



const RNBoolean R3Affine::
IsAffine(void) const
{
    // Return whether transform is affine
    return TRUE;
}



const RNBoolean R3Affine::
IsIsotropic(void) const
{
    // Return whether affine transformation is isotropic
    return matrix.IsIsotropic();
}



const RNScalar R3Affine::
ScaleFactor(void) const
{
    // Return scale factor imposed by transformation 
    R3Vector v(0.57735027, 0.57735027, 0.57735027);
    Apply(v);
    return v.Length();
}



void R3Affine::
Apply (R3Vector& vector) const
{
    // Transform vector
    vector = matrix * vector;
}



void R3Affine::
Apply (R3Point& point) const
{
    // Transform point
    point = matrix * point;
}



void R3Affine::
Apply (R3Transformation& transformation) const
{
    // Transform abstract transformation
    transformation.Transform(*this);
}



void R3Affine::
Apply (R3Affine& affine) const
{
    // Transform affine transformation
    affine.Transform(*this);
}



void R3Affine::
ApplyInverse (R3Vector& vector) const
{
    // Transform vector by inverse
    vector = InverseMatrix() * vector;
}



void R3Affine::
ApplyInverse (R3Point& point) const
{
    // Transform point by inverse
    point = InverseMatrix() * point;
}



void R3Affine::
ApplyInverse (R3Transformation& transformation) const
{
    // Transform abstract transformation by inverse
    transformation.InverseTransform(*this);
}



void R3Affine::
ApplyInverse (R3Affine& affine) const
{
    // Transform affine transformation by inverse
    affine.InverseTransform(*this);
}



void R3Affine::
ApplyInverseTranspose (R3Vector& vector) const
{
    // Transform (normal) vector by inverse transpose
    const R4Matrix& m = InverseMatrix();
    RNCoord x = m[0][0] * vector.X() + m[1][0] * vector.Y() + m[2][0] * vector.Z();
    RNCoord y = m[0][1] * vector.X() + m[1][1] * vector.Y() + m[2][1] * vector.Z();
    RNCoord z = m[0][2] * vector.X() + m[1][2] * vector.Y() + m[2][2] * vector.Z();
    vector = R3Vector(x, y, z);
}



void R3Affine::
ApplyTranspose (R3Vector& vector) const
{
    // Transform (normal) vector by inverse transpose
    const R4Matrix& m = Matrix();
    RNCoord x = m[0][0] * vector.X() + m[1][0] * vector.Y() + m[2][0] * vector.Z();
    RNCoord y = m[0][1] * vector.X() + m[1][1] * vector.Y() + m[2][1] * vector.Z();
    RNCoord z = m[0][2] * vector.X() + m[1][2] * vector.Y() + m[2][2] * vector.Z();
    vector = R3Vector(x, y, z);
}



void R3Affine::
XMirror(void) 
{
    // Mirror across YZ plane
    matrix.XScale(-1.0);
    InvalidateInverse();
    Mirror();
}



void R3Affine::
YMirror(void) 
{
    // Mirror across XZ plane
    matrix.YScale(-1.0);
    InvalidateInverse();
    Mirror();
}



void R3Affine::
ZMirror(void) 
{
    // Mirror across XY plane
    matrix.ZScale(-1.0);
    InvalidateInverse();
    Mirror();
}



void R3Affine::
Mirror(RNAxis axis)
{
    // Mirror along axis
    matrix.Scale(axis, -1.0);
    Mirror();
}



void R3Affine::
Mirror(const R3Plane& plane)
{
  // Compute translation
  R3Vector translation = plane.D() * plane.Normal();

  // Compute dimension to rotate to and reflect across
  RNDimension dim = plane.Normal().MinDimension();
  R3Vector target_vector = R3xyz_triad[dim];

  // Compute rotation axis
  R3Vector rotation_axis = plane.Normal() % target_vector;
  RNLength rotation_axis_length = rotation_axis.Length();
  if (RNIsZero(rotation_axis_length)) return;
  rotation_axis /= rotation_axis_length;

  // Compute rotation angle
  RNAngle rotation_angle = R3InteriorAngle(plane.Normal(), target_vector);
  
  // Transform
  Translate(-translation);
  Rotate(rotation_axis, -rotation_angle);
  Mirror(dim);
  Rotate(rotation_axis, rotation_angle);
  Translate(translation);
}  



void R3Affine::
Transform(const R3Transformation& transformation) 
{
    // Transform affine transformation
    transformation.Apply(*this);
}



void R3Affine::
Transform(const R3Affine& affine) 
{
    // Transform affine transformation
    matrix = matrix * affine.matrix;
    if (affine.HasMirror()) Mirror();
    InvalidateInverse();
}



void R3Affine::
InverseTransform(const R3Transformation& transformation) 
{
    // Inverse transform affine transformation
    transformation.ApplyInverse(*this);
}



void R3Affine::
InverseTransform(const R3Affine& affine) 
{
    // Transform affine transformation
    matrix = matrix * affine.InverseMatrix();
    if (affine.HasMirror()) Mirror();
    InvalidateInverse();
}



void R3Affine::
Reset(const R3Transformation& transformation) 
{
    // Reset transformation
    *this = R3identity_affine;
    Transform(transformation);
}



void R3Affine::
Reset (const R4Matrix& matrix, RNBoolean mirror)
{
    // Reset transformation
    this->matrix = matrix;
    this->flags = (mirror) ? R3_AFFINE_MIRROR_FLAG : 0;
}
















