/* Source file for the R2 affine transformation class */



/* Include files */

#include "R2Shapes/R2Shapes.h"



/* Public variables */

const R2Affine R2null_affine(
    R3Matrix(0.0, 0.0, 0.0, 
	     0.0, 0.0, 0.0, 
	     0.0, 0.0, 0.0)
);

const R2Affine R2identity_affine(
    R3Matrix(1.0, 0.0, 0.0, 
	     0.0, 1.0, 0.0, 
	     0.0, 0.0, 1.0)
);



/* Public functions */

int R2InitAffine()
{
    /* Return success */
    return TRUE;
}



void R2StopAffine()
{
}



R2Affine::
R2Affine(void)
{
}



R2Affine::
R2Affine(const R3Matrix& matrix, RNBoolean mirror)
    : matrix(matrix),
      mirror(mirror)
{
}



const RNBoolean R2Affine::
IsMirrored(void) const
{
    // Return whether affine transform is mirrored
    return mirror;
}



const RNBoolean R2Affine::
IsAffine(void) const
{
    // Return whether transform is affine
    return TRUE;
}



const RNBoolean R2Affine::
IsIsotropic(void) const
{
    // Return whether affine transformation is isotropic ???
    return TRUE;
}



const RNScalar R2Affine::
ScaleFactor(void) const
{
    // Return scale factor imposed by transformation 
    R2Vector v(0.707106781186, 0.707106781186);
    Apply(v);
    return v.Length();
}



void R2Affine::
Apply (R2Vector& vector) const
{
    // Transform vector
    vector = matrix * vector;
}



void R2Affine::
Apply (R2Point& point) const
{
    // Transform point
    point = matrix * point;
}



void R2Affine::
Apply (R2Transformation& transformation) const
{
    // Transform abstract transformation
    transformation.Transform(*this);
}



void R2Affine::
Apply (R2Affine& affine) const
{
    // Transform affine transformation
    affine.Transform(*this);
}



void R2Affine::
ApplyInverse (R2Vector& vector) const
{
    // Transform vector
    vector = InverseMatrix() * vector;
}



void R2Affine::
ApplyInverse (R2Point& point) const
{
    // Transform point
    point = InverseMatrix() * point;
}



void R2Affine::
ApplyInverse (R2Transformation& transformation) const
{
    // Transform abstract transformation
    transformation.InverseTransform(*this);
}



void R2Affine::
ApplyInverse (R2Affine& affine) const
{
    // Transform affine transformation
    affine.InverseTransform(*this);
}



void R2Affine::
Transform(const R2Transformation& transformation) 
{
    // Transform affine transformation
    transformation.Apply(*this);
}



void R2Affine::
Transform(const R2Affine& affine) 
{
    // Transform affine transformation
    if (affine.mirror) mirror = !mirror;
    matrix = matrix * affine.matrix;
}



void R2Affine::
InverseTransform(const R2Transformation& transformation) 
{
    // Inverse transform affine transformation
    transformation.ApplyInverse(*this);
}



void R2Affine::
InverseTransform(const R2Affine& affine) 
{
    // Transform affine transformation
    if (affine.mirror) mirror = !mirror;
    matrix = matrix * affine.InverseMatrix();
}



void R2Affine::
Reset(const R2Transformation& transformation) 
{
    // Reset transformation
    *this = R2identity_affine;
    Transform(transformation);
}



void R2Affine::
Reset (const R3Matrix& matrix, RNBoolean mirror)
{
    // Reset transformation
    this->matrix = matrix;
    this->mirror = mirror;
}






