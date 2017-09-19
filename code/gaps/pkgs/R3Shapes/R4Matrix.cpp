/* Source file for the GAPS matrix class */

/*

NOTES:

Matrix concatenation order is designed to be compatible with OpenGL.
The idea is that transformations can be applied (matrices multiplied)
in the same order as OpenGL calls would be made.  Make more global
transformations first (i.e., ones at the top of the scene hierachy).
Vectors are column vectors on the right of matrices.  Transformation
matrices are post-multiplied.

Matrix storage is NOT COMPATIBLE with OpenGL.  OpenGL represents
matrices in COLUMN MAJOR ORDER.  I find manipulating matrices in this
format both inconvenient (operator[]) and confusing.  Therefore, this
package represents matices in ROW MAJOR ORDER.  In order to pass
matrices to OpenGL, a transpose is required.  Look into this re
efficiency ???

*/



/* Include files */

#include "R3Shapes/R3Shapes.h"



/* Public variables */

const R4Matrix R4null_matrix (
    0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0
);

const R4Matrix R4identity_matrix (
    1.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 1.0
);



/* Public functions */

int 
R3InitMatrix()
{
    /* Return success */
    return TRUE;
}



void 
R3StopMatrix()
{
}



R4Matrix::
R4Matrix(void)
{
}



R4Matrix::
R4Matrix(const R4Matrix& matrix)
{
    // Assign matrix entries
    m[0][0] = matrix.m[0][0]; m[0][1] = matrix.m[0][1]; m[0][2] = matrix.m[0][2]; m[0][3] = matrix.m[0][3];
    m[1][0] = matrix.m[1][0]; m[1][1] = matrix.m[1][1]; m[1][2] = matrix.m[1][2]; m[1][3] = matrix.m[1][3];
    m[2][0] = matrix.m[2][0]; m[2][1] = matrix.m[2][1]; m[2][2] = matrix.m[2][2]; m[2][3] = matrix.m[2][3];
    m[3][0] = matrix.m[3][0]; m[3][1] = matrix.m[3][1]; m[3][2] = matrix.m[3][2]; m[3][3] = matrix.m[3][3];
}



R4Matrix::
R4Matrix(RNScalar a00, RNScalar a01, RNScalar a02, RNScalar a03,
         RNScalar a10, RNScalar a11, RNScalar a12, RNScalar a13,
	 RNScalar a20, RNScalar a21, RNScalar a22, RNScalar a23,
	 RNScalar a30, RNScalar a31, RNScalar a32, RNScalar a33)
{
    // Assign matrix entries
    m[0][0] = a00; m[0][1] = a01; m[0][2] = a02; m[0][3] = a03;
    m[1][0] = a10; m[1][1] = a11; m[1][2] = a12; m[1][3] = a13;
    m[2][0] = a20; m[2][1] = a21; m[2][2] = a22; m[2][3] = a23;
    m[3][0] = a30; m[3][1] = a31; m[3][2] = a32; m[3][3] = a33;
}



R4Matrix::
R4Matrix(const RNScalar *a)
{
    // Assign matrix entries
    m[0][0] = a[0]; m[0][1] = a[1]; m[0][2] = a[2]; m[0][3] = a[3];
    m[1][0] = a[4]; m[1][1] = a[5]; m[1][2] = a[6]; m[1][3] = a[7];
    m[2][0] = a[8]; m[2][1] = a[9]; m[2][2] = a[10]; m[2][3] = a[11];
    m[3][0] = a[12]; m[3][1] = a[13]; m[3][2] = a[14]; m[3][3] = a[15];
}



const RNBoolean R4Matrix::
IsZero (void) const
{
    // Return whether matrix is zero
    return (*this == R4null_matrix);
}



const RNBoolean R4Matrix::
IsIdentity (void) const
{
    // Return whether matrix is identity
    return (*this == R4identity_matrix);
}



const RNBoolean R4Matrix::
IsIsotropic(void) const
{
    // Return whether matrix transformation is isotropic
    RNScalar d0 = m[0][0]*m[0][0] + m[1][0]*m[1][0] + m[2][0]*m[2][0];
    RNScalar d1 = m[0][1]*m[0][1] + m[1][1]*m[1][1] + m[2][1]*m[2][1];
    if (RNIsNotEqual(d0, d1)) return FALSE;
    RNScalar d2 = m[0][2]*m[0][2] + m[1][2]*m[1][2] + m[2][2]*m[2][2];
    if (RNIsNotEqual(d0, d2)) return FALSE;
    return TRUE;
}



const RNBoolean R4Matrix::
HasTranslation(void) const
{
    // Return whether matrix transformation has translation
    if (RNIsNotZero(m[0][3])) return TRUE;
    if (RNIsNotZero(m[1][3])) return TRUE;
    if (RNIsNotZero(m[2][3])) return TRUE;
    return FALSE;
}



const RNBoolean R4Matrix::
HasScale(void) const
{
    // Return whether matrix transformation has scale
    RNScalar d0 = m[0][0]*m[0][0] + m[1][0]*m[1][0] + m[2][0]*m[2][0];
    if (RNIsNotEqual(d0, 1.0)) return TRUE;
    RNScalar d1 = m[0][1]*m[0][1] + m[1][1]*m[1][1] + m[2][1]*m[2][1];
    if (RNIsNotEqual(d1, 1.0)) return TRUE;
    RNScalar d2 = m[0][2]*m[0][2] + m[1][2]*m[1][2] + m[2][2]*m[2][2];
    if (RNIsNotEqual(d2, 1.0)) return TRUE;
    return FALSE;
}



const RNBoolean R4Matrix::
HasRotation(void) const
{
    // Return whether matrix transformation has rotation
    if (RNIsNotZero(m[0][1])) return TRUE;
    if (RNIsNotZero(m[0][2])) return TRUE;
    if (RNIsNotZero(m[1][0])) return TRUE;
    if (RNIsNotZero(m[1][2])) return TRUE;
    if (RNIsNotZero(m[2][0])) return TRUE;
    if (RNIsNotZero(m[2][1])) return TRUE;
    return FALSE;
}



const RNBoolean R4Matrix::
HasMirror(void) const
{
    // Return whether matrix transformation has mirror operator
    R3Vector vx(m[0][0], m[1][0], m[2][0]);
    R3Vector vy(m[0][1], m[1][1], m[2][1]);
    R3Vector vz(m[0][2], m[1][2], m[2][2]);
    return RNIsNegative(vz.Dot(vx % vy));
}



const RNScalar R4Matrix::
Determinant(void) const
{
    // Return matrix determinant
    return R4MatrixDet4(
        m[0][0], m[0][1], m[0][2], m[0][3],
        m[1][0], m[1][1], m[1][2], m[1][3],
        m[2][0], m[2][1], m[2][2], m[2][3],
        m[3][0], m[3][1], m[3][2], m[3][3]);
}



const R4Matrix R4Matrix::
Transpose(void) const
{
    // Return transpose of matrix
    return R4Matrix(
        m[0][0], m[1][0], m[2][0], m[3][0],
        m[0][1], m[1][1], m[2][1], m[3][1],
        m[0][2], m[1][2], m[2][2], m[3][2],
        m[0][3], m[1][3], m[2][3], m[3][3]);
}



const R4Matrix R4Matrix::
Inverse(void) const
{
    // Return inverse of matrix
    R4Matrix inverse(*this);
    inverse.Invert();
    return inverse;
}



void R4Matrix::
Flip(void)
{
    // Transpose matrix
    RNScalar tmp;
    tmp = m[1][0]; m[1][0] = m[0][1]; m[0][1] = tmp;
    tmp = m[2][0]; m[2][0] = m[0][2]; m[0][2] = tmp;
    tmp = m[3][0]; m[3][0] = m[0][3]; m[0][3] = tmp;
    tmp = m[1][2]; m[1][2] = m[2][1]; m[2][1] = tmp;
    tmp = m[1][3]; m[1][3] = m[3][1]; m[3][1] = tmp;
    tmp = m[2][3]; m[2][3] = m[3][2]; m[3][2] = tmp;
}



void R4Matrix::
Invert(void)
{
    // Copy matrix into local variables
    RNScalar Ma, Mb, Mc, Md, Me, Mf, Mg, Mh, Mi, Mj, Mk, Ml, Mm, Mn, Mo, Mp;
    Ma = m[0][0]; Mb = m[0][1]; Mc = m[0][2]; Md = m[0][3];
    Me = m[1][0]; Mf = m[1][1]; Mg = m[1][2]; Mh = m[1][3];
    Mi = m[2][0]; Mj = m[2][1]; Mk = m[2][2]; Ml = m[2][3];
    Mm = m[3][0]; Mn = m[3][1]; Mo = m[3][2]; Mp = m[3][3];

    // Compute sub-determinants and determinant
    RNScalar a1 = R4MatrixDet3(Mf, Mg, Mh, Mj, Mk, Ml, Mn, Mo, Mp);
    RNScalar a2 = R4MatrixDet3(Me, Mg, Mh, Mi, Mk, Ml, Mm, Mo, Mp);  
    RNScalar a3 = R4MatrixDet3(Me, Mf, Mh, Mi, Mj, Ml, Mm, Mn, Mp);
    RNScalar a4 = R4MatrixDet3(Me, Mf, Mg, Mi, Mj, Mk, Mm, Mn, Mo);
    RNScalar det = Ma*a1 - Mb*a2 + Mc*a3 - Md*a4;
    if (RNIsZero(det, 1E-8)) {
        RNWarning("Unable to invert matrix with zero determinant");
        return;
    }

    // Compute inverse matrix
    m[0][0] = (a1)/det;
    m[1][0] = -(a2)/det;
    m[2][0] = (a3)/det;
    m[3][0] = -(a4)/det;

    m[0][1] = -(R4MatrixDet3(Mb, Mc, Md, Mj, Mk, Ml, Mn, Mo, Mp))/det;
    m[1][1] = (R4MatrixDet3(Ma, Mc, Md, Mi, Mk, Ml, Mm, Mo, Mp))/det;
    m[2][1] = -(R4MatrixDet3(Ma, Mb, Md, Mi, Mj, Ml, Mm, Mn, Mp))/det;
    m[3][1] = (R4MatrixDet3(Ma, Mb, Mc, Mi, Mj, Mk, Mm, Mn, Mo))/det;

    m[0][2] = (R4MatrixDet3(Mb, Mc, Md, Mf, Mg, Mh, Mn, Mo, Mp))/det;
    m[1][2] = -(R4MatrixDet3(Ma, Mc, Md, Me, Mg, Mh, Mm, Mo, Mp))/det;
    m[2][2] = (R4MatrixDet3(Ma, Mb, Md, Me, Mf, Mh, Mm, Mn, Mp))/det;
    m[3][2] = -(R4MatrixDet3(Ma, Mb, Mc, Me, Mf, Mg, Mm, Mn, Mo))/det;

    m[0][3] = -(R4MatrixDet3(Mb, Mc, Md, Mf, Mg, Mh, Mj, Mk, Ml))/det;
    m[1][3] = (R4MatrixDet3(Ma, Mc, Md, Me, Mg, Mh, Mi, Mk, Ml))/det;
    m[2][3] = -(R4MatrixDet3(Ma, Mb, Md, Me, Mf, Mh, Mi, Mj, Ml))/det;
    m[3][3] = (R4MatrixDet3(Ma, Mb, Mc, Me, Mf, Mg, Mi, Mj, Mk))/det;
}



void R4Matrix:: 
XTranslate(RNScalar offset)
{
    // Translate matrix -- post-multiply by: 
    //   [ 1 0 0 tx ]
    //   [ 0 1 0 0  ]
    //   [ 0 0 1 0  ]
    //   [ 0 0 0 1  ]
    m[0][3] += m[0][0] * offset;
    m[1][3] += m[1][0] * offset;
    m[2][3] += m[2][0] * offset;
}



void R4Matrix:: 
YTranslate(RNScalar offset)
{
    // Translate matrix -- post-multiply by: 
    //   [ 1 0 0 0  ]
    //   [ 0 1 0 ty ]
    //   [ 0 0 1 0  ]
    //   [ 0 0 0 1  ]
    m[0][3] += m[0][1] * offset;
    m[1][3] += m[1][1] * offset;
    m[2][3] += m[2][1] * offset;
}



void R4Matrix:: 
ZTranslate(RNScalar offset)
{
    // Translate matrix -- post-multiply by: 
    //   [ 1 0 0 0  ]
    //   [ 0 1 0 0  ]
    //   [ 0 0 1 tz ]
    //   [ 0 0 0 1  ]
    m[0][3] += m[0][2] * offset;
    m[1][3] += m[1][2] * offset;
    m[2][3] += m[2][2] * offset;
}



void R4Matrix:: 
Translate(RNAxis axis, RNScalar offset)
{
    // Translate matrix along axis
    switch (axis) {
    case RN_XAXIS: 
	XTranslate(offset); 
	break;

    case RN_YAXIS: 
	YTranslate(offset); 
	break;

    case RN_ZAXIS: 
	ZTranslate(offset); 
	break;

    default: 
	RNWarning("Matrix translate along undefined axis");
	break;
    }
}



void R4Matrix::
Translate(RNScalar offset)
{
    // Translate matrix
    XTranslate(offset);
    YTranslate(offset);
    ZTranslate(offset);
}



void R4Matrix:: 
Translate(const R3Vector& offset)
{
    // Translate matrix
    XTranslate(offset.X());
    YTranslate(offset.Y());
    ZTranslate(offset.Z());
}



void R4Matrix:: 
XScale(RNScalar scale)
{
    // Scale matrix -- post-multiply by: 
    //   [ sx 0 0 0 ]
    //   [ 0  1 0 0 ]
    //   [ 0  0 1 0 ]
    //   [ 0  0 0 1 ]
    m[0][0] *= scale;
    m[1][0] *= scale;
    m[2][0] *= scale;
}



void R4Matrix:: 
YScale(RNScalar scale)
{
    // Scale matrix -- post-multiply by: 
    //   [ 1 0  0 0 ]
    //   [ 0 sy 0 0 ]
    //   [ 0 0  1 0 ]
    //   [ 0 0  0 1 ]
    m[0][1] *= scale;
    m[1][1] *= scale;
    m[2][1] *= scale;
}



void R4Matrix:: 
ZScale(RNScalar scale)
{
    // Scale matrix -- post-multiply by: 
    //   [ 1 0 0  0 ]
    //   [ 0 1 0  0 ]
    //   [ 0 0 sz 0 ]
    //   [ 0 0 0  1 ]
    m[0][2] *= scale;
    m[1][2] *= scale;
    m[2][2] *= scale;
}



void R4Matrix:: 
Scale(RNAxis axis, RNScalar scale)
{
    // Scale matrix along axis
    switch (axis) {
    case RN_XAXIS: 
	XScale(scale); 
	break;

    case RN_YAXIS: 
	YScale(scale); 
	break;

    case RN_ZAXIS: 
	ZScale(scale); 
	break;

    default: 
	RNWarning("Matrix scale along undefined axis");
	break;
    }
}



void R4Matrix:: 
Scale(RNScalar scale)
{
    // Scale matrix 
    XScale(scale);
    YScale(scale);
    ZScale(scale);
}



void R4Matrix:: 
Scale(const R3Vector& scale)
{
    // Scale matrix
    XScale(scale.X());
    YScale(scale.Y());
    ZScale(scale.Z());
}



void R4Matrix:: 
XRotate(RNAngle radians)
{
    // rotate matrix around X axis counterclockwise
    RNScalar c = cos(radians);
    RNScalar s = sin(radians);
    R4Matrix rotation(
        1.0, 0.0, 0.0, 0.0,
        0.0, c,   -s,  0.0,
        0.0, s,   c,   0.0,
        0.0, 0.0, 0.0, 1.0 );
    *this *= rotation;
}



void R4Matrix:: 
YRotate(RNAngle radians)
{
    // rotate matrix around Y axis counterclockwise
    RNScalar c = cos(radians);
    RNScalar s = sin(radians);
    R4Matrix rotation(
        c,   0.0, s,   0.0,
        0.0, 1.0, 0.0, 0.0,
        -s,  0.0, c,   0.0,
        0.0, 0.0, 0.0, 1.0 );
    *this *= rotation;
}



void R4Matrix:: 
ZRotate(RNAngle radians)
{
    // rotate matrix around Z axis counterclockwise 
    RNScalar c = cos(radians);
    RNScalar s = sin(radians);
    R4Matrix rotation(
        c,   -s,  0.0, 0.0,
        s,   c,   0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0 );
    *this *= rotation;
}



void R4Matrix:: 
Rotate(RNAxis axis, RNAngle radians)
{
    // rotate matrix around an axis counterclockwise
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



void R4Matrix:: 
Rotate(const R3Vector& radians)
{
#if 1
    // Rotate first around X, then around Y, and finally around Z
    R4Matrix rotation = R4identity_matrix;
    rotation.ZRotate(radians.Z());
    rotation.YRotate(radians.Y());
    rotation.XRotate(radians.X());
    *this *= rotation;
#else
    // This is faster way to do same thing (not tested)
    // Note: this matrix is for pre-multiplication -- need one for postmultiplication
    R4Matrix rotation(
      cos(radians.Y())*cos(radians.Z()),  
      sin(radians.X())*sin(radians.Y())*cos(radians.Z()) - cos(radians.X())*sin(radians.Z()), 
      cos(radians.X())*sin(radians.Y())*cos(radians.Z()) + sin(radians.X())*sin(radians.Z()), 0,
      cos(radians.Y())*sin(radians.Z()),
      sin(radians.X())*sin(radians.Y())*sin(radians.Z()) + cos(radians.X())*cos(radians.Z()),
      cos(radians.X())*sin(radians.Y())*sin(radians.Z()) - sin(radians.X())*cos(radians.Z()), 0,
      -sin(radians.Y()), 
      sin(radians.X())*cos(radians.Y()),
      cos(radians.X())*cos(radians.Y()), 0,
      0, 0, 0, 1);
    *this *= rotation;
#endif
}



void R4Matrix:: 
Rotate(const R3Vector& axis, RNAngle radians)
{
    // rotate matrix for arbitrary axis counterclockwise
    // From Graphics Gems I, p. 466
    RNCoord x = axis.X();
    RNCoord y = axis.Y();
    RNCoord z = axis.Z();
    RNScalar c = cos(radians);
    RNScalar s = sin(radians);
    RNScalar t = 1.0 - c;
    R4Matrix rotation(
        t*x*x + c, t*x*y - s*z, t*x*z + s*y, 0.0,
        t*x*y + s*z, t*y*y + c, t*y*z - s*x, 0.0,
        t*x*z - s*y, t*y*z + s*x, t*z*z + c, 0.0,
        0.0, 0.0, 0.0, 1.0);
    *this *= rotation;
}



void R4Matrix:: 
Rotate(const R3Vector& from, const R3Vector& to)
{
    // rotate matrix that takes direction of vector "from" -> "to"
    // This is a quickie hack -- there's got to be a better way
    RNAngle radians = R3InteriorAngle(from, to);
    R3Vector axis = from % to;
    axis.Normalize();
    Rotate(axis, radians);
}



void R4Matrix:: 
Rotate(const R3Quaternion& quaternion)
{
    // Rotate by quaternion
    Multiply(quaternion.Matrix());
}



void R4Matrix:: 
Add(const R4Matrix& a)
{
    // Add matrix entry-by-entry
    m[0][0] += a.m[0][0]; m[0][1] += a.m[0][1]; m[0][2] += a.m[0][2]; m[0][3] += a.m[0][3]; 
    m[1][0] += a.m[1][0]; m[1][1] += a.m[1][1]; m[1][2] += a.m[1][2]; m[1][3] += a.m[1][3]; 
    m[2][0] += a.m[2][0]; m[2][1] += a.m[2][1]; m[2][2] += a.m[2][2]; m[2][3] += a.m[2][3]; 
    m[3][0] += a.m[3][0]; m[3][1] += a.m[3][1]; m[3][2] += a.m[3][2]; m[3][3] += a.m[3][3];
}



void R4Matrix::
Subtract(const R4Matrix& a)
{
    // Subtract matrix entry-by-entry
    m[0][0] -= a.m[0][0]; m[0][1] -= a.m[0][1]; m[0][2] -= a.m[0][2]; m[0][3] -= a.m[0][3]; 
    m[1][0] -= a.m[1][0]; m[1][1] -= a.m[1][1]; m[1][2] -= a.m[1][2]; m[1][3] -= a.m[1][3]; 
    m[2][0] -= a.m[2][0]; m[2][1] -= a.m[2][1]; m[2][2] -= a.m[2][2]; m[2][3] -= a.m[2][3]; 
    m[3][0] -= a.m[3][0]; m[3][1] -= a.m[3][1]; m[3][2] -= a.m[3][2]; m[3][3] -= a.m[3][3];
}



R4Matrix& R4Matrix::
operator=(const R4Matrix& a)
{
    // Assign matrix entry-by-entry
    m[0][0] = a.m[0][0]; m[0][1] = a.m[0][1]; m[0][2] = a.m[0][2]; m[0][3] = a.m[0][3]; 
    m[1][0] = a.m[1][0]; m[1][1] = a.m[1][1]; m[1][2] = a.m[1][2]; m[1][3] = a.m[1][3]; 
    m[2][0] = a.m[2][0]; m[2][1] = a.m[2][1]; m[2][2] = a.m[2][2]; m[2][3] = a.m[2][3]; 
    m[3][0] = a.m[3][0]; m[3][1] = a.m[3][1]; m[3][2] = a.m[3][2]; m[3][3] = a.m[3][3];
    return *this;
}



R4Matrix& R4Matrix::
operator*=(RNScalar a)
{
    // Scale matrix entry-by-entry
    m[0][0] *= a; m[0][1] *= a; m[0][2] *= a; m[0][3] *= a; 
    m[1][0] *= a; m[1][1] *= a; m[1][2] *= a; m[1][3] *= a; 
    m[2][0] *= a; m[2][1] *= a; m[2][2] *= a; m[2][3] *= a; 
    m[3][0] *= a; m[3][1] *= a; m[3][2] *= a; m[3][3] *= a;
    return *this;
}



R4Matrix& R4Matrix::
operator/=(RNScalar a)
{
    // Scale matrix entry-by-entry
    m[0][0] /= a; m[0][1] /= a; m[0][2] /= a; m[0][3] /= a; 
    m[1][0] /= a; m[1][1] /= a; m[1][2] /= a; m[1][3] /= a; 
    m[2][0] /= a; m[2][1] /= a; m[2][2] /= a; m[2][3] /= a; 
    m[3][0] /= a; m[3][1] /= a; m[3][2] /= a; m[3][3] /= a;
    return *this;
}



R4Matrix 
operator-(const R4Matrix& a)
{
    // Negate matrix
    return R4Matrix(
   	-a.m[0][0], -a.m[0][1], -a.m[0][2], -a.m[0][3], 
	-a.m[1][0], -a.m[1][1], -a.m[1][2], -a.m[1][3], 
        -a.m[2][0], -a.m[2][1], -a.m[2][2], -a.m[2][3], 
        -a.m[3][0], -a.m[3][1], -a.m[3][2], -a.m[3][3]);
}



R4Matrix 
operator+(const R4Matrix& a, const R4Matrix& b) 
{
    // Sum matrix
    return R4Matrix(
  	a.m[0][0]+b.m[0][0], a.m[0][1]+b.m[0][1], a.m[0][2]+b.m[0][2], a.m[0][3]+b.m[0][3], 
  	a.m[1][0]+b.m[1][0], a.m[1][1]+b.m[1][1], a.m[1][2]+b.m[1][2], a.m[1][3]+b.m[1][3], 
  	a.m[2][0]+b.m[2][0], a.m[2][1]+b.m[2][1], a.m[2][2]+b.m[2][2], a.m[2][3]+b.m[2][3], 
  	a.m[3][0]+b.m[3][0], a.m[3][1]+b.m[3][1], a.m[3][2]+b.m[3][2], a.m[3][3]+b.m[3][3]);
}



R4Matrix 
operator-(const R4Matrix& a, const R4Matrix& b) 
{
    // Subtract matrix
    return R4Matrix(
  	a.m[0][0]-b.m[0][0], a.m[0][1]-b.m[0][1], a.m[0][2]-b.m[0][2], a.m[0][3]-b.m[0][3], 
  	a.m[1][0]-b.m[1][0], a.m[1][1]-b.m[1][1], a.m[1][2]-b.m[1][2], a.m[1][3]-b.m[1][3], 
  	a.m[2][0]-b.m[2][0], a.m[2][1]-b.m[2][1], a.m[2][2]-b.m[2][2], a.m[2][3]-b.m[2][3], 
  	a.m[3][0]-b.m[3][0], a.m[3][1]-b.m[3][1], a.m[3][2]-b.m[3][2], a.m[3][3]-b.m[3][3]);
}



R4Matrix 
operator*(const R4Matrix& a, RNScalar b)
{
    // Scale matrix
    return R4Matrix(
  	a.m[0][0]*b, a.m[0][1]*b, a.m[0][2]*b, a.m[0][3]*b, 
  	a.m[1][0]*b, a.m[1][1]*b, a.m[1][2]*b, a.m[1][3]*b, 
  	a.m[2][0]*b, a.m[2][1]*b, a.m[2][2]*b, a.m[2][3]*b, 
  	a.m[3][0]*b, a.m[3][1]*b, a.m[3][2]*b, a.m[3][3]*b);
}



inline R4Matrix 
operator/(const R4Matrix& a, RNScalar b) 
{
    // Scale matrix
    assert(b != 0.0);
    return R4Matrix(
  	a.m[0][0]*b, a.m[0][1]*b, a.m[0][2]*b, a.m[0][3]*b, 
  	a.m[1][0]*b, a.m[1][1]*b, a.m[1][2]*b, a.m[1][3]*b, 
  	a.m[2][0]*b, a.m[2][1]*b, a.m[2][2]*b, a.m[2][3]*b, 
  	a.m[3][0]*b, a.m[3][1]*b, a.m[3][2]*b, a.m[3][3]*b);
}



R4Matrix 
operator*(const R4Matrix& a, const R4Matrix& b) 
{
    R4Matrix result;
    int r, c;

    // Multiply matrix
    for (r=0; r<4; r++)
	for (c=0; c<4; c++)
	    result.m[r][c] = a.m[r][0] * b.m[0][c] + a.m[r][1] * b.m[1][c] +
		             a.m[r][2] * b.m[2][c] + a.m[r][3] * b.m[3][c];

    // Return result
    return result;
}



R3Vector
operator*(const R4Matrix& a, const R3Vector& v)
{
    // Multiply matrix by vector
    RNCoord x = a.m[0][0] * v.X() + a.m[0][1] * v.Y() + a.m[0][2] * v.Z();
    RNCoord y = a.m[1][0] * v.X() + a.m[1][1] * v.Y() + a.m[1][2] * v.Z();
    RNCoord z = a.m[2][0] * v.X() + a.m[2][1] * v.Y() + a.m[2][2] * v.Z();
    return R3Vector(x, y, z);
}



R3Point 
operator*(const R4Matrix& a, const R3Point& p)
{
    // Multiply matrix by point
    RNCoord x = a.m[0][0] * p.X() + a.m[0][1] * p.Y() + a.m[0][2] * p.Z() + a.m[0][3];
    RNCoord y = a.m[1][0] * p.X() + a.m[1][1] * p.Y() + a.m[1][2] * p.Z() + a.m[1][3];
    RNCoord z = a.m[2][0] * p.X() + a.m[2][1] * p.Y() + a.m[2][2] * p.Z() + a.m[2][3];
    return R3Point(x, y, z);
}



RNScalar R4MatrixDet2(
    RNScalar a, RNScalar b,
    RNScalar c, RNScalar d)
{
    /* Return determinant of 2x2 matrix */
    return (a * d - b * c);
}



RNScalar R4MatrixDet3 (
    RNScalar a, RNScalar b, RNScalar c, 
    RNScalar d, RNScalar e, RNScalar f, 
    RNScalar g, RNScalar h, RNScalar i)
{
    /* Return determinant of 3x3 matrix */
    return (a * (e * i - h * f) - b * (d * i - g * f) + c * (d * h - g * e));
}



RNScalar R4MatrixDet4 (
    RNScalar a, RNScalar b, RNScalar c, RNScalar d, 
    RNScalar e, RNScalar f, RNScalar g, RNScalar h, 
    RNScalar i, RNScalar j, RNScalar k, RNScalar l, 
    RNScalar m, RNScalar n, RNScalar o, RNScalar p)
{
    /* Return determinant of 4x4 matrix */
    RNScalar det = 0.0;
    det += a * R4MatrixDet3(f,g,h,j,k,l,n,o,p);
    det -= b * R4MatrixDet3(e,g,h,i,k,l,m,o,p);
    det += c * R4MatrixDet3(e,f,h,i,j,l,m,n,p);
    det -= d * R4MatrixDet3(e,f,g,i,j,k,m,n,o);
    return (det);
}




