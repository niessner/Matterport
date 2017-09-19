/* Source file for the GAPS matrix class */



/* Include files */

#include "R2Shapes/R2Shapes.h"



/* Public variables */

const R3Matrix R3null_matrix (
    0.0, 0.0, 0.0, 
    0.0, 0.0, 0.0, 
    0.0, 0.0, 0.0
);

const R3Matrix R3identity_matrix (
    1.0, 0.0, 0.0, 
    0.0, 1.0, 0.0, 
    0.0, 0.0, 1.0
);



/* Public functions */

int 
R2InitMatrix()
{
    /* Return success */
    return TRUE;
}



void 
R2StopMatrix()
{
}



R3Matrix::
R3Matrix(void)
{
}



R3Matrix::
R3Matrix(const R3Matrix& matrix)
{
    // Assign matrix entries
    m[0][0] = matrix.m[0][0]; m[0][1] = matrix.m[0][1]; m[0][2] = matrix.m[0][2]; 
    m[1][0] = matrix.m[1][0]; m[1][1] = matrix.m[1][1]; m[1][2] = matrix.m[1][2]; 
    m[2][0] = matrix.m[2][0]; m[2][1] = matrix.m[2][1]; m[2][2] = matrix.m[2][2]; 
}



R3Matrix::
R3Matrix(RNScalar a00, RNScalar a01, RNScalar a02, 
         RNScalar a10, RNScalar a11, RNScalar a12, 
	 RNScalar a20, RNScalar a21, RNScalar a22)
{
    // Assign matrix entries
    m[0][0] = a00; m[0][1] = a01; m[0][2] = a02; 
    m[1][0] = a10; m[1][1] = a11; m[1][2] = a12; 
    m[2][0] = a20; m[2][1] = a21; m[2][2] = a22; 
}



R3Matrix::
R3Matrix(const RNScalar *a)
{
    // Assign matrix entries
    m[0][0] = a[0]; m[0][1] = a[1]; m[0][2] = a[2]; 
    m[1][0] = a[3]; m[1][1] = a[4]; m[1][2] = a[5]; 
    m[2][0] = a[6]; m[2][1] = a[7]; m[2][2] = a[8]; 
}



const int R3Matrix::
IsZero (void) const
{
    // Return whether matrix is zero
    return (*this == R3null_matrix);
}



const int R3Matrix::
IsIdentity (void) const
{
    // Return whether matrix is identity
    return (*this == R3identity_matrix);
}



const RNBoolean R3Matrix::
IsIsotropic(void) const
{
    // Return whether matrix transformation is isotropic
    RNScalar d0 = m[0][0]*m[0][0] + m[1][0]*m[1][0] + m[2][0]*m[2][0];
    RNScalar d1 = m[0][1]*m[0][1] + m[1][1]*m[1][1] + m[2][1]*m[2][1];
    if (d0 != d1) return FALSE;
    return TRUE;
}



const RNBoolean R3Matrix::
HasTranslation(void) const
{
    // Return whether matrix transformation has translation
    if (RNIsNotZero(m[0][2])) return TRUE;
    if (RNIsNotZero(m[1][2])) return TRUE;
    return FALSE;
}



const RNBoolean R3Matrix::
HasScale(void) const
{
    // Return whether matrix transformation has scale
    RNScalar d0 = m[0][0]*m[0][0] + m[1][0]*m[1][0];
    if (RNIsNotEqual(d0, 1.0)) return TRUE;
    RNScalar d1 = m[0][1]*m[0][1] + m[1][1]*m[1][1];
    if (RNIsNotEqual(d1, 1.0)) return TRUE;
    return FALSE;
}



const RNBoolean R3Matrix::
HasRotation(void) const
{
    // Return whether matrix transformation has rotation
    if (RNIsNotZero(m[0][1])) return TRUE;
    if (RNIsNotZero(m[1][0])) return TRUE;
    return FALSE;
}



const RNBoolean R3Matrix::
HasMirror(void) const
{
    // Return whether matrix transformation has mirror operator
    R2Vector vx(m[0][0], m[1][0]);
    R2Vector vy(m[0][1], m[1][1]);
    return RNIsNegative(vx % vy);
}



const RNScalar R3Matrix::
Determinant(void) const
{
    // Return matrix determinant
    return R3MatrixDet3(
        m[0][0], m[0][1], m[0][2], 
        m[1][0], m[1][1], m[1][2], 
        m[2][0], m[2][1], m[2][2]);
}



const R3Matrix R3Matrix::
Transpose(void) const
{
    // Return transpose of matrix
    return R3Matrix(
        m[0][0], m[1][0], m[2][0], 
        m[0][1], m[1][1], m[2][1], 
        m[0][2], m[1][2], m[2][2]);
}



const R3Matrix R3Matrix::
Inverse(void) const
{
    // Return inverse of matrix
    R3Matrix inverse(*this);
    inverse.Invert();
    return inverse;
}



void R3Matrix::
Flip(void)
{
    // Transpose matrix
    RNScalar tmp;
    tmp = m[1][0]; m[1][0] = m[0][1]; m[0][1] = tmp;
    tmp = m[2][0]; m[2][0] = m[0][2]; m[0][2] = tmp;
    tmp = m[1][2]; m[1][2] = m[2][1]; m[2][1] = tmp;
}



void R3Matrix::
Invert(void)
{
    // Compute determinant
    RNScalar det = Determinant();
    if (RNIsZero(det, 1E-8)) {
        RNWarning("Unable to invert matrix with zero determinant");
        return;
    }

    // Copy values
    RNScalar m00 = m[0][0]; RNScalar m01 = m[0][1]; RNScalar m02 = m[0][2];
    RNScalar m10 = m[1][0]; RNScalar m11 = m[1][1]; RNScalar m12 = m[1][2];
    RNScalar m20 = m[2][0]; RNScalar m21 = m[2][1]; RNScalar m22 = m[2][2];

    // Compute inverse matrix
    m[0][0] = R3MatrixDet2(m11, m12, m21, m22)/det;
    m[0][1] = R3MatrixDet2(m02, m01, m22, m21)/det;
    m[0][2] = R3MatrixDet2(m01, m02, m11, m12)/det;
    m[1][0] = R3MatrixDet2(m12, m10, m22, m20)/det;
    m[1][1] = R3MatrixDet2(m00, m02, m20, m22)/det;
    m[1][2] = R3MatrixDet2(m02, m00, m12, m10)/det;
    m[2][0] = R3MatrixDet2(m10, m11, m20, m21)/det;
    m[2][1] = R3MatrixDet2(m01, m00, m21, m20)/det;
    m[2][2] = R3MatrixDet2(m00, m01, m10, m11)/det;
}



void R3Matrix:: 
XTranslate(RNScalar offset)
{
    // Translate matrix -- post-multiply by: 
    //   [ 1 0 tx ]
    //   [ 0 1 0  ]
    //   [ 0 0 1  ]
    m[0][2] += m[0][0] * offset;
    m[1][2] += m[1][0] * offset;
}



void R3Matrix:: 
YTranslate(RNScalar offset)
{
    // Translate matrix -- post-multiply by: 
    //   [ 1 0 0  ]
    //   [ 0 1 ty ]
    //   [ 0 0 1  ]
    m[0][2] += m[0][1] * offset;
    m[1][2] += m[1][1] * offset;
}



void R3Matrix:: 
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

    default: 
	RNWarning("Matrix translate along undefined axis");
	break;
    }
}



void R3Matrix::
Translate(RNScalar offset)
{
    // Translate matrix
    XTranslate(offset);
    YTranslate(offset);
}



void R3Matrix:: 
Translate(const R2Vector& offset)
{
    // Translate matrix
    XTranslate(offset.X());
    YTranslate(offset.Y());
}



void R3Matrix:: 
XScale(RNScalar scale)
{
    // Scale matrix -- post-multiply by: 
    //   [ sx 0 0 ]
    //   [ 0  1 0 ]
    //   [ 0  0 1 ]
    m[0][0] *= scale;
    m[1][0] *= scale;
}



void R3Matrix:: 
YScale(RNScalar scale)
{
    // Scale matrix -- post-multiply by: 
    //   [ 1 0  0 ]
    //   [ 0 sy 0 ]
    //   [ 0 0  1 ]
    m[0][1] *= scale;
    m[1][1] *= scale;
}



void R3Matrix:: 
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

    default: 
	RNWarning("Matrix scale along undefined axis");
	break;
    }
}



void R3Matrix:: 
Scale(RNScalar scale)
{
    // Scale matrix 
    XScale(scale);
    YScale(scale);
}



void R3Matrix:: 
Scale(const R2Vector& scale)
{
    // Scale matrix
    XScale(scale.X());
    YScale(scale.Y());
}



void R3Matrix:: 
Rotate(RNAngle radians)
{
    // rotate matrix counterclockwise 
    RNScalar c = cos(radians);
    RNScalar s = sin(radians);
    R3Matrix rotation(
        c,   -s,  0.0, 
        s,   c,   0.0, 
        0.0, 0.0, 1.0);
    *this *= rotation;
}



void R3Matrix:: 
Add(const R3Matrix& a)
{
    // Add matrix entry-by-entry
    m[0][0] += a.m[0][0]; m[0][1] += a.m[0][1]; m[0][2] += a.m[0][2]; 
    m[1][0] += a.m[1][0]; m[1][1] += a.m[1][1]; m[1][2] += a.m[1][2]; 
    m[2][0] += a.m[2][0]; m[2][1] += a.m[2][1]; m[2][2] += a.m[2][2]; 
}



void R3Matrix::
Subtract(const R3Matrix& a)
{
    // Subtract matrix entry-by-entry
    m[0][0] -= a.m[0][0]; m[0][1] -= a.m[0][1]; m[0][2] -= a.m[0][2]; 
    m[1][0] -= a.m[1][0]; m[1][1] -= a.m[1][1]; m[1][2] -= a.m[1][2]; 
    m[2][0] -= a.m[2][0]; m[2][1] -= a.m[2][1]; m[2][2] -= a.m[2][2]; 
}



R3Matrix& R3Matrix::
operator=(const R3Matrix& a)
{
    // Assign matrix entry-by-entry
    m[0][0] = a.m[0][0]; m[0][1] = a.m[0][1]; m[0][2] = a.m[0][2]; 
    m[1][0] = a.m[1][0]; m[1][1] = a.m[1][1]; m[1][2] = a.m[1][2]; 
    m[2][0] = a.m[2][0]; m[2][1] = a.m[2][1]; m[2][2] = a.m[2][2]; 
    return *this;
}



R3Matrix& R3Matrix::
operator*=(RNScalar a)
{
    // Scale matrix entry-by-entry
    m[0][0] *= a; m[0][1] *= a; m[0][2] *= a; 
    m[1][0] *= a; m[1][1] *= a; m[1][2] *= a; 
    m[2][0] *= a; m[2][1] *= a; m[2][2] *= a; 
    return *this;
}



R3Matrix& R3Matrix::
operator/=(RNScalar a)
{
    // Scale matrix entry-by-entry
    m[0][0] /= a; m[0][1] /= a; m[0][2] /= a; 
    m[1][0] /= a; m[1][1] /= a; m[1][2] /= a; 
    m[2][0] /= a; m[2][1] /= a; m[2][2] /= a; 
    return *this;
}



R3Matrix 
operator-(const R3Matrix& a)
{
    // Negate matrix
    return R3Matrix(
   	-a.m[0][0], -a.m[0][1], -a.m[0][2], 
	-a.m[1][0], -a.m[1][1], -a.m[1][2], 
        -a.m[2][0], -a.m[2][1], -a.m[2][2]);
}



R3Matrix 
operator+(const R3Matrix& a, const R3Matrix& b) 
{
    // Sum matrix
    return R3Matrix(
  	a.m[0][0]+b.m[0][0], a.m[0][1]+b.m[0][1], a.m[0][2]+b.m[0][2], 
  	a.m[1][0]+b.m[1][0], a.m[1][1]+b.m[1][1], a.m[1][2]+b.m[1][2], 
  	a.m[2][0]+b.m[2][0], a.m[2][1]+b.m[2][1], a.m[2][2]+b.m[2][2]);
}



R3Matrix 
operator-(const R3Matrix& a, const R3Matrix& b) 
{
    // Subtract matrix
    return R3Matrix(
  	a.m[0][0]-b.m[0][0], a.m[0][1]-b.m[0][1], a.m[0][2]-b.m[0][2], 
  	a.m[1][0]-b.m[1][0], a.m[1][1]-b.m[1][1], a.m[1][2]-b.m[1][2], 
  	a.m[2][0]-b.m[2][0], a.m[2][1]-b.m[2][1], a.m[2][2]-b.m[2][2]);
}



R3Matrix 
operator*(const R3Matrix& a, RNScalar b)
{
    // Scale matrix
    return R3Matrix(
  	a.m[0][0]*b, a.m[0][1]*b, a.m[0][2]*b, 
  	a.m[1][0]*b, a.m[1][1]*b, a.m[1][2]*b, 
  	a.m[2][0]*b, a.m[2][1]*b, a.m[2][2]*b);
}



inline R3Matrix 
operator/(const R3Matrix& a, RNScalar b) 
{
    // Scale matrix
    assert(b != 0.0);
    return R3Matrix(
  	a.m[0][0]*b, a.m[0][1]*b, a.m[0][2]*b, 
  	a.m[1][0]*b, a.m[1][1]*b, a.m[1][2]*b, 
  	a.m[2][0]*b, a.m[2][1]*b, a.m[2][2]*b);
}



R3Matrix 
operator*(const R3Matrix& a, const R3Matrix& b) 
{
    R3Matrix result;
    int r, c;

    // Multiply matrix
    for (r=0; r<3; r++)
	for (c=0; c<3; c++)
	    result.m[r][c] = a.m[r][0] * b.m[0][c] + 
		             a.m[r][1] * b.m[1][c] +
		             a.m[r][2] * b.m[2][c];

    // Return result
    return result;
}



R2Vector
operator*(const R3Matrix& a, const R2Vector& v)
{
    // Multiply matrix by vector
    RNCoord x = a.m[0][0] * v.X() + a.m[0][1] * v.Y();
    RNCoord y = a.m[1][0] * v.X() + a.m[1][1] * v.Y();
    return R2Vector(x, y);
}



R2Point 
operator*(const R3Matrix& a, const R2Point& p)
{
    // Multiply matrix by point
    RNCoord x = a.m[0][0] * p.X() + a.m[0][1] * p.Y() + a.m[0][2];
    RNCoord y = a.m[1][0] * p.X() + a.m[1][1] * p.Y() + a.m[1][2];
    return R2Point(x, y);
}



RNScalar R3MatrixDet2(
    RNScalar a, RNScalar b,
    RNScalar c, RNScalar d)
{
    /* Return determinant of 2x2 matrix */
    return (a * d - b * c);
}



RNScalar R3MatrixDet3 (
    RNScalar a, RNScalar b, RNScalar c, 
    RNScalar d, RNScalar e, RNScalar f, 
    RNScalar g, RNScalar h, RNScalar i)
{
    /* Return determinant of 3x3 matrix */
    return (a * (e * i - h * f) - b * (d * i - g * f) + c * (d * h - g * e));
}





