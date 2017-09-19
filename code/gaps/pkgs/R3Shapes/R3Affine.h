/* Include file for the R3 affine transformation class */



/* Initialization functions */

int R3InitAffine();
void R3StopAffine();



/* Class definition */

class R3Affine : public R3Transformation {
    public:
        // Constructor functions
	R3Affine(void);
	R3Affine(const R4Matrix& matrix, RNBoolean mirror = 0);
        virtual ~R3Affine(void);

        // Property functions/operators
        const R4Matrix& Matrix(void) const;
        const R4Matrix& InverseMatrix(void) const;
        const RNBoolean IsIdentity(void) const;
        const RNBoolean IsMirrored(void) const;
        const RNBoolean IsAffine(void) const;
        const RNBoolean IsIsotropic(void) const;
        const RNBoolean HasTranslation(void) const;
        const RNBoolean HasScale(void) const;
        const RNBoolean HasRotation(void) const;
        const RNBoolean HasMirror(void) const;
	const R3Affine Inverse(void) const;
	const RNBoolean operator==(const R3Affine& affine) const;
	const RNBoolean operator!=(const R3Affine& affine) const;

	// Application functions/operators
	void Apply(R3Vector& vector) const;
	void Apply(R3Point& point) const;
	void Apply(R3Transformation& transformation) const;
	void Apply(R3Affine& affine) const;
	void ApplyInverse(R3Vector& vector) const;
	void ApplyInverse(R3Point& point) const;
	void ApplyInverse(R3Transformation& transformation) const;
	void ApplyInverse(R3Affine& affine) const;
	void ApplyInverseTranspose(R3Vector& vector) const;
	void ApplyTranspose(R3Vector& vector) const;

	// Manipulation functions/operators
        void Invert(void);
        void Mirror(void);
        void XMirror(void);
        void YMirror(void);
        void ZMirror(void);
        void Mirror(RNAxis axis);
        void Mirror(const R3Plane& plane);
	void XTranslate(RNScalar offset);
	void YTranslate(RNScalar offset);
	void ZTranslate(RNScalar offset);
	void Translate(RNScalar offset);
	void Translate(RNAxis axis, RNScalar offset);
	void Translate(const R3Vector& offset);
	void XScale(RNScalar scale);
	void YScale(RNScalar scale);
	void ZScale(RNScalar scale);
	void Scale(RNScalar scale);
	void Scale(RNAxis axis, RNScalar scale);
	void Scale(const R3Vector& scale);
	void XRotate(RNAngle radians);
	void YRotate(RNAngle radians);
	void ZRotate(RNAngle radians);
	void Rotate(const R3Vector& radians);
	void Rotate(RNAxis axis, RNAngle radians);
	void Rotate(const R3Vector& vector, RNAngle radians);
	void Rotate(const R3Vector& from, const R3Vector& to);
	void Rotate(const R3Quaternion& quaternion);
	void Transform(const R3Transformation& transformation);
	void Transform(const R3Affine& affine);
	void InverseTransform(const R3Transformation& transformation);
	void InverseTransform(const R3Affine& affine);
	void Reset(const R3Transformation& transformation);
	void Reset(const R4Matrix& matrix, RNBoolean mirror = 0);

	// Draw functions/operators
        void Load(void) const;
        void Draw(void) const;
        void Push(void) const;
        void Pop(void) const;

	// Do not use these ???
        const RNScalar ScaleFactor(void) const;

    private:
	// Upkeep functions
	void InvalidateInverse(void);

    private:
	R4Matrix matrix;
	R4Matrix inverse;
        RNFlags flags;
};



/* Public variables */

extern const R3Affine R3null_affine;
extern const R3Affine R3identity_affine;
#define R3null_transformation R3null_affine
#define R3identity_transformation R3identity_affine



/* Affine transformation flags */

#define R3_AFFINE_MIRROR_FLAG           0x00000001
#define R3_AFFINE_INVERSE_UPTODATE_FLAG 0x00000002



/* Inline functions */

inline void R3Affine::
InvalidateInverse(void) 
{
    // Invalidate inverse
    flags.Remove(R3_AFFINE_INVERSE_UPTODATE_FLAG);
}



inline const R4Matrix& R3Affine::
Matrix(void) const
{
    // Return affine transform matrix
    return matrix;
}



inline const RNBoolean R3Affine::
HasTranslation(void) const
{
    // Return whether affine transformation has translation
    return matrix.HasTranslation();
}



inline const RNBoolean R3Affine::
HasScale(void) const
{
    // Return whether affine transformation has scale
    return matrix.HasScale();
}



inline const RNBoolean R3Affine::
HasRotation(void) const
{
    // Return whether affine transformation has rotation
    return matrix.HasRotation();
}



inline const RNBoolean R3Affine::
HasMirror(void) const
{
    // Return whether affine transformation has mirror operator
    return flags[R3_AFFINE_MIRROR_FLAG];
}



inline const R3Affine R3Affine::
Inverse(void) const
{
    // Return inverse affine transform 
    return R3Affine(InverseMatrix(), HasMirror());
}



inline void R3Affine::
Invert(void) 
{
    // Invert matrix
    matrix.Invert();
}



inline void R3Affine::
Mirror(void) 
{
    // Flip mirror flag
    if (HasMirror()) flags.Remove(R3_AFFINE_MIRROR_FLAG);
    else flags.Add(R3_AFFINE_MIRROR_FLAG);
}



inline void R3Affine::
XTranslate(RNScalar offset)
{
    // Translate matrix
    matrix.XTranslate(offset);
    InvalidateInverse();
}



inline void R3Affine::
YTranslate(RNScalar offset)
{
    // Translate matrix
    matrix.YTranslate(offset);
    InvalidateInverse();
}



inline void R3Affine::
ZTranslate(RNScalar offset)
{
    // Translate matrix
    matrix.ZTranslate(offset);
    InvalidateInverse();
}



inline void R3Affine::
Translate(RNScalar offset)
{
    // Translate matrix
    matrix.Translate(offset);
    InvalidateInverse();
}



inline void R3Affine::
Translate(RNAxis axis, RNScalar offset)
{
    // Translate matrix
    matrix.Translate(axis, offset);
    InvalidateInverse();
}



inline void R3Affine::
Translate(const R3Vector& offset)
{
    // Translate matrix
    matrix.Translate(offset);
    InvalidateInverse();
}



inline void R3Affine::
XScale(RNScalar scale)
{
    // Scale matrix
    matrix.XScale(scale);
    InvalidateInverse();
}



inline void R3Affine::
YScale(RNScalar scale)
{
    // Scale matrix
    matrix.YScale(scale);
    InvalidateInverse();
}



inline void R3Affine::
ZScale(RNScalar scale)
{
    // Scale matrix
    matrix.ZScale(scale);
    InvalidateInverse();
}



inline void R3Affine::
Scale(RNScalar scale)
{
    // Scale matrix
    matrix.Scale(scale);
    InvalidateInverse();
}



inline void R3Affine::
Scale(RNAxis axis, RNScalar scale)
{
    // Scale matrix
    matrix.Scale(axis, scale);
    InvalidateInverse();
}



inline void R3Affine::
Scale(const R3Vector& scale)
{
    // Scale matrix
    matrix.Scale(scale);
    InvalidateInverse();
}



inline void R3Affine::
XRotate(RNAngle radians)
{
    // Rotate matrix
    matrix.XRotate(radians);
    InvalidateInverse();
}



inline void R3Affine::
YRotate(RNAngle radians)
{
    // Rotate matrix
    matrix.YRotate(radians);
    InvalidateInverse();
}



inline void R3Affine::
ZRotate(RNAngle radians)
{
    // Rotate matrix
    matrix.ZRotate(radians);
    InvalidateInverse();
}



inline void R3Affine::
Rotate(const R3Vector& radians)
{
    // Rotate matrix
    matrix.Rotate(radians);
    InvalidateInverse();
}



inline void R3Affine::
Rotate(const R3Quaternion& quaternion)
{
    // Rotate matrix
    matrix.Rotate(quaternion);
    InvalidateInverse();
}



inline void R3Affine::
Rotate(RNAxis axis, RNAngle radians)
{
    // Rotate matrix
    matrix.Rotate(axis, radians);
    InvalidateInverse();
}



inline void R3Affine::
Rotate(const R3Vector& vector, RNAngle radians)
{
    // Rotate matrix
    matrix.Rotate(vector, radians);
    InvalidateInverse();
}



inline void R3Affine::
Rotate(const R3Vector& from, const R3Vector& to)
{
    // Rotate matrix
    matrix.Rotate(from, to);
    InvalidateInverse();
}



inline const RNBoolean R3Affine::
operator==(const R3Affine& affine) const
{
    // Return whether or not another affine is the same
    return (Matrix() == affine.Matrix()) && (HasMirror() == affine.HasMirror());
}



inline const RNBoolean R3Affine::
operator!=(const R3Affine& affine) const
{
    // Return whether affine is not equal
    return (!(*this == affine));
}



