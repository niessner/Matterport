/* Include file for the R2 affine transformation class */



/* Initialization functions */

int R2InitAffine();
void R2StopAffine();



/* Class definition */

class R2Affine : public R2Transformation {
    public:
        // Constructor functions
	R2Affine(void);
	R2Affine(const R3Matrix& matrix, RNBoolean mirror = 0);

        // Property functions/operators
        const R3Matrix& Matrix(void) const;
        const R3Matrix InverseMatrix(void) const;
        const RNBoolean IsMirrored(void) const;
        const RNBoolean IsAffine(void) const;
        const RNBoolean IsIsotropic(void) const;
        const RNBoolean HasTranslation(void) const;
        const RNBoolean HasScale(void) const;
        const RNBoolean HasRotation(void) const;
        const RNBoolean HasMirror(void) const;
	const R2Affine Inverse(void) const;
	const RNBoolean operator==(const R2Affine& affine) const;
	const RNBoolean operator!=(const R2Affine& affine) const;

	// Application functions/operators
	virtual void Apply(R2Vector& vector) const;
	virtual void Apply(R2Point& point) const;
	virtual void Apply(R2Transformation& transformation) const;
	virtual void Apply(R2Affine& affine) const;
	virtual void ApplyInverse(R2Vector& vector) const;
	virtual void ApplyInverse(R2Point& point) const;
	virtual void ApplyInverse(R2Transformation& transformation) const;
	virtual void ApplyInverse(R2Affine& affine) const;

	// Manipulation functions/operators
        void Invert(void);
        void XMirror(void);
        void YMirror(void);
        void Mirror(void);
	void XTranslate(RNScalar offset);
	void YTranslate(RNScalar offset);
	void Translate(RNScalar offset);
	void Translate(RNAxis axis, RNScalar offset);
	void Translate(const R2Vector& offset);
	void XScale(RNScalar scale);
	void YScale(RNScalar scale);
	void Scale(RNScalar scale);
	void Scale(RNAxis axis, RNScalar scale);
	void Scale(const R2Vector& scale);
	void Rotate(RNAngle radians);
	void Transform(const R2Transformation& transformation);
	void Transform(const R2Affine& affine);
	void InverseTransform(const R2Transformation& transformation);
	void InverseTransform(const R2Affine& affine);
	void Reset(const R2Transformation& transformation);
	void Reset(const R3Matrix& matrix, RNBoolean mirror = 0);

	// Draw functions/operators
        virtual void Load(void) const;
        virtual void Draw(void) const;
        virtual void Push(void) const;
        virtual void Pop(void) const;

	// Do not use these ???
        const RNScalar ScaleFactor(void) const;

    private:
	R3Matrix matrix;
        RNBoolean mirror;
};



/* Public variables */

extern const R2Affine R2null_affine;
extern const R2Affine R2identity_affine;
#define R2null_transformation R2null_affine
#define R2identity_transformation R2identity_affine



/* Inline functions */

inline const R3Matrix& R2Affine::
Matrix(void) const
{
    // Return affine transform matrix
    return matrix;
}



inline const R3Matrix R2Affine::
InverseMatrix(void) const
{
    // Return affine inverse transform matrix
    return matrix.Inverse();
}



inline const RNBoolean R2Affine::
HasTranslation(void) const
{
    // Return whether affine transformation has translation
    return matrix.HasTranslation();
}



inline const RNBoolean R2Affine::
HasScale(void) const
{
    // Return whether affine transformation has scale
    return matrix.HasScale();
}



inline const RNBoolean R2Affine::
HasRotation(void) const
{
    // Return whether affine transformation has rotation
    return matrix.HasRotation();
}



inline const RNBoolean R2Affine::
HasMirror(void) const
{
    // Return whether affine transformation has mirror operator
    return mirror;
}



inline const R2Affine R2Affine::
Inverse(void) const
{
    // Return inverse affine transform 
    return R2Affine(InverseMatrix(), HasMirror());
}



inline void R2Affine::
Invert(void) 
{
    // Invert matrix
    matrix.Invert();
}



inline void R2Affine::
XMirror(void) 
{
    // Mirror across YZ plane
    matrix.XScale(-1.0);
    mirror = !mirror;
}



inline void R2Affine::
YMirror(void) 
{
    // Mirror across XZ plane
    matrix.YScale(-1.0);
    mirror = !mirror;
}



inline void R2Affine::
Mirror(void) 
{
    // Flip mirror flag
    mirror = !mirror;
}



inline void R2Affine::
XTranslate(RNScalar offset)
{
    // Translate matrix
    matrix.XTranslate(offset);
}



inline void R2Affine::
YTranslate(RNScalar offset)
{
    // Translate matrix
    matrix.YTranslate(offset);
}



inline void R2Affine::
Translate(RNScalar offset)
{
    // Translate matrix
    matrix.Translate(offset);
}



inline void R2Affine::
Translate(RNAxis axis, RNScalar offset)
{
    // Translate matrix
    matrix.Translate(axis, offset);
}



inline void R2Affine::
Translate(const R2Vector& offset)
{
    // Translate matrix
    matrix.Translate(offset);
}



inline void R2Affine::
XScale(RNScalar scale)
{
    // Scale matrix
    matrix.XScale(scale);
}



inline void R2Affine::
YScale(RNScalar scale)
{
    // Scale matrix
    matrix.YScale(scale);
}



inline void R2Affine::
Scale(RNScalar scale)
{
    // Scale matrix
    matrix.Scale(scale);
}



inline void R2Affine::
Scale(RNAxis axis, RNScalar scale)
{
    // Scale matrix
    matrix.Scale(axis, scale);
}



inline void R2Affine::
Scale(const R2Vector& scale)
{
    // Scale matrix
    matrix.Scale(scale);
}



inline void R2Affine::
Rotate(RNAngle radians)
{
    // Rotate matrix
    matrix.Rotate(radians);
}



inline const RNBoolean R2Affine::
operator==(const R2Affine& affine) const
{
    // Return whether or not another affine is the same
    return (Matrix() == affine.Matrix()) && (HasMirror() == affine.HasMirror());
}



inline const RNBoolean R2Affine::
operator!=(const R2Affine& affine) const
{
    // Return whether affine is not equal
    return (!(*this == affine));
}



