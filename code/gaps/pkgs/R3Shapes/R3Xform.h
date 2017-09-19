/* Include file for the R3 transformation class */



/* Initialization functions */

int R3InitTransformation();
void R3StopTransformation();



/* Class definition */

class R3Transformation /* : public R3Base */ {
    public:
        // Destructor function
        virtual ~R3Transformation(void);

        // Property functions/operators
        virtual const RNBoolean IsIdentity(void) const;
        virtual const RNBoolean IsMirrored(void) const;
	virtual const RNBoolean IsLinear(void) const;
	virtual const RNBoolean IsAffine(void) const;
	virtual const RNBoolean IsIsotropic(void) const;

	// Application functions/operators
	virtual void Apply(R3Vector& vector) const = 0;
	virtual void Apply(R3Point& point) const = 0;
	virtual void Apply(R3Transformation& transformation) const = 0;
	virtual void Apply(R3Affine& affine) const = 0;
	virtual void ApplyInverse(R3Vector& vector) const = 0;
	virtual void ApplyInverse(R3Point& point) const = 0;
	virtual void ApplyInverse(R3Transformation& transformation) const = 0;
	virtual void ApplyInverse(R3Affine& affine) const = 0;
	virtual void ApplyInverseTranspose(R3Vector& vector) const = 0;
	virtual void ApplyTranspose(R3Vector& vector) const = 0;

	// Manipulation functions/operators
	virtual void Reset(const R3Transformation& transformation) = 0;
	virtual void Transform(const R3Transformation& transformation) = 0;
	virtual void InverseTransform(const R3Transformation& transformation) = 0;

        // Draw functions/operators
        virtual void Load(void) const = 0;
        virtual void Draw(void) const = 0;
        virtual void Push(void) const = 0;
        virtual void Pop(void) const = 0;

	// Do not use these ???
	virtual const RNScalar ScaleFactor(void) const;
};



