/* Include file for the R2 transformation class */



/* Initialization functions */

int R2InitTransformation();
void R2StopTransformation();



/* Class definition */

class R2Transformation /* : public R2Base */ {
    public:
        // Destructor function
        virtual ~R2Transformation(void);

        // Property functions/operators
	virtual const RNBoolean IsMirrored(void) const;
	virtual const RNBoolean IsLinear(void) const;
	virtual const RNBoolean IsAffine(void) const;
	virtual const RNBoolean IsIsotropic(void) const;

	// Application functions/operators
	virtual void Apply(R2Vector& vector) const = 0;
	virtual void Apply(R2Point& point) const = 0;
	virtual void Apply(R2Transformation& transformation) const = 0;
	virtual void Apply(R2Affine& affine) const = 0;
	virtual void ApplyInverse(R2Vector& vector) const = 0;
	virtual void ApplyInverse(R2Point& point) const = 0;
	virtual void ApplyInverse(R2Transformation& transformation) const = 0;
	virtual void ApplyInverse(R2Affine& affine) const = 0;

	// Manipulation functions/operators
	virtual void Reset(const R2Transformation& transformation) = 0;
	virtual void Transform(const R2Transformation& transformation) = 0;
	virtual void InverseTransform(const R2Transformation& transformation) = 0;

        // Draw functions/operators
        virtual void Load(void) const = 0;
        virtual void Draw(void) const = 0;
        virtual void Push(void) const = 0;
        virtual void Pop(void) const = 0;
};



