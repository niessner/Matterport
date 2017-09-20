/* Include file for the R3 cylinder class */



/* Initialization functions */

int R3InitCylinder();
void R3StopCylinder();



/* Class definition */

class R3Cylinder : public R3Solid {
    public:
        // Constructor functions
	R3Cylinder(void);
        R3Cylinder(const R3Cylinder& cylinder);
        R3Cylinder(const R3Span& axis, RNLength radius);
        R3Cylinder(const R3Point& p1, const R3Point& p2, RNLength radius);

        // Cylinder propetry functions/operators
        const R3Span& Axis(void) const;
        const R3Circle& Base(void) const;
        const R3Circle& Top(void) const;
        const RNLength Radius(void) const;
        const RNLength Height(void) const;
	const RNBoolean IsEmpty(void) const;
	const RNBoolean IsFinite(void) const;

        // Shape propetry functions/operators
	virtual const RNBoolean IsPoint(void) const;
	virtual const RNBoolean IsLinear(void) const ;
	virtual const RNBoolean IsPlanar(void) const ;
	virtual const RNBoolean IsConvex(void) const ;
        virtual const RNInterval NFacets(void) const;
	virtual const RNArea Area(void) const;
	virtual const RNVolume Volume(void) const;
	virtual const R3Point Centroid(void) const;
	virtual const R3Shape& BShape(void) const;
	virtual const R3Box BBox(void) const;
	virtual const R3Sphere BSphere(void) const;

        // Manipulation functions/operators
        virtual void Empty(void);
        virtual void Translate(const R3Vector& vector);
        virtual void Reposition(const R3Span& span);
        virtual void Resize(RNLength radius);
	virtual void Transform(const R3Transformation& transformation);

        // Draw functions/operators
        virtual void Draw(const R3DrawFlags draw_flags = R3_DEFAULT_DRAW_FLAGS) const;

        // Relationship functions/operators
	RNBoolean operator==(const R3Cylinder& cylinder) const;
	RNBoolean operator!=(const R3Cylinder& cylinder) const;

	// Standard shape definitions
	RN_CLASS_TYPE_DECLARATIONS(R3Cylinder);
        R3_SHAPE_RELATIONSHIP_DECLARATIONS(R3Cylinder);

    private:
	R3Span axis;
        R3Circle base;
        R3Circle top;
};



/* Public variables */

extern const R3Cylinder R3null_cylinder;
extern const R3Cylinder R3zero_cylinder;
extern const R3Cylinder R3unit_cylinder;
extern const R3Cylinder R3infinite_cylinder;



/* Inline functions */

inline const R3Span& R3Cylinder::
Axis(void) const
{
    // Return cylinder axis
    return axis;
}



inline const R3Circle& R3Cylinder::
Base(void) const
{
    // Return cylinder base circle
    return base;
}



inline const R3Circle& R3Cylinder::
Top(void) const
{
    // Return cylinder top circle
    return top;
}



inline const RNLength R3Cylinder::
Radius(void) const
{
    // Return cylinder radius
    return base.Radius();
}



inline const RNLength R3Cylinder::
Height(void) const
{
    // Return cylinder radius
    return axis.Length();
}



inline const RNBoolean R3Cylinder::
IsEmpty(void) const
{
    // Return whether cylinder is empty
    return (base.IsEmpty() || top.IsEmpty());
}



inline const RNBoolean R3Cylinder::
IsFinite(void) const
{
    // Return whether cylinder is finite
    return (base.IsFinite() && top.IsFinite());
}



inline RNBoolean R3Cylinder::
operator!=(const R3Cylinder& cylinder) const
{
    // Return whether cylinder is not equal
    return (!(*this == cylinder));
}



