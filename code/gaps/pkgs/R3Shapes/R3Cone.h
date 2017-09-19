/* Include file for the R3 cone class */



/* Initialization functions */

int R3InitCone();
void R3StopCone();



/* Class definition */

class R3Cone : public R3Solid {
    public:
        // Constructor functions
	R3Cone(void);
        R3Cone(const R3Cone& cone);
        R3Cone(const R3Span& axis, RNLength radius);
        R3Cone(const R3Point& p1, const R3Point& p2, RNLength radius);

        // Cone propetry functions/operators
        const R3Point& Apex(void) const;
        const R3Span& Axis(void) const;
        const R3Circle& Base(void) const;
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
	RNBoolean operator==(const R3Cone& cone) const;
	RNBoolean operator!=(const R3Cone& cone) const;

	// Standard shape definitions
	RN_CLASS_TYPE_DECLARATIONS(R3Cone);
        R3_SHAPE_RELATIONSHIP_DECLARATIONS(R3Cone);

    private:
	R3Span axis;
	R3Circle base;
};



/* Public variables */

extern const R3Cone R3null_cone;
extern const R3Cone R3zero_cone;
extern const R3Cone R3unit_cone;
extern const R3Cone R3infinite_cone;



/* Inline functions */

inline const R3Point& R3Cone::
Apex(void) const
{
    // Return cone apex
    return axis.End();
}



inline const R3Span& R3Cone::
Axis(void) const
{
    // Return cone axis
    return axis;
}



inline const R3Circle& R3Cone::
Base(void) const
{
    // Return cone base circle
    return base;
}



inline const RNLength R3Cone::
Radius(void) const
{
    // Return cone radius
    return base.Radius();
}



inline const RNLength R3Cone::
Height(void) const
{
    // Return cone height
    return axis.Length();
}



inline const RNBoolean R3Cone::
IsEmpty(void) const
{
    // Return whether cone is empty
    return (base.IsEmpty());
}



inline const RNBoolean R3Cone::
IsFinite(void) const
{
    // Return whether cone is empty
    return (base.IsFinite() && Apex().IsFinite());
}



inline RNBoolean R3Cone::
operator==(const R3Cone& cone) const
{
    // Return whether cone is equal
    return ((base == cone.base) && (axis == cone.axis));
}



inline RNBoolean R3Cone::
operator!=(const R3Cone& cone) const
{
    // Return whether cone is not equal
    return (!(*this == cone));
}



