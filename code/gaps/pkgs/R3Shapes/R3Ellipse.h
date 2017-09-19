/* Include file for the R3 ellipse class */



/* Initialization functions */

int R3InitEllipse();
void R3StopEllipse();



/* Class definition */

class R3Ellipse : public R3Surface {
    public:
        // Constructor functions
	R3Ellipse(void);
        R3Ellipse(const R3Ellipse& ellipse);
        R3Ellipse(const R3CoordSystem& cs, const R2Vector& radii);

        // Ellipse property functions/operators
        const R3Point& Center(void) const;
        const R3Triad& Axes(void) const;
        const R3Vector& Axis(RNDimension dim) const;
        const R2Vector& Radii(void) const;
	const R3Vector& Normal(void) const;
	const R3Plane Plane(void) const;
        const R3CoordSystem& CoordSystem(void) const;
        const RNLength Radius(int dim) const;
	const RNBoolean IsEmpty(void) const;
	const RNBoolean IsFinite(void) const;

        // Shape property functions/operators
	virtual const RNBoolean IsPoint(void) const;
	virtual const RNBoolean IsLinear(void) const;
	virtual const RNBoolean IsPlanar(void) const;
	virtual const RNBoolean IsConvex(void) const;
        virtual const RNInterval NFacets(void) const;
	virtual const RNArea Area(void) const;
	virtual const R3Point Centroid(void) const;
	virtual const R3Shape& BShape(void) const;
	virtual const R3Box BBox(void) const;
	virtual const R3Sphere BSphere(void) const;

        // Manipulation functions/operators
	virtual void Flip(void);
        virtual void Empty(void);
        virtual void Translate(const R3Vector& vector);
        virtual void Reposition(const R3Point& center);
        virtual void Resize(const R2Vector& radii);
	virtual void Transform(const R3Transformation& transformation);

        // Draw functions/operators
        virtual void Draw(const R3DrawFlags draw_flags = R3_DEFAULT_DRAW_FLAGS) const;

        // Relationship functions/operators
	RNBoolean operator==(const R3Ellipse& ellipse) const;
	RNBoolean operator!=(const R3Ellipse& ellipse) const;

	// Standard shape definitions
	RN_CLASS_TYPE_DECLARATIONS(R3Ellipse);
        R3_SHAPE_RELATIONSHIP_DECLARATIONS(R3Ellipse);

        // Upkeep functions
        void UpdateBBox(void);

    private:
        R3CoordSystem cs;
        R2Vector radii;
	R3Box bbox;
};



/* Inline functions */

inline const R3Point& R3Ellipse::
Center(void) const
{
    // Return ellipse center
    return cs.Origin();
}



inline const R3Triad& R3Ellipse::
Axes(void) const
{
    // Return triad
    return cs.Axes();
}



inline const R3CoordSystem& R3Ellipse::
CoordSystem(void) const
{
    // Return coordinate system
    return cs;
}



inline const R3Vector& R3Ellipse::
Axis(RNDimension dim) const
{
    // Return ellipse axis
    assert((dim >= 0) && (dim <= 1));
    return cs.Axes().Axis(dim);
}



inline const R3Vector& R3Ellipse::
Normal(void) const
{
    // Return normal to plane containing ellipse
    return cs.Axes().Axis(RN_Z);
}



inline const R3Plane R3Ellipse::
Plane(void) const
{
    // Return plane for coordinate system
    return R3Plane(Center(), Normal());
}



inline const R2Vector& R3Ellipse::
Radii(void) const
{
    // Return ellipse radii
    return radii;
}



inline const RNLength R3Ellipse::
Radius(int dim) const
{
    // Return ellipse radius in dimension
    return radii[dim];
}



inline const RNBoolean R3Ellipse::
IsEmpty(void) const
{
    // Return whether ellipse is empty
    return ((radii.X() < 0.0) || (radii.Y() < 0.0));
}



inline const RNBoolean R3Ellipse::
IsFinite(void) const
{
    // Return whether ellipse is finite
    return (!IsEmpty() && radii.IsFinite());
}



inline RNBoolean R3Ellipse::
operator==(const R3Ellipse& ellipse) const
{
    // Return whether ellipse is equal
    return ((cs == ellipse.cs) && (radii == ellipse.radii));
}



inline RNBoolean R3Ellipse::
operator!=(const R3Ellipse& ellipse) const
{
    // Return whether ellipse is not equal
    return (!(*this == ellipse));
}



