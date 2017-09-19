/* Include file for the R3 ellipsoid class */



/* Initialization functions */

int R3InitEllipsoid();
void R3StopEllipsoid();



/* Class definition */

class R3Ellipsoid : public R3Solid {
    public:
        // Constructor functions
	R3Ellipsoid(void);
        R3Ellipsoid(const R3Ellipsoid& ellipsoid);
        R3Ellipsoid(const R3CoordSystem& cs, const R3Vector& radii);

        // Ellipsoid property functions/operators
        const R3Point& Center(void) const;
        const R3CoordSystem& CoordSystem(void) const;
        const R3Vector& Radii(void) const;
	const RNBoolean IsEmpty(void) const;
	const RNBoolean IsFinite(void) const;

        // Shape property functions/operators
	virtual const RNBoolean IsPoint(void) const;
	virtual const RNBoolean IsLinear(void) const;
	virtual const RNBoolean IsPlanar(void) const;
	virtual const RNBoolean IsConvex(void) const;
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
        virtual void Reposition(const R3Point& center);
        virtual void Resize(const R3Vector& radii);
	virtual void Transform(const R3Transformation& transformation);

        // Draw functions/operators
        virtual void Draw(const R3DrawFlags draw_flags = R3_DEFAULT_DRAW_FLAGS) const;

        // Relationship functions/operators
	RNBoolean operator==(const R3Ellipsoid& ellipsoid) const;
	RNBoolean operator!=(const R3Ellipsoid& ellipsoid) const;

	// Standard shape definitions
	RN_CLASS_TYPE_DECLARATIONS(R3Ellipsoid);
        R3_SHAPE_RELATIONSHIP_DECLARATIONS(R3Ellipsoid);

        // Upkeep functions
        void UpdateBBox(void);

    private:
        R3CoordSystem cs;
        R3Vector radii;
	R3Box bbox;
};



/* Inline functions */

inline const R3Point& R3Ellipsoid::
Center(void) const
{
    // Return ellipsoid center
    return cs.Origin();
}



inline const R3CoordSystem& R3Ellipsoid::
CoordSystem(void) const
{
    // Return ellipsoid coordinate system
    return cs;
}



inline const R3Vector& R3Ellipsoid::
Radii(void) const
{
    // Return ellipsoid radii
    return radii;
}



inline const RNBoolean R3Ellipsoid::
IsEmpty(void) const
{
    // Return whether ellipsoid is empty
    return ((radii.X() < 0.0) || (radii.Y() < 0.0) || (radii.Z() < 0.0));
}



inline const RNBoolean R3Ellipsoid::
IsFinite(void) const
{
    // Return whether ellipsoid is finite
    return (!IsEmpty() && radii.IsFinite());
}



inline RNBoolean R3Ellipsoid::
operator==(const R3Ellipsoid& ellipsoid) const
{
    // Return whether ellipsoid is equal
    return ((cs == ellipsoid.cs) && (radii == ellipsoid.radii));
}



inline RNBoolean R3Ellipsoid::
operator!=(const R3Ellipsoid& ellipsoid) const
{
    // Return whether ellipsoid is not equal
    return (!(*this == ellipsoid));
}



