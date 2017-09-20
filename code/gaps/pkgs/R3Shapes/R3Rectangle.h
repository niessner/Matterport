/* Include file for the R3 rectangle class */



/* Initialization functions */

int R3InitRectangle();
void R3StopRectangle();



/* Class definition */

class R3Rectangle : public R3Surface {
    public:
        // Constructor functions
	R3Rectangle(void);
        R3Rectangle(const R3Rectangle& rectangle);
        R3Rectangle(const R3Point& center, const R3Vector& axis0, const R3Vector& axis1);
        R3Rectangle(const R3Point& center, const R3Vector& axis0, const  R3Vector& axis1, RNScalar radius0, RNScalar radius1); 

        // Rectangle propetry functions/operators
        const R3Point& Center(void) const;
        const R3Triad& Axes(void) const;
        const R3Vector& Axis(RNDimension dim) const;
        const R3Vector& Normal(void) const;
        const R3Plane Plane(void) const;
        const R3CoordSystem& CoordSystem(void) const;
        const RNLength Radius(RNDimension dim) const;
        const RNLength MinRadius(void) const;
        const RNLength MaxRadius(void) const;
        const RNLength DiagonalRadius(void) const;
        const R3Point Corner(RNQuadrant quadrant) const;
        const R3Point Corner(RNDirection dir0, RNDirection dir1) const;
	const RNBoolean IsEmpty(void) const;
	const RNBoolean IsFinite(void) const;
        
        // Shape propetry functions/operators
	virtual const RNBoolean IsPoint(void) const;
	virtual const RNBoolean IsLinear(void) const ;
	virtual const RNBoolean IsPlanar(void) const ;
	virtual const RNBoolean IsConvex(void) const ;
        virtual const RNInterval NFacets(void) const;
	virtual const RNArea Area(void) const;
	virtual const R3Point Centroid(void) const;
	virtual const R3Point ClosestPoint(const R3Point& point) const;
	virtual const R3Point FurthestPoint(const R3Point& point) const;
	virtual const R3Shape& BShape(void) const;
	virtual const R3Box BBox(void) const;
	virtual const R3Sphere BSphere(void) const;

        // Manipulation functions/operators
        virtual void Empty(void);
        virtual void Translate(const R3Vector& vector);
        virtual void Reposition(const R3Point& center);
        virtual void Resize(RNLength radius0, RNLength radius1);
        virtual void Reorient(const R3Vector& axis0, const R3Vector& axis1);
	virtual void Transform(const R3Transformation& transformation);
        virtual void Reset(const R3Point& center, const R3Vector& axis0, const  R3Vector& axis1, RNScalar radius0, RNScalar radius1); 

        // Draw functions/operators
        virtual void Draw(const R3DrawFlags draw_flags = R3_DEFAULT_DRAW_FLAGS) const;

        // Relationship functions/operators
	RNBoolean operator==(const R3Rectangle& rectangle) const;
	RNBoolean operator!=(const R3Rectangle& rectangle) const;

	// Standard shape definitions
	RN_CLASS_TYPE_DECLARATIONS(R3Rectangle);
        R3_SHAPE_RELATIONSHIP_DECLARATIONS(R3Rectangle);

    private:
        R3CoordSystem cs;
        RNLength radius[2];
};



/* Public variables */

extern const R3Rectangle R3null_rectangle;
extern const R3Rectangle R3zero_rectangle;



/* Inline functions */

inline const R3Point& R3Rectangle::
Center(void) const
{
    // Return center
    return cs.Origin();
}



inline const R3Triad& R3Rectangle::
Axes(void) const
{
    // Return triad
    return cs.Axes();
}



inline const R3CoordSystem& R3Rectangle::
CoordSystem(void) const
{
    // Return coordinate system
    return cs;
}



inline const R3Vector& R3Rectangle::
Axis(RNDimension dim) const
{
    // Return rectangle axis
    assert((dim >= 0) && (dim <= 1));
    return cs.Axes().Axis(dim);
}



inline const R3Vector& R3Rectangle::
Normal(void) const
{
  // Return normal to plane containing rectangle
  return cs.Axes().Axis(RN_Z);
}



inline const R3Plane R3Rectangle::
Plane(void) const
{
  // Return plane containing rectangle
  return R3Plane(Center(), Normal());
}



inline const RNLength R3Rectangle::
Radius(RNDimension dim) const
{
    // Return radius in dimension
    assert((dim >= 0) && (dim <= 1));
    return radius[dim];
}



inline const RNBoolean R3Rectangle::
IsEmpty(void) const
{
    // Return whether rectangle is empty
    return (MinRadius() < 0.0);
}



inline const RNLength R3Rectangle::
MinRadius(void) const
{
    // Return smallest radius
    return (Radius(0) < Radius(1)) ? Radius(0) : Radius(1);
}



inline const RNLength R3Rectangle::
MaxRadius(void) const
{
    // Check if rectangle is empty
    if (IsEmpty()) return -1.0;

    // Return smallest radius
    return (Radius(0) > Radius(1)) ? Radius(0) : Radius(1);
}



inline const RNLength R3Rectangle::
DiagonalRadius(void) const
{
    // Check if rectangle is empty
    if (IsEmpty()) return -1.0;

    // Return rectangle radius
    return sqrt(Radius(0)*Radius(0) + Radius(1)*Radius(1) + Radius(2)*Radius(2));
}



inline const R3Point R3Rectangle::
Corner(RNDirection dir0, RNDirection dir1) const
{
    // Check if rectangle is empty
    if (IsEmpty()) return R3zero_point;

    // Return corner position
    R3Point corner = Center();
    corner += (dir0 == RN_HI) ? Radius(0)*Axis(0) : -Radius(0)*Axis(0);
    corner += (dir1 == RN_HI) ? Radius(1)*Axis(1) : -Radius(1)*Axis(1);
    return corner;
}



inline const R3Point R3Rectangle::
Corner(RNQuadrant quadrant) const
{
    // Check if rectangle is empty
    if (IsEmpty()) return R3zero_point;

    // Return corner position
    R3Point corner = Center();
    corner += (quadrant & RN_PN_QUADRANT) ? Radius(0)*Axis(0) : -Radius(0)*Axis(0);
    corner += (quadrant & RN_NP_QUADRANT) ? Radius(1)*Axis(1) : -Radius(1)*Axis(1);
    return corner;
}



inline const RNBoolean R3Rectangle::
IsFinite(void) const
{
    // Check if rectangle is empty
    if (IsEmpty()) return FALSE;

    // Return whether rectangle is finite
    return (RNIsFinite(Radius(0)) && RNIsFinite(Radius(1)));
}



inline RNBoolean R3Rectangle::
operator!=(const R3Rectangle& rectangle) const
{
    // Return whether rectangle is not equal
    return (!(*this == rectangle));
}



