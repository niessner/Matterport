/* Include file for the R3 oriented box class */



/* Initialization functions */

int R3InitOrientedBox();
void R3StopOrientedBox();



/* Class definition */

class R3OrientedBox : public R3Solid {
    public:
        // Constructor functions
	R3OrientedBox(void);
        R3OrientedBox(const R3OrientedBox& box);
        R3OrientedBox(const R3Point& center, const R3Vector& axis0, const R3Vector& axis1, const R3Vector& axis2);
        R3OrientedBox(const R3Point& center, const R3Vector& axis0, const  R3Vector& axis1, RNScalar radius0, RNScalar radius1, RNScalar radius2); 
        R3OrientedBox(const RNArray<R3Point *>& points);

        // OrientedBox propetry functions/operators
        const R3Point& Center(void) const;
        const R3Triad& Axes(void) const;
        const R3CoordSystem& CoordSystem(void) const;
        const RNLength MinRadius(void) const;
        const RNLength MaxRadius(void) const;
        const RNLength DiagonalRadius(void) const;
        const R3Vector& Axis(RNDimension dim) const;
        const RNLength Radius(RNDimension dim) const;
        const R3Point Corner(RNOctant octant) const;
        const R3Point Corner(RNDirection dir0, RNDirection dir1, RNDirection dir2) const;
        const R3Plane Plane(RNSide side) const;
	const R3Plane Plane(RNDirection dir, RNDimension dim) const;
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
	virtual const R3Point ClosestPoint(const R3Point& point) const;
	virtual const R3Point FurthestPoint(const R3Point& point) const;
	virtual const R3Shape& BShape(void) const;
	virtual const R3Box BBox(void) const;
	virtual const R3Sphere BSphere(void) const;

        // Manipulation functions/operators
        virtual void Empty(void);
        virtual void Translate(const R3Vector& vector);
        virtual void Reposition(const R3Point& center);
        virtual void Resize(RNLength radius0, RNLength radius1, RNLength radius2);
        virtual void Reorient(const R3Vector& axis0, const R3Vector& axis1);
	virtual void Transform(const R3Transformation& transformation);
        virtual void Reset(const R3Point& center, const R3Vector& axis0, const  R3Vector& axis1, RNScalar radius0, RNScalar radius1, RNScalar radius2); 

        // Draw functions/operators
        virtual void Draw(const R3DrawFlags draw_flags = R3_DEFAULT_DRAW_FLAGS) const;

        // Relationship functions/operators
	RNBoolean operator==(const R3OrientedBox& box) const;
	RNBoolean operator!=(const R3OrientedBox& box) const;

	// Standard shape definitions
	RN_CLASS_TYPE_DECLARATIONS(R3OrientedBox);
        R3_SHAPE_RELATIONSHIP_DECLARATIONS(R3OrientedBox);

    private:
        R3CoordSystem cs;
        RNLength radius[3];
};



/* Public variables */

extern const R3OrientedBox R3null_oriented_box;
extern const R3OrientedBox R3zero_oriented_box;



/* Inline functions */

inline const R3Point& R3OrientedBox::
Center(void) const
{
    // Return center
    return cs.Origin();
}



inline const R3Triad& R3OrientedBox::
Axes(void) const
{
    // Return triad
    return cs.Axes();
}



inline const R3CoordSystem& R3OrientedBox::
CoordSystem(void) const
{
    // Return coordinate system
    return cs;
}



inline const R3Vector& R3OrientedBox::
Axis(RNDimension dim) const
{
    // Return box axis
    return cs.Axes().Axis(dim);
}



inline const RNLength R3OrientedBox::
Radius(RNDimension dim) const
{
    // Return radius in dimension
    return radius[dim];
}



inline const RNLength R3OrientedBox::
MinRadius(void) const
{
    // Return smallest radius
    RNLength result = Radius(0);
    if (Radius(1) < result) result = Radius(1);
    if (Radius(2) < result) result = Radius(2);
    return result;
}



inline const RNBoolean R3OrientedBox::
IsEmpty(void) const
{
    // Return whether box is empty
    return (MinRadius() < 0.0);
}



inline const RNLength R3OrientedBox::
MaxRadius(void) const
{
    // Check if box is empty
    if (IsEmpty()) return -1.0;

    // Return smallest radius
    RNLength result = Radius(0);
    if (Radius(1) > result) result = Radius(1);
    if (Radius(2) > result) result = Radius(2);
    return result;
}



inline const RNLength R3OrientedBox::
DiagonalRadius(void) const
{
    // Check if box is empty
    if (IsEmpty()) return -1.0;

    // Return box radius
    return sqrt(Radius(0)*Radius(0) + Radius(1)*Radius(1) + Radius(2)*Radius(2));
}



inline const R3Point R3OrientedBox::
Corner(RNDirection dir0, RNDirection dir1, RNDirection dir2) const
{
    // Check if box is empty
    if (IsEmpty()) return R3zero_point;

    // Return corner position
    R3Point corner = Center();
    corner += (dir0 == RN_HI) ? Radius(0)*Axis(0) : -Radius(0)*Axis(0);
    corner += (dir1 == RN_HI) ? Radius(1)*Axis(1) : -Radius(1)*Axis(1);
    corner += (dir2 == RN_HI) ? Radius(2)*Axis(2) : -Radius(2)*Axis(2);
    return corner;
}



inline const R3Point R3OrientedBox::
Corner(RNOctant octant) const
{
    // Check if box is empty
    if (IsEmpty()) return R3zero_point;

    // Return corner position
    R3Point corner = Center();
    corner += (octant & RN_PNN_OCTANT) ? Radius(0)*Axis(0) : -Radius(0)*Axis(0);
    corner += (octant & RN_NPN_OCTANT) ? Radius(1)*Axis(1) : -Radius(1)*Axis(1);
    corner += (octant & RN_NNP_OCTANT) ? Radius(2)*Axis(2) : -Radius(2)*Axis(2);
    return corner;
}



inline const RNBoolean R3OrientedBox::
IsFinite(void) const
{
    // Check if box is empty
    if (IsEmpty()) return FALSE;

    // Return whether box is finite
    return (RNIsFinite(Radius(0)) && RNIsFinite(Radius(1)) && RNIsFinite(Radius(2)));
}



inline RNBoolean R3OrientedBox::
operator!=(const R3OrientedBox& box) const
{
    // Return whether box is not equal
    return (!(*this == box));
}



