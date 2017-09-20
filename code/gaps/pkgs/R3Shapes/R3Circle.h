/* Include file for the R3 circle class */



/* Initialization functions */

int R3InitCircle();
void R3StopCircle();



/* Class definition */

class R3Circle : public R3Surface {
    public:
        // Constructor functions
	R3Circle(void);
        R3Circle(const R3Circle& circle);
        R3Circle(const R3Point& center, RNLength radius, const R3Vector& normal);
        R3Circle(const R3Point& p1, const R3Point& p2, const R3Point& p3);

        // Circle property functions/operators
        const R3Point& Center(void) const;
        const RNLength Radius(void) const;
	const R3Plane& Plane(void) const;
	const R3Vector& Normal(void) const;
	const RNBoolean IsEmpty(void) const;
	const RNBoolean IsFinite(void) const;

        // Shape property functions/operators
	virtual const RNBoolean IsPoint(void) const;
	virtual const RNBoolean IsLinear(void) const;
	virtual const RNBoolean IsPlanar(void) const;
	virtual const RNBoolean IsConvex(void) const;
        virtual const RNInterval NFacets(void) const;
	virtual const RNLength Length(void) const;
	virtual const RNArea Area(void) const;
	virtual const R3Point Centroid(void) const;
	virtual const R3Shape& BShape(void) const;
	virtual const R3Box BBox(void) const;
	virtual const R3Sphere BSphere(void) const;

        // Relationship functions/operators
	RNBoolean operator==(const R3Circle& circle) const;
	RNBoolean operator!=(const R3Circle& circle) const;

        // Manipulation functions/operators
	virtual void Flip(void);
	virtual void Reposition(const R3Point& center);
	virtual void Translate(const R3Vector& offset);
	virtual void Align(const R3Vector& normal);
	virtual void Resize(RNScalar radius);
        virtual void Transform(const R3Transformation& transformation);
	virtual void Reset(const R3Point& center, RNScalar radius, const R3Vector& normal);

        // Draw functions/operators
        virtual void Draw(const R3DrawFlags draw_flags = R3_DEFAULT_DRAW_FLAGS) const;

	// Standard shape definitions
	RN_CLASS_TYPE_DECLARATIONS(R3Circle);
        R3_SHAPE_RELATIONSHIP_DECLARATIONS(R3Circle);

    private:
        R3Point center;
        RNLength radius;
	R3Plane plane;
};



/* Public variables */

extern const R3Circle R3null_circle;
extern const R3Circle R3zero_circle;
extern const R3Circle R3unit_circle;
extern const R3Circle R3infinite_circle;

extern const int R3circle_npoints;
extern RNAngle R3circle_angles[];
extern R3Point R3circle_points[];
extern R2Point R3circle_texcoords[];



/* Inline functions */

inline const R3Point& R3Circle::
Center(void) const
{
    // Return circle center
    return center;
}



inline const RNLength R3Circle::
Radius(void) const
{
    // Return circle radius
    return radius;
}



inline const R3Plane& R3Circle::
Plane(void) const
{
    // Return plane containing circle
    return plane;
}



inline const R3Vector& R3Circle::
Normal(void) const
{
    // Return circle normal
    return plane.Normal();
}



inline const RNBoolean R3Circle::
IsEmpty(void) const
{
    // Return whether the circle is empty
    return (radius < 0.0);
}



inline const RNBoolean R3Circle::
IsFinite(void) const
{
    // Return whether the circle is finite
    return (center.IsFinite() && RNIsFinite(radius));
}



inline RNBoolean R3Circle::
operator==(const R3Circle& circle) const
{
    // Return whether circle is equal
    return ((center == circle.center) && 
	    (radius == circle.radius) && 
	    (plane == circle.plane));
}



inline RNBoolean R3Circle::
operator!=(const R3Circle& circle) const
{
    // Return whether circle is not equal
    return (!(*this == circle));
}



