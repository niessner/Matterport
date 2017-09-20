/* Include file for the GAPS halfspace class */



/* Initialization functions */

int R3InitHalfspace();
void R3StopHalfspace();



/* Class definition */

class R3Halfspace {
    public:
        // Constructor functions
	R3Halfspace(void);
        R3Halfspace(const R3Halfspace& halfspace);
	R3Halfspace(RNScalar a, RNScalar b, RNScalar c, RNScalar d);
	R3Halfspace(const RNScalar array[4]);
	R3Halfspace(const R3Vector& normal, RNScalar d);
	R3Halfspace(const R3Point& point, const R3Vector& normal);
	R3Halfspace(const R3Point& point, const R3Line& line);
	R3Halfspace(const R3Point& point, const R3Vector& vector1, const R3Vector& vector2);
	R3Halfspace(const R3Point& point1, const R3Point& point2, const R3Point& point3);
	R3Halfspace(const R3Plane& plane, int dummy /* avoid implicit conversion */);

        // Property functions/operators
	const R3Plane& Plane(void) const;
	const R3Vector& Normal(void) const;
	const RNBoolean IsZero(void) const;
	const RNBoolean operator==(const R3Halfspace& halfspace) const;
	const RNBoolean operator!=(const R3Halfspace& halfspace) const;

	// Manipulation functions/operators
	void Flip(void);
	void Mirror(const R3Plane& plane);
	void Translate(const R3Vector& vector);
	void Reposition(const R3Point& point);
	void Align(const R3Vector& vector);
        void Transform(const R3Transformation& transformation);
	void InverseTransform(const R3Transformation& transformation);
	void Reset(const R3Plane& plane);

        // Draw functions/operators
        void Draw(void) const;

	// Arithmetic functions/operators
	R3Halfspace operator-(void) const;
	
    private:
	R3Plane plane;
};



/* Public variables */

extern const R3Halfspace R3null_halfspace;
extern const R3Halfspace R3posx_halfspace;
extern const R3Halfspace R3posy_halfspace;
extern const R3Halfspace R3posz_halfspace;
extern const R3Halfspace R3negx_halfspace;
extern const R3Halfspace R3negy_halfspace;
extern const R3Halfspace R3negz_halfspace;



/* Inline functions */

inline const R3Plane& R3Halfspace::
Plane (void) const
{
    return plane;
}



inline const R3Vector& R3Halfspace::
Normal (void) const
{
    return plane.Normal();
}



inline const RNBoolean R3Halfspace::
IsZero (void) const
{
    // Return whether halfspace has zero normal vector
    return plane.IsZero();
}



inline const RNBoolean R3Halfspace::
operator==(const R3Halfspace& halfspace) const
{
    // Return whether halfspace is equal
    return (plane == halfspace.plane);
}



inline const RNBoolean R3Halfspace::
operator!=(const R3Halfspace& halfspace) const
{
    // Return whether halfspace is not equal
    return (!(*this == halfspace));
}



inline R3Halfspace R3Halfspace::
operator-(void) const
{
    return R3Halfspace(-plane.A(), -plane.B(), -plane.C(), -plane.D());
}



inline void R3Halfspace::
Flip(void)
{
    plane.Flip();
}



inline void R3Halfspace::
Translate(const R3Vector& vector) 
{
    plane.Translate(vector);
}



inline void R3Halfspace::
Reposition(const R3Point& point)
{
    plane.Reposition(point);
}



inline void R3Halfspace::
Align(const R3Vector& normal) 
{
    plane.Align(normal);
}



inline void R3Halfspace::
Reset(const R3Plane& plane)
{
    this->plane = plane;
}



