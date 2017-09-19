/* Include file for the GAPS coordinate system class */



/* Initialization functions */

int R3InitCoordinateSystem();
void R3StopCoordinateSystem();



/* Class definition */

class R3CoordSystem /* : public R3Base */ {
    public:
        // Constructor functions
	R3CoordSystem(void);
	R3CoordSystem(const R3CoordSystem& cs);
	R3CoordSystem(const R3Point& origin, const R3Triad& axes);

        // Property functions/operators
	const R3Point& Origin(void) const;
	const R3Triad& Axes(void) const;
        const R4Matrix Matrix(void) const;
        const R4Matrix InverseMatrix(void) const;
	const RNBoolean operator==(const R3CoordSystem& cs) const;
	const RNBoolean operator!=(const R3CoordSystem& cs) const;

	// Manipulation functions/operators
	void Translate(const R3Vector& offset);
	void Rotate(RNAxis axis, RNAngle radians);
	void Rotate(const R3Vector& axis, RNAngle radians);
	void Rotate(const R3Vector& from, const R3Vector& to);
	void Mirror(const R3Plane& plane);
        void Transform(const R3Transformation& transformation);
        void InverseTransform(const R3Transformation& transformation);
        void SetOrigin(const R3Point& origin);
        void SetAxes(const R3Triad& axes);
        void Reset(const R3Point& origin, const R3Triad& axes);

	// Draw functions/operators
        void Draw(void) const;
        void Outline(void) const;

    private:
	R3Point origin;
	R3Triad axes;
};



/* Public variables */

extern const R3CoordSystem R3xyz_coordinate_system;



/* Inline functions */

inline const R3Point& R3CoordSystem::
Origin(void) const
{
    // Return origin
    return origin;
}



inline const R3Triad& R3CoordSystem::
Axes(void) const
{
    // Return axes triad
    return axes;
}



inline const RNBoolean R3CoordSystem::
operator==(const R3CoordSystem& cs) const
{
    // Return whether coordinate system is equal
    return (origin == cs.origin) && (axes == cs.axes);
}



inline const RNBoolean R3CoordSystem::
operator!=(const R3CoordSystem& cs) const
{
    // Return whether coordinate system is not equal
    return !(*this == cs);
}



inline void R3CoordSystem::
Outline(void) const
{
    // Draw coordinate system
    Draw();
}



