/* Include file for the GAPS coordinate system class */



/* Initialization functions */

int R2InitCoordinateSystem();
void R2StopCoordinateSystem();



/* Class definition */

class R2CoordSystem /* : public R2Base */ {
    public:
        // Constructor functions
	R2CoordSystem(void);
	R2CoordSystem(const R2CoordSystem& cs);
	R2CoordSystem(const R2Point& origin, const R2Diad& axes);

        // Property functions/operators
	const R2Point& Origin(void) const;
	const R2Diad& Axes(void) const;
        const R3Matrix Matrix(void) const;
        const R3Matrix InverseMatrix(void) const;
	const RNBoolean operator==(const R2CoordSystem& cs) const;
	const RNBoolean operator!=(const R2CoordSystem& cs) const;

	// Manipulation functions/operators
	void Translate(const R2Vector& offset);
	void Rotate(RNAngle radians);
	void Mirror(const R2Line& line);
        void Transform(const R2Transformation& transformation);
        void InverseTransform(const R2Transformation& transformation);
        void SetOrigin(const R2Point& origin);
        void SetAxes(const R2Diad& axes);

	// Draw functions/operators
        void Draw(void) const;
        void Outline(void) const;

    private:
	R2Point origin;
	R2Diad axes;
};



/* Public variables */

extern const R2CoordSystem R2xy_coordinate_system;



/* Inline functions */

inline const R2Point& R2CoordSystem::
Origin(void) const
{
    // Return origin
    return origin;
}



inline const R2Diad& R2CoordSystem::
Axes(void) const
{
    // Return axes diad
    return axes;
}



inline const RNBoolean R2CoordSystem::
operator==(const R2CoordSystem& cs) const
{
    // Return whether coordinate system is equal
    return (origin == cs.origin) && (axes == cs.axes);
}



inline const RNBoolean R2CoordSystem::
operator!=(const R2CoordSystem& cs) const
{
    // Return whether coordinate system is not equal
    return !(*this == cs);
}



inline void R2CoordSystem::
Outline(void) const
{
    // Draw coordinate system
    Draw();
}



