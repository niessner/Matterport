/* Include file for the GAPS triad class */



/* Initialization functions */

int R3InitTriad();
void R3StopTriad();



/* Class definition */

class R3Triad /* : public R3Base */ {
    public:
        // Constructor functions
	R3Triad(void);
	R3Triad(const R3Triad& triad);
	R3Triad(const R3Vector& xaxis, const R3Vector& yaxis, const R3Vector& zaxis);
	R3Triad(const R3Vector& towards, const R3Vector& up);

        // Property functions/operators
        const R4Matrix Matrix(void) const;
        const R4Matrix InverseMatrix(void) const;
	const R3Vector& Axis(RNDimension dim) const;
	const R3Vector& operator[](RNDimension dim) const;
	const RNBoolean operator==(const R3Triad& triad) const;
	const RNBoolean operator!=(const R3Triad& triad) const;

	// Manipulation functions/operators
        void Normalize(void);
	void Rotate(RNAxis axis, RNAngle radians);
	void Rotate(const R3Vector& vector, RNAngle radians);
	void Rotate(const R3Vector& from, const R3Vector& to);
	void Mirror(const R3Plane& plane);
        void Transform(const R3Transformation& transformation);
        void InverseTransform(const R3Transformation& transformation);
	R3Triad& operator=(const R3Triad& triad);

	// Draw functions/operators
        void Draw(void) const;
        void Outline(void) const;

    private:
	R3Vector axis[3];
};



/* Public variables */

extern const R3Triad R3xyz_triad;



/* Inline functions */

inline const R3Vector& R3Triad::
Axis(RNDimension dim) const
{
    // Return axis vector
    assert((dim >= RN_X) && (dim <= RN_Z));
    return axis[dim];
}



inline const R3Vector& R3Triad::
operator[](RNDimension dim) const
{
    // Return axis vector
    return Axis(dim);
}



inline const RNBoolean R3Triad::
operator==(const R3Triad& triad) const
{
    // Return whether triad is equal
    // We can check only two of the three since everything is orthonormal
    return ((axis[0] == triad.axis[0]) && (axis[1] == triad.axis[1]));
}



inline const RNBoolean R3Triad::
operator!=(const R3Triad& triad) const
{
    // Return whether triad is not equal
    return !(*this == triad);
}



inline void R3Triad::
Outline(void) const
{
    // Draw triad
    Draw();
}



inline R3Triad& R3Triad::
operator=(const R3Triad& triad)
{
    axis[0] = triad.axis[0];
    axis[1] = triad.axis[1];
    axis[2] = triad.axis[2];
    return *this;
}



