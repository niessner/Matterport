/* Include file for the GAPS diad class */



/* Initialization functions */

int R2InitDiad();
void R2StopDiad();



/* Class definition */

class R2Diad /* : public R2Base */ {
    public:
        // Constructor functions
	R2Diad(void);
	R2Diad(const R2Diad& diad);
	R2Diad(const R2Vector& xaxis, const R2Vector& yaxis);

        // Property functions/operators
        const R3Matrix Matrix(void) const;
        const R3Matrix InverseMatrix(void) const;
	const R2Vector& Axis(RNDimension dim) const;
	const R2Vector& operator[](RNDimension dim) const;
	const RNBoolean operator==(const R2Diad& diad) const;
	const RNBoolean operator!=(const R2Diad& diad) const;

	// Manipulation functions/operators
        void Normalize(void);
	void Rotate(RNAngle radians);
	void Mirror(const R2Line& line);
        void Transform(const R2Transformation& transformation);
        void InverseTransform(const R2Transformation& transformation);
	R2Diad& operator=(const R2Diad& diad);

	// Draw functions/operators
        void Draw(void) const;
        void Outline(void) const;

    private:
	R2Vector axis[2];
};



/* Public variables */

extern const R2Diad R2xy_diad;



/* Inline functions */

inline const R2Vector& R2Diad::
Axis(RNDimension dim) const
{
    // Return axis vector
    assert((dim >= RN_X) && (dim <= RN_Y));
    return axis[dim];
}



inline const R2Vector& R2Diad::
operator[](RNDimension dim) const
{
    // Return axis vector
    return Axis(dim);
}



inline const RNBoolean R2Diad::
operator==(const R2Diad& diad) const
{
    // Return whether diad is equal
    // We can check only one of the three since everything is orthonormal
  return (axis[0] == diad.axis[0]);
}



inline const RNBoolean R2Diad::
operator!=(const R2Diad& diad) const
{
    // Return whether diad is not equal
    return !(*this == diad);
}



inline void R2Diad::
Outline(void) const
{
    // Draw diad
    Draw();
}



inline R2Diad& R2Diad::
operator=(const R2Diad& diad)
{
    axis[0] = diad.axis[0];
    axis[1] = diad.axis[1];
    return *this;
}



