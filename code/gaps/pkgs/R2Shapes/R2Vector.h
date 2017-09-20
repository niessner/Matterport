/* Include file for the R2 vector class */



/* Initialization functions */

int R2InitVector();
void R2StopVector();



/* Class definition */

class R2Vector {
    public:
        // Constructor functions
	R2Vector(void);
        R2Vector(const R2Vector& vector);
	R2Vector(RNCoord x, RNCoord y);
	R2Vector(const RNCoord array[2]);
	R2Vector(RNAngle angle);

        // Property functions/operators
	const RNCoord X(void) const;
	const RNCoord Y(void) const;
	const RNCoord Coord(RNDimension dim) const;
	const RNCoord operator[](RNDimension dim) const;
	const RNCoord *Coords(void) const;
	const RNBoolean IsZero(void) const;
	const RNBoolean IsFinite(void) const;
	const RNBoolean IsNormalized(void) const;
	const RNLength Length(void) const;
        const RNAngle Angle(void) const;
	const R2Point Point(void) const;
	const RNQuadrant Quadrant(void) const;
	const RNDimension MaxDimension(void) const;

	// Relationship functions/operators
	const RNBoolean operator==(const R2Vector& vector) const;
	const RNBoolean operator!=(const R2Vector& vector) const;
	const RNScalar Dot(const R2Vector& vector) const;
	const RNScalar Cross(const R2Vector& vector) const;

        // Manipulation functions/operators
	void X(RNCoord x);
	void Y(RNCoord y);
	void SetCoord(RNDimension dim, RNCoord coord);
	void Flip(void);
	void Normalize(void);
	void Scale(RNScalar a);
	void Rotate(RNAngle theta);
	void Project(const R2Vector& vector);
	void Mirror(const R2Line& line);
	void Transform(const R2Transformation& transformation);
	void InverseTransform(const R2Transformation& transformation);
	void Reset(RNCoord x, RNCoord y);

        // Draw functions/operators
        void Draw(void) const;

	// Assignment operators
	R2Vector& operator=(const R2Vector& vector);
	R2Vector& operator+=(const R2Vector& vector);
	R2Vector& operator-=(const R2Vector& vector);
	R2Vector& operator*=(const RNScalar a);
        R2Vector& operator*=(const R2Vector& vector);
        R2Vector& operator*=(const R3Matrix& matrix);
	R2Vector& operator/=(const RNScalar a);
	R2Vector& operator/=(const R2Vector& vector);

        // Arithmetic operators
	friend R2Vector operator+(const R2Vector& vector);
	friend R2Vector operator-(const R2Vector& vector);
	friend R2Vector operator+(const R2Vector& vector1, const R2Vector& vector2) ;
	friend R2Vector operator-(const R2Vector& vector1, const R2Vector& vector2);
	friend R2Vector operator*(const R2Vector& vector1, const R2Vector& vector2);
	friend R2Vector operator*(const R2Vector& vector, const RNScalar a);
	friend R2Vector operator*(const RNScalar a, const R2Vector& vector);
	friend R2Vector operator*(const R2Vector& vector, const R3Matrix& m);
	friend R2Vector operator/(const R2Vector& vector1, const R2Vector& vector2);
	friend R2Vector operator/(const R2Vector& vector, const RNScalar a);
	friend RNScalar operator%(const R2Vector& vector1, const R2Vector& vector2);

        // Undocumented functions/operators
  	RNCoord& operator[](RNDimension dim);

    private:
	RNCoord v[2];
};



/* Public variables */

extern const R2Vector R2null_vector;
extern const R2Vector R2ones_vector;
extern const R2Vector R2posx_vector;
extern const R2Vector R2posy_vector;
extern const R2Vector R2negx_vector;
extern const R2Vector R2negy_vector;
extern const R2Vector R2unknown_vector;
#define R2zero_vector R2null_vector
#define R2xaxis_vector R2posx_vector
#define R2yaxis_vector R2posy_vector



/* Inline functions */

inline const RNCoord R2Vector::
X (void) const
{
    return (v[0]);
}



inline const RNCoord R2Vector::
Y (void) const
{
    return (v[1]);
}



inline const RNCoord R2Vector::
Coord (RNDimension dim) const
{
    assert((dim>=RN_X) && (dim<=RN_Y));
    return (v[dim]);
}



inline const RNCoord R2Vector::
operator[](RNDimension dim) const
{
    assert((dim>=RN_X) && (dim<=RN_Y));
    return(v[dim]);
}



inline const RNCoord *R2Vector::
Coords(void) const
{
    // Return coords array
    return v;
}



inline const RNBoolean R2Vector::
IsZero (void) const
{
    // Return whether vector is zero
    return ((v[0] == 0.0) && (v[1] == 0.0));
}



inline const RNBoolean R2Vector::
IsNormalized (void) const
{
    // Return whether vector is normalized
    return (RNIsEqual(this->Length(), 1.0));
}



inline void R2Vector::
X (RNCoord x) 
{
    // Set X coord
    v[0] = x;
}



inline void R2Vector::
Y (RNCoord y) 
{
    // Set Y coord
    v[1] = y;
}



inline void R2Vector::
SetCoord (RNDimension dim, RNCoord coord) 
{
    // Set coord
    assert ((dim>=RN_X)&&(dim<=RN_Y));
    v[dim] = coord;
}



inline void R2Vector::
Reset(RNCoord x, RNCoord y) 
{
    // Set all coords
    v[0] = x;
    v[1] = y;
}



inline R2Vector 
operator*(const RNScalar a, const R2Vector& vector) 
{
    // Commute scaling
    return (vector * a);
}



inline RNCoord& R2Vector::
operator[] (RNDimension dim) 
{
    assert((dim>=RN_X) && (dim<=RN_Y));
    return(v[dim]);
}



