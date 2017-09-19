/* Include file for the GAPS vector class */



/* Initialization functions */

int R3InitVector();
void R3StopVector();



/* Class definition */

class R3Vector /* : public R3Base */ {
    public:
        // Constructor functions
	R3Vector(void);
        R3Vector(const R3Vector& vector);
	R3Vector(RNCoord x, RNCoord y, RNCoord z);
	R3Vector(const RNCoord array[3]);
	R3Vector(RNAngle pitch, RNAngle yaw);

        // Property functions/operators
	const RNCoord X(void) const;
	const RNCoord Y(void) const;
	const RNCoord Z(void) const;
	const RNCoord Coord(RNDimension dim) const;
	const RNCoord operator[](RNDimension dim) const;
	const RNCoord *Coords(void) const;
	const RNBoolean IsZero(void) const;
	const RNBoolean IsFinite(void) const;
	const RNBoolean IsNormalized(void) const;
	const RNLength Length(void) const;
	const R3Point Point(void) const;
	const RNSextant Sextant(void) const;
	const RNOctant Octant(void) const;
	const RNDimension MinDimension(void) const;
	const RNDimension MaxDimension(void) const;

	// Relationship functions/operators
	const RNBoolean operator==(const R3Vector& vector) const;
	const RNBoolean operator!=(const R3Vector& vector) const;
	const RNScalar Dot(const R3Vector& vector) const;

        // Manipulation functions/operators
	void SetX(RNCoord x);
	void SetY(RNCoord y);
	void SetZ(RNCoord z);
	void SetCoord(RNDimension dim, RNCoord coord);
	void Flip(void);
	void Normalize(void);
	void Cross(const R3Vector& vector);
	void XRotate(RNAngle radians);
	void YRotate(RNAngle radians);
	void ZRotate(RNAngle radians);
	void Rotate(const R3Vector& xyz_radians);
	void Rotate(const R3Quaternion& quaternion);
	void Rotate(RNAxis axis, RNAngle radians);
	void Rotate(const R3Vector& axis, RNAngle theta);
	void Project(const R3Vector& vector);
	void Project(const R3Plane& plane);
	void Mirror(const R3Plane& plane);
	void Transform(const R3Transformation& transformation);
	void InverseTransform(const R3Transformation& transformation);
	void Reset(RNCoord x, RNCoord y, RNCoord z);

        // Draw functions/operators
        void Draw(void) const;
        void Outline(void) const;

	// Assignment operators
	R3Vector& operator=(const R3Vector& vector);
	R3Vector& operator+=(const R3Vector& vector);
	R3Vector& operator-=(const R3Vector& vector);
	R3Vector& operator*=(const RNScalar a);
        R3Vector& operator*=(const R3Vector& vector);
	R3Vector& operator/=(const RNScalar a);
	R3Vector& operator/=(const R3Vector& vector);

        // Arithmetic operators
	friend R3Vector operator+(const R3Vector& vector);
	friend R3Vector operator-(const R3Vector& vector);
	friend R3Vector operator+(const R3Vector& vector1, const R3Vector& vector2);
	friend R3Vector operator-(const R3Vector& vector1, const R3Vector& vector2);
	friend R3Vector operator*(const R3Vector& vector1, const R3Vector& vector2);
	friend R3Vector operator*(const R3Vector& vector, const RNScalar a);
	friend R3Vector operator*(const RNScalar a, const R3Vector& vector);
	friend R3Vector operator/(const R3Vector& vector1, const R3Vector& vector2);
	friend R3Vector operator/(const R3Vector& vector, const RNScalar a);
	friend R3Vector operator%(const R3Vector& vector1, const R3Vector& vector2);

        // Undocumented functions/operators
  	RNCoord& operator[](RNDimension dim);

    private:
	RNCoord v[3];
};



/* Public variables */

extern const R3Vector R3null_vector;
extern const R3Vector R3ones_vector;
extern const R3Vector R3posx_vector;
extern const R3Vector R3posy_vector;
extern const R3Vector R3posz_vector;
extern const R3Vector R3negx_vector;
extern const R3Vector R3negy_vector;
extern const R3Vector R3negz_vector;
extern const R3Vector R3unknown_vector;
#define R3zero_vector R3null_vector
#define R3xaxis_vector R3posx_vector
#define R3yaxis_vector R3posy_vector
#define R3zaxis_vector R3posz_vector



/* Public function declarations */

R3Vector R3RandomDirection(void);



/* Inline functions */

inline const RNCoord R3Vector::
X (void) const
{
    return (v[0]);
}



inline const RNCoord R3Vector::
Y (void) const
{
    return (v[1]);
}



inline const RNCoord R3Vector::
Z (void) const
{
    return (v[2]);
}



inline const RNCoord R3Vector::
Coord (RNDimension dim) const
{
    assert((dim>=RN_X) && (dim<=RN_Z));
    return (v[dim]);
}



inline const RNCoord R3Vector::
operator[](RNDimension dim) const
{
    assert((dim>=RN_X) && (dim<=RN_Z));
    return(v[dim]);
}



inline const RNCoord *R3Vector::
Coords(void) const
{
    // Return coords array
    return v;
}



inline const RNBoolean R3Vector::
IsZero (void) const
{
    // Return whether vector is zero
    return ((v[0] == 0.0) && (v[1] == 0.0) && (v[2] == 0.0));
}



inline const RNBoolean R3Vector::
IsNormalized (void) const
{
    // Return whether vector is normalized
    return (RNIsEqual(this->Length(), 1.0));
}



inline void R3Vector::
SetX (RNCoord x) 
{
    // Set X coord
    v[0] = x;
}



inline void R3Vector::
SetY (RNCoord y) 
{
    // Set Y coord
    v[1] = y;
}



inline void R3Vector::
SetZ (RNCoord z) 
{
    // Set Z coord
    v[2] = z;
}



inline void R3Vector::
SetCoord (RNDimension dim, RNCoord coord) 
{
    // Set coord
    v[dim] = coord;
}



inline void R3Vector::
Reset(RNCoord x, RNCoord y, RNCoord z) 
{
    // Set all coords
    v[0] = x;
    v[1] = y;
    v[2] = z;
}



inline void R3Vector::
Outline (void) const
{
    // Draw vector
    Draw();
}



inline R3Vector
operator*(const RNScalar a, const R3Vector& vector)
{
    return vector * a;
}



inline RNCoord& R3Vector::
operator[] (RNDimension dim) 
{
    assert((dim>=RN_X) && (dim<=RN_Z));
    return(v[dim]);
}



