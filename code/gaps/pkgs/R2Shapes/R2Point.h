/* Include file for the R2 point class */



/* Initialization functions */

int R2InitPoint();
void R2StopPoint();



/* Class definition */

class R2Point {
    public:
        // Constructor/destructor functions
	R2Point(void);
	R2Point(const R2Point& point);
	R2Point(RNCoord x, RNCoord y);
	R2Point(const RNCoord array[2]);

        // Property functions/operators
	const RNCoord X(void) const;
	const RNCoord Y(void) const;
	const RNCoord Coord(RNDimension dim) const;
	const RNCoord operator[](RNDimension dim) const;
	const RNCoord *Coords(void) const;
	const RNBoolean IsZero(void) const;
	const RNBoolean IsFinite(void) const;
	const R2Vector Vector(void) const;

	// Relationship functions/operators
	const RNBoolean Collinear(const R2Point& point1, const R2Point& point2) const;
	const RNBoolean operator==(const R2Point& point) const;
	const RNBoolean operator!=(const R2Point& point) const;

        // Manipulation functions/operators
	void X(RNCoord x);
	void Y(RNCoord y);
	void SetCoord(RNDimension dim, RNCoord coord);
        void Translate(const R2Vector& vector);
	void Project(const R2Line& line);
	void Mirror(const R2Line& line);
	void Rotate(const R2Point& origin, RNAngle theta);
	void Transform(const R2Transformation& transformation);
	void InverseTransform(const R2Transformation& transformation);
	void Reset(RNCoord x, RNCoord y);

        // Draw functions/operators
        void Draw(void) const;

	// Assignment operators
	R2Point& operator=(const R2Point& point);
	R2Point& operator+=(const R2Point& point);
	R2Point& operator+=(const R2Vector& vector);
	R2Point& operator-=(const R2Vector& vector);
	R2Point& operator*=(const RNScalar a);
	R2Point& operator*=(const R3Matrix& m);
	R2Point& operator/=(const RNScalar a);

        // Arithmetic operators
	friend R2Point operator+(const R2Point& point);
	friend R2Point operator-(const R2Point& point);
	friend R2Point operator+(const R2Point& point1, const R2Point& point2) ;
	friend R2Point operator+(const R2Point& point, const R2Vector& vector) ;
	friend R2Point operator+(const R2Vector& vector, const R2Point& point) ;
	friend R2Vector operator-(const R2Point& point1, const R2Point& point2);
	friend R2Point operator-(const R2Point& point, const R2Vector& vector);
	friend R2Point operator*(const R2Point& point, const RNScalar a);
	friend R2Point operator*(const RNScalar a, const R2Point& point);
	friend R2Point operator*(const R2Point& point, const R3Matrix& m);
	friend R2Point operator/(const R2Point& point, const RNScalar a);

        // Undocumented functions/operators
  	RNCoord& operator[](RNDimension dim);

    private:
	RNCoord v[2];
};



/* Public variables */

extern const R2Point R2null_point;
extern const R2Point R2ones_point;
extern const R2Point R2posx_point;
extern const R2Point R2posy_point;
extern const R2Point R2negx_point;
extern const R2Point R2negy_point;
extern const R2Point R2infinite_point;
extern const R2Point R2unknown_point;
#define R2zero_point R2null_point



/* Inline functions */

inline const RNCoord R2Point::
X (void) const
{
    return(v[0]);
}



inline const RNCoord R2Point::
Y (void) const
{
    return(v[1]);
}



inline const RNCoord R2Point::
Coord (RNDimension dim) const
{
    assert ((dim>=RN_X)&&(dim<=RN_Y));
    return(v[dim]);
}



inline const RNCoord R2Point::
operator[](RNDimension dim) const
{
    assert((dim>=RN_X) && (dim<=RN_Y));
    return(v[dim]);
}



inline const RNCoord *R2Point::
Coords(void) const
{
    // Return coords array
    return v;
}



inline const RNBoolean R2Point::
IsZero(void) const
{
    // Return whether point is zero
    return ((v[0] == 0.0) && (v[1] == 0.0));
}



inline const R2Vector R2Point::
Vector(void) const
{
    // Return vector to point from origin
    return R2Vector(v[0], v[1]);
}



inline void R2Point::
X (RNCoord x) 
{
    // Set X coord
    v[0] = x;
}



inline void R2Point::
Y (RNCoord y) 
{
    // Set Y coord
    v[1] = y;
}



inline void R2Point::
SetCoord (RNDimension dim, RNCoord coord) 
{
    // Set coord
    assert ((dim>=RN_X)&&(dim<=RN_Y));
    v[dim] = coord;
}



inline void R2Point::
Reset(RNCoord x, RNCoord y) 
{
    // Set all coords
    v[0] = x;
    v[1] = y;
}



inline void R2Point::
Translate (const R2Vector& vector) 
{
    // Move point by vector
    *this += vector;
}



inline R2Point 
operator+(const R2Vector& vector, const R2Point& point)
{
    // Commute addition
    return (point + vector);
}



inline R2Point 
operator*(const RNScalar a, const R2Point& point)
{
    // Commute scale
    return (point * a);
}



inline RNCoord& R2Point::
operator[] (RNDimension dim) 
{
    assert((dim>=RN_X) && (dim<=RN_Y));
    return(v[dim]);
}



