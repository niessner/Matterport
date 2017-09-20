/* Include file for the GAPS point class */



/* Initialization functions */

int R3InitPoint();
void R3StopPoint();



/* Class definition */

class R3Point {
    public:
        // Constructor/destructor functions
	R3Point(void);
	R3Point(const R3Point& point);
	R3Point(RNCoord x, RNCoord y, RNCoord z);
	R3Point(const RNCoord array[3]);

        // Property functions/operators
	const RNCoord X(void) const;
	const RNCoord Y(void) const;
	const RNCoord Z(void) const;
	const RNCoord Coord(RNDimension dim) const;
	const RNCoord operator[](RNDimension dim) const;
	const RNCoord *Coords(void) const;
	const RNBoolean IsZero(void) const;
	const RNBoolean IsFinite(void) const;
	const R3Vector Vector(void) const;
	const R3Box BBox(void) const;
	const R3Sphere BSphere(void) const;

	// Relationship functions/operators
	const RNBoolean Collinear(const R3Point& point1, const R3Point& point2) const;
	const RNBoolean operator==(const R3Point& point) const;
	const RNBoolean operator!=(const R3Point& point) const;

        // Manipulation functions/operators
	void SetX(RNCoord x);
	void SetY(RNCoord y);
	void SetZ(RNCoord z);
	void SetCoord(RNDimension dim, RNCoord coord);
        void Translate(const R3Vector& vector);
	void Project(const R3Line& line);
	void Project(const R3Plane& plane);
	void Mirror(const R3Plane& plane);
	void XRotate(RNAngle radians);
	void YRotate(RNAngle radians);
	void ZRotate(RNAngle radians);
	void Rotate(const R3Vector& xyz_radians);
	void Rotate(const R3Quaternion& quaternion);
	void Rotate(RNAxis axis, RNAngle radians);
	void Rotate(const R3Vector& axis, RNAngle theta);
	void Rotate(const R3Line& axis, RNAngle theta);
	void Transform(const R3Transformation& transformation);
	void InverseTransform(const R3Transformation& transformation);
	void Reset(RNCoord x, RNCoord y, RNCoord z);

        // Draw functions/operators
        void Draw(void) const;

	// Assignment operators
	R3Point& operator=(const R3Point& point);
	R3Point& operator+=(const R3Point& point);
	R3Point& operator+=(const R3Vector& vector);
	R3Point& operator-=(const R3Vector& vector);
	R3Point& operator*=(const RNScalar a);
	R3Point& operator/=(const RNScalar a);

        // Arithmetic operators
	friend R3Point operator-(const R3Point& point);
	friend R3Point operator+(const R3Point& point1, const R3Point& point2);
	friend R3Point operator+(const R3Point& point, const R3Vector& vector);
	friend R3Point operator+(const R3Vector& vector, const R3Point& point);
	friend R3Vector operator-(const R3Point& point1, const R3Point& point2);
	friend R3Point operator-(const R3Point& point, const R3Vector& vector);
	friend R3Point operator*(const R3Point& point, const RNScalar a);
	friend R3Point operator*(const RNScalar a, const R3Point& point);
	friend R3Point operator/(const R3Point& point, const RNScalar a);

        // Undocumented functions/operators
  	RNCoord& operator[](RNDimension dim);

    private:
	RNCoord v[3];
};



/* Public variables */

extern const R3Point R3null_point;
extern const R3Point R3ones_point;
extern const R3Point R3posx_point;
extern const R3Point R3posy_point;
extern const R3Point R3posz_point;
extern const R3Point R3negx_point;
extern const R3Point R3negy_point;
extern const R3Point R3negz_point;
extern const R3Point R3infinity_point;
extern const R3Point R3unknown_point;
#define R3zero_point R3null_point



/* Inline functions */

inline const RNCoord R3Point::
X (void) const
{
    return(v[0]);
}



inline const RNCoord R3Point::
Y (void) const
{
    return(v[1]);
}



inline const RNCoord R3Point::
Z (void) const
{
    return(v[2]);
}



inline const RNCoord R3Point::
Coord (RNDimension dim) const
{
    assert ((dim>=RN_X)&&(dim<=RN_Z));
    return(v[dim]);
}



inline const RNCoord R3Point::
operator[](RNDimension dim) const
{
    assert((dim>=RN_X) && (dim<=RN_Z));
    return(v[dim]);
}



inline const RNCoord *R3Point::
Coords(void) const
{
    // Return coords array
    return v;
}



inline const RNBoolean R3Point::
IsZero(void) const
{
    // Return whether point is zero
    return ((v[0] == 0.0) && (v[1] == 0.0) && (v[2] == 0.0));
}



inline const R3Vector R3Point::
Vector(void) const
{
    // Return vector to point from origin
    return R3Vector(v[0], v[1], v[2]);
}



inline void R3Point::
SetX (RNCoord x) 
{
    // Set X coord
    v[0] = x;
}



inline void R3Point::
SetY (RNCoord y) 
{
    // Set Y coord
    v[1] = y;
}



inline void R3Point::
SetZ (RNCoord z) 
{
    // Set Z coord
    v[2] = z;
}



inline void R3Point::
SetCoord (RNDimension dim, RNCoord coord) 
{
    // Set coord
    v[dim] = coord;
}



inline void R3Point::
Reset(RNCoord x, RNCoord y, RNCoord z) 
{
    // Set all coords
    v[0] = x;
    v[1] = y;
    v[2] = z;
}



inline void R3Point::
Translate (const R3Vector& vector) 
{
    // Move point by vector
    *this += vector;
}



inline R3Point
operator+(const R3Vector& vector, const R3Point& point)
{
    // Commute vector addition
    return point + vector;
}



inline R3Point
operator*(const RNScalar a, const R3Point& point)
{
    // Commute scale
    return point * a;
}


inline RNCoord& R3Point::
operator[] (RNDimension dim) 
{
    assert((dim>=RN_X) && (dim<=RN_Z));
    return(v[dim]);
}




