/* Include file for the GAPS plane class */



/* Initialization functions */

int R3InitPlane();
void R3StopPlane();



/* Class definition */

class R3Plane {
    public:
        // Constructor functions
	R3Plane(void);
	R3Plane(const R3Plane& plane);
	R3Plane(RNScalar a, RNScalar b, RNScalar c, RNScalar d);
	R3Plane(const RNScalar array[4]);
	R3Plane(const R3Vector& normal, RNScalar d);
	R3Plane(const R3Point& point, const R3Vector& normal);
	R3Plane(const R3Point& point, const R3Line& line);
	R3Plane(const R3Point& point, const R3Vector& vector1, const R3Vector& vector2);
	R3Plane(const R3Point& point1, const R3Point& point2, const R3Point& point3);
	R3Plane(const RNArray<R3Point *>& points, RNBoolean polygon_vertices = TRUE);
	R3Plane(R3Point *points, int npoints, RNBoolean polygon_vertices = TRUE);

        // Property functions/operators
	const RNScalar A(void) const;
	const RNScalar B(void) const;
	const RNScalar C(void) const;
	const RNScalar D(void) const;
	const RNScalar operator[](int i) const;
	const R3Point Point(void) const;
	const R3Vector& Normal(void) const;
	const RNBoolean IsZero(void) const;
	const RNBoolean operator==(const R3Plane& plane) const;
	const RNBoolean operator!=(const R3Plane& plane) const;

	// Manipulation functions/operators
	void Flip(void);
	void Mirror(const R3Plane& plane);
	void Translate(const R3Vector& vector);
	void Reposition(const R3Point& point);
	void Align(const R3Vector& normal);
        void Transform(const R3Transformation& transformation);
	void InverseTransform(const R3Transformation& transformation);
	void Reset(const R3Point& point, const R3Vector& normal);

        // Draw functions/operators
        void Draw(void) const;

	// Arithmetic functions/operators
	R3Plane operator-(void) const;
	
	// Undocumented functions/operators
  	RNScalar& operator[](int i);

    private:
	R3Vector v;
	RNScalar d;
};



/* Public variables */

extern const R3Plane R3null_plane;
extern const R3Plane R3posxz_plane;
extern const R3Plane R3posxy_plane;
extern const R3Plane R3posyz_plane;
extern const R3Plane R3negxz_plane;
extern const R3Plane R3negxy_plane;
extern const R3Plane R3negyz_plane;
#define R3xz_plane R3posxz_plane
#define R3xy_plane R3posxy_plane
#define R3yz_plane R3posyz_plane



/* Inline functions */

inline const RNScalar R3Plane::
A (void) const
{
    return v[0];
}



inline const RNScalar R3Plane::
B (void) const
{
    return v[1];
}



inline const RNScalar R3Plane::
C (void) const
{
    return v[2];
}



inline const RNScalar R3Plane::
D (void) const
{
    return d;
}



inline const RNScalar R3Plane::
operator[](int i) const
{
    assert ((i>=0) && (i<=3));
    return ((i == 3) ? d : v[i]);
}



inline const R3Vector& R3Plane::
Normal (void) const
{
    return v;
}



inline const RNBoolean R3Plane::
IsZero (void) const
{
    // Return whether plane has zero normal vector
    return v.IsZero();
}



inline const RNBoolean R3Plane::
operator==(const R3Plane& plane) const
{
    // Return whether plane is equal
    return ((v == plane.v) && (d == plane.d));
}



inline const RNBoolean R3Plane::
operator!=(const R3Plane& plane) const
{
    // Return whether plane is not equal
    return (!(*this == plane));
}



inline R3Plane R3Plane::
operator-(void) const
{
    // Return plane with flipped orientation
    return R3Plane(-v, -d);
}



inline void R3Plane::
Flip(void) 
{
    v = -v;
    d = -d;
}



inline void R3Plane::
Align(const R3Vector& normal) 
{
    // Align plane normal - keep same distance to origin
    v = normal;
}



inline RNScalar& R3Plane::
operator[](int i) 
{
    assert ((i>=0) && (i<=3));
    return ((i == 3) ? d : v[i]);
}




