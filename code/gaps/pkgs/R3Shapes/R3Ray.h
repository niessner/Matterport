/* Include file for the GAPS ray class */



/* Initialization functions */

int R3InitRay();
void R3StopRay();



/* Class definition */

class R3Ray {
    public:
        // Constructor functions
	R3Ray(void);
	R3Ray(const R3Ray& ray);
	R3Ray(const R3Point& point, const R3Vector& vector, RNBoolean normalized = FALSE);
	R3Ray(const R3Point& point1, const R3Point& point2);
	R3Ray(RNCoord x1, RNCoord y1, RNCoord z1, RNCoord x2, RNCoord y2, RNCoord z2);

        // Property functions/operators
        const R3Point& Start(void) const;
        const R3Vector& Vector(void) const;
        const R3Line& Line(void) const;
	const R3Point Point(RNScalar t) const;
	const RNScalar T(const R3Point& point) const;
	const RNBoolean IsZero(void) const;
	const RNBoolean operator==(const R3Ray& ray) const;
	const RNBoolean operator!=(const R3Ray& ray) const;

	// Manipulation functions/operators
	void Flip(void);
	void Mirror(const R3Plane& plane);
        void Translate(const R3Vector& vector);
        void Reposition(const R3Point& point);
        void Align(const R3Vector& vector, RNBoolean normalized = FALSE);
        void Transform(const R3Transformation& transformation);
	void InverseTransform(const R3Transformation& transformation);
	void Reset(const R3Point& point, const R3Vector& vector, RNBoolean normalized = FALSE);

        // Draw functions/operators
        void Draw(void) const;

	// Arithmetic functions/operators
	R3Ray operator-(void) const;
	
    private:
	R3Line line;
};



/* Public variables */

extern const R3Ray R3null_ray;
extern const R3Ray R3posx_ray;
extern const R3Ray R3posy_ray;
extern const R3Ray R3posz_ray;
extern const R3Ray R3negx_ray;
extern const R3Ray R3negy_ray;
extern const R3Ray R3negz_ray;
#define R3xaxis_ray R3posx_ray
#define R3yaxis_ray R3posy_ray
#define R3zaxis_ray R3posz_ray



/* Inline functions */

inline const R3Point& R3Ray::
Start(void) const
{
    // Return source point of ray
    return line.Point();
}



inline const R3Vector& R3Ray::
Vector(void) const
{
    // Return direction vector of ray
    return line.Vector();
}



inline const R3Line& R3Ray::
Line(void) const
{
    // Return line containing ray
    return line;
}



inline const RNBoolean R3Ray::
IsZero (void) const
{
    // Return whether ray has zero vector
    return line.IsZero();
}



inline const RNBoolean R3Ray::
operator==(const R3Ray& ray) const
{
    // Return whether ray is equal
    return (line == ray.line);
}



inline const RNBoolean R3Ray::
operator!=(const R3Ray& ray) const
{
    // Return whether ray is not equal
    return (!(*this == ray));
}



inline R3Ray R3Ray::
operator-(void) const
{
    // Return ray with flipped orientation
    return R3Ray(line.Point(), -(line.Vector()));
}



inline void R3Ray::
Flip(void)
{
    // Flip direction of ray
    line.Flip();
}



inline void R3Ray::
Mirror(const R3Plane& plane)
{
    // Mirror ray over plane
    line.Mirror(plane);
}



inline void R3Ray::
Translate(const R3Vector& vector)
{
    // Move endpoint of ray
    line.Translate(vector);
}



inline void R3Ray::
Reposition(const R3Point& point)
{
    // Set endpoint of ray
    line.Reposition(point);
}



inline void R3Ray::
Align(const R3Vector& vector, RNBoolean normalized)
{
    // Set vector of ray
    line.Align(vector, normalized);
}



inline void R3Ray::
Reset(const R3Point& point, const R3Vector& vector, RNBoolean normalized)
{
    // Reset ray
    line.Reset(point, vector, normalized);
}



