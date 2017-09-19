/* Include file for the GAPS ray class */



/* Initialization functions */

int R2InitRay();
void R2StopRay();



/* Class definition */

class R2Ray {
    public:
        // Constructor functions
	R2Ray(void);
	R2Ray(const R2Ray& ray);
	R2Ray(const R2Point& point, const R2Vector& vector, RNBoolean normalized = FALSE);
	R2Ray(const R2Point& point1, const R2Point& point2);
	R2Ray(RNCoord x1, RNCoord y1, RNCoord x2, RNCoord y2);

        // Property functions/operators
        const R2Point& Start(void) const;
        const R2Vector& Vector(void) const;
        const R2Vector& Normal(void) const;
        const R2Line& Line(void) const;
	const R2Point Point(RNScalar t) const;
	const RNScalar T(const R2Point& point) const;
	const RNBoolean IsZero(void) const;
	const RNBoolean operator==(const R2Ray& ray) const;
	const RNBoolean operator!=(const R2Ray& ray) const;

	// Manipulation functions/operators
	void Flip(void);
	void Project(const R2Line& line);
	void Mirror(const R2Line& line);
        void Translate(const R2Vector& vector);
        void Reposition(const R2Point& point);
        void Align(const R2Vector& vector);
        void Transform(const R2Transformation& transformation);
	void InverseTransform(const R2Transformation& transformation);
	void Reset(const R2Point& point, const R2Vector& vector, RNBoolean normalized = FALSE);

        // Draw functions/operators
        void Draw(void) const;

	// Arithmetic functions/operators
	R2Ray operator-(void) const;
	
    private:
	R2Line line;
	R2Point start;
};



/* Public variables */

extern const R2Ray R2null_ray;
extern const R2Ray R2posx_ray;
extern const R2Ray R2posy_ray;
extern const R2Ray R2negx_ray;
extern const R2Ray R2negy_ray;
#define R2xaxis_ray R2posx_ray
#define R2yaxis_ray R2posy_ray



/* Inline functions */

inline const R2Point& R2Ray::
Start(void) const
{
    // Return source point of ray
    return start;
}



inline const R2Vector& R2Ray::
Vector(void) const
{
    // Return direction vector of ray
    return line.Vector();
}



inline const R2Vector& R2Ray::
Normal(void) const
{
    // Return normal vector of ray
    return line.Normal();
}



inline const R2Line& R2Ray::
Line(void) const
{
    // Return line containing ray
    return line;
}



inline const RNBoolean R2Ray::
IsZero (void) const
{
    // Return whether ray has zero vector
    return line.IsZero();
}



inline const RNBoolean R2Ray::
operator==(const R2Ray& ray) const
{
    // Return whether ray is equal
    return ((line == ray.line) && 
	    (start == ray.start));
}



inline const RNBoolean R2Ray::
operator!=(const R2Ray& ray) const
{
    // Return whether ray is not equal
    return (!(*this == ray));
}



inline R2Ray R2Ray::
operator-(void) const
{
    // Return ray with flipped orientation
    return R2Ray(start, -(line.Vector()));
}



inline void R2Ray::
Flip(void)
{
    // Flip direction of ray
    line.Flip();
}



inline void R2Ray::
Translate(const R2Vector& vector)
{
    // Move ray
    line.Translate(vector);
    start.Translate(vector);
}



inline void R2Ray::
Reposition(const R2Point& point)
{
    // Set endpoint of ray
    line.Reposition(point);
    start = point;
}



inline void R2Ray::
Align(const R2Vector& vector)
{
    // Set vector of ray
    line.Align(vector);
}



inline void R2Ray::
Reset(const R2Point& point, const R2Vector& vector, RNBoolean normalized)
{
    // Reset ray
    line.Reset(point, vector, normalized);
    start = point;
}













