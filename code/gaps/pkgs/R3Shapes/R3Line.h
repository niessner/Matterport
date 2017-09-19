/* Include file for the GAPS line class */



/* Initialization functions */

int R3InitLine();
void R3StopLine();



/* Class definition */

class R3Line {
    public:
        // Constructor functions
	R3Line(void);
	R3Line(const R3Line& line);
	R3Line(const R3Point& point, const R3Vector& vector, RNBoolean normalized = FALSE);
	R3Line(const R3Point& point1, const R3Point& point2);
	R3Line(RNCoord x1, RNCoord y1, RNCoord z1, RNCoord x2, RNCoord y2, RNCoord z2);
	R3Line(const RNArray<R3Point *>& points);
	R3Line(R3Point *points, int npoints);

        // Property functions/operators
        const R3Point& Point(void) const;
        const R3Vector& Vector(void) const;
	const RNBoolean IsZero(void) const;
	const RNBoolean operator==(const R3Line& line) const;
	const RNBoolean operator!=(const R3Line& line) const;

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
	R3Line operator-(void) const;
	
    private:
	R3Point point;
        R3Vector vector;
};



/* Public variables */

extern const R3Line R3null_line;
extern const R3Line R3posx_line;
extern const R3Line R3posy_line;
extern const R3Line R3posz_line;
extern const R3Line R3negx_line;
extern const R3Line R3negy_line;
extern const R3Line R3negz_line;
#define R3xaxis_line R3posx_line
#define R3yaxis_line R3posy_line
#define R3zaxis_line R3posz_line



/* Inline functions */

inline const R3Point& R3Line::
Point(void) const
{
    // Return point on line
    return point;
}



inline const R3Vector& R3Line::
Vector(void) const
{
    // Return direction vector of line
    return vector;
}



inline const RNBoolean R3Line::
IsZero (void) const
{
    // Return whether line has zero vector
    return vector.IsZero();
}



inline const RNBoolean R3Line::
operator!=(const R3Line& line) const
{
    // Return whether line is not equal
    return (!(*this == line));
}



inline R3Line R3Line::
operator-(void) const
{
    // Return line with flipped orientation
    return R3Line(point, -vector);
}



inline void R3Line::
Flip(void)
{
    // Flip direction of line
    vector.Flip();
}



inline void R3Line::
Mirror(const R3Plane& plane)
{
    // Mirror line over plane
    point.Mirror(plane);
    vector.Mirror(plane);
}



inline void R3Line::
Translate(const R3Vector& vector)
{
    // Move point on line
    this->point += vector;
}



inline void R3Line::
Reposition(const R3Point& point)
{
    // Set point on line
    this->point = point;
}



inline void R3Line::
Align(const R3Vector& vector, RNBoolean normalized)
{
    // Set vector of line
    this->vector = vector;
    if (!normalized) this->vector.Normalize();
}



inline void R3Line::
Reset(const R3Point& point, const R3Vector& vector, RNBoolean normalized)
{
    // Reset line
    this->point = point;
    this->vector = vector;
    if (!normalized) this->vector.Normalize();
}



