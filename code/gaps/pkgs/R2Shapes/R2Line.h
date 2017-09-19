/* Include file for the GAPS line class */



/* Initialization functions */

int R2InitLine();
void R2StopLine();



/* Class definition */

class R2Line {
    public:
        // Constructor functions
	R2Line(void);
	R2Line(const R2Line& line);
	R2Line(const RNScalar a, const RNScalar b, const RNScalar c);
	R2Line(const RNScalar array[3]);
	R2Line(const R2Point& point, const R2Vector& vector, RNBoolean normalized = FALSE);
	R2Line(const R2Point& point1, const R2Point& point2);
	R2Line(RNCoord x1, RNCoord y1, RNCoord x2, RNCoord y2);

        // Property functions/operators
	const RNScalar A(void) const;
	const RNScalar B(void) const;
	const RNScalar C(void) const;
	const RNScalar operator[](int i) const;
        const R2Point Point(void) const;
        const R2Vector& Vector(void) const;
        const R2Vector& Normal(void) const;
	const RNBoolean IsZero(void) const;
	const RNBoolean operator==(const R2Line& line) const;
	const RNBoolean operator!=(const R2Line& line) const;

	// Manipulation functions/operators
	void Flip(void);
	void Mirror(const R2Line& line);
	void Project(const R2Line& line);
        void Translate(const R2Vector& vector);
        void Reposition(const R2Point& point);
        void Align(const R2Vector& vector, RNBoolean normalized = FALSE);
        void Transform(const R2Transformation& transformation);
	void InverseTransform(const R2Transformation& transformation);
	void Reset(const R2Point& point, const R2Vector& vector, RNBoolean normalized = FALSE);

        // Draw functions/operators
        void Draw(void) const;

	// Arithmetic functions/operators
	R2Line operator-(void) const;
	
    private:
	R2Vector vector;
	R2Vector normal;
	RNScalar c;
};



/* Public variables */

extern const R2Line R2null_line;
extern const R2Line R2posx_line;
extern const R2Line R2posy_line;
extern const R2Line R2negx_line;
extern const R2Line R2negy_line;
#define R2xaxis_line R2posx_line
#define R2yaxis_line R2posy_line



/* Inline functions */

inline const RNScalar R2Line::
A(void) const
{
    // Return A coefficient of AX+BY+C=0
    return normal.X();
}



inline const RNScalar R2Line::
B(void) const
{
    // Return B coefficient of AX+BY+C=0
    return normal.Y();
}



inline const RNScalar R2Line::
C(void) const
{
    // Return C coefficient of AX+BY+C=0
    return c;
}



inline const RNScalar R2Line::
operator[](int i) const
{
    assert ((i>=0) && (i<=2));
    return ((i == 2) ? c : normal[i]);
}



inline const R2Point R2Line::
Point(void) const
{
    // Return point on line
    return R2zero_point + normal * -c;
}



inline const R2Vector& R2Line::
Vector(void) const
{
    // Return direction vector of line
    return vector;
}



inline const R2Vector& R2Line::
Normal(void) const
{
    // Return normal
    return normal;
}



inline const RNBoolean R2Line::
IsZero (void) const
{
    // Return whether line has zero vector
    return vector.IsZero();
}



inline const RNBoolean R2Line::
operator!=(const R2Line& line) const
{
    // Return whether line is not equal
    return (!(*this == line));
}



inline R2Line R2Line::
operator-(void) const
{
    // Return line with flipped orientation
    return R2Line(-A(), -B(), -C());
}



inline void R2Line::
Flip(void)
{
    // Flip direction of line
    vector.Flip();
    normal.Flip();
    c = -c;
}



inline void R2Line::
Translate(const R2Vector& vector)
{
    // Move line by vector - there's got to be a better way ???
    Reposition(Point() + vector);
}



inline void R2Line::
Reset(const R2Point& point, const R2Vector& vector, RNBoolean normalized)
{
    // Reset line
    Align(vector, normalized);
    Reposition(point);
}





