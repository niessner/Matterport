/* Include file for the GAPS halfspace class */



/* Initialization functions */

int R2InitHalfspace();
void R2StopHalfspace();



/* Class definition */

class R2Halfspace {
    public:
        // Constructor functions
	R2Halfspace(void);
        R2Halfspace(const R2Halfspace& halfspace);
	R2Halfspace(RNScalar a, RNScalar b, RNScalar c);
	R2Halfspace(const RNScalar array[3]);
	R2Halfspace(const R2Point& point, const R2Vector& normal);
	R2Halfspace(const R2Point& point1, const R2Point& point2);
	R2Halfspace(RNCoord x1, RNCoord y1, RNCoord x2, RNCoord y2);
	R2Halfspace(const R2Line& line, int dummy);

        // Property functions/operators
	const R2Line& Line(void) const;
	const R2Vector& Normal(void) const;
	const RNBoolean IsZero(void) const;
	const RNBoolean operator==(const R2Halfspace& halfspace) const;
	const RNBoolean operator!=(const R2Halfspace& halfspace) const;

	// Manipulation functions/operators
	void Flip(void);
	void Mirror(const R2Line& line);
	void Translate(const R2Vector& vector);
	void Reposition(const R2Point& point);
	void Align(const R2Vector& vector);
        void Transform(const R2Transformation& transformation);
	void InverseTransform(const R2Transformation& transformation);
	void Reset(const R2Line& line);

        // Draw functions/operators
        void Draw(void) const;

	// Arithmetic functions/operators
	R2Halfspace operator-(void) const;
	
    private:
	R2Line line;
};



/* Public variables */

extern const R2Halfspace R2null_halfspace;
extern const R2Halfspace R2posx_halfspace;
extern const R2Halfspace R2posy_halfspace;
extern const R2Halfspace R2posz_halfspace;
extern const R2Halfspace R2negx_halfspace;
extern const R2Halfspace R2negy_halfspace;
extern const R2Halfspace R2negz_halfspace;



/* Inline functions */

inline const R2Line& R2Halfspace::
Line (void) const
{
    return line;
}



inline const R2Vector& R2Halfspace::
Normal (void) const
{
    return line.Normal();
}



inline const RNBoolean R2Halfspace::
IsZero (void) const
{
    // Return whether halfspace has zero normal vector
    return line.IsZero();
}



inline const RNBoolean R2Halfspace::
operator==(const R2Halfspace& halfspace) const
{
    // Return whether halfspace is equal
    return (line == halfspace.line);
}



inline const RNBoolean R2Halfspace::
operator!=(const R2Halfspace& halfspace) const
{
    // Return whether halfspace is not equal
    return (!(*this == halfspace));
}



inline R2Halfspace R2Halfspace::
operator-(void) const
{
    return R2Halfspace(-line.A(), -line.B(), -line.C());
}



inline void R2Halfspace::
Flip(void)
{
    line.Flip();
}



inline void R2Halfspace::
Translate(const R2Vector& vector) 
{
    line.Translate(vector);
}



inline void R2Halfspace::
Reposition(const R2Point& point)
{
    line.Reposition(point);
}



inline void R2Halfspace::
Align(const R2Vector& vector) 
{
    line.Align(vector);
}



inline void R2Halfspace::
Reset(const R2Line& line)
{
    this->line = line;
}















