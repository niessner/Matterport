/* Include file for the GAPS span class */



/* Initialization functions */

int R2InitSpan();
void R2StopSpan();



/* Class definition */

class R2Span {
    public:
        // Constructor functions
	R2Span(void);
	R2Span(const R2Span& span);
	R2Span(const R2Point& point, const R2Vector& vector);
	R2Span(const R2Point& point1, const R2Point& point2);
	R2Span(RNCoord x1, RNCoord y1, RNCoord x2, RNCoord y2);

        // Property functions/operators
        const R2Point& Start(void) const;
        const R2Point& End(void) const;
        const R2Vector& Vector(void) const;
        const R2Vector& Normal(void) const;
        const R2Point Point(int k) const;
	const R2Point Point(RNScalar t) const;
        const R2Point& operator[](int k) const;
        const R2Ray& Ray(void) const;
        const R2Line& Line(void) const;
        const R2Point Midpoint(void) const;
        const R2Point Centroid(void) const;
        const RNLength Length(void) const;
	const RNScalar T(const R2Point& point) const;
        const RNBoolean IsPoint(void) const;
	const RNBoolean operator==(const R2Span& span) const;
	const RNBoolean operator!=(const R2Span& span) const;
        const R2Box BBox(void) const;
        const R2Circle BCircle(void) const;

	// Manipulation functions/operators
	void Flip(void);
	void Project(const R2Line& line);
	void Mirror(const R2Line& line);
        void Translate(const R2Vector& vector);
        void Reposition(int k, const R2Point& point);
        void Align(const R2Vector& vector);
        void Transform(const R2Transformation& transformation);
        void InverseTransform(const R2Transformation& transformation);
        void Reset(const R2Point& point1, const R2Point& point2);
	RNClassID Clip(const R2Line& line);

        // Draw functions/operators
        void Draw(void) const;

	// Arithmetic functions/operators
	R2Span operator-(void) const;
	
    private:
	R2Ray ray;
	R2Point end;
	RNLength length;
};



/* Inline functions */

inline const R2Point& R2Span::
Start (void) const
{
    // Return start point of span
    return ray.Start();
}



inline const R2Point& R2Span::
End (void) const
{
    // Return end point of span
    return end;
}



inline const R2Vector& R2Span::
Vector(void) const
{
    // Return direction vector of span 
    return ray.Vector();
}



inline const R2Vector& R2Span::
Normal(void) const
{
    // Return normal vector of span 
    return ray.Normal();
}



inline const R2Point& R2Span::
operator[] (int k) const
{
    // Return kth endpoint of span
    assert((k>=0)&&(k<=1));
    return (k==0) ? Start() : End();
}



inline const R2Point R2Span::
Point (int k) const
{
    // Return kth endpoint of span
    return (*this)[k];
}



inline const R2Ray& R2Span::
Ray(void) const
{
    // Return ray along span
    return ray;
}



inline const R2Line& R2Span::
Line(void) const
{
    // Return line containing span
    return ray.Line();
}



inline const R2Point R2Span::
Centroid(void) const
{
    // Return midpoint of span
    return Midpoint();
}



inline const RNLength R2Span::
Length(void) const
{
    // Return length of span
    return length;
}



inline const RNScalar R2Span::
T(const R2Point& point) const
{
    // Return parametric value of closest point on span
    return ray.T(point);
}



inline const RNBoolean R2Span::
operator!=(const R2Span& span) const
{
    // Return whether span is not equal
    return (!(*this == span));
}



inline R2Span R2Span::
operator-(void) const
{
    // Return span with flipped orientation
    return R2Span(End(), Start());
}



