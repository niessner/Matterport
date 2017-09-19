/* Include file for the GAPS span class */



/* Initialization functions */

int R3InitSpan();
void R3StopSpan();



/* Class definition */

class R3Span {
    public:
        // Constructor functions
	R3Span(void);
	R3Span(const R3Span& span);
	R3Span(const R3Point& point, const R3Vector& vector);
	R3Span(const R3Point& point1, const R3Point& point2);
	R3Span(RNCoord x1, RNCoord y1, RNCoord z1, RNCoord x2, RNCoord y2, RNCoord z2);

        // Property functions/operators
        const R3Point& Start(void) const;
        const R3Point& End(void) const;
        const R3Vector& Vector(void) const;
        const R3Point Point(int k) const;
	const R3Point Point(RNScalar t) const;
        const R3Point& operator[](int k) const;
        const R3Ray& Ray(void) const;
        const R3Line& Line(void) const;
        const R3Point Midpoint(void) const;
        const R3Point Centroid(void) const;
        const R3Box BBox(void) const;
        const R3Sphere BSphere(void) const;
        const RNLength Length(void) const;
	const RNScalar T(const R3Point& point) const;
        const RNBoolean IsPoint(void) const;
	const RNBoolean operator==(const R3Span& span) const;
	const RNBoolean operator!=(const R3Span& span) const;

	// Manipulation functions/operators
	void Flip(void);
	void Mirror(const R3Plane& plane);
        void Translate(const R3Vector& vector);
        void Reposition(int k, const R3Point& point);
        void Align(const R3Vector& vector);
        void Transform(const R3Transformation& transformation);
        void InverseTransform(const R3Transformation& transformation);
        void Reset(const R3Point& point1, const R3Point& point2);
	RNClassID Clip(const R3Plane& plane);

        // Draw functions/operators
        void Draw(void) const;

	// Arithmetic functions/operators
	R3Span operator-(void) const;

    private:
	R3Ray ray;
	R3Point end;
	RNLength length;
};



/* Public variables */

extern R3Span R3null_span;




/* Inline functions */

inline const R3Point& R3Span::
Start (void) const
{
    // Return start point of span
    return ray.Start();
}



inline const R3Point& R3Span::
End (void) const
{
    // Return end point of span
    return end;
}



inline const R3Vector& R3Span::
Vector(void) const
{
    // Return direction vector of span 
    return ray.Vector();
}



inline const R3Point& R3Span::
operator[] (int k) const
{
    // Return kth endpoint of span
    assert((k>=0)&&(k<=1));
    return (k==0) ? Start() : End();
}



inline const R3Point R3Span::
Point (int k) const
{
    // Return kth endpoint of span
    return (*this)[k];
}



inline const R3Ray& R3Span::
Ray(void) const
{
    // Return ray along span
    return ray;
}



inline const R3Line& R3Span::
Line(void) const
{
    // Return line containing span
    return ray.Line();
}



inline const R3Point R3Span::
Centroid(void) const
{
    // Return midpoint of span
    return Midpoint();
}



inline const RNLength R3Span::
Length(void) const
{
    // Return length of span
    return length;
}



inline const RNBoolean R3Span::
operator!=(const R3Span& span) const
{
    // Return whether span is not equal
    return (!(*this == span));
}



inline R3Span R3Span::
operator-(void) const
{
    // Return span with flipped orientation
    return R3Span(End(), Start());
}



inline void R3Span::
Mirror(const R3Plane& plane)
{
    // Mirror span over plane
    ray.Mirror(plane);
    end.Mirror(plane);
}





