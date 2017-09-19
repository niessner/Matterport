/* Include file for the R2 Box class */



/* Initialization functions */

int R2InitBox();
void R2StopBox();



/* Class definition */

class R2Box : public R2Solid {
    public:
        // Constructor functions
	R2Box(void);
        R2Box(const R2Box& box);
	R2Box(const R2Point& min, const R2Point& max);
	R2Box(const R2Point& center, RNLength xradius, RNLength yradius);
	R2Box(RNCoord xmin, RNCoord ymin, RNCoord xmax, RNCoord ymax);

        // Box property functions/operators
	const RNCoord XMin(void) const;
	const RNCoord YMin(void) const;
	const RNCoord XMax(void) const;
	const RNCoord YMax(void) const;
	const RNCoord Coord(RNDirection dir, RNDimension dim) const;
	const R2Point& Min(void) const;
	const R2Point& Max(void) const;
	const R2Point Corner(RNQuadrant quadrant) const;
	const R2Point Corner(RNDirection xdir, RNDirection ydir) const;
	const R2Point Centroid(void) const;
	const R2Point ClosestPoint(const R2Point& point) const;
	const R2Box Quadrant(RNQuadrant octant) const;
	const R2Box Quadrant(RNDirection xdir, RNDirection ydir) const;
	const RNBoolean IsEmpty(void) const;
	const RNBoolean IsFinite(void) const;
        const int NDimensions(void) const;
	const RNBoolean IsAxisNull(const RNAxis axis) const;
	const RNLength XLength(void) const;
	const RNLength YLength(void) const;
	const RNLength AxisLength(const RNAxis axis) const;
	const RNLength DiagonalLength(void) const;
	const RNLength XRadius(void) const;
	const RNLength YRadius(void) const;
	const RNLength AxisRadius(const RNAxis axis) const;
	const RNLength DiagonalRadius(void) const;
	const RNCoord XCenter(void) const;
	const RNCoord YCenter(void) const;
	const RNCoord AxisCenter(const RNAxis axis) const;
	const RNAxis ShortestAxis(void) const;
	const RNAxis LongestAxis(void) const;
	const RNLength ShortestAxisLength(void) const;
	const RNLength LongestAxisLength(void) const;

	// Shape property functions
	virtual const RNBoolean IsPoint(void) const;
	virtual const RNBoolean IsLinear(void) const;
	virtual const RNBoolean IsConvex(void) const;
	virtual const RNArea Area(void) const;
	virtual const R2Shape& BShape(void) const;
	virtual const R2Box BBox(void) const;
	virtual const R2Circle BCircle(void) const;

        // Manipulation functions/operators
        virtual void Empty(void);
        virtual void Inflate(RNScalar fraction);
        virtual void Translate(const R2Vector& vector);
	virtual void Union(const R2Point& point);
	virtual void Union(const R2Box& box);
	virtual void Union(const R2Circle& circle);
	virtual void Intersect(const R2Box& box);
	virtual void Transform(const R2Transformation& transformation);
	virtual void Reset(const R2Point& min, const R2Point& max);

        // Draw functions/operators
        virtual void Draw(const R2DrawFlags draw_flags = R2_DEFAULT_DRAW_FLAGS) const;

        // Relationship functions/operators
	RNBoolean operator==(const R2Box& box) const;
	RNBoolean operator!=(const R2Box& box) const;

        // Undocumented functions/operators
	const R2Point& operator[](RNDirection dir) const;
	R2Point& operator[](RNDirection dir);

    private:
	R2Point minpt;
	R2Point maxpt;
};



/* Public variables */

extern const R2Box R2null_box;
extern const R2Box R2zero_box;
extern const R2Box R2unit_box;
extern const R2Box R2infinite_box;



/* Inline functions */

inline const RNCoord R2Box::
XMin (void) const
{
    // Return X coordinate of low corner
    return minpt.X();
}



inline const RNCoord R2Box::
YMin (void) const
{
    // Return Y coordinate of low corner
    return minpt.Y();
}



inline const RNCoord R2Box::
XMax (void) const
{
    // Return X coordinate of high corner
    return maxpt.X();
}



inline const RNCoord R2Box::
YMax (void) const
{
    // Return Y coordinate of high corner
    return maxpt.Y();
}



inline const RNCoord R2Box::
Coord (RNDirection dir, RNDimension dim) const
{
    // Return requested coordinate 
    return (dir == RN_LO) ? minpt[dim] : maxpt[dim];
}



inline const R2Point& R2Box::
Min (void) const
{
    // Return point (min, min, min)
    return minpt;
}



inline const R2Point& R2Box::
Max (void) const
{
    // Return point (max, max, max)
    return maxpt;
}



inline const R2Box R2Box::
Quadrant (RNQuadrant quadrant) const
{
    // Return box in quadrant (xdir, ydir)
    return Quadrant((quadrant & RN_PN_QUADRANT) ? RN_HI : RN_LO, 
		    (quadrant & RN_NP_QUADRANT) ? RN_HI : RN_LO);
}



inline const RNBoolean R2Box::
IsEmpty (void) const
{
    // Return whether bounding box contains no space
    return ((minpt.X() > maxpt.X()) || (minpt.Y() > maxpt.Y()));
}



inline const RNBoolean R2Box::
IsAxisNull (const RNAxis axis) const
{
    // Return whether bounding box contains just one point along axis
    return (minpt[axis] == maxpt[axis]);
}




inline const RNLength R2Box::
AxisLength (const RNAxis axis) const
{
    // Return length of Box along axis
    return (maxpt[axis] - minpt[axis]);
}



inline const RNLength R2Box::
XLength (void) const
{
    // Return length in X dimension
    return this->AxisLength(RN_XAXIS);
}



inline const RNLength R2Box::
YLength (void) const
{
    // Return length in Y dimension
    return this->AxisLength(RN_YAXIS);
}



inline const RNLength R2Box::
AxisRadius (const RNAxis axis) const
{
    // Return radius of Box along axis
    return 0.5 * (maxpt[axis] - minpt[axis]);
}



inline const RNLength R2Box::
XRadius (void) const
{
    // Return radius in X dimension
    return this->AxisRadius(RN_XAXIS);
}



inline const RNLength R2Box::
YRadius (void) const
{
    // Return radius in Y dimension
    return this->AxisRadius(RN_YAXIS);
}



inline const RNLength R2Box::
DiagonalRadius (void) const
{
    // Return radius of Box along diagonal
    return (0.5 * DiagonalLength());
}



inline const RNCoord R2Box::
AxisCenter (const RNAxis axis) const
{
    // Return radius of Box along axis
    return 0.5 * (maxpt[axis] + minpt[axis]);
}



inline const RNCoord R2Box::
XCenter (void) const
{
    // Return center in X dimension
    return this->AxisCenter(RN_XAXIS);
}



inline const RNCoord R2Box::
YCenter (void) const
{
    // Return center in Y dimension
    return this->AxisCenter(RN_YAXIS);
}



inline RNBoolean R2Box::
operator==(const R2Box& box) const
{
    // Return whether box is equal
    return ((minpt == box.minpt) && (maxpt == box.maxpt));
}



inline RNBoolean R2Box::
operator!=(const R2Box& box) const
{
    // Return whether box is not equal
    return (!(*this == box));
}



inline const R2Point& R2Box::
operator[] (RNDirection dir) const
{
    // Return min or max point 
    return (dir == RN_LO) ? minpt : maxpt;
}



inline R2Point& R2Box::
operator[] (RNDirection dir) 
{
    // Return min or max point 
    return (dir == RN_LO) ? minpt : maxpt;
}



