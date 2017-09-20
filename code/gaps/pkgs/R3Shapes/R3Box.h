/* Include file for the R3 Box class */



/* Initialization functions */

int R3InitBox();
void R3StopBox();



/* Class definition */

class R3Box : public R3Solid {
    public:
        // Constructor functions
	R3Box(void);
        R3Box(const R3Box& box);
	R3Box(const R3Point& min, const R3Point& max);
	R3Box(RNCoord xmin, RNCoord ymin, RNCoord zmin, 
	      RNCoord xmax, RNCoord ymax, RNCoord zmax);

        // Box property functions/operators
	const RNCoord XMin(void) const;
	const RNCoord YMin(void) const;
	const RNCoord ZMin(void) const;
	const RNCoord XMax(void) const;
	const RNCoord YMax(void) const;
	const RNCoord ZMax(void) const;
	const RNCoord Coord(RNDirection dir, RNDimension dim) const;
	const R3Point& Min(void) const;
	const R3Point& Max(void) const;
	const R3Point Corner(RNOctant octant) const;
	const R3Point Corner(RNDirection xdir, RNDirection ydir, RNDirection zdir) const;
	const R3Plane Plane(RNSide side) const;
	const R3Plane Plane(RNDirection dir, RNDimension dim) const;
	const R3Box Side(RNSide side) const;
	const R3Box Side(RNDirection dir, RNDimension dim) const;
	const R3Box Octant(RNOctant octant) const;
	const R3Box Octant(RNDirection xdir, RNDirection ydir, RNDirection zdir) const;
	const RNBoolean IsEmpty(void) const;
	const RNBoolean IsFinite(void) const;
        const int NDimensions(void) const;
        const int NDimensionsAlongSide(RNDirection dir, RNDimension dim) const;
	const RNBoolean IsSideNull(RNDirection dir, RNDimension dim) const;
	const RNBoolean IsSideLinear(RNDirection dir, RNDimension dim) const;	
        const int NDimensionsAlongAxis(const RNAxis axis) const;
	const RNBoolean IsAxisNull(const RNAxis axis) const;
	const RNLength XLength(void) const;
	const RNLength YLength(void) const;
	const RNLength ZLength(void) const;
	const RNLength AxisLength(const RNAxis axis) const;
	const RNLength DiagonalLength(void) const;
	const RNLength XRadius(void) const;
	const RNLength YRadius(void) const;
	const RNLength ZRadius(void) const;
	const RNLength AxisRadius(const RNAxis axis) const;
	const RNLength DiagonalRadius(void) const;
	const RNCoord XCenter(void) const;
	const RNCoord YCenter(void) const;
	const RNCoord ZCenter(void) const;
	const RNCoord AxisCenter(const RNAxis axis) const;
	const RNAxis ShortestAxis(void) const;
	const RNAxis LongestAxis(void) const;
	const RNLength ShortestAxisLength(void) const;
	const RNLength LongestAxisLength(void) const;

	// Shape property functions
	virtual const RNBoolean IsPoint(void) const;
	virtual const RNBoolean IsLinear(void) const;
	virtual const RNBoolean IsPlanar(void) const;
	virtual const RNBoolean IsConvex(void) const;
        virtual const RNInterval NFacets(void) const;
	virtual const RNArea Area(void) const;
	virtual const RNVolume Volume(void) const;
	virtual const R3Point Centroid(void) const;
	virtual const R3Point ClosestPoint(const R3Point& point) const;
	virtual const R3Point FurthestPoint(const R3Point& point) const;
	virtual const R3Shape& BShape(void) const;
	virtual const R3Box BBox(void) const;
	virtual const R3Sphere BSphere(void) const;

        // Manipulation functions/operators
        virtual void Empty(void);
        virtual void Inflate(RNScalar fraction);
        virtual void Translate(const R3Vector& vector);
	virtual void Union(const R3Point& point);
	virtual void Union(const R3Box& box);
	virtual void Union(const R3Sphere& sphere);
	virtual void Intersect(const R3Box& box);
	virtual void Transform(const R3Transformation& transformation);
	virtual void Reset(const R3Point& min, const R3Point& max);
	
        // Draw functions/operators
        virtual void Draw(const R3DrawFlags draw_flags = R3_DEFAULT_DRAW_FLAGS) const;

        // Relationship functions/operators
	RNBoolean operator==(const R3Box& box) const;
	RNBoolean operator!=(const R3Box& box) const;

	// Standard shape definitions
	RN_CLASS_TYPE_DECLARATIONS(R3Box);
        R3_SHAPE_RELATIONSHIP_DECLARATIONS(R3Box);

        // Undocumented functions/operators
	const R3Point& operator[](RNDirection dir) const;
	R3Point& operator[](RNDirection dir);

    private:
	R3Point minpt;
	R3Point maxpt;
};



/* Public variables */

extern const R3Box R3null_box;
extern const R3Box R3zero_box;
extern const R3Box R3unit_box;
extern const R3Box R3infinite_box;



/* Inline functions */

inline const RNCoord R3Box::
XMin (void) const
{
    // Return X coordinate of low corner
    return minpt.X();
}



inline const RNCoord R3Box::
YMin (void) const
{
    // Return Y coordinate of low corner
    return minpt.Y();
}



inline const RNCoord R3Box::
ZMin (void) const
{
    // Return Z coordinate of low corner
    return minpt.Z();
}



inline const RNCoord R3Box::
XMax (void) const
{
    // Return X coordinate of high corner
    return maxpt.X();
}



inline const RNCoord R3Box::
YMax (void) const
{
    // Return Y coordinate of high corner
    return maxpt.Y();
}



inline const RNCoord R3Box::
ZMax (void) const
{
    // Return Z coordinate of high corner
    return maxpt.Z();
}



inline const RNCoord R3Box::
Coord (RNDirection dir, RNDimension dim) const
{
    // Return requested coordinate 
    return (dir == RN_LO) ? minpt[dim] : maxpt[dim];
}



inline const R3Point& R3Box::
Min (void) const
{
    // Return point (min, min, min)
    return minpt;
}



inline const R3Point& R3Box::
Max (void) const
{
    // Return point (max, max, max)
    return maxpt;
}



inline const R3Box R3Box::
Side (RNSide side) const
{
    // Return box along side
    return Side(side & 0x1, side >> 1);
}



inline const R3Plane R3Box::
Plane (RNSide side) const
{
    // Return box along side
    return Plane(side & 0x1, side >> 1);
}



inline const R3Box R3Box::
Octant (RNOctant octant) const
{
    // Return box in octant (xdir, ydir, zdir)
    return Octant((octant & RN_PNN_OCTANT) ? RN_HI : RN_LO, 
		  (octant & RN_NPN_OCTANT) ? RN_HI : RN_LO, 
		  (octant & RN_NNP_OCTANT) ? RN_HI : RN_LO);
}



inline const RNBoolean R3Box::
IsEmpty (void) const
{
    // Return whether bounding box contains no space
    return ((minpt.X() > maxpt.X()) || (minpt.Y() > maxpt.Y()) || (minpt.Z() > maxpt.Z()));
}



inline const RNBoolean R3Box::
IsSideNull (RNDirection dir, RNDimension dim) const
{
    // Return whether bounding box contains just one point along side
    return (this->NDimensionsAlongSide(dir, dim) == 0);
}



inline const RNBoolean R3Box::
IsSideLinear (RNDirection dir, RNDimension dim) const
{
    // Return whether bounding box contains just one line segment along side
    return (this->NDimensionsAlongSide(dir, dim) <= 1);
}



inline const RNBoolean R3Box::
IsAxisNull (const RNAxis axis) const
{
    // Return whether bounding box contains just one point along axis
    return (minpt[axis] == maxpt[axis]);
}




inline const RNLength R3Box::
AxisLength (const RNAxis axis) const
{
    // Return length of Box along axis
    return (maxpt[axis] - minpt[axis]);
}



inline const RNLength R3Box::
XLength (void) const
{
    // Return length in X dimension
    return this->AxisLength(RN_XAXIS);
}



inline const RNLength R3Box::
YLength (void) const
{
    // Return length in Y dimension
    return this->AxisLength(RN_YAXIS);
}



inline const RNLength R3Box::
ZLength (void) const
{
    // Return length in Z dimension
    return this->AxisLength(RN_ZAXIS);
}



inline const RNLength R3Box::
AxisRadius (const RNAxis axis) const
{
    // Return radius of Box along axis
    return 0.5 * (maxpt[axis] - minpt[axis]);
}



inline const RNLength R3Box::
XRadius (void) const
{
    // Return radius in X dimension
    return this->AxisRadius(RN_XAXIS);
}



inline const RNLength R3Box::
YRadius (void) const
{
    // Return radius in Y dimension
    return this->AxisRadius(RN_YAXIS);
}



inline const RNLength R3Box::
ZRadius (void) const
{
    // Return radius in Z dimension
    return this->AxisRadius(RN_ZAXIS);
}



inline const RNLength R3Box::
DiagonalRadius (void) const
{
    // Return radius of Box along diagonal
    return (0.5 * DiagonalLength());
}



inline const RNCoord R3Box::
AxisCenter (const RNAxis axis) const
{
    // Return radius of Box along axis
    return 0.5 * (maxpt[axis] + minpt[axis]);
}



inline const RNCoord R3Box::
XCenter (void) const
{
    // Return center in X dimension
    return this->AxisCenter(RN_XAXIS);
}



inline const RNCoord R3Box::
YCenter (void) const
{
    // Return center in Y dimension
    return this->AxisCenter(RN_YAXIS);
}



inline const RNCoord R3Box::
ZCenter (void) const
{
    // Return center in Z dimension
    return this->AxisCenter(RN_ZAXIS);
}



inline RNBoolean R3Box::
operator==(const R3Box& box) const
{
    // Return whether box is equal
    return ((minpt == box.minpt) && (maxpt == box.maxpt));
}



inline RNBoolean R3Box::
operator!=(const R3Box& box) const
{
    // Return whether box is not equal
    return (!(*this == box));
}



inline const R3Point& R3Box::
operator[] (RNDirection dir) const
{
    // Return min or max point 
    return (dir == RN_LO) ? minpt : maxpt;
}



inline R3Point& R3Box::
operator[] (RNDirection dir) 
{
    // Return min or max point 
    return (dir == RN_LO) ? minpt : maxpt;
}



