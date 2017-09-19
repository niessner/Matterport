/* Include file for the R2 polygon class */



/* Initialization functions */

int R2InitPolygon();
void R2StopPolygon();



/* Class definition */

class R2Polygon : public R2Solid {
    public:
        // Constructor functions
	R2Polygon(void);
        R2Polygon(const R2Polygon& polygon);
        R2Polygon(const RNArray<R2Point *>& points, RNBoolean clockwise = FALSE);
        R2Polygon(const R2Point *points, int npoints, RNBoolean clockwise = FALSE);
        ~R2Polygon(void);

        // Polygon property functions/operators
        const int NPoints(void) const;
        const R2Point& Point(int k) const;
	const RNBoolean IsEmpty(void) const;
	const RNBoolean IsFinite(void) const;
	const RNBoolean IsClockwise(void) const;
        const R2Point& operator[](int k) const;
        const R2Point ClosestPoint(const R2Point& point) const;

        // Shape property functions/operators
	virtual const RNBoolean IsPoint(void) const;
	virtual const RNBoolean IsLinear(void) const ;
	virtual const RNBoolean IsConvex(void) const ;
	virtual const RNArea Area(void) const;
	virtual const RNScalar Convexity(void) const;
	virtual const RNLength Perimeter(void) const;
	virtual const R2Point Centroid(void) const;
	virtual const R2Shape& BShape(void) const;
	virtual const R2Box BBox(void) const;
	virtual const R2Circle BCircle(void) const;
	virtual const R2Polygon ConvexHull(void) const;

        // Point properties
        R2Vector Normal(int k, RNLength radius = 0) const;
        R2Line Tangent(int k, RNLength radius = 0) const;
        RNAngle InteriorAngle(int k, RNLength radius = 0) const;
        RNScalar Curvature(int k, RNLength radius = 0) const;

        // Manipulation functions/operators
        virtual void Empty(void);
        virtual void Clip(const R2Line& line);
        virtual void Clip(const R2Box& box);
	virtual void Transform(const R2Transformation& transformation);

        // Draw functions/operators
        virtual void Draw(const R2DrawFlags draw_flags = R2_DEFAULT_DRAW_FLAGS) const;
        virtual void Print(FILE *fp = stdout) const;

        // Input/output functions
        virtual int ReadTheraFile(const char *filename);

    private:
        R2Point *points;
        int npoints;
        R2Box bbox;
        RNBoolean clockwise;
};



/* Inline functions */

inline const int R2Polygon::
NPoints(void) const
{
  // Return number of points
  return npoints;
}



inline const R2Point& R2Polygon::
Point(int k) const
{
  // Return Kth point
  return points[k];
}



inline const RNBoolean R2Polygon::
IsEmpty(void) const
{
    // Return whether the polygon is null
    return (npoints == 0);
}



inline const RNBoolean R2Polygon::
IsFinite(void) const
{
    // Return whether the polygon has finite extent
    return bbox.IsFinite();
}




inline const RNBoolean R2Polygon::
IsClockwise(void) const
{
    // Return whether the polygon has points in clockwise order
    return clockwise;
}




inline const R2Point& R2Polygon::
operator[](int k) const
{
  // Return Kth point
  return Point(k);
}





