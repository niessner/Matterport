/* Include file for the R2 polyline class */



/* Initialization functions */

int R2InitPolyline();
void R2StopPolyline();



/* Class definition */

class R2Polyline : public R2Curve {
    public:
        // Constructor functions
	R2Polyline(void);
        R2Polyline(const R2Polyline& polyline);
        R2Polyline(const RNArray<R2Point *>& points);
        R2Polyline(const R2Point *points, int npoints);
        ~R2Polyline(void);

        // Polyline property functions/operators
        const int NPoints(void) const;
        const R2Point& Point(int k) const;
	const RNBoolean IsEmpty(void) const;
	const RNBoolean IsFinite(void) const;
        const R2Point& operator[](int k) const;
        const R2Point ClosestPoint(const R2Point& point) const;

        // Shape property functions/operators
	virtual const RNBoolean IsPoint(void) const;
	virtual const RNBoolean IsLinear(void) const ;
	virtual const RNLength Length(void) const;
	virtual const R2Point Centroid(void) const;
	virtual const R2Shape& BShape(void) const;
	virtual const R2Box BBox(void) const;
	virtual const R2Circle BCircle(void) const;

        // Point properties
        R2Vector Normal(int k, RNLength radius = 0) const;
        R2Line Tangent(int k, RNLength radius = 0) const;
        RNAngle InteriorAngle(int k, RNLength radius = 0) const;
        RNScalar Curvature(int k, RNLength radius = 0) const;

        // Manipulation functions/operators
        virtual void Empty(void);
	virtual void Transform(const R2Transformation& transformation);

        // Draw functions/operators
        virtual void Draw(const R2DrawFlags draw_flags = R2_DEFAULT_DRAW_FLAGS) const;
        virtual void Print(FILE *fp = stdout) const;

    private:
        R2Point *points;
        int npoints;
        R2Box bbox;
};



/* Inline functions */

inline const int R2Polyline::
NPoints(void) const
{
  // Return number of points
  return npoints;
}



inline const R2Point& R2Polyline::
Point(int k) const
{
  // Return Kth point
  return points[k];
}



inline const RNBoolean R2Polyline::
IsEmpty(void) const
{
    // Return whether the polyline is null
    return (npoints == 0);
}



inline const RNBoolean R2Polyline::
IsFinite(void) const
{
    // Return whether the polyline has finite extent
    return bbox.IsFinite();
}




inline const R2Point& R2Polyline::
operator[](int k) const
{
  // Return Kth point
  return Point(k);
}





