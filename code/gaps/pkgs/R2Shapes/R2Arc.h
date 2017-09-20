/* Include file for the R2 arc class */



/* Initialization functions */

int R2InitArc();
void R2StopArc();



/* Class definition */

class R2Arc : public R2Curve {
    public:
        // Constructor functions
	R2Arc(void);
        R2Arc(const R2Arc& arc);
        R2Arc(const R2Point& center, RNLength radius, RNAngle start, RNAngle stop);

        // Arc property functions/operators
        const R2Circle& Circle(void) const;
        const R2Point& Center(void) const;
        const RNLength Radius(void) const;
        const RNAngle StartAngle(void) const;
        const RNAngle StopAngle(void) const;
        const RNAngle SweepAngle(void) const;
        const R2Point Point(RNAngle angle) const;
        const R2Point MidPoint(void) const;
        const R2Point StartPoint(void) const;
        const R2Point StopPoint(void) const;
	const RNBoolean IsEmpty(void) const;
	const RNBoolean IsFinite(void) const;

        // Shape property functions/operators
	virtual const RNBoolean IsPoint(void) const;
	virtual const RNBoolean IsLinear(void) const ;
	virtual const RNBoolean IsConvex(void) const ;
	virtual const R2Point Centroid(void) const;
	virtual const R2Shape& BShape(void) const;
	virtual const R2Box BBox(void) const;
	virtual const R2Circle BCircle(void) const;

        // Manipulation functions/operators
        virtual void Empty(void);
	virtual void SetStartAngle(RNAngle theta);
	virtual void SetStopAngle(RNAngle theta);
        virtual void Translate(const R2Vector& vector);
        virtual void Reposition(const R2Point& center);
        virtual void Resize(RNLength radius);
	virtual void Transform(const R2Transformation& transformation);

        // Draw functions/operators
        virtual void Draw(const R2DrawFlags draw_flags = R2_DEFAULT_DRAW_FLAGS) const;

    private:
        R2Circle circle;
	RNAngle start;
	RNAngle stop;
};



/* Public variables */

extern const R2Arc R2null_arc;
extern const R2Arc R2zero_arc;



/* Inline functions */

inline const R2Circle& R2Arc::
Circle(void) const
{
    // Return arc circle
    return circle;
}



inline const R2Point& R2Arc::
Center(void) const
{
    // Return arc center
    return circle.Center();
}



inline const RNLength R2Arc::
Radius(void) const
{
    // Return arc radius
    return circle.Radius();
}



inline const RNAngle R2Arc::
StartAngle(void) const
{
    // Return arc start angle
    return start;
}



inline const RNAngle R2Arc::
StopAngle(void) const
{
    // Return arc start angle
    return stop;
}



inline const RNAngle R2Arc::
SweepAngle(void) const
{
    // Return angle within arc sector
    return stop - start;
}



inline const R2Point R2Arc::
MidPoint(void) const
{
    // Return midpoint on arc
    return Point(0.5 * (start + stop));
}



inline const R2Point R2Arc::
StartPoint(void) const
{
    // Return start point on arc
    return Point(start);
}



inline const R2Point R2Arc::
StopPoint(void) const
{
    // Return start point on arc
    return Point(stop);
}



inline const RNBoolean R2Arc::
IsEmpty(void) const
{
    // Return whether the arc is null
    return circle.IsEmpty();
}



inline const RNBoolean R2Arc::
IsFinite(void) const
{
    // Return whether the arc has finite radius
    return circle.IsFinite();
}






