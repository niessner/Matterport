/* Include file for the R2 circle class */



/* Initialization functions */

int R2InitCircle();
void R2StopCircle();



/* Class definition */

class R2Circle : public R2Solid {
    public:
        // Constructor functions
	R2Circle(void);
        R2Circle(const R2Circle& circle);
        R2Circle(const R2Point& center, RNLength radius);

        // Circle property functions/operators
        const R2Point& Center(void) const;
        const RNLength Radius(void) const;
	const RNBoolean IsEmpty(void) const;
	const RNBoolean IsFinite(void) const;

        // Shape property functions/operators
	virtual const RNBoolean IsPoint(void) const;
	virtual const RNBoolean IsLinear(void) const ;
	virtual const RNBoolean IsConvex(void) const ;
	virtual const RNArea Area(void) const;
	virtual const R2Point Centroid(void) const;
	virtual const R2Shape& BShape(void) const;
	virtual const R2Box BBox(void) const;
	virtual const R2Circle BCircle(void) const;

        // Manipulation functions/operators
        virtual void Empty(void);
        virtual void Translate(const R2Vector& vector);
        virtual void Reposition(const R2Point& center);
        virtual void Resize(RNLength radius);
	virtual void Transform(const R2Transformation& transformation);

        // Draw functions/operators
        virtual void Draw(const R2DrawFlags draw_flags = R2_DEFAULT_DRAW_FLAGS) const;

    private:
        R2Point center;
        RNLength radius;
};



/* Public variables */

extern const R2Circle R2null_circle;
extern const R2Circle R2zero_circle;
extern const R2Circle R2unit_circle;
extern const R2Circle R2infinite_circle;

extern const int R2circle_npoints;
extern RNAngle R2circle_angles[];
extern R2Point R2circle_points[];



/* Inline functions */

inline const R2Point& R2Circle::
Center(void) const
{
    // Return circle center
    return center;
}



inline const RNLength R2Circle::
Radius(void) const
{
    // Return circle radius
    return radius;
}



inline const RNBoolean R2Circle::
IsEmpty(void) const
{
    // Return whether the circle is null
    return (radius < 0.0);
}



inline const RNBoolean R2Circle::
IsFinite(void) const
{
    // Return whether the circle has finite radius
    return RNIsFinite(radius);
}






