/* Include file for the R3 sphere class */



/* Initialization functions */

int R3InitSphere();
void R3StopSphere();



/* Class definition */

class R3Sphere : public R3Solid {
    public:
        // Constructor functions
	R3Sphere(void);
        R3Sphere(const R3Sphere& sphere);
        R3Sphere(const R3Point& center, RNLength radius);

        // Sphere property functions/operators
        const R3Point& Center(void) const;
        const RNLength Radius(void) const;
	const RNBoolean IsEmpty(void) const;
	const RNBoolean IsFinite(void) const;

        // Shape property functions/operators
	virtual const RNBoolean IsPoint(void) const;
	virtual const RNBoolean IsLinear(void) const ;
	virtual const RNBoolean IsPlanar(void) const ;
	virtual const RNBoolean IsConvex(void) const ;
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
        virtual void Translate(const R3Vector& vector);
        virtual void Reposition(const R3Point& center);
        virtual void Resize(RNLength radius);
	virtual void Transform(const R3Transformation& transformation);

        // Draw functions/operators
        virtual void Draw(const R3DrawFlags draw_flags = R3_DEFAULT_DRAW_FLAGS) const;

	// Standard shape definitions
	RN_CLASS_TYPE_DECLARATIONS(R3Sphere);
        R3_SHAPE_RELATIONSHIP_DECLARATIONS(R3Sphere);


    private:
        R3Point center;
        RNLength radius;
};



/* Public variables */

extern const R3Sphere R3null_sphere;
extern const R3Sphere R3zero_sphere;
extern const R3Sphere R3unit_sphere;
extern const R3Sphere R3infinite_sphere;



/* Inline functions */

inline const R3Point& R3Sphere::
Center(void) const
{
    // Return sphere center
    return center;
}



inline const RNLength R3Sphere::
Radius(void) const
{
    // Return sphere radius
    return radius;
}



inline const RNBoolean R3Sphere::
IsEmpty(void) const
{
    // Return whether the sphere is null
    return (radius < 0.0);
}



inline const RNBoolean R3Sphere::
IsFinite(void) const
{
    // Return whether the sphere has finite radius
    return RNIsFinite(radius);
}






