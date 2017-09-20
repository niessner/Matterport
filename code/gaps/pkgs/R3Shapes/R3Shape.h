/* Include file for the R3 shape class */



/* Initialization functions */

int R3InitShape();
void R3StopShape();



/* Class definition */

class R3Shape /* : public R3Base */ {
    public:
        // Constructors/destructors
        virtual ~R3Shape(void) {};

        // Property functions/operators
	virtual const RNBoolean IsPoint(void) const;
	virtual const RNBoolean IsCurve(void) const;
	virtual const RNBoolean IsLinear(void) const;
	virtual const RNBoolean IsSurface(void) const;
	virtual const RNBoolean IsPlanar(void) const;
	virtual const RNBoolean IsSolid(void) const;
	virtual const RNBoolean IsConvex(void) const;
	virtual const RNInterval NFacets(void) const;
	virtual const RNLength Length(void) const;
	virtual const RNArea Area(void) const;
	virtual const RNVolume Volume(void) const;
	// virtual const R3Line *Line(void) const;
	// virtual const R3Plane *Plane(void) const;
	virtual const R3Point Centroid(void) const;
	virtual const R3Point ClosestPoint(const R3Point& point) const;
	virtual const R3Point FurthestPoint(const R3Point& point) const;
	virtual const R3Shape& BShape(void) const;
	virtual const R3Box BBox(void) const;
	virtual const R3Sphere BSphere(void) const;
        virtual R3Shape *Copy(void) const;

	// Manipulation functions/operators
	virtual void Transform(const R3Transformation& transformation);

	// Draw functions/operations
	virtual void Draw(const R3DrawFlags draw_flags = R3_DEFAULT_DRAW_FLAGS) const;
	virtual void Outline(const R3DrawFlags draw_flags = R3_EDGES_DRAW_FLAG) const;

	// Standard shape type functions/operators
	RN_CLASS_TYPE_DECLARATIONS(R3Shape);
#define R3_SHAPE_RELATIONSHIP_DECLARATIONS(__type) \
	virtual RNLength Distance(const R3Point& point) const; \
	virtual RNLength Distance(const R3Line& line) const; \
    	virtual RNLength Distance(const R3Ray& ray) const; \
    	virtual RNLength Distance(const R3Span& span) const; \
    	virtual RNLength Distance(const R3Plane& plane) const; \
    	virtual RNLength Distance(const R3Halfspace& halfspace) const; \
    	virtual RNLength Distance(const R3Box& box) const; \
    	virtual RNLength Distance(const R3OrientedBox& box) const; \
    	virtual RNLength Distance(const R3Sphere& sphere) const; \
    	virtual RNLength Distance(const R3Shape& shape) const; \
    	virtual RNLength SignedDistance(const R3Plane& plane) const; \
    	virtual RNBoolean Contains(const R3Point& point) const; \
    	virtual RNBoolean Contains(const R3Line& line) const; \
    	virtual RNBoolean Contains(const R3Ray& ray) const; \
    	virtual RNBoolean Contains(const R3Span& span) const; \
    	virtual RNBoolean Contains(const R3Plane& plane) const; \
    	virtual RNBoolean Contains(const R3Halfspace& halfspace) const; \
    	virtual RNBoolean Contains(const R3Box& box) const; \
    	virtual RNBoolean Contains(const R3OrientedBox& box) const; \
    	virtual RNBoolean Contains(const R3Sphere& sphere) const; \
    	virtual RNBoolean Contains(const R3Shape& shape) const; \
    	virtual RNBoolean Inside(const R3Point& point) const; \
    	virtual RNBoolean Inside(const R3Line& line) const; \
    	virtual RNBoolean Inside(const R3Ray& ray) const; \
    	virtual RNBoolean Inside(const R3Span& span) const; \
    	virtual RNBoolean Inside(const R3Plane& plane) const; \
    	virtual RNBoolean Inside(const R3Halfspace& halfspace) const; \
    	virtual RNBoolean Inside(const R3Box& box) const; \
    	virtual RNBoolean Inside(const R3OrientedBox& box) const; \
    	virtual RNBoolean Inside(const R3Sphere& sphere) const; \
    	virtual RNBoolean Inside(const R3Shape& shape) const; \
    	virtual RNClassID Intersects(const R3Point& point) const; \
    	virtual RNClassID Intersects(const R3Line& line) const; \
    	virtual RNClassID Intersects(const R3Ray& ray, \
	    R3Point *hit_point = NULL, R3Vector *hit_normal = NULL, RNScalar *hit_t = NULL) const; \
    	virtual RNClassID Intersects(const R3Span& span) const; \
    	virtual RNClassID Intersects(const R3Plane& plane) const; \
    	virtual RNClassID Intersects(const R3Halfspace& halfspace) const; \
    	virtual RNClassID Intersects(const R3Box& box) const; \
    	virtual RNClassID Intersects(const R3OrientedBox& box) const; \
    	virtual RNClassID Intersects(const R3Sphere& sphere) const; \
    	virtual RNClassID Intersects(const R3Shape& shape) const
	R3_SHAPE_RELATIONSHIP_DECLARATIONS(R3Shape);
};



/* Shape type enumeration ??? */

enum {
    R3_POINT_CLASS_ID = 1,
    R3_LINE_CLASS_ID = 2,
    R3_RAY_CLASS_ID = 3,
    R3_SPAN_CLASS_ID = 4,
    R3_PLANE_CLASS_ID = 5,
    R3_POLYGON_CLASS_ID = 6,
    R3_CIRCLE_CLASS_ID = 7,
    R3_ELLIPSE_CLASS_ID = 8
};



/* Other function declarations */

const R3Box R3ShapeBBox(const void *data, void *appl);




