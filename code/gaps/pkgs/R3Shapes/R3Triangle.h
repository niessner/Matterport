/* Include file for the R3 triangle class */



/* Initialization functions */

int R3InitTriangle();
void R3StopTriangle();



/* Triangle vertex class definition */

class R3TriangleVertex {
    public:
        // Constructors
	R3TriangleVertex(void);
        R3TriangleVertex(const R3TriangleVertex& vertex);
        R3TriangleVertex(const R3Point& position);
        R3TriangleVertex(const R3Point& position, const R3Vector& normal);
        R3TriangleVertex(const R3Point& position, const RNRgb& color);
        R3TriangleVertex(const R3Point& position, const R2Point& texcoords);
        R3TriangleVertex(const R3Point& position, const R3Vector& normal, const RNRgb& color);
        R3TriangleVertex(const R3Point& position, const R3Vector& normal, const RNRgb& color, const R2Point& texcoords);
        R3TriangleVertex(const R3Point& position, const R3Vector& normal, const R2Point& texcoords);

	// Property functions/operators
	const R3Point& Position(void) const;
	const R3Vector& Normal(void) const;
        const RNRgb& Color(void) const;
	const R2Point& TextureCoords(void) const;
        const RNFlags Flags(void) const;
        const RNMark Mark(void) const;

        // Query functions/operations
        RNBoolean HasNormal(void) const;
        RNBoolean HasColor(void) const;
        RNBoolean HasTextureCoords(void) const;

	// Manipulation functions/operators
        void Mirror(const R3Plane& plane);
        void Transform(const R3Transformation& transformation);
	void SetPosition(const R3Point& position);
	void SetNormal(const R3Vector& normal);
        void SetColor(const RNRgb& color);
	void SetTextureCoords(const R2Point& texcoords);
        void SetSharedFlag(void);
        void SetMark(RNMark mark);

    public:
        R3Point position;
	R3Vector normal;
        RNRgb color;
	R2Point texcoords;
        RNFlags flags;
        RNMark mark;
};



/* Triangle class definition */

class R3Triangle : public R3Surface {
    public:
        // Constructor functions
        R3Triangle(void);
        R3Triangle(const R3Triangle& triangle);
        R3Triangle(R3TriangleVertex *v0, R3TriangleVertex *v1, R3TriangleVertex *v2);
        R3Triangle(R3TriangleVertex *vertices[3]);
        virtual ~R3Triangle(void);
  
	// Vertex access functions/operators
	R3TriangleVertex *Vertex(int index) const;
	R3TriangleVertex *V0(void) const;
	R3TriangleVertex *V1(void) const;
	R3TriangleVertex *V2(void) const;

	// Triangle property functions/operators
	const R3Plane& Plane(void) const;
	const R3Vector& Normal(void) const;
        const R3Box& Box(void) const;
        const RNFlags Flags(void) const;
        const RNMark Mark(void) const;
	const RNBoolean IsFinite(void) const;
	const RNBoolean IsDegenerate(void) const;

        // Shape property functions/operators
	virtual const RNBoolean IsPoint(void) const;
	virtual const RNBoolean IsLinear(void) const;
	virtual const RNBoolean IsPlanar(void) const;
	virtual const RNBoolean IsConvex(void) const;
        virtual const RNInterval NFacets(void) const;
        virtual const RNLength Length(void) const;
        virtual const RNArea Area(void) const;
        virtual const R3Point Centroid(void) const;
        virtual const R3Point RandomPoint(void) const;
        virtual const R3Point ClosestPoint(const R3Point& point) const;
        virtual const R3Point FurthestPoint(const R3Point& point) const;
        virtual const R3Shape& BShape(void) const;
        virtual const R3Box BBox(void) const;
        virtual const R3Sphere BSphere(void) const;

        // Manipulation functions/operators
	virtual void Flip(void);
	virtual void Transform(const R3Transformation& transformation);
	virtual void Mirror(const R3Plane& plane);
	virtual void Reset(R3TriangleVertex *v1, R3TriangleVertex *v2, R3TriangleVertex *v3);
	virtual void Reset(R3TriangleVertex *vertices[3]);
        virtual void SetMark(RNMark mark);
	virtual void Update(void);  

        // Draw functions/operators
        virtual void Draw(const R3DrawFlags draw_flags = R3_DEFAULT_DRAW_FLAGS) const;

	// Arithmetic functions/operators
	R3Triangle operator-(void) const;
	
	// Standard shape definitions
	RN_CLASS_TYPE_DECLARATIONS(R3Triangle);
        R3_SHAPE_RELATIONSHIP_DECLARATIONS(R3Triangle);

    private:
	R3TriangleVertex *v[3];
        R3Plane plane;
        R3Box bbox;
        RNFlags flags;
        RNMark mark;
};



/* Public functions -- for debugging */

void R3PrintTriangle(const R3Triangle *triangle);



/* Triangle vertex inline functions **************/

inline const R3Point& R3TriangleVertex::
Position(void) const
{
    // Return vertex position
    return position;
}



inline const R3Vector& R3TriangleVertex::
Normal(void) const
{
    // Return vertex normal
    return normal;
}



inline const RNRgb& R3TriangleVertex::
Color(void) const
{
    // Return color
    return color;
}



inline const R2Point& R3TriangleVertex::
TextureCoords(void) const
{
    // Return vertex texture coordinates
    return texcoords;
}



inline const RNFlags R3TriangleVertex::
Flags(void) const
{
    // Return vertex flags
    return flags;
}



inline const RNMark R3TriangleVertex::
Mark(void) const
{
    // Return vertex mark
    return mark;
}



inline RNBoolean R3TriangleVertex::
HasNormal(void) const
{
  // Return whether triangle vertex has explicit normal
  return flags[R3_VERTEX_NORMALS_DRAW_FLAG];
}



inline RNBoolean R3TriangleVertex::
HasColor(void) const
{
  // Return whether triangle vertex has explicit color
  return flags[R3_VERTEX_COLORS_DRAW_FLAG];
}



inline RNBoolean R3TriangleVertex::
HasTextureCoords(void) const
{
  // Return whether triangle vertex has explicit texture coordinates
  return flags[R3_VERTEX_TEXTURE_COORDS_DRAW_FLAG];
}



inline void R3TriangleVertex::
SetPosition(const R3Point& point)
{
    // Reset vertex position
    position = point;
}



inline void R3TriangleVertex::
SetNormal(const R3Vector& vector)
{
    // Reset vertex normal
    normal = vector;

    // Update flags
    if (!normal.IsZero()) flags.Add(R3_VERTEX_NORMALS_DRAW_FLAG);
}



inline void R3TriangleVertex::
SetColor(const RNRgb& rgb)
{
    // Reset vertex color
    color = rgb;

    // Update flags
    flags.Add(R3_VERTEX_COLORS_DRAW_FLAG);
}



inline void R3TriangleVertex::
SetTextureCoords(const R2Point& point)
{
    // Reset vertex texture coordinates
    texcoords = point;

    // Update flags
    flags.Add(R3_VERTEX_TEXTURE_COORDS_DRAW_FLAG);
}



inline void R3TriangleVertex::
SetSharedFlag(void)
{
    // Set mark
    this->flags.Add(R3_VERTEX_SHARED_FLAG);
}



inline void R3TriangleVertex::
SetMark(RNMark mark)
{
    // Reset vertex mark
    this->mark = mark;
}



/* Triangle inline functions **************/

inline R3TriangleVertex *R3Triangle::
Vertex(int k) const
{
    // Return kth vertex
    assert((k >= 0) && (k < 3));
    return v[k];
}



inline R3TriangleVertex *R3Triangle::
V0(void) const
{
    // Return first vertex
    return v[0];
}



inline R3TriangleVertex *R3Triangle::
V1(void) const
{
    // Return second vertex
    return v[1];
}



inline R3TriangleVertex *R3Triangle::
V2(void) const
{
    // Return third vertex
    return v[2];
}



inline const R3Plane& R3Triangle::
Plane(void) const
{
    // Return plane containing triangle
    return plane;
}



inline const R3Vector& R3Triangle::
Normal(void) const
{
    // Return triangle normal
    return plane.Normal();
}



inline const R3Box& R3Triangle::
Box(void) const
{
    // Return bounding box 
    return bbox;
}



inline const RNFlags R3Triangle::
Flags(void) const
{
    // Return flags
    return flags;
}



inline const RNMark R3Triangle::
Mark(void) const
{
    // Return mark
    return mark;
}



inline const RNBoolean R3Triangle::
IsDegenerate(void) const
{
    // Return whether triangle is degenerate (zero area?)
    return IsLinear();
}



