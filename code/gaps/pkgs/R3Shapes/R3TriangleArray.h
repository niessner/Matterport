/* Include file for the R3 triangle class */



/* Initialization functions */

int R3InitTriangleArray();
void R3StopTriangleArray();



/* Triangle class definition */

class R3TriangleArray : public R3Surface {
    public:
        // Constructor functions
        R3TriangleArray(void);
        R3TriangleArray(const R3TriangleArray& array);
        R3TriangleArray(const RNArray<R3TriangleVertex *>& vertices, const RNArray<R3Triangle *>& triangles);
        virtual ~R3TriangleArray(void);
  
        // Triangle array properties
        const R3Box& Box(void) const;

	// Vertex access functions/operators
        int NVertices(void) const;
	R3TriangleVertex *Vertex(int index) const;

	// Triangle access functions/operators
        int NTriangles(void) const;
	R3Triangle *Triangle(int index) const;

        // Shape property functions/operators
	virtual const RNBoolean IsPoint(void) const;
	virtual const RNBoolean IsLinear(void) const;
	virtual const RNBoolean IsPlanar(void) const;
	virtual const RNBoolean IsConvex(void) const;
        virtual const RNInterval NFacets(void) const;
        virtual const RNLength Length(void) const;
        virtual const RNArea Area(void) const;
        virtual const R3Point Centroid(void) const;
        virtual const R3Point ClosestPoint(const R3Point& point) const;
        virtual const R3Point FurthestPoint(const R3Point& point) const;
        virtual const R3Shape& BShape(void) const;
        virtual const R3Box BBox(void) const;
        virtual const R3Sphere BSphere(void) const;

        // Manipulation functions/operators
	virtual void Flip(void);
	virtual void Mirror(const R3Plane& plane);
	virtual void Transform(const R3Transformation& transformation);
        virtual void Subdivide(RNLength max_edge_length);
	virtual void MoveVertex(R3TriangleVertex *vertex, const R3Point& position);
	virtual void Update(void);  

        // Draw functions/operators
        virtual void Draw(const R3DrawFlags draw_flags = R3_DEFAULT_DRAW_FLAGS) const;

	// Standard shape definitions
	RN_CLASS_TYPE_DECLARATIONS(R3TriangleArray);
        R3_SHAPE_RELATIONSHIP_DECLARATIONS(R3TriangleArray);

    private:
	RNArray<R3TriangleVertex *> vertices;
	RNArray<R3Triangle *> triangles;
        R3Box bbox;
};



/* Inline functions */

inline const R3Box& R3TriangleArray::
Box(void) const
{
    // Return bounding box
    return bbox;
}



inline int R3TriangleArray::
NVertices(void) const
{
    // Return number of vertices
    return vertices.NEntries();
}



inline R3TriangleVertex *R3TriangleArray::
Vertex(int k) const
{
    // Return kth vertex
    return vertices[k];
}



inline int R3TriangleArray::
NTriangles(void) const
{
    // Return number of triangles
    return triangles.NEntries();
}



inline R3Triangle *R3TriangleArray::
Triangle(int k) const
{
    // Return kth triangle
    return triangles[k];
}



