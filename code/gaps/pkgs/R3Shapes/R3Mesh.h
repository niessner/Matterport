// Include file for the R3 mesh class 



// Mesh element types

typedef enum {
  R3_MESH_NULL_TYPE,
  R3_MESH_VERTEX_TYPE,
  R3_MESH_EDGE_TYPE,
  R3_MESH_FACE_TYPE,
} R3MeshType;



// Mesh class declarations

class R3MeshVertex;
class R3MeshEdge;
class R3MeshFace;



// Mesh vertex definition

class R3MeshVertex {
  friend class R3Mesh;
  public:
    R3MeshVertex(void);
    virtual ~R3MeshVertex(void);
  protected:
    RNArray<R3MeshEdge *> edges;
    R3Point position;
    R3Vector normal;
    R2Point texcoords;
    RNRgb color;
    RNScalar curvature;
    int id;
    RNFlags flags;
    RNScalar value;
    RNMark mark;
    void *data;
};
  
  
  
// Mesh edge definition

class R3MeshEdge {
  friend class R3Mesh;
  public:
    R3MeshEdge(void);
    virtual ~R3MeshEdge(void);
  protected:
    class R3MeshVertex *vertex[2];
    class R3MeshFace *face[2];
    RNLength length;
    int id;
    RNFlags flags;
    RNScalar value;
    RNMark mark;
    void *data;
};
  


// Mesh face definition

class R3MeshFace {
  friend class R3Mesh;
  public:
    R3MeshFace(void);
    virtual ~R3MeshFace(void);
  protected:
    class R3MeshVertex *vertex[3];
    class R3MeshEdge *edge[3];
    R3Plane plane;
    R3Box bbox;
    RNArea area;
    int material;
    int segment;
    int category;
    int id;
    RNFlags flags;
    RNScalar value;
    RNMark mark;
    void *data;
};



// Mesh interesection definition

struct R3MeshIntersection {
  R3MeshType type;
  R3MeshVertex *vertex;
  R3MeshEdge *edge;
  R3MeshFace *face;
  R3Point point;
  RNScalar t;
};
  
  

// Useful constant definitions 

#define R3_MESH_NAME_LENGTH 128
#define RN_CCW 0
#define RN_CW  1 



// Mesh definition 

class R3Mesh {
  public:
    // CONSTRUCTOR FUNCTIONS
    R3Mesh(void);
      // Create an empty mesh
    R3Mesh(const R3Mesh& mesh);
      // Create a copy of a mesh
    virtual ~R3Mesh(void);
      // Deallocate all data associated with mesh

    // VERTEX/EDGE/FACE ENUMERATION FUNCTIONS
    int NVertices(void) const;
      // Return number of vertices in mesh
    R3MeshVertex *Vertex(int k) const;
      // Return kth vertex of mesh
    int NEdges(void) const;
      // Return number of edges in mesh
    R3MeshEdge *Edge(int k) const;
      // Return kth edge of mesh
    int NFaces(void) const;
      // Return number of faces in mesh
    R3MeshFace *Face(int k) const;
      // Return kth face of mesh

    // MESH PROPERTIES
    const char *Name(void) const;
      // Return name of mesh
    const R3Box& BBox(void) const;
      // Return bounding box of all vertices in mesh
    R3Point Centroid(void) const;
      // Return centroid of mesh
    RNArea Area(void) const;
      // Return surface area of mesh
    RNLength AverageEdgeLength(void) const;
      // Return average edge length
    RNLength AverageRadius(const R3Point *centroid = NULL) const;
      // Return average distance from point on surface to centroid
    R3Triad PrincipleAxes(const R3Point *centroid = NULL, RNScalar *variances = NULL) const;
      // Return principle axes of mesh
    R3Affine PCANormalizationTransformation(RNBoolean translate = TRUE, RNBoolean rotate = TRUE, RNBoolean scale = TRUE) const;
      // Return transformation that aligns centroid, principal axes, and scale with R3xyz_triad
    void *Data(void) const;
     // Return user data associated with mesh

    // VERTEX PROPERTIES
    const R3Point& VertexPosition(const R3MeshVertex *vertex) const;
      // Returns the position of the vertex
    const R3Vector& VertexNormal(const R3MeshVertex *vertex) const;
      // Returns the normal of the surface at the vertex
    const RNRgb& VertexColor(const R3MeshVertex *vertex) const;
      // Returns the color of the surface at the vertex
    const R2Point& VertexTextureCoords(const R3MeshVertex *vertex) const;
      // Returns the texture coordinates at the vertex
    RNArea VertexArea(const R3MeshVertex *vertex) const;
      // Returns the area of the surface associated with this vertex
    R3Plane VertexTangentPlane(const R3MeshVertex *vertex) const;
      // Returns the tangent plane at vertex 
    RNScalar VertexCurvature(const R3MeshVertex *vertex) const;
      // Returns the (Gauss) curvature of the surface at the vertex
    RNScalar VertexGaussCurvature(const R3MeshVertex *vertex) const;
      // Returns the Gauss curvature of the surface at the vertex
    RNScalar VertexMeanCurvature(const R3MeshVertex *vertex) const;
      // Returns the mean curvature of the surface at the vertex
    RNLength VertexAverageEdgeLength(const R3MeshVertex *vertex) const;
      // Returns the average length of edges attached to vertex
    R3Vector VertexLaplacianVector(const R3MeshVertex *vertex) const;
      // Returns the vector from a cotan weighted average of neighbors to the vertex
    int VertexValence(const R3MeshVertex *vertex) const;
      // Returns the number of edges attached to vertex
    int VertexID(const R3MeshVertex *vertex) const;
      // Returns the ID stored with a vertex
    RNFlags VertexFlags(const R3MeshVertex *vertex) const;
      // Returns the flags stored with a vertex
    RNScalar VertexValue(const R3MeshVertex *vertex) const;
      // Returns the value stored with a vertex
    RNMark VertexMark(const R3MeshVertex *vertex) const;
      // Returns the mark stored with a vertex
    void *VertexData(const R3MeshVertex *vertex) const;
      // Returns the user data stored with a vertex
  
    // EDGE PROPERTIES
    RNLength EdgeLength(const R3MeshEdge *edge) const;
      // Returns the length of the edge
    R3Point EdgeMidpoint(const R3MeshEdge *edge) const;
      // Returns the midpoint of the edge
    RNAngle EdgeInteriorAngle(const R3MeshEdge *edge) const;
      // Returns interior angle between adjacent faces in range 0 to 2*PI (returns 0 if boundary edge)
    R3Vector EdgeNormal(const R3MeshEdge *edge) const;
      // Returns the normal of the surface at the edge
    RNScalar EdgeAspect(const R3MeshEdge *edge) const;
      // Returns ratio of "distance between vertices across faces" over edge length
    R3Vector EdgeVector(const R3MeshEdge *edge) const;
      // Returns a vector along the edge (v2 - v1)
    R3Vector EdgeDirection(const R3MeshEdge *edge) const;
      // Returns a normalized vector along edge (v2 - v1)/||v2 - v1||
    R3Span EdgeSpan(const R3MeshEdge *edge) const;
      // Returns the span covering the edge
    R3Line EdgeLine(const R3MeshEdge *edge) const;
      // Returns the line supporting the edge
    R3Box EdgeBBox(const R3MeshEdge *edge) const;
      // Returns the bounding box of the edge
    int EdgeID(const R3MeshEdge *edge) const;
      // Returns the ID stored with a edge
    RNFlags EdgeFlags(const R3MeshEdge *edge) const;
      // Returns the flags stored with the edge
    RNScalar EdgeValue(const R3MeshEdge *edge) const;
      // Returns the value stored with a edge
    RNMark EdgeMark(const R3MeshEdge *edge) const;
      // Returns the mark stored with the edge
    void *EdgeData(const R3MeshEdge *edge) const;
      // Returns the user data stored with the edge

    // FACE PROPERTIES
    const R3Vector& FaceNormal(const R3MeshFace *face) const;
      // Returns the normal of the face
    const R3Plane& FacePlane(const R3MeshFace *face) const;
      // Returns the plane equation of the face
    R3Point FaceCentroid(const R3MeshFace *face) const;
      // Returns the centroid of the face
    R3Point FacePoint(const R3MeshFace *face, RNMagnitude barycentrics[3]) const;
      // Returns the point on the face with the given barycentric coordinates
    R3Point FaceBarycentric(const R3MeshFace *face, const R3Point& point) const;
      // Returns the barycentric coordinates of point on face
    const R3Box& FaceBBox(const R3MeshFace *face) const;
      // Returns the bounding box of the face
    RNArea FaceArea(const R3MeshFace *face) const;
      // Returns the area of the face
    RNScalar FaceAspect(const R3MeshFace *face) const;
      // Returns the ratio of the circumradius to twice its inradius
    int FaceMaterial(const R3MeshFace *face) const;
      // Returns the material ID associated with a face
    int FaceSegment(const R3MeshFace *face) const;
      // Returns the segment ID associated with a face
    int FaceCategory(const R3MeshFace *face) const;
      // Returns the category ID associated with a face
    int FaceID(const R3MeshFace *face) const;
      // Returns the ID stored with a face
    RNFlags FaceFlags(const R3MeshFace *face) const;
      // Returns the flags stored with a face
    RNScalar FaceValue(const R3MeshFace *face) const;
      // Returns the value stored with a face
    RNMark FaceMark(const R3MeshFace *face) const;
      // Returns the mark stored with a face
    void *FaceData(const R3MeshFace *face) const;
      // Returns the user data stored with a face
  
    // MESH PROPOERTY SETTING FUNCTIONS
    void SetName(const char *name);
      // Set mesh name

    // VERTEX PROPOERTY SETTING FUNCTIONS
    void SetVertexPosition(R3MeshVertex *vertex, const R3Point& position);
      // Set new position for vertex
    void SetVertexNormal(R3MeshVertex *vertex, const R3Vector& normal);
      // Set new normal for vertex
    void SetVertexColor(R3MeshVertex *vertex, const RNRgb& color);
      // Set new color for vertex
    void SetVertexTextureCoords(R3MeshVertex *vertex, const R2Point& texcoords);
      // Set new texture coordinates for vertex
    void SetVertexValue(R3MeshVertex *vertex, RNScalar value);
      // Set the value stored with a vertex
    void SetVertexMark(R3MeshVertex *vertex, RNMark mark);
      // Set the mark stored with a vertex
    void SetVertexData(R3MeshVertex *vertex, void *data);
      // Set the data stored with a vertex
    void SetVertexFlags(R3MeshVertex *vertex, RNFlags flags);
      // Set the flags stored with a vertex
  
    // EDGE PROPOERTY SETTING FUNCTIONS
    void SetEdgeMark(R3MeshEdge *edge, RNMark mark);
      // Set the mark stored with a edge
    void SetEdgeValue(R3MeshEdge *edge, RNScalar value);
      // Set the value stored with a edge
    void SetEdgeData(R3MeshEdge *edge, void *data);
      // Set the data stored with a edge
    void SetEdgeFlags(R3MeshEdge *edge, RNFlags flags);
      // Set the flags stored with a edge
  
    // FACE PROPOERTY SETTING FUNCTIONS
    void SetFacePlane(R3MeshFace *face, const R3Plane& plane);
      // Set the plane equation for a face
    void SetFaceMaterial(R3MeshFace *face, int material);
      // Set the material associated with a face
    void SetFaceSegment(R3MeshFace *face, int segment);
      // Set the segment associated with a face
    void SetFaceCategory(R3MeshFace *face, int category);
      // Set the category associated with a face
    void SetFaceMark(R3MeshFace *face, RNMark mark);
      // Set the mark stored with a face
    void SetFaceValue(R3MeshFace *face, RNScalar value);
      // Set the value stored with a face
    void SetFaceData(R3MeshFace *face, void *data);
      // Set the data stored with a face
    void SetFaceFlags(R3MeshFace *face, RNFlags flags);
      // Set the flags stored with a face
  
    // VERTEX TRAVERSAL FUNCTIONS
    R3MeshVertex *VertexOnVertex(const R3MeshVertex *vertex) const;
      // Returns any vertex adjacent to vertex
    R3MeshVertex *VertexOnVertex(const R3MeshVertex *vertex, int k) const;
      // Returns kth vertex adjacent to vertex (neighbor vertices are not sorted!)
    R3MeshEdge *EdgeOnVertex(const R3MeshVertex *vertex) const;
      // Returns any edge connected to vertex
    R3MeshEdge *EdgeOnVertex(const R3MeshVertex *vertex, int k) const;
      // Returns kth edge connected to vertex (edges are not sorted!)
    R3MeshEdge *EdgeBetweenVertices(const R3MeshVertex *vertex1, const R3MeshVertex *vertex2) const;
      // Returns edge spanning between two vertices (or NULL if none)
    R3MeshEdge *EdgeOnVertex(const R3MeshVertex *vertex, const R3MeshEdge *edge, RNDirection dir = RN_CCW) const;
      // Returns edge on vertex in dir with respect to edge
    R3MeshEdge *EdgeOnVertex(const R3MeshVertex *vertex, const R3MeshFace *face, RNDirection dir = RN_CCW) const;
      // Returns edge on vertex in dir with respect to face
    R3MeshEdge *EdgeAcrossVertex(const R3MeshVertex *vertex, const R3MeshEdge *edge, const R3MeshFace *face) const;
      // Returns edge on the other side of a vertex from an edge on the same face 
    R3MeshFace *FaceOnVertex(const R3MeshVertex *vertex) const;
      // Returns any face connected to vertex
    R3MeshFace *FaceBetweenVertices(const R3MeshVertex *vertex1, const R3MeshVertex *vertex2, RNDirection dir = RN_CCW) const;
      // Returns face between two vertices (in dir with respect to v1) or NULL if there is none
    R3MeshFace *FaceOnVertex(const R3MeshVertex *vertex, const R3MeshEdge *edge, RNDirection dir = RN_CCW) const;
      // Returns face on vertex in dir with respect to edge
    R3MeshFace *FaceOnVertex(const R3MeshVertex *vertex, const R3MeshFace *face, RNDirection dir = RN_CCW) const;
      // Returns face on vertex in dir with respect to face
  
    // EDGE TRAVERSAL FUNCTIONS
    R3MeshVertex *VertexOnEdge(const R3MeshEdge *edge) const;
      // Returns any vertex on edge 
    R3MeshVertex *VertexOnEdge(const R3MeshEdge *edge, int k) const;
      // Returns kth vertex on edge 
    R3MeshVertex *VertexAcrossEdge(const R3MeshEdge *edge, const R3MeshVertex *vertex) const;
      // Returns vertex across edge 
    R3MeshVertex *VertexBetweenEdges(const R3MeshEdge *edge1, const R3MeshEdge *edge2) const;
      // Returns vertex between two edges (or NULL if there is none)
    R3MeshVertex *VertexOnEdge(const R3MeshEdge *edge, const R3MeshFace *face, RNDirection dir = RN_CCW) const;
      // Returns vertex of edge in dir with respect to face
    R3MeshEdge *EdgeOnEdge(const R3MeshEdge *edge) const;
      // Returns any edge on edge 
    R3MeshEdge *EdgeOnEdge(const R3MeshEdge *edge, const R3MeshVertex *vertex, RNDirection dir = RN_CCW) const;
      // Returns edge on edge in dir with respect to vertex
    R3MeshFace *FaceOnEdge(const R3MeshEdge *edge) const;
      // Returns any face on edge 
    R3MeshFace *FaceOnEdge(const R3MeshEdge *edge, int k) const;
      // Returns kth face on edge 
    R3MeshFace *FaceOnEdge(const R3MeshEdge *edge, const R3MeshVertex *vertex, RNDirection dir = RN_CCW) const;
      // Returns face on edge in dir with respect to vertex
    R3MeshFace *FaceBetweenEdges(const R3MeshEdge *edge1, const R3MeshEdge *edge2) const;
    // Returns the face between two edges (or NULL if there is none)
    R3MeshFace *FaceAcrossEdge(const R3MeshEdge *edge, const R3MeshFace *face) const;
      // Returns face across edge 
  
    // FACE TRAVERSAL FUNCTIONS
    R3MeshVertex *VertexOnFace(const R3MeshFace *face) const;
      // Returns any vertex on face
    R3MeshVertex *VertexOnFace(const R3MeshFace *face, int k) const;
      // Returns kth vertex on face
    R3MeshVertex *VertexOnFace(const R3MeshFace *face, const R3MeshVertex *vertex, RNDirection dir = RN_CCW) const;
      // Returns vertex on face in dir with respect to vertex
    R3MeshVertex *VertexOnFace(const R3MeshFace *face, const R3MeshEdge *edge, RNDirection dir = RN_CCW) const;
      // Returns vertex on face in dir with respect to edge
    R3MeshVertex *VertexAcrossFace(const R3MeshFace *face, const R3MeshEdge *edge) const;
      // Returns the vertex across a triangle from an edge (the edge is the base and the vertex is the apex)
    R3MeshVertex *VertexBetweenFaces(const R3MeshFace *face1, const R3MeshFace *face2, RNDirection dir = RN_CCW) const;
      // Returns the vertex between two faces (in dir with respect to face1) or NULL if there is none
    R3MeshEdge *EdgeOnFace(const R3MeshFace *face) const;
      // Returns any edge connected to face
    R3MeshEdge *EdgeOnFace(const R3MeshFace *face, int k) const;
      // Returns kth edge connected to face
    R3MeshEdge *EdgeOnFace(const R3MeshFace *face, const R3MeshVertex *vertex, RNDirection dir = RN_CCW) const;
      // Returns edge on face in dir with respect to vertex
    R3MeshEdge *EdgeOnFace(const R3MeshFace *face, const R3MeshEdge *edge, RNDirection dir = RN_CCW) const;
      // Return edge on face in dir with respect to edge
    R3MeshEdge *EdgeAcrossFace(const R3MeshFace *face, const R3MeshVertex *vertex) const;
      // Returns the edge across a triangle from a vertex (the edge is the base and the vertex is the apex)
    R3MeshEdge *EdgeBetweenFaces(const R3MeshFace *face1, const R3MeshFace *face2) const;
      // Returns the edge between two faces (or NULL if there is none)
    R3MeshFace *FaceOnFace(const R3MeshFace *face) const;
      // Returns any face connected to face
    R3MeshFace *FaceOnFace(const R3MeshFace *face, int k) const;
      // Returns kth face connected to face
  
    // TOPOLOGY QUERY FUNCTIONS
    RNBoolean IsVertexOnEdge(const R3MeshVertex *vertex, const R3MeshEdge *edge) const;
      // Returns whether vertex lies on edge
    RNBoolean IsVertexOnFace(const R3MeshVertex *vertex, const R3MeshFace *face) const;
      // Returns whether vertex lies on face
    RNBoolean IsVertexOnBoundary(const R3MeshVertex *vertex) const;
      // Returns whether vertex lies on mesh boundary
    RNBoolean IsVertexOnMesh(const R3MeshVertex *vertex) const;
      // Returns whether vertex is part of mesh
    RNBoolean IsEdgeOnFace(const R3MeshEdge *edge, const R3MeshFace *face) const;
      // Returns whether edge lies on face
    RNBoolean IsEdgeOnBoundary(const R3MeshEdge *edge) const;
      // Returns whether edge lies on mesh boundary
    RNBoolean IsEdgeOnMesh(const R3MeshEdge *edge) const;
      // Returns whether edge is part of mesh
    RNBoolean IsFaceOnBoundary(const R3MeshFace *face) const;
      // Returns whether face has any edge on mesh boundary
    RNBoolean IsFaceOnMesh(const R3MeshFace *face) const;
      // Returns whether face is part of mesh
    void FindConnectedVertices(R3MeshVertex *seed, RNArray<R3MeshVertex *>& vertices);
      // Fill array with vertices in same connected component as seed
    void FindConnectedFaces(R3MeshFace *seed, RNArray<R3MeshFace *>& faces);
      // Fill array with faces in same connected component as seed
  
    // TOPOLOGY MANIPULATION FUNCTIONS
    virtual void Empty(void);
      // Delete all vertices, edges, faces
    virtual R3MeshVertex *CreateVertex(const R3MeshVertex& source, R3MeshVertex *vertex = NULL);
      // Create a new vertex with position, normal, color and texture coordinates of source
    virtual R3MeshVertex *CreateVertex(const R3Point& position, R3MeshVertex *vertex = NULL);
      // Create a new vertex at a position (compute normal)
    virtual R3MeshVertex *CreateVertex(const R3Point& position, const R3Vector& normal, R3MeshVertex *vertex = NULL);
      // Create a new vertex at a position with a given normal
    virtual R3MeshVertex *CreateVertex(const R3Point& position, const R3Vector& normal, const RNRgb& color, R3MeshVertex *vertex = NULL);
      // Create a new vertex at a position with a given normal and color
    virtual R3MeshVertex *CreateVertex(const R3Point& position, const R3Vector& normal, const RNRgb& color, const R2Point& texcoords, R3MeshVertex *vertex = NULL);
      // Create a new vertex at a position with a given normal and color and texture coordinates
    virtual R3MeshEdge *CreateEdge(R3MeshVertex *v1, R3MeshVertex *v2, R3MeshEdge *edge = NULL);
      // Create a new edge between two vertices
    virtual R3MeshFace *CreateFace(R3MeshVertex *v1, R3MeshVertex *v2, R3MeshVertex *v3, R3MeshFace *face = NULL);
      // Create a face with vertices (v1, v2, v3)
    virtual R3MeshFace *CreateFace(R3MeshEdge *e1, R3MeshEdge *e2, R3MeshEdge *e3, R3MeshFace *face = NULL);
      // Create a face with edges (e1, e2, e3)
    virtual R3MeshFace *CreateFace(R3MeshVertex *v1, R3MeshVertex *v2, R3MeshVertex *v3, R3MeshEdge *e1, R3MeshEdge *e2, R3MeshEdge *e3, R3MeshFace *face = NULL);
      // Create a face with vertices (v1, v2, v3) and edges (e1, e2, e3)
    virtual void DeleteVertex(R3MeshVertex *v);
      // Delete a vertex and all edges/faces attached to it
    virtual void DeleteEdge(R3MeshEdge *e);
      // Delete an edge (attached vertices are unaffected, attached faces are deleted)
    virtual void DeleteFace(R3MeshFace *f);
      // Delete a face (attached vertices and edges are unaffected)
    virtual void DeleteUnusedVertices(void);
      // Delete vertices not attached to any edge
    virtual void DeleteUnusedEdges(void);
      // Delete edges not attached to any face
    virtual R3MeshVertex *MergeVertex(R3MeshVertex *v1, R3MeshVertex *v2);
      // Merge vertex v2 into vertex v1 (v2 is deleted)
    virtual void MergeCoincidentVertices(RNLength epsilon = -1.0);
      // Merge all vertices within epsilon of each other (negative epsilon says "select epsilon automatically")
    virtual R3MeshVertex *CollapseEdge(R3MeshEdge *edge, const R3Point& point);
      // Collapses an edge into a vertex at point
    virtual R3MeshVertex *CollapseFace(R3MeshFace *face, const R3Point& point);
      // Collapses a face into a vertex at point
    virtual R3MeshVertex *CollapseEdge(R3MeshEdge *edge);
      // Collapses an edge into a vertex at point
    virtual R3MeshVertex *CollapseFace(R3MeshFace *face);
      // Collapses a face into a vertex at point
    virtual R3MeshVertex *SplitEdge(R3MeshEdge *edge, const R3Point& point, R3MeshEdge **e0 = NULL, R3MeshEdge **e1 = NULL);
      // Creates a vertex at point, and inserts it into an edge, returns new vertex and edges (in e0-e1)
    virtual R3MeshVertex *SplitEdge(R3MeshEdge *edge, const R3Plane& plane);
      // Creates a vertex at intersection with plane, and inserts it into an edge
    virtual R3MeshVertex *SubdivideEdge(R3MeshEdge *edge);
      // Splits edge at midpoint
    virtual int SwapEdge(R3MeshEdge *edge);
      // Swap edge to use vertices across from adjacent faces
    virtual void FlipEdge(R3MeshEdge *edge);
      // Reverse direction of edge (not important for most applications)
    virtual void FlipFace(R3MeshFace *face);
      // Reverses orientation of face
    virtual R3MeshVertex *SplitFace(R3MeshFace *face, const R3Point& point, 
      R3MeshFace **f0 = NULL, R3MeshFace **f1 = NULL, R3MeshFace **f2 = NULL);
      // Creates a new vertex at point, removes face, and builds new faces between edges and new vertices, returns new vertex and faces (in f0-f2)
    virtual R3MeshFace *SubdivideFace(R3MeshFace *face);
      // Splits face into four by subdividing each edge at midpoint (returns middle face)
    void CollapseShortEdges(RNLength min_edge_length);
      // Collapse edges shorter than min_edge_length
    void SubdivideLongEdges(RNLength max_edge_length);
      // Subdivide edges longer than max_edge_length
    void FlipFaces(void);
      // Flip orientation of all faces
    void SubdivideFaces(void);
      // Subdivide all faces, replacing each triangle by four with vertices at edge midpoints
    void SwapEdges(RNAngle min_angle_improvement = 0.1);
     // Swap edges to improve increase interior angles of faces
    void FillHole(R3MeshEdge *boundary_edge);
     // Fill hole bounded by the given edge
    void FillHoles(void);
     // Fill all holes in mesh

    // GEOMETRY MANIPULATION FUNCTIONS
    virtual void Smooth(RNScalar factor = 1.0);
      // Smooth mesh
    virtual void Sharpen(RNScalar factor = 1.0);
      // Sharpen mesh 
    virtual void AddRandomNoise(RNScalar factor);
      // Add noise to vertex positions (factor is relative to average edge length at each vertex)
    virtual void Inflate(RNScalar factor);
      // Move every vertex along normal vector by factor times average edge length at vertex
    virtual void Transform(const R3Transformation& transformation);
      // Transform mesh
 
    // SHAPE CREATION FUNCTIONS
    void CreateBox(const R3Box& box);
      // Create mesh elements for a box and add to mesh
    void CreateOrientedBox(const R3OrientedBox& box);
      // Create mesh elements for an oriented box and add to mesh
    void CreateSphere(const R3Sphere& sphere, RNLength vertex_spacing = 0);
      // Create mesh elements for a sphere and add to mesh
    void CreateEllipsoid(const R3Ellipsoid& ellipsoid, RNLength vertex_spacing = 0);
      // Create mesh elements for a ellipsoid and add to mesh
    void CreateCone(const R3Cone& cone, RNLength vertex_spacing = 0);
      // Create mesh elements for a cone and add to mesh
    void CreateCylinder(const R3Cylinder& cylinder, RNLength vertex_spacing = 0);
      // Create mesh elements for a cylinder and add to mesh
    void CreateCircle(const R3Circle& circle, RNLength vertex_spacing = 0);
      // Create mesh elements for a circle and add to mesh
    void CreateEllipse(const R3Ellipse& ellipse, RNLength vertex_spacing = 0);
      // Create mesh elements for a ellipse and add to mesh
    void CreateRectangle(const R3Rectangle& rectangle);
      // Create mesh elements for a rectangle and add to mesh
    void CreateTriangle(const R3Triangle& triangle);
      // Create mesh elements for a triangle and add to mesh
    void CreateTriangleArray(const R3TriangleArray& triangles);
      // Create mesh elements for a triangle array and add to mesh
    void CreateCopy(const R3Mesh& mesh);
      // Create mesh elements for a copy of mesh

    // DRAW FUNCTIONS
    virtual void Draw(void) const;
      // Draws the faces
    virtual void DrawVertices(void) const;
      // Draws the vertices
    virtual void DrawEdges(void) const;
      // Draws the edges
    virtual void DrawFaces(void) const;
      // Draws the faces
    virtual void DrawVertexIDs(void) const;
      // Draws the vertex IDs into color buffer
    virtual void DrawEdgeIDs(void) const;
      // Draws the edge IDs into color buffer
    virtual void DrawFaceIDs(void) const;
      // Draws the face IDs into color buffer
    virtual void DrawVertex(R3MeshVertex *vertex) const;
      // Draws one vertex
    virtual void DrawEdge(R3MeshEdge *edge) const;
      // Draws one edge
    virtual void DrawFace(R3MeshFace *face) const;
      // Draws one face

    // INTERSECTION FUNCTIONS
    virtual R3MeshType Intersection(const R3Ray& ray, R3MeshIntersection *intersection = NULL) const;
      // Returns first intersection along ray 
    virtual R3MeshType Intersection(const R3Ray& ray, R3MeshFace *face, R3MeshIntersection *intersection = NULL) const;
      // Returns intersection with face along ray 
  
    // CLOSEST POINT FUNCTIONS
    R3Point ClosestPoint(const R3Point& point, R3MeshIntersection *closest_point = NULL) const;
      // Returns closest point on mesh
    R3Point ClosestPointOnEdge(const R3MeshEdge *edge, const R3Point& point, R3MeshIntersection *closest_point = NULL) const;
      // Returns closest point on edge
    R3Point ClosestPointOnFace(const R3MeshFace *face, const R3Point& point, R3MeshIntersection *closest_point = NULL) const;
      // Returns closest point on face

    // SURFACE DISTANCE FUNCTIONS
    RNLength *GeodesicDistances(const R3MeshVertex *source_vertex, RNLength max_distance = 0) const;
      // Returns array of geodesic distances from source vertex to all other vertices (return array is indexed by VertexID).
    RNLength DijkstraDistance(const R3MeshVertex *source_vertex, const R3MeshVertex *destination_vertex, RNArray<R3MeshEdge *> *edges = NULL) const;
      // Returns approximate geodesic distance from source_vertex to destination_vertex.  
      // If "edges" is non-NULL (it should have been created as an empty RNArray<R3MeshEdge *>), 
      // then it will be filled in with the edges in the shortest path from the destination_vertex back to the source_vertex
    RNLength *DijkstraDistances(const R3MeshVertex *source_vertex, RNLength max_distance = 0, R3MeshEdge **edges = NULL) const;
      // Returns array of approximate geodesic distances from source vertex to all other vertices (return array is indexed by VertexID).
    RNLength *DijkstraDistances(const RNArray<R3MeshVertex *>& source_vertices, RNLength max_distance = 0, R3MeshEdge **edges = NULL) const;
      // Returns array of approximate geodesic distances from closest source vertex to all other vertices (return array is indexed by VertexID).
      // If max_distance is non-zero, then accurate distances will be computed only up to that threshold
      // If "edges" is non-NULL (it should have already been allocated with size NVertices), 
      // then it will be filled in with the edge from each vertex (indexed by VertexID) to its ancestor in the 
      // shortest spanning tree back to the source_vertex
    RNLength *DijkstraDistances(const R3MeshVertex *source_vertex, RNLength max_distance, RNArray<R3MeshVertex *>& vertices, R3MeshEdge **edges = NULL) const;
      // Returns array of approximate geodesic distances from source vertex to other vertices within max_distance (return array is indexed by VertexID).
      // The "vertices" RNArray is filled with vertices withing the max_distance threshold, and only those vertices have a valid distance in the return array
    RNLength *DijkstraDistances(const RNArray<R3MeshVertex *>& source_vertices, RNLength max_distance, RNArray<R3MeshVertex *>& vertices, R3MeshEdge **edges = NULL) const;
      // Returns array of approximate geodesic distances from closest source vertex to other vertices within max_distance (return array is indexed by VertexID).
      // The "vertices" RNArray is filled with vertices withing the max_distance threshold, and only those vertices have a valid distance in the return array
      // If "edges" is non-NULL (it should have already been allocated with size NVertices),
      // then it will be filled in with the edge from each vertex (indexed by VertexID) to its ancestor in the
      // shortest spanning tree back to the source_vertex
    RNLength TracePath(const R3MeshVertex *source_vertex, const R3Vector& source_direction, RNLength max_distance,
      R3Point *position = NULL, R3MeshFace **face = NULL, R3MeshIntersection *intersections = NULL, int *nintersections = NULL) const;
      // Traces a path from the source_vertex starting in the source_direction along the mesh surface for max_distance.
      // The position and its face found at the end of the path are returned if the given addresses are non-NULL
      // Also, the sequence and number of edge intersections is filled into intersections and nintersections 
      // (sufficient memory must be allocated for intersections by the caller if this feature is used)

    // POINT SAMPLING FUNCTIONS
    R3Point RandomPointOnFace(const R3MeshFace *face) const;
      // Returns random point on face

    // I/O FUNCTIONS
    virtual int ReadFile(const char *filename);
      // Loads data structure from a file, returns 0 if error
    virtual int ReadObjFile(const char *filename);
      // Loads data structure from Wavefront file (.obj), returns 0 if error
    virtual int ReadOffFile(const char *filename);
      // Loads data structure from OFF (.off) file, returns 0 if error
    virtual int ReadRayFile(const char *filename);
      // Loads data structure from ray file (.ray), returns 0 if error
    virtual int ReadPlyFile(const char *filename);
      // Loads data structure from ply file (.ply), returns 0 if error
    virtual int ReadCattFile(const char *filename);
      // Loads data structure from a Catt acoustic file (.cat), returns 0 if error
    virtual int ReadHoppeFile(const char *filename);
      // Loads data structure from a file in Hugues' mesh format (.m), returns 0 if error
    virtual int ReadIfsFile(const char *filename);
      // Loads data structure from ifs file (.ifs), returns 0 if error
    virtual int ReadSTLFile(const char *filename);
      // Loads data structure from STL file (.stl), returns 0 if error
    virtual int ReadVRMLFile(const char *filename);
      // Loads data structure from VRML file (.wrl), returns 0 if error
    virtual int WriteFile(const char *filename) const;
     // Writes a file, returns number of faces written (0 is error)
    virtual int WriteRayFile(const char *filename) const;
     // Writes a ray file, returns number of faces written (0 is error)
    virtual int WritePlyFile(const char *filename, RNBoolean binary = TRUE) const;
     // Writes a PLY file, returns number of faces written (0 is error)
    virtual int WriteObjFile(const char *filename) const;
     // Writes a Wavefront OBJ file, returns number of faces written (0 is error)
    virtual int WriteOffFile(const char *filename) const;
     // Writes a OFF (.off) file, returns number of faces written (0 is error)
    virtual int WriteCattFile(const char *filename) const;
     // Writes a Catt Acoustics (.cat) file, returns number of faces written (0 is error)
    virtual int WriteIfsFile(const char *filename) const;
     // Writes a IFS (.ifs) file, returns number of faces written (0 is error)
    virtual int WriteSTLFile(const char *filename) const;
     // Writes a STL (.stl) file, returns number of faces written (0 is error)

    // DEBUG FUNCTIONS
    virtual RNBoolean IsValid(void) const;
      // Returns whether data structure is valid
    virtual RNBoolean IsValid(R3MeshVertex *vertex) const;
      // Returns whether vertex data structure is valid
    virtual RNBoolean IsValid(R3MeshEdge *edge) const;
      // Returns whether edge data structure is valid
    virtual RNBoolean IsValid(R3MeshFace *face) const;
      // Returns whether face data structure is valid
    void SetData(void *data);
      // Set the user data stored with mesh
  
  protected:
    // INTERNAL DELETE FUNCTIONS
    virtual void DeallocateVertex(R3MeshVertex *v);
    virtual void DeallocateEdge(R3MeshEdge *e);
    virtual void DeallocateFace(R3MeshFace *f);

    // INTERNAL UPDATE FUNCTIONS
    virtual void UpdateVertexNormal(R3MeshVertex *v) const;  
    virtual void UpdateVertexCurvature(R3MeshVertex *v) const;  
    virtual void UpdateEdgeLength(R3MeshEdge *e) const;  
    virtual void UpdateFaceArea(R3MeshFace *f) const;
    virtual void UpdateFacePlane(R3MeshFace *f) const;  
    virtual void UpdateFaceBBox(R3MeshFace *f) const;  
    virtual void UpdateFaceRefs(R3MeshFace *f, R3MeshVertex *v1, R3MeshVertex *v2, R3MeshVertex *v3,
                                R3MeshEdge *e1, R3MeshEdge *e2, R3MeshEdge *e3);
  
  protected:
    // Arrays of all vertices, edges, faces
    RNArray<R3MeshVertex *> vertices;
    RNArray<R3MeshEdge *> edges;
    RNArray<R3MeshFace *> faces;

    // Storage (if allocated stuff in blocks)
    R3MeshVertex *vertex_block;
    R3MeshEdge *edge_block;
    R3MeshFace *face_block;

    // Other attributes
    char name[R3_MESH_NAME_LENGTH];
    R3Box bbox;
    void *data;
};



////////////////////////////////////////////////////////////////////////
// Flag definitions
////////////////////////////////////////////////////////////////////////

#define R3_MESH_BBOX_UPTODATE             1
#define R3_MESH_VERTEX_ALLOCATED          1
#define R3_MESH_VERTEX_NORMAL_UPTODATE    2
#define R3_MESH_VERTEX_CURVATURE_UPTODATE 4
#define R3_MESH_VERTEX_USER_FLAG          8
#define R3_MESH_EDGE_ALLOCATED            1
#define R3_MESH_EDGE_LENGTH_UPTODATE      2
#define R3_MESH_EDGE_USER_FLAG            4
#define R3_MESH_FACE_ALLOCATED            1
#define R3_MESH_FACE_AREA_UPTODATE        2
#define R3_MESH_FACE_PLANE_UPTODATE       4
#define R3_MESH_FACE_BBOX_UPTODATE        8
#define R3_MESH_FACE_USER_FLAG           16



////////////////////////////////////////////////////////////////////////
// Public variables
////////////////////////////////////////////////////////////////////////

extern RNMark R3mesh_mark;



////////////////////////////////////////////////////////////////////////
// VERTEX/EDGE/FACE ENUMERATION FUNCTIONS
////////////////////////////////////////////////////////////////////////

inline int R3Mesh::
NVertices(void) const
{
  // Return number of vertices
  return vertices.NEntries();
}



inline R3MeshVertex *R3Mesh::
Vertex(int k) const
{
  // Return kth vertex
  return vertices[k];
}



inline int R3Mesh::
NEdges(void) const
{
  // Return number of edges
  return edges.NEntries();
}



inline R3MeshEdge *R3Mesh::
Edge(int k) const
{
  // Return kth edge
  return edges[k];
}



inline int R3Mesh::
NFaces(void) const
{
  // Return number of faces
  return faces.NEntries();
}



inline R3MeshFace *R3Mesh::
Face(int k) const
{
  // Return kth face
  return faces[k];
}



////////////////////////////////////////////////////////////////////////
// MESH PROPERTY FUNCTIONS
////////////////////////////////////////////////////////////////////////

inline const char *R3Mesh::
Name(void) const
{
  // Return name
  return name;
}



inline const R3Box& R3Mesh::
BBox(void) const
{
  // Return bounding box
  return bbox;
}



inline void *R3Mesh:: 
Data(void) const
{
  // Returns the user data stored with the mesh
  return data; 
}



////////////////////////////////////////////////////////////////////////
// VERTEX PROPERTY FUNCTIONS
////////////////////////////////////////////////////////////////////////

inline const R3Point& R3Mesh:: 
VertexPosition(const R3MeshVertex *v) const
{ 
  // Return the position of the vertex
  return v->position;
}



inline const R3Vector& R3Mesh::
VertexNormal(const R3MeshVertex *v) const
{
  // Update the vertex normal
  if (!(v->flags[R3_MESH_VERTEX_NORMAL_UPTODATE]))
    UpdateVertexNormal((R3MeshVertex *) v);

  // Return the normal at the vertex
  return v->normal;
}



inline const RNRgb& R3Mesh::
VertexColor(const R3MeshVertex *v) const
{
  // Return the color at the vertex
  return v->color;
}



inline const R2Point& R3Mesh::
VertexTextureCoords(const R3MeshVertex *v) const
{
  // Return the texture coordinates at the vertex
  return v->texcoords;
}



inline R3Plane R3Mesh::
VertexTangentPlane(const R3MeshVertex *v) const
{
  // Return tangent plane to surface at vertex
  return R3Plane(VertexPosition(v), VertexNormal(v));
}



inline RNScalar R3Mesh::
VertexCurvature(const R3MeshVertex *v) const
{
  // Update the vertex curvature
  if (!(v->flags[R3_MESH_VERTEX_CURVATURE_UPTODATE]))
    UpdateVertexCurvature((R3MeshVertex *) v);

  // Return the curvature at the vertex
  return v->curvature;
}



inline int R3Mesh:: 
VertexValence(const R3MeshVertex *v) const
{
  // Returns the number of edges attached to vertex
  return v->edges.NEntries(); 
}



inline int R3Mesh:: 
VertexID(const R3MeshVertex *v) const
{
  // Returns the ID stored with a vertex
  return v->id; 
}



inline RNFlags R3Mesh:: 
VertexFlags(const R3MeshVertex *v) const
{
  // Returns the flags stored with a vertex
  return v->flags; 
}



inline RNScalar R3Mesh:: 
VertexValue(const R3MeshVertex *v) const
{
  // Returns the value stored with a vertex
  return v->value; 
}



inline RNMark R3Mesh:: 
VertexMark(const R3MeshVertex *v) const
{
  // Returns the mark stored with a vertex
  return v->mark; 
}



inline void *R3Mesh:: 
VertexData(const R3MeshVertex *v) const
{
  // Returns the data stored with a vertex
  return v->data; 
}



////////////////////////////////////////////////////////////////////////
// EDGE PROPERTY FUNCTIONS
////////////////////////////////////////////////////////////////////////

inline RNLength R3Mesh::
EdgeLength(const R3MeshEdge *e) const
{
  // Update the edge length
  if (!(e->flags[R3_MESH_EDGE_LENGTH_UPTODATE]))
    UpdateEdgeLength((R3MeshEdge *) e);

  // Return the edge length
  return e->length;
}



inline R3Point R3Mesh:: 
EdgeMidpoint(const R3MeshEdge *e) const
{
  // Returns the midpoint of the edge
  assert(e->vertex[0] && e->vertex[1]);
  return 0.5 * (e->vertex[0]->position + e->vertex[1]->position);
}



inline R3Vector R3Mesh:: 
EdgeVector(const R3MeshEdge *e) const
{
  // Returns the vector pointing in direction of the edge
  assert(e->vertex[0] && e->vertex[1]);
  return e->vertex[1]->position - e->vertex[0]->position;
}



inline R3Span R3Mesh:: 
EdgeSpan(const R3MeshEdge *e) const
{
  // Returns the span covering the edge
  assert(e->vertex[0] && e->vertex[1]);
  return R3Span(e->vertex[0]->position, e->vertex[1]->position);
}



inline R3Line R3Mesh:: 
EdgeLine(const R3MeshEdge *e) const
{
  // Returns the line supporting the edge
  assert(e->vertex[0] && e->vertex[1]);
  return R3Line(e->vertex[0]->position, e->vertex[1]->position);
}



inline R3Box R3Mesh:: 
EdgeBBox(const R3MeshEdge *e) const
{
  // Returns the bounding box of the edge
  return R3Box(e->vertex[0]->position, e->vertex[1]->position);
}



inline int R3Mesh:: 
EdgeID(const R3MeshEdge *e) const
{
  // Returns the ID stored with a edge
  return e->id; 
}



inline RNFlags R3Mesh:: 
EdgeFlags(const R3MeshEdge *e) const
{
  // Returns the flags stored with a edge
  return e->flags; 
}



inline RNScalar R3Mesh:: 
EdgeValue(const R3MeshEdge *e) const
{
  // Returns the value stored with a edge
  return e->value; 
}



inline RNMark R3Mesh:: 
EdgeMark(const R3MeshEdge *e) const
{
  // Returns the mark stored with a edge
  return e->mark; 
}



inline void *R3Mesh:: 
EdgeData(const R3MeshEdge *e) const
{
  // Returns the data stored with a edge
  return e->data; 
}



////////////////////////////////////////////////////////////////////////
// FACE PROPERTY FUNCTIONS
////////////////////////////////////////////////////////////////////////

inline RNArea R3Mesh::
FaceArea(const R3MeshFace *f) const
{
  // Update the face area
  if (!(f->flags[R3_MESH_FACE_AREA_UPTODATE]))
    UpdateFaceArea((R3MeshFace *) f);

  // Return the area of the face
  return f->area;
}



inline const R3Plane& R3Mesh::
FacePlane(const R3MeshFace *f) const
{
  // Update the face plane
  if (!(f->flags[R3_MESH_FACE_PLANE_UPTODATE]))
    UpdateFacePlane((R3MeshFace *) f);

    // Return the plane containing the face
  return f->plane;
}



inline const R3Vector& R3Mesh:: 
FaceNormal(const R3MeshFace *f) const
{
  // Returns the normal of the face
  return FacePlane(f).Normal(); 
}



inline int R3Mesh:: 
FaceMaterial(const R3MeshFace *f) const
{
  // Returns the material associated with a face
  return f->material; 
}



inline int R3Mesh:: 
FaceSegment(const R3MeshFace *f) const
{
  // Returns the segment associated with a face
  return f->segment; 
}



inline int R3Mesh:: 
FaceCategory(const R3MeshFace *f) const
{
  // Returns the category associated with a face
  return f->category; 
}



inline int R3Mesh:: 
FaceID(const R3MeshFace *f) const
{
  // Returns the ID stored with a face
  return f->id; 
}



inline RNFlags R3Mesh:: 
FaceFlags(const R3MeshFace *f) const
{
  // Returns the flags stored with a face
  return f->flags; 
}



inline RNScalar R3Mesh:: 
FaceValue(const R3MeshFace *f) const
{
  // Returns the value stored with a face
  return f->value; 
}



inline RNMark R3Mesh:: 
FaceMark(const R3MeshFace *f) const
{
  // Returns the mark stored with a face
  return f->mark; 
}



inline void *R3Mesh:: 
FaceData(const R3MeshFace *f) const
{
  // Returns the data stored with a face
  return f->data; 
}



////////////////////////////////////////////////////////////////////////
// MANIPULATION FUNCTIONS
////////////////////////////////////////////////////////////////////////

inline void R3Mesh::
SetName(const char *name)
{
  // Copy name
  strncpy(this->name, name, R3_MESH_NAME_LENGTH);
}



////////////////////////////////////////////////////////////////////////
// VERTEX MANIPULATION FUNCTIONS
////////////////////////////////////////////////////////////////////////

inline void R3Mesh:: 
SetVertexNormal(R3MeshVertex *v, const R3Vector& normal)
{
  // Set new normal for vertex
  v->normal = normal;

  // Remember that vertex normal is up-to-date
  v->flags.Add(R3_MESH_VERTEX_NORMAL_UPTODATE);
}



inline void R3Mesh:: 
SetVertexColor(R3MeshVertex *v, const RNRgb& color)
{
  // Set new color for vertex
  v->color = color;
}



inline void R3Mesh:: 
SetVertexTextureCoords(R3MeshVertex *v, const R2Point& texcoords)
{
  // Set new texture coordinates for vertex
  v->texcoords = texcoords;
}



inline void R3Mesh:: 
SetVertexValue(R3MeshVertex *v, RNScalar value)
{
  // Set the value stored with a vertex
  v->value = value; 
}



inline void R3Mesh:: 
SetVertexMark(R3MeshVertex *v, RNMark mark)
{
  // Set the mark stored with a vertex
  v->mark = mark; 
}



inline void R3Mesh:: 
SetVertexData(R3MeshVertex *v, void *data)
{
  // Set the data stored with a vertex
  v->data = data; 
}



inline void R3Mesh:: 
SetVertexFlags(R3MeshVertex *v, RNFlags flags)
{
  // Set the flags stored with a vertex
  v->flags = flags; 
}



////////////////////////////////////////////////////////////////////////
// EDGE MANIPULATION FUNCTIONS
////////////////////////////////////////////////////////////////////////

inline void R3Mesh:: 
SetEdgeValue(R3MeshEdge *e, RNScalar value)
{
  // Set the value stored with a edge
  e->value = value; 
}



inline void R3Mesh:: 
SetEdgeMark(R3MeshEdge *e, RNMark mark)
{
  // Set the mark stored with a edge
  e->mark = mark; 
}



inline void R3Mesh:: 
SetEdgeData(R3MeshEdge *e, void *data)
{
  // Set the data stored with a edge
  e->data = data; 
}



inline void R3Mesh:: 
SetEdgeFlags(R3MeshEdge *e, RNFlags flags)
{
  // Set the flags stored with a edge
  e->flags = flags; 
}



////////////////////////////////////////////////////////////////////////
// FACE MANIPULATION FUNCTIONS
////////////////////////////////////////////////////////////////////////

inline void R3Mesh:: 
SetFacePlane(R3MeshFace *f, const R3Plane& plane)
{
  // Set the plane stored with a face
  f->plane = plane; 

  // Remember that face plane is up-to-date
  f->flags.Add(R3_MESH_FACE_PLANE_UPTODATE);
}



inline void R3Mesh:: 
SetFaceMaterial(R3MeshFace *f, int material)
{
  // Set the material stored with a face
  f->material = material; 
}



inline void R3Mesh:: 
SetFaceSegment(R3MeshFace *f, int segment)
{
  // Set the segment stored with a face
  f->segment = segment; 
}



inline void R3Mesh:: 
SetFaceCategory(R3MeshFace *f, int category)
{
  // Set the category stored with a face
  f->category = category; 
}



inline void R3Mesh:: 
SetFaceValue(R3MeshFace *f, RNScalar value)
{
  // Set the value stored with a face
  f->value = value; 
}



inline void R3Mesh:: 
SetFaceMark(R3MeshFace *f, RNMark mark)
{
  // Set the mark stored with a face
  f->mark = mark; 
}



inline void R3Mesh:: 
SetFaceData(R3MeshFace *f, void *data)
{
  // Set the data stored with a face
  f->data = data; 
}



inline void R3Mesh:: 
SetFaceFlags(R3MeshFace *f, RNFlags flags)
{
  // Set the flags stored with a face
  f->flags = flags; 
}



////////////////////////////////////////////////////////////////////////
// TOPOLOGY TRAVERSAL FUNCTIONS
////////////////////////////////////////////////////////////////////////

inline R3MeshEdge *R3Mesh:: 
EdgeOnVertex(const R3MeshVertex *v) const
{
  // Returns any edge connected to vertex
  if (v->edges.IsEmpty()) return NULL;
  else return v->edges[0];
}



inline R3MeshEdge *R3Mesh:: 
EdgeOnVertex(const R3MeshVertex *v, int k) const
{
  // Returns any face connected to vertex
  assert((0 <= k) && (k < edges.NEntries()));
  return v->edges[k];
}



inline R3MeshEdge *R3Mesh::
EdgeOnVertex(const R3MeshVertex *v, const R3MeshEdge *e, RNDirection dir) const
{
  // Go through the face on v that is in dir with respect to e
  R3MeshFace *f = FaceOnVertex(v, e, dir);
  return (f) ? EdgeOnVertex(v, f, dir) : NULL;
}



inline R3MeshVertex *R3Mesh:: 
VertexOnVertex(const R3MeshVertex *v) const
{
  // Returns any vertex adjacent to vertex
  R3MeshEdge *e = EdgeOnVertex(v);
  if (e) return VertexAcrossEdge(e, v);
  else return NULL;
}



inline R3MeshVertex *R3Mesh:: 
VertexOnVertex(const R3MeshVertex *v, int k) const
{
  // Returns kth vertex adjacent to vertex (adjacent vertices are not sorted!)
  R3MeshEdge *e = EdgeOnVertex(v, k);
  if (e) return VertexAcrossEdge(e, v);
  else return NULL;
}



inline R3MeshFace *R3Mesh:: 
FaceOnVertex(const R3MeshVertex *v) const
{
  // Returns any face connected to vertex
  for (int i = 0; i < VertexValence(v); i++) {
    R3MeshEdge *e = EdgeOnVertex(v, i);
    R3MeshFace *f = FaceOnEdge(e);
    if (f) return f;
  }

  // No face found
  return NULL;
}



inline R3MeshFace *R3Mesh:: 
FaceOnVertex(const R3MeshVertex *v, const R3MeshEdge *e, RNDirection dir) const
{
  // Returns face on vertex in dir with respect to edge
  return FaceOnEdge(e, v, dir); 
}



inline R3MeshFace *R3Mesh:: 
FaceOnVertex(const R3MeshVertex *v, const R3MeshFace *f, RNDirection dir) const
{
  // Returns face on vertex in dir with respect to face
  return FaceOnVertex(v, EdgeOnVertex(v, f, dir), dir); 
}



inline R3MeshVertex *R3Mesh:: 
VertexOnEdge(const R3MeshEdge *e) const
{
  // Returns any vertex on edge 
  return e->vertex[1]; 
}



inline R3MeshVertex *R3Mesh:: 
VertexOnEdge(const R3MeshEdge *e, int k) const
{
  // Returns kth vertex on edge 
  assert((k >= 0) && (k <= 1));
  return e->vertex[k]; 
}



inline R3MeshVertex *R3Mesh:: 
VertexAcrossEdge(const R3MeshEdge *e, const R3MeshVertex *v) const
{
  // Returns vertex across edge 
  assert((v == e->vertex[0]) || (v == e->vertex[1]));
  return (v == e->vertex[0]) ? e->vertex[1] : e->vertex[0];
}



inline R3MeshVertex *R3Mesh:: 
VertexOnEdge(const R3MeshEdge *e, const R3MeshFace *f, RNDirection dir) const
{
  // Returns vertex of edge in dir with respect to face
  assert((f == e->face[0]) || (f == e->face[1]));
  return (f == e->face[0]) ? e->vertex[1-dir] : e->vertex[dir];
}



inline R3MeshFace *R3Mesh:: 
FaceOnEdge(const R3MeshEdge *e) const
{
  // Returns any face on edge 
  return (e->face[1]) ? e->face[1] : e->face[0]; 
}



inline R3MeshFace *R3Mesh:: 
FaceOnEdge(const R3MeshEdge *e, int k) const
{
  // Returns kth face on edge 
  return e->face[k]; 
  
}



inline R3MeshFace *R3Mesh:: 
FaceOnEdge(const R3MeshEdge *e, const R3MeshVertex *v, RNDirection dir) const
{
  // Returns face on edge in dir with respect to vertex
  assert((v == e->vertex[0]) || (v == e->vertex[1]));
  return (v == e->vertex[0]) ? e->face[dir] : e->face[1-dir];
}



inline R3MeshFace *R3Mesh:: 
FaceAcrossEdge(const R3MeshEdge *e, const R3MeshFace *f) const
{
  // Returns face across edge 
  // assert((f == e->face[0]) || (f == e->face[1]));
  return (f == e->face[0]) ? e->face[1] : e->face[0];
}



inline R3MeshEdge *R3Mesh:: 
EdgeOnEdge(const R3MeshEdge *e) const
{
  // Returns any edge adjacent to e 
  R3MeshFace *f = FaceOnEdge(e);
  return (f) ? EdgeOnFace(f, e, RN_CCW) : NULL;
}



inline R3MeshEdge *R3Mesh:: 
EdgeOnEdge(const R3MeshEdge *e, const R3MeshVertex *v, RNDirection dir) const
{
  // Returns edge adjacent to e in dir with respect to vertex
  R3MeshFace *f = FaceOnEdge(e, v, dir);
  return (f) ? EdgeAcrossVertex(v, e, f) : NULL;
}



inline R3MeshVertex *R3Mesh:: 
VertexOnFace(const R3MeshFace *f) const
{
  // Returns any vertex on face
  return f->vertex[0];
}



inline R3MeshVertex *R3Mesh:: 
VertexOnFace(const R3MeshFace *f, int k) const
{
  // Returns kth vertex on face
  return f->vertex[k];
}



inline R3MeshVertex *R3Mesh:: 
VertexOnFace(const R3MeshFace *f, const R3MeshVertex *v, RNDirection dir) const
{
  // Returns vertex on face in dir with respect to vertex
  return VertexOnFace(f, EdgeOnFace(f, v, dir), dir); 
}



inline R3MeshVertex *R3Mesh:: 
VertexOnFace(const R3MeshFace *f, const R3MeshEdge *e, RNDirection dir) const
{
  // Returns vertex on face in dir with respect to edge
  return VertexOnEdge(e, f, dir); 
}



inline R3MeshVertex *R3Mesh:: 
VertexAcrossFace(const R3MeshFace *f, const R3MeshEdge *e) const
{
  // Returns the vertex across a triangle from an edge (the edge is the base and the vertex is the apex)
  assert(IsEdgeOnFace(e, f));
  if (!IsVertexOnEdge(f->vertex[0], e)) return f->vertex[0];
  if (!IsVertexOnEdge(f->vertex[1], e)) return f->vertex[1];
  if (!IsVertexOnEdge(f->vertex[2], e)) return f->vertex[2];
  return NULL;
}



inline R3MeshEdge *R3Mesh:: 
EdgeOnFace(const R3MeshFace *f) const
{
  // Returns any edge connected to face
  return f->edge[0]; 
}



inline R3MeshEdge *R3Mesh:: 
EdgeOnFace(const R3MeshFace *f, int k) const
{
  // Returns kth edge connected to face
  return f->edge[k]; 
}



inline R3MeshEdge *R3Mesh:: 
EdgeOnFace(const R3MeshFace *f, const R3MeshVertex *v, RNDirection dir) const
{
  // Returns edge on face in dir with respect to vertex
  if (f->vertex[0] == v) {
    if (dir == RN_CCW) return f->edge[0];
    else return f->edge[2];
  }
  else if (f->vertex[1] == v) {
    if (dir == RN_CCW) return f->edge[1];
    else return f->edge[0];
  }
  else if (f->vertex[2] == v) {
    if (dir == RN_CCW) return f->edge[2];
    else return f->edge[1];
  }
  else {
    // Vertex not found
    return NULL;
  }
}



inline R3MeshEdge *R3Mesh:: 
EdgeOnFace(const R3MeshFace *f, const R3MeshEdge *e, RNDirection dir) const
{
  // Return edge on face in dir with respect to edge
  if (f->edge[0] == e) {
    if (dir == RN_CCW) return f->edge[1];
    else return f->edge[2];
  }
  else if (f->edge[1] == e) {
    if (dir == RN_CCW) return f->edge[2];
    else return f->edge[0];
  }
  else if (f->edge[2] == e) {
    if (dir == RN_CCW) return f->edge[0];
    else return f->edge[1];
  }
  else {
    // Edge not found
    return NULL;
  }
}



inline R3MeshEdge *R3Mesh:: 
EdgeAcrossFace(const R3MeshFace *f, const R3MeshVertex *v) const
{
  // Returns the edge across a triangle from a vertex (the edge is the base and the vertex is the apex)
  assert(IsVertexOnFace(v, f));
  if (!IsVertexOnEdge(v, f->edge[0]))return f->edge[0];
  if (!IsVertexOnEdge(v, f->edge[1]))return f->edge[1];
  if (!IsVertexOnEdge(v, f->edge[2]))return f->edge[2];
  return NULL;
}



inline R3MeshFace *R3Mesh:: 
FaceOnFace(const R3MeshFace *f) const
{
  // Returns any face connected to face
  R3MeshFace *face;
  face = FaceAcrossEdge(f->edge[0], f);
  if (face) return face;
  face = FaceAcrossEdge(f->edge[1], f);
  if (face) return face;
  face = FaceAcrossEdge(f->edge[2], f);
  if (face) return face;
  return NULL;
}



inline R3MeshFace *R3Mesh:: 
FaceOnFace(const R3MeshFace *f, int k) const
{
  // Returns kth face connected to face
  return FaceAcrossEdge(f->edge[k], f);
}



////////////////////////////////////////////////////////////////////////
// TOPOLOGY QUERY FUNCTIONS
////////////////////////////////////////////////////////////////////////

inline RNBoolean R3Mesh:: 
IsVertexOnEdge(const R3MeshVertex *v, const R3MeshEdge *e) const
{
  // Returns whether vertex lies on edge
  return ((v == e->vertex[0]) || (v == e->vertex[1])); 
}



inline RNBoolean R3Mesh::
IsVertexOnFace(const R3MeshVertex *v, const R3MeshFace *f) const
{
  // Return whether vertex lies on face
  return (f->vertex[0] == v) || (f->vertex[1] == v) || (f->vertex[2] == v);
}



inline RNBoolean R3Mesh:: 
IsEdgeOnFace(const R3MeshEdge *e, const R3MeshFace *f) const
{
  // Returns whether edge lies on face
  return ((f == e->face[0]) || (f == e->face[1]));
} 



inline RNBoolean R3Mesh:: 
IsEdgeOnBoundary(const R3MeshEdge *e) const
{
  // Returns whether edge lies on boundary 
  return ((e->face[0] == NULL) || (e->face[1] == NULL));
} 



inline RNBoolean R3Mesh:: 
IsVertexOnMesh(const R3MeshVertex *v) const
{
  // Returns whether vertex is part of mesh
  if (!v) return FALSE;
  if ((v->id < 0) || (v->id >= NVertices())) return FALSE;
  if (vertices[v->id] != v) return FALSE;
  return TRUE;
} 



inline RNBoolean R3Mesh:: 
IsEdgeOnMesh(const R3MeshEdge *e) const
{
  // Returns whether edge is part of mesh
  if (!e) return FALSE;
  if ((e->id < 0) || (e->id >= NEdges())) return FALSE;
  if (edges[e->id] != e) return FALSE;
  return TRUE;
} 



inline RNBoolean R3Mesh:: 
IsFaceOnMesh(const R3MeshFace *f) const
{
  // Returns whether face is part of mesh
  if (!f) return FALSE;
  if ((f->id < 0) || (f->id >= NFaces())) return FALSE;
  if (faces[f->id] != f) return FALSE;
  return TRUE;
} 



////////////////////////////////////////////////////////////////////////
// DRAW FUNCTIONS
////////////////////////////////////////////////////////////////////////

inline void R3Mesh:: 
Draw(void) const
{
  // Draw the mesh
  DrawFaces();
} 


////////////////////////////////////////////////////////////////////////
// OTHER FUNCTIONS
////////////////////////////////////////////////////////////////////////

inline void R3Mesh:: 
SetData(void *data)
{
  // Set the user data stored with the mesh
  this->data = data; 
}



