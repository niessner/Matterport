/* Source file for the R3 mesh class */




// Include files

#include "R3Shapes/R3Shapes.h"
#include "ply.h"



// Public variables

RNMark R3mesh_mark = 1;



////////////////////////////////////////////////////////////////////////
// CONSTRUCTORS, DESTRUCTORS
////////////////////////////////////////////////////////////////////////

R3Mesh::
R3Mesh(void) 
  : vertex_block(NULL),
    edge_block(NULL),
    face_block(NULL),
    bbox(R3null_box),
    data(NULL)
{
  // Initialize name
  name[0] = '\0';

  // Just checking
  assert(IsValid());
}



R3Mesh::
R3Mesh(const R3Mesh& mesh)
  : vertex_block(NULL),
    edge_block(NULL),
    face_block(NULL),
    bbox(mesh.bbox),
    data(NULL)
{
  // Copy vertices 
  vertex_block = new R3MeshVertex [ mesh.NVertices() ];
  for (int i = 0; i < mesh.NVertices(); i++) {
    R3MeshVertex *vertex = mesh.Vertex(i);
    const R3Point& position = mesh.VertexPosition(vertex);
    const R3Vector& normal = mesh.VertexNormal(vertex);
    const RNRgb& color = mesh.VertexColor(vertex);
    const R2Point& texcoords = mesh.VertexTextureCoords(vertex);
    R3MeshVertex *copy_vertex = this->CreateVertex(position, normal, color, texcoords, &vertex_block[i] );
    if (this->VertexID(copy_vertex) != i) RNAbort("Mismatching vertex id"); 
  }

  // Copy edges
  edge_block = new R3MeshEdge [ mesh.NEdges() ];
  for (int i = 0; i < mesh.NEdges(); i++) {
    R3MeshEdge *edge = mesh.Edge(i);
    R3MeshVertex *v0 = mesh.VertexOnEdge(edge, 0);
    R3MeshVertex *v1 = mesh.VertexOnEdge(edge, 1);
    int i0 = mesh.VertexID(v0);
    int i1 = mesh.VertexID(v1);
    R3MeshVertex *copy_v0 = this->Vertex(i0);
    R3MeshVertex *copy_v1 = this->Vertex(i1);
    R3MeshEdge *copy_edge = this->CreateEdge(copy_v0, copy_v1, &edge_block[i] );
    if (this->EdgeID(copy_edge) != i) RNAbort("Mismatching edge id"); 
  }

  // Copy faces
  face_block = new R3MeshFace [ mesh.NFaces() ];
  for (int i = 0; i < mesh.NFaces(); i++) {
    R3MeshFace *face = mesh.Face(i);
    R3MeshVertex *v0 = mesh.VertexOnFace(face, 0);
    R3MeshVertex *v1 = mesh.VertexOnFace(face, 1);
    R3MeshVertex *v2 = mesh.VertexOnFace(face, 2);
    int i0 = mesh.VertexID(v0);
    int i1 = mesh.VertexID(v1);
    int i2 = mesh.VertexID(v2);
    R3MeshVertex *copy_v0 = this->Vertex(i0);
    R3MeshVertex *copy_v1 = this->Vertex(i1);
    R3MeshVertex *copy_v2 = this->Vertex(i2);
    R3MeshFace *copy_face = this->CreateFace(copy_v0, copy_v1, copy_v2, &face_block[i]);
    if (this->FaceID(copy_face) != i) RNAbort("Mismatching face id"); 
    this->SetFaceMaterial(copy_face, mesh.FaceMaterial(face));
    this->SetFaceSegment(copy_face, mesh.FaceSegment(face));
    this->SetFaceCategory(copy_face, mesh.FaceCategory(face));
  }
}



R3Mesh::
~R3Mesh(void) 
{
  // Delete everything
  Empty();
}



////////////////////////////////////////////////////////////////////////
// MESH PROPERTY FUNCTIONS
////////////////////////////////////////////////////////////////////////

R3Point R3Mesh::
Centroid(void) const
{
  // Return centroid of mesh
  RNArea area = 0;
  R3Point centroid(0,0,0);
  if (NFaces() > 0) {
    for (int i = 0; i < NFaces(); i++) {
      R3MeshFace *face = Face(i);
      RNArea face_area = FaceArea(face);
      R3Point face_centroid = FaceCentroid(face);
      centroid += face_area * face_centroid;
      area += face_area;
    }
  }
  else {
    for (int i = 0; i < NVertices(); i++) {
      R3MeshVertex *vertex = Vertex(i);
      centroid += VertexPosition(vertex);
      area += 1.0;
    }
  }

  // Return centroid
  if (area == 0) return centroid;
  else return centroid / area;
}



RNArea R3Mesh::
Area(void) const
{
  // Return area of mesh
  RNArea area = 0;
  for (int i = 0; i < NFaces(); i++) {
    R3MeshFace *face = Face(i);
    RNArea face_area = FaceArea(face);
    area += face_area;
  }

  // Return total area
  return area;
}



RNLength R3Mesh::
AverageEdgeLength(void) const
{
  // Check number of edges
  if (NEdges() == 0) return 0;

  // Sum edge lengths
  RNLength sum = 0;
  for (int i = 0; i < NEdges(); i++) {
    R3MeshEdge *edge = Edge(i);
    sum += EdgeLength(edge);
  }

  // Return average
  return sum / NEdges();
}



RNLength R3Mesh::
AverageRadius(const R3Point *center) const
{
  // Get centroid
  R3Point centroid = (center) ? *center : Centroid();

  // Check if mesh has faces
  if (NFaces() > 0) {
    // Compute average distance between a position on the surface and a center point
    RNArea area = 0;
    RNLength distance = 0;
    for (int i = 0; i < NFaces(); i++) {
      R3MeshFace *face = Face(i);
      RNArea face_area = FaceArea(face);
      R3Point face_centroid = FaceCentroid(face);
      distance += face_area * R3Distance(face_centroid, centroid);
      area += face_area;
    }

    // Return weighted average
    if (area == 0) return 0;
    else return distance / area;
  }
  else if (NVertices() > 0) {
    // Sum distances between vertices and center point
    RNLength sum = 0;
    for (int i = 0; i < NVertices(); i++) {
      R3MeshVertex *vertex = Vertex(i);
      const R3Point& position = VertexPosition(vertex);
      sum += R3Distance(position, centroid);
    }

    // Return average
    return sum / NVertices();
  }

  // Should not get here unless mesh is empty
  return 0.0;
}



R3Triad R3Mesh::
PrincipleAxes(const R3Point *center, RNScalar *variances) const
{
  // Get centroid
  R3Point centroid = (center) ? *center : Centroid();

  // Compute covariance matrix
  RNArea area = 0;
  RNScalar m[9] = { 0 };
  if (NFaces() > 0) {
    for (int i = 0; i < NFaces(); i++) {
      R3MeshFace *face = Face(i);
      RNArea face_area = FaceArea(face);
      R3Point face_centroid = FaceCentroid(face);
      RNScalar x = face_centroid[0] - centroid[0];
      RNScalar y = face_centroid[1] - centroid[1];
      RNScalar z = face_centroid[2] - centroid[2];
      m[0] += face_area * x*x;
      m[4] += face_area * y*y;
      m[8] += face_area * z*z;
      m[1] += face_area * x*y;
      m[3] += face_area * x*y;
      m[2] += face_area * x*z;
      m[6] += face_area * x*z;
      m[5] += face_area * y*z;
      m[7] += face_area * y*z;
      area += face_area;
    }
  }
  else {
    for (int i = 0; i < NVertices(); i++) {
      R3MeshVertex *vertex = Vertex(i);
      R3Point position = VertexPosition(vertex);
      RNScalar x = position[0] - centroid[0];
      RNScalar y = position[1] - centroid[1];
      RNScalar z = position[2] - centroid[2];
      m[0] += x*x;
      m[4] += y*y;
      m[8] += z*z;
      m[1] += x*y;
      m[3] += x*y;
      m[2] += x*z;
      m[6] += x*z;
      m[5] += y*z;
      m[7] += y*z;
      area += 1.0;
    }
  }

  // Normalize covariance matrix
  if (area == 0) return R3xyz_triad;
  for (int i = 0; i < 9; i++) m[i] /= area;

  // Compute eigenvalues and eigenvectors
  RNScalar U[9];
  RNScalar W[3];
  RNScalar Vt[9];
  RNSvdDecompose(3, 3, m, U, W, Vt);  // m == U . DiagonalMatrix(W) . Vt

  // Copy principle axes into more convenient form
  // W has eigenvalues (greatest to smallest) and Vt has eigenvectors (normalized)
  R3Vector axes[3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      axes[i][j] = Vt[3*i+j];
    }
  }

  // Flip axes so that "heavier" on positive side for first two dimensions
  int positive_count[3] = { 0, 0, 0 };
  int negative_count[3] = { 0, 0, 0 };
  for (int i = 0; i < NVertices(); i++) {
    R3MeshVertex *vertex = Vertex(i);
    R3Vector vertex_vector = VertexPosition(vertex) - centroid;
    for (int j = 0; j < 2; j++) {
      RNScalar dot = axes[j].Dot(vertex_vector);
      if (dot > 0.0) positive_count[j]++;
      else negative_count[j]++;
    }
  }
  for (int j =0; j < 2; j++) {
    if (positive_count[j] < negative_count[j]) {
      axes[j].Flip();
    }
  }

  // Set third axis to form orthonormal triad with other two
  axes[2] = axes[0] % axes[1];

  // Just checking
  assert(RNIsEqual(axes[0].Length(), 1.0, RN_BIG_EPSILON));
  assert(RNIsEqual(axes[1].Length(), 1.0, RN_BIG_EPSILON));
  assert(RNIsEqual(axes[2].Length(), 1.0, RN_BIG_EPSILON));
  assert(RNIsZero(axes[0].Dot(axes[1]), RN_BIG_EPSILON));
  assert(RNIsZero(axes[1].Dot(axes[2]), RN_BIG_EPSILON));
  assert(RNIsZero(axes[0].Dot(axes[2]), RN_BIG_EPSILON));

  // Return variances (eigenvalues)
  if (variances) {
    variances[0] = W[0];
    variances[1] = W[1];
    variances[2] = W[2];
  }

  // Return triad
  return R3Triad(axes[0], axes[1], axes[2]);
}



R3Affine R3Mesh::
PCANormalizationTransformation(RNBoolean translate, RNBoolean rotate, int scale) const
{
  // Initialize transformation
  R3Affine affine(R3identity_affine);

  // Compute center of mass
  R3Point centroid = Centroid();

  // Translate center of mass back to original (if not translating)
  if (!translate) {
    affine.Translate(centroid.Vector());
  }

  // Scale by inverse of radius
  if ((scale != 0) && (scale != 2)) {
    RNScalar radius = AverageRadius(&centroid);
    if (RNIsPositive(radius)) affine.Scale(1.0 / radius);
  }

  // Rotate to align principal axes with XYZ
  if (rotate || (scale == 2)) {
    RNScalar variances[3];
    R3Triad triad = PrincipleAxes(&centroid, variances);
    if (!rotate) affine.Transform(R3Affine(triad.InverseMatrix()));
    if (scale == 2) {
      if (variances[0] > 0) affine.XScale(1.0 / variances[0]);
      if (variances[1] > 0) affine.XScale(1.0 / variances[1]);
      if (variances[2] > 0) affine.XScale(1.0 / variances[2]);
    }
    affine.Transform(R3Affine(triad.InverseMatrix()));
  }

  // Translate center of mass to origin
  affine.Translate(-(centroid.Vector()));
  
  // Return PCA normalization transformation
  return affine;
}



////////////////////////////////////////////////////////////////////////
// VERTEX, EDGE, FACE PROPERTY FUNCTIONS
////////////////////////////////////////////////////////////////////////

RNArea R3Mesh::
VertexArea(const R3MeshVertex *v) const
{
  // Return 1/3 the sum of areas of attached faces
  RNLength sum = 0;
  for (int i = 0; i < VertexValence(v); i++) {
    R3MeshEdge *e = EdgeOnVertex(v, i);
    R3MeshFace *f = FaceOnEdge(e, v, RN_CCW);
    if (!f) continue;
    sum += FaceArea(f);
  }
  return sum / 3;
}



RNLength R3Mesh:: 
VertexAverageEdgeLength(const R3MeshVertex *v) const
{
  // Check number of edges
  if (VertexValence(v) == 0) return 0;

  // Returns the average length of edges attached to vertex
  RNLength sum = 0;
  for (int i = 0; i < VertexValence(v); i++) {
    R3MeshEdge *e = EdgeOnVertex(v, i);
    RNLength length = EdgeLength(e);
    sum += length;
  }
  return sum / VertexValence(v);
}



RNScalar R3Mesh::
VertexGaussCurvature(const R3MeshVertex *v) const
{
  // Sum areas and interior angles of adjacent faces
  RNArea area = 0;
  RNAngle angle = 0;
  const R3Point& p = VertexPosition(v);
  for (int i = 0; i < VertexValence(v); i++) {
    R3MeshEdge *e1 = EdgeOnVertex(v, i);
    if (!e1) { angle = 0; area = 0; break; }
    R3MeshVertex *v1 = VertexAcrossEdge(e1, v);
    R3MeshEdge *e2 = EdgeOnVertex(v, e1, RN_CCW);
    if (!e2) { angle = 0; area = 0; break; }
    R3MeshVertex *v2 = VertexAcrossEdge(e2, v);
    R3MeshFace *f = FaceOnVertex(v, e1, RN_CCW);
    if (!f) { angle = 0; area = 0; break; }
    const R3Point& p1 = VertexPosition(v1);
    const R3Point& p2 = VertexPosition(v2);
    angle += R3InteriorAngle(p1 - p, p2 - p);
    area += FaceArea(f);
  }

  // Compute curvature using Gauss-Bonet
  RNScalar curvature = (area > 0) ? 3  * (RN_TWO_PI - angle) / area : 0;

  // Return Gauss curvature
  return curvature;
}



RNScalar R3Mesh::
VertexMeanCurvature(const R3MeshVertex *v) const
{
  // Compute signed length of laplacian vector projected onto normal
  RNArea area = VertexArea(v);
  if (area == 0) return 0;
  R3Vector laplacian = VertexLaplacianVector(v);
  RNScalar curvature = -1 * laplacian.Dot(VertexNormal(v));
  return curvature;
}



R3Vector R3Mesh::
VertexLaplacianVector(const R3MeshVertex *v1) const
{
  // Compute vector from vertex to cotan weighed average of neighbors 
  // From http://www.mpi-inf.mpg.de/~ag4-gm/handouts/06gm_surf3.pdf
  R3Vector laplacian = R3zero_vector;
  const R3Point& p1 = VertexPosition(v1);
  for (int j = 0; j < VertexValence(v1); j++) {
    R3MeshEdge *e = EdgeOnVertex(v1, j);
    R3MeshVertex *v2 = VertexAcrossEdge(e, v1);
    const R3Point& p2 = VertexPosition(v2);

    // Compute cotan weight
    double weight = 0;
    for (int k = 0; k < 2; k++) {
      R3MeshFace *f = FaceOnEdge(e, k);
      if (!f) continue;
      R3MeshVertex *v3 = VertexAcrossFace(f, e);
      const R3Point& p3 = VertexPosition(v3);
      R3Vector vec1 = p1 - p3; vec1.Normalize();
      R3Vector vec2 = p2 - p3; vec2.Normalize();
      RNAngle angle = R3InteriorAngle(vec1, vec2);
      if (angle == 0) continue; 
      double tan_angle = tan(angle);
      if (tan_angle == 0) continue; 
      weight += 1.0 / tan_angle;
    }

    // Add weighted position
    laplacian += weight * (p2 - p1);
  }

  // Normalize by area (maybe should be 3X -- i.e., for face areas, not vertex area)
  laplacian /= 4 * VertexArea(v1);

  // Return laplacian vector
  return laplacian;
}



R3Vector R3Mesh:: 
EdgeDirection(const R3MeshEdge *e) const
{
  // Returns the normalized vector pointing in direction of the edge
  R3Vector vector = EdgeVector(e);
  RNLength d = EdgeLength(e);
  if (RNIsPositive(d)) vector /= d;
  return vector;
}



RNAngle R3Mesh:: 
EdgeInteriorAngle(const R3MeshEdge *e) const
{
  // Get attached faces
  R3MeshFace *f0 = FaceOnEdge(e, 0);
  R3MeshFace *f1 = FaceOnEdge(e, 1);
  if (!f0 || !f1) return 0.0;

  // Returns the dihedral angle between the attached faces (in range 0 to 2*PI)
  R3Vector cross_normals = FaceNormal(f0) % FaceNormal(f1);
  RNScalar cos_angle = FaceNormal(f0).Dot(FaceNormal(f1));
  RNScalar acos_angle = (RNIsEqual(cos_angle * cos_angle, 1.0, 1.0E-6)) ? 0.0 : acos(cos_angle);
  RNBoolean convex = (cross_normals.Dot(EdgeVector(e)) > 0);
  RNScalar angle = (convex) ? (RN_PI - acos_angle) : (RN_PI + acos_angle);
  assert((angle >= 0.0) && (angle <= RN_TWO_PI));
  return angle;
}



R3Vector R3Mesh::
EdgeNormal(const R3MeshEdge *e) const
{
  // Return the normal to the surface at the edge
  R3Vector normal = R3zero_vector;
  R3MeshFace *f0 = FaceOnEdge(e, 0);
  R3MeshFace *f1 = FaceOnEdge(e, 1);
  if (f0) normal += FaceNormal(f0);
  if (f1) normal += FaceNormal(f1);
  normal.Normalize();
  return normal;
}



RNScalar R3Mesh:: 
EdgeAspect(const R3MeshEdge *e) const
{
  // Get attached faces
  R3MeshFace *f0 = FaceOnEdge(e, 0);
  R3MeshFace *f1 = FaceOnEdge(e, 1);
  if (!f0 || !f1) return 0.0;

  // Get edge length
  RNScalar edge_length = EdgeLength(e);
  if (RNIsZero(edge_length)) return RN_INFINITY;

  // Get distances between opposite vertices
  R3MeshVertex *v0 = VertexAcrossFace(f0, e);
  R3MeshVertex *v1 = VertexAcrossFace(f1, e);
  RNScalar opposite_length = R3Distance(v0->position, v1->position);

  // Return ratio
  return opposite_length / edge_length;
}



RNScalar R3Mesh:: 
FaceAspect(const R3MeshFace *f) const
{
  R3MeshEdge *e0 = EdgeOnFace(f, 0);
  R3MeshEdge *e1 = EdgeOnFace(f, 1);
  R3MeshEdge *e2 = EdgeOnFace(f, 2);
  RNLength a = EdgeLength(e0);
  if (RNIsZero(a)) return RN_INFINITY;
  RNLength b = EdgeLength(e1);
  if (RNIsZero(b)) return RN_INFINITY;
  RNLength c = EdgeLength(e2);
  if (RNIsZero(c)) return RN_INFINITY;
  RNScalar s = (a + b + c) / 2.0;
  return (a * b * c) / (8.0 * (s - a) * (s - b) * (s - c));
}



R3Point R3Mesh::
FaceCentroid(const R3MeshFace *f) const
{
  // Return face centroid
  R3Point centroid = R3zero_point;
  centroid += f->vertex[0]->position;
  centroid += f->vertex[1]->position;
  centroid += f->vertex[2]->position;
  return centroid / 3.0;
}



R3Point R3Mesh::
FacePoint(const R3MeshFace *f, RNMagnitude barycentrics[3]) const
{
  // Return point on the face with the given barycentric coordinates
  R3Point point = R3zero_point;
  point += barycentrics[0] * f->vertex[0]->position;
  point += barycentrics[1] * f->vertex[1]->position;
  point += barycentrics[2] * f->vertex[2]->position;
  return point;
}



R3Point R3Mesh::
FaceBarycentric(const R3MeshFace *f, const R3Point &point) const
{
  // From http://www.devmaster.net/wiki/Ray-triangle_intersection
  R3Point p0 = f->vertex[0]->position;
  R3Point p1 = f->vertex[1]->position;
  R3Point p2 = f->vertex[2]->position;
  R3Vector b = p1 - p0;
  R3Vector c = p2 - p0;
  R3Vector p = point - p0;
  RNDimension dim = f->plane.Normal().MaxDimension();
  RNDimension dim0 = (dim+1)%3;
  RNDimension dim1 = (dim+2)%3;
  RNScalar denom = b[dim1]*c[dim0] - b[dim0]*c[dim1];
  if (denom == 0) return R3zero_point;
  RNScalar b1 = (p[dim1]*c[dim0] - p[dim0]*c[dim1]) / denom;
  RNScalar b2 = (p[dim1]*b[dim0] - p[dim0]*b[dim1]) / -denom;
  RNScalar b0 = 1 - b1 - b2;
  return R3Point(b0, b1, b2);
}



const R3Box& R3Mesh::
FaceBBox(const R3MeshFace *f) const
{
  // Update the face bbox
  if (!(f->flags[R3_MESH_FACE_BBOX_UPTODATE]))
    UpdateFaceBBox((R3MeshFace *) f);

  // Return the bbox of the face
  return f->bbox;
}



////////////////////////////////////////////////////////////////////////
// TOPOLOGY QUERY FUNCTIONS
////////////////////////////////////////////////////////////////////////

RNBoolean R3Mesh::
IsVertexOnBoundary(const R3MeshVertex *v) const
{
  // Return whether vertex lies on boundary
  for (int i = 0; i < VertexValence(v); i++) {
    R3MeshEdge *e = EdgeOnVertex(v, i);
    if (IsEdgeOnBoundary(e)) return TRUE;
  }

  // No adjacent edge is on boundary
  return FALSE;
}



RNBoolean R3Mesh::
IsFaceOnBoundary(const R3MeshFace *f) const
{
  // Return whether face has any edge on boundary
  if (IsEdgeOnBoundary(EdgeOnFace(f, 0))) return TRUE;
  if (IsEdgeOnBoundary(EdgeOnFace(f, 1))) return TRUE;
  if (IsEdgeOnBoundary(EdgeOnFace(f, 2))) return TRUE;
  return FALSE;
}



////////////////////////////////////////////////////////////////////////
// GEOMETRY MANIPULATION FUNCTIONS
////////////////////////////////////////////////////////////////////////

void R3Mesh::
SetVertexPosition(R3MeshVertex *v, const R3Point& position)
{
  // Set vertex position
  v->position = position;

  // Mark vertex in need of update to normal and curvature
  v->flags.Remove(R3_MESH_VERTEX_NORMAL_UPTODATE);
  v->flags.Remove(R3_MESH_VERTEX_CURVATURE_UPTODATE);

  // Mark edges/faces in need of update
  for (int i = 0; i < v->edges.NEntries(); i++) {
    R3MeshEdge *e = v->edges[i];
    e->flags.Remove(R3_MESH_EDGE_LENGTH_UPTODATE);
    R3MeshVertex *neighbor = VertexAcrossEdge(e, v);
    neighbor->flags.Remove(R3_MESH_VERTEX_NORMAL_UPTODATE);
    neighbor->flags.Remove(R3_MESH_VERTEX_CURVATURE_UPTODATE);
    R3MeshFace *f0 = e->face[0];
    if (f0) f0->flags.Remove(R3_MESH_FACE_PLANE_UPTODATE | R3_MESH_FACE_BBOX_UPTODATE);
    R3MeshFace *f1 = e->face[1];
    if (f1) f1->flags.Remove(R3_MESH_FACE_PLANE_UPTODATE | R3_MESH_FACE_BBOX_UPTODATE);
  }

  // Update bounding box
  bbox.Union(position);
}



////////////////////////////////////////////////////////////////////////
// TOPOLOGY TRAVERSAL FUNCTIONS
////////////////////////////////////////////////////////////////////////

R3MeshEdge *R3Mesh::
EdgeBetweenVertices(const R3MeshVertex *v1, const R3MeshVertex *v2) const
{
  // Check if vertices are same
  if (v1 == v2) return NULL;

  // Search for edge in v1's list
  for (int i = 0; i < v1->edges.NEntries(); i++) {
    R3MeshEdge *e = v1->edges[i];
    if (IsVertexOnEdge(v2, e)) return e;
  }

  // Not found
  return NULL;
}



R3MeshEdge *R3Mesh::
EdgeOnVertex(const R3MeshVertex *v, const R3MeshFace *f, RNDirection dir) const
{
  // Return edge on v in dir with respect to f
  if (f->vertex[0] == v) {
    if (dir == RN_CCW) return f->edge[2];
    else return f->edge[0];
  }
  else if (f->vertex[1] == v) {
    if (dir == RN_CCW) return f->edge[0];
    else return f->edge[1];
  }
  else if (f->vertex[2] == v) {
    if (dir == RN_CCW) return f->edge[1];
    else return f->edge[2];
  }
  else {
    // Vertex is not on face
    return NULL;
  }
}



R3MeshEdge *R3Mesh:: 
EdgeAcrossVertex(const R3MeshVertex *v, const R3MeshEdge *e, const R3MeshFace *f) const
{
  // Returns edge on the other side of a vertex from an edge on the same face 
  if (v == f->vertex[0]) {
    if (e == f->edge[0]) return f->edge[2];
    else if (e == f->edge[2]) return f->edge[0];
  }
  else if (v == f->vertex[1]) {
    if (e == f->edge[1]) return f->edge[0];
    else if (e == f->edge[0]) return f->edge[1];
  }
  else if (v == f->vertex[2]) {
    if (e == f->edge[2]) return f->edge[1];
    else if (e == f->edge[1]) return f->edge[2];
  }
  return NULL;
}



R3MeshFace *R3Mesh::
FaceBetweenVertices(const R3MeshVertex *v1, const R3MeshVertex *v2, RNDirection dir) const
{
  RNAbort("Not implemented");
  return NULL;
}



R3MeshVertex *R3Mesh::
VertexBetweenEdges(const R3MeshEdge *e1, const R3MeshEdge *e2) const
{
  // Return vertex shared by edges
  if (e1 == e2) return NULL;
  if ((e1->vertex[0] == e2->vertex[0]) || (e1->vertex[0] == e2->vertex[1])) return e1->vertex[0];
  else if ((e1->vertex[1] == e2->vertex[0]) || (e1->vertex[1] == e2->vertex[1])) return e1->vertex[1];
  else return NULL;
}



R3MeshFace *R3Mesh::
FaceBetweenEdges(const R3MeshEdge *e1, const R3MeshEdge *e2) const
{
  // Return face shared by edges
  if (e1 == e2) return NULL;
  if ((e1->face[0]) && ((e1->face[0] == e2->face[0]) || (e1->face[0] == e2->face[1]))) return e1->face[0];
  else if ((e1->face[1]) && ((e1->face[1] == e2->face[0]) || (e1->face[1] == e2->face[1]))) return e1->face[1];
  else return NULL;
}



R3MeshVertex *R3Mesh::
VertexBetweenFaces(const R3MeshFace *f1, const R3MeshFace *f2, RNDirection dir) const
{
  // Check if faces are same
  if (f1 == f2) return NULL;

  // Return vertex connected to both f1 and f2
  R3MeshEdge *e = EdgeBetweenFaces(f1, f2);
  if (e) {
    // Faces lie on same edge, return vertex in dir with respect to f1
    return VertexOnFace(f1, e, dir);
  }
  else {
    // Faces do not lie on the same edge, check if they share a vertex
    if (f1->vertex[0] == f2->vertex[0]) return f1->vertex[0];
    if (f1->vertex[0] == f2->vertex[1]) return f1->vertex[0];
    if (f1->vertex[0] == f2->vertex[2]) return f1->vertex[0];
    if (f1->vertex[1] == f2->vertex[0]) return f1->vertex[1];
    if (f1->vertex[1] == f2->vertex[1]) return f1->vertex[1];
    if (f1->vertex[1] == f2->vertex[2]) return f1->vertex[1];
    if (f1->vertex[2] == f2->vertex[0]) return f1->vertex[2];
    if (f1->vertex[2] == f2->vertex[1]) return f1->vertex[2];
    if (f1->vertex[2] == f2->vertex[2]) return f1->vertex[2];
    return NULL;
  }
}



R3MeshEdge *R3Mesh::
EdgeBetweenFaces(const R3MeshFace *f1, const R3MeshFace *f2) const
{
  if (f1 == f2) return NULL;
  if (FaceAcrossEdge(f1->edge[0], f1) == f2) return f1->edge[0];
  if (FaceAcrossEdge(f1->edge[1], f1) == f2) return f1->edge[1];
  if (FaceAcrossEdge(f1->edge[2], f1) == f2) return f1->edge[2];
  return NULL;
}



////////////////////////////////////////////////////////////////////////
// TOPOLOGY QUERY FUNCTIONS
////////////////////////////////////////////////////////////////////////

void R3Mesh::
FindConnectedVertices(R3MeshVertex *seed, RNArray<R3MeshVertex *>& vertices)
{
  // Fill array with vertices in same connected component as seed
  R3mesh_mark++;
  RNArray<R3MeshVertex *> stack;
  stack.Insert(seed);
  SetVertexMark(seed, R3mesh_mark);
  while (!stack.IsEmpty()) {
    R3MeshVertex *vertex = stack.Tail();
    stack.RemoveTail();
    vertices.Insert(vertex);
    for (int i = 0; i < VertexValence(vertex); i++) {
      R3MeshEdge *edge = EdgeOnVertex(vertex, i);
      R3MeshVertex *neighbor = VertexAcrossEdge(edge, vertex);
      if (VertexMark(neighbor) != R3mesh_mark) {
        SetVertexMark(neighbor, R3mesh_mark);
        stack.Insert(neighbor);
      }
    }
  }
}


void R3Mesh::
FindConnectedFaces(R3MeshFace *seed, RNArray<R3MeshFace *>& faces)
{
  // Fill array with faces in same connected component as seed
  R3mesh_mark++;
  RNArray<R3MeshFace *> stack;
  stack.Insert(seed);
  while (!stack.IsEmpty()) {
    R3MeshFace *face = stack.Tail();
    stack.RemoveTail();
    if (FaceMark(face) != R3mesh_mark) {
      faces.Insert(face);
      if (FaceOnFace(face, 0)) stack.Insert(FaceOnFace(face, 0));
      if (FaceOnFace(face, 1)) stack.Insert(FaceOnFace(face, 1));
      if (FaceOnFace(face, 2)) stack.Insert(FaceOnFace(face, 2));
      SetFaceMark(face, R3mesh_mark);
    }
  }
}


////////////////////////////////////////////////////////////////////////
// MESH MANIPULATION FUNCTIONS
////////////////////////////////////////////////////////////////////////

void R3Mesh::
Empty(void)
{
  // Delete all faces, edges, vertices
  while (NFaces() > 0) DeleteFace(Face(0));
  while (NEdges() > 0) DeleteEdge(Edge(0));
  while (NVertices() > 0) DeleteVertex(Vertex(0));

  // Delete the blocks of data
  if (vertex_block) { delete [] vertex_block; vertex_block = NULL; }
  if (edge_block) { delete [] edge_block; edge_block = NULL; }
  if (face_block) { delete [] face_block; face_block = NULL; }

  // Reset bounding box
  bbox = R3null_box;
}



void R3Mesh::
Smooth(RNScalar factor)
{
  // Copy vertex positions
  R3Point *positions = new R3Point [ NVertices() ];
  for (int i = 0; i < NVertices(); i++) {
    R3MeshVertex *vertex = Vertex(i);
    positions[i] = VertexPosition(vertex);
  }

  // Update vertex positions
  for (int i = 0; i < NVertices(); i++) {
    R3MeshVertex *vertex = Vertex(i);
    RNScalar weight = 1;
    R3Point position = positions[vertex->id];
    RNLength radius = VertexAverageEdgeLength(vertex);
    if (radius > 0) {
      for (int j = 0; j < VertexValence(vertex); j++) {
        R3MeshEdge *edge = EdgeOnVertex(vertex, j);
        R3MeshVertex *neighbor_vertex = VertexAcrossEdge(edge, vertex);
        const R3Point& neighbor_position = positions[neighbor_vertex->id];
        RNLength length = EdgeLength(edge);
        RNLength sigma = (length > 0) ? length / radius : 1;
        RNScalar w = factor * exp((length * length) / (-2.0 * sigma * sigma));
        position += w * neighbor_position;
        weight += w;
      }
      
      // Update vertex position
      R3Point smooth_position = position / weight;
      SetVertexPosition(vertex, smooth_position);
    }
  }

  // Delete copy of vertex positions
  delete [] positions;
}



void R3Mesh::
Sharpen(RNScalar factor)
{
  // Copy vertex positions
  R3Point *positions = new R3Point [ NVertices() ];
  for (int i = 0; i < NVertices(); i++) {
    R3MeshVertex *vertex = Vertex(i);
    positions[i] = VertexPosition(vertex);
  }

  // Update vertex positions
  for (int i = 0; i < NVertices(); i++) {
    R3MeshVertex *vertex = Vertex(i);
    RNScalar weight = 1;
    R3Point position = positions[vertex->id];
    RNLength radius = VertexAverageEdgeLength(vertex);
    for (int j = 0; j < VertexValence(vertex); j++) {
      R3MeshEdge *edge = EdgeOnVertex(vertex, j);
      R3MeshVertex *neighbor_vertex = VertexAcrossEdge(edge, vertex);
      const R3Point& neighbor_position = positions[neighbor_vertex->id];
      RNLength length = EdgeLength(edge);
      RNLength sigma = length / radius;
      RNScalar w = factor * exp((length * length) / (-2.0 * sigma * sigma));
      position += w * neighbor_position;
      weight += w;
    }

    // Update vertex position
    R3Point smooth_position = position / weight;
    R3Point sharpen_position = 2 * positions[vertex->id] - smooth_position.Vector();
    SetVertexPosition(vertex, sharpen_position);
  }

  // Delete copy of vertex positions
  delete [] positions;
}



void R3Mesh::
AddRandomNoise(RNScalar factor)
{
  // Change every vertex position
  bbox = R3null_box;
  for (int i = 0; i < NVertices(); i++) {
    R3MeshVertex *vertex = Vertex(i);
    RNLength max_distance = factor * VertexAverageEdgeLength(vertex);
    R3Vector random_vector(RNRandomScalar(), RNRandomScalar(), RNRandomScalar());
    vertex->position += max_distance * random_vector;
    vertex->flags.Remove(R3_MESH_VERTEX_NORMAL_UPTODATE);
    bbox.Union(vertex->position);
  }

  // Mark every edge out of date
  for (int i = 0; i < NEdges(); i++) {
    R3MeshEdge *edge = Edge(i);
    edge->flags.Remove(R3_MESH_EDGE_LENGTH_UPTODATE);
  }

  // Mark every face out of date
  for (int i = 0; i < NFaces(); i++) {
    R3MeshFace *face = Face(i);
    face->flags.Remove(R3_MESH_FACE_PLANE_UPTODATE | R3_MESH_FACE_BBOX_UPTODATE);
  }
}



void R3Mesh::
Inflate(RNScalar factor)
{
  // Change every vertex position
  bbox = R3null_box;
  for (int i = 0; i < NVertices(); i++) {
    R3MeshVertex *vertex = Vertex(i);
    RNLength distance = factor * VertexAverageEdgeLength(vertex);
    R3Vector normal_vector = VertexNormal(vertex);
    vertex->position += distance * normal_vector;
    vertex->flags.Remove(R3_MESH_VERTEX_NORMAL_UPTODATE);
    bbox.Union(vertex->position);
  }

  // Mark every edge out of date
  for (int i = 0; i < NEdges(); i++) {
    R3MeshEdge *edge = Edge(i);
    edge->flags.Remove(R3_MESH_EDGE_LENGTH_UPTODATE);
  }

  // Mark every face out of date
  for (int i = 0; i < NFaces(); i++) {
    R3MeshFace *face = Face(i);
    face->flags.Remove(R3_MESH_FACE_PLANE_UPTODATE | R3_MESH_FACE_BBOX_UPTODATE);
  }
}



void R3Mesh::
Transform(const R3Transformation& transformation)
{
  // Transform every vertex position
  bbox = R3null_box;
  for (int i = 0; i < NVertices(); i++) {
    R3MeshVertex *vertex = Vertex(i);
    vertex->position.Transform(transformation);
    vertex->flags.Remove(R3_MESH_VERTEX_NORMAL_UPTODATE);
    bbox.Union(vertex->position);
  }

  // Mark every edge out of date
  for (int i = 0; i < NEdges(); i++) {
    R3MeshEdge *edge = Edge(i);
    edge->flags.Remove(R3_MESH_EDGE_LENGTH_UPTODATE);
  }

  // Mark every face out of date
  for (int i = 0; i < NFaces(); i++) {
    R3MeshFace *face = Face(i);
    face->flags.Remove(R3_MESH_FACE_PLANE_UPTODATE | R3_MESH_FACE_BBOX_UPTODATE);
  }
}




////////////////////////////////////////////////////////////////////////
// TOPOLOGY MANIPULATION FUNCTIONS
////////////////////////////////////////////////////////////////////////

R3MeshVertex *R3Mesh::
CreateVertex(const R3MeshVertex& source_vertex, R3MeshVertex *v)
{
  // Create vertex
  if (!v) {
    v = new R3MeshVertex();
    v->flags.Add(R3_MESH_VERTEX_ALLOCATED);
  }

  // Set position/normal of new vertex
  SetVertexPosition(v, source_vertex.position);
  if (!source_vertex.normal.IsZero()) SetVertexNormal(v, source_vertex.normal);
  SetVertexColor(v, source_vertex.color);
  SetVertexTextureCoords(v, source_vertex.texcoords);

  // Set ID of new vertex
  v->id = vertices.NEntries();

  // Insert vertex into array
  vertices.Insert(v);

  // Return vertex
  return v;
}



R3MeshVertex *R3Mesh::
CreateVertex(const R3Point& position, R3MeshVertex *v)
{
  // Create vertex
  if (!v) {
    v = new R3MeshVertex();
    v->flags.Add(R3_MESH_VERTEX_ALLOCATED);
  }

  // Set position of new vertex
  SetVertexPosition(v, position);

  // Set ID of new vertex
  v->id = vertices.NEntries();

  // Insert vertex into array
  vertices.Insert(v);

  // Return vertex
  return v;
}



R3MeshVertex *R3Mesh::
CreateVertex(const R3Point& position, const R3Vector& normal, R3MeshVertex *v)
{
  // Create vertex
  if (!v) {
    v = new R3MeshVertex();
    v->flags.Add(R3_MESH_VERTEX_ALLOCATED);
  }

  // Set position/normal of new vertex
  SetVertexPosition(v, position);
  if (!normal.IsZero()) SetVertexNormal(v, normal);

  // Set ID of new vertex
  v->id = vertices.NEntries();

  // Insert vertex into array
  vertices.Insert(v);

  // Return vertex
  return v;
}



R3MeshVertex *R3Mesh::
CreateVertex(const R3Point& position, const R3Vector& normal, const RNRgb& color, R3MeshVertex *v)
{
  // Create vertex
  if (!v) {
    v = new R3MeshVertex();
    v->flags.Add(R3_MESH_VERTEX_ALLOCATED);
  }

  // Set position/normal of new vertex
  SetVertexPosition(v, position);
  if (!normal.IsZero()) SetVertexNormal(v, normal);
  SetVertexColor(v, color);

  // Set ID of new vertex
  v->id = vertices.NEntries();

  // Insert vertex into array
  vertices.Insert(v);

  // Return vertex
  return v;
}



R3MeshVertex *R3Mesh::
CreateVertex(const R3Point& position, const R3Vector& normal, const RNRgb& color, const R2Point& texcoords, R3MeshVertex *v)
{
  // Create vertex
  if (!v) {
    v = new R3MeshVertex();
    v->flags.Add(R3_MESH_VERTEX_ALLOCATED);
  }

  // Set position/normal of new vertex
  SetVertexPosition(v, position);
  if (!normal.IsZero()) SetVertexNormal(v, normal);
  SetVertexColor(v, color);
  SetVertexTextureCoords(v, texcoords);

  // Set ID of new vertex
  v->id = vertices.NEntries();

  // Insert vertex into array
  vertices.Insert(v);

  // Return vertex
  return v;
}



R3MeshEdge *R3Mesh::
CreateEdge(R3MeshVertex *v1, R3MeshVertex *v2, R3MeshEdge *e)
{
  // Create edge
  if (!e) {
    e = new R3MeshEdge();
    e->flags.Add(R3_MESH_EDGE_ALLOCATED);
  }

  // Update edge-vertex relations
  e->vertex[0] = v1;
  e->vertex[1] = v2;

  // Set ID of new edge
  e->id = edges.NEntries();

  // Insert edge into vertex lists
  v1->edges.Insert(e);
  v2->edges.Insert(e);

  // Insert edge into array
  edges.Insert(e);

  // Return edge
  return e;
}



R3MeshFace *R3Mesh::
CreateFace(R3MeshVertex *v1, R3MeshVertex *v2, R3MeshVertex *v3, R3MeshFace *f)
{
  // Get/create edges
  R3MeshEdge *e1 = EdgeBetweenVertices(v1, v2);
  if (!e1) e1 = CreateEdge(v1, v2);
  R3MeshEdge *e2 = EdgeBetweenVertices(v2, v3);
  if (!e2) e2 = CreateEdge(v2, v3);
  R3MeshEdge *e3 = EdgeBetweenVertices(v3, v1);
  if (!e3) e3 = CreateEdge(v3, v1);
  assert(e1 && e2 && e3);

  // Create face
  return CreateFace(v1, v2, v3, e1, e2, e3, f);
}



R3MeshFace *R3Mesh::
CreateFace(R3MeshEdge *e1, R3MeshEdge *e2, R3MeshEdge *e3, R3MeshFace *f)
{
  // Get vertices
  R3MeshVertex *v1 = VertexBetweenEdges(e3, e1);
  R3MeshVertex *v2 = VertexBetweenEdges(e1, e2);
  R3MeshVertex *v3 = VertexBetweenEdges(e2, e3);
  assert(v1 && v2 && v3);

  // Create face
  return CreateFace(v1, v2, v3, e1, e2, e3, f);
}



R3MeshFace *R3Mesh::
CreateFace(R3MeshVertex *v1, R3MeshVertex *v2, R3MeshVertex *v3,
           R3MeshEdge *e1, R3MeshEdge *e2, R3MeshEdge *e3, R3MeshFace *f)
{
  // Check if two faces share same side of same edge
  if ((e1->vertex[0] == v1) && e1->face[0]) return NULL;
  if ((e1->vertex[0] == v2) && e1->face[1]) return NULL;
  if ((e2->vertex[0] == v2) && e2->face[0]) return NULL;
  if ((e2->vertex[0] == v3) && e2->face[1]) return NULL;
  if ((e3->vertex[0] == v3) && e3->face[0]) return NULL;
  if ((e3->vertex[0] == v1) && e3->face[1]) return NULL;

  // Create face
  if (!f) {
    f = new R3MeshFace();
    f->flags.Add(R3_MESH_FACE_ALLOCATED);
  }

  // Update face pointers
  UpdateFaceRefs(f, v1, v2, v3, e1, e2, e3);

  // Set face ID
  f->id = faces.NEntries();

  // Insert face into array
  faces.Insert(f);

  // Return face
  return f;
}



void R3Mesh::
DeallocateVertex(R3MeshVertex *v)
{
  // Remove vertex from mesh array by moving last vertex
  assert(!vertices.IsEmpty());
  R3MeshVertex *tail = vertices.Tail();
  int index = v->id;
  RNArrayEntry *entry = vertices.KthEntry(index);
  assert(vertices.EntryContents(entry) == v);
  vertices.EntryContents(entry) = tail;
  tail->id = index;
  vertices.RemoveTail();

  // Reset ID to ease debugging
  v->id = -1;

  // Deallocate vertex
  if (v->flags[R3_MESH_VERTEX_ALLOCATED]) delete v;
}



void R3Mesh::
DeallocateEdge(R3MeshEdge *e)
{
  // Remove edge from mesh array by moving last edge
  assert(!edges.IsEmpty());
  R3MeshEdge *tail = edges.Tail();
  int index = e->id;
  RNArrayEntry *entry = edges.KthEntry(index);
  assert(edges.EntryContents(entry) == e);
  edges.EntryContents(entry) = tail;
  tail->id = index;
  edges.RemoveTail();

  // Reset ID to ease debugging
  e->id = -1;

  // Deallocate edge
  if (e->flags[R3_MESH_EDGE_ALLOCATED]) delete e;
}



void R3Mesh::
DeallocateFace(R3MeshFace *f)
{
  // Remove face from mesh array by moving last face
  assert(!faces.IsEmpty());
  R3MeshFace *tail = faces.Tail();
  int index = f->id;
  RNArrayEntry *entry = faces.KthEntry(index);
  assert(faces.EntryContents(entry) == f);
  faces.EntryContents(entry) = tail;
  tail->id = index;
  faces.RemoveTail();

  // Reset ID to ease debugging
  f->id = -1;

  // Deallocate face
  if (f->flags[R3_MESH_FACE_ALLOCATED]) delete f;
}



void R3Mesh::
DeleteVertex(R3MeshVertex *v)
{
  // Delete edges attached to vertex
  RNArray<R3MeshEdge *> vertex_edges = v->edges;
  for (int i = 0; i < vertex_edges.NEntries(); i++) 
    DeleteEdge(vertex_edges.Kth(i));

  // Deallocate vertex
  DeallocateVertex(v);
}



void R3Mesh::
DeleteEdge(R3MeshEdge *e)
{
  // Delete faces attached to edge
  if (e->face[0]) DeleteFace(e->face[0]);
  if (e->face[1]) DeleteFace(e->face[1]);

  // Remove edge from vertex arrays
  e->vertex[0]->edges.Remove(e);
  e->vertex[1]->edges.Remove(e);

  // Deallocate edge
  DeallocateEdge(e);
}



void R3Mesh::
DeleteFace(R3MeshFace *f)
{
  // Update edge-face relations
  if (f->edge[0]->face[0] == f) f->edge[0]->face[0] = NULL;
  if (f->edge[0]->face[1] == f) f->edge[0]->face[1] = NULL;
  if (f->edge[1]->face[0] == f) f->edge[1]->face[0] = NULL;
  if (f->edge[1]->face[1] == f) f->edge[1]->face[1] = NULL;
  if (f->edge[2]->face[0] == f) f->edge[2]->face[0] = NULL;
  if (f->edge[2]->face[1] == f) f->edge[2]->face[1] = NULL;

  // Deallocate face
  DeallocateFace(f);
}



void R3Mesh::
DeleteUnusedVertices(void)
{
  // Find unused vertices
  RNArray<R3MeshVertex *> unused_vertices;
  for (int i = 0; i < NVertices(); i++) {
    R3MeshVertex *vertex = Vertex(i);
    if (VertexValence(vertex) == 0) unused_vertices.Insert(vertex);
  }

  // Delete unused vertices
  for (int i = 0; i < unused_vertices.NEntries(); i++) {
    R3MeshVertex *vertex = unused_vertices.Kth(i);
    DeleteVertex(vertex);
  }
}



void R3Mesh::
DeleteUnusedEdges(void)
{
  // Find unused edges
  RNArray<R3MeshEdge *> unused_edges;
  for (int i = 0; i < NEdges(); i++) {
    R3MeshEdge *edge = Edge(i);
    if (!FaceOnEdge(edge, 0) && !FaceOnEdge(edge, 1)) unused_edges.Insert(edge);
  }

  // Delete unused edges
  for (int i = 0; i < unused_edges.NEntries(); i++) {
    R3MeshEdge *edge = unused_edges.Kth(i);
    DeleteEdge(edge);
  }
}



R3MeshVertex* R3Mesh::
MergeVertex(R3MeshVertex *v1, R3MeshVertex *v2)
{
#if 0
  // Check if there is an edge between the vertices
  // R3MeshEdge *e = EdgeBetweenVertices(v1, v2);
  // if (e) return CollapseEdge(e, v1->position);

  while (TRUE) {
    // Find a vertex v3 connected by edges to both v1 and v2
    R3MeshEdge *e1 = NULL;
    R3MeshEdge *e2 = NULL;
    R3MeshVertex *v3 = NULL;
    for (int i = 0; i < VertexValence(v1); i++) {
      e1 = EdgeOnVertex(v1, i);
      R3MeshVertex *ve1 = VertexAcrossEdge(e1, v1);
      for (int j = 0; j < VertexValence(v2); j++) {
        e2 = EdgeOnVertex(v2, j);
        R3MeshVertex *ve2 = VertexAcrossEdge(e2, v2);
        if (ve1 == ve2) { 
          v3 = ve1; 
          break; 
        }
      }
    }

    if (!v3) break;


    // Find all vertices in the triangle between v1, v2, and v3
    R3mesh_mark++;
    SetVertexMark(v1, R3mesh_mark);
    SetVertexMark(v2, R3mesh_mark);
    SetVertexMark(v3, R3mesh_mark);
    RNArray<R3MeshFace *> faces;
    R3MeshFace *f1 = FaceOnEdge(e1, v1, RN_CW);
    if (f1) faces.Insert(f1);
    R3MeshFace *f2 = FaceOnEdge(e2, v2, RN_CCW);
    if (f2) faces.Insert(f2);
    while (!faces.IsEmpty()) {
    }

  }

  // Get vertices on edge
  R3MeshVertex *v[2];
  v[0] = edge->vertex[0];
  v[1] = edge->vertex[1];

  // Get faces on edge
  R3MeshFace *f[2];
  f[0] = edge->face[0];
  f[1] = edge->face[1];

  // Get other edges and vertex on faces
  R3MeshEdge *ef[2][2] = { { NULL, NULL }, { NULL, NULL } };
  R3MeshFace *ff[2][2] = { { NULL, NULL }, { NULL, NULL } };
  R3MeshVertex *vf[2] = { NULL, NULL };
  if (f[0]) {
    ef[0][0] = EdgeAcrossVertex(v[0], edge, f[0]);
    ef[0][1] = EdgeAcrossVertex(v[1], edge, f[0]);
    vf[0] = VertexAcrossEdge(ef[0][1], v[1]);
    assert(vf[0] == VertexAcrossEdge(ef[0][0], v[0]));
    ff[0][0] = FaceAcrossEdge(ef[0][0], f[0]);
    ff[0][1] = FaceAcrossEdge(ef[0][1], f[0]);
  }
  if (f[1]) {
    ef[1][0] = EdgeAcrossVertex(v[0], edge, f[1]);
    ef[1][1] = EdgeAcrossVertex(v[1], edge, f[1]);
    vf[1] = VertexAcrossEdge(ef[1][1], v[1]);
    assert(vf[1] == VertexAcrossEdge(ef[1][0], v[0]));
    ff[1][0] = FaceAcrossEdge(ef[1][0], f[1]);
    ff[1][1] = FaceAcrossEdge(ef[1][1], f[1]);
  }

  // Make an array of vertices attached by edges to both v1 and v2
  RNArray<R3MeshVertices *> adjacent_vertices;
  for (int i = 0; i < VertexValence(v1); i++) {
    R3MeshEdge *e1 = EdgeOnVertex(v1, i);
    R3MeshVertex *ve1 = VertexAcrossEdge(e1, v1);
    for (int j = 0; j < VertexValence(v2); j++) {
      R3MeshEdge *e2 = EdgeOnVertex(v2, j);
      R3MeshVertex *ve2 = VertexAcrossEdge(e2, v2);
      if (ve1 == ve2) {
        // Check for fin
        for (int k1 = 0; k1 < 2; k1++) {
          R3MeshFace *f1 = FaceOnEdge(e1, k1);
          if (f1) {
            R3MeshVertex *vf1 = VertexAcrossFace(f1, e1);
            for (int k2 = 0; k2 < 2; k2++) {
              R3MeshFace *f2 = FaceOnEdge(e2, k2);
              if (f2) {
                R3MeshVertex *vf2 = VertexAcrossFace(f2, e2);
                if (vf1 == vf2) return NULL;
              }
            }
          }
        }

        // Remember adjacent vertex
        adjacent_vertices.Insert(ve1);
      }
    }
  }

  // Check if vertices are connected
  if (adjacent_vertices.IsEmpty()) {
 }

  // v1 and v2 are not connected, merge references to/from v2 into v1
  for (int i = 0; i < v2->edges.NEntries(); i++) {
    R3MeshEdge *e2 = v2->edges.Kth(i);

    // Update edge-vertex relations 
    for (int j = 0; j < 2; j++) {
      if (e2->vertex[j] == v2) e2->vertex[j] = v1;
    }

    // Update face-vertex relations 
    for (int j = 0; j < 2; j++) {
      R3MeshEdge *f2 = e2->face[j];
      for (int k = 0; k < 3; k++) {
        if (f2->vertex[k] == v2) f2->vertex[k] = v1;
      }
    }

    // Update vertex-edge relations 
    v1->edges.Insert(e2);
  }
 
  // Deallocate vertex
  DeallocateVertex(v2);

  // Return remaining vertex
  return v1;
#else
  RNAbort("Not implemented");
  return NULL;
#endif
}



void R3Mesh::
MergeCoincidentVertices(RNLength epsilon)
{
  // Compute epsilon
  if (epsilon < 0.0) epsilon = 0.0001 * bbox.DiagonalLength();

  // Load all vertices into a regular grid
  const int NUM_GRID_CELLS = 32;
  RNArray<R3MeshVertex *> grid[NUM_GRID_CELLS][NUM_GRID_CELLS][NUM_GRID_CELLS];
  for (int i = 0; i < vertices.NEntries(); i++) {
    R3MeshVertex *v = vertices[i];
    const R3Point& p = v->position;
    int ix = (int) (NUM_GRID_CELLS * (p.X() - bbox.XMin()) / bbox.XLength());
    int iy = (int) (NUM_GRID_CELLS * (p.Y() - bbox.YMin()) / bbox.YLength());
    int iz = (int) (NUM_GRID_CELLS * (p.Z() - bbox.ZMin()) / bbox.ZLength());
    if (ix < 0) ix = 0; else if (ix > NUM_GRID_CELLS-1) ix = NUM_GRID_CELLS-1;
    if (iy < 0) iy = 0; else if (iy > NUM_GRID_CELLS-1) iy = NUM_GRID_CELLS-1;
    if (iz < 0) iz = 0; else if (iz > NUM_GRID_CELLS-1) iz = NUM_GRID_CELLS-1;
    grid[ix][iy][iz].Insert(v);
  }

  // Merge coincident vertices in same grid cell
  for (int ix = 0; ix < NUM_GRID_CELLS; ix++) {
    for (int iy = 0; iy < NUM_GRID_CELLS; iy++) {
      for (int iz = 0; iz < NUM_GRID_CELLS; iz++) {
        // Consider all pairs of vertices in same cell
        for (int j = 0; j < grid[ix][iy][iz].NEntries(); j++) {
          R3MeshVertex *v1 = grid[ix][iy][iz].Kth(j);
          for (int k = j+1; k < grid[ix][iy][iz].NEntries(); k++) {
            R3MeshVertex *v2 = grid[ix][iy][iz].Kth(k);
            // Check if vertices are coincident
            if (R3Contains(v1->position, v2->position)) {
              // Merge coincident vertices
              MergeVertex(v1, v2);
            }
          }
        }
      }
    }
  }
}



R3MeshVertex *R3Mesh::
CollapseEdge(R3MeshEdge *edge, const R3Point& point)
{
  // Just to be safe
  if (NFaces() <= 4) return NULL;

  // Get vertices on edge
  R3MeshVertex *v[2];
  v[0] = edge->vertex[0];
  v[1] = edge->vertex[1];

  // Get faces on edge
  R3MeshFace *f[2];
  f[0] = edge->face[0];
  f[1] = edge->face[1];

#if 1
  // Get other edges and vertex on faces
  R3MeshEdge *ef[2][2] = { { NULL, NULL }, { NULL, NULL } };
  R3MeshFace *ff[2][2] = { { NULL, NULL }, { NULL, NULL } };
  R3MeshVertex *vf[2] = { NULL, NULL };
  if (f[0]) {
    ef[0][0] = EdgeAcrossVertex(v[0], edge, f[0]);
    ef[0][1] = EdgeAcrossVertex(v[1], edge, f[0]);
    vf[0] = VertexAcrossEdge(ef[0][1], v[1]);
    assert(vf[0] == VertexAcrossEdge(ef[0][0], v[0]));
    ff[0][0] = FaceAcrossEdge(ef[0][0], f[0]);
    ff[0][1] = FaceAcrossEdge(ef[0][1], f[0]);
  }
  if (f[1]) {
    ef[1][0] = EdgeAcrossVertex(v[0], edge, f[1]);
    ef[1][1] = EdgeAcrossVertex(v[1], edge, f[1]);
    vf[1] = VertexAcrossEdge(ef[1][1], v[1]);
    assert(vf[1] == VertexAcrossEdge(ef[1][0], v[0]));
    ff[1][0] = FaceAcrossEdge(ef[1][0], f[1]);
    ff[1][1] = FaceAcrossEdge(ef[1][1], f[1]);
  }

  // Check if will create a fin
  for (int i = 0; i < VertexValence(v[0]); i++) {
    R3MeshEdge *e0 = EdgeOnVertex(v[0], i);
    if ((e0 == ef[0][0]) || (e0 == ef[1][0])) continue;
    R3MeshVertex *ve0 = VertexAcrossEdge(e0, v[0]);
    for (int j = 0; j < VertexValence(v[1]); j++) {
      R3MeshEdge *e1 = EdgeOnVertex(v[1], j);
      if ((e1 == ef[0][1]) || (e1 == ef[1][1])) continue;
      R3MeshVertex *ve1 = VertexAcrossEdge(e1, v[1]);
      if (ve0 == ve1) return NULL;
    }
  }

  // Update vertex-edge relations
  v[0]->edges.Remove(edge);
  for (int i = 0; i < v[1]->edges.NEntries(); i++) {
    R3MeshEdge *ev1 = v[1]->edges[i];
    if (ev1 == edge) continue;
    if (ev1 == ef[0][1]) continue;
    if (ev1 == ef[1][1]) continue;
    v[0]->edges.Insert(ev1);
  }
  if (f[0]) {
    assert(vf[0] && ef[0][1]);
    vf[0]->edges.Remove(ef[0][1]);
  }
  if (f[1]) {
    assert(vf[1] && ef[1][1]);
    vf[1]->edges.Remove(ef[1][1]);
  }

  // Update edge-vertex relations (v[1]->v[0]) 
  for (int i = 0; i < v[1]->edges.NEntries(); i++) {
    R3MeshEdge *ev = v[1]->edges[i];
    for (int j = 0; j < 2; j++) {
      if (ev->vertex[j] == v[1]) ev->vertex[j] = v[0];
    }
  }
  
  // Update edge-face relations (ef[i][1]->ef[i][0])
  for (int i = 0; i < 2; i++) {
    if (!f[i]) continue;
    assert(ef[i][0]);
    for (int j = 0; j < 2; j++) {
      if (ef[i][0]->face[j] == f[i]) ef[i][0]->face[j] = ff[i][1];
    }
  }

  // Update face-vertex relations (v[1]->v[0])
  for (int i = 0; i < v[1]->edges.NEntries(); i++) {
    R3MeshEdge *ev = v[1]->edges[i];
    for (int j = 0; j < 2; j++) {
      R3MeshFace *fv = ev->face[j];
      if (fv) {
        for (int k = 0; k < 3; k++) {
          if (fv->vertex[k] == v[1]) fv->vertex[k] = v[0];
        }
      }
    }
  }
  
  // Update face-edge relations (ef[i][1]->ef[i][0])
  for (int i = 0; i < 2; i++) {
    if (!f[i]) continue;
    if (!ff[i][1]) continue;
    for (int j = 0; j < 3; j++) {
      if (ff[i][1]->edge[j] == ef[i][1]) {
        ff[i][1]->edge[j] = ef[i][0];
      }
    }
  }

  // Update vertex properties
  RNScalar d0 = R3Distance(point, VertexPosition(v[0]));
  RNScalar d1 = R3Distance(point, VertexPosition(v[1]));
  RNScalar dsum = d0 + d1;
  if (RNIsPositive(dsum)) {
    d0 /= dsum;  d1 /= dsum;
    SetVertexTextureCoords(v[0], d1*VertexTextureCoords(v[0]) + d0*VertexTextureCoords(v[1])); 
    SetVertexColor(v[0], d1*VertexColor(v[0]) + d0*VertexColor(v[1])); 
  }

  // Deallocate faces 
  if (f[0]) DeallocateFace(f[0]);
  if (f[1]) DeallocateFace(f[1]);

  // Deallocate edges
  if (f[0]) DeallocateEdge(ef[0][1]);
  if (f[1]) DeallocateEdge(ef[1][1]);
  DeallocateEdge(edge);

  // Deallocate vertex
  DeallocateVertex(v[1]);

  // Set position of remaining vertex
  SetVertexPosition(v[0], point);

  // Return remaining vertex
  return v[0];
#else
  // Set vertex position 
  SetVertexPosition(v[0], point);

  // Make arrays of everything attached to v[1]
  R3mesh_mark++;
  RNArray<R3MeshFace *> faces_to_remove;
  RNArray<R3MeshEdge *> edges_to_remove;
  RNArray<R3MeshVertex *> cw_vertices;
  RNArray<R3MeshVertex *> ccw_vertices;
  for (int i = 0; i < v[1]->edges.NEntries(); i++) {
    R3MeshEdge *v1_e = v[1]->edges[i];
    edges_to_remove.Insert(v1_e);
    R3MeshFace *v1_f = FaceOnEdge(v1_e, v[1], RN_CCW);
    if (v1_f) {
      faces_to_remove.Insert(v1_f);
      if ((v1_f != f[0]) && (v1_f != f[1])) {
        cw_vertices.Insert(VertexOnFace(v1_f, v[1], RN_CW));
        ccw_vertices.Insert(VertexOnFace(v1_f, v[1], RN_CCW));
      }
    }
  }

  // Remove everything attached to v1
  for (int i = 0; i < faces_to_remove.NEntries(); i++) 
    DeleteFace(faces_to_remove[i]);
  for (int i = 0; i < edges_to_remove.NEntries(); i++) 
    DeleteEdge(edges_to_remove[i]);
  DeleteVertex(v[1]);

  // Recreate everything by replacing v1 with v0
  for (int i = 0; i < cw_vertices.NEntries(); i++) {
    CreateFace(v[0], cw_vertices[i], ccw_vertices[i]);
  }

  // Return remaining vertex
  return v[0];
#endif
}



R3MeshVertex *R3Mesh::
CollapseEdge(R3MeshEdge *edge)
{
  // Collapse edge and put new vertex at midpoint
  return CollapseEdge(edge, EdgeMidpoint(edge));
}



R3MeshVertex *R3Mesh::
CollapseFace(R3MeshFace *f, const R3Point& point)
{
  R3MeshVertex *v0 = VertexOnFace(f, 0);
  R3MeshVertex *v1 = VertexOnFace(f, 1);
  R3MeshVertex *v2 = VertexOnFace(f, 2);
  R3MeshEdge *e01 = EdgeBetweenVertices(v0, v1);
  R3MeshVertex *v01 = CollapseEdge(e01, point);
  if (!v01) return NULL;
  R3MeshEdge *e012 = EdgeBetweenVertices(v01, v2);
  assert(e012);
  R3MeshVertex *v012 = CollapseEdge(e012, point);
  return v012;
}



R3MeshVertex *R3Mesh::
CollapseFace(R3MeshFace *face)
{
  // Collapse face and put new vertex at centroid
  return CollapseFace(face, FaceCentroid(face));
}



R3MeshVertex *R3Mesh::
SplitEdge(R3MeshEdge *edge, const R3Point& point, R3MeshEdge **e0, R3MeshEdge **e1)
{
  // Create vertex at split point
  R3MeshVertex *vertex = CreateVertex(R3zero_point);

  // Get vertices on edge
  R3MeshVertex *v[2];
  v[0] = edge->vertex[0];
  v[1] = edge->vertex[1];

  // Get faces on edge
  R3MeshFace *f[2];
  f[0] = edge->face[0];
  f[1] = edge->face[1];

  // Get other edges and vertex on f[0]
  R3MeshEdge *ef[2][2] = { { NULL, NULL }, { NULL, NULL } };
  R3MeshVertex *vf[2] = { NULL, NULL };
  int m[2] = { 0, 0 };
  int s[2] = { 0, 0 };
  int c[2] = { 0, 0 };
  if (f[0]) {
    ef[0][0] = EdgeAcrossVertex(v[0], edge, f[0]);
    ef[0][1] = EdgeAcrossVertex(v[1], edge, f[0]);
    vf[0] = VertexAcrossEdge(ef[0][1], v[1]);
    assert(vf[0] == VertexAcrossEdge(ef[0][0], v[0]));
    m[0] = FaceMaterial(f[0]);
    s[0] = FaceSegment(f[0]);
    c[0] = FaceCategory(f[0]);
  }
  if (f[1]) {
    ef[1][0] = EdgeAcrossVertex(v[0], edge, f[1]);
    ef[1][1] = EdgeAcrossVertex(v[1], edge, f[1]);
    vf[1] = VertexAcrossEdge(ef[1][1], v[1]);
    assert(vf[1] == VertexAcrossEdge(ef[1][0], v[0]));
    m[1] = FaceMaterial(f[1]);
    s[1] = FaceSegment(f[1]);
    c[1] = FaceCategory(f[1]);
  }

  // Delete the edge and the two adjacent faces
  if (f[0]) DeleteFace(f[0]);
  if (f[1]) DeleteFace(f[1]);
  DeleteEdge(edge);

  // Create two new edges
  R3MeshEdge *e[2];
  e[0] = CreateEdge(v[0], vertex);
  e[1] = CreateEdge(vertex, v[1]);

  // Create two new faces on one side
  if (f[0]) {
    // Create edge splitting f[0]
    R3MeshEdge *es = CreateEdge(vertex, vf[0]);

    // Create new face on v[0] side of split
    assert((vf[0] != v[0]) && (vf[0] != vertex) && (v[0] != vertex));
    assert((ef[0][0] != e[0]) && (ef[0][0] != es) && (e[0] != es));
    R3MeshFace *f00 = CreateFace(vf[0], v[0], vertex, ef[0][0], e[0], es);
    SetFaceMaterial(f00, m[0]);
    SetFaceSegment(f00, s[0]);
    SetFaceCategory(f00, c[0]);

    // Create new face on v[1] side of split
    assert((es != e[1]) && (es != ef[0][1]) && (e[1] != ef[0][1]));
    R3MeshFace *f01 = CreateFace(vf[0], vertex, v[1], es, e[1], ef[0][1]);
    SetFaceMaterial(f01, m[0]);
    SetFaceSegment(f01, s[0]);
    SetFaceCategory(f01, c[0]);
  }

  // Create two new faces on other side
  if (f[1]) {
    // Create edge splitting f[1]
    R3MeshEdge *es = CreateEdge(vertex, vf[1]);
    
    // Create new face on v[1] side of split
    assert((vf[1] != v[1]) && (vf[1] != vertex) && (v[1] != vertex));
    assert((ef[1][1] != e[1]) && (ef[1][1] != es) && (e[1] != es));
    R3MeshFace *f10 = CreateFace(vf[1], v[1], vertex, ef[1][1], e[1], es);
    SetFaceMaterial(f10, m[1]);
    SetFaceSegment(f10, s[1]);
    SetFaceCategory(f10, c[1]);

    // Create new face on v[0] side of split
    assert((es != e[0]) && (es != ef[1][0]) && (e[0] != ef[1][0]));
    R3MeshFace *f11 = CreateFace(vf[1], vertex, v[0], es, e[0], ef[1][0]);
    SetFaceMaterial(f11, m[1]);
    SetFaceSegment(f11, s[1]);
    SetFaceCategory(f11, c[1]);
  }

  // Interpolate vertex properties
  RNScalar d0 = R3Distance(point, VertexPosition(v[0]));
  RNScalar d1 = R3Distance(point, VertexPosition(v[1]));
  RNScalar dsum = d0 + d1;
  if (RNIsPositive(dsum)) {
    d0 /= dsum;  d1 /= dsum;
    SetVertexTextureCoords(vertex, d1*VertexTextureCoords(v[0]) + d0*VertexTextureCoords(v[1])); 
    SetVertexColor(vertex, d1*VertexColor(v[0]) + d0*VertexColor(v[1])); 
  }

  // Set vertex position (also marks things in need of update)
  SetVertexPosition(vertex, point);

  // Return edges
  if (e0) *e0 = e[0];
  if (e1) *e1 = e[1];

  // Return vertex
  return vertex;
}



R3MeshVertex *R3Mesh::
SplitEdge(R3MeshEdge *e, const R3Plane& plane)
{
  // Compute intersection point with plane
  RNScalar t;
  R3Point point;
  R3Span edge_span = EdgeSpan(e);
  if (!R3Intersects(edge_span, plane, &point, &t)) return NULL;
  if ((RNIsEqual(t, 0.0)) || (RNIsEqual(t, edge_span.Length()))) return NULL;

  // Split edge
  return SplitEdge(e, point);
}



R3MeshVertex *R3Mesh::
SubdivideEdge(R3MeshEdge *e)
{
  // Split edge near midpoint
  return SplitEdge(e, EdgeMidpoint(e));
}



R3MeshVertex *R3Mesh::
SplitFace(R3MeshFace *f, const R3Point& point, R3MeshFace **f0, R3MeshFace **f1, R3MeshFace **f2)
{
  // Find vertices/edges bounding face
  R3MeshVertex *v0 = VertexOnFace(f, 0);
  R3MeshEdge *e0 = EdgeOnFace(f, v0, RN_CCW);
  R3MeshVertex *v1 = VertexAcrossEdge(e0, v0);
  R3MeshEdge *e1 = EdgeOnFace(f, v1, RN_CCW);
  R3MeshVertex *v2 = VertexAcrossEdge(e1, v1);
  R3MeshEdge *e2 = EdgeOnFace(f, v2, RN_CCW);
  int m = FaceMaterial(f);
  int s = FaceSegment(f);
  int c = FaceCategory(f);

  // Delete face
  DeleteFace(f);

  // Create new vertex at point
  R3MeshVertex *vertex = CreateVertex(point);

  // Create three new edges
  R3MeshEdge *s0 = CreateEdge(v0, vertex);
  R3MeshEdge *s1 = CreateEdge(v1, vertex);
  R3MeshEdge *s2 = CreateEdge(v2, vertex);

  // Create three new faces
  R3MeshFace *t0 = CreateFace(vertex, v0, v1, s0, e0, s1);
  R3MeshFace *t1 = CreateFace(vertex, v1, v2, s1, e1, s2);
  R3MeshFace *t2 = CreateFace(vertex, v2, v0, s2, e2, s0);
  
  // Set face materials
  SetFaceMaterial(t0, m);
  SetFaceMaterial(t1, m);
  SetFaceMaterial(t2, m);

  // Set face segments
  SetFaceSegment(t0, s);
  SetFaceSegment(t1, s);
  SetFaceSegment(t2, s);

  // Set face categories
  SetFaceCategory(t0, c);
  SetFaceCategory(t1, c);
  SetFaceCategory(t2, c);

  // Interpolate vertex properties
  RNScalar d0 = R3Distance(point, VertexPosition(v0));
  RNScalar d1 = R3Distance(point, VertexPosition(v1));
  RNScalar d2 = R3Distance(point, VertexPosition(v2));
  RNScalar dsum = d0 + d1 + d2;
  if (RNIsPositive(dsum)) {
    RNScalar t0 = (d1+d2)/dsum;   RNScalar t1 = (d0+d2)/dsum;  RNScalar t2 = (d0+d1)/dsum;
    SetVertexTextureCoords(vertex, t0*VertexTextureCoords(v0) + t1*VertexTextureCoords(v1) + t2*VertexTextureCoords(v2)); 
    SetVertexColor(vertex, t0*VertexColor(v0) + t1*VertexColor(v1) + t2*VertexColor(v2)); 
  }

  // Return created faces
  if (f0) *f0 = t0;
  if (f1) *f1 = t1;
  if (f2) *f2 = t2;

  // Return new vertex
  return vertex;
}



R3MeshFace *R3Mesh::
SubdivideFace(R3MeshFace *f)
{
  // Find vertices/edges bounding face
  R3MeshVertex *v0 = VertexOnFace(f, 0);
  R3MeshEdge *e0 = EdgeOnFace(f, v0, RN_CCW);
  R3MeshVertex *v1 = VertexAcrossEdge(e0, v0);
  R3MeshEdge *e1 = EdgeOnFace(f, v1, RN_CCW);
  R3MeshVertex *v2 = VertexAcrossEdge(e1, v1);
  R3MeshEdge *e2 = EdgeOnFace(f, v2, RN_CCW);
  int m = FaceMaterial(f);
  int s = FaceSegment(f);
  int c = FaceCategory(f);

  // Delete face
  DeleteFace(f);

  // Subdivide edges
  R3MeshVertex *ve0 = SubdivideEdge(e0);
  R3MeshVertex *ve1 = SubdivideEdge(e1);
  R3MeshVertex *ve2 = SubdivideEdge(e2);

  // Create new faces
  R3MeshFace *f1 = CreateFace(v0, ve0, ve2);  
  R3MeshFace *f2 = CreateFace(v1, ve1, ve0);  
  R3MeshFace *f3 = CreateFace(v2, ve2, ve1);  
  R3MeshFace *f4 = CreateFace(ve0, ve1, ve2);  

  // Set face materials
  SetFaceMaterial(f1, m);
  SetFaceMaterial(f2, m);
  SetFaceMaterial(f3, m);
  SetFaceMaterial(f4, m);

  // Set face segments
  SetFaceSegment(f1, s);
  SetFaceSegment(f2, s);
  SetFaceSegment(f3, s);
  SetFaceSegment(f4, s);

  // Set face categorys
  SetFaceCategory(f1, c);
  SetFaceCategory(f2, c);
  SetFaceCategory(f3, c);
  SetFaceCategory(f4, c);

  // Return interior face
  return f4;
}



int R3Mesh::
SwapEdge(R3MeshEdge *edge)
{
  // Get vertices on edge
  R3MeshVertex *v[2];
  v[0] = edge->vertex[0];
  v[1] = edge->vertex[1];

  // Get faces on edge
  R3MeshFace *f[2];
  f[0] = edge->face[0];
  f[1] = edge->face[1];
  if (!f[0] || !f[1]) return 0;

  // Get edges and vertices across faces
  R3MeshEdge *ef[2][2] = { { NULL, NULL }, { NULL, NULL } };
  R3MeshVertex *vf[2] = { NULL, NULL };
  ef[0][0] = EdgeAcrossVertex(v[0], edge, f[0]);
  ef[0][1] = EdgeAcrossVertex(v[1], edge, f[0]);
  vf[0] = VertexAcrossEdge(ef[0][1], v[1]);
  assert(vf[0] == VertexAcrossEdge(ef[0][0], v[0]));
  ef[1][0] = EdgeAcrossVertex(v[0], edge, f[1]);
  ef[1][1] = EdgeAcrossVertex(v[1], edge, f[1]);
  vf[1] = VertexAcrossEdge(ef[1][1], v[1]);
  assert(vf[1] == VertexAcrossEdge(ef[1][0], v[0]));
  if (EdgeBetweenVertices(vf[0], vf[1])) return 0;

  // Update vertices
  v[0]->edges.Remove(edge);
  v[1]->edges.Remove(edge);
  vf[0]->edges.Insert(edge);
  vf[1]->edges.Insert(edge);

  // Update edge
  edge->vertex[0] = vf[0];
  edge->vertex[1] = vf[1];
  edge->length = 0;
  edge->flags = 0;

  // Update faces
  UpdateFaceRefs(f[0], vf[0], vf[1], v[1], edge, ef[1][1], ef[0][1]);
  UpdateFacePlane(f[0]);
  UpdateFaceBBox(f[0]);
  UpdateFaceRefs(f[1], vf[1], vf[0], v[0], edge, ef[0][0], ef[1][0]);
  UpdateFacePlane(f[1]);
  UpdateFaceBBox(f[1]);

  // Return success
  return 1;
}



void R3Mesh::
FlipEdge(R3MeshEdge *e)
{
  // Reverse order of vertices
  R3MeshVertex *v = e->vertex[0];
  e->vertex[0] = e->vertex[1];
  e->vertex[1] = v;

  // Reverse order of faces
  R3MeshFace *f = e->face[0];
  e->face[0] = e->face[1];
  e->face[1] = f;
}



void R3Mesh::
FlipFace(R3MeshFace *f)
{
  // Reverse orientation of plane
  f->plane.Flip();

  // Reverse order of vertices
  R3MeshVertex *vswap = f->vertex[0];
  f->vertex[0] = f->vertex[2];
  f->vertex[2] = vswap;

  // Reverse order of edges
  R3MeshEdge *eswap = f->edge[0];
  f->edge[0] = f->edge[1];
  f->edge[1] = eswap;
}



void R3Mesh::
FlipFaces(void)
{
  // Reverse orientation of all faces
  for (int i = 0; i < NFaces(); i++) {
    R3MeshFace *f = Face(i);
    FlipFace(f);
  }
}



static RNScalar 
EdgeLengthCallback(R3MeshEdge *edge, void *data)
{
  // Return value associated with edge for heap sorting
  const R3Mesh *mesh = (R3Mesh *) data;
  return mesh->EdgeLength(edge);
}



void R3Mesh::
CollapseShortEdges(RNLength min_edge_length)
{
  // Create priority queue for short edges
  RNHeap<R3MeshEdge *> heap(EdgeLengthCallback, NULL, this, TRUE);
  R3mesh_mark++;

  // Add short edges to priority queue
  for (int i = 0; i < NEdges(); i++) {
    R3MeshEdge *edge = Edge(i);
    RNLength length = EdgeLength(edge);
    if (length < min_edge_length) {
      SetEdgeMark(edge, R3mesh_mark);
      heap.Push(edge);
    }
  }

  // Collapse edges in shortest to longest order
  while (!heap.IsEmpty()) {
    R3MeshEdge *edge = heap.Pop();
    assert(EdgeMark(edge) == R3mesh_mark);
    assert(IsEdgeOnMesh(edge));

    // Remove edges to be deleted from queue (e01, e11)
    R3MeshVertex *v1 = VertexOnEdge(edge, 1);
    R3MeshFace *f0 = FaceOnEdge(edge, 0);  
    R3MeshFace *f1 = FaceOnEdge(edge, 1);  
    R3MeshEdge *e01 = (f0) ? EdgeAcrossVertex(v1, edge, f0) : NULL;
    R3MeshEdge *e11 = (f1) ? EdgeAcrossVertex(v1, edge, f1) : NULL;
    if (e01 && (EdgeMark(e01) == R3mesh_mark)) { heap.Remove(e01); SetEdgeMark(e01, 0); }
    if (e11 && (EdgeMark(e11) == R3mesh_mark)) { heap.Remove(e11); SetEdgeMark(e11, 0); }

    // Collapse edge
    R3MeshVertex *vertex = CollapseEdge(edge);
    if (!vertex) continue;

    // Update adjacent edges
    for (int i = 0; i < VertexValence(vertex); i++) {
      R3MeshEdge *edge = EdgeOnVertex(vertex, i);
      if (EdgeMark(edge) == R3mesh_mark) {
        heap.Update(edge);
      }
      else {
        RNLength length = EdgeLength(edge);
        if (length < min_edge_length) {
          SetEdgeMark(edge, R3mesh_mark);
          heap.Push(edge);
        }
      }
    }
  }
}



void R3Mesh::
SubdivideLongEdges(RNLength max_edge_length)
{
  // Create priority queue for long edges
  RNHeap<R3MeshEdge *> heap(EdgeLengthCallback, NULL, this, FALSE);

  // Add long edges to priority queue
  for (int i = 0; i < NEdges(); i++) {
    R3MeshEdge *edge = Edge(i);
    RNLength length = EdgeLength(edge);
    if (length > max_edge_length) heap.Push(edge);
  }

  // Split edges in longest to shortest order
  while (!heap.IsEmpty()) {
    R3MeshEdge *edge = heap.Pop();
    R3MeshVertex *vertex = SubdivideEdge(edge);
    for (int i = 0; i < VertexValence(vertex); i++) {
      R3MeshEdge *edge = EdgeOnVertex(vertex, i);
      RNLength length = EdgeLength(edge);
      if (length > max_edge_length) heap.Push(edge);
    }
  }
}



void R3Mesh::
SubdivideFaces(void)
{
  // Remember cardinalities
  int nfaces = NFaces();
  int nedges = NEdges();
  int nvertices = NVertices();

  // Allocate temporary memory for vertex indices
  int *vertex_indices = new int [ nfaces * 4 * 3 ];

  // Compute vertex indices for new faces
  for (int i = 0; i < nfaces; i++) {
    R3MeshFace *face = Face(i);
    R3MeshVertex *v0 = VertexOnFace(face, 0);
    R3MeshVertex *v1 = VertexOnFace(face, 1);
    R3MeshVertex *v2 = VertexOnFace(face, 2);
    R3MeshEdge *e0 = EdgeOnFace(face, 0);
    R3MeshEdge *e1 = EdgeOnFace(face, 1);
    R3MeshEdge *e2 = EdgeOnFace(face, 2);
    vertex_indices[i * 12 + 0] = VertexID(v0);
    vertex_indices[i * 12 + 1] = nvertices + EdgeID(e0);
    vertex_indices[i * 12 + 2] = nvertices + EdgeID(e2);
    vertex_indices[i * 12 + 3] = VertexID(v1);
    vertex_indices[i * 12 + 4] = nvertices + EdgeID(e1);
    vertex_indices[i * 12 + 5] = nvertices + EdgeID(e0);
    vertex_indices[i * 12 + 6] = VertexID(v2);
    vertex_indices[i * 12 + 7] = nvertices + EdgeID(e2);
    vertex_indices[i * 12 + 8] = nvertices + EdgeID(e1);
    vertex_indices[i * 12 + 9] = nvertices + EdgeID(e0);
    vertex_indices[i * 12 + 10] = nvertices + EdgeID(e1);
    vertex_indices[i * 12 + 11] = nvertices + EdgeID(e2);
  }

  // Create vertex at midpoint of every edge
  for (int i = 0; i < nedges; i++) {
    R3MeshEdge *edge = Edge(i);
    R3MeshVertex *v0 = VertexOnEdge(edge, 0);
    R3MeshVertex *v1 = VertexOnEdge(edge, 1);
    R3MeshVertex *vertex = CreateVertex(EdgeMidpoint(edge));
    SetVertexTextureCoords(vertex, 0.5*VertexTextureCoords(v0) + 0.5*VertexTextureCoords(v1)); 
    SetVertexColor(vertex, 0.5*VertexColor(v0) + 0.5*VertexColor(v1)); 
  }

  // Delete all edges and faces
  while (NFaces() > 0) DeleteFace(Face(0));
  while (NEdges() > 0) DeleteEdge(Edge(0));

  // Create new edges and faces
  for (int i = 0; i < nfaces; i++) {
    for (int j = 0; j < 4; j++) {
      R3MeshVertex *v0 = Vertex(vertex_indices[i*12 + j*3 + 0]);
      R3MeshVertex *v1 = Vertex(vertex_indices[i*12 + j*3 + 1]);
      R3MeshVertex *v2 = Vertex(vertex_indices[i*12 + j*3 + 2]);
      CreateFace(v0, v1, v2);
    }
  }

  // Delete temporary memory for vertex indices
  delete [] vertex_indices;
}



////////////////////////////////////////////////////////////////////////
// EDGE SWAPPING FUNCTIONS
////////////////////////////////////////////////////////////////////////

static RNScalar 
R3MeshEdgeSwapValue(R3Mesh *mesh, R3MeshEdge *e)
{
  // Get/check faces
  R3MeshFace *f0 = mesh->FaceOnEdge(e, 0);
  R3MeshFace *f1 = mesh->FaceOnEdge(e, 1);
  if (!f0 || !f1) return -RN_INFINITY;

  // Get/check vertices
  R3MeshVertex *v0 = mesh->VertexOnEdge(e, 0);
  R3MeshVertex *v1 = mesh->VertexOnEdge(e, 1);
  R3MeshVertex *vf0 = mesh->VertexAcrossFace(f0, e);
  R3MeshVertex *vf1 = mesh->VertexAcrossFace(f1, e);
  if (mesh->EdgeBetweenVertices(vf0, vf1)) return -RN_INFINITY;

  // Check if swap would flip face normals
  R3Vector n0 = mesh->FaceNormal(f0);
  R3Vector n1 = mesh->FaceNormal(f1);
  const R3Point& p0 = mesh->VertexPosition(v0);
  const R3Point& p1 = mesh->VertexPosition(v1);
  const R3Point& pf0 = mesh->VertexPosition(vf0);
  const R3Point& pf1 = mesh->VertexPosition(vf1);
  R3Vector vec00 = pf0 - p0;
  R3Vector vec01 = pf1 - p0;
  R3Vector vec10 = pf0 - p1;
  R3Vector vec11 = pf1 - p1;
  R3Vector nf0 = vec01 % vec00;
  R3Vector nf1 = vec10 % vec11;
  if (RNIsNegativeOrZero(nf0.Dot(n0))) return -RN_INFINITY;
  if (RNIsNegativeOrZero(nf0.Dot(n1))) return -RN_INFINITY;
  if (RNIsNegativeOrZero(nf1.Dot(n0))) return -RN_INFINITY;
  if (RNIsNegativeOrZero(nf1.Dot(n1))) return -RN_INFINITY;

  // Compute interior angles at vertices
  R3Vector vecA = p1 - p0;
  R3Vector vecB = pf1 - pf0;
  RNAngle A00 = R3InteriorAngle(vecA, vec00);
  RNAngle A01 = R3InteriorAngle(vecA, vec01);
  RNAngle A10 = R3InteriorAngle(-vecA, vec10);
  RNAngle A11 = R3InteriorAngle(-vecA, vec11);
  RNAngle B00 = R3InteriorAngle(-vecB, vec00);
  RNAngle B01 = R3InteriorAngle(vecB, vec01);
  RNAngle B10 = R3InteriorAngle(-vecB, vec10);
  RNAngle B11 = R3InteriorAngle(vecB, vec11);

  // Determine minimum interior angles
  RNAngle A = A00;
  if (A01 < A) A = A01;
  if (A10 < A) A = A10;
  if (A11 < A) A = A11;
  RNAngle B = B00;
  if (B01 < B) B = B01;
  if (B10 < B) B = B10;
  if (B11 < B) B = B11;
  if (RNIsZero(B)) return -RN_INFINITY;

  // Compute difference in minimum interior angles
  RNAngle delta = B - A;

  // Return difference in interior angles
  return delta;
}



void R3Mesh::
SwapEdges(RNAngle min_angle_improvement)
{
  // Initialize edge values
  RNHeap<RNScalar *> heap(0, -1, FALSE);
  RNScalar *edge_values = new RNScalar [ NEdges() ];
  for (int i = 0; i < NEdges(); i++) {
    R3MeshEdge *edge = Edge(i);
    edge_values[i] = R3MeshEdgeSwapValue(this, edge);
    if (edge_values[i] < min_angle_improvement) continue;
    heap.Push(&edge_values[i]);
  }

  // Swap edges
  while (!heap.IsEmpty()) {
    // Pop edge
    RNScalar *edge_value = heap.Pop();
    int edge_index = edge_value - edge_values;
    R3MeshEdge *edge = Edge(edge_index);
    if (*edge_value < min_angle_improvement) break;

    // Check current edge value
    RNScalar current_value = R3MeshEdgeSwapValue(this, edge);
    if (current_value != edge_values[edge_index]) {
      edge_values[edge_index] = current_value;
      if (current_value >= min_angle_improvement) heap.Push(&edge_values[edge_index]);
      continue;
    }

    // Swap the edge
    if (!SwapEdge(edge)) continue;

    // Update edge value
    edge_values[edge_index] = -current_value; // R3MeshEdgeSwapValue(this, edge);

    // Update neighbors
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        R3MeshFace *face = FaceOnEdge(edge, i);
        if (!face) continue;
        R3MeshEdge *neighbor_edge = EdgeOnFace(face, edge, j);
        int neighbor_index = EdgeID(neighbor_edge);
        RNScalar old_neighbor_value = edge_values[neighbor_index];
        edge_values[neighbor_index] = R3MeshEdgeSwapValue(this, neighbor_edge);
        if (old_neighbor_value >= min_angle_improvement) {
          if (edge_values[neighbor_index] >= min_angle_improvement) heap.Update(&edge_values[neighbor_index]);
          else heap.Remove(&edge_values[neighbor_index]);
        }
        else {
          if (edge_values[neighbor_index] >= min_angle_improvement) heap.Push(&edge_values[neighbor_index]);
        }
      }
    }
  }

  // Delete edge values
  delete [] edge_values;
}



////////////////////////////////////////////////////////////////////////
// HOLE FILLING FUNCTIONS
////////////////////////////////////////////////////////////////////////

struct R3MeshBoundaryVertexData {
  R3MeshVertex *vertex;
  R3MeshEdge *edge;
  R3MeshBoundaryVertexData *prev;
  R3MeshBoundaryVertexData *next;
  R3MeshBoundaryVertexData **heappointer;
  RNScalar error;
};


#if 0

static void
AssertBoundaryVertexData(R3Mesh *mesh, R3MeshBoundaryVertexData *head)
{
#ifndef NDEBUG
  int counter = 1;
  R3MeshBoundaryVertexData *vdata = head;
  do {
    if (!(vdata->vertex)) RNAbort("%d", counter++);
    if (!(vdata->edge)) RNAbort("%d", counter++);
    if (!(mesh->IsVertexOnEdge(vdata->vertex, vdata->edge))) RNAbort("%d", counter++);
    R3MeshFace *face = mesh->FaceOnVertex(vdata->vertex, vdata->edge, RN_CW);
    if (face) {
      if (!(mesh->IsVertexOnFace(vdata->vertex, face))) RNAbort("%d", counter++);
      if (!(mesh->IsEdgeOnFace(vdata->edge, face))) RNAbort("%d", counter++);
      if ((vdata->edge != mesh->EdgeOnVertex(vdata->vertex, face, RN_CCW))) RNAbort("%d", counter++);
      if (!(mesh->VertexOnFace(face, 0))) RNAbort("%d", counter++);
      if (!(mesh->VertexOnFace(face, 1))) RNAbort("%d", counter++);
      if (!(mesh->VertexOnFace(face, 2))) RNAbort("%d", counter++);
    }
    if (!(vdata->prev)) RNAbort("%d", counter++);
    if ((vdata->prev->next != vdata)) RNAbort("%d", counter++);
    if (!(mesh->IsVertexOnEdge(vdata->vertex, vdata->prev->edge))) RNAbort("%d", counter++);
    R3MeshFace *prev_face = mesh->FaceOnVertex(vdata->prev->vertex, vdata->prev->edge, RN_CW);
    if (prev_face) {
      if (!(mesh->IsVertexOnFace(vdata->vertex, prev_face))) RNAbort("%d", counter++);
    }
    if (!(vdata->next)) RNAbort("%d", counter++);
    if ((vdata->next->prev != vdata)) RNAbort("%d", counter++);
    if (!(mesh->IsVertexOnEdge(vdata->next->vertex, vdata->edge))) RNAbort("%d", counter++);
    R3MeshFace *next_face = mesh->FaceOnVertex(vdata->next->vertex, vdata->next->edge, RN_CW);
    if (next_face) {
      if (!(mesh->IsVertexOnFace(vdata->next->vertex, next_face))) RNAbort("%d", counter++);
    }
    vdata = vdata->next;
  } while (vdata != head);
#endif
}

#endif



static void
FillHoleBoundary(R3Mesh *mesh, R3MeshBoundaryVertexData *head)
{
  // Just checking
  // AssertBoundaryVertexData(mesh, head);

  // Initialize heap
  R3MeshBoundaryVertexData tmp;
  RNHeap<R3MeshBoundaryVertexData *> heap(&tmp, &(tmp.error), &(tmp.heappointer));
  R3MeshBoundaryVertexData *vdata = head;
  if (vdata) do {
    heap.Push(vdata);
    vdata = vdata->next;
  } while (vdata != head);

  // Iteratively select boundary vertex to form apex of new triangle
  while (heap.NEntries() > 3) {
    R3MeshBoundaryVertexData *vdata = heap.Pop();
    R3MeshBoundaryVertexData *prev_vdata = vdata->prev;
    R3MeshBoundaryVertexData *next_vdata = vdata->next;
    // AssertBoundaryVertexData(mesh, vdata);

    // Create face
    R3MeshVertex *vertex = vdata->vertex;
    R3MeshVertex *prev_vertex = vdata->prev->vertex;
    R3MeshVertex *next_vertex = vdata->next->vertex;
    R3MeshFace *face = mesh->CreateFace(prev_vertex, vertex, next_vertex);
    if (!face) { /* ERROR */ return; }

    // Update edge 
    vdata->prev->edge = mesh->EdgeAcrossFace(face, vertex);

    // Remove vdata from linked list
    vdata->prev->next = vdata->next;
    vdata->next->prev = vdata->prev;
    delete vdata;

    // Update error of prev_vdata and next_vdata
    for (int i = 0; i < 2; i++) {
      R3MeshBoundaryVertexData *vdata = (i == 0) ? prev_vdata : next_vdata;
      if (!vdata) continue;
      R3MeshVertex *prev_vertex = vdata->prev->vertex;
      R3MeshVertex *next_vertex = vdata->next->vertex;
      const R3Point& position = mesh->VertexPosition(vertex);
      const R3Point& prev_position = mesh->VertexPosition(prev_vertex);
      const R3Point& next_position = mesh->VertexPosition(next_vertex);
      const R3Vector& prev_normal = mesh->VertexNormal(prev_vertex);
      const R3Vector& next_normal = mesh->VertexNormal(next_vertex);
      R3Vector prev_vector = prev_position - position;
      R3Vector next_vector = next_position - position;
      R3Vector cross = prev_vector % next_vector;
      RNAngle interior_angle = R3InteriorAngle(prev_vector, next_vector);
      if (cross.Dot(mesh->VertexNormal(vdata->vertex)) > 0) interior_angle = RN_TWO_PI - interior_angle;
      RNAngle normal_angle = R3InteriorAngle(prev_normal, next_normal);
      vdata->error = interior_angle + normal_angle;
      heap.Update(vdata);
    }
  }

  // Create last triangle
  if (heap.NEntries() == 3) {
    R3MeshBoundaryVertexData *vdata = heap.Pop();
    R3MeshVertex *vertex = vdata->vertex;
    R3MeshVertex *prev_vertex = vdata->prev->vertex;
    R3MeshVertex *next_vertex = vdata->next->vertex;
    mesh->CreateFace(prev_vertex, vertex, next_vertex);
    delete vdata;
    vdata = heap.Pop(); delete vdata;
    vdata = heap.Pop(); delete vdata;
  }
}



void R3Mesh::
FillHole(R3MeshEdge *seed_edge)
{
  // Check if seed edge is on boundary
  R3MeshFace *seed_face = FaceOnEdge(seed_edge);
  if (!seed_face) return; 
  if (FaceAcrossEdge(seed_edge, seed_face)) return;
  R3MeshVertex *seed_vertex = VertexOnEdge(seed_edge, seed_face, RN_CCW);

  // Create counterclockwise linked list of boundary vertices
  R3mesh_mark++;
  R3MeshBoundaryVertexData *head = NULL;
  R3MeshBoundaryVertexData *tail = NULL;
  R3MeshEdge *edge = seed_edge;
  R3MeshVertex *vertex = seed_vertex;
  // R3MeshFace *face = seed_face;
  do {
    // Create boundary vertex 
    R3MeshBoundaryVertexData *vdata = new R3MeshBoundaryVertexData();
    vdata->vertex = vertex;
    vdata->edge = edge;
    vdata->prev = NULL;
    vdata->next = NULL;
    vdata->heappointer = NULL;
    vdata->error = 0;

    // Mark vertex (so can detect figure 8)
    SetVertexMark(vertex, R3mesh_mark);

    // Add boundary vertex to list
    if (!head) head = vdata;
    if (tail) tail->next = vdata;
    vdata->prev = tail;
    tail = vdata;

    // Update vdata error 
    if (vdata->prev) {
      R3MeshVertex *prev_vertex = vdata->prev->vertex;
      R3MeshVertex *next_vertex = VertexAcrossEdge(vdata->edge, vdata->vertex);
      const R3Point& position = VertexPosition(vertex);
      const R3Point& prev_position = VertexPosition(prev_vertex);
      const R3Point& next_position = VertexPosition(next_vertex);
      const R3Vector& prev_normal = VertexNormal(prev_vertex);
      const R3Vector& next_normal = VertexNormal(next_vertex);
      R3Vector prev_vector = prev_position - position;
      R3Vector next_vector = next_position - position;
      R3Vector cross = prev_vector % next_vector;
      RNAngle interior_angle = R3InteriorAngle(prev_vector, next_vector);
      if (cross.Dot(VertexNormal(vdata->vertex)) > 0) interior_angle = RN_TWO_PI - interior_angle;
      RNAngle normal_angle = R3InteriorAngle(prev_normal, next_normal);
      vdata->error = interior_angle + normal_angle;
    }

    // Find next vertex
    vertex = VertexAcrossEdge(edge, vertex);

    // Check for repeated vertex (figure 8)
    if (VertexMark(vertex) == R3mesh_mark) {
      if (vertex != seed_vertex) {
        // Find previous instance of repeated vertex
        R3MeshBoundaryVertexData *vdata = tail->prev;
        while (vdata->vertex != vertex) vdata = vdata->prev;

        // Remove part of list before repeated vertex
        head = vdata;
        vdata = head->prev;
        head->prev = NULL;
        while (vdata) { 
          R3MeshBoundaryVertexData *prev = vdata->prev; 
          delete vdata; 
          vdata = prev; 
        }

        // Now have a complete loop
        break;
      }
    }

    // Find next edge
    R3MeshEdge *e = edge;
    while (e) {
      R3MeshFace *f = FaceOnVertex(vertex, e, RN_CCW);
      if (!f) { assert(e != edge); edge = e; break; }
      e = EdgeOnFace(f, vertex, RN_CW);
      if (!e || (e == edge)) return; 
    } 

    // Find next face
    // face = FaceOnVertex(vertex, edge, RN_CW);
  } while (vertex != seed_vertex);

  // Complete loop
  if (head && tail) {
    // Add link between head and tail
    head->prev = tail;
    tail->next = head;

    // Update head error 
    if (head) {
      R3MeshBoundaryVertexData *vdata = head;
      R3MeshVertex *prev_vertex = vdata->prev->vertex;
      R3MeshVertex *next_vertex = vdata->next->vertex;
      const R3Point& position = VertexPosition(vertex);
      const R3Point& prev_position = VertexPosition(prev_vertex);
      const R3Point& next_position = VertexPosition(next_vertex);
      const R3Vector& prev_normal = VertexNormal(prev_vertex);
      const R3Vector& next_normal = VertexNormal(next_vertex);
      R3Vector prev_vector = prev_position - position;
      R3Vector next_vector = next_position - position;
      R3Vector cross = prev_vector % next_vector;
      RNAngle interior_angle = R3InteriorAngle(prev_vector, next_vector);
      if (cross.Dot(VertexNormal(vdata->vertex)) > 0) interior_angle = RN_TWO_PI - interior_angle;
      RNAngle normal_angle = R3InteriorAngle(prev_normal, next_normal);
      vdata->error = interior_angle + normal_angle;
    }
  }

  // Split any edges spanning (back side of) boundary
  R3MeshBoundaryVertexData *vdata = head;
  do {
    R3MeshVertex *vertex = vdata->vertex;
    for (int i = 0; i < VertexValence(vertex); i++) {
      R3MeshEdge *edge = EdgeOnVertex(vertex, i);
      if (edge == vdata->edge) continue;
      if (edge == vdata->prev->edge) continue;
      R3MeshVertex *neighbor = VertexAcrossEdge(edge, vertex);
      if (VertexMark(neighbor)) SubdivideEdge(edge);
    }
    vdata = vdata->next;
  } while (vdata != head);

  // Fill hole
  FillHoleBoundary(this, head);
}



void R3Mesh::
FillHoles(void)
{
  // Iteratively find boundary edge and fill hole
  for (int i = 0; i < NEdges(); i++) {
    R3MeshEdge *edge = Edge(i);
    if (!IsEdgeOnBoundary(edge)) continue;
    FillHole(edge);
  }
}




////////////////////////////////////////////////////////////////////////
// DRAW FUNCTIONS
////////////////////////////////////////////////////////////////////////

void R3Mesh::
DrawVertices(void) const
{
  // Draw all vertices
  for (int i = 0; i < vertices.NEntries(); i++)
    DrawVertex(vertices[i]);
}



void R3Mesh::
DrawEdges(void) const
{
  // Draw all edges
  for (int i = 0; i < edges.NEntries(); i++)
    DrawEdge(edges[i]);
}



void R3Mesh::
DrawFaces(void) const
{
  // Draw all faces
  for (int i = 0; i < faces.NEntries(); i++) 
    DrawFace(faces[i]);
}



void R3Mesh::
DrawVertexIDs(void) const
{
  // Draw all vertex IDs
  glDisable(GL_LIGHTING);
  for (int i = 0; i < vertices.NEntries(); i++) {
    assert(vertices[i]->id == i);
    unsigned char r = ((i << 16) & 0xFF);
    unsigned char g = (i << 8) & 0xFF;
    unsigned char b = (i << 0) & 0xFF;
    glColor3ub(r, g, b);
    DrawVertex(vertices[i]);
  }
  glEnable(GL_LIGHTING);
}



void R3Mesh::
DrawEdgeIDs(void) const
{
  // Draw all edge IDs
  glDisable(GL_LIGHTING);
  for (int i = 0; i < edges.NEntries(); i++) {
    assert(edges[i]->id == i);
    unsigned char r = ((i << 16) & 0xFF);
    unsigned char g = (i << 8) & 0xFF;
    unsigned char b = (i << 0) & 0xFF;
    glColor3ub(r, g, b);
    DrawEdge(edges[i]);
  }
  glEnable(GL_LIGHTING);
}



void R3Mesh::
DrawFaceIDs(void) const
{
  // Draw all face IDs
  glDisable(GL_LIGHTING);
  for (int i = 0; i < faces.NEntries(); i++) {
    assert(faces[i]->id == i);
    unsigned char r = (i << 16) & 0xFF;
    unsigned char g = (i << 8) & 0xFF;
    unsigned char b = (i << 0) & 0xFF;
    glColor3ub(r, g, b);
    DrawFace(faces[i]);
  }
  glEnable(GL_LIGHTING);
}



void R3Mesh::
DrawVertex(R3MeshVertex *v) const
{
  // Draw box around vertex 
  RNScalar d = 0.001 * BBox().LongestAxisLength();
  R3Sphere(v->position, d).Draw();
}



void R3Mesh::
DrawEdge(R3MeshEdge *e) const
{
  // Draw edge
  R3BeginLine();
  R3LoadPoint(e->vertex[0]->position);
  R3LoadPoint(e->vertex[1]->position);
  R3EndLine();
}



void R3Mesh::
DrawFace(R3MeshFace *f) const
{
  // Draw polygon
  R3BeginPolygon();
  R3LoadNormal(FaceNormal(f));
  R3LoadPoint(f->vertex[0]->position);
  R3LoadPoint(f->vertex[1]->position);
  R3LoadPoint(f->vertex[2]->position);
  R3EndPolygon();
}



////////////////////////////////////////////////////////////////////////
// INTERSECTION FUNCTIONS
////////////////////////////////////////////////////////////////////////

R3MeshType R3Mesh::
Intersection(const R3Ray& ray, R3MeshIntersection *intersection) const
{
  // Initialize intersection variables
  R3MeshType type = R3_MESH_NULL_TYPE;
  if (intersection) {
    intersection->type = R3_MESH_NULL_TYPE;
    intersection->vertex = NULL;
    intersection->edge = NULL;
    intersection->face = NULL;
    intersection->t = RN_INFINITY;
  }

  // Check bounding box for intersection 
  if (R3Intersects(ray, bbox)) {
    // Check each face to find closest intersection
    RNScalar min_t = FLT_MAX;
    for (int i = 0; i < faces.NEntries(); i++) {
      // Get ith face
      R3MeshFace *f = faces[i];

      // Get intersection with face
      R3MeshIntersection face_intersection;
      if (Intersection(ray, f, &face_intersection)) {
        if (face_intersection.t < min_t) {
          if (intersection) *intersection = face_intersection;
          type = face_intersection.type;
          min_t = face_intersection.t;
        }
      }
    }
  }

  // Return intersection type
  return type;
}



R3MeshType R3Mesh::
Intersection(const R3Ray& ray, R3MeshFace *f, R3MeshIntersection *intersection) const
{
  // Initialize pick variables
  R3MeshType type = R3_MESH_NULL_TYPE;
  if (intersection) {
    intersection->type = R3_MESH_NULL_TYPE;
    intersection->vertex = NULL;
    intersection->edge = NULL;
    intersection->face = NULL;
    intersection->t = RN_INFINITY;
  }

  // Check if face intersects ray
  RNScalar t;
  R3Point p;
  R3Plane plane = FacePlane(f);
  if (R3Intersects(ray, plane, &p, &t) || R3Intersects(ray, -plane, &p, &t)) {
    // Check if face bbox contains p
    if (R3Contains(FaceBBox(f), p)) {
#if 1
      // Get convenient variables
      R3Point& p0 = f->vertex[0]->position;
      R3Point& p1 = f->vertex[1]->position;
      R3Point& p2 = f->vertex[2]->position;

      // Check side of first edge
      R3Vector e0 = p1 - p0;
      e0.Normalize();
      R3Vector n0 = FaceNormal(f) % e0;
      R3Plane s0(p0, n0);
      RNScalar b0 = R3SignedDistance(s0, p);
      if (RNIsNegative(b0)) return R3_MESH_NULL_TYPE;

      // Check side of second edge
      R3Vector e1 = p2 - p1;
      e1.Normalize();
      R3Vector n1 = FaceNormal(f) % e1;
      R3Plane s1(p1, n1);
      RNScalar b1 = R3SignedDistance(s1, p);
      if (RNIsNegative(b1)) return R3_MESH_NULL_TYPE;

      // Check side of third edge
      R3Vector e2 = p0 - p2;
      e2.Normalize();
      R3Vector n2 = FaceNormal(f) % e2;
      R3Plane s2(p2, n2);
      RNScalar b2 = R3SignedDistance(s2, p);
      if (RNIsNegative(b2)) return R3_MESH_NULL_TYPE;

      // Fill in interesection info
      // Could easily return intersection type, etc. but too lazy now
      type = R3_MESH_FACE_TYPE;
      if (intersection) {
        intersection->type = type;
        intersection->point = p;
        intersection->vertex = NULL;
        intersection->edge = NULL;
        intersection->face = f;
        intersection->t = t;
      }
#else
      // This does not always work -- don't know why
      // Compute interpolation parameters for intersection point (from Graphics Gems I, page 393)
      RNScalar a, b;
      RNDimension dim = FaceNormal(f).MaxDimension();
      RNDimension dim1 = (dim + 1) % 3;
      RNDimension dim2 = (dim + 2) % 3;
      R3Point& p0 = f->vertex[0]->position;
      R3Point& p1 = f->vertex[1]->position;
      R3Point& p2 = f->vertex[2]->position;
      RNScalar u0 = p[dim1] - p0[dim1];
      RNScalar v0 = p[dim2] - p0[dim2];
      RNScalar u1 = p1[dim1] - p0[dim1];
      RNScalar u2 = p2[dim1] - p0[dim1];
      RNScalar v1 = p1[dim2] - p0[dim2];
      RNScalar v2 = p2[dim2] - p0[dim2];
      if (RNIsZero(u1)) {
        if (RNIsZero(u2)) {
          a = RN_INFINITY;
          b = RN_INFINITY;
        }
        else {
          b = u0/u2;
          if (RNIsLess(b, 0.0) || RNIsGreater(b, 1.0) || RNIsZero(v1)) a = b = RN_INFINITY;
          else a = (v0 - b*v2) / v1;
        }
      }
      else {
        RNScalar denom = v2*u1 - u2*v1;
        if (RNIsZero(denom)) a = b = RN_INFINITY;
        else {
          b = (v0*u1 - u0*v1) / denom;
          if (RNIsLess(b, 0.0) || RNIsGreater(b, 1.0) || RNIsZero(u1)) a = b = RN_INFINITY;
          else a = (u0 - b*u2) / u1;
        }
      }

      // Check if intersection point is inside triangle
      if (RNIsPositiveOrZero(a, 0.001) && RNIsPositiveOrZero(b, 0.001) && RNIsLess(a+b, 1.0, 0.001)) {
        if (RNIsZero(a, 0.01)) {
          if (RNIsZero(b, 0.01)) {
            type = R3_MESH_VERTEX_TYPE;
            if (intersection) {
              intersection->type = type;
              intersection->point = p;
              intersection->vertex = f->vertex[0];
              intersection->edge = f->edge[0];
              intersection->face = f;
              intersection->t = t;
            }
          }
          else if (RNIsEqual(b, 1.0, 0.01)) {
            type = R3_MESH_VERTEX_TYPE;
            if (intersection) {
              intersection->type = type;
              intersection->point = p;
              intersection->vertex = f->vertex[2];
              intersection->edge = f->edge[2];
              intersection->face = f;
              intersection->t = t;
            }
          }
          else {
            type = R3_MESH_EDGE_TYPE;
            if (intersection) {
              intersection->type = type;
              intersection->point = p;
              intersection->vertex = NULL;
              intersection->edge = f->edge[2];
              intersection->face = f;
              intersection->t = t;
            }
          }
        }
        else if (RNIsEqual(a, 1.0, 0.01)) {
          type = R3_MESH_VERTEX_TYPE;
          if (intersection) {
            intersection->type = type;
            intersection->point = p;
            intersection->vertex = f->vertex[1];
            intersection->edge = f->edge[1];
            intersection->face = f;
            intersection->t = t;
          }
        }                      
        else {
          if (RNIsZero(b, 0.01)) {
            type = R3_MESH_EDGE_TYPE;
            if (intersection) {
              intersection->type = type;
              intersection->point = p;
              intersection->vertex = NULL;
              intersection->edge = f->edge[0];
              intersection->face = f;
              intersection->t = t;
            }
          }
          else if (RNIsEqual(a+b, 1.0, 0.01)) {
            type = R3_MESH_EDGE_TYPE;
            if (intersection) {
              intersection->type = type;
              intersection->point = p;
              intersection->vertex = NULL;
              intersection->edge = f->edge[1];
              intersection->face = f;
              intersection->t = t;
            }
          }
          else {
            type = R3_MESH_FACE_TYPE;
            if (intersection) {
              intersection->type = type;
              intersection->point = p;
              intersection->vertex = NULL;
              intersection->edge = NULL;
              intersection->face = f;
              intersection->t = t;
            }
          }
        }
      }
#endif
    }
  }

  // Return intersection type
  return type;
}



////////////////////////////////////////////////////////////////////////
// CLOSEST POINT FUNCTIONS
////////////////////////////////////////////////////////////////////////

R3Point R3Mesh::
ClosestPointOnEdge(const R3MeshEdge *e, const R3Point& point, R3MeshIntersection *closest_point) const
{
  // Compute parametric value of closest point on edge
  R3MeshVertex *v0 = VertexOnEdge(e, 0);
  R3MeshVertex *v1 = VertexOnEdge(e, 1);
  const R3Point& p0 = VertexPosition(v0);
  const R3Point& p1 = VertexPosition(v1);
  R3Vector edge_vector = p1 - p0;
  RNScalar edge_length = edge_vector.Length();
  if (edge_length == 0) return p0;
  edge_vector /= edge_length;
  R3Vector point_vector = point - p0;
  RNScalar t = edge_vector.Dot(point_vector);

  // Compute closest point position
  R3Point p;
  if (t <= 0) p = p0;
  else if (t >= edge_length) p = p1;
  else p = p0 + t * edge_vector;
  
  // Fill in closest point info
  if (closest_point) {
    closest_point->edge = (R3MeshEdge *) e;
    closest_point->face = NULL;
    closest_point->point = p;
    closest_point->t = R3Distance(p, point);
    if (t <= 0) {
      closest_point->type = R3_MESH_VERTEX_TYPE;
      closest_point->vertex = v0;
    }
    else if (t >= edge_length) {
      closest_point->type = R3_MESH_VERTEX_TYPE;
      closest_point->vertex = v1;
    }
    else {
      closest_point->type = R3_MESH_EDGE_TYPE;
      closest_point->vertex = NULL;
    }
  }

  // Return closest point
  return p;
}



R3Point R3Mesh::
ClosestPointOnFace(const R3MeshFace *face, const R3Point& point, R3MeshIntersection *closest_point) const
{
  // Project point point onto face plane
  const R3Plane& plane = FacePlane(face);
  const R3Vector& face_normal = plane.Normal();
  RNScalar plane_signed_distance = R3SignedDistance(plane, point);
  R3Point plane_point = point - plane_signed_distance * face_normal;

  // Check if point is outside each edge
  R3MeshEdge *closest_edge = NULL;
  R3Point closest_edge_point = R3zero_point;
  RNScalar closest_edge_distance_squared = FLT_MAX;
  for (int i = 0; i < 3; i++) {
    R3MeshEdge *edge = EdgeOnFace(face, i);
    R3MeshVertex *v0 = VertexOnEdge(edge, face, RN_CW);
    R3MeshVertex *v1 = VertexOnEdge(edge, face, RN_CCW);
    R3Point p0 = v0->position;
    R3Point p1 = v1->position;
    R3Vector edge_vector = p1 - p0;
    edge_vector.Normalize();
    R3Vector edge_normal = face_normal % edge_vector;
    R3Plane edge_plane(p0, edge_normal);
    RNScalar b = R3SignedDistance(edge_plane, plane_point);
    if (b < 0) {
      R3Point edge_point = ClosestPointOnEdge(edge, point);
      RNScalar distance_squared = R3SquaredDistance(edge_point, point);
      if (distance_squared < closest_edge_distance_squared) {
        closest_edge_distance_squared = distance_squared;
        closest_edge_point = edge_point;
        closest_edge = edge;
      }
    }
  }

  // Compute closest position
  R3Point p = (closest_edge) ? closest_edge_point : plane_point;

  // Fill in closest point info
  if (closest_point) {
    closest_point->face = (R3MeshFace *) face;
    closest_point->point = p;
    closest_point->t = R3Distance(p, point);
    if (closest_edge) {
      closest_point->edge = closest_edge;
      R3Span edge_span = EdgeSpan(closest_edge);
      RNScalar t = edge_span.T(point);
      if (t <= 0) {
        closest_point->type = R3_MESH_VERTEX_TYPE;
        closest_point->vertex = VertexOnEdge(closest_edge, 0);
      }
      else if (t >= edge_span.Length()) {
        closest_point->type = R3_MESH_VERTEX_TYPE;
        closest_point->vertex = VertexOnEdge(closest_edge, 1);
      }
      else {
        closest_point->type = R3_MESH_EDGE_TYPE;
        closest_point->vertex = NULL;
      }
    }
    else {
      closest_point->type = R3_MESH_FACE_TYPE;
      closest_point->vertex = NULL;
      closest_point->edge = NULL;
    }
  }

  // Return closest point
  return p;
}



R3Point R3Mesh::
ClosestPoint(const R3Point& point, R3MeshIntersection *closest_point) const
{
  // Check each face to find closest point
  R3Point closest_position = R3zero_point;
  RNScalar closest_distance_squared = FLT_MAX;
  for (int i = 0; i < faces.NEntries(); i++) {
    R3MeshFace *f = faces[i];
    R3MeshIntersection intersection;
    R3Point position = ClosestPointOnFace(f, point, (closest_point) ? &intersection : NULL);
    RNLength distance_squared = R3SquaredDistance(position, point);
    if (distance_squared < closest_distance_squared) {
      if (closest_point) *closest_point = intersection;
      closest_distance_squared = distance_squared;
      closest_position = position;
    }
  }

  // Return closest position on mesh
  return closest_position;
}



////////////////////////////////////////////////////////////////////////
// POINT SAMPLING STUFF
////////////////////////////////////////////////////////////////////////

R3Point R3Mesh::
RandomPointOnFace(const R3MeshFace *face) const
{
  // Seed random number generator
  static RNBoolean seed = 0;
  if (!seed) { seed = 1; RNSeedRandomScalar(); }

  // Get vertex positions
  R3MeshVertex *v0 = VertexOnFace(face, 0);
  R3MeshVertex *v1 = VertexOnFace(face, 1);
  R3MeshVertex *v2 = VertexOnFace(face, 2);
  const R3Point& p0 = VertexPosition(v0);
  const R3Point& p1 = VertexPosition(v1);
  const R3Point& p2 = VertexPosition(v2);


  // Return random point on face
  RNScalar r1 = sqrt(RNRandomScalar());
  RNScalar r2 = RNRandomScalar();
  R3Point p = p0 * (1.0 - r1) + p1 * r1 * (1.0 - r2) + p2 * r1 * r2;
  return p;
}



////////////////////////////////////////////////////////////////////////
// EXACT GEODESIC DISTANCE STUFF
////////////////////////////////////////////////////////////////////////

struct R3MeshGeodesicBeam {
  const R3MeshVertex *source; // vertex from which beam originated
  const R3MeshFace *face; // face from which beam arrived at edge
  const R3MeshEdge *edge; // edge reached by beam
  R3Affine affine; // transformation from 3D to 2D beam space
  R3Span span; // cross-section of beam in 2D beam space
  RNLength distance; // lower bound on total distance to span
  RNLength distance_to_source; // may be non-zero if beam was spawned from saddle vertex
  struct R3MeshGeodesicBeam *parent;
};



static void
R3MeshCreateGeodesicBeamsForVertex(const R3Mesh *mesh, const R3MeshVertex *source, 
  RNLength max_distance, RNLength distance_to_source, R3MeshGeodesicBeam *parent,
  RNLength *vertex_distances, RNLength *edge_distances, RNLength *face_distances,
  RNHeap<R3MeshGeodesicBeam *>& heap, RNArray<R3MeshGeodesicBeam *>& beams)
{
  // Create beams through edges across faces adjacent to source
  for (int i = 0; i < mesh->VertexValence(source); i++) {
    R3MeshEdge *e = mesh->EdgeOnVertex(source, i);
    R3MeshFace *face = mesh->FaceOnEdge(e, source, RN_CCW);
    if (!face) continue;

    // Get edge across face from  source
    R3MeshEdge *edge = mesh->EdgeAcrossFace(face, source);
    if (!mesh->FaceAcrossEdge(edge, face)) continue;

    // Construct beam coordinate system
    R3MeshVertex *v_ccw = mesh->VertexOnFace(face, edge, RN_CCW);
    R3MeshVertex *v_cw = mesh->VertexOnFace(face, edge, RN_CW);
    R3Vector face_normal = mesh->FaceNormal(face);
    R3Span edge_span(mesh->VertexPosition(v_ccw), mesh->VertexPosition(v_cw));
    R3Vector edge_vector = edge_span.Vector();
    R3Vector edge_normal = face_normal % edge_vector; edge_normal.Normalize();
    R3Triad edge_triad(edge_vector, edge_normal, face_normal);
    R3Point source_position = mesh->VertexPosition(source);
    R3CoordSystem face_cs(source_position, edge_triad);
    R3Affine affine(face_cs.InverseMatrix(), 0);

    // Get edge span in beam coordinate system
    edge_span.Transform(affine);

    // Get lower bound on distance to edge along beam
    RNLength beam_distance = distance_to_source + R3Distance(R3zero_point, edge_span);
    if ((max_distance > 0) && (beam_distance > max_distance)) continue;

    // Check distance to edge
    if (beam_distance >= edge_distances[mesh->EdgeID(edge)]) continue;

    // Check distance to neighbor face
    R3MeshFace *neighbor_face = mesh->FaceAcrossEdge(edge, face);
    if ((neighbor_face) && (beam_distance >= face_distances[mesh->FaceID(neighbor_face)])) continue;

    // Create beam
    R3MeshGeodesicBeam *beam = new R3MeshGeodesicBeam();
    beam->source = source;
    beam->face = face;
    beam->edge = edge;
    beam->affine = affine;
    beam->span = edge_span;
    beam->distance = beam_distance;
    beam->distance_to_source = distance_to_source;
    beam->parent = parent;

    // Add beam to heap
    heap.Push(beam);
    beams.Insert(beam);
  }
}
 


static RNBoolean
R3MeshVertexIsSaddleOrBoundary(const R3Mesh *mesh, const R3MeshVertex *vertex)
{
  // Return whether vertex is a saddle or boundary
  int nconvex = 0;
  int nconcave = 0;
  RNScalar angle_epsilon = 0.1;
  for (int i = 0; i < mesh->VertexValence(vertex); i++) {
    R3MeshEdge *e = mesh->EdgeOnVertex(vertex, i);
    if (!mesh->FaceOnEdge(e, 0)) return TRUE;
    if (!mesh->FaceOnEdge(e, 1)) return TRUE;
    RNAngle angle = mesh->EdgeInteriorAngle(e);
    if (RNIsGreater(angle, RN_PI, angle_epsilon)) nconcave++;
    else if (RNIsLess(angle, RN_PI, angle_epsilon)) nconvex++;
    if (nconvex * nconcave != 0) return TRUE;
  }

  // Not a saddle
  return FALSE;
}



RNLength *R3Mesh::
GeodesicDistances(const R3MeshVertex *source_vertex, RNLength max_distance) const
{
  // Create arrays of upper bounds on distances
  RNLength *vertex_distances = new RNLength [ NVertices() ];
  RNLength *edge_distances = new RNLength [ NEdges() ];
  RNLength *face_distances = new RNLength [ NFaces() ];
  assert(vertex_distances && edge_distances && face_distances);
  for (int i = 0; i < NVertices(); i++) vertex_distances[i] = FLT_MAX;
  for (int i = 0; i < NEdges(); i++) edge_distances[i] = FLT_MAX;
  for (int i = 0; i < NFaces(); i++) face_distances[i] = FLT_MAX;
  vertex_distances [VertexID(source_vertex)] = 0;

  // Initialize heap of beams originating at source vertex
  R3MeshGeodesicBeam tmp;
  RNHeap<R3MeshGeodesicBeam *> heap(&tmp, &(tmp.distance));
  RNArray<R3MeshGeodesicBeam *> beams;
  R3MeshCreateGeodesicBeamsForVertex(this, source_vertex, max_distance, 0, NULL, 
    vertex_distances, edge_distances, face_distances, heap, beams);
  for (int i = 0; i < VertexValence(source_vertex); i++) {
    R3MeshEdge *e = EdgeOnVertex(source_vertex, i);
    R3MeshVertex *vertex = VertexAcrossEdge(e, source_vertex);
    edge_distances [EdgeID(e)] = EdgeLength(e);
    vertex_distances [VertexID(vertex)] = EdgeLength(e);
    if (R3MeshVertexIsSaddleOrBoundary(this, vertex)) {
      R3MeshCreateGeodesicBeamsForVertex(this, vertex, max_distance, EdgeLength(e), NULL, 
        vertex_distances, edge_distances, face_distances, heap, beams);
    }
  }

  // Pop beams off heap
  while (!heap.IsEmpty()) {
    // Get beam
    R3MeshGeodesicBeam *beam = heap.Pop();
    assert(beam);

    // Get face across edge
    R3MeshFace *face = FaceAcrossEdge(beam->edge, beam->face);
    assert(face);

    // Check beam distances versus upper bounds on distances found already
    if (beam->distance > beam->distance - beam->distance_to_source + vertex_distances[VertexID(beam->source)]) continue;
    if (beam->distance >= edge_distances[EdgeID(beam->edge)]) continue;
    if (beam->distance >= face_distances[FaceID(face)]) continue;

    // Construct beam extent
    R3Point left_span_position = beam->span.Start();
    R3Vector left_vector =  left_span_position.Vector();
    R3Vector left_normal = left_vector % R3posz_vector; left_normal.Normalize();
    R3Plane left_plane(left_span_position, left_normal);
    R3Point right_span_position = beam->span.End();
    R3Vector right_vector =  right_span_position.Vector();
    R3Vector right_normal = R3posz_vector % right_vector; right_normal.Normalize();
    R3Plane right_plane(right_span_position, right_normal);

    // Find affine that takes points on face into beam coordinate system
    R3Vector face_normal = FaceNormal(face);
    face_normal.Transform(beam->affine);
    R3Affine affine = R3identity_affine;
    affine.Translate(beam->span.Midpoint().Vector());
    affine.Rotate(face_normal, R3posz_vector);
    affine.Translate(-(beam->span.Midpoint().Vector()));
    affine.Transform(beam->affine);

    // Spawn beams across edges
    for (int i = 0; i < 2; i++) {
      R3MeshEdge *edge = EdgeOnFace(face, beam->edge, i);
      R3MeshFace *neighbor_face = FaceAcrossEdge(edge, face);
      if (!neighbor_face) continue;

      // Check versus upper-bound distances
      if (beam->distance >= edge_distances[EdgeID(edge)]) continue;
      if (beam->distance >= face_distances[FaceID(neighbor_face)]) continue;

      // Get edge span in unfolded beam coordinate system
      R3MeshVertex *v_ccw = VertexOnFace(face, edge, RN_CCW);
      R3MeshVertex *v_cw = VertexOnFace(face, edge, RN_CW);
      R3Point p_ccw = VertexPosition(v_ccw);
      R3Point p_cw = VertexPosition(v_cw);
      p_ccw.Transform(affine);
      p_cw.Transform(affine);
      R3Span edge_span(p_ccw, p_cw);
      // assert(RNIsZero(p_ccw[2]));
      // assert(RNIsZero(p_cw[2]));

      // Check if source is on correct side of edge span
      R3Vector span_vector = edge_span.Vector();
      if (RNIsZero(span_vector.Length())) continue;
      R3Vector span_normal = R3posz_vector % span_vector; span_normal.Normalize();
      R3Plane span_plane(edge_span.Midpoint(), span_normal);
      if (RNIsPositive(R3SignedDistance(span_plane, R3zero_point))) continue;
      
      // Check for edge-beam intersection
      // RNScalar left_dcw = R3SignedDistance(left_plane, edge_span.End());
      // if (RNIsNegative(left_dcw)) continue;
      // RNScalar right_dccw = R3SignedDistance(right_plane, edge_span.Start());
      // if (RNIsNegative(right_dccw)) continue;
        
      // Clip edge span to beam
      if (!edge_span.Clip(left_plane)) continue;
      if (!edge_span.Clip(right_plane)) continue;

      // Get lower bound on distance to edge span
      RNLength beam_distance = beam->distance_to_source + R3Distance(R3zero_point, edge_span);
      if ((max_distance > 0) && (beam_distance > max_distance)) continue;

      // Check versus upper-bound distances
      if (beam_distance >= edge_distances[EdgeID(edge)]) continue;
      if (beam_distance >= face_distances[FaceID(neighbor_face)]) continue;

      // Create beam
      R3MeshGeodesicBeam *child_beam = new R3MeshGeodesicBeam();
      child_beam->source = beam->source;
      child_beam->face = face;
      child_beam->edge = edge;
      child_beam->affine = affine;
      child_beam->span = edge_span;
      child_beam->distance = beam_distance;
      child_beam->distance_to_source = beam->distance_to_source;
      child_beam->parent = beam;

      // Add beam to heap
      heap.Push(child_beam);
      beams.Insert(child_beam);
    }

    // Update upper bound distances and spawn beams from vertices
    for (int i = 0; i < 3; i++) {
      R3MeshVertex *vertex = VertexOnFace(face, i);
      R3Point vertex_position = VertexPosition(vertex);
      vertex_position.Transform(affine);

      // Compute distance to vertex along this path
      RNLength vertex_distance = 0;
      RNBoolean vertex_visible = FALSE;
      if (R3SignedDistance(left_plane, vertex_position) < 0) {
        // Vertex is to the left of the beam
        RNLength left_distance1 = R3Distance(R3zero_point, left_span_position);
        RNLength left_distance2 = R3Distance(left_span_position, vertex_position);
        vertex_distance = beam->distance_to_source + left_distance1 + left_distance2;
      }
      else if (R3SignedDistance(right_plane, vertex_position) < 0) {
        // Vertex is to the right of the beam
        RNLength right_distance1 = R3Distance(R3zero_point, right_span_position);
        RNLength right_distance2 = R3Distance(right_span_position, vertex_position);
        vertex_distance = beam->distance_to_source + right_distance1 + right_distance2;
      }
      else {
        // Vertex is inside the beam
        vertex_distance = beam->distance_to_source + R3Distance(R3zero_point, vertex_position);
        vertex_visible = TRUE;
      }

      // Update upper bound distances
      if (vertex_distance < vertex_distances[VertexID(vertex)]) {
        // Remember vertex distance
        RNLength old_vertex_distance = vertex_distances[VertexID(vertex)];
        vertex_distances[VertexID(vertex)] = vertex_distance;

        // Update edge distances
        for (int j = 0; j < 2; j++) {
          R3MeshEdge *edge = EdgeOnFace(face, beam->edge, j);
          if (vertex_distance < edge_distances[EdgeID(edge)]) {
            R3MeshVertex *v = VertexAcrossEdge(edge, vertex);
            if (vertex_distance > vertex_distances[VertexID(v)]) {
              edge_distances[EdgeID(edge)] = vertex_distance;
            }
          }
        }

        // Update face distance
        if (vertex_distance < face_distances[FaceID(face)]) {
          R3MeshVertex *v0 = VertexOnFace(face, vertex, RN_CCW);
          R3MeshVertex *v1 = VertexOnFace(face, vertex, RN_CW);
          if ((vertex_distance > vertex_distances[VertexID(v0)]) && (vertex_distance > vertex_distances[VertexID(v1)])) {
            face_distances[FaceID(face)] = vertex_distance;
          }
        }

        // Spawn beam from "across" vertex if it is a visible saddle or boundary
        if (vertex_visible && (vertex == VertexAcrossFace(face, beam->edge))) { 
          if (RNIsLess(vertex_distance, old_vertex_distance)) {
            if (R3MeshVertexIsSaddleOrBoundary(this, vertex)) {
              R3MeshCreateGeodesicBeamsForVertex(this, vertex, max_distance, vertex_distance, beam, 
                vertex_distances, edge_distances, face_distances, heap, beams);
            }
          }
        }
      }
    }
  }

  // Delete all beams
  for (int i = 0; i < beams.NEntries(); i++) {
    delete beams[i];
  }

  // Delete edge and face distances
  delete [] edge_distances;
  delete [] face_distances;

  // Return vertex distances
  return vertex_distances;
}



////////////////////////////////////////////////////////////////////////
// DIJKSTRA DISTANCE STUFF
////////////////////////////////////////////////////////////////////////

struct DijkstraData {
  R3MeshVertex *vertex;
  R3MeshEdge *edge;
  double distance_to_destination;
  DijkstraData **heappointer;
  double distance;
};



RNScalar DijkstraDataValue(DijkstraData *data, void *)
{
  return data->distance;
}



RNLength R3Mesh::
DijkstraDistance(const R3MeshVertex *source_vertex, const R3MeshVertex *destination_vertex, RNArray<R3MeshEdge *> *edges) const
{
  // Initialize return value
  RNLength destination_distance = RN_INFINITY;
  R3Point destination_position = VertexPosition(destination_vertex);

  // Allocate temporary data
  DijkstraData *vertex_data = new DijkstraData [ NVertices() ];
  if (!vertex_data) {
    fprintf(stderr, "Unable to allocate temporary data for geodesic distances.\n");
    return RN_INFINITY;
  }

  // Initialize all data
  for (int i = 0; i < NVertices(); i++) {
    vertex_data[i].vertex = Vertex(i);
    vertex_data[i].edge = NULL;
    vertex_data[i].distance_to_destination = RN_UNKNOWN;
    vertex_data[i].heappointer = NULL;
    vertex_data[i].distance = FLT_MAX;
  }

  // Initialize priority queue
  DijkstraData tmp;
  static RNHeap<DijkstraData *> heap(&tmp, &(tmp.distance), &(tmp.heappointer));
  DijkstraData *data = &vertex_data[ VertexID(source_vertex) ];
  data->distance_to_destination = R3Distance(VertexPosition(source_vertex), destination_position);
  data->distance = data->distance_to_destination;
  heap.Push(data);

  // Visit other vertices computing shortest distance
  while (!heap.IsEmpty()) {
    DijkstraData *data = heap.Pop();
    R3MeshVertex *vertex = data->vertex;
    if (vertex == destination_vertex) { destination_distance = data->distance; break; }
    for (int i = 0; i < VertexValence(vertex); i++) {
      R3MeshEdge *edge = EdgeOnVertex(vertex, i);
      R3MeshVertex *neighbor_vertex = VertexAcrossEdge(edge, vertex);
      DijkstraData *neighbor_data = &vertex_data [ VertexID(neighbor_vertex) ];
      if (neighbor_data->distance_to_destination == RN_UNKNOWN) {
        const R3Point& neighbor_position = VertexPosition(neighbor_vertex);
        neighbor_data->distance_to_destination = R3Distance(neighbor_position, destination_position);
      }
      RNScalar old_distance = neighbor_data->distance;
      RNScalar new_distance = EdgeLength(edge) + data->distance - data->distance_to_destination + neighbor_data->distance_to_destination;
      if (new_distance < old_distance) {
        neighbor_data->edge = edge;
        neighbor_data->distance = new_distance;
        if (old_distance < FLT_MAX) heap.Update(neighbor_data);
        else heap.Push(neighbor_data);
      }
    }
  }

  // Fill in path of edges
  if (edges) {
    const R3MeshVertex *vertex = destination_vertex;
    while (vertex != source_vertex) {
      DijkstraData *data = &vertex_data [ VertexID(vertex) ];
      if (!data->edge) break;
      edges->Insert(data->edge);
      vertex = VertexAcrossEdge(data->edge, vertex);
    }
  }

  // Delete temporary data
  heap.Empty();
  delete [] vertex_data;

  // Return distance
  return destination_distance;
}



RNLength *R3Mesh::
DijkstraDistances(const R3MeshVertex *source_vertex, RNLength max_distance, R3MeshEdge **edges) const
{
  // Compute dijkstra distances from source vertex
  RNArray<R3MeshVertex *> source_vertices;
  source_vertices.Insert((R3MeshVertex *) source_vertex);
  return DijkstraDistances(source_vertices, max_distance, edges);
}



RNLength *R3Mesh::
DijkstraDistances(const RNArray<R3MeshVertex *>& source_vertices, RNLength max_distance, R3MeshEdge **edges) const
{
  // Allocate array of distances (to return)
  RNLength *distances = new RNLength [ NVertices() ];
  if (!distances) {
    fprintf(stderr, "Unable to allocate array of distances.\n");
    return NULL;
  }

  // Allocate temporary data
  DijkstraData *vertex_data = new DijkstraData [ NVertices() ];
  if (!vertex_data) {
    fprintf(stderr, "Unable to allocate temporary data for geodesic distances.\n");
    delete [] distances;
    return NULL;
  }

  // Initialize all data
  for (int i = 0; i < NVertices(); i++) {
    vertex_data[i].vertex = Vertex(i);
    vertex_data[i].edge = NULL;
    vertex_data[i].distance_to_destination = RN_UNKNOWN;
    vertex_data[i].heappointer = NULL;
    vertex_data[i].distance = FLT_MAX;
    distances[i] = FLT_MAX;
  }

  // Initialize priority queue
  DijkstraData tmp;
  static RNHeap<DijkstraData *> heap(&tmp, &(tmp.distance), &(tmp.heappointer));
  for (int i = 0; i < source_vertices.NEntries(); i++) {
    R3MeshVertex *source_vertex = source_vertices[i];
    DijkstraData *data = &vertex_data[ VertexID(source_vertex) ];
    data->distance = 0;
    heap.Push(data);
  }

  // Visit vertices computing shortest distance to closest source vertex
  while (!heap.IsEmpty()) {
    DijkstraData *data = heap.Pop();
    R3MeshVertex *vertex = data->vertex;
    if ((max_distance > 0) && (data->distance > max_distance)) break;
    distances[ VertexID(vertex) ] = data->distance;
    if (edges) edges[ VertexID(vertex) ] = data->edge;
    for (int i = 0; i < VertexValence(vertex); i++) {
      R3MeshEdge *edge = EdgeOnVertex(vertex, i);
      R3MeshVertex *neighbor_vertex = VertexAcrossEdge(edge, vertex);
      DijkstraData *neighbor_data = &vertex_data [ VertexID(neighbor_vertex) ];
      RNScalar old_distance = neighbor_data->distance;
      RNScalar new_distance = EdgeLength(edge) + data->distance;
      if (new_distance < old_distance) {
        neighbor_data->edge = edge;
        neighbor_data->distance = new_distance;
        if (old_distance < FLT_MAX) heap.Update(neighbor_data);
        else heap.Push(neighbor_data);
      }
    }
  }

  // Delete temporary data
  heap.Empty();
  delete [] vertex_data;

  // Return distances
  return distances;
}



RNLength *R3Mesh::
DijkstraDistances(const R3MeshVertex *source_vertex, RNLength max_distance, RNArray<R3MeshVertex *>& neighbor_vertices, R3MeshEdge **edges) const
{
  // Compute dijkstra distances from source vertex
  RNArray<R3MeshVertex *> source_vertices;
  source_vertices.Insert((R3MeshVertex *) source_vertex);
  return DijkstraDistances(source_vertices, max_distance, neighbor_vertices, edges);
}



RNLength *R3Mesh::
DijkstraDistances(const RNArray<R3MeshVertex *>& source_vertices, RNLength max_distance, RNArray<R3MeshVertex *>& neighbor_vertices, R3MeshEdge **edges) const
{
  // Empty array of neighbor vertices to be returned
  neighbor_vertices.Empty();

  // Allocate array of distances (to return)
  RNLength *distances = new RNLength [ NVertices() ];
  if (!distances) {
    fprintf(stderr, "Unable to allocate array of distances.\n");
    return NULL;
  }

  // Allocate temporary data
  DijkstraData *vertex_data = new DijkstraData [ NVertices() ];
  if (!vertex_data) {
    fprintf(stderr, "Unable to allocate temporary data for geodesic distances.\n");
    delete [] distances;
    return NULL;
  }

  // Increment mesh mark
  R3mesh_mark++;

  // Initialize priority queue
  DijkstraData tmp;
  static RNHeap<DijkstraData *> heap(&tmp, &(tmp.distance), &(tmp.heappointer));
  for (int i = 0; i < source_vertices.NEntries(); i++) {
    R3MeshVertex *source_vertex = source_vertices[i];
    source_vertex->mark = R3mesh_mark;
    DijkstraData *data = &vertex_data[ VertexID(source_vertex) ];
    data->vertex = source_vertex;
    data->edge = NULL;
    data->distance_to_destination = 0;
    data->heappointer = NULL;
    data->distance = 0;
    heap.Push(data);
  }

  // Visit vertices computing shortest distance to closest source vertex
  while (!heap.IsEmpty()) {
    DijkstraData *data = heap.Pop();
    R3MeshVertex *vertex = data->vertex;
    if ((max_distance > 0) && (data->distance > max_distance)) break;
    neighbor_vertices.Insert(vertex);
    distances[ VertexID(vertex) ] = data->distance;
    if (edges) edges[ VertexID(vertex) ] = data->edge;
    for (int i = 0; i < VertexValence(vertex); i++) {
      R3MeshEdge *edge = EdgeOnVertex(vertex, i);
      R3MeshVertex *neighbor_vertex = VertexAcrossEdge(edge, vertex);
      DijkstraData *neighbor_data = &vertex_data [ VertexID(neighbor_vertex) ];
      RNScalar new_distance = EdgeLength(edge) + data->distance;
      if (VertexMark(neighbor_vertex) != R3mesh_mark) {
        neighbor_vertex->mark = R3mesh_mark;
        neighbor_data->vertex = neighbor_vertex;
        neighbor_data->edge = edge;
        neighbor_data->distance_to_destination = RN_UNKNOWN;
        neighbor_data->heappointer = NULL;
        neighbor_data->distance = new_distance;
        heap.Push(neighbor_data);
      }
      else {
        RNScalar old_distance = neighbor_data->distance;
        if (new_distance < old_distance) {
          neighbor_data->edge = edge;
          neighbor_data->distance = new_distance;
          heap.Update(neighbor_data);
        }
      }
    }
  }

  // Delete temporary data
  heap.Empty();
  delete [] vertex_data;

  // Return distances
  return distances;
}



////////////////////////////////////////////////////////////////////////
// PATH STUFF
////////////////////////////////////////////////////////////////////////

RNLength R3Mesh::
TracePath(const R3MeshVertex *source_vertex, const R3Vector& tangent_direction, RNLength max_distance,
  R3Point *final_position, R3MeshFace **final_face, R3MeshIntersection *intersections, int *nintersections) const
{
  // Increment mesh mark (used to avoid visiting same edge/face twice in recursive function)
  R3Mesh *mesh = (R3Mesh *) this;
  R3mesh_mark++;

  // Initialize result
  if (final_position) *final_position = VertexPosition(source_vertex);
  if (final_face) *final_face = NULL;
  if (nintersections) *nintersections = 0;
  if (max_distance <= 0) return 0;

  // Get convenient variables
  const R3Point& source_position = VertexPosition(source_vertex);
  R3Vector source_normal = VertexNormal(source_vertex);
  R3Vector split_direction = tangent_direction % source_normal;
  RNLength split_direction_length = split_direction.Length();
  if (split_direction_length == 0) return 0;
  split_direction /= split_direction_length;
  R3Plane split_plane(source_position, split_direction);

  // Find first face, edge, and position
  R3MeshFace *cur_face = NULL;
  R3MeshEdge *cur_edge = NULL;
  R3Point cur_position = R3zero_point;
  for (int i = 0; i < VertexValence(source_vertex); i++) {
    R3MeshEdge *cw_edge = EdgeOnVertex(source_vertex, i);
    R3MeshEdge *ccw_edge = EdgeOnVertex(source_vertex, cw_edge, RN_CCW);
    R3MeshFace *face = FaceOnVertex(source_vertex, cw_edge, RN_CCW);
    if (!face) continue;
    R3MeshEdge *edge = EdgeAcrossFace(face, source_vertex);

    // Check position of cw_vertex
    R3MeshVertex *cw_vertex = VertexAcrossEdge(cw_edge, source_vertex);
    const R3Point& cw_position = VertexPosition(cw_vertex);
    RNScalar cw_d = R3SignedDistance(split_plane, cw_position);
    if (cw_d == 0) {
      cur_position = cw_position;
      cur_face = face;
      cur_edge = edge;
      break;
    }

    // Check position of ccw_vertex
    R3MeshVertex *ccw_vertex = VertexAcrossEdge(ccw_edge, source_vertex);
    const R3Point& ccw_position = VertexPosition(ccw_vertex);
    RNScalar ccw_d = R3SignedDistance(split_plane, ccw_position);
    if (ccw_d == 0) {
      cur_position = ccw_position;
      cur_face = face;
      cur_edge = edge;
      break;
    }

    // Check if source direction splits cw and ccw vertices
    if ((cw_d > 0)  && (ccw_d < 0)) {
      RNScalar t = cw_d / (cw_d - ccw_d);
      cur_position = cw_position + t * (ccw_position - cw_position);
      cur_face = face;
      cur_edge = edge;
      break;
    }
  }      

  // Get next distance and direction
  RNLength cur_distance = R3Distance(cur_position, source_position);
  R3Vector cur_direction = cur_position - source_position; cur_direction.Normalize();

  // Check if found next face
  if (!cur_face) return 0;

  // Check if max_distance lies within next face
  if (cur_distance >= max_distance) {
    RNScalar t = max_distance / cur_distance;
    if (final_position) *final_position = (1-t)*source_position + t*cur_position;
    if (final_face) *final_face = cur_face;
    return max_distance;
  }

  // Trace path
  while (TRUE) {
    // Mark cur edge and face
    mesh->SetEdgeMark(cur_edge, R3mesh_mark);
    mesh->SetFaceMark(cur_face, R3mesh_mark);

    // Add intersection
    if (intersections && nintersections) {
      intersections[*nintersections].type = R3_MESH_EDGE_TYPE;
      intersections[*nintersections].vertex = NULL;
      intersections[*nintersections].edge = cur_edge;
      intersections[*nintersections].face = cur_face;
      intersections[*nintersections].point = cur_position;
      intersections[*nintersections].t = cur_distance;
      (*nintersections)++;
    }

    // Get next face
    R3MeshFace *next_face = FaceAcrossEdge(cur_edge, cur_face);
    if (!next_face) {
      if (final_position) *final_position = cur_position;
      if (final_face) *final_face = cur_face;
      return cur_distance;
    }
    
    // Get next direction 
    R3Vector cur_face_normal = FaceNormal(cur_face);
    R3Vector next_face_normal = FaceNormal(next_face);
    R4Matrix matrix = R4identity_matrix;
    matrix.Rotate(cur_face_normal, next_face_normal);
    R3Vector next_direction = matrix * cur_direction;
    next_direction.Normalize();
    
    // Get plane splitting next face along next_direction
    R3Vector split_plane_normal = next_face_normal % next_direction;  
    split_plane_normal.Normalize();
    R3Plane split_plane(cur_position, split_plane_normal);
    
    // Get next edge and position
    R3MeshVertex *vertex = VertexAcrossFace(next_face, cur_edge);
    R3Point next_position = VertexPosition(vertex);
    RNScalar d = R3SignedDistance(split_plane, next_position);
    R3MeshEdge *next_edge = EdgeOnFace(next_face, cur_edge, (d < 0) ? RN_CW : RN_CCW);
    if (d != 0) {
      // Split does not hit vertex
      R3Vector next_edge_vector = EdgeVector(next_edge);
      RNScalar denom = split_plane.Normal().Dot(next_edge_vector);
      if (denom != 0) {
        R3Point next_edge_start = VertexPosition(VertexOnEdge(next_edge, 0));
        RNScalar t = -(R3SignedDistance(split_plane, next_edge_start)) / denom;
        next_position = next_edge_start + next_edge_vector * t;
      }
    }

    // Get next distance
    RNLength next_distance = cur_distance + R3Distance(cur_position, next_position);
    if (next_distance >= max_distance) {
      assert(next_distance > cur_distance);
      RNScalar t = (max_distance - cur_distance) / (next_distance - cur_distance);
      if (final_position) *final_position = (1-t)*cur_position + t*next_position;
      if (final_face) *final_face = next_face;
      return max_distance;
    }

    // Go to next face
    cur_edge = next_edge;
    cur_face = next_face;
    cur_position = next_position;
    cur_direction = next_direction;
    cur_distance = next_distance;
  }

  // Return error
  return -1;
}



////////////////////////////////////////////////////////////////////////
// SHAPE CREATION FUNCTIONS
////////////////////////////////////////////////////////////////////////

void R3Mesh::
CreateBox(const R3Box& box)
{
  // Create box 
  R3MeshVertex *v1 = CreateVertex(box.Corner(RN_NNP_OCTANT));
  R3MeshVertex *v2 = CreateVertex(box.Corner(RN_PNP_OCTANT));
  R3MeshVertex *v3 = CreateVertex(box.Corner(RN_NPP_OCTANT));
  R3MeshVertex *v4 = CreateVertex(box.Corner(RN_PPP_OCTANT));
  R3MeshVertex *v5 = CreateVertex(box.Corner(RN_NPN_OCTANT));
  R3MeshVertex *v6 = CreateVertex(box.Corner(RN_PPN_OCTANT));
  R3MeshVertex *v7 = CreateVertex(box.Corner(RN_NNN_OCTANT));
  R3MeshVertex *v8 = CreateVertex(box.Corner(RN_PNN_OCTANT));
  CreateFace(v1, v2, v4);
  CreateFace(v1, v4, v3);
  CreateFace(v3, v4, v6);
  CreateFace(v3, v6, v5);
  CreateFace(v5, v6, v8);
  CreateFace(v5, v8, v7);
  CreateFace(v7, v8, v2);
  CreateFace(v7, v2, v1);
  CreateFace(v2, v8, v6);
  CreateFace(v2, v6, v4);
  CreateFace(v7, v1, v3);
  CreateFace(v7, v3, v5);
}



void R3Mesh::
CreateOrientedBox(const R3OrientedBox& box)
{
  // Create mesh elements for an oriented box and add to mesh
  RNAbort("Not implemented");
}



void R3Mesh::
CreateSphere(const R3Sphere& sphere, RNLength vertex_spacing)
{
  // Create mesh elements for a sphere and add to mesh
  RNAbort("Not implemented");
}



void R3Mesh::
CreateEllipsoid(const R3Ellipsoid& ellipsoid, RNLength vertex_spacing)
{
  // Create mesh elements for a ellipsoid and add to mesh
  RNAbort("Not implemented");
}



void R3Mesh::
CreateCone(const R3Cone& cone, RNLength vertex_spacing)
{
  // Create mesh elements for a cone and add to mesh
  RNAbort("Not implemented");
}



void R3Mesh::
CreateCylinder(const R3Cylinder& cylinder, RNLength vertex_spacing)
{
  // Determine number of vertices per circle
  int n = 32;
  if (vertex_spacing > 0) {
    n = (int) (2.0 * RN_PI * cylinder.Radius() / vertex_spacing);
    if (n < 8) n = 8;
  }

  // Allocate vertex pointers
  R3MeshVertex ***vertices = new R3MeshVertex **[2];
  vertices[0] = new R3MeshVertex * [ n ];
  vertices[1] = new R3MeshVertex * [ n ];

  // Create vertices
  for (int i = 0; i < n; i++) {
    RNAngle a = i*RN_TWO_PI/n;
    RNScalar x = cylinder.Radius()*cos(a);
    RNScalar y = cylinder.Radius()*sin(a);
    vertices[0][i] = CreateVertex(R3Point(x, y, -0.5*cylinder.Height()));
    vertices[1][i] = CreateVertex(R3Point(x, y,  0.5*cylinder.Height()));
  }

  // Create bottom faces
  R3MeshVertex *v0c = CreateVertex(R3Point(0, 0, -0.5*cylinder.Height()));
  for (int i = 0; i < n; i++) {
    R3MeshVertex *v00 = vertices[0][i];
    R3MeshVertex *v01 = vertices[0][(i+1)%n];
    CreateFace(v0c, v01, v00);
  }
  
  // Create top faces
  R3MeshVertex *v1c = CreateVertex(R3Point(0, 0, 0.5*cylinder.Height()));
  for (int i = 0; i < n; i++) {
    R3MeshVertex *v10 = vertices[1][i];
    R3MeshVertex *v11 = vertices[1][(i+1)%n];
    CreateFace(v1c, v10, v11);
  }
  
  // Create side faces
  for (int i = 0; i < n; i++) {
    R3MeshVertex *v00 = vertices[0][i];
    R3MeshVertex *v01 = vertices[0][(i+1)%n];
    R3MeshVertex *v10 = vertices[1][i];
    R3MeshVertex *v11 = vertices[1][(i+1)%n];
    CreateFace(v00, v01, v11);
    CreateFace(v00, v11, v10);
  }

  // Delete vertex pointers
  delete [] vertices[0];
  delete [] vertices[1];
  delete [] vertices;
}



void R3Mesh::
CreateCircle(const R3Circle& circle, RNLength vertex_spacing)
{
  // Create mesh elements for a circle and add to mesh
  RNAbort("Not implemented");
}



void R3Mesh::
CreateEllipse(const R3Ellipse& ellipse, RNLength vertex_spacing)
{
  // Create mesh elements for a ellipse and add to mesh
  RNAbort("Not implemented");
}



void R3Mesh::
CreateRectangle(const R3Rectangle& rectangle)
{
  // Create mesh elements for a rectangle and add to mesh
  RNAbort("Not implemented");
}



void R3Mesh::
CreateTriangle(const R3Triangle& triangle)
{
  // Create mesh elements for a triangle and add to mesh
  RNAbort("Not implemented");
}



void R3Mesh::
CreateTriangleArray(const R3TriangleArray& triangles)
{
  // Create mesh elements for a triangle array and add to mesh
  RNAbort("Not implemented");
}



void R3Mesh::
CreateCopy(const R3Mesh& mesh)
{
  // Copy vertices
  RNArray<R3MeshVertex *> verts;
  for (int i = 0; i < mesh.NVertices(); i++) {
    R3MeshVertex *vertex = mesh.Vertex(i);
    const R3Point& position = mesh.VertexPosition(vertex);
    const R3Vector& normal = mesh.VertexNormal(vertex);
    const RNRgb& color = mesh.VertexColor(vertex);
    const R2Point& texcoords = mesh.VertexTextureCoords(vertex);
    R3MeshVertex *copy_vertex = this->CreateVertex(position, normal, color, texcoords);
    if (!copy_vertex) return;
    verts.Insert(copy_vertex);
  }

  // Copy edges
  for (int i = 0; i < mesh.NEdges(); i++) {
    R3MeshEdge *edge = mesh.Edge(i);
    R3MeshVertex *v0 = mesh.VertexOnEdge(edge, 0);
    R3MeshVertex *v1 = mesh.VertexOnEdge(edge, 1);
    int i0 = mesh.VertexID(v0);
    int i1 = mesh.VertexID(v1);
    R3MeshVertex *copy_v0 = verts.Kth(i0);
    R3MeshVertex *copy_v1 = verts.Kth(i1);
    R3MeshEdge *copy_edge = this->CreateEdge(copy_v0, copy_v1);
    if (!copy_edge) return;
  }

  // Copy faces
  for (int i = 0; i < mesh.NFaces(); i++) {
    R3MeshFace *face = mesh.Face(i);
    R3MeshVertex *v0 = mesh.VertexOnFace(face, 0);
    R3MeshVertex *v1 = mesh.VertexOnFace(face, 1);
    R3MeshVertex *v2 = mesh.VertexOnFace(face, 2);
    int i0 = mesh.VertexID(v0);
    int i1 = mesh.VertexID(v1);
    int i2 = mesh.VertexID(v2);
    R3MeshVertex *copy_v0 = verts.Kth(i0);
    R3MeshVertex *copy_v1 = verts.Kth(i1);
    R3MeshVertex *copy_v2 = verts.Kth(i2);
    R3MeshFace *copy_face = this->CreateFace(copy_v0, copy_v1, copy_v2);
    if (!copy_face) return;
    this->SetFaceMaterial(copy_face, mesh.FaceMaterial(face));
    this->SetFaceSegment(copy_face, mesh.FaceSegment(face));
    this->SetFaceCategory(copy_face, mesh.FaceCategory(face));
  }
}



////////////////////////////////////////////////////////////////////////
// I/O FUNCTIONS
////////////////////////////////////////////////////////////////////////

int R3Mesh::
ReadFile(const char *filename)
{
  // Parse input filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    printf("Filename %s has no extension (e.g., .ply)\n", filename);
    return 0;
  }

  // Read file of appropriate type
  if (!strncmp(extension, ".obj", 4)) {
    if (!ReadObjFile(filename)) return 0;
  }
  else if (!strncmp(extension, ".off", 4)) {
    if (!ReadOffFile(filename)) return 0;
  }
  else if (!strncmp(extension, ".ray", 4))  {
    if (!ReadRayFile(filename)) return 0;
  }
  else if (!strncmp(extension, ".ply", 4)) { 
    if (!ReadPlyFile(filename)) return 0;
  }
  else if (!strncmp(extension, ".cat", 4)) { 
    if (!ReadCattFile(filename)) return 0;
  }
  else if (!strncmp(extension, ".m", 4)) { 
    if (!ReadHoppeFile(filename)) return 0;
  }  
  else if (!strncmp(extension, ".ifs", 4)) { 
    if (!ReadIfsFile(filename)) return 0;
  }
  else if (!strncmp(extension, ".stl", 4))  {
    if (!ReadSTLFile(filename)) return 0;
  }
  else if (!strncmp(extension, ".wrl", 4))  {
    if (!ReadVRMLFile(filename)) return 0;
  }
  else {
    RNFail("Unable to read file %s (unrecognized extension: %s)\n", filename, extension);
    return 0;
  }

  // Set mesh name, if it doesn't already have one
  if (strlen(Name()) == 0) SetName(filename);

  // Return success
  return 1;
}



int R3Mesh::
ReadObjFile(const char *filename)
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "r"))) {
    RNFail("Unable to open file %s", filename);
    return 0;
  }

  // Read body
  char buffer[1024];
  int line_count = 0;
  RNArray<R2Point *> texture_coords;
  RNArray<R3Vector *> normals;
  RNArray<R3MeshVertex *> verts;
  RNArray<R3MeshVertex *> degenerate_triangle_vertices;
  while (fgets(buffer, 1023, fp)) {
    // Increment line counter
    line_count++;

    // Skip white space
    char *bufferp = buffer;
    while (isspace(*bufferp)) bufferp++;

    // Skip blank lines and comments
    if (*bufferp == '#') continue;
    if (*bufferp == '\0') continue;

    // Get keyword
    char keyword[80];
    if (sscanf(bufferp, "%s", keyword) != 1) {
      RNFail("Syntax error on line %d in file %s", line_count, filename);
      return 0;
    }

    // Check keyword
    if (!strcmp(keyword, "v")) {
      // Read vertex coordinates
      double x, y, z;
      if (sscanf(bufferp, "%s%lf%lf%lf", keyword, &x, &y, &z) != 4) {
        RNFail("Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Create vertex
      R3MeshVertex *v = CreateVertex(R3Point(x, y, z));
      verts.Insert(v);
    }
    else if (!strcmp(keyword, "vt")) {
      // Read texture coordinates
      double u, v;
      if (sscanf(bufferp, "%s%lf%lf", keyword, &u, &v) != 3) {
        fprintf(stderr, "Syntax error on line %d in OBJ file", line_count);
        return 0;
      }

      // Create texture coordinates
      R2Point *vt = new R2Point(u, v);
      texture_coords.Insert(vt);
    }
    else if (!strcmp(keyword, "vn")) {
      // Read texture coordinates
      double x, y, z;
      if (sscanf(bufferp, "%s%lf%lf%lf", keyword, &x, &y, &z) != 4) {
        fprintf(stderr, "Syntax error on line %d in OBJ file", line_count);
        return 0;
      }

      // Create normal
      R3Vector *vn = new R3Vector(x, y, z);
      normals.Insert(vn);
    }
    else if (!strcmp(keyword, "f")) {
      // Read vertex indices
      int quad = 1;
      char s[4][128] = { { '\0' }, { '\0' }, { '\0' },{ '\0' } }; 
      if (sscanf(bufferp, "%s%s%s%s%s", keyword, s[0], s[1], s[2], s[3]) != 5) {
        quad = 0;;
        if (sscanf(bufferp, "%s%s%s%s", keyword, s[0], s[1], s[2]) != 4) {
          fprintf(stderr, "Syntax error on line %d in OBJ file", line_count);
          return 0;
        }
      }

      // Parse vertex indices
      int n = (quad) ? 4 : 3;
      R3MeshVertex *v[4] = { NULL, NULL, NULL, NULL };
      for (int i = 0; i < n; i++) {
        char *sv = s[i];
        char *st = strchr(sv, '/');
        if (st) *(st++) = 0;
        char *sn = (st) ? strchr(st, '/') : NULL;
        if (sn) *(sn++) = 0;
        if (sn && (strlen(sn) == 0)) sn = NULL; 
        if (st && (strlen(st) == 0)) st = NULL;
        int vi = (sv) ? atoi(sv) : 0;
        int ti = (st) ? atoi(st) : 0;
        int ni = (sn) ? atoi(sn) : 0;
        v[i] = verts.Kth(vi-1);
        R2Point texcoords(0,0);
        R3Vector normal(0, 0, 0);
        if ((ti > 0) && ((ti-1) < texture_coords.NEntries())) texcoords = *(texture_coords.Kth(ti-1));
        if ((ni > 0) && ((ni-1) < normals.NEntries())) normal = *(normals.Kth(ni-1));
        if (((ti > 0) && !R2Contains(texcoords, VertexTextureCoords(v[i]))) ||
            ((ni > 0) && !R3Contains(normal, VertexNormal(v[i])))) {
          v[i] = CreateVertex(VertexPosition(v[i]), normal, RNgray_rgb, texcoords);
        }
      }

      // Check vertices
      if ((v[0] == v[1]) || (v[1] == v[2]) || (v[0] == v[2])) continue;
      if ((quad) && ((v[3] == v[0]) || (v[3] == v[1]) || (v[3] == v[2]))) quad = 0;

      // Create first triangle
      if (RNIsPositive(R3Distance(VertexPosition(v[0]), VertexPosition(v[1]))) &&
          RNIsPositive(R3Distance(VertexPosition(v[1]), VertexPosition(v[2]))) &&
          RNIsPositive(R3Distance(VertexPosition(v[2]), VertexPosition(v[0])))) {
        if (!CreateFace(v[0], v[1], v[2])) {
          // Must have been degeneracy (e.g., flips or three faces sharing an edge)
          // Remember for later processing (to preserve vertex indices)
          degenerate_triangle_vertices.Insert(v[0]);
          degenerate_triangle_vertices.Insert(v[1]);
          degenerate_triangle_vertices.Insert(v[2]);
        }
      }

      // Create second triangle
      if (quad) {
        if (RNIsPositive(R3Distance(VertexPosition(v[0]), VertexPosition(v[2]))) &&
            RNIsPositive(R3Distance(VertexPosition(v[2]), VertexPosition(v[3]))) &&
            RNIsPositive(R3Distance(VertexPosition(v[0]), VertexPosition(v[3])))) {
          if (!CreateFace(v[0], v[2], v[3])) {
            // Must have been degeneracy (e.g., flips or three faces sharing an edge)
            // Remember for later processing (to preserve vertex indices)
            degenerate_triangle_vertices.Insert(v[0]);
            degenerate_triangle_vertices.Insert(v[2]);
            degenerate_triangle_vertices.Insert(v[3]);
          }
        }
      }
    }
  }

  // Create degenerate triangles
  for (int i = 0; i <= degenerate_triangle_vertices.NEntries()-3; i+=3) {
    R3MeshVertex *v1 = degenerate_triangle_vertices.Kth(i+0);
    R3MeshVertex *v2 = degenerate_triangle_vertices.Kth(i+1);
    R3MeshVertex *v3 = degenerate_triangle_vertices.Kth(i+2);
    if (!CreateFace(v1, v2, v3)) {
      if (!CreateFace(v1, v3, v2)) {
        // Note: these vertices are allocated separately, and so they will not be deleted (memory leak)
        R3MeshVertex *v1a = CreateVertex(VertexPosition(v1));
        R3MeshVertex *v2a = CreateVertex(VertexPosition(v2));
        R3MeshVertex *v3a = CreateVertex(VertexPosition(v3));
        CreateFace(v1a, v2a, v3a);
      }
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R3Mesh::
ReadOffFile(const char *filename)
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "r"))) {
    fprintf(stderr, "Unable to open file %s\n", filename);
    return 0;
  }

  // Read file
  int nverts = 0;
  int nfaces = 0;
  int nedges = 0;
  int line_count = 0;
  int vertex_count = 0;
  int face_count = 0;
  char buffer[1024];
  char header[64];
  RNArray<R3MeshVertex *> degenerate_triangle_vertices;
  while (fgets(buffer, 1023, fp)) {
    // Increment line counter
    line_count++;

    // Skip white space
    char *bufferp = buffer;
    while (isspace(*bufferp)) bufferp++;

    // Skip blank lines and comments
    if (*bufferp == '#') continue;
    if (*bufferp == '\0') continue;

    // Check section
    if (nverts == 0) {
      // Read header keyword
      if (strstr(bufferp, "OFF")) {
        // Check if counts are on first line
        int tmp;
        if (sscanf(bufferp, "%s%d%d%d", header, &tmp, &nfaces, &nedges) == 4) {
          nverts = tmp;
        }
      }
      else {
        // Read counts from second line
        if ((sscanf(bufferp, "%d%d%d", &nverts, &nfaces, &nedges) != 3) || (nverts == 0)) {
          RNFail("Syntax error reading header on line %d in file %s\n", line_count, filename);
          fclose(fp);
          return 0;
        }
      }
    }
    else if (vertex_count < nverts) {
      // Read vertex coordinates
      double x, y, z;
      if (sscanf(bufferp, "%lf%lf%lf", &x, &y, &z) != 3) {
        RNFail("Syntax error with vertex coordinates on line %d in file %s\n", line_count, filename);
        fclose(fp);
        return 0;
      }

      // Create vertex
      CreateVertex(R3Point(x, y, z));

      // Increment counter
      vertex_count++;
    }
    else if (face_count < nfaces) {
      // Read number of vertices in face 
      int face_nverts = 0;
      bufferp = strtok(bufferp, " \t");
      if (bufferp) face_nverts = atoi(bufferp);
      else {
        RNFail("Syntax error with face on line %d in file %s\n", line_count, filename);
        fclose(fp);
        return 0;
      }

      // Read vertex indices for face
      R3MeshVertex *v1 = NULL;
      R3MeshVertex *v2 = NULL;
      R3MeshVertex *v3 = NULL;
      for (int i = 0; i < face_nverts; i++) {
        bufferp = strtok(NULL, " \t");
        if (bufferp) {
          R3MeshVertex *v = Vertex(atoi(bufferp));
          if (!v1) v1 = v;
          else v3 = v;
        }
        else {
          RNFail("Syntax error with face on line %d in file %s\n", line_count, filename);
          fclose(fp);
          return 0;
        }

        // Create triangle
        if (v1 && v2 && v3 && (v1 != v2) && (v2 != v3) && (v1 != v3)) {
          if (!CreateFace(v1, v2, v3)) {
            // Must have been degeneracy (e.g., flips or three faces sharing an edge)
            // Remember for later processing (to preserve vertex indices)
            degenerate_triangle_vertices.Insert(v1);
            degenerate_triangle_vertices.Insert(v2);
            degenerate_triangle_vertices.Insert(v3);
          }
        }

        // Move to next triangle
        v2 = v3;
      }

      // Increment counter
      face_count++;
    }
    else {
      // Should never get here
      RNFail("Found extra text starting at line %d in file %s\n", line_count, filename);
      break;
    }
  }

  // Create degenerate triangles
  for (int i = 0; i <= degenerate_triangle_vertices.NEntries()-3; i+=3) {
    R3MeshVertex *v1 = degenerate_triangle_vertices.Kth(i+0);
    R3MeshVertex *v2 = degenerate_triangle_vertices.Kth(i+1);
    R3MeshVertex *v3 = degenerate_triangle_vertices.Kth(i+2);
    if (!CreateFace(v1, v2, v3)) {
      if (!CreateFace(v1, v3, v2)) {
        // Note: these vertices are allocated separately, and so they will not be deleted (memory leak)
        R3MeshVertex *v1a = CreateVertex(VertexPosition(v1));
        R3MeshVertex *v2a = CreateVertex(VertexPosition(v2));
        R3MeshVertex *v3a = CreateVertex(VertexPosition(v3));
        CreateFace(v1a, v2a, v3a);
      }
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R3Mesh::
ReadRayFile(const char *filename)
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "r"))) {
    RNFail("Unable to open file %s", filename);
    return 0;
  }

  // Read body
  char cmd[128];
  int triangle_count = 0;
  int command_number = 1;
  RNArray<R3MeshVertex *> degenerate_triangle_vertices;
  while (fscanf(fp, "%s", cmd) == 1) {
    if (!strcmp(cmd, "#vertex")) {
      // Read data
      double px, py, pz;
      double nx, ny, nz;
      double ts, tt;
      if (fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf%lf", &px, &py, &pz, &nx, &ny, &nz, &ts, &tt) != 8) {
        RNFail("Unable to read vertex at command %d in file %s", command_number, filename);
        return 0;
      }

      // Create vertex
      R3Point point(px, py, pz);
      R3Vector normal(nx, ny, nz);
      if (normal.IsZero()) CreateVertex(point);
      else CreateVertex(point, normal);
    }
    else if (!strcmp(cmd, "#shape_triangle")) {
      // Read data
      int m;
      int i1, i2, i3;
      if (fscanf(fp, "%d%d%d%d", &m, &i1, &i2, &i3) != 4) {
        RNFail("Unable to read triangle at command %d in file %s", command_number, filename);
        return 0;
      }

      // Get vertices
      R3MeshVertex *v1 = vertices.Kth(i1);
      R3MeshVertex *v2 = vertices.Kth(i2);
      R3MeshVertex *v3 = vertices.Kth(i3);

      // Check vertices
      if ((v1 == v2) || (v2 == v3) || (v1 == v3)) continue;

      // Create face
      R3MeshFace *face = CreateFace(v1, v2, v3);
      if (!face) {
        // Must have been degeneracy (e.g., flips or three faces sharing an edge)
        // Remember for later processing (to preserve vertex indices)
        degenerate_triangle_vertices.Insert(v1);
        degenerate_triangle_vertices.Insert(v2);
        degenerate_triangle_vertices.Insert(v3);
      }

      // Set face material
      if (face) SetFaceMaterial(face, m);

      // Increment triangle counter
      triangle_count++;
    }
	
    // Increment command number
    command_number++;
  }

  // Create degenerate triangles
  for (int i = 0; i <= degenerate_triangle_vertices.NEntries()-3; i+=3) {
    R3MeshVertex *v1 = degenerate_triangle_vertices.Kth(i+0);
    R3MeshVertex *v2 = degenerate_triangle_vertices.Kth(i+1);
    R3MeshVertex *v3 = degenerate_triangle_vertices.Kth(i+2);
    if (!CreateFace(v1, v2, v3)) {
      if (!CreateFace(v1, v3, v2)) {
        // Note: these vertices are allocated separately, and so they will not be deleted (memory leak)
        R3MeshVertex *v1a = CreateVertex(VertexPosition(v1));
        R3MeshVertex *v2a = CreateVertex(VertexPosition(v2));
        R3MeshVertex *v3a = CreateVertex(VertexPosition(v3));
        CreateFace(v1a, v2a, v3a);
      }
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



static void
CreatePlyRangeGridFace(R3Mesh *mesh, 
  R3MeshVertex *v0, R3MeshVertex *v1, R3MeshVertex *v2, 
  RNScalar max_aspect_ratio = 8)
{
  // Get vertex positions
  const R3Point& p0 = mesh->VertexPosition(v0);
  const R3Point& p1 = mesh->VertexPosition(v1);
  const R3Point& p2 = mesh->VertexPosition(v2);

  // Get edge lengths
  RNScalar a = R3Distance(p0, p1);
  if (RNIsZero(a)) return;
  RNScalar b = R3Distance(p1, p2);
  if (RNIsZero(b)) return;
  RNScalar c = R3Distance(p2, p0);
  if (RNIsZero(c)) return;

  // Compute aspect ratio
  RNScalar s = (a + b + c) / 2.0;
  RNScalar da = s - a;
  if (RNIsZero(da)) return;
  RNScalar db = s - b;
  if (RNIsZero(db)) return;
  RNScalar dc = s - c;
  if (RNIsZero(dc)) return;
  RNScalar aspect = (a * b * c) / (8.0 * da * db * dc);
  if (aspect > max_aspect_ratio) return;

  // Check edge length ratios
  if (a < b) { RNScalar swap = a; a = b; b = swap; }
  if (a < c) { RNScalar swap = a; a = c; c = swap; }
  if (b < c) { RNScalar swap = b; b = c; c = swap; }
  if (a / c > max_aspect_ratio) return;

  // Create face
  mesh->CreateFace(v0, v1, v2);
}



int R3Mesh::
ReadPlyFile(const char *filename)
{
  FILE *fp;
  int i,j;
  PlyFile *ply;
  int nelems;
  PlyProperty **plist;
  char **elist;
  int file_type;
  int nprops;
  int num_elems;
  char *elem_name;
  float version;
  int num_rows = 0;
  int num_cols = 0;

  typedef struct PlyVertex {
    float x, y, z;
    float nx, ny, nz;
    float tx, ty;
    unsigned char red, green, blue;
  } PlyVertex;

  typedef struct PlyFace {
    unsigned char nverts;
    int *verts;
    int material;
    int segment;
    int category;
  } PlyFace;

  typedef struct PlyRangeGrid {
    unsigned char nverts;
    int *verts;
  } PlyRangeGrid;

  // List of property information for a vertex 
  static PlyProperty vert_props[] = { 
    {(char *) "x", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex,x), 0, 0, 0, 0},
    {(char *) "y", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex,y), 0, 0, 0, 0},
    {(char *) "z", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex,z), 0, 0, 0, 0},
    {(char *) "nx", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex,nx), 0, 0, 0, 0},
    {(char *) "ny", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex,ny), 0, 0, 0, 0},
    {(char *) "nz", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex,nz), 0, 0, 0, 0},
    {(char *) "tx", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex,tx), 0, 0, 0, 0},
    {(char *) "ty", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex,ty), 0, 0, 0, 0},
    {(char *) "red", PLY_UCHAR, PLY_UCHAR, offsetof(PlyVertex,red), 0, 0, 0, 0},
    {(char *) "green", PLY_UCHAR, PLY_UCHAR, offsetof(PlyVertex,green), 0, 0, 0, 0},
    {(char *) "blue", PLY_UCHAR, PLY_UCHAR, offsetof(PlyVertex,blue), 0, 0, 0, 0}
  };

  // List of property information for a face 
  static PlyProperty face_props[] = { 
    {(char *) "vertex_indices", PLY_INT, PLY_INT, offsetof(PlyFace,verts), 1, PLY_UCHAR, PLY_UCHAR, offsetof(PlyFace,nverts)},
    {(char *) "vertex_index", PLY_INT, PLY_INT, offsetof(PlyFace,verts), 1, PLY_UCHAR, PLY_UCHAR, offsetof(PlyFace,nverts)},
    {(char *) "material_id", PLY_INT, PLY_INT, offsetof(PlyFace,material), 0, 0, 0, 0},
    {(char *) "segment_id", PLY_INT, PLY_INT, offsetof(PlyFace,segment), 0, 0, 0, 0},
    {(char *) "category_id", PLY_INT, PLY_INT, offsetof(PlyFace,category), 0, 0, 0, 0}
  };

  // List of property information for a range_grid
  static PlyProperty range_grid_props[] = { 
    {(char *) "vertex_indices", PLY_INT, PLY_INT, offsetof(PlyFace,verts), 1, PLY_UCHAR, PLY_UCHAR, offsetof(PlyFace,nverts)},
  };

  // Open file 
  fp = fopen(filename, "rb");
  if (!fp) {
    RNFail("Unable to open file: %s", filename);
    return 0;
  }

  // Read PLY header
  ply = ply_read (fp, &nelems, &elist);
  if (!ply) {
    RNFail("Unable to read ply file: %s", filename);
    fclose(fp);
    return 0;
  }
  
  // Get header info
  ply_get_info (ply, &version, &file_type);

  // Read all elements
  for (i = 0; i < nelems; i++) {
    // Get the description of the element 
    elem_name = elist[i];
    plist = ply_get_element_description (ply, elem_name, &num_elems, &nprops);

    // Check element type
    if (equal_strings ("vertex", elem_name)) {
      // Allocate block of vertices
      vertex_block = new R3MeshVertex [num_elems];

      // Resize array of vertices
      vertices.Resize(num_elems);

      // set up for getting vertex elements 
      RNBoolean has_normals = 0;
      RNBoolean has_texcoords = 0;
      RNBoolean has_colors = 0;
      for (j = 0; j < nprops; j++) {
	if (equal_strings("x", plist[j]->name)) ply_get_property (ply, elem_name, &vert_props[0]);
	else if (equal_strings("y", plist[j]->name)) ply_get_property (ply, elem_name, &vert_props[1]);
	else if (equal_strings("z", plist[j]->name)) ply_get_property (ply, elem_name, &vert_props[2]);
	else if (equal_strings("nx", plist[j]->name)) ply_get_property (ply, elem_name, &vert_props[3]); 
	else if (equal_strings("ny", plist[j]->name)) ply_get_property (ply, elem_name, &vert_props[4]); 
	else if (equal_strings("nz", plist[j]->name)) ply_get_property (ply, elem_name, &vert_props[5]);
	else if (equal_strings("tx", plist[j]->name)) ply_get_property (ply, elem_name, &vert_props[6]); 
	else if (equal_strings("ty", plist[j]->name)) ply_get_property (ply, elem_name, &vert_props[7]); 
	else if (equal_strings("red", plist[j]->name)) ply_get_property (ply, elem_name, &vert_props[8]); 
	else if (equal_strings("green", plist[j]->name)) ply_get_property (ply, elem_name, &vert_props[9]); 
	else if (equal_strings("blue", plist[j]->name)) ply_get_property (ply, elem_name, &vert_props[10]);
	if (equal_strings("nx", plist[j]->name)) has_normals = 1;
	else if (equal_strings("tx", plist[j]->name)) has_texcoords = 1;
	else if (equal_strings("red", plist[j]->name)) has_colors = 1;
      }

      // grab all the vertex elements 
      for (j = 0; j < num_elems; j++) {
        // Read vertex into local struct
        PlyVertex plyvertex;
        ply_get_element(ply, (void *) &plyvertex);

        // Create mesh vertex
        R3Point position(plyvertex.x, plyvertex.y, plyvertex.z);
        R3MeshVertex *v = CreateVertex(position, &vertex_block[j]);
        if (has_normals) {
          R3Vector normal(plyvertex.nx, plyvertex.ny, plyvertex.nz);
          SetVertexNormal(v, normal);
        }
        if (has_texcoords) {
          R2Point texcoords(plyvertex.tx, plyvertex.ty);
          SetVertexTextureCoords(v, texcoords);
        }
        if (has_colors) {
          RNRgb color(plyvertex.red/255.0, plyvertex.green/255.0, plyvertex.blue/255.0);
          SetVertexColor(v, color);
        }
      }
    }
    else if (equal_strings ("face", elem_name)) {
      // Resize array of faces
      faces.Resize(num_elems);

      // set up for getting face elements 
      for (j = 0; j < nprops; j++) {
	if (equal_strings("vertex_indices", plist[j]->name)) ply_get_property (ply, elem_name, &face_props[0]);
	else if (equal_strings("vertex_index", plist[j]->name)) ply_get_property (ply, elem_name, &face_props[1]);
	else if (equal_strings("material_id", plist[j]->name)) ply_get_property (ply, elem_name, &face_props[2]);
	else if (equal_strings("segment_id", plist[j]->name)) ply_get_property (ply, elem_name, &face_props[3]);
	else if (equal_strings("category_id", plist[j]->name)) ply_get_property (ply, elem_name, &face_props[4]);
      }

      // Create stuff for degenerate triangles
      RNArray<R3MeshVertex *> degenerate_triangle_vertices;
      int *degenerate_triangle_materials = new int [ num_elems ];
      int *degenerate_triangle_segments = new int [ num_elems ];
      int *degenerate_triangle_categories = new int [ num_elems ];
      int num_degenerate_triangles = 0;

      // grab all the face elements 
      for (j = 0; j < num_elems; j++) {
        // Read face into local struct
        PlyFace plyface;
        plyface.nverts = 0;
        plyface.verts = NULL;
        plyface.material = -1;
        plyface.segment = -1;
        plyface.category = -1;
        ply_get_element(ply, (void *) &plyface);

        // Create mesh face(s)
        R3MeshVertex *v1 = vertices[plyface.verts[0]];
        for (int k = 2; k < plyface.nverts; k++) {
          // Get vertices
          R3MeshVertex *v2 = vertices[plyface.verts[k-1]];
          R3MeshVertex *v3 = vertices[plyface.verts[k]];

          // Check plyface
          assert(plyface.verts[0] >= 0);
          assert(plyface.verts[k-1] >= 0);
          assert(plyface.verts[k] >= 0);
          assert(plyface.verts[0] < vertices.NEntries());
          assert(plyface.verts[k-1] < vertices.NEntries());
          assert(plyface.verts[k] < vertices.NEntries());

          // Check vertices
          if ((v1 == v2) || (v2 == v3) || (v1 == v3)) continue;

          // Create face
          R3MeshFace *f = CreateFace(v1, v2, v3);
          if (f) {
            // Set material/segment/category
            SetFaceMaterial(f, plyface.material);
            SetFaceSegment(f, plyface.segment);
            SetFaceCategory(f, plyface.category);
          }
          else {
            // Must have been degeneracy (e.g., three faces sharing an edge)
            // Remember for later processing (to preserve vertex indices)
            degenerate_triangle_vertices.Insert(v1);
            degenerate_triangle_vertices.Insert(v2);
            degenerate_triangle_vertices.Insert(v3);
            degenerate_triangle_materials[num_degenerate_triangles] = plyface.material;
            degenerate_triangle_segments[num_degenerate_triangles] = plyface.segment;
            degenerate_triangle_categories[num_degenerate_triangles] = plyface.category;
            num_degenerate_triangles++;
          }
        }

        // Free face data allocated by ply
        if (plyface.verts) free(plyface.verts);
      }
      
      // Create degenerate triangles (do this at end to preserve face ordering)
      for (int i = 0; i < num_degenerate_triangles; i++) {
        R3MeshVertex *v1 = degenerate_triangle_vertices.Kth(3*i + 0);
        R3MeshVertex *v2 = degenerate_triangle_vertices.Kth(3*i + 1);
        R3MeshVertex *v3 = degenerate_triangle_vertices.Kth(3*i + 2);
        R3MeshFace *f = CreateFace(v1, v2, v3);
        if (!f) {
          f = CreateFace(v1, v3, v2);
          if (!f) {
            // Note: these vertices are allocated separately, and so they will not be deleted (memory leak)
            R3MeshVertex *v1a = CreateVertex(VertexPosition(v1));
            R3MeshVertex *v2a = CreateVertex(VertexPosition(v2));
            R3MeshVertex *v3a = CreateVertex(VertexPosition(v3));
            f = CreateFace(v1a, v2a, v3a);
          }
        }
        if (f) {
          SetFaceMaterial(f, degenerate_triangle_materials[i]);
          SetFaceSegment(f, degenerate_triangle_segments[i]);
          SetFaceCategory(f, degenerate_triangle_categories[i]);
        }
      }

      // Delete memory for degenerate triangles
      delete [] degenerate_triangle_materials;
      delete [] degenerate_triangle_segments;
      delete [] degenerate_triangle_categories;
    }
    else if (equal_strings ("range_grid", elem_name)) {
      // Get num_cols and num_rows for range_grid
      int ninfo = 0;
      char **info = ply_get_obj_info(ply, &ninfo);
      for (int i = 0; i < ninfo; i++) {
        char *str = strdup(info[i]);
        char *strp = strtok(str, " \n\t");
        if (strp) {
          if (!strcmp(strp, "num_rows")) {
            strp = strtok(NULL, " \n\t");
            if (strp) num_rows = atoi(strp);
          }
          else if (!strcmp(strp, "num_cols")) {
            strp = strtok(NULL, " \n\t");
            if (strp) num_cols = atoi(strp);
          }
        }
        free(str);
      }

      // check num_cols and num_rows
      if ((num_cols > 0) && (num_rows > 0) && (num_cols * num_rows == num_elems)) {
        // Allocate vertex grid
        R3MeshVertex **vertex_grid = new R3MeshVertex * [num_elems];
        for (j = 0; j < num_elems; j++) vertex_grid[j] = NULL;

        // set up for getting range grid elements 
        for (j = 0; j < nprops; j++) {
          if (equal_strings("vertex_indices", plist[j]->name)) ply_get_property (ply, elem_name, &range_grid_props[0]);
        }

        // grab all the range grid vertices
        for (j = 0; j < num_elems; j++) {
          // Read face into local struct
          PlyRangeGrid ply_range_grid;
          ply_range_grid.nverts = 0;
          ply_range_grid.verts = NULL;
          ply_get_element(ply, (void *) &ply_range_grid);
          if (!ply_range_grid.verts) continue;
          if (ply_range_grid.nverts == 0) continue;
          vertex_grid[j] = vertices[ply_range_grid.verts[0]];
          free(ply_range_grid.verts);
        }

        // create faces
        for (int row = 0; row < num_rows-1; row++) {
          for (int col = 0; col < num_cols-1; col++) {
            R3MeshVertex *v00 = vertex_grid[(row)*num_cols + (col)];
            R3MeshVertex *v01 = vertex_grid[(row+1)*num_cols + (col)];
            R3MeshVertex *v10 = vertex_grid[(row)*num_cols + (col+1)];
            R3MeshVertex *v11 = vertex_grid[(row+1)*num_cols + (col+1)];
            if (v00) {
              if (v01) {
                if (v10) {
                  if (v11) {
                    CreatePlyRangeGridFace(this, v00, v10, v11);
                    CreatePlyRangeGridFace(this, v00, v11, v01);
                  }
                  else {
                    CreatePlyRangeGridFace(this, v00, v10, v01);
                  }
                }
                else {
                  if (v11) {
                    CreatePlyRangeGridFace(this, v00, v11, v01);
                  }
                }
              }
              else {
                if (v10) {
                  if (v11) {
                    CreatePlyRangeGridFace(this, v00, v10, v11);
                  }
                }
              }
            }
            else {
              if (v01) {
                if (v10) {
                  if (v11) {
                    CreatePlyRangeGridFace(this, v10, v11, v01);
                  }
                }
              }
            }
          }
        }

        // Delete vertex grid
        delete [] vertex_grid;
      }
    }
    else {
      ply_get_other_element (ply, elem_name, num_elems);
    }
  }

  // Close the file 
  ply_close (ply);

  // Return success
  return 1;
}



int R3Mesh::
ReadCattFile(const char *filename)
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "r"))) {
    RNFail("Unable to open file %s", filename);
    return 0;
  }

  // Read body
  char buffer[1024];
  int line_count = 0;
  RNBoolean corners = FALSE;
  RNBoolean planes = FALSE;
  RNArray<R3MeshVertex *> degenerate_triangle_vertices;
  while (fgets(buffer, 1023, fp)) {
    // Increment line counter
    line_count++;

    // Skip white space
    char *bufferp = buffer;
    while (isspace(*bufferp)) bufferp++;

    // Skip blank lines
    if (*bufferp == '\0') continue;

    // Get keyword
    char keyword[80];
    if (sscanf(bufferp, "%s", keyword) != 1) {
      RNFail("Syntax error on line %d in file %s", line_count, filename);
      return 0;
    }

    // Check keyword
    if (keyword[0] == '%') {
      corners = planes = FALSE;
      if (!strcmp(keyword, "%CORNERS")) corners = TRUE;
      else if (!strcmp(keyword, "%PLANES")) planes = TRUE;
      else if (!strcmp(keyword, "%EOF")) break;
      continue;
    }

    // Read data
    if (corners) {
      // Read vertex coordinates
      int id;
      double x, y, z;
      if (sscanf(bufferp, "%d%lf%lf%lf", &id, &x, &y, &z) != 4) {
        RNFail("Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Create vertex
      CreateVertex(R3Point(x, y, z));
    }
    else if (planes) {
      // Read plane header
      int id;
      char buffer1[256];
      char buffer2[256];
      if ((sscanf(bufferp, "%d%s%s", &id, buffer1, buffer2) != 3) ||
          (strcmp(buffer1, "/")) || (strcmp(buffer2, "/RIGID"))) {
        RNFail("Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Read vertices
      RNArray<R3MeshVertex *> vertices;
      if (fgets(buffer, 1023, fp)) {
        line_count++;
        bufferp = strtok(buffer, "\t ");
        while (bufferp) {
          int id = atoi(bufferp);
          if ((id <= 0) || (id > NVertices())) {
            RNFail("Bogus vertex index on line %d in file %s", line_count, filename);
            return 0;
          }
          R3MeshVertex *v = Vertex(id-1);
          assert(v);
          vertices.Insert(v);
          bufferp = strtok(NULL, "\t ");
        }
      }

      // Create face(s)
      for (int i = 2; i < vertices.NEntries(); i++) {
        if (vertices[0] == vertices[i-1]) continue;
        if (vertices[i-1] == vertices[i]) continue;
        if (vertices[0] == vertices[i]) continue;
        if (!CreateFace(vertices[0], vertices[i-1], vertices[i])) {
          // Must have been degeneracy (e.g., flips or three faces sharing an edge)
          // Remember for later processing (to preserve vertex indices)
          degenerate_triangle_vertices.Insert(vertices[0]);
          degenerate_triangle_vertices.Insert(vertices[i-1]);
          degenerate_triangle_vertices.Insert(vertices[i]);
        }
      }
    }
  }
    
  // Create degenerate triangles
  for (int i = 0; i <= degenerate_triangle_vertices.NEntries()-3; i+=3) {
    R3MeshVertex *v1 = degenerate_triangle_vertices.Kth(i+0);
    R3MeshVertex *v2 = degenerate_triangle_vertices.Kth(i+1);
    R3MeshVertex *v3 = degenerate_triangle_vertices.Kth(i+2);
    if (!CreateFace(v1, v2, v3)) {
      if (!CreateFace(v1, v3, v2)) {
        // Note: these vertices are allocated separately, and so they will not be deleted (memory leak)
        R3MeshVertex *v1a = CreateVertex(VertexPosition(v1));
        R3MeshVertex *v2a = CreateVertex(VertexPosition(v2));
        R3MeshVertex *v3a = CreateVertex(VertexPosition(v3));
        CreateFace(v1a, v2a, v3a);
      }
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}    



int R3Mesh::
ReadHoppeFile(const char *filename)
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "r"))) {
    fprintf(stderr, "Unable to open file %s\n", filename);
    return 0;
  }

  // Read file
  int line_count = 0;
  char buffer[1024];
  RNArray<R3MeshVertex *> vertices;
  RNArray<R3MeshVertex *> degenerate_triangle_vertices;
  while (fgets(buffer, 1023, fp)) {
    // Increment line counter
    line_count++;

    // Skip white space
    char *bufferp = buffer;
    while (isspace(*bufferp)) bufferp++;

    // Skip blank lines and comments
    if (*bufferp == '#') continue;
    if (*bufferp == '\0') continue;

    // Read command
    char cmd[1024];
    if (sscanf(buffer, "%s", cmd) != (unsigned int) 1) {
      fprintf(fp, "Syntax error line %d in file %s\n", line_count, filename);
      fclose(fp);
      return 0;
    }

    // Check command
    if (!strcmp(cmd, "Vertex")) {
      // Read vertex info
      int index;
      double x, y, z;
      if (sscanf(bufferp, "%s%d%lf%lf%lf", cmd, &index, &x, &y, &z) != (unsigned int) 5) {
        RNFail("Syntax error with vertex coordinates on line %d in file %s\n", line_count, filename);
        fclose(fp);
        return 0;
      }

      // Create vertex
      R3MeshVertex *vertex = CreateVertex(R3Point(x, y, z));
      vertices.Insert(vertex);
    }
    else if (!strcmp(cmd, "Face")) {
      // Read face info
      int index, i0, i1, i2;
      if (sscanf(bufferp, "%s%d%d%d%d", cmd, &index, &i0, &i1, &i2) != (unsigned int) 5) {
        RNFail("Syntax error with face info on line %d in file %s\n", line_count, filename);
        fclose(fp);
        return 0;
      }

      // Find vertices
      assert((i0 >= 0) && (i0 < vertices.NEntries()));
      assert((i1 >= 0) && (i1 < vertices.NEntries()));
      assert((i2 >= 0) && (i2 < vertices.NEntries()));
      R3MeshVertex *v0 = vertices.Kth(i0-1);
      R3MeshVertex *v1 = vertices.Kth(i1-1);
      R3MeshVertex *v2 = vertices.Kth(i2-1);

      // Create triangle
      if ((v0 != v1) && (v0 != v2) && (v1 != v2)) {
        if (!CreateFace(v0, v1, v2)) {
          // Must have been degeneracy (e.g., flips or three faces sharing an edge)
          // Remember for later processing (to preserve vertex indices)
          degenerate_triangle_vertices.Insert(v0);
          degenerate_triangle_vertices.Insert(v1);
          degenerate_triangle_vertices.Insert(v2);
        }
      }
    }
  }

  // Create degenerate triangles
  for (int i = 0; i <= degenerate_triangle_vertices.NEntries()-3; i+=3) {
    R3MeshVertex *v1 = degenerate_triangle_vertices.Kth(i+0);
    R3MeshVertex *v2 = degenerate_triangle_vertices.Kth(i+1);
    R3MeshVertex *v3 = degenerate_triangle_vertices.Kth(i+2);
    if (!CreateFace(v1, v2, v3)) {
      if (!CreateFace(v1, v3, v2)) {
        // Note: these vertices are allocated separately, and so they will not be deleted (memory leak)
        R3MeshVertex *v1a = CreateVertex(VertexPosition(v1));
        R3MeshVertex *v2a = CreateVertex(VertexPosition(v2));
        R3MeshVertex *v3a = CreateVertex(VertexPosition(v3));
        CreateFace(v1a, v2a, v3a);
      }
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



static int 
ReadIfsString(FILE *fp, char *buffer, int maxlength)
{
  // Read length
  unsigned int length;
  if (fread(&length, sizeof(unsigned int), 1, fp) != 1) {
    RNFail("Unable to read string length");
    return -1;
  }

  // Read characters
  if ((int) length > maxlength) length = maxlength;
  if (fread(buffer, sizeof(unsigned char), length, fp) != length) {
    RNFail("Unable to read string characters");
    return -1;
  }

  // Return length of string
  return length;
}



int R3Mesh::
ReadIfsFile(const char *filename)
{
  char buffer[1024];

  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "rb"))) {
    fprintf(stderr, "Unable to open file %s\n", filename);
    return 0;
  }

  // Read magic string
  if (ReadIfsString(fp, buffer, 1024) < 0) {
    RNFail("Unable to read header of %s", filename);
    return 0;
  }
  if (strcmp(buffer, "IFS")) {
    RNFail("Bad magic string in IFS file header of %s", filename);
    return 0;
  }

  // Read version number
  float version;
  if (fread(&version, sizeof(float), 1, fp) != 1) {
    RNFail("Unable to read version of %s", filename);
    return 0;
  }
  if (version != 1.0) {
    RNFail("Bad version number in file header of %s", filename);
    return 0;
  }

  // Read model name
  if (ReadIfsString(fp, buffer, 1024) < 0) {
    RNFail("Unable to read model name of %s", filename);
    return 0;
  }

  // Read vertex header
  if (ReadIfsString(fp, buffer, 1024) < 0) {
    RNFail("Unable to read vertex header of %s", filename);
    return 0;
  }
  if (strcmp(buffer, "VERTICES")) {
    RNFail("Bad vertex header in %s", filename);
    return 0;
  }

  // Read number of vertices
  unsigned int nverts;
  if (fread(&nverts, sizeof(unsigned int), 1, fp) != 1) {
    RNFail("Unable to read number of vertices in %s", filename);
    return 0;
  }

  // Allocate block of vertices
  vertex_block = new R3MeshVertex [nverts];
  
  // Resize array of vertices
  vertices.Resize(nverts);

  // Read vertices
  for (unsigned int i = 0; i < nverts; i++) {
    float p[3];
    if (fread(p, sizeof(float), 3, fp) != 3) {
      RNFail("Unable to read vertex %d in %s", i, filename);
      return 0;
    }

    // Create mesh vertex
    if (!CreateVertex(R3Point(p[0], p[1], p[2]), &vertex_block[i])) {
      RNFail("Unable to create vertex %d in %s", i, filename);
      return 0;
    }
  }

  // Read triangle header
  if (ReadIfsString(fp, buffer, 1024) < 0) {
    RNFail("Unable to read vertex header of %s", filename);
    return 0;
  }
  if (strcmp(buffer, "TRIANGLES")) {
    RNFail("Bad triangle header in %s", filename);
    return 0;
  }

  // Read number of faces
  unsigned int nfaces;
  if (fread(&nfaces, sizeof(unsigned int), 1, fp) != 1) {
    RNFail("Unable to read number of faces in %s", filename);
    return 0;
  }

  // Allocate block of faces
  face_block = new R3MeshFace [nfaces];

  // Resize array of faces
  faces.Resize(nfaces);

  // Allocate array for degenerate triangles
  RNArray<R3MeshVertex *> degenerate_triangle_vertices;

  // Read triangles
  for (unsigned int i = 0; i < nfaces; i++) {
    int v[3];
    if (fread(v, sizeof(unsigned int), 3, fp) != 3) {
      RNFail("Unable to read vertex index in triangle %d of %s", i, filename);
      return 0;
    }

    // Get vertices
    R3MeshVertex *v0 = vertices[v[0]];
    R3MeshVertex *v1 = vertices[v[1]];
    R3MeshVertex *v2 = vertices[v[2]];

    // Check vertices
    if ((v0 == v1) || (v1 == v2) || (v0 == v2)) continue;

    // Create mesh face
    if (!CreateFace(v0, v1, v2, &face_block[i])) {
      // Must have been degeneracy (e.g., flips or three faces sharing an edge)
      // Remember for later processing (to preserve vertex indices)
      degenerate_triangle_vertices.Insert(v0);
      degenerate_triangle_vertices.Insert(v1);
      degenerate_triangle_vertices.Insert(v2);
    }
  }

  // Create degenerate triangles
  for (int i = 0; i <= degenerate_triangle_vertices.NEntries()-3; i+=3) {
    R3MeshVertex *v1 = degenerate_triangle_vertices.Kth(i+0);
    R3MeshVertex *v2 = degenerate_triangle_vertices.Kth(i+1);
    R3MeshVertex *v3 = degenerate_triangle_vertices.Kth(i+2);
    if (!CreateFace(v1, v2, v3)) {
      if (!CreateFace(v1, v3, v2)) {
        // Note: these vertices are allocated separately, and so they will not be deleted (memory leak)
        R3MeshVertex *v1a = CreateVertex(VertexPosition(v1));
        R3MeshVertex *v2a = CreateVertex(VertexPosition(v2));
        R3MeshVertex *v3a = CreateVertex(VertexPosition(v3));
        CreateFace(v1a, v2a, v3a);
      }
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R3Mesh::
ReadSTLFile(const char *filename)
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "r"))) {
    RNFail("Unable to open file %s", filename);
    return 0;
  }

  // Read body
  int line_number = 0;
  char buffer[1024], cmd[1024];
  RNArray<R3MeshVertex *> vertices;
  RNArray<R3MeshVertex *> degenerate_triangle_vertices;
  while (fgets(buffer, 1024, fp)) {
    line_number++;
    if (sscanf(buffer, "%s", cmd) != 1) continue;
    if (!strcmp(cmd, "vertex")) {
      RNCoord x, y, z;
      if (sscanf(buffer, "%s%lf%lf%lf", cmd, &x, &y, &z) != 4) {
        RNFail("Error in vertex of STL file %s on line %d", filename, line_number);
        return 0;
      }

      // Make vertex
      const R3Point position(x, y, z);
      R3MeshVertex *vertex = CreateVertex(position);
      vertices.Insert(vertex);
    }
    else if (!strcmp(cmd, "outer")) {
      // Just checking
      assert(vertices.IsEmpty());
    }
    else if (!strcmp(cmd, "endloop")) {
      // Create face(s)
      for (int i = 2; i < vertices.NEntries(); i++) {
        if (vertices[0] == vertices[i-1]) continue;
        if (vertices[i-1] == vertices[i]) continue;
        if (vertices[0] == vertices[i]) continue;
        if (!CreateFace(vertices[0], vertices[i-1], vertices[i])) {
          // Must have been degeneracy (e.g., flips or three faces sharing an edge)
          // Remember for later processing (to preserve vertex indices)
          degenerate_triangle_vertices.Insert(vertices[0]);
          degenerate_triangle_vertices.Insert(vertices[i-1]);
          degenerate_triangle_vertices.Insert(vertices[i]);
        }
        vertices.Empty();
      }
    }
  }

  // Create degenerate triangles
  for (int i = 0; i <= degenerate_triangle_vertices.NEntries()-3; i+=3) {
    R3MeshVertex *v1 = degenerate_triangle_vertices.Kth(i+0);
    R3MeshVertex *v2 = degenerate_triangle_vertices.Kth(i+1);
    R3MeshVertex *v3 = degenerate_triangle_vertices.Kth(i+2);
    if (!CreateFace(v1, v2, v3)) {
      if (!CreateFace(v1, v3, v2)) {
        // Note: these vertices are allocated separately, and so they will not be deleted (memory leak)
        R3MeshVertex *v1a = CreateVertex(VertexPosition(v1));
        R3MeshVertex *v2a = CreateVertex(VertexPosition(v2));
        R3MeshVertex *v3a = CreateVertex(VertexPosition(v3));
        CreateFace(v1a, v2a, v3a);
      }
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R3Mesh::
ReadVRMLFile(const char *filename)
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "r"))) {
    fprintf(stderr, "Unable to open file %s\n", filename);
    return 0;
  }

  // Read file
  char buffer[1024];
  int line_count = 0;
  RNBoolean coordinate3_active = FALSE;
  RNBoolean coordinate3_point_active = FALSE;
  RNBoolean indexedfaceset_active = FALSE;
  RNBoolean indexedfaceset_coordindex_active = FALSE;
  RNArray<R3MeshVertex *> vertices;
  RNArray<R3MeshVertex *> degenerate_triangle_vertices;
  while (fgets(buffer, 1023, fp)) {
    // Increment line counter
    line_count++;

    // Skip white space
    char *bufferp = buffer;
    while (isspace(*bufferp)) bufferp++;

    // Skip blank lines and comments
    if (*bufferp == '#') continue;
    if (*bufferp == '\0') continue;

    // Get first token
    bufferp = strtok(bufferp, " ,");
    if (!bufferp) break;

    // Check node being parsed
    if (coordinate3_point_active) {
      if (*bufferp == ']') {
        coordinate3_point_active = FALSE;
      } 
      else {
        do {
          if (*bufferp && !isspace(*bufferp)) {
            // Read next coordinate
            static int coordinate_index = 0;
            static RNScalar coordinates[3];
            coordinates[coordinate_index] = atof(bufferp);
            coordinate_index++;

            // Create vertex if have three coordinates
            if (coordinate_index == 3) {
              R3Point position(coordinates[0], coordinates[1], coordinates[2]);
              R3MeshVertex *vertex = CreateVertex(position);
              vertices.Insert(vertex);
              coordinate_index = 0;
            }
          }
          bufferp = strtok(NULL, " ,");
        } while (bufferp);
      }
    }
    else if (indexedfaceset_coordindex_active) {
      if (*bufferp == ']') {
        indexedfaceset_active = FALSE;
        indexedfaceset_coordindex_active = FALSE;
      } 
      else {
        do {
          // Read next coordinate index
          if (*bufferp && !isspace(*bufferp)) {
            static RNArray<R3MeshVertex *> face_vertices;
            int vertex_index = atoi(bufferp);
            if (vertex_index >= 0) {
              R3MeshVertex *vertex = vertices[vertex_index];
              face_vertices.Insert(vertex);
            }
            else {
              // Create face(s)
              R3MeshVertex *v1 = face_vertices[0];
              for (int k = 2; k < face_vertices.NEntries(); k++) {
                // Get vertices
                R3MeshVertex *v2 = face_vertices[k-1];
                R3MeshVertex *v3 = face_vertices[k];

                // Check vertices
                if ((v1 == v2) || (v2 == v3) || (v1 == v3)) continue;
 
                // Create face
                if (!CreateFace(v1, v2, v3)) {
                  // Must have been degeneracy (e.g., flips or three faces sharing an edge)
                  // Remember for later processing (to preserve vertex indices)
                  degenerate_triangle_vertices.Insert(v1);
                  degenerate_triangle_vertices.Insert(v2);
                  degenerate_triangle_vertices.Insert(v3);
                }
              }
              face_vertices.Empty();
            }
          }
          bufferp = strtok(NULL, " ,");
        } while (bufferp);
      }
    }
    else {
      if (!strcmp(bufferp, "Coordinate3")) coordinate3_active = TRUE;
      else if (!strcmp(bufferp, "IndexedFaceSet")) indexedfaceset_active = TRUE;
      else if ((coordinate3_active) && !strcmp(bufferp, "point")) coordinate3_point_active = TRUE;
      else if ((indexedfaceset_active) && !strcmp(bufferp, "coordIndex")) indexedfaceset_coordindex_active = TRUE;
    }
  }

  // Create degenerate triangles
  for (int i = 0; i <= degenerate_triangle_vertices.NEntries()-3; i+=3) {
    R3MeshVertex *v1 = degenerate_triangle_vertices.Kth(i+0);
    R3MeshVertex *v2 = degenerate_triangle_vertices.Kth(i+1);
    R3MeshVertex *v3 = degenerate_triangle_vertices.Kth(i+2);
    if (!CreateFace(v1, v2, v3)) {
      if (!CreateFace(v1, v3, v2)) {
        // Note: these vertices are allocated separately, and so they will not be deleted (memory leak)
        R3MeshVertex *v1a = CreateVertex(VertexPosition(v1));
        R3MeshVertex *v2a = CreateVertex(VertexPosition(v2));
        R3MeshVertex *v3a = CreateVertex(VertexPosition(v3));
        CreateFace(v1a, v2a, v3a);
      }
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R3Mesh::
WriteFile(const char *filename) const
{
  // Parse input filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    printf("Filename %s has no extension (e.g., .ply)", filename);
    return 0;
  }

  // Write file of appropriate type
  if (!strncmp(extension, ".ray", 4)) 
    return WriteRayFile(filename);
  else if (!strncmp(extension, ".ply", 4)) 
    return WritePlyFile(filename);
  else if (!strncmp(extension, ".obj", 4)) 
    return WriteObjFile(filename);
  else if (!strncmp(extension, ".off", 4)) 
    return WriteOffFile(filename);
  else if (!strncmp(extension, ".cat", 4)) 
    return WriteCattFile(filename);
  else if (!strncmp(extension, ".ifs", 4)) 
    return WriteIfsFile(filename);
  else if (!strncmp(extension, ".stl", 4)) 
    return WriteSTLFile(filename);
  else {
    RNFail("Unable to write file %s (unrecognized extension: %s)", filename, extension);
    return 0;
  }
}



int R3Mesh::
WriteRayFile(const char *filename) const
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "w"))) {
    RNFail("Unable to open file %s", filename);
    return 0;
  }

  // Write vertices
  for (int i = 0; i < NVertices(); i++) {
    R3MeshVertex *vertex = Vertex(i);
    const R3Point& p = VertexPosition(vertex);
    const R3Vector& n = VertexNormal(vertex);
    fprintf(fp, "#vertex %g %g %g %g %g %g 0 0\n", p.X(), p.Y(), p.Z(), n.X(), n.Y(), n.Z());
  }

  // Write Faces
  for (int i = 0; i < NFaces(); i++) {
    R3MeshFace *face = Face(i);
    R3MeshVertex *v0 = VertexOnFace(face, 0);
    R3MeshVertex *v1 = VertexOnFace(face, 1);
    R3MeshVertex *v2 = VertexOnFace(face, 2);
    int m = FaceMaterial(face);
    fprintf(fp, "#shape_triangle %d %d %d %d\n", m, v0->id, v1->id, v2->id);
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R3Mesh::
WritePlyFile(const char *filename, RNBoolean binary)  const
{
  typedef struct PlyVertex {
    float x, y, z;
    float nx, ny, nz;
    float tx, ty;
    unsigned char red, green, blue;
  } PlyVertex;

  typedef struct PlyFace {
    unsigned char nverts;
    int *verts;
    int material;
    int segment;
    int category;
  } PlyFace;

  // Element names
  char *elem_names[] = { (char *) "vertex", (char *) "face" };

  // List of property information for a vertex 
  static PlyProperty vert_props[] = { 
    {(char *) "x", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex,x), 0, 0, 0, 0},
    {(char *) "y", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex,y), 0, 0, 0, 0},
    {(char *) "z", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex,z), 0, 0, 0, 0},
    {(char *) "nx", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex,nx), 0, 0, 0, 0},
    {(char *) "ny", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex,ny), 0, 0, 0, 0},
    {(char *) "nz", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex,nz), 0, 0, 0, 0},
    {(char *) "tx", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex,tx), 0, 0, 0, 0},
    {(char *) "ty", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex,ty), 0, 0, 0, 0},
    {(char *) "red", PLY_UCHAR, PLY_UCHAR, offsetof(PlyVertex,red), 0, 0, 0, 0},
    {(char *) "green", PLY_UCHAR, PLY_UCHAR, offsetof(PlyVertex,green), 0, 0, 0, 0},
    {(char *) "blue", PLY_UCHAR, PLY_UCHAR, offsetof(PlyVertex,blue), 0, 0, 0, 0}
  };

  // List of property information for a vertex 
  static PlyProperty face_props[] = { 
    {(char *) "vertex_indices", PLY_INT, PLY_INT, offsetof(PlyFace,verts), 1, PLY_UCHAR, PLY_UCHAR, offsetof(PlyFace,nverts)},
    {(char *) "material_id", PLY_INT, PLY_INT, offsetof(PlyFace,material), 0, 0, 0, 0},
    {(char *) "segment_id", PLY_INT, PLY_INT, offsetof(PlyFace,segment), 0, 0, 0, 0},
    {(char *) "category_id", PLY_INT, PLY_INT, offsetof(PlyFace,category), 0, 0, 0, 0}
  };

  // Open ply file
  float version;
  int file_type = (binary) ? PLY_BINARY_NATIVE : PLY_ASCII;
  PlyFile *ply = ply_open_for_writing((char *) filename, 2, elem_names, file_type, &version);
  if (!ply) return -1;

  // Describe vertex properties
  ply_element_count(ply, (char *) "vertex", NVertices());
  ply_describe_property(ply, (char *) "vertex", &vert_props[0]);
  ply_describe_property(ply, (char *) "vertex", &vert_props[1]);
  ply_describe_property(ply, (char *) "vertex", &vert_props[2]);
  ply_describe_property(ply, (char *) "vertex", &vert_props[3]);
  ply_describe_property(ply, (char *) "vertex", &vert_props[4]);
  ply_describe_property(ply, (char *) "vertex", &vert_props[5]);
  ply_describe_property(ply, (char *) "vertex", &vert_props[6]);
  ply_describe_property(ply, (char *) "vertex", &vert_props[7]);
  ply_describe_property(ply, (char *) "vertex", &vert_props[8]);
  ply_describe_property(ply, (char *) "vertex", &vert_props[9]);
  ply_describe_property(ply, (char *) "vertex", &vert_props[10]);

  // Describe face properties
  ply_element_count(ply, (char *) "face", NFaces());
  ply_describe_property(ply, (char *) "face", &face_props[0]);
  ply_describe_property(ply, (char *) "face", &face_props[1]);
  ply_describe_property(ply, (char *) "face", &face_props[2]);
  ply_describe_property(ply, (char *) "face", &face_props[3]);

  // Complete header
  ply_header_complete(ply);

  // Write vertices
  ply_put_element_setup(ply, (char *) "vertex");
  for (int i = 0; i < NVertices(); i++) {
    R3MeshVertex *v = Vertex(i);
    const R3Point& p = VertexPosition(v);
    const R3Vector& n = VertexNormal(v);
    const R2Point& t = VertexTextureCoords(v);
    const RNRgb& c = VertexColor(v);
    PlyVertex ply_vertex;
    ply_vertex.x = p.X();
    ply_vertex.y = p.Y();
    ply_vertex.z = p.Z();
    ply_vertex.nx = n.X();
    ply_vertex.ny = n.Y();
    ply_vertex.nz = n.Z();
    ply_vertex.tx = t.X();
    ply_vertex.ty = t.Y();
    ply_vertex.red = n.X();
    ply_vertex.red = (unsigned char) (255.0 * c.R());
    ply_vertex.green = (unsigned char) (255.0 * c.G());
    ply_vertex.blue = (unsigned char) (255.0 * c.B());
    ply_put_element(ply, (void *) &ply_vertex);
  }

  // Write faces
  ply_put_element_setup(ply, (char *) "face");
  for (int i = 0; i < NFaces(); i++) {
    R3MeshFace *f = Face(i);
    static int verts[3];
    static PlyFace ply_face = { 3, verts, 0 };
    ply_face.verts[0] = VertexID(VertexOnFace(f, 0));
    ply_face.verts[1] = VertexID(VertexOnFace(f, 1));
    ply_face.verts[2] = VertexID(VertexOnFace(f, 2));
    ply_face.material = FaceMaterial(f);
    ply_face.segment = FaceSegment(f);
    ply_face.category = FaceCategory(f);
    ply_put_element(ply, (void *) &ply_face);
  }

  // Close the file 
  ply_close(ply);

  // Return number of faces written
  return 1;
}



int R3Mesh::
WriteObjFile(const char *filename) const
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "w"))) {
    RNFail("Unable to open file %s", filename);
    return 0;
  }

  // Check if have texcoords
  RNBoolean has_texcoords = FALSE;
  for (int i = 0; i < NVertices(); i++) {
    R3MeshVertex *vertex = Vertex(i);
    const R2Point& t = VertexTextureCoords(vertex);
    if (!t.IsZero()) { has_texcoords = TRUE; break; }
  }

  // Write vertices
  for (int i = 0; i < NVertices(); i++) {
    R3MeshVertex *vertex = Vertex(i);
    const R3Point& p = VertexPosition(vertex);
    const R2Point& t = VertexTextureCoords(vertex);
    if (has_texcoords) fprintf(fp, "vt %g %g\n", t.X(), t.Y());
    fprintf(fp, "v %g %g %g\n", p.X(), p.Y(), p.Z());
  }

  // Write Faces
  for (int i = 0; i < NFaces(); i++) {
    R3MeshFace *face = Face(i);
    R3MeshVertex *v0 = VertexOnFace(face, 0);
    R3MeshVertex *v1 = VertexOnFace(face, 1);
    R3MeshVertex *v2 = VertexOnFace(face, 2);
    if (has_texcoords) {
      fprintf(fp, "f %d/%d %d/%d %d/%d\n",
        v0->id + 1, v1->id + 1, v2->id + 1,
        v0->id + 1, v1->id + 1, v2->id + 1);
    }
    else {
      fprintf(fp, "f %d %d %d\n",
        v0->id + 1, v1->id + 1, v2->id + 1);
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R3Mesh::
WriteOffFile(const char *filename) const
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "w"))) {
    RNFail("Unable to open file %s", filename);
    return 0;
  }

  // Write header
  fprintf(fp, "OFF\n");
  fprintf(fp, "%d %d %d\n", NVertices(), NFaces(), NEdges());

  // Write vertices
  for (int i = 0; i < NVertices(); i++) {
    R3MeshVertex *vertex = Vertex(i);
    const R3Point& p = VertexPosition(vertex);
    fprintf(fp, "%g %g %g\n", p.X(), p.Y(), p.Z());
  }

  // Write Faces
  for (int i = 0; i < NFaces(); i++) {
    R3MeshFace *face = Face(i);
    R3MeshVertex *v0 = VertexOnFace(face, 0);
    R3MeshVertex *v1 = VertexOnFace(face, 1);
    R3MeshVertex *v2 = VertexOnFace(face, 2);
    fprintf(fp, "3 %d %d %d\n", v0->id, v1->id, v2->id);
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R3Mesh::
WriteCattFile(const char *filename) const
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "w"))) {
    RNFail("Unable to open file %s", filename);
    return 0;
  }

  // Write vertex header
  fprintf(fp, "%%CORNERS\n\n");

  // Write vertices
  for (int i = 0; i < NVertices(); i++) {
    R3MeshVertex *vertex = Vertex(i);
    int id = VertexID(vertex);
    const R3Point& p = VertexPosition(vertex);
    fprintf(fp, "%d %g %g %g\n", id+1, p.X(), p.Y(), p.Z());
  }

  // Write face header
  fprintf(fp, "\n%%PLANES\n\n");

  // Write Faces
  for (int i = 0; i < NFaces(); i++) {
    R3MeshFace *face = Face(i);
    R3MeshVertex *v0 = VertexOnFace(face, 0);
    R3MeshVertex *v1 = VertexOnFace(face, 1);
    R3MeshVertex *v2 = VertexOnFace(face, 2);
    fprintf(fp, "%d / /RIGID\n%d %d %d\n\n", i, VertexID(v0) + 1, VertexID(v1) + 1, VertexID(v2) + 1);
  }

  // Write trailer
  fprintf(fp, "%%EOF\n");

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



static int 
WriteIfsString(FILE *fp, const char *buffer)
{
  // Write length
  unsigned int length = strlen(buffer) + 1;
  if (fwrite(&length, sizeof(unsigned int), 1, fp) != 1) {
    RNFail("Unable to write string length");
    return -1;
  }

  // Write characters
  if (fwrite(buffer, sizeof(unsigned char), length, fp) != length) {
    RNFail("Unable to write string characters");
    return -1;
  }

  // Return length of string
  return length;
}



int R3Mesh::
WriteIfsFile(const char *filename) const
{
  float version = 1.0;
  unsigned int nverts = NVertices();
  unsigned int nfaces = NFaces();

  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "wb"))) {
    RNFail("Unable to open file %s", filename);
    return 0;
  }

  // Write header
  WriteIfsString(fp, "IFS");
  fwrite(&version, sizeof(float), 1, fp);
  WriteIfsString(fp, "filname");

  // Write vertices
  WriteIfsString(fp, "VERTICES");
  fwrite(&nverts, sizeof(unsigned int), 1, fp);
  for (unsigned int i = 0; i < nverts; i++) {
    R3MeshVertex *vertex = Vertex(i);
    const R3Point& position = VertexPosition(vertex);
    float p[3];
    p[0] = position.X();
    p[1] = position.Y();
    p[2] = position.Z();
    fwrite(&p, sizeof(float), 3, fp);
  }

  // Write faces
  WriteIfsString(fp, "TRIANGLES");
  fwrite(&nfaces, sizeof(unsigned int), 1, fp);
  for (int i = 0; i < NFaces(); i++) {
    R3MeshFace *face = Face(i);
    R3MeshVertex *v0 = VertexOnFace(face, 0);
    R3MeshVertex *v1 = VertexOnFace(face, 1);
    R3MeshVertex *v2 = VertexOnFace(face, 2);
    unsigned int v[3];
    v[0] = v0->id;
    v[1] = v1->id;
    v[2] = v2->id;
    fwrite(&v, sizeof(unsigned int), 3, fp);
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R3Mesh::
WriteSTLFile(const char *filename) const
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "w"))) {
    RNFail("Unable to open file %s", filename);
    return 0;
  }

  // Write header
  fprintf(fp, "solid %s\n", filename);

  // Write Faces
  for (int i = 0; i < NFaces(); i++) {
    R3MeshFace *face = Face(i);
    const R3Vector& n = FaceNormal(face);
    fprintf(fp, "facet normal %g %g %g\n", n[0], n[1], n[2]);
    fprintf(fp, "  outer loop\n");
    for (int j = 0; j < 3; j++) {
      R3MeshVertex *vertex = VertexOnFace(face, j);
      const R3Point& p = VertexPosition(vertex);
      fprintf(fp, "   vertex %g %g %g\n", p[0], p[1], p[2]);
    }
    fprintf(fp, "  endloop\n");
    fprintf(fp, "endfacet\n");
  }

  // Write trailer
  fprintf(fp, "endsolid\n");

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// USEFUL DEBUG FUNCTIONS
////////////////////////////////////////////////////////////////////////

RNBoolean R3Mesh::
IsValid(R3MeshVertex *v) const
{
#ifndef NDEBUG
  // Check if vertex is on mesh
  assert(IsVertexOnMesh(v));

  // Check each vertex-edge and vertex-face relation
  for (int j = 0; j < VertexValence(v); j++) {
    R3MeshEdge *e = EdgeOnVertex(v, j);
    assert(IsVertexOnEdge(v, e));
    assert((!e->face[0]) || IsVertexOnFace(v, e->face[0]));
    assert((!e->face[1]) || IsVertexOnFace(v, e->face[1]));
  }
#endif

  // Return OK status
  return TRUE;
}



RNBoolean R3Mesh::
IsValid(R3MeshEdge *e) const
{
#ifndef NDEBUG
  // Check if edge is on mesh
  assert(IsEdgeOnMesh(e));

  // Check vertices
  assert(e->vertex[0]);
  assert(e->vertex[1]);
  assert(e->vertex[0] != e->vertex[1]);

  // Check each edge-vertex relation
  for (int j = 0; j < 2; j++) {
    R3MeshVertex *v = e->vertex[j];
    assert(v && v->edges.FindEntry(e));
  }

  // Check each edge-face relation
  for (int j = 0; j < 2; j++) {
    R3MeshFace *f = e->face[j];
    assert(!f || (f->edge[0] == e) || (f->edge[1] == e) || (f->edge[2] == e));
  }
#endif

  // Return OK status
  return TRUE;
}



RNBoolean R3Mesh::
IsValid(R3MeshFace *f) const
{
#ifndef NDEBUG
  // Check if face is on mesh
  assert(IsFaceOnMesh(f));

  // Check vertices
  assert(f->vertex[0]);
  assert(f->vertex[1]);
  assert(f->vertex[1]);
  assert(f->vertex[0] != f->vertex[1]);
  assert(f->vertex[1] != f->vertex[2]);
  assert(f->vertex[0] != f->vertex[2]);

  // Check edges
  assert(f->edge[0]);
  assert(f->edge[1]);
  assert(f->edge[1]);
  assert(f->edge[0] != f->edge[1]);
  assert(f->edge[1] != f->edge[2]);
  assert(f->edge[0] != f->edge[2]);

  // Check each face-vertex and face-edge relation
  for (int j = 0; j < 3; j++) {
    R3MeshVertex *v = f->vertex[j];
    R3MeshEdge *e = f->edge[j];
    assert(v && ((e->vertex[0] == v) || (e->vertex[1] == v)));
    assert(e && ((e->face[0] == f) || (e->face[1] == f)));
    assert(v == VertexOnFace(f, e, RN_CW));
    assert(e == EdgeOnFace(f, v, RN_CCW));
  }
#endif

  // Return OK status
  return TRUE;
}



RNBoolean R3Mesh::
IsValid(void) const
{
#ifndef NDEBUG
  // Check each vertex, edge, and face
  for (int i = 0; i < NVertices(); i++) assert(IsValid(Vertex(i)));
  for (int i = 0; i < NEdges(); i++) assert(IsValid(Edge(i)));
  for (int i = 0; i < NFaces(); i++) assert(IsValid(Face(i)));
#endif

  // Return success
  return TRUE;
}



////////////////////////////////////////////////////////////////////////
// INTERNAL UPDATE FUNCTIONS
////////////////////////////////////////////////////////////////////////

void R3Mesh::
UpdateVertexNormal(R3MeshVertex *v) const
{
  // Average face normals
  v->normal = R3zero_vector;
  for (int i = 0; i < v->edges.NEntries(); i++) {
    R3MeshEdge *e = v->edges[i];
    if (e->vertex[0] == v) {
      if (e->face[0]) 
        v->normal += FaceNormal(e->face[0]);
    }
    else {
      if (e->face[1]) 
        v->normal += FaceNormal(e->face[1]);
    }
  }
  v->normal.Normalize();

  // Update flags
  v->flags.Add(R3_MESH_VERTEX_NORMAL_UPTODATE);
}



void R3Mesh::
UpdateVertexCurvature(R3MeshVertex *v) const
{
  // Compute and remember Gauss curvature
  v->curvature = VertexGaussCurvature(v);

  // Update flags
  v->flags.Add(R3_MESH_VERTEX_CURVATURE_UPTODATE);
}



void R3Mesh::
UpdateEdgeLength(R3MeshEdge *e) const
{
  // Reset length
  e->length = R3Distance(e->vertex[0]->position, e->vertex[1]->position);

  // Update flags
  e->flags.Add(R3_MESH_EDGE_LENGTH_UPTODATE);
}



void R3Mesh::
UpdateFaceArea(R3MeshFace *f) const
{
  // Compute area of face
  R3Vector v1 = f->vertex[1]->position - f->vertex[0]->position;
  R3Vector v2 = f->vertex[2]->position - f->vertex[0]->position;
  R3Vector v3 = v1 % v2;
  f->area = 0.5 * v3.Length();

  // Update flags
  f->flags.Add(R3_MESH_FACE_AREA_UPTODATE);
}



void R3Mesh::
UpdateFacePlane(R3MeshFace *f) const
{
  // Reset plane 
  f->plane = R3Plane(f->vertex[0]->position, f->vertex[1]->position, f->vertex[2]->position);

  // Update flags
  f->flags.Add(R3_MESH_FACE_PLANE_UPTODATE);
}



void R3Mesh::
UpdateFaceBBox(R3MeshFace *f) const
{
  // Update face bbox
  f->bbox = R3null_box;
  f->bbox.Union(f->vertex[0]->position);
  f->bbox.Union(f->vertex[1]->position);
  f->bbox.Union(f->vertex[2]->position);

  // Update flags
  f->flags.Add(R3_MESH_FACE_BBOX_UPTODATE);
}



void R3Mesh::
UpdateFaceRefs(R3MeshFace *f,
               R3MeshVertex *v1, R3MeshVertex *v2, R3MeshVertex *v3,
               R3MeshEdge *e1, R3MeshEdge *e2, R3MeshEdge *e3)
{
  // Update face-vertex relations
  f->vertex[0] = v1;
  f->vertex[1] = v2;
  f->vertex[2] = v3;

  // Update face-edge relations
  f->edge[0] = e1;
  f->edge[1] = e2;
  f->edge[2] = e3;

  // Update edge-face relations
  if (e1->vertex[0] == v1) e1->face[0] = f;
  else { assert(e1->vertex[1] == v1); e1->face[1] = f; }
  if (e2->vertex[0] == v2) e2->face[0] = f;
  else { assert(e2->vertex[1] == v2); e2->face[1] = f; }
  if (e3->vertex[0] == v3) e3->face[0] = f;
  else { assert(e3->vertex[1] == v3); e3->face[1] = f; }

  // Invalidate face properties
  f->flags.Remove(R3_MESH_FACE_AREA_UPTODATE | R3_MESH_FACE_PLANE_UPTODATE | R3_MESH_FACE_BBOX_UPTODATE);

  // Invalidate vertex properties
  v1->flags.Remove(R3_MESH_VERTEX_NORMAL_UPTODATE | R3_MESH_VERTEX_CURVATURE_UPTODATE);
  v2->flags.Remove(R3_MESH_VERTEX_NORMAL_UPTODATE | R3_MESH_VERTEX_CURVATURE_UPTODATE);
  v3->flags.Remove(R3_MESH_VERTEX_NORMAL_UPTODATE | R3_MESH_VERTEX_CURVATURE_UPTODATE);
}
 
   

////////////////////////////////////////////////////////////////////////
// EUCLIDEAN DISTANCE FUNCTIONS
////////////////////////////////////////////////////////////////////////

RNLength
R3Distance(R3Mesh *mesh, R3MeshVertex *vertex, const R3Point& point, R3MeshIntersection* closest_point)
{
  // Compute squared distance from vertex to point
  const R3Point& vertex_position = mesh->VertexPosition(vertex);
  RNLength distance = R3Distance(point, vertex_position);

  // Fill in closest_point info
  if (closest_point) {
    closest_point->type = R3_MESH_VERTEX_TYPE;
    closest_point->vertex = vertex;
    closest_point->edge = NULL;
    closest_point->face = NULL;
    closest_point->point = vertex_position;
    closest_point->t = distance; 
  }

  // Return distance
  return distance;
}



RNLength 
R3Distance(R3Mesh *mesh, R3MeshEdge *edge, const R3Point& point, R3MeshIntersection* closest_point)
{
  // Compute distance from edge to point
  R3Span edge_span = mesh->EdgeSpan(edge);
  RNScalar distance = R3Distance(point, edge_span);

  // Fill in closest point info
  if (closest_point) {
    closest_point->edge = edge;
    closest_point->face = NULL;
    closest_point->t = distance; 
    RNScalar edge_t = edge_span.T(point);
    if (RNIsEqual(edge_t, 0)) {
      R3MeshVertex *vertex = mesh->VertexOnEdge(edge, 0);
      closest_point->type = R3_MESH_VERTEX_TYPE;
      closest_point->vertex = vertex;
      closest_point->point = mesh->VertexPosition(vertex);
    }
    else if (RNIsEqual(edge_t, edge_span.Length())) {
      R3MeshVertex *vertex = mesh->VertexOnEdge(edge, 1);
      closest_point->type = R3_MESH_VERTEX_TYPE;
      closest_point->vertex = vertex;
      closest_point->point = mesh->VertexPosition(vertex);
    }
    else {
      closest_point->type = R3_MESH_EDGE_TYPE;
      closest_point->vertex = NULL;
      closest_point->point = edge_span.Point(edge_t);
    }
  }

  // Return squared distance
  return distance;
}



RNLength
R3Distance(R3Mesh *mesh, R3MeshFace *face, const R3Point& point, R3MeshIntersection *closest_point) 
{
  // Compute projection of point onto face plane
  const R3Plane& plane = mesh->FacePlane(face);
  const R3Vector& face_normal = plane.Normal();
  RNScalar plane_signed_distance = R3SignedDistance(plane, point);
  R3Point plane_point = point - plane_signed_distance * face_normal;
  RNLength closest_distance = fabs(plane_signed_distance);

  // Check if point is outside any of the edges
  R3MeshEdge *closest_edge = NULL;
  for (int i = 0; i < 3; i++) {
    R3MeshEdge *edge = mesh->EdgeOnFace(face, i);
    R3MeshVertex *v0 = mesh->VertexOnEdge(edge, face, RN_CW);
    R3MeshVertex *v1 = mesh->VertexOnEdge(edge, face, RN_CCW);
    R3Point p0 = mesh->VertexPosition(v0);
    R3Point p1 = mesh->VertexPosition(v1);
    R3Vector edge_vector = p1 - p0;
    edge_vector.Normalize();
    R3Vector edge_normal = face_normal % edge_vector;
    R3Plane edge_plane(p0, edge_normal);
    RNScalar b = R3SignedDistance(edge_plane, plane_point);
    if (b < 0) {
      R3MeshIntersection intersection;
      RNLength distance = R3Distance(mesh, edge, point, &intersection);
      if (!closest_edge || (distance < closest_distance)) {
        if (closest_point) *closest_point = intersection;
        closest_distance = distance;
        closest_edge = edge;
      }
    }
  }

  // Fill in closest point info
  if (closest_point) {
    closest_point->face = face;
    if (!closest_edge) {
      closest_point->type = R3_MESH_FACE_TYPE;
      closest_point->vertex = NULL;
      closest_point->edge = NULL;
      closest_point->point = plane_point;
      closest_point->t = closest_distance;
    }
  }


  // Return distance
  return closest_distance;
}



////////////////////////////////////////////////////////////////////////
// CONSTRUCTORS FOR VERTEX, EDGE, FACE
////////////////////////////////////////////////////////////////////////

R3MeshVertex::
R3MeshVertex(void) 
  : position(0.0, 0.0, 0.0),
    normal(0.0, 0.0, 0.0),
    texcoords(0.0, 0.0),
    color(0.0, 0.0, 0.0),
    curvature(0),
    id(-1),
    flags(0),
    value(0.0),
    mark(0),
    data(NULL)
{
}



R3MeshVertex::
~R3MeshVertex(void) 
{
}



R3MeshEdge::
R3MeshEdge(void) 
  : length(0),
    id(-1),
    flags(0),
    value(0.0),
    mark(0),
    data(NULL)
{
  // Initialize edge fields
  vertex[0] = vertex[1] = NULL;
  face[0] = face[1] = NULL;
}



R3MeshEdge::
~R3MeshEdge(void) 
{
}



R3MeshFace::
R3MeshFace(void) 
  : plane(0.0, 0.0, 0.0, 0.0),
    bbox(1.0, 1.0, 1.0, -1.0, -1.0, -1.0),
    material(-1),
    segment(-1),
    category(-1),
    id(-1),
    flags(0),
    value(0.0),
    mark(0),
    data(NULL)
{
  // Initialize face fields
  vertex[0] = vertex[1] = vertex[2] = NULL;
  edge[0] = edge[1] = edge[2] = NULL;
}



R3MeshFace::
~R3MeshFace(void) 
{
}




