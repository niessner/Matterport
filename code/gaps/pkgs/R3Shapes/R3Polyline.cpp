/* Source file for the Catmull Rom Splineclass */



/* Include files */

#include "R3Shapes/R3Shapes.h"



R3Polyline::
R3Polyline(void)
  : vertex_positions(NULL),
    vertex_parameters(NULL),
    vertex_datas(NULL),
    nvertices(0),
    bbox(FLT_MAX, FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX)
{
}



R3Polyline::
R3Polyline(const R3Polyline& curve)
  : vertex_positions(NULL),
    vertex_parameters(NULL),
    vertex_datas(NULL),
    nvertices(curve.nvertices),
    bbox(curve.BBox())
{
  // Copy vertex_positions and parameters
  if (nvertices > 0) {
    this->vertex_positions = new R3Point [ nvertices ];
    this->vertex_parameters = new RNScalar [ nvertices ];
    for (int i = 0; i < nvertices; i++) {
      this->vertex_positions[i] = curve.vertex_positions[i];
      this->vertex_parameters[i] = curve.vertex_parameters[i];
    }
  }
}



R3Polyline::
R3Polyline(const RNArray<R3Point *>& points)
  : vertex_positions(NULL),
    vertex_parameters(NULL),
    vertex_datas(NULL),
    nvertices(points.NEntries()),
    bbox(FLT_MAX, FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX)
{
  // Copy vertex_positions
  if (points.NEntries() > 0) {
    this->vertex_positions = new R3Point [ points.NEntries() ];
    for (int i = 0; i < points.NEntries(); i++) {
      this->vertex_positions[i] = *(points[i]);
    }
  }
}



R3Polyline::
R3Polyline(const R3Point *points, int npoints)
  : vertex_positions(NULL),
    vertex_parameters(NULL),
    vertex_datas(NULL),
    nvertices(npoints),
    bbox(FLT_MAX, FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX)
{
  // Copy vertex_positions
  if (npoints > 0) {
    this->vertex_positions = new R3Point [ npoints ];
    for (int i = 0; i < npoints; i++) {
      this->vertex_positions[i] = points[i];
    }
  }
}



R3Polyline::
R3Polyline(const R3Point *points, const RNScalar *parameters, int npoints)
  : vertex_positions(NULL),
    vertex_parameters(NULL),
    vertex_datas(NULL),
    nvertices(npoints),
    bbox(FLT_MAX, FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX)
{
  // Copy vertex_positions and parameters
  if (npoints > 0) {
    this->vertex_positions = new R3Point [ npoints ];
    this->vertex_parameters = new RNScalar [ npoints ];
    for (int i = 0; i < npoints; i++) {
      this->vertex_positions[i] = points[i];
      this->vertex_parameters[i] = parameters[i];
    }
  }
}



R3Polyline::
~R3Polyline(void)
{
  // Delete stuff
  if (vertex_positions) delete [] vertex_positions;
  if (vertex_parameters) delete [] vertex_parameters;
  if (vertex_datas) delete [] vertex_datas;
}



const RNScalar R3Polyline::
StartParameter(void) const
{
  // Return parameter at start of spline
  return 0;
}



const RNScalar R3Polyline::
EndParameter(void) const
{
  // Return parameter at end of spline
  if (vertex_parameters) return vertex_parameters[nvertices-1];
  else if (nvertices > 0) return nvertices-1;
  else return 0;
}



const R3Point R3Polyline::
StartPosition(void) const
{
  // Return position at start of spline
  if (nvertices == 0) return R3zero_point;
  return vertex_positions[0];
}



const R3Point R3Polyline::
EndPosition(void) const
{
  // Return position at end of spline
  if (nvertices == 0) return R3zero_point;
  return vertex_positions[nvertices-1];
}



R3Vector R3Polyline::
VertexDerivative(int i) const
{
  // Return derivative of polyline at ith vertex
  if (nvertices < 2) return R3zero_vector;
  if (i == 0) {
    RNScalar denom = VertexParameter(1) - VertexParameter(0);
    if (RNIsZero(denom)) return R3zero_vector;
    else return (VertexPosition(1) - VertexPosition(0)) / denom;
  }
  else if (i == nvertices-1) {
    RNScalar denom = VertexParameter(nvertices-1) - VertexParameter(nvertices-2);
    if (RNIsZero(denom)) return R3zero_vector;
    else return (VertexPosition(nvertices-1) - VertexPosition(nvertices-2)) / denom;
  }
  else {
    assert(nvertices > 3);
    RNScalar denom = VertexParameter(i+1) - VertexParameter(i-1);
    if (RNIsZero(denom)) return R3zero_vector;
    else return (VertexPosition(i+1) - VertexPosition(i-1)) / denom;
  }

  // Should never get here
  return R3zero_vector;
}



R3Vector R3Polyline::
VertexDirection(int i) const
{
  // Return tangent direction of polyline at ith vertex
  R3Vector direction = VertexDerivative(i);
  direction.Normalize();
  return direction;
}



void *R3Polyline::
VertexData(int i) const
{
  // Return user data associated with vertex
  if (!vertex_datas) return NULL;
  else return vertex_datas[i];
}



const RNBoolean R3Polyline::
IsPoint(void) const
{
  // Return whether spline covers one point
  if (nvertices == 0) return FALSE;
  if (nvertices == 1) return TRUE;

  // Look for vertices at different positions
  for (int i = 1; i < nvertices; i++) {
    if (!R3Contains(vertex_positions[0], vertex_positions[i])) {
      return FALSE;
    }
  }

  // All vertices are at same position
  return TRUE;
}



const RNBoolean R3Polyline::
IsLinear(void) const
{
  // Return whether spline lies on one line
  if (nvertices == 0) return FALSE;
  if (nvertices <= 2) return TRUE;

  // Check if vertices all line on same line
  R3Line line;
  int line_npositions = 0;
  R3Point line_positions[2];
  for (int i = 0; i < nvertices; i++) {
    if (line_npositions == 0) {
      line_positions[line_npositions++] = vertex_positions[i];
    }
    else if ((line_npositions == 1) && !R3Contains(line_positions[0], vertex_positions[i])) {
      line_positions[line_npositions++] = vertex_positions[i];
      line.Reset(line_positions[0], line_positions[1] - line_positions[0]);
    }
    else if (!R3Contains(line, vertex_positions[i])) {
      return FALSE;
    }
  }

  // All vertices are on same line
  return TRUE;
}



const RNBoolean R3Polyline::
IsPlanar(void) const
{
  // Return whether spline lies on one line
  if (nvertices == 0) return FALSE;
  if (nvertices <= 3) return TRUE;

  // Check if vertices all line on same plane
  R3Plane plane;
  int plane_npositions = 0;
  R3Point plane_positions[3];
  for (int i = 0; i < nvertices; i++) {
    if (plane_npositions == 0) {
      plane_positions[plane_npositions++] = vertex_positions[i];
    }
    else if ((plane_npositions == 1) &&
             !R3Contains(plane_positions[0], vertex_positions[i])) {
      plane_positions[plane_npositions++] = vertex_positions[i];
    }
    else if ((plane_npositions == 2) &&
             !R3Contains(plane_positions[0], vertex_positions[i]) &&
             !R3Contains(plane_positions[1], vertex_positions[i])) {
      plane_positions[plane_npositions++] = vertex_positions[i];
      R3Vector v1 = plane_positions[1] - plane_positions[0];
      R3Vector v2 = plane_positions[2] - plane_positions[0];
      R3Vector normal = v1 % v2; normal.Normalize();
      plane.Reset(plane_positions[0], normal);
    }
    else if (!R3Contains(plane, vertex_positions[i])) {
      return FALSE;
    }
  }

  // All vertices are on same plane
  return TRUE;
}



const RNLength R3Polyline::
Length(void) const
{
  // For now, return length of control polygon
  RNLength length = 0;
  for (int i = 1; i < nvertices; i++) {
    length += R3Distance(vertex_positions[i-1], vertex_positions[i]);
  }
  return length;
}



const R3Point R3Polyline::
Centroid(void) const
{
  // Return centroid of vertices weighted by parameters
  RNScalar weight = 0;
  R3Point sum = R3zero_point;
  for (int i = 0; i < nvertices; i++) {
    RNScalar weight = 0;
    if (i > 0) weight += VertexParameter(i) - VertexParameter(i-1);
    if (i < nvertices-1) weight += VertexParameter(i+1) - VertexParameter(i);
    sum += vertex_positions[i];
  }

  // Return weighted average
  if (weight == 0) return R3zero_point;
  return sum / weight;
}



const R3Shape& R3Polyline::
BShape(void) const
{
  // Return bounding box
  return bbox;
}



const R3Box R3Polyline::
BBox(void) const
{
  // Return bounding box
  if (bbox[0][0] == FLT_MAX) 
    ((R3Polyline *) this)->UpdateBBox();
  return bbox;
}



static int
BinarySearchForVertexIndex(const R3Polyline *spline, int min_index, int max_index, RNScalar u)
{
  // Binary search for vertex index
  if (min_index == max_index) return min_index;
  if (spline->VertexParameter(min_index+1) > u) return min_index;
  if (spline->VertexParameter(max_index) < u) return max_index;
  int mid_index = (min_index + max_index)/2;
  RNScalar mid_u = spline->VertexParameter(mid_index);
  if (u == mid_u) return mid_index;
  else if (u < mid_u) return BinarySearchForVertexIndex(spline, min_index, mid_index-1, u);
  else return BinarySearchForVertexIndex(spline, mid_index, max_index, u);
}



int R3Polyline::
VertexIndex(RNScalar u) const
{
  // Return index at start of segment containing u
  if (u < 0) return 0;
  else if (nvertices < 3) return 0;
  else if (u >= VertexParameter(nvertices-2)) return nvertices-2;
  else if (!vertex_parameters) return (int) u;
  else return BinarySearchForVertexIndex(this, 0, nvertices-2, u);
}



R3Point R3Polyline::
PointPosition(RNScalar u) const
{
  // Return position at parametric value u

  // Get index variables
  int i = VertexIndex(u);
  double denom = VertexParameter(i+1) - VertexParameter(i);
  if (RNIsZero(denom)) return VertexPosition(i);
  double t = (u - VertexParameter(i)) / denom;

  // Get control vertex_positions
  R3Point P1 = VertexPosition(i);
  R3Point P2 = VertexPosition(i+1);

  // Return interpolated position
  return (1-t)*P1 + t*P2;
}



R3Vector R3Polyline::
PointDerivative(RNScalar u) const
{
  // Return tangent vector at parametric value u

  // Get index variables
  int i = VertexIndex(u);
  R3Point p1 = VertexPosition(i);
  R3Point p2 = VertexPosition(i+1);
  RNScalar u1 = VertexParameter(i);
  RNScalar u2 = VertexParameter(i+1);
  if (u2 <= u1) return R3zero_vector;
  else return (p2 - p1) / (u2 - u1);
}



void R3Polyline::
Transform(const R3Transformation& transformation)
{
  // Transform vertices
  for (int i = 0; i < nvertices; i++) {
    vertex_positions[i].Transform(transformation);
  }

  // Mark bbox as out of date
  bbox[0][0] = bbox[0][1] = bbox[0][2] = FLT_MAX;
  bbox[1][0] = bbox[1][1] = bbox[1][2] = -FLT_MAX;
}



void R3Polyline::
SetVertexPosition(int k, const R3Point& position)
{
  // Set vertex position
  vertex_positions[k] = position;

  // Mark bbox as out of date
  bbox[0][0] = bbox[0][1] = bbox[0][2] = FLT_MAX;
  bbox[1][0] = bbox[1][1] = bbox[1][2] = -FLT_MAX;
}



void R3Polyline::
SetVertexParameter(int k, RNScalar u)
{
  // Check number of points
  if (nvertices == 0) return;
  
  // Allocate vertex parameters
  if (!vertex_parameters) {
    vertex_parameters = new RNScalar [ nvertices ];
    for (int i = 0; i < k; i++) vertex_parameters[i] = (u*i)/k;
    vertex_parameters[k] = u;
    for (int i = k+1; i < nvertices; i++) vertex_parameters[i] = (i-k)+u;
  }
  else {
    // Set vertex parameter
    vertex_parameters[k] = u;
  }
}



void R3Polyline::
SetVertexData(int k, void *data)
{
  // Allocate vertex datas
  if (!vertex_datas) {
    vertex_datas = new void * [ NVertices() ];
    if (!vertex_datas) {
      fprintf(stderr, "Unable to allocate vertex datas\n");
      return;
    }
  }

  // Set vertex data
  vertex_datas[k] = data;
}



void R3Polyline::
UpdateBBox(void) 
{
  // Initialize bounding box
  bbox = R3null_box;

  // Union vertex_positions
  for (int i = 0; i < nvertices; i++) {
    bbox.Union(VertexPosition(i));
  }
}



void R3Polyline::
Draw(const RNFlags draw_flags) const
{
  // Draw polyline
  Outline();
}



void R3Polyline::
Outline(const R3DrawFlags flags) const
{
  // Draw polyline
  glBegin(GL_LINE_STRIP);
  for (int i = 0; i < nvertices; i++) 
    R3LoadPoint(VertexPosition(i));
  glEnd();
}
