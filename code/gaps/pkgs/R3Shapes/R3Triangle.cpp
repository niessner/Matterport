/* Source file for the R3 triangle class */



/* Include files */

#include "R3Shapes/R3Shapes.h"



/* Class type definitions */

RN_CLASS_TYPE_DEFINITIONS(R3Triangle);



/* Initialization functions ********************/

int 
R3InitTriangle()
{
    /* Return success */
    return TRUE;
}



void 
R3StopTriangle()
{
}



/* Triangle vertex functions ********************/

R3TriangleVertex::
R3TriangleVertex(void)
    : position(0,0,0),
      normal(0,0,0),
      color(0, 0, 0),
      texcoords(0,0),
      flags(0),
      mark(0)
{
}



R3TriangleVertex::
R3TriangleVertex(const R3TriangleVertex& vertex)
    : position(vertex.position), 
      normal(vertex.normal),
      color(vertex.color),
      texcoords(vertex.texcoords),
      flags(vertex.flags),
      mark(0)
{
}



R3TriangleVertex::
R3TriangleVertex(const R3Point& position)
    : position(position),
      normal(0,0,0),
      color(0, 0, 0),
      texcoords(0,0),
      flags(0),
      mark(0)
{
}



R3TriangleVertex::
R3TriangleVertex(const R3Point& position, const R3Vector& normal)
    : position(position),
      normal(normal),
      color(0,0,0),
      texcoords(0,0),
      flags(0),
      mark(0)
{
  // Check normal
  if (!normal.IsZero()) flags.Add(R3_VERTEX_NORMALS_DRAW_FLAG);
}



R3TriangleVertex::
R3TriangleVertex(const R3Point& position, const RNRgb& color)
    : position(position),
      normal(0,0,0),
      color(color),
      texcoords(0,0),
      flags(0),
      mark(0)
{
  // Remember that color was set
  flags.Add(R3_VERTEX_COLORS_DRAW_FLAG);
}



R3TriangleVertex::
R3TriangleVertex(const R3Point& position, const R2Point& texcoords)
    : position(position),
      normal(0,0,0),
      color(0, 0, 0),
      texcoords(texcoords),
      flags(R3_VERTEX_TEXTURE_COORDS_DRAW_FLAG),
      mark(0)
{
}



R3TriangleVertex::
R3TriangleVertex(const R3Point& position, const R3Vector& normal, const RNRgb& color)
    : position(position),
      normal(normal),
      color(color),
      texcoords(0, 0),
      flags(R3_VERTEX_COLORS_DRAW_FLAG),
      mark(0)
{
  // Check normal
  if (!normal.IsZero()) flags.Add(R3_VERTEX_NORMALS_DRAW_FLAG);
}



R3TriangleVertex::
R3TriangleVertex(const R3Point& position, const R3Vector& normal, const RNRgb& color, const R2Point& texcoords)
    : position(position),
      normal(normal),
      color(color),
      texcoords(texcoords),
      flags(R3_VERTEX_COLORS_DRAW_FLAG | R3_VERTEX_TEXTURE_COORDS_DRAW_FLAG),
      mark(0)
{
  // Check normal
  if (!normal.IsZero()) flags.Add(R3_VERTEX_NORMALS_DRAW_FLAG);
}



R3TriangleVertex::
R3TriangleVertex(const R3Point& position, const R3Vector& normal, const R2Point& texcoords)
    : position(position),
      normal(normal),
      color(0, 0, 0),
      texcoords(texcoords),
      flags(R3_VERTEX_TEXTURE_COORDS_DRAW_FLAG),
      mark(0)
{
  // Check normal
  if (!normal.IsZero()) flags.Add(R3_VERTEX_NORMALS_DRAW_FLAG);
}



void R3TriangleVertex::
Mirror (const R3Plane& plane)
{
  // Mirror position and normal
  position.Mirror(plane);
  normal.Mirror(plane);
}



void R3TriangleVertex::
Transform (const R3Transformation& transformation)
{
  // Transform position and normal
  transformation.Apply(position);
  transformation.ApplyInverseTranspose(normal);
}



/* Triangle functions ********************/

R3Triangle::
R3Triangle(void)
    : plane(R3null_plane),
      bbox(R3null_box),
      flags(0),
      mark(0)
{
    // Initialize vertices
    v[0] = v[1] = v[2] = NULL;
}



R3Triangle::
R3Triangle(const R3Triangle& triangle)
  : plane(triangle.plane),
    bbox(triangle.bbox),
    flags(triangle.flags),
    mark(triangle.mark)
{
    // Copy vertices 
    if (!triangle.v[0] || triangle.v[0]->Flags()[R3_VERTEX_SHARED_FLAG]) v[0] = triangle.v[0];
    else v[0] = new R3TriangleVertex(*(triangle.v[0]));
    if (!triangle.v[1] || triangle.v[1]->Flags()[R3_VERTEX_SHARED_FLAG]) v[1] = triangle.v[1];
    else v[1] = new R3TriangleVertex(*(triangle.v[1]));
    if (!triangle.v[2] || triangle.v[2]->Flags()[R3_VERTEX_SHARED_FLAG]) v[2] = triangle.v[2];
    else v[2] = new R3TriangleVertex(*(triangle.v[2]));

    // Update the triangle
    Update();
}



R3Triangle::
R3Triangle(R3TriangleVertex *v0, R3TriangleVertex *v1, R3TriangleVertex *v2)
  : plane(R3null_plane),
    bbox(R3null_box),
    flags(0),
    mark(0)
{
    // NOTE: The vertices will be deleted when the triangle is deleted if they do not have not been SetSharedFlag
  
    // Reset triangle
    Reset(v0, v1, v2);
}



R3Triangle::
R3Triangle(R3TriangleVertex *vertices[3])
  : plane(R3null_plane),
    bbox(R3null_box),
    flags(0),
    mark(0)
{
    // NOTE: The vertices will be deleted when the triangle is deleted if they do not have not been SetSharedFlag

    // Reset triangle
    Reset(vertices[0], vertices[1], vertices[2]);
}



R3Triangle::
~R3Triangle(void)
{
  // Delete vertices
  if (v[0] && !v[0]->Flags()[R3_VERTEX_SHARED_FLAG]) delete v[0];
  if (v[1] && !v[1]->Flags()[R3_VERTEX_SHARED_FLAG]) delete v[1];
  if (v[2] && !v[2]->Flags()[R3_VERTEX_SHARED_FLAG]) delete v[2];
}



const RNBoolean R3Triangle::
IsPoint (void) const
{
    // Return whether triangle covers only one vertex
    if (!R3Contains(v[0]->Position(), v[1]->Position())) return FALSE;
    if (!R3Contains(v[0]->Position(), v[2]->Position())) return FALSE;
    return TRUE;
}



const RNBoolean R3Triangle::
IsLinear (void) const
{
    // Return whether triangle lies on a line
    R3Line line(v[0]->Position(), v[1]->Position());
    if (R3Contains(line, v[2]->Position())) return TRUE;
    return FALSE;
}



const RNBoolean R3Triangle::
IsPlanar(void) const
{
    // All triangles are planar
    return TRUE;
}



const RNBoolean R3Triangle::
IsConvex(void) const
{
    // Assume all triangles are convex
    return TRUE;
}



const RNBoolean R3Triangle::
IsFinite (void) const
{
    // Check each point
    if (!(v[0]->Position().IsFinite())) return FALSE;
    if (!(v[1]->Position().IsFinite())) return FALSE;
    if (!(v[2]->Position().IsFinite())) return FALSE;
    return TRUE;
}



const RNInterval R3Triangle::
NFacets(void) const
{
    // Return number of trianglets (triangles)
    return RNInterval(1.0, 1.0);
}



const RNLength R3Triangle::
Length (void) const
{
    // Return perimeter of closed triangle
    RNLength length = 0.0;
    length += R3Distance(v[0]->Position(), v[1]->Position());
    length += R3Distance(v[1]->Position(), v[2]->Position());
    length += R3Distance(v[2]->Position(), v[0]->Position());
    return length;
}



const RNArea R3Triangle::
Area(void) const
{
    // Compute area using Newell's method
    R3Vector sum = R3null_vector;
    sum += v[0]->Position().Vector() % v[2]->Position().Vector();
    sum += v[1]->Position().Vector() % v[0]->Position().Vector();
    sum += v[2]->Position().Vector() % v[1]->Position().Vector();
    return 0.5 * sum.Length();
}



const R3Point R3Triangle::
Centroid(void) const
{
    // Return centroid of vertex set
    R3Point centroid = R3zero_point;
    centroid += v[0]->Position();
    centroid += v[1]->Position();
    centroid += v[2]->Position();
    centroid /= 3;
    return centroid;
}



const R3Point R3Triangle::
RandomPoint(void) const
{
  // Get vertex positions
  const R3Point& p0 = v[0]->Position();
  const R3Point& p1 = v[1]->Position();
  const R3Point& p2 = v[2]->Position();
  
  // Return random point on triangle
  RNScalar r1 = sqrt(RNRandomScalar());
  RNScalar r2 = RNRandomScalar();
  R3Point p = p0 * (1.0 - r1) + p1 * r1 * (1.0 - r2) + p2 * r1 * r2;
  return p;
}



const R3Point R3Triangle::
ClosestPoint(const R3Point& point) const
{
    // Start with projection onto plane
    R3Point closest_point = point;
    closest_point.Project(Plane());

    // Project onto edge plane if point is on other side
    for (int i = 0; i < 3; i++) {
      R3TriangleVertex *v0 = v[i];
      R3TriangleVertex *v1 = v[(i+1)%3];
      R3Point p0 = v0->Position();
      R3Point p1 = v1->Position();
      R3Span span(p0, p1);
      if (RNIsZero(span.Length())) return Centroid();
      R3Vector n = span.Vector() % Normal();
      R3Plane edge_plane(p0, n);
      RNScalar d1 = R3SignedDistance(edge_plane, closest_point);
      if (RNIsPositive(d1)) {
        RNScalar t = span.T(closest_point);
        if (t < 0) t = 0;
        if (t > span.Length()) t = span.Length();
        closest_point = span.Point(t);
      }
    }
   
    // Return closest point;
    return closest_point;
}



const R3Point R3Triangle::
FurthestPoint(const R3Point& point) const
{
    // Find furthest point
    RNLength furthest_dd = 0;
    R3Point furthest_point = R3zero_point;
    for (int i = 0; i < 3; i++) {
      R3TriangleVertex *vertex = v[i];
      R3Point p = vertex->Position();
      RNLength dd = R3SquaredDistance(p, point);
      if (dd > furthest_dd) {
        furthest_dd = dd;
        furthest_point = p;
      }
    }

    // Return furthest point
    return furthest_point;
}



const R3Shape& R3Triangle::
BShape(void) const
{
    // Return bounding shape
    return bbox;
}



const R3Box R3Triangle::
BBox(void) const
{
    // Return bounding box
    return bbox;
}



const R3Sphere R3Triangle::
BSphere(void) const
{
    // Return bounding sphere (is this right???)
    R3Point centroid = Centroid();
    RNLength radius = R3Distance(centroid, v[0]->Position());
    return R3Sphere(centroid, radius);
}



R3Triangle R3Triangle::
operator-(void) const
{
    // Return flipped copy of triangle
    R3Triangle triangle(*this);
    triangle.Flip();
    return triangle;
}



void R3Triangle::
Flip (void) 
{
    // Flip orientation of plane
    plane.Flip();

    // Reverse order of vertices
    R3TriangleVertex *swap = v[0];
    v[0] = v[2];
    v[2] = swap;
}



void R3Triangle::
Mirror (const R3Plane& plane)
{
    // Mirror vertices (what if they are shared?)
    v[0]->Mirror(plane);
    v[1]->Mirror(plane);
    v[2]->Mirror(plane);

    // Reverse order of vertices
    R3TriangleVertex *swap = v[0];
    v[0] = v[2];
    v[2] = swap;

    // Update
    Update();
}



void R3Triangle::
Transform (const R3Transformation& transformation) 
{
    // Transform vertices
    v[0]->Transform(transformation);
    v[1]->Transform(transformation);
    v[2]->Transform(transformation);

    // Check if need to flip orientation
    if (transformation.IsMirrored()) {
      R3TriangleVertex *swap = v[0];
      v[0] = v[2];
      v[2] = swap;
    }

    // Update
    Update();
}



void R3Triangle::
Reset(R3TriangleVertex *v0, R3TriangleVertex *v1, R3TriangleVertex *v2)
{
    // Reset vertices
    v[0] = v0;
    v[1] = v1;
    v[2] = v2;
    
    // Update
    Update();
}



void R3Triangle::
Reset(R3TriangleVertex *vertices[3])
{
    // Reset vertices
    v[0] = vertices[0];
    v[1] = vertices[1];
    v[2] = vertices[2];
    
    // Update
    Update();
}



void R3Triangle::
SetMark(RNMark mark)
{
    // Set mark
    this->mark = mark;
}



void R3Triangle::
Update(void)
{
    // Recompute plane equation
    plane = R3Plane(v[0]->Position(), v[1]->Position(), v[2]->Position()); 

    // Recompute bounding box
    bbox = R3null_box;
    bbox.Union(v[0]->Position());
    bbox.Union(v[1]->Position());
    bbox.Union(v[2]->Position());

    // Recompute flags
    flags.Add(R3_EVERYTHING_DRAW_FLAGS);
    for (int i = 0; i < 3; i++) {
	RNFlags vertex_flags = v[i]->Flags();
	if (!(vertex_flags[R3_VERTEX_NORMALS_DRAW_FLAG]))
	    flags.Remove(R3_VERTEX_NORMALS_DRAW_FLAG);
	if (!(vertex_flags[R3_VERTEX_COLORS_DRAW_FLAG]))
	    flags.Remove(R3_VERTEX_COLORS_DRAW_FLAG);
	if (!(vertex_flags[R3_VERTEX_TEXTURE_COORDS_DRAW_FLAG]))
	    flags.Remove(R3_VERTEX_TEXTURE_COORDS_DRAW_FLAG);
    }
}




