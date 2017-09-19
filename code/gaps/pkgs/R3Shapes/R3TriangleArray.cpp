/* Source file for the R3 triangle class */



/* Include files */

#include "R3Shapes/R3Shapes.h"



/* Class type definitions */

RN_CLASS_TYPE_DEFINITIONS(R3TriangleArray);



/* Public functions */

int 
R3InitTriangleArray()
{
    /* Return success */
    return TRUE;
}



void 
R3StopTriangleArray()
{
}



R3TriangleArray::
R3TriangleArray(void)
    : bbox(R3null_box)
{
}



R3TriangleArray::
R3TriangleArray(const R3TriangleArray& array)
  : vertices(),
    triangles(),
    bbox(array.bbox)
{
    // Copy vertices
    for (int i = 0; i < array.vertices.NEntries(); i++) {
      R3TriangleVertex *oldv = array.vertices.Kth(i);
      R3TriangleVertex *newv = new R3TriangleVertex(*oldv);
      vertices.Insert(newv);
      newv->SetSharedFlag();
      oldv->SetMark(i);
    }

    // Copy triangles
    for (int i = 0; i < array.triangles.NEntries(); i++) {
      R3Triangle *oldt = array.triangles.Kth(i);
      R3TriangleVertex *v0 = vertices.Kth(oldt->V0()->Mark());
      R3TriangleVertex *v1 = vertices.Kth(oldt->V1()->Mark());
      R3TriangleVertex *v2 = vertices.Kth(oldt->V2()->Mark());
      R3Triangle *newt = new R3Triangle(v0, v1, v2);
      triangles.Insert(newt);
    }

    // Update everything
    Update();
}



R3TriangleArray::
R3TriangleArray(const RNArray<R3TriangleVertex *>& vertices, const RNArray<R3Triangle *>& triangles)
  : vertices(vertices),
    triangles(triangles),
    bbox(R3null_box)
{
    // Update everything
    Update();
}



R3TriangleArray::
~R3TriangleArray(void)
{
    // Delete triangles and vertices
    for (int i = 0; i < triangles.NEntries(); i++) delete triangles[i];
    for (int i = 0; i < vertices.NEntries(); i++) delete vertices[i];
}



const RNBoolean R3TriangleArray::
IsPoint (void) const
{
  // Check if all vertices are at same position
  if (vertices.NEntries() > 1) {
    const R3Point& p0 = vertices[0]->Position();
    for (int i = 1; i < vertices.NEntries(); i++) 
      if (!R3Contains(p0, vertices[i]->Position())) return FALSE;
  }
  return TRUE;
}



const RNBoolean R3TriangleArray::
IsLinear (void) const
{
    // Assume not ???
    return FALSE;
}



const RNBoolean R3TriangleArray::
IsPlanar(void) const
{
  // Check if all triangles are on same plane
  if (triangles.NEntries() > 1) {
    const R3Plane& p0 = triangles[0]->Plane();
    for (int i = 1; i < triangles.NEntries(); i++) 
      if (!R3Contains(p0, triangles[i]->Plane())) return FALSE;
  }
  return TRUE;
}



const RNBoolean R3TriangleArray::
IsConvex(void) const
{
    // If planar - may still be concave ???
    return IsPlanar();
}



const RNInterval R3TriangleArray::
NFacets(void) const
{
    // Return number of trianglets (triangles)
    return RNInterval(triangles.NEntries(), triangles.NEntries());
}



const RNLength R3TriangleArray::
Length (void) const
{
    // Return cumulative perimeter of triangles
    RNLength length = 0.0;
    for (int i = 0; i < triangles.NEntries(); i++)
      length += triangles[i]->Length();
    return length;
}



const RNArea R3TriangleArray::
Area(void) const
{
    // Return cumulative area of triangles
    RNArea area = 0.0;
    for (int i = 0; i < triangles.NEntries(); i++)
      area += triangles[i]->Area();
    return area;
}



const R3Point R3TriangleArray::
Centroid(void) const
{
    // Return centroid of vertices
    R3Point centroid = R3zero_point;
    for (int i = 0; i < vertices.NEntries(); i++) 
      centroid += vertices[i]->Position();
    centroid /= (RNScalar) vertices.NEntries();
    return centroid;
}



const R3Point R3TriangleArray::
ClosestPoint(const R3Point& point) const
{
    // Find closest point
    R3Point closest_point = R3zero_point;
    RNLength closest_dd = RN_INFINITY * RN_INFINITY;
    for (int i = 0; i < triangles.NEntries(); i++) {
      R3Triangle *triangle = triangles.Kth(i);
      R3Point p = triangle->ClosestPoint(point);
      RNLength dd = R3SquaredDistance(p, point);
      if (dd < closest_dd) {
        closest_point = p;
        closest_dd = dd;
      }
    }
        
    // Return closest point
    return closest_point;
}



const R3Point R3TriangleArray::
FurthestPoint(const R3Point& point) const
{
    // Find furthest point
    R3Point furthest_point = R3zero_point;
    RNLength furthest_dd = 0;
    for (int i = 0; i < vertices.NEntries(); i++) {
      R3TriangleVertex *vertex = vertices.Kth(i);
      R3Point p = vertex->Position();
      RNLength dd = R3SquaredDistance(p, point);
      if (dd > furthest_dd) {
        furthest_point = p;
        furthest_dd = dd;
      }
    }

    // Return furthest point
    return furthest_point;
}



const R3Shape& R3TriangleArray::
BShape(void) const
{
    // Return bounding shape
    return bbox;
}



const R3Box R3TriangleArray::
BBox(void) const
{
    // Return bounding box
    return bbox;
}



const R3Sphere R3TriangleArray::
BSphere(void) const
{
    // Return bounding sphere
    return bbox.BSphere();
}



void R3TriangleArray::
Flip (void)
{
    // Flip the triangles
    for (int i = 0; i < triangles.NEntries(); i++)
      triangles[i]->Flip();
}



void R3TriangleArray::
Mirror (const R3Plane& plane)
{
    // Mirror vertices
    for (int i = 0; i < vertices.NEntries(); i++) 
      vertices[i]->Mirror(plane);

    // Update the triangles
    for (int i = 0; i < triangles.NEntries(); i++)
      triangles[i]->Update();

    // Update the bounding box
    Update();
}



void R3TriangleArray::
Transform (const R3Transformation& transformation) 
{
    // Transform vertices
    for (int i = 0; i < vertices.NEntries(); i++) 
      vertices[i]->Transform(transformation);

    // Update the triangles
    for (int i = 0; i < triangles.NEntries(); i++)
      triangles[i]->Update();

    // Check if need to flip triangles
    if (transformation.IsMirrored()) Flip();

    // Update the bounding box
    Update();
}



void R3TriangleArray::
Subdivide(RNLength max_edge_length) 
{
  // Get convenient variables
  RNLength max_edge_length_squared = max_edge_length * max_edge_length;
  
  // Subdivide triangles with edges that are too long
  RNArray<R3Triangle *> stack = triangles;
  while (!stack.IsEmpty()) {
    // Get triangle from stack
    R3Triangle *triangle = stack.Tail();
    stack.RemoveTail();

    // Get vertices
    R3TriangleVertex *v0 = triangle->V0();
    R3TriangleVertex *v1 = triangle->V1();
    R3TriangleVertex *v2 = triangle->V2();
    const R3Point& p0 = v0->Position();
    const R3Point& p1 = v1->Position();
    const R3Point& p2 = v2->Position();
    RNLength dd01 = R3SquaredDistance(p0, p1);
    RNLength dd12 = R3SquaredDistance(p1, p2);
    RNLength dd20 = R3SquaredDistance(p2, p0);

    // Check if edge is too long
    if ((dd01 > max_edge_length_squared) && (dd01 >= dd12) && (dd01 >= dd20)) {
      // Create new vertex at midpoint of edge between v0 and v1
      R3Point position = 0.5*(p0 + p1);
      R3Vector normal = 0.5*(v0->Normal() + v1->Normal()); normal.Normalize();
      R3TriangleVertex *v = new R3TriangleVertex(position, normal);
      if (v0->flags[R3_VERTEX_TEXTURE_COORDS_DRAW_FLAG] && v1->flags[R3_VERTEX_TEXTURE_COORDS_DRAW_FLAG])
        v->SetTextureCoords(0.5*(v0->TextureCoords() + v1->TextureCoords()));
      vertices.Insert(v);

      // Adjust existing triangle
      triangle->Reset(v0, v, v2);
      stack.Insert(triangle);

      // Create new triangle
      R3Triangle *t = new R3Triangle(v1, v2, v);
      triangles.Insert(t);
      stack.Insert(t);
    }
    else if ((dd12 > max_edge_length_squared) && (dd12 >= dd01) && (dd12 >= dd20)) {
      // Create new vertex at midpoint of edge between v1 and v2
      R3Point position = 0.5*(p1 + p2);
      R3Vector normal = 0.5*(v1->Normal() + v2->Normal()); normal.Normalize();
      R3TriangleVertex *v = new R3TriangleVertex(position, normal);
      if (v1->flags[R3_VERTEX_TEXTURE_COORDS_DRAW_FLAG] && v2->flags[R3_VERTEX_TEXTURE_COORDS_DRAW_FLAG])
        v->SetTextureCoords(0.5*(v1->TextureCoords() + v2->TextureCoords()));
      vertices.Insert(v);

      // Adjust existing triangle
      triangle->Reset(v1, v, v0);
      stack.Insert(triangle);

      // Create new triangle
      R3Triangle *t = new R3Triangle(v2, v0, v);
      triangles.Insert(t);
      stack.Insert(t);
    }
    else if ((dd20 > max_edge_length_squared) && (dd20 >= dd01) && (dd20 >= dd12)) {
      // Create new vertex at midpoint of edge between v2 and v0
      R3Point position = 0.5*(p2 + p0);
      R3Vector normal = 0.5*(v2->Normal() + v0->Normal()); normal.Normalize();
      R3TriangleVertex *v = new R3TriangleVertex(position, normal);
      if (v2->flags[R3_VERTEX_TEXTURE_COORDS_DRAW_FLAG] && v0->flags[R3_VERTEX_TEXTURE_COORDS_DRAW_FLAG])
        v->SetTextureCoords(0.5*(v2->TextureCoords() + v0->TextureCoords()));
      vertices.Insert(v);

      // Adjust existing triangle
      triangle->Reset(v2, v, v1);
      stack.Insert(triangle);

      // Create new triangle
      R3Triangle *t = new R3Triangle(v0, v1, v);
      triangles.Insert(t);
      stack.Insert(t);
    }
  }

  // Update
  Update();
}



void R3TriangleArray::
MoveVertex(R3TriangleVertex *vertex, const R3Point& position)
{
    // Move vertex
    vertex->SetPosition(position);

    // Update triangles
    for (int i = 0; i < triangles.NEntries(); i++) {
      R3Triangle *triangle = triangles.Kth(i);
      for (int j = 0; j < 3; j++) {
        if (triangle->Vertex(i) == vertex) {
          triangle->Update();
          break;
        }
      }
    }

    // Update
    Update();
}



void R3TriangleArray::
Update(void)
{
    // Mark vertices as shared
    for (int i = 0; i < vertices.NEntries(); i++) {
      R3TriangleVertex *v = vertices.Kth(i);
      v->SetSharedFlag();
    }
    
    // Recompute bounding box
    bbox = R3null_box;
    for (int i = 0; i < vertices.NEntries(); i++) {
      R3TriangleVertex *v = vertices.Kth(i);
      bbox.Union(v->Position());
    }
}




