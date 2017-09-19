/* Source file for the R3 scene element class */



/* Include files */

#include "R3Graphics.h"



/* Member functions */

R3SceneElement::
R3SceneElement(R3Material *material)
  : node(NULL),
    material(material),
    shapes(),
    opengl_id(0),
    bbox(FLT_MAX, FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX)
{
}



R3SceneElement::
R3SceneElement(const R3SceneElement& element)
  : node(NULL),
    material(element.material),
    shapes(),
    opengl_id(0),
    bbox(element.bbox)
{
}



R3SceneElement::
~R3SceneElement(void)
{
  // Delete display list
  if (opengl_id > 0) glDeleteLists(opengl_id, 1); 

  // Remove from node
  if (node) node->RemoveElement(this);
}



const R3Point R3SceneElement::
Centroid(void) const
{
  // Return centroid
  return BBox().Centroid();
}



const RNInterval R3SceneElement::
NFacets(void) const
{
  // Initialize nfacets
  RNInterval nfacets(0,0);
  
  // Add facets from shapes
  for (int i = 0; i < NShapes(); i++) {
    R3Shape *shape = Shape(i);
    nfacets += shape->NFacets();
  }
  
  // Return nfacets
  return nfacets;
}



const RNLength R3SceneElement::
Length(void) const
{
  // Initialize length
  RNLength length = 0;
  
  // Add length from shapes
  for (int i = 0; i < NShapes(); i++) {
    R3Shape *shape = Shape(i);
    length += shape->Length();
  }
  
  // Return length
  return length;
}



const RNArea R3SceneElement::
Area(void) const
{
  // Initialize area
  RNArea area = 0;
  
  // Add area from shapes
  for (int i = 0; i < NShapes(); i++) {
    R3Shape *shape = Shape(i);
    area += shape->Area();
  }
  
  // Return area
  return area;
}



const RNVolume R3SceneElement::
Volume(void) const
{
  // Initialize volume
  RNVolume volume = 0;
  
  // Add volume from shapes
  for (int i = 0; i < NShapes(); i++) {
    R3Shape *shape = Shape(i);
    volume += shape->Volume();
  }

  // Return volume
  return volume;
}



void R3SceneElement::
SetMaterial(R3Material *material) 
{
  // Set material
  this->material = material;
}



void R3SceneElement::
InsertShape(R3Shape *shape) 
{
  // Insert shape
  shapes.Insert(shape);

  // Invalidate bounding box
  InvalidateBBox();
}



void R3SceneElement::
RemoveShape(R3Shape *shape) 
{
  // Remove shape
  shapes.Remove(shape);

  // Invalidate bounding box
  InvalidateBBox();
}



void R3SceneElement::
Transform(const R3Transformation& transformation)
{
  // Transform all shapes
  for (int i = 0; i < NShapes(); i++) {
    R3Shape *shape = Shape(i);
    shape->Transform(transformation);
  }


  // Invalidate bounding box
  InvalidateBBox();
}



RNLength R3SceneElement::
Distance(const R3Point& point) const
{
  // Return distance from point to closest point in any shape
  RNLength distance = RN_INFINITY;
  if (!FindClosest(point, NULL, NULL, NULL, &distance)) return RN_INFINITY;
  return distance;
}



RNBoolean R3SceneElement::
FindClosest(const R3Point& point, R3Shape **hit_shape,
  R3Point *hit_point, R3Vector *hit_normal, RNLength *hit_d,
  RNLength min_d, RNLength max_d) const
{
  // Check if bounding box is within max_d
  RNScalar bbox_d = R3Distance(point, BBox());
  if (bbox_d > max_d) return FALSE;

  // Find distance to closest shape
  RNBoolean found = FALSE;
  for (int i = 0; i < NShapes(); i++) {
    R3Shape *shape = Shape(i);
    R3Point closest = shape->ClosestPoint(point);
    RNLength d = R3Distance(closest, point);
    if ((d >= min_d) && (d <= max_d)) {
      if (hit_shape) *hit_shape = shape;
      if (hit_point) *hit_point = closest;
      if (hit_normal) *hit_normal = R3zero_vector; // ???
      if (hit_d) *hit_d = d;
      found = TRUE;
      max_d = d;
    }
  }

  // Return whether hit any shape
  return found;
}



RNBoolean R3SceneElement::
Intersects(const R3Ray& ray, R3Shape **hit_shape,
  R3Point *hit_point, R3Vector *hit_normal, RNScalar *hit_t,
  RNScalar min_t, RNScalar max_t) const
{
  // Variables
  RNScalar bbox_t;
  RNScalar closest_t = max_t;
  R3Point point;
  R3Vector normal;
  RNScalar t;

  // Check if ray intersects bounding box
  if (!R3Contains(BBox(), ray.Start())) {
    if (!R3Intersects(ray, BBox(), NULL, NULL, &bbox_t)) return FALSE;
    if (RNIsGreater(bbox_t, max_t)) return FALSE;
  }

  // Intersect with shapes
  for (int i = 0; i < NShapes(); i++) {
    R3Shape *shape = Shape(i);
    if (shape->Intersects(ray, &point, &normal, &t)) {
      if ((t >= min_t) && (t <= closest_t)) {
        if (hit_shape) *hit_shape = shape;
        if (hit_point) *hit_point = point;
        if (hit_normal) *hit_normal = normal;
        if (hit_t) *hit_t = t;
        closest_t = t;
      }
    }
  }

  // Return whether hit any shape
  return (closest_t == max_t) ? FALSE : TRUE;
}



void R3SceneElement::
Draw(const R3DrawFlags draw_flags, const RNArray<R3Material *> *materials) const
{
  // Check shapes
  if (NShapes() == 0) return;

  // Draw material
  if (material && draw_flags[R3_SURFACE_MATERIAL_DRAW_FLAG]) {
    int material_index = material->SceneIndex();
    if (materials && (material_index >= 0) && (material_index < materials->NEntries()) && materials->Kth(material_index)) materials->Kth(material_index)->Draw();
    else material->Draw();
  }

  // Draw shapes
  if (0 && (draw_flags == R3_DEFAULT_DRAW_FLAGS)) {
    // Create display list
    if (opengl_id == 0) {
      // Begin display list
      R3SceneElement *element = (R3SceneElement *) this;
      element->opengl_id = glGenLists(1);
      glNewList(opengl_id, GL_COMPILE);
      
      // Draw shapes
      for (int i = 0; i < NShapes(); i++) {
        R3Shape *shape = Shape(i);
        shape->Draw(draw_flags);
      }

      // End display list
      glEndList();
    }

    // Call display list
    glCallList(opengl_id);
  }
  else {
    // Draw shapes with non-default draw flags
    for (int i = 0; i < NShapes(); i++) {
      R3Shape *shape = Shape(i);
      shape->Draw(draw_flags);
    }
  }
}



void R3SceneElement::
UpdateBBox(void)
{
  // Update bounding box
  bbox = R3null_box;
  for (int i = 0; i < NShapes(); i++) {
    R3Shape *shape = Shape(i);
    bbox.Union(shape->BBox());
  }
}



void R3SceneElement::
InvalidateBBox(void)
{
  // Invalidate bounding box
  bbox[0][0] = FLT_MAX;

  // Invalidate node bounding box
  if (node) node->InvalidateBBox();
}

