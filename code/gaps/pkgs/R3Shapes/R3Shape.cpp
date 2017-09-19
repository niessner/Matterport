/* Source file for the R3 shape class */



/* Include files */

#include "R3Shapes/R3Shapes.h"



/* Class type definitions */

RN_CLASS_TYPE_DEFINITIONS(R3Shape);



/* Public functions */

int 
R3InitShape()
{
    /* Return success */
    return TRUE;
}



void 
R3StopShape()
{
}



R3Shape *R3Shape::
Copy(void) const
{
  // Return copy of appropriate type
  if (ClassID() == R3Box::CLASS_ID()) return new R3Box(*((R3Box *) this));
  else if (ClassID() == R3OrientedBox::CLASS_ID()) return new R3OrientedBox(*((R3OrientedBox *) this));
  else if (ClassID() == R3Cylinder::CLASS_ID()) return new R3Cylinder(*((R3Cylinder *) this));
  else if (ClassID() == R3Cone::CLASS_ID()) return new R3Cone(*((R3Cone *) this));
  else if (ClassID() == R3Ellipsoid::CLASS_ID()) return new R3Ellipsoid(*((R3Ellipsoid *) this));
  else if (ClassID() == R3Sphere::CLASS_ID()) return new R3Sphere(*((R3Sphere *) this));
  else if (ClassID() == R3Triangle::CLASS_ID()) return new R3Triangle(*((R3Triangle *) this));
  else if (ClassID() == R3TriangleArray::CLASS_ID()) return new R3TriangleArray(*((R3TriangleArray *) this));
  else if (ClassID() == R3Circle::CLASS_ID()) return new R3Circle(*((R3Circle *) this));
  else if (ClassID() == R3Ellipse::CLASS_ID()) return new R3Ellipse(*((R3Ellipse *) this));
  else if (ClassID() == R3Rectangle::CLASS_ID()) return new R3Rectangle(*((R3Rectangle *) this));
  else if (ClassID() == R3Polyline::CLASS_ID()) return new R3Polyline(*((R3Polyline *) this));
  else if (ClassID() == R3CatmullRomSpline::CLASS_ID()) return new R3CatmullRomSpline(*((R3CatmullRomSpline *) this));

  // Should never get here
  RNAbort("Unrecognized shape class id: %d\n", ClassID());
  return NULL;
}



const RNBoolean R3Shape::
IsPoint(void) const
{
    // Assume the worst - this may be overridden
    return FALSE;
}



const RNBoolean R3Shape::
IsCurve(void) const
{
    // Assume the worst - this may be overridden
    return FALSE;
}



const RNBoolean R3Shape::
IsLinear(void) const
{
    // Assume the worst - this may be overridden
    return FALSE;
}



const RNBoolean R3Shape::
IsSurface(void) const
{
    // Assume the worst - this may be overridden
    return FALSE;
}



const RNBoolean R3Shape::
IsPlanar(void) const
{
    // Assume the worst - this may be overridden
    return FALSE;
}



const RNBoolean R3Shape::
IsSolid(void) const
{
    // Assume the worst - this may be overridden
    return FALSE;
}



const RNBoolean R3Shape::
IsConvex(void) const
{
    // Assume the worst - this may be overridden
    return FALSE;
}



const RNInterval R3Shape::
NFacets(void) const
{
    // This should be overridden
    return RNzero_interval;
}



const RNArea R3Shape::
Length(void) const
{
    // Assume the worst - this may be overridden
    return 0.0;
}



const RNArea R3Shape::
Area(void) const
{
    // Assume the worst - this may be overridden
    return 0.0;
}



const RNVolume R3Shape::
Volume(void) const
{
    // Assume the worst - this may be overridden
    return 0.0;
}



#if FALSE

const R3Line *R3Shape::
Line(void) const
{
    // May be overridden if the shape is linear
    return NULL;
}



const R3Plane *R3Shape::
Plane(void) const
{
    // May be overridden if the shape is planar
    return NULL;
}

#endif



const R3Point R3Shape::
Centroid(void) const
{
    // Default is centroid of bbox
    return BBox().Centroid();
}



const R3Point R3Shape::
ClosestPoint(const R3Point& point) const
{
    // Assume the worst - this should be overridden
    RNAbort("Not Implemented");
    return R3zero_point;
}



const R3Point R3Shape::
FurthestPoint(const R3Point& point) const
{
    // Assume the worst - this should be overridden
    RNAbort("Not Implemented");
    return R3infinity_point;
}



const R3Shape& R3Shape::
BShape(void) const
{
    // Assume the worst - this may be overridden
    RNAbort("Not Implemented");
    return R3infinite_box;
}



const R3Box R3Shape::
BBox(void) const
{
    // Return bounding shape's bounding box - this may be overridden
    return BShape().BBox();
}



const R3Sphere R3Shape::
BSphere(void) const
{
    // Return bounding shape's bounding sphere - this may be overridden
    return BShape().BSphere();
}



void R3Shape::
Transform(const R3Transformation& /* transformation */)
{
    RNAbort("Not Implemented");
}



void R3Shape::
Draw(const R3DrawFlags ) const
{
    RNAbort("Not Implemented");
}



void R3Shape::
Outline(const R3DrawFlags flags) const
{
    // Draw edges only
    if (flags[R3_EDGES_DRAW_FLAG]) Draw(R3_EDGES_DRAW_FLAG);
}



const R3Box 
R3ShapeBBox(const void *shape, void *)
{
    // R3Shape bounding box function -- useful callback function
    return (*((R3Shape * const *) shape))->BBox();
}



// Shape relationship functions

RNLength R3Shape:: Distance(const R3Point& ) const { RNAbort("Not implemented"); return 0.0; }
RNLength R3Shape:: Distance(const R3Line& ) const { RNAbort("Not implemented"); return 0.0; }
RNLength R3Shape:: Distance(const R3Ray& ) const { RNAbort("Not implemented"); return 0.0; }
RNLength R3Shape:: Distance(const R3Span& ) const { RNAbort("Not implemented"); return 0.0; }
RNLength R3Shape:: Distance(const R3Plane& ) const { RNAbort("Not implemented"); return 0.0; }
RNLength R3Shape:: Distance(const R3Halfspace& ) const { RNAbort("Not implemented"); return 0.0; }
RNLength R3Shape:: Distance(const R3Box& ) const { RNAbort("Not implemented"); return 0.0; }
RNLength R3Shape:: Distance(const R3OrientedBox& ) const { RNAbort("Not implemented"); return 0.0; }
RNLength R3Shape:: Distance(const R3Sphere& ) const { RNAbort("Not implemented"); return 0.0; }
RNLength R3Shape:: Distance(const R3Shape& ) const { RNAbort("Not implemented"); return 0.0; }
RNLength R3Shape:: SignedDistance(const R3Plane& ) const { RNAbort("Not implemented"); return 0.0; }
RNBoolean R3Shape:: Contains(const R3Point& ) const { RNAbort("Not implemented"); return FALSE; }
RNBoolean R3Shape:: Contains(const R3Line& ) const { RNAbort("Not implemented"); return FALSE; }
RNBoolean R3Shape:: Contains(const R3Ray& ) const { RNAbort("Not implemented"); return FALSE; }
RNBoolean R3Shape:: Contains(const R3Span& ) const { RNAbort("Not implemented"); return FALSE; }
RNBoolean R3Shape:: Contains(const R3Plane& ) const { RNAbort("Not implemented"); return FALSE; }
RNBoolean R3Shape:: Contains(const R3Halfspace& ) const { RNAbort("Not implemented"); return FALSE; }
RNBoolean R3Shape:: Contains(const R3Box& ) const { RNAbort("Not implemented"); return FALSE; }
RNBoolean R3Shape:: Contains(const R3OrientedBox& ) const { RNAbort("Not implemented"); return FALSE; }
RNBoolean R3Shape:: Contains(const R3Sphere& ) const { RNAbort("Not implemented"); return FALSE; }
RNBoolean R3Shape:: Contains(const R3Shape& ) const { RNAbort("Not implemented"); return FALSE; }
RNBoolean R3Shape:: Inside(const R3Point& ) const { RNAbort("Not implemented"); return FALSE; }
RNBoolean R3Shape:: Inside(const R3Line& ) const { RNAbort("Not implemented"); return FALSE; }
RNBoolean R3Shape:: Inside(const R3Ray& ) const { RNAbort("Not implemented"); return FALSE; }
RNBoolean R3Shape:: Inside(const R3Span& ) const { RNAbort("Not implemented"); return FALSE; }
RNBoolean R3Shape:: Inside(const R3Plane& ) const { RNAbort("Not implemented"); return FALSE; }
RNBoolean R3Shape:: Inside(const R3Halfspace& ) const { RNAbort("Not implemented"); return FALSE; }
RNBoolean R3Shape:: Inside(const R3Box& ) const { RNAbort("Not implemented"); return FALSE; }
RNBoolean R3Shape:: Inside(const R3OrientedBox& ) const { RNAbort("Not implemented"); return FALSE; }
RNBoolean R3Shape:: Inside(const R3Sphere& ) const { RNAbort("Not implemented"); return FALSE; }
RNBoolean R3Shape:: Inside(const R3Shape& ) const { RNAbort("Not implemented"); return FALSE; }
RNClassID R3Shape:: Intersects(const R3Point& ) const { RNAbort("Not implemented"); return RN_NULL_CLASS_ID; }
RNClassID R3Shape:: Intersects(const R3Line& ) const { RNAbort("Not implemented"); return RN_NULL_CLASS_ID; }
RNClassID R3Shape:: Intersects(const R3Ray& , R3Point *, R3Vector *, RNScalar *) const { RNAbort("Not implemented"); return RN_NULL_CLASS_ID; }
RNClassID R3Shape:: Intersects(const R3Span& ) const { RNAbort("Not implemented"); return RN_NULL_CLASS_ID; }
RNClassID R3Shape:: Intersects(const R3Plane& ) const { RNAbort("Not implemented"); return RN_NULL_CLASS_ID; }
RNClassID R3Shape:: Intersects(const R3Halfspace& ) const { RNAbort("Not implemented"); return RN_NULL_CLASS_ID; }
RNClassID R3Shape:: Intersects(const R3Box& ) const { RNAbort("Not implemented"); return RN_NULL_CLASS_ID; }
RNClassID R3Shape:: Intersects(const R3OrientedBox& ) const { RNAbort("Not implemented"); return RN_NULL_CLASS_ID; }
RNClassID R3Shape:: Intersects(const R3Sphere& ) const { RNAbort("Not implemented"); return RN_NULL_CLASS_ID; }
RNClassID R3Shape:: Intersects(const R3Shape& ) const { RNAbort("Not implemented"); return RN_NULL_CLASS_ID; }


// Shape relationship function definition macro 

#define R3_SHAPE_RELATIONSHIP_DEFINITIONS(__type) \
    RNLength __type:: Distance(const R3Point& point) const { return R3Distance(*this, point); } \
    RNLength __type:: Distance(const R3Line& line) const { return R3Distance(*this, line); } \
    RNLength __type:: Distance(const R3Ray& ray) const { return R3Distance(*this, ray); } \
    RNLength __type:: Distance(const R3Span& span) const { return R3Distance(*this, span); } \
    RNLength __type:: Distance(const R3Plane& plane) const { return R3Distance(*this, plane); } \
    RNLength __type:: Distance(const R3Halfspace& halfspace) const { return R3Distance(*this, halfspace); } \
    RNLength __type:: Distance(const R3Box& box) const { return R3Distance(*this, box); } \
    RNLength __type:: Distance(const R3OrientedBox& box) const { return R3Distance(*this, box); } \
    RNLength __type:: Distance(const R3Sphere& sphere) const { return R3Distance(*this, sphere); } \
    RNLength __type:: Distance(const R3Shape& shape) const { return R3Distance(*this, shape); } \
    RNLength __type:: SignedDistance(const R3Plane& plane) const { return R3SignedDistance(plane, *this); } \
    RNBoolean __type:: Contains(const R3Point& point) const { return R3Contains(*this, point); } \
    RNBoolean __type:: Contains(const R3Line& line) const { return R3Contains(*this, line); } \
    RNBoolean __type:: Contains(const R3Ray& ray) const { return R3Contains(*this, ray); } \
    RNBoolean __type:: Contains(const R3Span& span) const { return R3Contains(*this, span); } \
    RNBoolean __type:: Contains(const R3Plane& plane) const { return R3Contains(*this, plane); } \
    RNBoolean __type:: Contains(const R3Halfspace& halfspace) const { return R3Contains(*this, halfspace); } \
    RNBoolean __type:: Contains(const R3Box& box) const { return R3Contains(*this, box); } \
    RNBoolean __type:: Contains(const R3OrientedBox& box) const { return R3Contains(*this, box); } \
    RNBoolean __type:: Contains(const R3Sphere& sphere) const { return R3Contains(*this, sphere); } \
    RNBoolean __type:: Contains(const R3Shape& shape) const { return R3Contains(*this, shape); } \
    RNBoolean __type:: Inside(const R3Point& point) const { return R3Inside(*this, point); } \
    RNBoolean __type:: Inside(const R3Line& line) const { return R3Inside(*this, line); } \
    RNBoolean __type:: Inside(const R3Ray& ray) const { return R3Inside(*this, ray); } \
    RNBoolean __type:: Inside(const R3Span& span) const { return R3Inside(*this, span); } \
    RNBoolean __type:: Inside(const R3Plane& plane) const { return R3Inside(*this, plane); } \
    RNBoolean __type:: Inside(const R3Halfspace& halfspace) const { return R3Inside(*this, halfspace); } \
    RNBoolean __type:: Inside(const R3Box& box) const { return R3Inside(*this, box); } \
    RNBoolean __type:: Inside(const R3OrientedBox& box) const { return R3Inside(*this, box); } \
    RNBoolean __type:: Inside(const R3Sphere& sphere) const { return R3Inside(*this, sphere); } \
    RNBoolean __type:: Inside(const R3Shape& shape) const { return R3Inside(*this, shape); } \
    RNClassID __type:: Intersects(const R3Point& point) const { return R3Intersects(*this, point); } \
    RNClassID __type:: Intersects(const R3Line& line) const { return R3Intersects(*this, line); } \
    RNClassID __type:: Intersects(const R3Ray& ray, R3Point *hit_point, R3Vector *hit_normal, RNScalar *hit_t) const \
        { return R3Intersects(*this, ray, hit_point, hit_normal, hit_t); } \
    RNClassID __type:: Intersects(const R3Span& span) const { return R3Intersects(*this, span); } \
    RNClassID __type:: Intersects(const R3Plane& plane) const { return R3Intersects(*this, plane); } \
    RNClassID __type:: Intersects(const R3Halfspace& halfspace) const { return R3Intersects(*this, halfspace); } \
    RNClassID __type:: Intersects(const R3Box& box) const { return R3Intersects(*this, box); } \
    RNClassID __type:: Intersects(const R3OrientedBox& box) const { return R3Intersects(*this, box); } \
    RNClassID __type:: Intersects(const R3Sphere& sphere) const { return R3Intersects(*this, sphere); } \
    RNClassID __type:: Intersects(const R3Shape& shape) const { return R3Intersects(*this, shape); }



// Relationship function definitions

R3_SHAPE_RELATIONSHIP_DEFINITIONS(R3Triangle);
R3_SHAPE_RELATIONSHIP_DEFINITIONS(R3TriangleArray);
R3_SHAPE_RELATIONSHIP_DEFINITIONS(R3Circle);
R3_SHAPE_RELATIONSHIP_DEFINITIONS(R3Ellipse);
R3_SHAPE_RELATIONSHIP_DEFINITIONS(R3Rectangle);
R3_SHAPE_RELATIONSHIP_DEFINITIONS(R3Box);
R3_SHAPE_RELATIONSHIP_DEFINITIONS(R3OrientedBox);
R3_SHAPE_RELATIONSHIP_DEFINITIONS(R3Sphere);
R3_SHAPE_RELATIONSHIP_DEFINITIONS(R3Ellipsoid);
R3_SHAPE_RELATIONSHIP_DEFINITIONS(R3Cylinder);
R3_SHAPE_RELATIONSHIP_DEFINITIONS(R3Cone);





