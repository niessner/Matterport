/* Source file for the R3 cylinder class */



/* Include files */

#include "R3Shapes/R3Shapes.h"



/* Public variables */

const R3Cylinder R3null_cylinder(R3Point(0.0, 0.0, 0.0), R3Point(0.0, 0.0, 0.0), -1.0);
const R3Cylinder R3zero_cylinder(R3Point(0.0, 0.0, 0.0), R3Point(0.0, 0.0, 0.0), 0.0);
const R3Cylinder R3unit_cylinder(R3Point(0.0, 0.0, -1.0), R3Point(0.0, 0.0, 1.0), 1.0);
const R3Cylinder R3infinite_cylinder(R3Point(0.0, 0.0, -RN_INFINITY), R3Point(0.0, 0.0, RN_INFINITY), RN_INFINITY);



/* Class type definitions */

RN_CLASS_TYPE_DEFINITIONS(R3Cylinder);



/* Public functions */

int 
R3InitCylinder()
{
    /* Return success */
    return TRUE;
}



void 
R3StopCylinder()
{
}



R3Cylinder::
R3Cylinder(void)
{
}



R3Cylinder::
R3Cylinder(const R3Cylinder& cylinder)
    : axis(cylinder.axis),
      base(cylinder.base),
      top(cylinder.top)
{
}



R3Cylinder::
R3Cylinder(const R3Span& axis, RNLength radius)
    : axis(axis),
      base(axis.Start(), radius, -(axis.Vector())),
      top(axis.End(), radius, axis.Vector())
{
}



R3Cylinder::
R3Cylinder(const R3Point& p1, const R3Point& p2, RNLength radius)
    : axis(p1, p2),
      base(axis.Start(), radius, -(axis.Vector())),
      top(axis.End(), radius, axis.Vector())
{
}



const RNBoolean R3Cylinder::
IsPoint(void) const
{
    // A cylinder lies on a single point if axis is point and it has radius zero
    return ((Radius() == 0.0) && (Height() == 0.0));
}



const RNBoolean R3Cylinder::
IsLinear(void) const
{
    // A cylinder only lies within a line if it has zero radius
    return (Radius() == 0.0);
}



const RNBoolean R3Cylinder::
IsPlanar(void) const
{
    // A cylinder only lies within a plane if it has zero height or zero radius
    return ((Radius() == 0.0) || (Height() == 0.0));
}



const RNBoolean R3Cylinder::
IsConvex(void) const
{
    // All cylinders are convex
    return TRUE;
}



RNBoolean R3Cylinder::
operator==(const R3Cylinder& cylinder) const
{
    // Return whether cylinder is equal
    return ((axis == cylinder.axis) && 
	    (base == cylinder.base) &&
	    (top == cylinder.top));
}



const RNInterval R3Cylinder::
NFacets(void) const
{
    // Return number of facets (polygons)
    // Min is 8 sides, max is 32 sides (plus top and base)
    return RNInterval(10.0, 34.0);
}



const RNArea R3Cylinder::
Area(void) const
{
    // Return surface area of cylinder
    return (2.0 * RN_PI * Radius() * (Radius() + Height()));
}



const RNVolume R3Cylinder::
Volume(void) const
{
    // Return volume of cylinder
    return (RN_PI * Radius() * Radius() * Height());
}



const R3Point R3Cylinder::
Centroid(void) const
{
    // Return centroid of cylinder
    return axis.Centroid();
}



const R3Shape& R3Cylinder::
BShape(void) const
{
    // Return self ???
    return *this;
}



const R3Box R3Cylinder::
BBox(void) const
{
    // Return bounding box of cylinder ???
    R3Box bbox(R3null_box);
    bbox.Union(base.BBox());
    bbox.Union(top.BBox());
    return bbox;
}



const R3Sphere R3Cylinder::
BSphere(void) const
{
    // Return BSphere
    RNScalar a = 0.5 * Height();
    RNScalar r = sqrt(a * a + Radius() * Radius());
    return R3Sphere(Centroid(), r);
}



void R3Cylinder::
Translate(const R3Vector& vector)
{
    // Move cylinder axis
    axis.Translate(vector);
    base.Translate(vector);
    top.Translate(vector);
}



void R3Cylinder::
Reposition(const R3Span& axis)
{
    // Set cylinder axis
    this->axis = axis;
    base.Reposition(axis.Start());
    base.Align(-(axis.Vector()));
    top.Reposition(axis.End());
    top.Align(axis.Vector());
}



void R3Cylinder::
Resize(RNLength radius) 
{
    // Set cylinder radius
    base.Resize(radius);
    top.Resize(radius);
}



void R3Cylinder::
Transform (const R3Transformation& transformation)
{
    // Transform 
    axis.Transform(transformation);
    base.Transform(transformation);
    top.Transform(transformation);
}



void R3Cylinder::
Empty(void)
{
    // Empty cylinder
    *this = R3null_cylinder;
}




