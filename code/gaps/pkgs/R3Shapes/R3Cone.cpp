/* Source file for the R3 cone class */



/* Include files */

#include "R3Shapes/R3Shapes.h"



/* Public variables */

const R3Cone R3null_cone(R3Point(0.0, 0.0, 0.0), R3Point(0.0, 0.0, 0.0), -1.0);
const R3Cone R3zero_cone(R3Point(0.0, 0.0, 0.0), R3Point(0.0, 0.0, 0.0), 0.0);
const R3Cone R3unit_cone(R3Point(0.0, 0.0, -1.0), R3Point(0.0, 0.0, 1.0), 1.0);
const R3Cone R3infinite_cone(R3Point(0.0, 0.0, -RN_INFINITY), R3Point(0.0, 0.0, RN_INFINITY), RN_INFINITY);



/* Class type definitions */

RN_CLASS_TYPE_DEFINITIONS(R3Cone);



/* Public functions */

int 
R3InitCone()
{
    /* Return success */
    return TRUE;
}



void 
R3StopCone()
{
}



R3Cone::
R3Cone(void)
{
}



R3Cone::
R3Cone(const R3Cone& cone)
    : axis(cone.axis),
      base(cone.base)
{
}



R3Cone::
R3Cone(const R3Span& axis, RNLength radius)
    : axis(axis),
      base(axis.Start(), radius, -(axis.Vector()))
{
}



R3Cone::
R3Cone(const R3Point& p1, const R3Point& p2, RNLength radius)
    : axis(p1, p2),
      base(p1, radius, -(axis.Vector()))
{
}



const RNBoolean R3Cone::
IsPoint(void) const
{
    // A cone lies on a single point if axis is point and it has radius zero
    return ((Height() == 0.0) && (Radius() == 0.0));
}



const RNBoolean R3Cone::
IsLinear(void) const
{
    // A cone only lies within a line if it has zero radius
    return (Radius() == 0.0);
}



const RNBoolean R3Cone::
IsPlanar(void) const
{
    // A cone only lies within a plane if it has zero height or zero radius
    return ((Height() == 0.0) || (Radius() == 0.0));
}



const RNBoolean R3Cone::
IsConvex(void) const
{
    // All cones are convex
    return TRUE;
}



const RNInterval R3Cone::
NFacets(void) const
{
    // Return number of facets (polygons)
    // Min is 8 sides, max is 32 sides (plus base)
    return RNInterval(9.0, 33.0);
}



const RNArea R3Cone::
Area(void) const
{
    // Return surface area of cone ???
    RNLength s = sqrt((Height() * Height()) + (Radius() * Radius()));
    return (RN_PI * Radius() * (Radius() + s));
}



const RNVolume R3Cone::
Volume(void) const
{
    // Return volume of cone
    return (0.33333333 * RN_PI * Radius() * Radius() * Height());
}



const R3Point R3Cone::
Centroid(void) const
{
    // Return centroid of cone ???
    return axis.Point((RNScalar) 0.33333333);
}



const R3Shape& R3Cone::
BShape(void) const
{
    // Return self ???
    return *this;
}



const R3Box R3Cone::
BBox(void) const
{
    // Return bounding box 
    R3Box bbox(Apex(), Apex());
    bbox.Union(base.BBox());
    return bbox;
}



const R3Sphere R3Cone::
BSphere(void) const
{
    // Compute distance up axis for sphere center using quadradic formula
    if (RNIsLessOrEqual(Height(), Radius())) return base.BSphere();
    RNScalar d = 1.0 - 4.0 * (Radius() * Radius() - Height());
    if (RNIsNegativeOrZero(d)) return base.BSphere();
    d = sqrt(d);
    if (d <= 1.0) d = 0.5 * (1.0 + d);
    else  d = 0.5 * (1.0 + d);
    R3Point center = base.Center() + d * axis.Vector();
    return R3Sphere(center, Height() - d); 
}



void R3Cone::
Translate(const R3Vector& vector)
{
    // Move cone axis
    axis.Translate(vector);
    base.Translate(vector);
}



void R3Cone::
Reposition(const R3Span& axis)
{
    // Set cone axis
    this->axis = axis;
    base.Reposition(axis.Start());
    base.Align(-axis.Vector());
}



void R3Cone::
Resize(RNLength radius) 
{
    // Set cone radius
    base.Resize(radius);
}




void R3Cone::
Transform (const R3Transformation& transformation)
{
    // Transform 
    axis.Transform(transformation);
    base.Transform(transformation);
}



void R3Cone::
Empty(void)
{
    // Empty cone
    *this = R3null_cone;
}




