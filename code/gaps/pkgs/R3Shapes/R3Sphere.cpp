/* Source file for the R3 sphere class */



/* Include files */

#include "R3Shapes/R3Shapes.h"



/* Public variables */

const R3Sphere R3null_sphere(R3Point(0.0, 0.0, 0.0), -1.0);
const R3Sphere R3zero_sphere(R3Point(0.0, 0.0, 0.0), 0.0);
const R3Sphere R3unit_sphere(R3Point(0.0, 0.0, 0.0), 1.0);
const R3Sphere R3infinite_sphere(R3Point(0.0, 0.0, 0.0), RN_INFINITY);



/* Class type definitions */

RN_CLASS_TYPE_DEFINITIONS(R3Sphere);



/* Public functions */

int 
R3InitSphere()
{
    /* Return success */
    return TRUE;
}



void 
R3StopSphere()
{
}



R3Sphere::
R3Sphere(void)
{
}



R3Sphere::
R3Sphere(const R3Sphere& sphere)
    : center(sphere.center),
      radius(sphere.radius)
{
}



R3Sphere::
R3Sphere(const R3Point& center, RNLength radius)
    : center(center),
      radius(radius)
{
}



const RNBoolean R3Sphere::
IsPoint(void) const
{
    // A sphere only lies on a single point if it has radius zero
    return RNIsZero(radius);
}



const RNBoolean R3Sphere::
IsLinear(void) const
{
    // A sphere only lies within a line if it is a point
    return this->IsPoint();
}



const RNBoolean R3Sphere::
IsPlanar(void) const
{
    // A sphere only lies within a plane if it is a point
    return this->IsPoint();
}



const RNBoolean R3Sphere::
IsConvex(void) const
{
    // All spheres are convex
    return TRUE;
}



const RNInterval R3Sphere::
NFacets(void) const
{
    // Return number of facets (polygons)
    return RNInterval(128.0, 128.0);
}



const RNArea R3Sphere::
Area(void) const
{
    // Return surface area of sphere
    return (4.0 * RN_PI * radius * radius);
}



const RNVolume R3Sphere::
Volume(void) const
{
    // Return volume of sphere
    return (1.3333333333333 * RN_PI * radius * radius * radius);
}



const R3Point R3Sphere::
Centroid(void) const
{
    // Return centroid of sphere
    return center;
}



const R3Point R3Sphere::
ClosestPoint(const R3Point& point) const
{
    // Return closest point in sphere
    if (radius <= 0) return center;
    R3Vector v = point - Center();
    RNLength d = v.Length();
    if (d < radius) return point;
    else return center + radius * v / d;
}



const R3Point R3Sphere::
FurthestPoint(const R3Point& point) const
{
    // Return furthest point in sphere
    if (radius <= 0) return center;
    R3Vector v = point - Center();
    RNLength d = v.Length();
    if (RNIsZero(d)) return R3Point(center[0] + radius, center[1], center[2]);
    else return center - radius * v / d;
}



const R3Shape& R3Sphere::
BShape(void) const
{
    // Return self
    return *this;
}



const R3Box R3Sphere::
BBox(void) const
{
    // Return bounding box of sphere
    return R3Box(center.X() - radius, center.Y() - radius, center.Z() - radius,
                 center.X() + radius, center.Y() + radius, center.Z() + radius);
}



const R3Sphere R3Sphere::
BSphere(void) const
{
    // Return self
    return *this;
}



void R3Sphere::
Empty(void)
{
    // Empty sphere
    *this = R3null_sphere;
}



void R3Sphere::
Translate(const R3Vector& vector)
{
    // Move sphere center
    center.Translate(vector);
}



void R3Sphere::
Reposition(const R3Point& center)
{
    // Set sphere center
    this->center = center;
}



void R3Sphere::
Resize(RNLength radius) 
{
    // Set sphere radius
    this->radius = radius;
}




void R3Sphere::
Transform (const R3Transformation& transformation)
{
    // Transform center 
    center.Transform(transformation);

    // Scale radius 
    if (!transformation.IsIsotropic()) RNWarning("Sphere transformed by anisotropic transformation");
    radius *= transformation.ScaleFactor();
}




