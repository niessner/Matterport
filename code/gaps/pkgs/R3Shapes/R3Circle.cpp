/* Source file for the R3 circle class */



/* Include files */

#include "R3Shapes/R3Shapes.h"



/* Public variables */

const R3Circle R3null_circle(R3Point(0.0, 0.0, 0.0), -1.0, R3Vector(0.0, 0.0, 0.0));
const R3Circle R3zero_circle(R3Point(0.0, 0.0, 0.0), 0.0, R3Vector(0.0, 0.0, 1.0));
const R3Circle R3unit_circle(R3Point(0.0, 0.0, 0.0), 1.0, R3Vector(0.0, 0.0, 1.0));
const R3Circle R3infinite_circle(R3Point(0.0, 0.0, 0.0), RN_INFINITY, R3Vector(0.0, 0.0, 1.0));



/* Unit circle coordinates */

const int R3circle_npoints = 32;
RNAngle R3circle_angles[R3circle_npoints];
R3Point R3circle_points[R3circle_npoints];
R2Point R3circle_texcoords[R3circle_npoints];



/* Class type definitions */

RN_CLASS_TYPE_DEFINITIONS(R3Circle);



/* Public functions */

int 
R3InitCircle()
{
    // Initialize circle angles and points 
    for (int i = 0; i < R3circle_npoints; i++) {
        R3circle_angles[i] = ((RNScalar) i * RN_TWO_PI) / (RNScalar) R3circle_npoints;
        R3circle_points[i] = R3Point(cos(R3circle_angles[i]), sin(R3circle_angles[i]), 0.0);
        R3circle_texcoords[i] = R2Point(0.5 * (R3circle_points[i].X() + 1.0), 0.5 * (R3circle_points[i].Y() + 1.0));
    }

    /* Return success */
    return TRUE;
}



void 
R3StopCircle()
{
}



R3Circle::
R3Circle(void)
{
}



R3Circle::
R3Circle(const R3Circle& circle)
    : center(circle.center),
      radius(circle.radius),
      plane(circle.plane)
{
}



R3Circle::
R3Circle(const R3Point& center, RNLength radius, const R3Vector& normal)
    : center(center),
      radius(radius),
      plane(center, normal)
{
}



R3Circle::
R3Circle(const R3Point& p1, const R3Point& p2, const R3Point& p3)
    : center(R3null_point),
      radius(0.0),
      plane(p1, p2, p3)
{
    RNAbort("Not Implemented");
}



const RNBoolean R3Circle::
IsPoint(void) const
{
    // A circle only lies on a single point if it has radius zero
    return (radius == 0.0);
}



const RNBoolean R3Circle::
IsLinear(void) const
{
    // A circle only lies within a line if it is a point
    return IsPoint();
}



const RNBoolean R3Circle::
IsPlanar(void) const
{
    // All circles are planar
    return TRUE;
}



const RNBoolean R3Circle::
IsConvex(void) const
{
    // All circles are convex
    return TRUE;
}



const RNInterval R3Circle::
NFacets(void) const
{
    // Return number of facets (polygons)
    return RNInterval(1.0, 1.0);
}



const RNLength R3Circle::
Length(void) const
{
    // Return circumference of circle
    return (2.0 * RN_PI * radius);
}



const RNArea R3Circle::
Area(void) const
{
    // Return area of circle
    return (RN_PI * radius * radius);
}



const R3Point R3Circle::
Centroid(void) const
{
    // Return centroid of circle
    return center;
}


const R3Shape& R3Circle::
BShape(void) const
{
    // Return self ???
    return *this;
}



const R3Box R3Circle::
BBox(void) const
{
    // Return bounding box of circle ???
    return R3Box(center.X() - radius, center.Y() - radius, center.Z() - radius,
                 center.X() + radius, center.Y() + radius, center.Z() + radius);
}



const R3Sphere R3Circle::
BSphere(void) const
{
    // Return bounding sphere of circle
    return R3Sphere(center, radius);
}



void R3Circle::
Flip (void) 
{
    // Flip orientation of plane
    plane.Flip();
}



void R3Circle::
Reposition (const R3Point& center) 
{
    // Reposition circle's center
    this->center = center;
    plane.Reposition(center);
}



void R3Circle::
Translate (const R3Vector& offset) 
{
    // Translate circle
    center.Translate(offset);
    plane.Translate(offset);
}



void R3Circle::
Align (const R3Vector& normal) 
{
    // Rotate around center to align with normal
    plane.Reset(center, normal);
}



void R3Circle::
Resize (RNScalar radius)
{
    // Set radius
    this->radius = radius;
}



void R3Circle::
Transform (const R3Transformation& transformation)
{
    // Transform circle ???
    center.Transform(transformation);
    plane.Transform(transformation);
    if (!transformation.IsIsotropic()) RNWarning("Circle transformed by anisotropic transformation");
    radius *= transformation.ScaleFactor();
}



void R3Circle::
Reset (const R3Point& center, RNScalar radius, const R3Vector& normal) 
{
    // Rotate around center to align with normal
    this->center = center;
    this->radius = radius;
    plane.Reset(center, normal);
}





