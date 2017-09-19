/* Source file for the R2 circle class */



/* Include files */

#include "R2Shapes/R2Shapes.h"



/* Public variables */

const R2Circle R2null_circle(R2Point(0.0, 0.0), -1.0);
const R2Circle R2zero_circle(R2Point(0.0, 0.0), 0.0);
const R2Circle R2unit_circle(R2Point(0.0, 0.0), 1.0);
const R2Circle R2infinite_circle(R2Point(0.0, 0.0), RN_INFINITY);



/* Unit circle coordinates */

const int R2circle_npoints = 32;
RNAngle R2circle_angles[R2circle_npoints];
R2Point R2circle_points[R2circle_npoints];



/* Public functions */

int 
R2InitCircle()
{
    // Initialize circle angles and points 
    for (int i = 0; i < R2circle_npoints; i++) {
        R2circle_angles[i] = ((RNScalar) i * RN_TWO_PI) / (RNScalar) R2circle_npoints;
        R2circle_points[i] = R2Point(cos(R2circle_angles[i]), sin(R2circle_angles[i]));
    }

    /* Return success */
    return TRUE;
}



void 
R2StopCircle()
{
}



R2Circle::
R2Circle(void)
{
}



R2Circle::
R2Circle(const R2Circle& circle)
    : center(circle.center),
      radius(circle.radius)
{
}



R2Circle::
R2Circle(const R2Point& center, RNLength radius)
    : center(center),
      radius(radius)
{
}



const RNBoolean R2Circle::
IsPoint(void) const
{
    // A circle only lies on a single point if it has radius zero
    return RNIsZero(radius);
}



const RNBoolean R2Circle::
IsLinear(void) const
{
    // A circle only lies within a line if it is a point
    return IsPoint();
}



const RNBoolean R2Circle::
IsConvex(void) const
{
    // All circles are convex
    return TRUE;
}



const RNArea R2Circle::
Area(void) const
{
    // Return area of circle
    return (RN_PI * radius * radius);
}



const R2Point R2Circle::
Centroid(void) const
{
    // Return centroid of circle
    return center;
}



const R2Shape& R2Circle::
BShape(void) const
{
    // Return self
    return *this;
}



const R2Box R2Circle::
BBox(void) const
{
    // Return bounding box of circle
    return R2Box(center.X() - radius, center.Y() - radius,
                 center.X() + radius, center.Y() + radius);
}



const R2Circle R2Circle::
BCircle(void) const
{
    // Return self
    return *this;
}



void R2Circle::
Empty(void)
{
    // Empty circle
    *this = R2null_circle;
}



void R2Circle::
Translate(const R2Vector& vector)
{
    // Move circle center
    center.Translate(vector);
}



void R2Circle::
Reposition(const R2Point& center)
{
    // Set circle center
    this->center = center;
}



void R2Circle::
Resize(RNLength radius) 
{
    // Set circle radius
    this->radius = radius;
}



void R2Circle::
Transform (const R2Transformation& transformation)
{
    // Transform center 
    center.Transform(transformation);

    // Scale radius 
    if (!transformation.IsIsotropic()) RNWarning("Circle transformed by anisotropic transformation");
    R2Vector v(R2ones_vector);
    v.Transform(transformation);
    radius *= v.Length() / RN_SQRT_TWO;
}





