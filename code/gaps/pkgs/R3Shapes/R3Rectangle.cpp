/* Source file for the R3 rectangle class */



/* Include files */

#include "R3Shapes/R3Shapes.h"



/* Public variables */

const R3Rectangle R3null_rectangle(R3Point(0.0, 0.0, 0.0), R3Vector(1.0, 0.0, 0.0), R3Vector(0.0, 1.0, 0.0), -1.0, -1.0);
const R3Rectangle R3zero_rectangle(R3Point(0.0, 0.0, 0.0), R3Vector(1.0, 0.0, 0.0), R3Vector(0.0, 1.0, 0.0),  0.0,  0.0);



/* Class type definitions */

RN_CLASS_TYPE_DEFINITIONS(R3Rectangle);



/* Public functions */

int 
R3InitRectangle()
{
    /* Return success */
    return TRUE;
}



void 
R3StopRectangle()
{
}



R3Rectangle::
R3Rectangle(void)
{
}



R3Rectangle::
R3Rectangle(const R3Rectangle& rectangle)
    : cs(rectangle.cs)
{
    // Copy radii from rectangle
    this->radius[0] = rectangle.radius[0];
    this->radius[1] = rectangle.radius[1];
}



R3Rectangle::
R3Rectangle(const R3Point& center, const R3Vector& axis0, const R3Vector& axis1)
{
    // Initialize everything
    Reset(center, axis0, axis1, axis0.Length(), axis1.Length());
}



R3Rectangle::
R3Rectangle(const R3Point& center, 
    const R3Vector& axis0, const R3Vector& axis1, 
    RNLength radius0, RNLength radius1)
{
    // Initialize everything
    Reset(center, axis0, axis1, radius0, radius1);
}



const R3Point R3Rectangle::
ClosestPoint(const R3Point& point) const
{
    // Check if rectangle is empty
    if (IsEmpty()) return R3zero_point;

    // Start with original point
    R3Point closest_point = point;

    // Project onto plane
    closest_point.Project(Plane());

    // Iteratively project onto each edge, if outside
    for (int sign = -1; sign <= 1; sign += 2) {
      for (int dim = 0; dim <= 1; dim++) {
        R3Plane plane(Center() + sign * Radius(dim) * Axis(dim), sign * Axis(dim));
        RNScalar d = R3SignedDistance(plane, closest_point);
        if (d > 0) closest_point -= d * plane.Normal();
      }
    }

    // Return closest point
    return closest_point;
}



const R3Point R3Rectangle::
FurthestPoint(const R3Point& point) const
{
    // Check if rectangle is empty
    if (IsEmpty()) return R3zero_point;

    // Start with original point
    R3Point closest_point = point;

    // Project onto plane
    closest_point.Project(Plane());

    // Iteratively project onto each edge, if inside
    for (int sign = -1; sign <= 1; sign += 2) {
      for (int dim = 0; dim <= 1; dim++) {
        R3Plane plane(Center() + sign * Radius(dim) * Axis(dim), sign * Axis(dim));
        RNScalar d = R3SignedDistance(plane, closest_point);
        if (d < 0) closest_point -= d * plane.Normal();
      }
    }

    // Return closest point
    return closest_point;
}



const RNBoolean R3Rectangle::
IsPoint(void) const
{
    // A rectangle lies on a single point if its maximum radius is zero
    return (RNIsZero(MaxRadius()));
}



const RNBoolean R3Rectangle::
IsLinear(void) const
{
    // Check if rectangle is empty
    if (IsEmpty()) return FALSE;

    // A rectangle lies on a line if its smallest radius is zero
    return RNIsZero(MinRadius());
}



const RNBoolean R3Rectangle::
IsPlanar(void) const
{
    // Check if rectangle is empty
    if (IsEmpty()) return FALSE;

    // All rectangles are planar
    return TRUE;
}



const RNBoolean R3Rectangle::
IsConvex(void) const
{
    // Check if rectangle is empty
    if (IsEmpty()) return FALSE;

    // All rectangles are convex
    return TRUE;
}



RNBoolean R3Rectangle::
operator==(const R3Rectangle& rectangle) const
{
    // Return whether rectangle is equal
    return ((Center() == Center()) &&
            (Axis(0) == rectangle.Axis(0)) && 
	    (Axis(1) == rectangle.Axis(1)) && 
	    (Radius(0) == rectangle.Radius(0)) && 
	    (Radius(1) == rectangle.Radius(1)));
}



const RNInterval R3Rectangle::
NFacets(void) const
{
    // Check if rectangle is empty
    if (IsEmpty()) return RNInterval(0, 0);

    // Return number of facets (polygons)
    return RNInterval(1, 1);
}



const RNArea R3Rectangle::
Area(void) const
{
    // Check if rectangle is empty
    if (IsEmpty()) return 0.0;

    // Return surface area of rectangle
    return 4*Radius(0)*Radius(1);
}



const R3Point R3Rectangle::
Centroid(void) const
{
    // Check if rectangle is empty
    if (IsEmpty()) return R3zero_point;

    // Return centroid of rectangle
    return Center();
}



const R3Shape& R3Rectangle::
BShape(void) const
{
    // Return this rectangle
    return *this;
}



const R3Box R3Rectangle::
BBox(void) const
{
    // Check if rectangle is empty
    if (IsEmpty()) return R3null_box;

    // Return bounding box 
    R3Box bbox = R3null_box;
    for (int i = 0; i < RN_NUM_QUADRANTS; i++) 
      bbox.Union(Corner(i));
    return bbox;
}



const R3Sphere R3Rectangle::
BSphere(void) const
{
    // Check if rectangle is empty
    if (IsEmpty()) return R3null_sphere;

    // Return BSphere
    return R3Sphere(Center(), DiagonalRadius());
}



void R3Rectangle::
Translate(const R3Vector& vector)
{
    // Check if rectangle is empty
    if (IsEmpty()) return;

    // Move rectangle coordinate system
    cs.Translate(vector);
}



void R3Rectangle::
Reposition(const R3Point& center)
{
    // Check if rectangle is empty
    if (IsEmpty()) return;

    // Set rectangle center
    cs.SetOrigin(center);
}



void R3Rectangle::
Resize(RNLength radius0, RNLength radius1) 
{
    // Check if rectangle is empty
    if (IsEmpty()) return;

    // Set rectangle radius
    this->radius[0] = radius0;
    this->radius[1] = radius1;
}



void R3Rectangle::
Reorient(const R3Vector& axis0, const R3Vector& axis1) 
{
    // Check if rectangle is empty
    if (IsEmpty()) return;

    // Set rectangle orientation
    Reset(cs.Origin(), axis0, axis1, radius[0], radius[1]);
}



void R3Rectangle::
Transform (const R3Transformation& transformation)
{
    // Check if rectangle is empty
    if (IsEmpty()) return;

    // Transform 
    cs.Transform(transformation);
    radius[0] *= transformation.ScaleFactor();
    radius[1] *= transformation.ScaleFactor();
}



void R3Rectangle::
Empty(void)
{
    // Empty rectangle
    *this = R3null_rectangle;
}



void R3Rectangle::
Reset(const R3Point& center, 
    const R3Vector& axis0, const R3Vector& axis1, 
    RNLength radius0, RNLength radius1)
{
    // Check axes
    assert(RNIsPositive(axis0.Length()));
    assert(RNIsPositive(axis1.Length()));

    // Normalize and orthogonalize axes
    R3Vector axis[3];
    axis[0] = axis0;
    axis[0].Normalize();
    axis[1] = axis1;
    axis[1].Normalize();
    axis[1] = axis[0] % axis[1];
    axis[1] = axis[1] % axis[0];
    axis[1].Normalize();
    axis[2] = axis[0] % axis[1];
    axis[2].Normalize();

    // Assign radii
    radius[0] = radius0;
    radius[1] = radius1;

    // Sort axes
    R3Vector swap_axis;
    RNScalar swap_radius;
    if (radius[1] > radius[0]) {
      swap_axis = axis[0]; swap_radius = radius[0];
      axis[0] = axis[1]; radius[0] = radius[1];
      axis[1] = swap_axis; radius[1] = swap_radius;
    }

    // Ensure right-handed coordinate system
    axis[2] = axis[0] % axis[1];
    axis[2].Normalize();

    // Set coordinate system
    cs.Reset(center, R3Triad(axis[0], axis[1], axis[2]));
}


    
    










