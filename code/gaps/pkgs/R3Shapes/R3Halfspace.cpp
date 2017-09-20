/* Source file for the GAPS halfspace class */



/* Include files */

#include "R3Shapes/R3Shapes.h"



/* Public variables */

const R3Halfspace R3null_halfspace(0.0, 0.0, 0.0, 0.0);
const R3Halfspace R3posx_halfspace(1.0, 0.0, 0.0, 0.0);
const R3Halfspace R3posy_halfspace(0.0, 1.0, 0.0, 0.0);
const R3Halfspace R3posz_halfspace(0.0, 0.0, 1.0, 0.0);
const R3Halfspace R3negx_halfspace(-1.0, 0.0, 0.0, 0.0);
const R3Halfspace R3negy_halfspace(0.0, -1.0, 0.0, 0.0);
const R3Halfspace R3negz_halfspace(0.0, 0.0, -1.0, 0.0);



/* Public functions */

int 
R3InitHalfspace()
{
    /* Return success */
    return TRUE;
}



void 
R3StopHalfspace()
{
}



R3Halfspace::
R3Halfspace (void)
{
}



R3Halfspace::
R3Halfspace (const R3Halfspace& halfspace)
    : plane(halfspace.plane)
{
}



R3Halfspace::
R3Halfspace(RNScalar a, RNScalar b, RNScalar c, RNScalar d)
    : plane(a, b, c, d)
{
}



R3Halfspace::
R3Halfspace(const RNScalar a[4])
    : plane(a)
{
}



R3Halfspace::
R3Halfspace(const R3Vector& normal, RNScalar d)
    : plane(normal, d)
{
}



R3Halfspace::
R3Halfspace(const R3Point& point, const R3Vector& normal)
    : plane(point, normal)
{
}



R3Halfspace::
R3Halfspace(const R3Point& point, const R3Line& line)
    : plane(point, line)
{
}



R3Halfspace::
R3Halfspace(const R3Point& point, const R3Vector& vector1, const R3Vector& vector2)
    : plane(point, vector1, vector2)
{
}



R3Halfspace::
R3Halfspace(const R3Point& point1, const R3Point& point2, const R3Point& point3)
    : plane(point1, point2, point3)
{
}



R3Halfspace::
R3Halfspace(const R3Plane& plane, int )
    : plane(plane)
{
}



void R3Halfspace::
Mirror (const R3Plane& plane)
{
    // Mirror halfspace with respect to plane
    this->plane.Mirror(plane);
}



void R3Halfspace::
Transform (const R3Transformation& transformation)
{
    // Transform halfspace
    plane.Transform(transformation);
}



void R3Halfspace::
InverseTransform (const R3Transformation& transformation)
{
    // Transform halfspace
    plane.InverseTransform(transformation);
}




