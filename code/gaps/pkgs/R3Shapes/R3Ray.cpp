/* Source file for the GAPS ray class */



/* Include files */

#include "R3Shapes/R3Shapes.h"



/* Public variables */

const R3Ray R3null_ray(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
const R3Ray R3posx_ray(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
const R3Ray R3posy_ray(0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
const R3Ray R3posz_ray(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
const R3Ray R3negx_ray(0.0, 0.0, 0.0, -1.0, 0.0, 0.0);
const R3Ray R3negy_ray(0.0, 0.0, 0.0, 0.0, -1.0, 0.0);
const R3Ray R3negz_ray(0.0, 0.0, 0.0, 0.0, 0.0, -1.0);



/* Public functions */

int R3InitRay()
{
    /* Return success */
    return TRUE;
}



void R3StopRay()
{
}



R3Ray::
R3Ray(void)
{
}



R3Ray::
R3Ray(const R3Ray& ray)
    : line(ray.line)
{
}



R3Ray::
R3Ray(const R3Point& point, const R3Vector& vector, RNBoolean normalized)
    : line(point, vector, normalized)
{
}



R3Ray::
R3Ray(const R3Point& point1, const R3Point& point2)
    : line(point1, point2)
{
}



R3Ray::
R3Ray(RNCoord x1, RNCoord y1, RNCoord z1, RNCoord x2, RNCoord y2, RNCoord z2)
    : line(x1, y1, z1, x2, y2, z2)
{
}



const R3Point R3Ray::
Point(RNScalar t) const
{
    // Return point along span
    return (Start() + Vector() * t);
}



const RNScalar R3Ray::
T(const R3Point& point) const
{
    // Return parametric value of closest point on ray
    if (IsZero()) return 0.0;
    RNScalar denom = Vector().Dot(Vector());
    assert(RNIsNotZero(denom));
    R3Vector topoint = point - Start();
    return (Vector().Dot(topoint) / denom);
}



void R3Ray::
Transform (const R3Transformation& transformation)
{
    // Transform line
    line.Transform(transformation);
}



void R3Ray::
InverseTransform (const R3Transformation& transformation)
{
    // Transform line
    line.InverseTransform(transformation);
}





