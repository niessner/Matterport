/* Source file for the GAPS ray class */



/* Include files */

#include "R2Shapes/R2Shapes.h"



/* Public variables */

const R2Ray R2null_ray(0.0, 0.0, 0.0, 0.0);
const R2Ray R2posx_ray(0.0, 0.0, 1.0, 0.0);
const R2Ray R2posy_ray(0.0, 0.0, 0.0, 1.0);
const R2Ray R2negx_ray(0.0, 0.0, -1.0, 0.0);
const R2Ray R2negy_ray(0.0, 0.0, 0.0, -1.0);



/* Public functions */

int R2InitRay()
{
    /* Return success */
    return TRUE;
}



void R2StopRay()
{
}



R2Ray::
R2Ray(void)
{
}



R2Ray::
R2Ray(const R2Ray& ray)
    : line(ray.line),
      start(ray.start)
{
}



R2Ray::
R2Ray(const R2Point& point, const R2Vector& vector, RNBoolean normalized)
    : line(point, vector, normalized),
      start(point)
{
}



R2Ray::
R2Ray(const R2Point& point1, const R2Point& point2)
    : line(point1, point2),
      start(point1)
{
}



R2Ray::
R2Ray(RNCoord x1, RNCoord y1, RNCoord x2, RNCoord y2)
    : line(x1, y1, x2, y2),
      start(x1, y1)
{
}



const R2Point R2Ray::
Point(RNScalar t) const
{
    // Return point along span
    return (Start() + Vector() * t);
}



const RNScalar R2Ray::
T(const R2Point& point) const
{
    // Return parametric value of closest point on ray
    RNScalar denom = Vector().Dot(Vector());
    if (RNIsZero(denom)) return 0.0;
    R2Vector topoint = point - Start();
    return (Vector().Dot(topoint) / denom);
}



void R2Ray::
Mirror (const R2Line& mirror)
{
    // Mirror ray across line
    line.Mirror(mirror);
    start.Mirror(mirror);
}



void R2Ray::
Project (const R2Line& target)
{
    // Project ray onto target line
    line.Project(target);
    start.Project(target);
}



void R2Ray::
Transform (const R2Transformation& transformation)
{
    // Transform line and point
    line.Transform(transformation);
    start.Transform(transformation);
}



void R2Ray::
InverseTransform (const R2Transformation& transformation)
{
    // Transform line and point
    line.InverseTransform(transformation);
    start.InverseTransform(transformation);
}














