/* Source file for GAPS parallel utility */



/* Include files */

#include "R3Shapes/R3Shapes.h"



RNBoolean R3Parallel(const R3Vector& vector1, const R3Vector& vector2)
{
    // Return whether two vectors are coincident (or anti-coincident)
    if ((RNIsEqual(vector1.X(), vector2.X())) &&
        (RNIsEqual(vector1.Y(), vector2.Y())) && 
        (RNIsEqual(vector1.Z(), vector2.Z()))) return TRUE;
    if ((RNIsEqual(vector1.X(), -vector2.X())) &&
        (RNIsEqual(vector1.Y(), -vector2.Y())) && 
        (RNIsEqual(vector1.Z(), -vector2.Z()))) return TRUE;
    return FALSE;
}



RNBoolean R3Parallel(const R3Vector& vector, const R3Line& line)
{
    // Return whether vector and line are parallel
    return R3Parallel(vector, line.Vector());
}



RNBoolean R3Parallel(const R3Vector& vector, const R3Ray& ray)
{
    // Return whether vector and ray are parallel
    return R3Parallel(vector, ray.Vector());
}



RNBoolean R3Parallel(const R3Vector& vector, const R3Span& span)
{
    // Return whether vector and span are parallel
    return R3Parallel(vector, span.Vector());
}



RNBoolean R3Parallel(const R3Vector& vector, const R3Plane& plane)
{
    // Return whether vector and plane are parallel
    return R3Perpendicular(vector, plane.Normal());
}



RNBoolean R3Parallel(const R3Line& line1, const R3Line& line2)
{
    // Return whether line1 and line2 are parallel
    return R3Parallel(line1.Vector(), line2.Vector());
}



RNBoolean R3Parallel(const R3Line& line, const R3Ray& ray)
{
    // Return whether line and ray are parallel
    return R3Parallel(line.Vector(), ray.Vector());
}



RNBoolean R3Parallel(const R3Line& line, const R3Span& span)
{
    // Return whether line and span are parallel
    return R3Parallel(line.Vector(), span.Vector());
}



RNBoolean R3Parallel(const R3Line& line, const R3Plane& plane)
{
    // Return whether line and plane are parallel
    return R3Perpendicular(line.Vector(), plane.Normal());
}



RNBoolean R3Parallel(const R3Ray& ray1, const R3Ray& ray2)
{
    // Return whether ray1 and ray2 are parallel
    return R3Parallel(ray1.Vector(), ray2.Vector());
}



RNBoolean R3Parallel(const R3Ray& ray, const R3Span& span)
{
    // Return whether ray and span are parallel
    return R3Parallel(ray.Vector(), span.Vector());
}



RNBoolean R3Parallel(const R3Ray& ray, const R3Plane& plane)
{
    // Return whether ray and plane are parallel
    return R3Perpendicular(ray.Vector(), plane.Normal());
}



RNBoolean R3Parallel(const R3Span& span1, const R3Span& span2)
{
    // Return whether span1 and span2 are parallel
    return R3Parallel(span1.Vector(), span2.Vector());
}



RNBoolean R3Parallel(const R3Span& span, const R3Plane& plane)
{
    // Return whether span and plane are parallel
    return R3Perpendicular(span.Vector(), plane.Normal());
}



RNBoolean R3Parallel(const R3Plane& plane1, const R3Plane& plane2)
{
    // Return whether plane1 and plane2 are parallel
    return R3Parallel(plane1.Normal(), plane2.Normal());
}




