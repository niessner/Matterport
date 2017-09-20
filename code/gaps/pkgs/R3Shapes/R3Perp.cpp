/* Source file for GAPS perpendicular utility */



/* Include files */

#include "R3Shapes/R3Shapes.h"



RNBoolean R3Perpendicular(const R3Vector& vector1, const R3Vector& vector2)
{
    // Normalized vectors ???
    // Return whether vector1 and vector2 are perpendicular
    return RNIsZero(vector1.Dot(vector2));
}



RNBoolean R3Perpendicular(const R3Vector& vector, const R3Line& line)
{
    // Return whether vector and line are perpendicular
    return R3Perpendicular(vector, line.Vector());
}



RNBoolean R3Perpendicular(const R3Vector& vector, const R3Ray& ray)
{
    // Return whether vector and ray are perpendicular
    return R3Perpendicular(vector, ray.Vector());
}



RNBoolean R3Perpendicular(const R3Vector& vector, const R3Span& span)
{
    // Return whether vector and span are perpendicular
    return R3Perpendicular(vector, span.Vector());
}



RNBoolean R3Perpendicular(const R3Vector& vector, const R3Plane& plane)
{
    // Return whether vector and plane are perpendicular
    return R3Parallel(vector, plane.Normal());
}



RNBoolean R3Perpendicular(const R3Line& line1, const R3Line& line2)
{
    // Return whether line1 and line2 are perpendicular
    return R3Perpendicular(line1.Vector(), line2.Vector());
}



RNBoolean R3Perpendicular(const R3Line& line, const R3Ray& ray)
{
    // Return whether line and ray are perpendicular
    return R3Perpendicular(line.Vector(), ray.Vector());
}



RNBoolean R3Perpendicular(const R3Line& line, const R3Span& span)
{
    // Return whether line and span are perpendicular
    return R3Perpendicular(line.Vector(), span.Vector());
}



RNBoolean R3Perpendicular(const R3Line& line, const R3Plane& plane)
{
    // Return whether line and plane are perpendicular
    return R3Parallel(line.Vector(), plane.Normal());
}



RNBoolean R3Perpendicular(const R3Ray& ray1, const R3Ray& ray2)
{
    // Return whether ray1 and ray2 are perpendicular
    return R3Perpendicular(ray1.Vector(), ray2.Vector());
}



RNBoolean R3Perpendicular(const R3Ray& ray, const R3Span& span)
{
    // Return whether ray and span are perpendicular
    return R3Perpendicular(ray.Vector(), span.Vector());
}



RNBoolean R3Perpendicular(const R3Ray& ray, const R3Plane& plane)
{
    // Return whether ray and plane are perpendicular
    return R3Parallel(ray.Vector(), plane.Normal());
}



RNBoolean R3Perpendicular(const R3Span& span1, const R3Span& span2)
{
    // Return whether span1 and span2 are perpendicular
    return R3Perpendicular(span1.Vector(), span2.Vector());
}



RNBoolean R3Perpendicular(const R3Span& span, const R3Plane& plane)
{
    // Return whether span and plane are perpendicular
    return R3Parallel(span.Vector(), plane.Normal());
}



RNBoolean R3Perpendicular(const R3Plane& plane1, const R3Plane& plane2)
{
    // Return whether plane1 and plane2 are perpendicular
    return R3Perpendicular(plane1.Normal(), plane2.Normal());
}







