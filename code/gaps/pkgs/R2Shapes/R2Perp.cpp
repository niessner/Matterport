/* Source file for GAPS perpendicular utility */



/* Include files */

#include "R2Shapes/R2Shapes.h"



RNBoolean R2Perpendicular(const R2Vector& vector1, const R2Vector& vector2)
{
    // Normalized vectors ???
    // Return whether vector1 and vector2 are perpendicular
    return RNIsZero(vector1.Dot(vector2));
}



RNBoolean R2Perpendicular(const R2Vector& vector, const R2Line& line)
{
    // Return whether vector and line are perpendicular
    return R2Perpendicular(vector, line.Vector());
}



RNBoolean R2Perpendicular(const R2Vector& vector, const R2Ray& ray)
{
    // Return whether vector and ray are perpendicular
    return R2Perpendicular(vector, ray.Vector());
}



RNBoolean R2Perpendicular(const R2Vector& vector, const R2Span& span)
{
    // Return whether vector and span are perpendicular
    return R2Perpendicular(vector, span.Vector());
}



RNBoolean R2Perpendicular(const R2Line& line1, const R2Line& line2)
{
    // Return whether line1 and line2 are perpendicular
    return R2Perpendicular(line1.Vector(), line2.Vector());
}



RNBoolean R2Perpendicular(const R2Line& line, const R2Ray& ray)
{
    // Return whether line and ray are perpendicular
    return R2Perpendicular(line.Vector(), ray.Vector());
}



RNBoolean R2Perpendicular(const R2Line& line, const R2Span& span)
{
    // Return whether line and span are perpendicular
    return R2Perpendicular(line.Vector(), span.Vector());
}



RNBoolean R2Perpendicular(const R2Ray& ray1, const R2Ray& ray2)
{
    // Return whether ray1 and ray2 are perpendicular
    return R2Perpendicular(ray1.Vector(), ray2.Vector());
}



RNBoolean R2Perpendicular(const R2Ray& ray, const R2Span& span)
{
    // Return whether ray and span are perpendicular
    return R2Perpendicular(ray.Vector(), span.Vector());
}



RNBoolean R2Perpendicular(const R2Span& span1, const R2Span& span2)
{
    // Return whether span1 and span2 are perpendicular
    return R2Perpendicular(span1.Vector(), span2.Vector());
}




