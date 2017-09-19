/* Source file for GAPS parallel utility */



/* Include files */

#include "R2Shapes/R2Shapes.h"



RNBoolean R2Parallel(const R2Vector& vector1, const R2Vector& vector2)
{
    // Return whether two vectors are coincident (or anti-coincident)
    if ((RNIsEqual(vector1.X(), vector2.X())) &&
        (RNIsEqual(vector1.Y(), vector2.Y()))) return TRUE;
    if ((RNIsEqual(vector1.X(), -vector2.X())) &&
        (RNIsEqual(vector1.Y(), -vector2.Y()))) return TRUE;
    return FALSE;
}



RNBoolean R2Parallel(const R2Vector& vector, const R2Line& line)
{
    // Return whether vector and line are parallel
    return R2Parallel(vector, line.Vector());
}



RNBoolean R2Parallel(const R2Vector& vector, const R2Ray& ray)
{
    // Return whether vector and ray are parallel
    return R2Parallel(vector, ray.Vector());
}



RNBoolean R2Parallel(const R2Vector& vector, const R2Span& span)
{
    // Return whether vector and span are parallel
    return R2Parallel(vector, span.Vector());
}



RNBoolean R2Parallel(const R2Line& line1, const R2Line& line2)
{
    // Return whether line1 and line2 are parallel
    return R2Parallel(line1.Vector(), line2.Vector());
}



RNBoolean R2Parallel(const R2Line& line, const R2Ray& ray)
{
    // Return whether line and ray are parallel
    return R2Parallel(line.Vector(), ray.Vector());
}



RNBoolean R2Parallel(const R2Line& line, const R2Span& span)
{
    // Return whether line and span are parallel
    return R2Parallel(line.Vector(), span.Vector());
}



RNBoolean R2Parallel(const R2Ray& ray1, const R2Ray& ray2)
{
    // Return whether ray1 and ray2 are parallel
    return R2Parallel(ray1.Vector(), ray2.Vector());
}



RNBoolean R2Parallel(const R2Ray& ray, const R2Span& span)
{
    // Return whether ray and span are parallel
    return R2Parallel(ray.Vector(), span.Vector());
}



RNBoolean R2Parallel(const R2Span& span1, const R2Span& span2)
{
    // Return whether span1 and span2 are parallel
    return R2Parallel(span1.Vector(), span2.Vector());
}




