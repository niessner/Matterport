/* Include file for GAPS perpendicular utility */



/* Function declarations */

RNBoolean R3Perpendicular(const R3Vector& vector1, const R3Vector& vector2);
RNBoolean R3Perpendicular(const R3Vector& vector, const R3Line& line);
RNBoolean R3Perpendicular(const R3Vector& vector, const R3Ray& ray);
RNBoolean R3Perpendicular(const R3Vector& vector, const R3Span& span);
RNBoolean R3Perpendicular(const R3Vector& vector, const R3Plane& plane);

RNBoolean R3Perpendicular(const R3Line& line, const R3Vector& vector);
RNBoolean R3Perpendicular(const R3Line& line1, const R3Line& line2);
RNBoolean R3Perpendicular(const R3Line& line, const R3Ray& ray);
RNBoolean R3Perpendicular(const R3Line& line, const R3Span& span);
RNBoolean R3Perpendicular(const R3Line& line, const R3Plane& plane);

RNBoolean R3Perpendicular(const R3Ray& ray, const R3Vector& vector);
RNBoolean R3Perpendicular(const R3Ray& ray, const R3Line& line);
RNBoolean R3Perpendicular(const R3Ray& ray1, const R3Ray& ray2);
RNBoolean R3Perpendicular(const R3Ray& ray, const R3Span& span);
RNBoolean R3Perpendicular(const R3Ray& ray, const R3Plane& plane);

RNBoolean R3Perpendicular(const R3Span& span, const R3Vector& vector);
RNBoolean R3Perpendicular(const R3Span& span, const R3Line& line);
RNBoolean R3Perpendicular(const R3Span& span, const R3Ray& ray);
RNBoolean R3Perpendicular(const R3Span& span1, const R3Span& span2);
RNBoolean R3Perpendicular(const R3Span& span, const R3Plane& plane);

RNBoolean R3Perpendicular(const R3Plane& plane, const R3Vector& vector);
RNBoolean R3Perpendicular(const R3Plane& plane, const R3Line& line);
RNBoolean R3Perpendicular(const R3Plane& plane, const R3Ray& ray);
RNBoolean R3Perpendicular(const R3Plane& plane, const R3Span& span);
RNBoolean R3Perpendicular(const R3Plane& plane1, const R3Plane& plane2);



/* Inline functions */

inline RNBoolean R3Perpendicular(const R3Line& line, const R3Vector& vector)
{
    // Perpendicular is commutative
    return R3Perpendicular(vector, line);
}



inline RNBoolean R3Perpendicular(const R3Ray& ray, const R3Vector& vector)
{
    // Perpendicular is commutative
    return R3Perpendicular(vector, ray);
}



inline RNBoolean R3Perpendicular(const R3Ray& ray, const R3Line& line)
{
    // Perpendicular is commutative
    return R3Perpendicular(line, ray);
}



inline RNBoolean R3Perpendicular(const R3Span& span, const R3Vector& vector)
{
    // Perpendicular is commutative
    return R3Perpendicular(vector, span);
}



inline RNBoolean R3Perpendicular(const R3Span& span, const R3Line& line)
{
    // Perpendicular is commutative
    return R3Perpendicular(line, span);
}



inline RNBoolean R3Perpendicular(const R3Span& span, const R3Ray& ray)
{
    // Perpendicular is commutative
    return R3Perpendicular(ray, span);
}



inline RNBoolean R3Perpendicular(const R3Plane& plane, const R3Vector& vector)
{
    // Perpendicular is commutative
    return R3Perpendicular(vector, plane);
}



inline RNBoolean R3Perpendicular(const R3Plane& plane, const R3Line& line)
{
    // Perpendicular is commutative
    return R3Perpendicular(line, plane);
}



inline RNBoolean R3Perpendicular(const R3Plane& plane, const R3Ray& ray)
{
    // Perpendicular is commutative
    return R3Perpendicular(ray, plane);
}



inline RNBoolean R3Perpendicular(const R3Plane& plane, const R3Span& span)
{
    // Perpendicular is commutative
    return R3Perpendicular(span, plane);
}








