/* Include file for GAPS parallel utility */



/* Function declarations */

RNBoolean R3Parallel(const R3Vector& vector1, const R3Vector& vector2);
RNBoolean R3Parallel(const R3Vector& vector, const R3Line& line);
RNBoolean R3Parallel(const R3Vector& vector, const R3Ray& ray);
RNBoolean R3Parallel(const R3Vector& vector, const R3Span& span);
RNBoolean R3Parallel(const R3Vector& vector, const R3Plane& plane);

RNBoolean R3Parallel(const R3Line& line, const R3Vector& vector);
RNBoolean R3Parallel(const R3Line& line1, const R3Line& line2);
RNBoolean R3Parallel(const R3Line& line, const R3Ray& ray);
RNBoolean R3Parallel(const R3Line& line, const R3Span& span);
RNBoolean R3Parallel(const R3Line& line, const R3Plane& plane);

RNBoolean R3Parallel(const R3Ray& ray, const R3Vector& vector);
RNBoolean R3Parallel(const R3Ray& ray, const R3Line& line);
RNBoolean R3Parallel(const R3Ray& ray1, const R3Ray& ray2);
RNBoolean R3Parallel(const R3Ray& ray, const R3Span& span);
RNBoolean R3Parallel(const R3Ray& ray, const R3Plane& plane);

RNBoolean R3Parallel(const R3Span& span, const R3Vector& vector);
RNBoolean R3Parallel(const R3Span& span, const R3Line& line);
RNBoolean R3Parallel(const R3Span& span, const R3Ray& ray);
RNBoolean R3Parallel(const R3Span& span1, const R3Span& span2);
RNBoolean R3Parallel(const R3Span& span, const R3Plane& plane);

RNBoolean R3Parallel(const R3Plane& plane, const R3Vector& vector);
RNBoolean R3Parallel(const R3Plane& plane, const R3Line& line);
RNBoolean R3Parallel(const R3Plane& plane, const R3Ray& ray);
RNBoolean R3Parallel(const R3Plane& plane, const R3Span& span);
RNBoolean R3Parallel(const R3Plane& plane1, const R3Plane& plane2);



/* Inline functions */

inline RNBoolean R3Parallel(const R3Line& line, const R3Vector& vector)
{
    // Parallel is commutative
    return R3Parallel(vector, line);
}



inline RNBoolean R3Parallel(const R3Ray& ray, const R3Vector& vector)
{
    // Parallel is commutative
    return R3Parallel(vector, ray);
}



inline RNBoolean R3Parallel(const R3Ray& ray, const R3Line& line)
{
    // Parallel is commutative
    return R3Parallel(line, ray);
}



inline RNBoolean R3Parallel(const R3Span& span, const R3Vector& vector)
{
    // Parallel is commutative
    return R3Parallel(vector, span);
}



inline RNBoolean R3Parallel(const R3Span& span, const R3Line& line)
{
    // Parallel is commutative
    return R3Parallel(line, span);
}



inline RNBoolean R3Parallel(const R3Span& span, const R3Ray& ray)
{
    // Parallel is commutative
    return R3Parallel(ray, span);
}



inline RNBoolean R3Parallel(const R3Plane& plane, const R3Vector& vector)
{
    // Parallel is commutative
    return R3Parallel(vector, plane);
}



inline RNBoolean R3Parallel(const R3Plane& plane, const R3Line& line)
{
    // Parallel is commutative
    return R3Parallel(line, plane);
}



inline RNBoolean R3Parallel(const R3Plane& plane, const R3Ray& ray)
{
    // Parallel is commutative
    return R3Parallel(ray, plane);
}



inline RNBoolean R3Parallel(const R3Plane& plane, const R3Span& span)
{
    // Parallel is commutative
    return R3Parallel(span, plane);
}





