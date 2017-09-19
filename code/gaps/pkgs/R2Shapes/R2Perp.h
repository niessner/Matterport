/* Include file for GAPS perpendicular utility */



/* Function declarations */

RNBoolean R2Perpendicular(const R2Vector& vector1, const R2Vector& vector2);
RNBoolean R2Perpendicular(const R2Vector& vector, const R2Line& line);
RNBoolean R2Perpendicular(const R2Vector& vector, const R2Ray& ray);
RNBoolean R2Perpendicular(const R2Vector& vector, const R2Span& span);

RNBoolean R2Perpendicular(const R2Line& line, const R2Vector& vector);
RNBoolean R2Perpendicular(const R2Line& line1, const R2Line& line2);
RNBoolean R2Perpendicular(const R2Line& line, const R2Ray& ray);
RNBoolean R2Perpendicular(const R2Line& line, const R2Span& span);

RNBoolean R2Perpendicular(const R2Ray& ray, const R2Vector& vector);
RNBoolean R2Perpendicular(const R2Ray& ray, const R2Line& line);
RNBoolean R2Perpendicular(const R2Ray& ray1, const R2Ray& ray2);
RNBoolean R2Perpendicular(const R2Ray& ray, const R2Span& span);

RNBoolean R2Perpendicular(const R2Span& span, const R2Vector& vector);
RNBoolean R2Perpendicular(const R2Span& span, const R2Line& line);
RNBoolean R2Perpendicular(const R2Span& span, const R2Ray& ray);
RNBoolean R2Perpendicular(const R2Span& span1, const R2Span& span2);



/* Inline functions */

inline RNBoolean R2Perpendicular(const R2Line& line, const R2Vector& vector)
{
    // Perpendicular is commutative
    return R2Perpendicular(vector, line);
}



inline RNBoolean R2Perpendicular(const R2Ray& ray, const R2Vector& vector)
{
    // Perpendicular is commutative
    return R2Perpendicular(vector, ray);
}



inline RNBoolean R2Perpendicular(const R2Ray& ray, const R2Line& line)
{
    // Perpendicular is commutative
    return R2Perpendicular(line, ray);
}



inline RNBoolean R2Perpendicular(const R2Span& span, const R2Vector& vector)
{
    // Perpendicular is commutative
    return R2Perpendicular(vector, span);
}



inline RNBoolean R2Perpendicular(const R2Span& span, const R2Line& line)
{
    // Perpendicular is commutative
    return R2Perpendicular(line, span);
}



inline RNBoolean R2Perpendicular(const R2Span& span, const R2Ray& ray)
{
    // Perpendicular is commutative
    return R2Perpendicular(ray, span);
}



