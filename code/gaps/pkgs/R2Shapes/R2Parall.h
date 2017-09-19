/* Include file for GAPS parallel utility */



/* Function declarations */

RNBoolean R2Parallel(const R2Vector& vector1, const R2Vector& vector2);
RNBoolean R2Parallel(const R2Vector& vector, const R2Line& line);
RNBoolean R2Parallel(const R2Vector& vector, const R2Ray& ray);
RNBoolean R2Parallel(const R2Vector& vector, const R2Span& span);

RNBoolean R2Parallel(const R2Line& line, const R2Vector& vector);
RNBoolean R2Parallel(const R2Line& line1, const R2Line& line2);
RNBoolean R2Parallel(const R2Line& line, const R2Ray& ray);
RNBoolean R2Parallel(const R2Line& line, const R2Span& span);

RNBoolean R2Parallel(const R2Ray& ray, const R2Vector& vector);
RNBoolean R2Parallel(const R2Ray& ray, const R2Line& line);
RNBoolean R2Parallel(const R2Ray& ray1, const R2Ray& ray2);
RNBoolean R2Parallel(const R2Ray& ray, const R2Span& span);

RNBoolean R2Parallel(const R2Span& span, const R2Vector& vector);
RNBoolean R2Parallel(const R2Span& span, const R2Line& line);
RNBoolean R2Parallel(const R2Span& span, const R2Ray& ray);
RNBoolean R2Parallel(const R2Span& span1, const R2Span& span2);



/* Inline functions */

inline RNBoolean R2Parallel(const R2Line& line, const R2Vector& vector)
{
    // Parallel is commutative
    return R2Parallel(vector, line);
}



inline RNBoolean R2Parallel(const R2Ray& ray, const R2Vector& vector)
{
    // Parallel is commutative
    return R2Parallel(vector, ray);
}



inline RNBoolean R2Parallel(const R2Ray& ray, const R2Line& line)
{
    // Parallel is commutative
    return R2Parallel(line, ray);
}



inline RNBoolean R2Parallel(const R2Span& span, const R2Vector& vector)
{
    // Parallel is commutative
    return R2Parallel(vector, span);
}



inline RNBoolean R2Parallel(const R2Span& span, const R2Line& line)
{
    // Parallel is commutative
    return R2Parallel(line, span);
}



inline RNBoolean R2Parallel(const R2Span& span, const R2Ray& ray)
{
    // Parallel is commutative
    return R2Parallel(ray, span);
}



