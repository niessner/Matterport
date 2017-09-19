/* Include file for R2 miscelleaneous relationship utility */



/* Miscellaneous relationship function declarations */

RNAngle R2InteriorAngle(const R2Vector& vector1, const R2Vector& vector2);
RNAngle R2InteriorAngle(const R2Vector& vector, const R2Line& line);
RNAngle R2InteriorAngle(const R2Vector& vector, const R2Ray& ray);
RNAngle R2InteriorAngle(const R2Vector& vector, const R2Span& span);

RNArea R2IntersectionArea(const R2Box& box1, const R2Box& box2);
RNArea R2IntersectionArea(const R2Circle& circle1, const R2Circle& circle2);

RNBoolean R2Abuts(const R2Box& box1, const R2Box& box2);

int R2Splits(const R2Line& line, const R2Span& span);
int R2Splits(const R2Line& line, const R2Span& span, R2Span *below_result, R2Span *above_result = NULL);
