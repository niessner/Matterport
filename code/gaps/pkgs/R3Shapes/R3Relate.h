/* Include file for R3 miscelleaneous relationship utility */



/* Miscellaneous relationship function declarations */

R3Point R3Interpolate(const R3Point& point1, const R3Point& point2, RNScalar t1, RNScalar t2, RNScalar t);

RNAngle R3InteriorAngle(const R3Vector& vector1, const R3Vector& vector2);
RNAngle R3InteriorAngle(const R3Vector& vector, const R3Line& line);
RNAngle R3InteriorAngle(const R3Vector& vector, const R3Ray& ray);
RNAngle R3InteriorAngle(const R3Vector& vector, const R3Span& span);
RNAngle R3InteriorAngle(const R3Vector& vector, const R3Plane& plane);

RNBoolean R3Abuts(const R3Box& box1, const R3Box& box2);

int R3Splits(const R3Plane& plane, const R3Span& span);
int R3Splits(const R3Plane& plane, const R3Span& span, R3Span *below_result, R3Span *above_result = NULL);
int R3Splits(const R3Plane& plane, const R3Triangle& triangle);
