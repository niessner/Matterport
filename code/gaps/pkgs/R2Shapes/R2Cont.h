/* Include file for GAPS inside/containment utility */



/* Useful macros */

#define R2Inside(__primitive1, __primitive2) \
    R2Contains(__primitive2, __primitive1)



/* Function declarations */

RNBoolean R2Contains(const R2Vector& vector1, const R2Vector& vector2);

RNBoolean R2Contains(const R2Point& point1, const R2Point& point2);
RNBoolean R2Contains(const R2Point& point, const R2Line& line);
RNBoolean R2Contains(const R2Point& point, const R2Ray& ray);
RNBoolean R2Contains(const R2Point& point, const R2Span& span);
RNBoolean R2Contains(const R2Point& point, const R2Halfspace& halfspace);
RNBoolean R2Contains(const R2Point& point, const R2Box& box);
RNBoolean R2Contains(const R2Point& point, const R2Circle& circle);
RNBoolean R2Contains(const R2Point& point, const R2Shape& shape);

RNBoolean R2Contains(const R2Line& line, const R2Point& point);
RNBoolean R2Contains(const R2Line& line1, const R2Line& line2);
RNBoolean R2Contains(const R2Line& line, const R2Ray& ray);
RNBoolean R2Contains(const R2Line& line, const R2Span& span);
RNBoolean R2Contains(const R2Line& line, const R2Halfspace& halfspace);
RNBoolean R2Contains(const R2Line& line, const R2Box& box);
RNBoolean R2Contains(const R2Line& line, const R2Circle& circle);
RNBoolean R2Contains(const R2Line& line, const R2Shape& shape);

RNBoolean R2Contains(const R2Ray& ray, const R2Point& point);
RNBoolean R2Contains(const R2Ray& ray, const R2Line& line);
RNBoolean R2Contains(const R2Ray& ray1, const R2Ray& ray2);
RNBoolean R2Contains(const R2Ray& ray, const R2Span& span);
RNBoolean R2Contains(const R2Ray& ray, const R2Halfspace& halfspace);
RNBoolean R2Contains(const R2Ray& ray, const R2Box& box);
RNBoolean R2Contains(const R2Ray& ray, const R2Circle& circle);
RNBoolean R2Contains(const R2Ray& ray, const R2Shape& shape);

RNBoolean R2Contains(const R2Span& span, const R2Point& point);
RNBoolean R2Contains(const R2Span& span, const R2Line& line);
RNBoolean R2Contains(const R2Span& span, const R2Ray& ray);
RNBoolean R2Contains(const R2Span& span1, const R2Span& span2);
RNBoolean R2Contains(const R2Span& span, const R2Halfspace& halfspace);
RNBoolean R2Contains(const R2Span& span, const R2Box& box);
RNBoolean R2Contains(const R2Span& span, const R2Circle& circle);
RNBoolean R2Contains(const R2Span& span, const R2Shape& shape);

RNBoolean R2Contains(const R2Halfspace& halfspace, const R2Point& point);
RNBoolean R2Contains(const R2Halfspace& halfspace, const R2Line& line);
RNBoolean R2Contains(const R2Halfspace& halfspace, const R2Ray& ray);
RNBoolean R2Contains(const R2Halfspace& halfspace, const R2Span& span);
RNBoolean R2Contains(const R2Halfspace& halfspace1, const R2Halfspace& halfspace2);
RNBoolean R2Contains(const R2Halfspace& halfspace, const R2Box& box);
RNBoolean R2Contains(const R2Halfspace& halfspace, const R2Circle& circle);
RNBoolean R2Contains(const R2Halfspace& halfspace, const R2Shape& shape);

RNBoolean R2Contains(const R2Box& box, const R2Point& point);
RNBoolean R2Contains(const R2Box& box, const R2Line& line);
RNBoolean R2Contains(const R2Box& box, const R2Ray& ray);
RNBoolean R2Contains(const R2Box& box, const R2Span& span);
RNBoolean R2Contains(const R2Box& box, const R2Halfspace& halfspace);
RNBoolean R2Contains(const R2Box& box1, const R2Box& box2);
RNBoolean R2Contains(const R2Box& box, const R2Circle& circle);
RNBoolean R2Contains(const R2Box& box, const R2Shape& shape);

RNBoolean R2Contains(const R2Circle& circle, const R2Point& point);
RNBoolean R2Contains(const R2Circle& circle, const R2Line& line);
RNBoolean R2Contains(const R2Circle& circle, const R2Ray& ray);
RNBoolean R2Contains(const R2Circle& circle, const R2Span& span);
RNBoolean R2Contains(const R2Circle& circle, const R2Halfspace& halfspace);
RNBoolean R2Contains(const R2Circle& circle, const R2Box& box);
RNBoolean R2Contains(const R2Circle& circle1, const R2Circle& circle2);
RNBoolean R2Contains(const R2Circle& circle, const R2Shape& shape);

RNBoolean R2Contains(const R2Polygon& polygon, const R2Point& point);
RNBoolean R2Contains(const R2Polygon& polygon, const R2Line& line);
RNBoolean R2Contains(const R2Polygon& polygon, const R2Ray& ray);
RNBoolean R2Contains(const R2Polygon& polygon, const R2Span& span);
RNBoolean R2Contains(const R2Polygon& polygon, const R2Halfspace& halfspace);
RNBoolean R2Contains(const R2Polygon& polygon, const R2Box& box);
RNBoolean R2Contains(const R2Polygon& polygon, const R2Circle& circle);
RNBoolean R2Contains(const R2Polygon& polygon1, const R2Polygon& polygon2);
RNBoolean R2Contains(const R2Polygon& polygon, const R2Shape& shape);

RNBoolean R2Contains(const R2Shape& shape, const R2Point& point);
RNBoolean R2Contains(const R2Shape& shape, const R2Line& line);
RNBoolean R2Contains(const R2Shape& shape, const R2Ray& ray);
RNBoolean R2Contains(const R2Shape& shape, const R2Span& span);
RNBoolean R2Contains(const R2Shape& shape, const R2Halfspace& halfspace);
RNBoolean R2Contains(const R2Shape& shape, const R2Box& box);
RNBoolean R2Contains(const R2Shape& shape, const R2Circle& circle);
RNBoolean R2Contains(const R2Shape& shape1, const R2Shape& shape2);











