/* Include file for GAPS distance utility */



/* Function declarations */

RNLength R2Distance(const R2Point& point1, const R2Point& point2);
RNLength R2Distance(const R2Point& point, const R2Line& line);
RNLength R2Distance(const R2Point& point, const R2Ray& ray);
RNLength R2Distance(const R2Point& point, const R2Span& span);
RNLength R2Distance(const R2Point& point, const R2Halfspace& halfspace);
RNLength R2Distance(const R2Point& point, const R2Arc& arc);
RNLength R2Distance(const R2Point& point, const R2Arc& arc);
RNLength R2Distance(const R2Point& point, const R2Polyline& polyline);
RNLength R2Distance(const R2Point& point, const R2Polygon& polygon);
RNLength R2Distance(const R2Point& point, const R2Box& box);
RNLength R2Distance(const R2Point& point, const R2Circle& circle);
RNLength R2SquaredDistance(const R2Point& point1, const R2Point& point2);

RNLength R2Distance(const R2Line& line, const R2Point& point);
RNLength R2Distance(const R2Line& line1, const R2Line& line2);
RNLength R2Distance(const R2Line& line, const R2Ray& ray);
RNLength R2Distance(const R2Line& line, const R2Span& span);
RNLength R2Distance(const R2Line& line, const R2Halfspace& halfspace);
RNLength R2Distance(const R2Line& line, const R2Box& box);
RNLength R2Distance(const R2Line& line, const R2Circle& circle);

RNLength R2SignedDistance(const R2Line& line, const R2Point& point);
RNLength R2SignedDistance(const R2Line& line1, const R2Line& line2);
RNLength R2SignedDistance(const R2Line& line, const R2Ray& ray);
RNLength R2SignedDistance(const R2Line& line, const R2Span& span);
RNLength R2SignedDistance(const R2Line& line, const R2Halfspace& halfspace);
RNLength R2SignedDistance(const R2Line& line, const R2Box& box);
RNLength R2SignedDistance(const R2Line& line, const R2Circle& circle);

RNLength R2Distance(const R2Ray& ray, const R2Point& point);
RNLength R2Distance(const R2Ray& ray, const R2Line& line);
RNLength R2Distance(const R2Ray& ray1, const R2Ray& ray2);
RNLength R2Distance(const R2Ray& ray, const R2Span& span);
RNLength R2Distance(const R2Ray& ray, const R2Halfspace& halfspace);
RNLength R2Distance(const R2Ray& ray, const R2Box& box);
RNLength R2Distance(const R2Ray& ray, const R2Circle& circle);

RNLength R2Distance(const R2Span& span, const R2Point& point);
RNLength R2Distance(const R2Span& span, const R2Line& line);
RNLength R2Distance(const R2Span& span, const R2Ray& ray);
RNLength R2Distance(const R2Span& span1, const R2Span& span2);
RNLength R2Distance(const R2Span& span, const R2Halfspace& halfspace);
RNLength R2Distance(const R2Span& span, const R2Box& box);
RNLength R2Distance(const R2Span& span, const R2Circle& circle);

RNLength R2Distance(const R2Halfspace& halfspace, const R2Point& point);
RNLength R2Distance(const R2Halfspace& halfspace, const R2Line& line);
RNLength R2Distance(const R2Halfspace& halfspace, const R2Ray& ray);
RNLength R2Distance(const R2Halfspace& halfspace, const R2Span& span);
RNLength R2Distance(const R2Halfspace& halfspace1, const R2Halfspace& halfspace2);
RNLength R2Distance(const R2Halfspace& halfspace, const R2Box& box);
RNLength R2Distance(const R2Halfspace& halfspace, const R2Circle& circle);

RNLength R2Distance(const R2Arc& arc, const R2Point& point);

RNLength R2Distance(const R2Polyline& polyline, const R2Point& point);

RNLength R2Distance(const R2Polygon& polygon, const R2Point& point);

RNLength R2Distance(const R2Box& box, const R2Point& point);
RNLength R2Distance(const R2Box& box, const R2Line& line);
RNLength R2Distance(const R2Box& box, const R2Ray& ray);
RNLength R2Distance(const R2Box& box, const R2Span& span);
RNLength R2Distance(const R2Box& box, const R2Halfspace& halfspace);
RNLength R2Distance(const R2Box& box1, const R2Box& box2);
RNLength R2Distance(const R2Box& box, const R2Circle& circle);

RNLength R2Distance(const R2Circle& circle, const R2Point& point);
RNLength R2Distance(const R2Circle& circle, const R2Line& line);
RNLength R2Distance(const R2Circle& circle, const R2Ray& ray);
RNLength R2Distance(const R2Circle& circle, const R2Span& span);
RNLength R2Distance(const R2Circle& circle, const R2Halfspace& halfspace);
RNLength R2Distance(const R2Circle& circle, const R2Box& box);
RNLength R2Distance(const R2Circle& circle1, const R2Circle& circle2);



/* Inline functions */

inline RNLength R2Distance(const R2Line& line, const R2Point& point)
{
    // Distance is commutative
    return R2Distance(point, line);
}



inline RNLength R2Distance(const R2Ray& ray, const R2Point& point)
{
    // Distance is commutative
    return R2Distance(point, ray);
}



inline RNLength R2Distance(const R2Ray& ray, const R2Line& line)
{
    // Distance is commutative
    return R2Distance(line, ray);
}



inline RNLength R2Distance(const R2Span& span, const R2Point& point)
{
    // Distance is commutative
    return R2Distance(point, span);
}



inline RNLength R2Distance(const R2Span& span, const R2Line& line)
{
    // Distance is commutative
    return R2Distance(line, span);
}



inline RNLength R2Distance(const R2Span& span, const R2Ray& ray)
{
    // Distance is commutative
    return R2Distance(ray, span);
}



inline RNLength R2Distance(const R2Halfspace& halfspace, const R2Point& point)
{
    // Distance is commutative
    return R2Distance(point, halfspace);
}



inline RNLength R2Distance(const R2Halfspace& halfspace, const R2Line& line)
{
    // Distance is commutative
    return R2Distance(line, halfspace);
}



inline RNLength R2Distance(const R2Halfspace& halfspace, const R2Ray& ray)
{
    // Distance is commutative
    return R2Distance(ray, halfspace);
}



inline RNLength R2Distance(const R2Halfspace& halfspace, const R2Span& span)
{
    // Distance is commutative
    return R2Distance(span, halfspace);
}



inline RNLength R2Distance(const R2Arc& arc, const R2Point& point)
{
    // Distance is commutative
    return R2Distance(point, arc);
}



inline RNLength R2Distance(const R2Polyline& polyline, const R2Point& point)
{
    // Distance is commutative
    return R2Distance(point, polyline);
}



inline RNLength R2Distance(const R2Polygon& polygon, const R2Point& point)
{
    // Distance is commutative
    return R2Distance(point, polygon);
}



inline RNLength R2Distance(const R2Box& box, const R2Point& point)
{
    // Distance is commutative
    return R2Distance(point, box);
}



inline RNLength R2Distance(const R2Box& box, const R2Line& line)
{
    // Distance is commutative
    return R2Distance(line, box);
}



inline RNLength R2Distance(const R2Box& box, const R2Ray& ray)
{
    // Distance is commutative
    return R2Distance(ray, box);
}



inline RNLength R2Distance(const R2Box& box, const R2Span& span)
{
    // Distance is commutative
    return R2Distance(span, box);
}



inline RNLength R2Distance(const R2Box& box, const R2Halfspace& halfspace)
{
    // Distance is commutative
    return R2Distance(halfspace, box);
}



inline RNLength R2Distance(const R2Circle& circle, const R2Point& point)
{
    // Distance is commutative
    return R2Distance(point, circle);
}



inline RNLength R2Distance(const R2Circle& circle, const R2Line& line)
{
    // Distance is commutative
    return R2Distance(line, circle);
}



inline RNLength R2Distance(const R2Circle& circle, const R2Ray& ray)
{
    // Distance is commutative
    return R2Distance(ray, circle);
}



inline RNLength R2Distance(const R2Circle& circle, const R2Span& span)
{
    // Distance is commutative
    return R2Distance(span, circle);
}



inline RNLength R2Distance(const R2Circle& circle, const R2Halfspace& halfspace)
{
    // Distance is commutative
    return R2Distance(halfspace, circle);
}



inline RNLength R2Distance(const R2Circle& circle, const R2Box& box)
{
    // Distance is commutative
    return R2Distance(box, circle);
}




