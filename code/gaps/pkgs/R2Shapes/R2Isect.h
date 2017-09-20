/* Include file for GAPS intersection utility */



/* Function declarations */

RNClassID R2Intersects(const R2Point& point1, const R2Point& point2);
RNClassID R2Intersects(const R2Point& point, const R2Line& line);
RNClassID R2Intersects(const R2Point& point, const R2Ray& ray, 
    RNScalar *hit_t = NULL); 
RNClassID R2Intersects(const R2Point& point, const R2Span& span, 
    RNScalar *hit_t = NULL); 
RNClassID R2Intersects(const R2Point& point, const R2Halfspace& halfspace); 
RNClassID R2Intersects(const R2Point& point, const R2Box& box); 
RNClassID R2Intersects(const R2Point& point, const R2Circle& circle); 

RNClassID R2Intersects(const R2Line& line, const R2Point& point);
RNClassID R2Intersects(const R2Line& line1, const R2Line& line2, 
    R2Point *hit_point = NULL);
RNClassID R2Intersects(const R2Line& line, const R2Ray& ray, 
    R2Point *hit_point = NULL, RNScalar *hit_t = NULL);
RNClassID R2Intersects(const R2Line& line, const R2Span& span, 
    R2Point *hit_point = NULL, RNScalar *hit_t = NULL);
RNClassID R2Intersects(const R2Line& line, const R2Halfspace& halfspace, 
    R2Ray *result = NULL);
RNClassID R2Intersects(const R2Line& line, const R2Box& box, 
    R2Point *hit_point = NULL);
RNClassID R2Intersects(const R2Line& line, const R2Circle& circle, 
    R2Point *hit_point = NULL);

RNClassID R2Intersects(const R2Ray& ray, const R2Point& point, 
    RNScalar *hit_t = NULL); 
RNClassID R2Intersects(const R2Ray& ray, const R2Line& line, 
    R2Point *hit_point = NULL, RNScalar *hit_t = NULL);
RNClassID R2Intersects(const R2Ray& ray1, const R2Ray& ray2, 
    R2Point *hit_point = NULL, RNScalar *hit_t1 = NULL, RNScalar *hit_t2 = NULL);
RNClassID R2Intersects(const R2Ray& ray, const R2Span& span, 
    R2Point *hit_point = NULL, RNScalar *hit_tray = NULL, RNScalar *hit_tspan = NULL);
RNClassID R2Intersects(const R2Ray& ray, const R2Halfspace& halfspace, 
    R2Point *hit_point = NULL, RNScalar *hit_t = NULL);
RNClassID R2Intersects(const R2Ray& ray, const R2Box& box, 
    R2Point *hit_point1 = NULL, R2Vector *hit_normal1 = NULL, RNScalar *hit_t1 = NULL);
RNClassID R2Intersects(const R2Ray& ray, const R2Circle& circle, 
    R2Point *hit_point1 = NULL, R2Vector *hit_normal1 = NULL, RNScalar *hit_t1 = NULL);

RNClassID R2Intersects(const R2Span& span, const R2Point& point, 
    RNScalar *hit_t = NULL);
RNClassID R2Intersects(const R2Span& span, const R2Line& line, 
    R2Point *hit_point = NULL, RNScalar *hit_t = NULL);
RNClassID R2Intersects(const R2Span& span, const R2Ray& ray, 
    R2Point *hit_point = NULL, RNScalar *hit_tspan = NULL, RNScalar *hit_tray = NULL);
RNClassID R2Intersects(const R2Span& span1, const R2Span& span2, 
    R2Point *hit_point = NULL, RNScalar *hit_t1 = NULL, RNScalar *hit_t2 = NULL);
RNClassID R2Intersects(const R2Span& span, const R2Halfspace& halfspace, 
    R2Point *hit_point = NULL, RNScalar *hit_t = NULL);
RNClassID R2Intersects(const R2Span& span, const R2Box& box, 
    R2Point *hit_point1 = NULL, RNScalar *hit_t1 = NULL);
RNClassID R2Intersects(const R2Span& span, const R2Circle& circle, 
    R2Point *hit_point1 = NULL, RNScalar *hit_t1 = NULL);

RNClassID R2Intersects(const R2Halfspace& halfspace, const R2Point& point);
RNClassID R2Intersects(const R2Halfspace& halfspace, const R2Line& line, 
    R2Ray *result = NULL);
RNClassID R2Intersects(const R2Halfspace& halfspace, const R2Ray& ray, 
    R2Point *hit_point = NULL, R2Vector *hit_normal = NULL, RNScalar *hit_t = NULL);
RNClassID R2Intersects(const R2Halfspace& halfspace, const R2Span& span, 
    R2Point *hit_point = NULL, RNScalar *hit_t = NULL);
RNClassID R2Intersects(const R2Halfspace& halfspace1, const R2Halfspace& halfspace2);
RNClassID R2Intersects(const R2Halfspace& halfspace, const R2Box& box);
RNClassID R2Intersects(const R2Halfspace& halfspace, const R2Circle& circle);


RNClassID R2Intersects(const R2Box& box, const R2Point& point);
RNClassID R2Intersects(const R2Box& box, const R2Line& line, 
    R2Point *hit_point = NULL);
RNClassID R2Intersects(const R2Box& box, const R2Ray& ray, 
    R2Point *hit_point1 = NULL, R2Vector *hit_normal = NULL, RNScalar *hit_t1 = NULL);
RNClassID R2Intersects(const R2Box& box, const R2Span& span, 
    R2Point *hit_point1 = NULL, RNScalar *hit_t1 = NULL);
RNClassID R2Intersects(const R2Box& box, const R2Halfspace& halfspace);
RNClassID R2Intersects(const R2Box& box1, const R2Box& box2, 
    R2Box *result = NULL);
RNClassID R2Intersects(const R2Box& box, const R2Circle& circle);

RNClassID R2Intersects(const R2Circle& circle, const R2Point& point); 
RNClassID R2Intersects(const R2Circle& circle, const R2Line& line, 
    R2Point *hit_point = NULL);
RNClassID R2Intersects(const R2Circle& circle, const R2Ray& ray, 
    R2Point *hit_point1 = NULL, R2Vector *hit_normal = NULL, RNScalar *hit_t1 = NULL);
RNClassID R2Intersects(const R2Circle& circle, const R2Span& span, 
    R2Point *hit_point1 = NULL, RNScalar *hit_t1 = NULL);
RNClassID R2Intersects(const R2Circle& circle, const R2Halfspace& halfspace);
RNClassID R2Intersects(const R2Circle& circle, const R2Box& box);
RNClassID R2Intersects(const R2Circle& circle1, const R2Circle& circle2);



/* Inline functions */

inline RNClassID R2Intersects(const R2Line& line, const R2Point& point)
{
    // Intersection is commutative
    return R2Intersects(point, line);
}



inline RNClassID R2Intersects(const R2Ray& ray, const R2Point& point, RNScalar *hit_t)
{
    // Intersection is commutative
    return R2Intersects(point, ray, hit_t);
}



inline RNClassID R2Intersects(const R2Ray& ray, const R2Line& line, R2Point *hit_point, RNScalar *hit_t)
{
    // Intersection is commutative
    return R2Intersects(line, ray, hit_point, hit_t);
}



inline RNClassID R2Intersects(const R2Span& span, const R2Point& point, RNScalar *hit_t)
{
    // Intersection is commutative
    return R2Intersects(point, span, hit_t);
}



inline RNClassID R2Intersects(const R2Span& span, const R2Line& line, R2Point *hit_point, RNScalar *hit_t)
{
    // Intersection is commutative
    return R2Intersects(line, span, hit_point, hit_t);
}



inline RNClassID R2Intersects(const R2Span& span, const R2Ray& ray, R2Point *hit_point, RNScalar *hit_tspan, RNScalar *hit_tray)
{
    // Intersection is commutative
    return R2Intersects(ray, span, hit_point, hit_tray, hit_tspan);
}



inline RNClassID R2Intersects(const R2Halfspace& halfspace, const R2Point& point)
{
    // Intersection is commutative
    return R2Intersects(point, halfspace);
}



inline RNClassID R2Intersects(const R2Halfspace& halfspace, const R2Line& line, R2Ray *result)
{
    // Intersection is commutative
    return R2Intersects(line, halfspace, result);
}



inline RNClassID R2Intersects(const R2Halfspace& halfspace, const R2Ray& ray, R2Point *hit_point, RNScalar *hit_t)
{
    // Intersection is commutative
    return R2Intersects(ray, halfspace, hit_point, hit_t);
}



inline RNClassID R2Intersects(const R2Halfspace& halfspace, const R2Span& span, R2Point *hit_point, RNScalar *hit_t)
{
    // Intersection is commutative
    return R2Intersects(span, halfspace, hit_point, hit_t);
}



inline RNClassID R2Intersects(const R2Box& box, const R2Point& point)
{
    // Intersection is commutative
    return R2Intersects(point, box);
}



inline RNClassID R2Intersects(const R2Box& box, const R2Line& line, R2Point *hit_point)
{
    // Intersection is commutative
    return R2Intersects(line, box, hit_point);
}



inline RNClassID R2Intersects(const R2Box& box, const R2Ray& ray, 
    R2Point *hit_point1, R2Vector *hit_normal1, RNScalar *hit_t1)
{
    // Intersection is commutative
    return R2Intersects(ray, box, hit_point1, hit_normal1, hit_t1);
}



inline RNClassID R2Intersects(const R2Box& box, const R2Span& span, R2Point *hit_point1, RNScalar *hit_t1)
{
    // Intersection is commutative
    return R2Intersects(span, box, hit_point1, hit_t1);
}



inline RNClassID R2Intersects(const R2Box& box, const R2Halfspace& halfspace)
{
    // Intersection is commutative
    return R2Intersects(halfspace, box);
}



inline RNClassID R2Intersects(const R2Circle& circle, const R2Point& point)
{
    // Intersection is commutative
    return R2Intersects(point, circle);
}



inline RNClassID R2Intersects(const R2Circle& circle, const R2Line& line, R2Point *hit_point)
{
    // Intersection is commutative
    return R2Intersects(line, circle, hit_point);
}



inline RNClassID R2Intersects(const R2Circle& circle, const R2Ray& ray, 
    R2Point *hit_point1, R2Vector *hit_normal1, RNScalar *hit_t1)
{
    // Intersection is commutative
    return R2Intersects(ray, circle, hit_point1, hit_normal1, hit_t1);
}



inline RNClassID R2Intersects(const R2Circle& circle, const R2Span& span, R2Point *hit_point1, RNScalar *hit_t1)
{
    // Intersection is commutative
    return R2Intersects(span, circle, hit_point1, hit_t1);
}



inline RNClassID R2Intersects(const R2Circle& circle, const R2Halfspace& halfspace)
{
    // Intersection is commutative
    return R2Intersects(halfspace, circle);
}



inline RNClassID R2Intersects(const R2Circle& circle, const R2Box& box)
{
    // Intersection is commutative
    return R2Intersects(box, circle);
}




