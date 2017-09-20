/* Include file for GAPS intersection utility */



/* Function declarations */

RNClassID R3Intersects(const R3Point& point1, const R3Point& point2);
RNClassID R3Intersects(const R3Point& point, const R3Line& line);
RNClassID R3Intersects(const R3Point& point, const R3Ray& ray, 
    RNScalar *hit_t = NULL); 
RNClassID R3Intersects(const R3Point& point, const R3Span& span, 
    RNScalar *hit_t = NULL); 
RNClassID R3Intersects(const R3Point& point, const R3Plane& plane); 
RNClassID R3Intersects(const R3Point& point, const R3Triangle& triangle); 
RNClassID R3Intersects(const R3Point& point, const R3Circle& circle); 
RNClassID R3Intersects(const R3Point& point, const R3Ellipse& ellipse); 
RNClassID R3Intersects(const R3Point& point, const R3Rectangle& rectangle); 
RNClassID R3Intersects(const R3Point& point, const R3Halfspace& halfspace); 
RNClassID R3Intersects(const R3Point& point, const R3Box& box); 
RNClassID R3Intersects(const R3Point& point, const R3OrientedBox& box); 
RNClassID R3Intersects(const R3Point& point, const R3Sphere& sphere); 
RNClassID R3Intersects(const R3Point& point, const R3Cylinder& cylinder); 
RNClassID R3Intersects(const R3Point& point, const R3Cone& cone); 
RNClassID R3Intersects(const R3Point& point, const R3Shape& shape);

RNClassID R3Intersects(const R3Line& line, const R3Point& point);
RNClassID R3Intersects(const R3Line& line1, const R3Line& line2, 
    R3Point *hit_point = NULL);
RNClassID R3Intersects(const R3Line& line, const R3Ray& ray, 
    R3Point *hit_point = NULL, RNScalar *hit_t = NULL);
RNClassID R3Intersects(const R3Line& line, const R3Span& span, 
    R3Point *hit_point = NULL, RNScalar *hit_t = NULL);
RNClassID R3Intersects(const R3Line& line, const R3Plane& plane, 
    R3Point *hit_point = NULL);
RNClassID R3Intersects(const R3Line& line, const R3Triangle& triangle, 
    R3Point *hit_point = NULL);
RNClassID R3Intersects(const R3Line& line, const R3Halfspace& halfspace, 
    R3Ray *result = NULL);
RNClassID R3Intersects(const R3Line& line, const R3Box& box, 
    R3Point *hit_point = NULL);
RNClassID R3Intersects(const R3Line& line, const R3OrientedBox& box, 
    R3Point *hit_point = NULL);
RNClassID R3Intersects(const R3Line& line, const R3Sphere& sphere, 
    R3Point *hit_point = NULL);
RNClassID R3Intersects(const R3Line& line, const R3Shape& shape);

RNClassID R3Intersects(const R3Ray& ray, const R3Point& point, 
    RNScalar *hit_t = NULL); 
RNClassID R3Intersects(const R3Ray& ray, const R3Line& line, 
    R3Point *hit_point = NULL, RNScalar *hit_t = NULL);
RNClassID R3Intersects(const R3Ray& ray1, const R3Ray& ray2, 
    R3Point *hit_point = NULL, RNScalar *hit_t1 = NULL, RNScalar *hit_t2 = NULL);
RNClassID R3Intersects(const R3Ray& ray, const R3Span& span, 
    R3Point *hit_point = NULL, RNScalar *hit_tray = NULL, RNScalar *hit_tspan = NULL);
RNClassID R3Intersects(const R3Ray& ray, const R3Plane& plane, 
    R3Point *hit_point = NULL, RNScalar *hit_t = NULL);
RNClassID R3Intersects(const R3Ray& ray, const R3Halfspace& halfspace, 
    R3Point *hit_point = NULL, RNScalar *hit_t = NULL);
RNClassID R3Intersects(const R3Ray& ray, const R3Triangle& triangle, 
    R3Point *hit_point1 = NULL, R3Vector *hit_normal1 = NULL, RNScalar *hit_t1 = NULL);
RNClassID R3Intersects(const R3Ray& ray, const R3TriangleArray& array, 
    R3Point *hit_point1 = NULL, R3Vector *hit_normal1 = NULL, RNScalar *hit_t1 = NULL);
RNClassID R3Intersects(const R3Ray& ray, const R3Circle& circle, 
    R3Point *hit_point1 = NULL, R3Vector *hit_normal1 = NULL, RNScalar *hit_t1 = NULL);
RNClassID R3Intersects(const R3Ray& ray, const R3Ellipse& ellipse, 
    R3Point *hit_point1 = NULL, R3Vector *hit_normal1 = NULL, RNScalar *hit_t1 = NULL);
RNClassID R3Intersects(const R3Ray& ray, const R3Rectangle& rectangle, 
    R3Point *hit_point1 = NULL, R3Vector *hit_normal1 = NULL, RNScalar *hit_t1 = NULL);
RNClassID R3Intersects(const R3Ray& ray, const R3Box& box, 
    R3Point *hit_point1 = NULL, R3Vector *hit_normal1 = NULL, RNScalar *hit_t1 = NULL);
RNClassID R3Intersects(const R3Ray& ray, const R3OrientedBox& box, 
    R3Point *hit_point1 = NULL, R3Vector *hit_normal1 = NULL, RNScalar *hit_t1 = NULL);
RNClassID R3Intersects(const R3Ray& ray, const R3Sphere& sphere, 
    R3Point *hit_point1 = NULL, R3Vector *hit_normal1 = NULL, RNScalar *hit_t1 = NULL);
RNClassID R3Intersects(const R3Ray& ray, const R3Cylinder& cylinder, 
    R3Point *hit_point1 = NULL, R3Vector *hit_norma1l = NULL, RNScalar *hit_t1 = NULL);
RNClassID R3Intersects(const R3Ray& ray, const R3Cone& cone, 
    R3Point *hit_point1 = NULL, R3Vector *hit_normal1 = NULL, RNScalar *hit_t1 = NULL);
RNClassID R3Intersects(const R3Ray& ray, const R3Shape& shape, 
    R3Point *hit_point = NULL, R3Vector *hit_normal = NULL, RNScalar *hit_t = NULL);

RNClassID R3Intersects(const R3Span& span, const R3Point& point, 
    RNScalar *hit_t = NULL);
RNClassID R3Intersects(const R3Span& span, const R3Line& line, 
    R3Point *hit_point = NULL, RNScalar *hit_t = NULL);
RNClassID R3Intersects(const R3Span& span, const R3Ray& ray, 
    R3Point *hit_point = NULL, RNScalar *hit_tspan = NULL, RNScalar *hit_tray = NULL);
RNClassID R3Intersects(const R3Span& span1, const R3Span& span2, 
    R3Point *hit_point = NULL, RNScalar *hit_t1 = NULL, RNScalar *hit_t2 = NULL);
RNClassID R3Intersects(const R3Span& span, const R3Plane& plane, 
    R3Point *hit_point = NULL, RNScalar *hit_t = NULL);
RNClassID R3Intersects(const R3Span& span, const R3Triangle& triangle, 
    R3Point *hit_point1 = NULL, RNScalar *hit_t1 = NULL);
RNClassID R3Intersects(const R3Span& span, const R3Halfspace& halfspace, 
    R3Point *hit_point = NULL, RNScalar *hit_t = NULL);
RNClassID R3Intersects(const R3Span& span, const R3Shape& shape, 
    R3Point *hit_point = NULL, RNScalar *hit_t = NULL);

RNClassID R3Intersects(const R3Plane& plane, const R3Point& point); 
RNClassID R3Intersects(const R3Plane& plane, const R3Line& line, 
    R3Point *hit_point = NULL);
RNClassID R3Intersects(const R3Plane& plane, const R3Ray& ray, 
    R3Point *hit_point = NULL, RNScalar *hit_t = NULL);
RNClassID R3Intersects(const R3Plane& plane, const R3Span& span, 
    R3Point *hit_point = NULL, RNScalar *hit_t = NULL);
RNClassID R3Intersects(const R3Plane& plane1, const R3Plane& plane2, 
    R3Line *result = NULL);
RNClassID R3Intersects(const R3Plane& plane1, const R3Plane& plane2, const R3Plane& plane3, 
    R3Point *hit_point = NULL);
RNClassID R3Intersects(const R3Plane& plane, const R3Halfspace& halfspace);
RNClassID R3Intersects(const R3Plane& plane, const R3Box& box);
RNClassID R3Intersects(const R3Plane& plane, const R3OrientedBox& box);
RNClassID R3Intersects(const R3Plane& plane, const R3Sphere& sphere);
RNClassID R3Intersects(const R3Plane& plane, const R3Ellipsoid& ellipsoid, R3Ellipse *result = NULL);
RNClassID R3Intersects(const R3Plane& plane, const R3Shape& shape);

RNClassID R3Intersects(const R3Triangle& triangle, const R3Point& point); 
RNClassID R3Intersects(const R3Triangle& triangle, const R3Line& line, 
    R3Point *hit_point = NULL);
RNClassID R3Intersects(const R3Triangle& triangle, const R3Ray& ray, 
    R3Point *hit_point = NULL, R3Vector *hit_normal = NULL, RNScalar *hit_t = NULL);
RNClassID R3Intersects(const R3Triangle& triangle, const R3Span& span, 
    R3Point *hit_point1 = NULL, RNScalar *hit_t1 = NULL);
RNClassID R3Intersects(const R3Triangle& triangle, const R3Halfspace& halfspace); 
RNClassID R3Intersects(const R3Triangle& triangle, const R3Box& box);

RNClassID R3Intersects(const R3TriangleArray& array, const R3Ray& ray, 
    R3Point *hit_point1 = NULL, R3Vector *hit_normal1 = NULL, RNScalar *hit_t1 = NULL);

RNClassID R3Intersects(const R3Circle& circle, const R3Ray& ray, 
    R3Point *hit_point1 = NULL, R3Vector *hit_normal = NULL, RNScalar *hit_t1 = NULL);
RNClassID R3Intersects(const R3Circle& circle, const R3Halfspace& halfspace);

RNClassID R3Intersects(const R3Halfspace& halfspace, const R3Point& point);
RNClassID R3Intersects(const R3Halfspace& halfspace, const R3Line& line, 
    R3Ray *result = NULL);
RNClassID R3Intersects(const R3Halfspace& halfspace, const R3Ray& ray, 
    R3Point *hit_point = NULL, R3Vector *hit_normal = NULL, RNScalar *hit_t = NULL);
RNClassID R3Intersects(const R3Halfspace& halfspace, const R3Span& span, 
    R3Point *hit_point = NULL, RNScalar *hit_t = NULL);
RNClassID R3Intersects(const R3Halfspace& halfspace, const R3Plane& plane);
RNClassID R3Intersects(const R3Halfspace& halfspace, const R3Triangle& triangle);
RNClassID R3Intersects(const R3Halfspace& halfspace, const R3Circle& circle);
RNClassID R3Intersects(const R3Halfspace& halfspace1, const R3Halfspace& halfspace2);
RNClassID R3Intersects(const R3Halfspace& halfspace, const R3Box& box);
RNClassID R3Intersects(const R3Halfspace& halfspace, const R3OrientedBox& box);
RNClassID R3Intersects(const R3Halfspace& halfspace, const R3Sphere& sphere);
RNClassID R3Intersects(const R3Halfspace& halfspace, const R3Cylinder& cylinder);
RNClassID R3Intersects(const R3Halfspace& halfspace, const R3Cone& cone);
RNClassID R3Intersects(const R3Halfspace& halfspace, const R3Shape& shape);

RNClassID R3Intersects(const R3Box& box, const R3Point& point);
RNClassID R3Intersects(const R3Box& box, const R3Line& line, 
    R3Point *hit_point = NULL);
RNClassID R3Intersects(const R3Box& box, const R3Ray& ray, 
    R3Point *hit_point1 = NULL, R3Vector *hit_normal = NULL, RNScalar *hit_t1 = NULL);
RNClassID R3Intersects(const R3Box& box, const R3Span& span, 
    R3Point *hit_point1 = NULL, RNScalar *hit_t1 = NULL);
RNClassID R3Intersects(const R3Box& box, const R3Plane& plane);
RNClassID R3Intersects(const R3Box& box, const R3Halfspace& halfspace);
RNClassID R3Intersects(const R3Box& box1, const R3Box& box2, 
    R3Box *result = NULL);
RNClassID R3Intersects(const R3Box& box1, const R3OrientedBox& box2);
RNClassID R3Intersects(const R3Box& box, const R3Sphere& sphere);
RNClassID R3Intersects(const R3Box& box, const R3Shape& shape);

RNClassID R3Intersects(const R3OrientedBox& box, const R3Point& point);
RNClassID R3Intersects(const R3OrientedBox& box, const R3Ray& ray, 
    R3Point *hit_point1 = NULL, R3Vector *hit_normal1 = NULL, RNScalar *hit_t1 = NULL);
RNClassID R3Intersects(const R3OrientedBox& box, const R3Halfspace& halfspace);
RNClassID R3Intersects(const R3OrientedBox& box1, const R3Box& box2);
RNClassID R3Intersects(const R3OrientedBox& box1, const R3OrientedBox& box2);
RNClassID R3Intersects(const R3OrientedBox& box, const R3Sphere& sphere);
RNClassID R3Intersects(const R3OrientedBox& box, const R3Shape& shape);

RNClassID R3Intersects(const R3Sphere& sphere, const R3Point& point); 
RNClassID R3Intersects(const R3Sphere& sphere, const R3Line& line, 
    R3Point *hit_point = NULL);
RNClassID R3Intersects(const R3Sphere& sphere, const R3Ray& ray, 
    R3Point *hit_point1 = NULL, R3Vector *hit_normal = NULL, RNScalar *hit_t1 = NULL);
RNClassID R3Intersects(const R3Sphere& sphere, const R3Span& span, 
    R3Point *hit_point1 = NULL, RNScalar *hit_t1 = NULL);
RNClassID R3Intersects(const R3Sphere& sphere, const R3Plane& plane);
RNClassID R3Intersects(const R3Sphere& sphere, const R3Halfspace& halfspace);
RNClassID R3Intersects(const R3Sphere& sphere, const R3Box& box);
RNClassID R3Intersects(const R3Sphere& sphere, const R3OrientedBox& box);
RNClassID R3Intersects(const R3Sphere& sphere1, const R3Sphere& sphere2);
RNClassID R3Intersects(const R3Sphere& sphere, const R3Shape& shape);

RNClassID R3Intersects(const R3Ellipsoid& ellipsoid, const R3Plane& plane, R3Ellipse *result = NULL);

RNClassID R3Intersects(const R3Cylinder& cylinder, const R3Point& point);
RNClassID R3Intersects(const R3Cylinder& cylinder, const R3Ray& ray, 
    R3Point *hit_point1 = NULL, R3Vector *hit_normal1 = NULL, RNScalar *hit_t1 = NULL);
RNClassID R3Intersects(const R3Cylinder& cylinder, const R3Halfspace& halfspace);

RNClassID R3Intersects(const R3Cone& cone, const R3Point& point);
RNClassID R3Intersects(const R3Cone& cone, const R3Ray& ray, 
    R3Point *hit_point1 = NULL, R3Vector *hit_normal1 = NULL, RNScalar *hit_t1 = NULL);
RNClassID R3Intersects(const R3Cone& cone, const R3Halfspace& halfspace);

RNClassID R3Intersects(const R3Shape& shape, const R3Point& point);
RNClassID R3Intersects(const R3Shape& shape, const R3Line& line);
RNClassID R3Intersects(const R3Shape& shape, const R3Ray& ray, 
    R3Point *hit_point = NULL, R3Vector *hit_normal = NULL, RNScalar *hit_t = NULL);
RNClassID R3Intersects(const R3Shape& shape, const R3Span& span); 
RNClassID R3Intersects(const R3Shape& shape, const R3Plane& plane);
RNClassID R3Intersects(const R3Shape& shape, const R3Halfspace& halfspace);
RNClassID R3Intersects(const R3Shape& shape, const R3Box& box);
RNClassID R3Intersects(const R3Shape& shape, const R3OrientedBox& box);
RNClassID R3Intersects(const R3Shape& shape, const R3Sphere& sphere);
RNClassID R3Intersects(const R3Shape& shape1, const R3Shape& shape2);



/* Inline functions */

inline RNClassID R3Intersects(const R3Line& line, const R3Point& point)
{
    // Intersection is commutative
    return R3Intersects(point, line);
}



inline RNClassID R3Intersects(const R3Ray& ray, const R3Point& point, RNScalar *hit_t)
{
    // Intersection is commutative
    return R3Intersects(point, ray, hit_t);
}



inline RNClassID R3Intersects(const R3Ray& ray, const R3Line& line, R3Point *hit_point, RNScalar *hit_t)
{
    // Intersection is commutative
    return R3Intersects(line, ray, hit_point, hit_t);
}



inline RNClassID R3Intersects(const R3Span& span, const R3Point& point, RNScalar *hit_t)
{
    // Intersection is commutative
    return R3Intersects(point, span, hit_t);
}



inline RNClassID R3Intersects(const R3Span& span, const R3Line& line, R3Point *hit_point, RNScalar *hit_t)
{
    // Intersection is commutative
    return R3Intersects(line, span, hit_point, hit_t);
}



inline RNClassID R3Intersects(const R3Span& span, const R3Ray& ray, R3Point *hit_point, RNScalar *hit_tspan, RNScalar *hit_tray)
{
    // Intersection is commutative
    return R3Intersects(ray, span, hit_point, hit_tray, hit_tspan);
}



inline RNClassID R3Intersects(const R3Plane& plane, const R3Point& point)
{
    // Intersection is commutative
    return R3Intersects(point, plane);
}



inline RNClassID R3Intersects(const R3Plane& plane, const R3Line& line, R3Point *hit_point)
{
    // Intersection is commutative
    return R3Intersects(line, plane, hit_point);
}



inline RNClassID R3Intersects(const R3Plane& plane, const R3Ray& ray, R3Point *hit_point, RNScalar *hit_t)
{
    // Intersection is commutative
    return R3Intersects(ray, plane, hit_point, hit_t);
}



inline RNClassID R3Intersects(const R3Plane& plane, const R3Span& span, R3Point *hit_point, RNScalar *hit_t)
{
    // Intersection is commutative
    return R3Intersects(span, plane, hit_point, hit_t);
}



inline RNClassID R3Intersects(const R3Triangle& triangle, const R3Point& point)
{
    // Intersection is commutative
    return R3Intersects(point, triangle);
}



inline RNClassID R3Intersects(const R3Triangle& triangle, const R3Line& line, 
    R3Point *hit_point)
{
    // Intersection is commutative
    return R3Intersects(line, triangle, hit_point);
}



inline RNClassID R3Intersects(const R3Triangle& triangle, const R3Ray& ray, 
    R3Point *hit_point, R3Vector *hit_normal, RNScalar *hit_t)
{
    // Intersection is commutative
    return R3Intersects(ray, triangle, hit_point, hit_normal, hit_t);
}



inline RNClassID R3Intersects(const R3Triangle& triangle, const R3Span& span, 
    R3Point *hit_point, RNScalar *hit_t)
{
    // Intersection is commutative
    return R3Intersects(span, triangle, hit_point, hit_t);
}



inline RNClassID R3Intersects(const R3TriangleArray& array, const R3Ray& ray, 
    R3Point *hit_point1, R3Vector *hit_normal1, RNScalar *hit_t1)
{
    // Intersection is commutative
    return R3Intersects(ray, array, hit_point1, hit_normal1, hit_t1);
}



inline RNClassID R3Intersects(const R3Circle& circle, const R3Ray& ray, 
    R3Point *hit_point, R3Vector *hit_normal, RNScalar *hit_t)
{
    // Intersection is commutative
    return R3Intersects(ray, circle, hit_point, hit_normal, hit_t);
}



inline RNClassID R3Intersects(const R3Circle& circle, const R3Halfspace& halfspace)
{
    // Intersection is commutative
    return R3Intersects(halfspace, circle);
}



inline RNClassID R3Intersects(const R3Halfspace& halfspace, const R3Point& point)
{
    // Intersection is commutative
    return R3Intersects(point, halfspace);
}



inline RNClassID R3Intersects(const R3Halfspace& halfspace, const R3Line& line, R3Ray *result)
{
    // Intersection is commutative
    return R3Intersects(line, halfspace, result);
}



inline RNClassID R3Intersects(const R3Halfspace& halfspace, const R3Ray& ray, R3Point *hit_point, RNScalar *hit_t)
{
    // Intersection is commutative
    return R3Intersects(ray, halfspace, hit_point, hit_t);
}



inline RNClassID R3Intersects(const R3Halfspace& halfspace, const R3Span& span, R3Point *hit_point, RNScalar *hit_t)
{
    // Intersection is commutative
    return R3Intersects(span, halfspace, hit_point, hit_t);
}



inline RNClassID R3Intersects(const R3Halfspace& halfspace, const R3Plane& plane)
{
    // Intersection is commutative
    return R3Intersects(plane, halfspace);
}



inline RNClassID R3Intersects(const R3Halfspace& halfspace, const R3Triangle& triangle)
{
    // Intersection is commutative
    return R3Intersects(triangle, halfspace);
}



inline RNClassID R3Intersects(const R3Box& box, const R3Point& point)
{
    // Intersection is commutative
    return R3Intersects(point, box);
}



inline RNClassID R3Intersects(const R3Box& box, const R3Line& line, R3Point *hit_point)
{
    // Intersection is commutative
    return R3Intersects(line, box, hit_point);
}



inline RNClassID R3Intersects(const R3Box& box, const R3Ray& ray, 
    R3Point *hit_point1, R3Vector *hit_normal1, RNScalar *hit_t1)
{
    // Intersection is commutative
    return R3Intersects(ray, box, hit_point1, hit_normal1, hit_t1);
}



inline RNClassID R3Intersects(const R3Box& box, const R3Span& span, R3Point *hit_point1, RNScalar *hit_t1)
{
    // Intersection is commutative
    return R3Intersects(span, box, hit_point1, hit_t1);
}



inline RNClassID R3Intersects(const R3Box& box, const R3Plane& plane)
{
    // Intersection is commutative
    return R3Intersects(plane, box);
}



inline RNClassID R3Intersects(const R3Box& box, const R3Halfspace& halfspace)
{
    // Intersection is commutative
    return R3Intersects(halfspace, box);
}



inline RNClassID R3Intersects(const R3Box& box, const R3Triangle& triangle)
{
    // Intersection is commutative
    return R3Intersects(triangle, box);
}



inline RNClassID R3Intersects(const R3OrientedBox& box, const R3Point& point)
{
    // Intersection is commutative
    return R3Intersects(point, box);
}



inline RNClassID R3Intersects(const R3OrientedBox& box, const R3Line& line)
{
    // Intersection is commutative
    return R3Intersects(line, box);
}



inline RNClassID R3Intersects(const R3OrientedBox& box, const R3Ray& ray, R3Point *hit_point1, R3Vector *hit_normal1, RNScalar *hit_t1)
{
    // Intersection is commutative
    return R3Intersects(ray, box, hit_point1, hit_normal1, hit_t1);
}



inline RNClassID R3Intersects(const R3OrientedBox& box, const R3Plane& plane)
{
    // Intersection is commutative
    return R3Intersects(plane, box);
}



inline RNClassID R3Intersects(const R3OrientedBox& box, const R3Halfspace& halfspace)
{
    // Intersection is commutative
    return R3Intersects(halfspace, box);
}



inline RNClassID R3Intersects(const R3OrientedBox& box1, const R3Box& box2)
{
    // Intersection is commutative
    return R3Intersects(box2, box1);
}



inline RNClassID R3Intersects(const R3Sphere& sphere, const R3Point& point)
{
    // Intersection is commutative
    return R3Intersects(point, sphere);
}



inline RNClassID R3Intersects(const R3Sphere& sphere, const R3Line& line, R3Point *hit_point)
{
    // Intersection is commutative
    return R3Intersects(line, sphere, hit_point);
}



inline RNClassID R3Intersects(const R3Sphere& sphere, const R3Ray& ray, 
    R3Point *hit_point1, R3Vector *hit_normal1, RNScalar *hit_t1)
{
    // Intersection is commutative
    return R3Intersects(ray, sphere, hit_point1, hit_normal1, hit_t1);
}



inline RNClassID R3Intersects(const R3Sphere& sphere, const R3Span& span, R3Point *hit_point1, RNScalar *hit_t1)
{
    // Intersection is commutative
    return R3Intersects(span, sphere, hit_point1, hit_t1);
}



inline RNClassID R3Intersects(const R3Sphere& sphere, const R3Plane& plane)
{
    // Intersection is commutative
    return R3Intersects(plane, sphere);
}



inline RNClassID R3Intersects(const R3Sphere& sphere, const R3Halfspace& halfspace)
{
    // Intersection is commutative
    return R3Intersects(halfspace, sphere);
}



inline RNClassID R3Intersects(const R3Sphere& sphere, const R3Box& box)
{
    // Intersection is commutative
    return R3Intersects(box, sphere);
}



inline RNClassID R3Intersects(const R3Sphere& sphere, const R3OrientedBox& box)
{
    // Intersection is commutative
    return R3Intersects(box, sphere);
}



inline RNClassID R3Intersects(const R3Ellipsoid& ellipsoid, const R3Plane& plane, R3Ellipse *result)
{
    // Intersection is commutative
    return R3Intersects(plane, ellipsoid, result);
}



inline RNClassID R3Intersects(const R3Cylinder& cylinder, const R3Point& point)
{
    // Intersection is commutative
    return R3Intersects(point, cylinder);
}



inline RNClassID R3Intersects(const R3Cylinder& cylinder, const R3Ray& ray, 
    R3Point *hit_point1, R3Vector *hit_normal1, RNScalar *hit_t1)
{
    // Intersection is commutative
    return R3Intersects(ray, cylinder, hit_point1, hit_normal1, hit_t1);
}



inline RNClassID R3Intersects(const R3Cylinder& cylinder, const R3Halfspace& halfspace)
{
    // Intersection is commutative
    return R3Intersects(halfspace, cylinder);
}



inline RNClassID R3Intersects(const R3Cone& cone, const R3Point& point)
{
    // Intersection is commutative
    return R3Intersects(point, cone);
}



inline RNClassID R3Intersects(const R3Cone& cone, const R3Ray& ray, 
    R3Point *hit_point1, R3Vector *hit_normal1, RNScalar *hit_t1)
{
    // Intersection is commutative
    return R3Intersects(ray, cone, hit_point1, hit_normal1, hit_t1);
}



inline RNClassID R3Intersects(const R3Cone& cone, const R3Halfspace& halfspace)
{
    // Intersection is commutative
    return R3Intersects(halfspace, cone);
}



inline RNClassID R3Intersects(const R3Shape& shape, const R3Point& point)
{
    // Intersection is commutative
    return R3Intersects(point, shape);
}



inline RNClassID R3Intersects(const R3Shape& shape, const R3Line& line)
{
    // Intersection is commutative
    return R3Intersects(line, shape);
}



inline RNClassID R3Intersects(const R3Shape& shape, const R3Ray& ray, 
    R3Point *hit_point, R3Vector *hit_normal, RNScalar *hit_t)
{
    // Intersection is commutative
    return R3Intersects(ray, shape, hit_point, hit_normal, hit_t);
}



inline RNClassID R3Intersects(const R3Shape& shape, const R3Span& span)
{
    // Intersection is commutative
    return R3Intersects(span, shape);
}



inline RNClassID R3Intersects(const R3Shape& shape, const R3Plane& plane)
{
    // Intersection is commutative
    return R3Intersects(plane, shape);
}



inline RNClassID R3Intersects(const R3Shape& shape, const R3Halfspace& halfspace)
{
    // Intersection is commutative
    return R3Intersects(halfspace, shape);
}



inline RNClassID R3Intersects(const R3Shape& shape, const R3Box& box)
{
    // Intersection is commutative
    return R3Intersects(box, shape);
}



inline RNClassID R3Intersects(const R3Shape& shape, const R3OrientedBox& box)
{
    // Intersection is commutative
    return R3Intersects(box, shape);
}



inline RNClassID R3Intersects(const R3Shape& shape, const R3Sphere& sphere)
{
    // Intersection is commutative
    return R3Intersects(sphere, shape);
}







