/* Include file for GAPS distance utility */



/* Function declarations */

RNLength R3Distance(const R3Point& point1, const R3Point& point2);
RNLength R3Distance(const R3Point& point, const R3Line& line);
RNLength R3Distance(const R3Point& point, const R3Ray& ray);
RNLength R3Distance(const R3Point& point, const R3Span& span);
RNLength R3Distance(const R3Point& point, const R3Plane& plane);
RNLength R3Distance(const R3Point& point, const R3Halfspace& halfspace);
RNLength R3Distance(const R3Point& point, const R3Box& box);
RNLength R3Distance(const R3Point& point, const R3OrientedBox& box);
RNLength R3Distance(const R3Point& point, const R3Sphere& sphere);
RNLength R3Distance(const R3Point& point, const R3Cone& cone);
RNLength R3Distance(const R3Point& point, const R3Shape& shape);
RNLength R3SquaredDistance(const R3Point& point1, const R3Point& point2);

RNLength R3Distance(const R3Line& line, const R3Point& point);
RNLength R3Distance(const R3Line& line1, const R3Line& line2);
RNLength R3Distance(const R3Line& line, const R3Ray& ray);
RNLength R3Distance(const R3Line& line, const R3Span& span);
RNLength R3Distance(const R3Line& line, const R3Plane& plane);
RNLength R3Distance(const R3Line& line, const R3Halfspace& halfspace);
RNLength R3Distance(const R3Line& line, const R3Box& box);
RNLength R3Distance(const R3Line& line, const R3OrientedBox& box);
RNLength R3Distance(const R3Line& line, const R3Sphere& sphere);
RNLength R3Distance(const R3Line& line, const R3Shape& shape);

RNLength R3Distance(const R3Ray& ray, const R3Point& point);
RNLength R3Distance(const R3Ray& ray, const R3Line& line);
RNLength R3Distance(const R3Ray& ray1, const R3Ray& ray2);
RNLength R3Distance(const R3Ray& ray, const R3Span& span);
RNLength R3Distance(const R3Ray& ray, const R3Plane& plane);
RNLength R3Distance(const R3Ray& ray, const R3Halfspace& halfspace);
RNLength R3Distance(const R3Ray& ray, const R3Box& box);
RNLength R3Distance(const R3Ray& ray, const R3OrientedBox& box);
RNLength R3Distance(const R3Ray& ray, const R3Sphere& sphere);
RNLength R3Distance(const R3Ray& ray, const R3Shape& shape);

RNLength R3Distance(const R3Span& span, const R3Point& point);
RNLength R3Distance(const R3Span& span, const R3Line& line);
RNLength R3Distance(const R3Span& span, const R3Ray& ray);
RNLength R3Distance(const R3Span& span1, const R3Span& span2);
RNLength R3Distance(const R3Span& span, const R3Plane& plane);
RNLength R3Distance(const R3Span& span, const R3Halfspace& halfspace);
RNLength R3Distance(const R3Span& span, const R3Box& box);
RNLength R3Distance(const R3Span& span, const R3OrientedBox& box);
RNLength R3Distance(const R3Span& span, const R3Sphere& sphere);
RNLength R3Distance(const R3Span& span, const R3Shape& shape);

RNLength R3Distance(const R3Plane& plane, const R3Point& point);
RNLength R3Distance(const R3Plane& plane, const R3Line& line);
RNLength R3Distance(const R3Plane& plane, const R3Ray& ray);
RNLength R3Distance(const R3Plane& plane, const R3Span& span);
RNLength R3Distance(const R3Plane& plane1, const R3Plane& plane2);
RNLength R3Distance(const R3Plane& plane, const R3Halfspace& halfspace);
RNLength R3Distance(const R3Plane& plane, const R3Box& box);
RNLength R3Distance(const R3Plane& plane, const R3OrientedBox& box);
RNLength R3Distance(const R3Plane& plane, const R3Sphere& sphere);
RNLength R3Distance(const R3Plane& plane, const R3Shape& shape);

RNLength R3SignedDistance(const R3Plane& plane, const R3Point& point);
RNLength R3SignedDistance(const R3Plane& plane, const R3Line& line);
RNLength R3SignedDistance(const R3Plane& plane, const R3Ray& ray);
RNLength R3SignedDistance(const R3Plane& plane, const R3Span& span);
RNLength R3SignedDistance(const R3Plane& plane1, const R3Plane& plane2);
RNLength R3SignedDistance(const R3Plane& plane, const R3Triangle& triangle);
RNLength R3SignedDistance(const R3Plane& plane, const R3Halfspace& halfspace);
RNLength R3SignedDistance(const R3Plane& plane, const R3Box& box);
RNLength R3SignedDistance(const R3Plane& plane, const R3OrientedBox& box);
RNLength R3SignedDistance(const R3Plane& plane, const R3Sphere& sphere);
RNLength R3SignedDistance(const R3Plane& plane, const R3Shape& shape);

RNLength R3Distance(const R3Halfspace& halfspace, const R3Point& point);
RNLength R3Distance(const R3Halfspace& halfspace, const R3Line& line);
RNLength R3Distance(const R3Halfspace& halfspace, const R3Ray& ray);
RNLength R3Distance(const R3Halfspace& halfspace, const R3Span& span);
RNLength R3Distance(const R3Halfspace& halfspace, const R3Plane& plane);
RNLength R3Distance(const R3Halfspace& halfspace1, const R3Halfspace& halfspace2);
RNLength R3Distance(const R3Halfspace& halfspace, const R3Box& box);
RNLength R3Distance(const R3Halfspace& halfspace, const R3OrientedBox& box);
RNLength R3Distance(const R3Halfspace& halfspace, const R3Sphere& sphere);
RNLength R3Distance(const R3Halfspace& halfspace, const R3Shape& shape);

RNLength R3Distance(const R3Box& box, const R3Point& point);
RNLength R3Distance(const R3Box& box, const R3Line& line);
RNLength R3Distance(const R3Box& box, const R3Ray& ray);
RNLength R3Distance(const R3Box& box, const R3Span& span);
RNLength R3Distance(const R3Box& box, const R3Plane& plane);
RNLength R3Distance(const R3Box& box, const R3Halfspace& halfspace);
RNLength R3Distance(const R3Box& box1, const R3Box& box2);
RNLength R3Distance(const R3Box& box1, const R3OrientedBox& box2);
RNLength R3Distance(const R3Box& box, const R3Sphere& sphere);
RNLength R3Distance(const R3Box& box, const R3Shape& shape);

RNLength R3Distance(const R3OrientedBox& box, const R3Point& point);
RNLength R3Distance(const R3OrientedBox& box, const R3Line& line);
RNLength R3Distance(const R3OrientedBox& box, const R3Ray& ray);
RNLength R3Distance(const R3OrientedBox& box, const R3Span& span);
RNLength R3Distance(const R3OrientedBox& box, const R3Plane& plane);
RNLength R3Distance(const R3OrientedBox& box, const R3Halfspace& halfspace);
RNLength R3Distance(const R3OrientedBox& box1, const R3Box& box2);
RNLength R3Distance(const R3OrientedBox& box1, const R3OrientedBox& box2);
RNLength R3Distance(const R3OrientedBox& box, const R3Sphere& sphere);
RNLength R3Distance(const R3OrientedBox& box, const R3Shape& shape);

RNLength R3Distance(const R3Sphere& sphere, const R3Point& point);
RNLength R3Distance(const R3Sphere& sphere, const R3Line& line);
RNLength R3Distance(const R3Sphere& sphere, const R3Ray& ray);
RNLength R3Distance(const R3Sphere& sphere, const R3Span& span);
RNLength R3Distance(const R3Sphere& sphere, const R3Plane& plane);
RNLength R3Distance(const R3Sphere& sphere, const R3Halfspace& halfspace);
RNLength R3Distance(const R3Sphere& sphere, const R3Box& box);
RNLength R3Distance(const R3Sphere& sphere, const R3OrientedBox& box);
RNLength R3Distance(const R3Sphere& sphere1, const R3Sphere& sphere2);
RNLength R3Distance(const R3Sphere& sphere, const R3Shape& shape);

RNLength R3Distance(const R3Cone& cone, const R3Point& point);

RNLength R3Distance(const R3Shape& shape, const R3Point& point);
RNLength R3Distance(const R3Shape& shape, const R3Line& line);
RNLength R3Distance(const R3Shape& shape, const R3Ray& ray);
RNLength R3Distance(const R3Shape& shape, const R3Span& span);
RNLength R3Distance(const R3Shape& shape, const R3Plane& plane);
RNLength R3Distance(const R3Shape& shape, const R3Halfspace& halfspace);
RNLength R3Distance(const R3Shape& shape, const R3Box& box);
RNLength R3Distance(const R3Shape& shape, const R3OrientedBox& box);
RNLength R3Distance(const R3Shape& shape, const R3Sphere& sphere);
RNLength R3Distance(const R3Shape& shape1, const R3Shape& shape2);



/* Abstract shape function type definition */

#ifdef TYPES
#define R3_NUM_SHAPE_TYPES 32
typedef RNLength R3DistanceFunction(const R3Shape& shape1, const R3Shape& shape2);
#endif



/* Inline functions */

inline RNLength R3Distance(const R3Line& line, const R3Point& point)
{
    // Distance is commutative
    return R3Distance(point, line);
}



inline RNLength R3Distance(const R3Ray& ray, const R3Point& point)
{
    // Distance is commutative
    return R3Distance(point, ray);
}



inline RNLength R3Distance(const R3Ray& ray, const R3Line& line)
{
    // Distance is commutative
    return R3Distance(line, ray);
}



inline RNLength R3Distance(const R3Span& span, const R3Point& point)
{
    // Distance is commutative
    return R3Distance(point, span);
}



inline RNLength R3Distance(const R3Span& span, const R3Line& line)
{
    // Distance is commutative
    return R3Distance(line, span);
}



inline RNLength R3Distance(const R3Span& span, const R3Ray& ray)
{
    // Distance is commutative
    return R3Distance(ray, span);
}



inline RNLength R3Distance(const R3Plane& plane, const R3Point& point)
{
    // Distance is commutative
    return R3Distance(point, plane);
}



inline RNLength R3Distance(const R3Plane& plane, const R3Line& line)
{
    // Distance is commutative
    return R3Distance(line, plane);
}



inline RNLength R3Distance(const R3Plane& plane, const R3Ray& ray)
{
    // Distance is commutative
    return R3Distance(ray, plane);
}



inline RNLength R3Distance(const R3Plane& plane, const R3Span& span)
{
    // Distance is commutative
    return R3Distance(span, plane);
}



inline RNLength R3Distance(const R3Halfspace& halfspace, const R3Point& point)
{
    // Distance is commutative
    return R3Distance(point, halfspace);
}



inline RNLength R3Distance(const R3Halfspace& halfspace, const R3Line& line)
{
    // Distance is commutative
    return R3Distance(line, halfspace);
}



inline RNLength R3Distance(const R3Halfspace& halfspace, const R3Ray& ray)
{
    // Distance is commutative
    return R3Distance(ray, halfspace);
}



inline RNLength R3Distance(const R3Halfspace& halfspace, const R3Span& span)
{
    // Distance is commutative
    return R3Distance(span, halfspace);
}



inline RNLength R3Distance(const R3Halfspace& halfspace, const R3Plane& plane)
{
    // Distance is commutative
    return R3Distance(plane, halfspace);
}



inline RNLength R3Distance(const R3Box& box, const R3Point& point)
{
    // Distance is commutative
    return R3Distance(point, box);
}



inline RNLength R3Distance(const R3Box& box, const R3Line& line)
{
    // Distance is commutative
    return R3Distance(line, box);
}



inline RNLength R3Distance(const R3Box& box, const R3Ray& ray)
{
    // Distance is commutative
    return R3Distance(ray, box);
}



inline RNLength R3Distance(const R3Box& box, const R3Span& span)
{
    // Distance is commutative
    return R3Distance(span, box);
}



inline RNLength R3Distance(const R3Box& box, const R3Plane& plane)
{
    // Distance is commutative
    return R3Distance(plane, box);
}



inline RNLength R3Distance(const R3Box& box, const R3Halfspace& halfspace)
{
    // Distance is commutative
    return R3Distance(halfspace, box);
}



inline RNLength R3Distance(const R3OrientedBox& box, const R3Point& point)
{
    // Distance is commutative
    return R3Distance(point, box);
}



inline RNLength R3Distance(const R3OrientedBox& box, const R3Line& line)
{
    // Distance is commutative
    return R3Distance(line, box);
}



inline RNLength R3Distance(const R3OrientedBox& box, const R3Ray& ray)
{
    // Distance is commutative
    return R3Distance(ray, box);
}



inline RNLength R3Distance(const R3OrientedBox& box, const R3Span& span)
{
    // Distance is commutative
    return R3Distance(span, box);
}



inline RNLength R3Distance(const R3OrientedBox& box, const R3Plane& plane)
{
    // Distance is commutative
    return R3Distance(plane, box);
}



inline RNLength R3Distance(const R3OrientedBox& box, const R3Halfspace& halfspace)
{
    // Distance is commutative
    return R3Distance(halfspace, box);
}



inline RNLength R3Distance(const R3OrientedBox& box1, const R3Box& box2)
{
    // Distance is commutative
    return R3Distance(box2, box1);
}



inline RNLength R3Distance(const R3Sphere& sphere, const R3Point& point)
{
    // Distance is commutative
    return R3Distance(point, sphere);
}



inline RNLength R3Distance(const R3Sphere& sphere, const R3Line& line)
{
    // Distance is commutative
    return R3Distance(line, sphere);
}



inline RNLength R3Distance(const R3Sphere& sphere, const R3Ray& ray)
{
    // Distance is commutative
    return R3Distance(ray, sphere);
}



inline RNLength R3Distance(const R3Sphere& sphere, const R3Span& span)
{
    // Distance is commutative
    return R3Distance(span, sphere);
}



inline RNLength R3Distance(const R3Sphere& sphere, const R3Plane& plane)
{
    // Distance is commutative
    return R3Distance(plane, sphere);
}



inline RNLength R3Distance(const R3Sphere& sphere, const R3Halfspace& halfspace)
{
    // Distance is commutative
    return R3Distance(halfspace, sphere);
}



inline RNLength R3Distance(const R3Sphere& sphere, const R3Box& box)
{
    // Distance is commutative
    return R3Distance(box, sphere);
}



inline RNLength R3Distance(const R3Sphere& sphere, const R3OrientedBox& box)
{
    // Distance is commutative
    return R3Distance(box, sphere);
}



inline RNLength R3Distance(const R3Cone& cone, const R3Point& point)
{
    // Distance is commutative
    return R3Distance(point, cone);
}



inline RNLength R3Distance(const R3Shape& shape, const R3Point& point)
{
    // Distance is commutative
    return R3Distance(point, shape);
}



inline RNLength R3Distance(const R3Shape& shape, const R3Line& line)
{
    // Distance is commutative
    return R3Distance(line, shape);
}



inline RNLength R3Distance(const R3Shape& shape, const R3Ray& ray)
{
    // Distance is commutative
    return R3Distance(ray, shape);
}



inline RNLength R3Distance(const R3Shape& shape, const R3Span& span)
{
    // Distance is commutative
    return R3Distance(span, shape);
}



inline RNLength R3Distance(const R3Shape& shape, const R3Plane& plane)
{
    // Distance is commutative
    return R3Distance(plane, shape);
}



inline RNLength R3Distance(const R3Shape& shape, const R3Halfspace& halfspace)
{
    // Distance is commutative
    return R3Distance(halfspace, shape);
}



inline RNLength R3Distance(const R3Shape& shape, const R3Box& box)
{
    // Distance is commutative
    return R3Distance(box, shape);
}



inline RNLength R3Distance(const R3Shape& shape, const R3OrientedBox& box)
{
    // Distance is commutative
    return R3Distance(box, shape);
}



inline RNLength R3Distance(const R3Shape& shape, const R3Sphere& sphere)
{
    // Distance is commutative
    return R3Distance(sphere, shape);
}





