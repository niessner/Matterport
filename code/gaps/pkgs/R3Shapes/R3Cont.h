/* Include file for GAPS inside/containment utility */



/* Useful macros */

#define R3Inside(__primitive1, __primitive2) \
    R3Contains(__primitive2, __primitive1)



/* Function declarations */

RNBoolean R3Contains(const R3Vector& vector1, const R3Vector& vector2);

RNBoolean R3Contains(const R3Point& point1, const R3Point& point2);
RNBoolean R3Contains(const R3Point& point, const R3Line& line);
RNBoolean R3Contains(const R3Point& point, const R3Ray& ray);
RNBoolean R3Contains(const R3Point& point, const R3Span& span);
RNBoolean R3Contains(const R3Point& point, const R3Plane& plane);
RNBoolean R3Contains(const R3Point& point, const R3Triangle& triangle);
RNBoolean R3Contains(const R3Point& point, const R3Circle& circle);
RNBoolean R3Contains(const R3Point& point, const R3Ellipse& ellipse);
RNBoolean R3Contains(const R3Point& point, const R3Rectangle& rectangle);
RNBoolean R3Contains(const R3Point& point, const R3Halfspace& halfspace);
RNBoolean R3Contains(const R3Point& point, const R3Box& box);
RNBoolean R3Contains(const R3Point& point, const R3OrientedBox& box);
RNBoolean R3Contains(const R3Point& point, const R3Sphere& sphere);
RNBoolean R3Contains(const R3Point& point, const R3Shape& shape);

RNBoolean R3Contains(const R3Line& line, const R3Point& point);
RNBoolean R3Contains(const R3Line& line1, const R3Line& line2);
RNBoolean R3Contains(const R3Line& line, const R3Ray& ray);
RNBoolean R3Contains(const R3Line& line, const R3Span& span);
RNBoolean R3Contains(const R3Line& line, const R3Plane& plane);
RNBoolean R3Contains(const R3Line& line, const R3Halfspace& halfspace);
RNBoolean R3Contains(const R3Line& line, const R3Box& box);
RNBoolean R3Contains(const R3Line& line, const R3OrientedBox& box);
RNBoolean R3Contains(const R3Line& line, const R3Sphere& sphere);
RNBoolean R3Contains(const R3Line& line, const R3Shape& shape);

RNBoolean R3Contains(const R3Ray& ray, const R3Point& point);
RNBoolean R3Contains(const R3Ray& ray, const R3Line& line);
RNBoolean R3Contains(const R3Ray& ray1, const R3Ray& ray2);
RNBoolean R3Contains(const R3Ray& ray, const R3Span& span);
RNBoolean R3Contains(const R3Ray& ray, const R3Plane& plane);
RNBoolean R3Contains(const R3Ray& ray, const R3Halfspace& halfspace);
RNBoolean R3Contains(const R3Ray& ray, const R3Box& box);
RNBoolean R3Contains(const R3Ray& ray, const R3OrientedBox& box);
RNBoolean R3Contains(const R3Ray& ray, const R3Sphere& sphere);
RNBoolean R3Contains(const R3Ray& ray, const R3Shape& shape);

RNBoolean R3Contains(const R3Span& span, const R3Point& point);
RNBoolean R3Contains(const R3Span& span, const R3Line& line);
RNBoolean R3Contains(const R3Span& span, const R3Ray& ray);
RNBoolean R3Contains(const R3Span& span1, const R3Span& span2);
RNBoolean R3Contains(const R3Span& span, const R3Plane& plane);
RNBoolean R3Contains(const R3Span& span, const R3Halfspace& halfspace);
RNBoolean R3Contains(const R3Span& span, const R3Box& box);
RNBoolean R3Contains(const R3Span& span, const R3OrientedBox& box);
RNBoolean R3Contains(const R3Span& span, const R3Sphere& sphere);
RNBoolean R3Contains(const R3Span& span, const R3Shape& shape);

RNBoolean R3Contains(const R3Plane& plane, const R3Point& point);
RNBoolean R3Contains(const R3Plane& plane, const R3Line& line);
RNBoolean R3Contains(const R3Plane& plane, const R3Ray& ray);
RNBoolean R3Contains(const R3Plane& plane, const R3Span& span);
RNBoolean R3Contains(const R3Plane& plane, const R3Triangle& triangle);
RNBoolean R3Contains(const R3Plane& plane1, const R3Plane& plane2);
RNBoolean R3Contains(const R3Plane& plane, const R3Halfspace& halfspace);
RNBoolean R3Contains(const R3Plane& plane, const R3Box& box);
RNBoolean R3Contains(const R3Plane& plane, const R3OrientedBox& box);
RNBoolean R3Contains(const R3Plane& plane, const R3Sphere& sphere);
RNBoolean R3Contains(const R3Plane& plane, const R3Shape& shape);

RNBoolean R3Contains(const R3Triangle& triangle, const R3Point& point);
RNBoolean R3Contains(const R3Triangle& triangle, const R3Line& line);
RNBoolean R3Contains(const R3Triangle& triangle, const R3Ray& ray);
RNBoolean R3Contains(const R3Triangle& triangle, const R3Span& span);
RNBoolean R3Contains(const R3Triangle& triangle, const R3Plane& plane);
RNBoolean R3Contains(const R3Triangle& triangle1, const R3Triangle& triangle2);
RNBoolean R3Contains(const R3Triangle& triangle, const R3Halfspace& halfspace);
RNBoolean R3Contains(const R3Triangle& triangle, const R3Box& box);
RNBoolean R3Contains(const R3Triangle& triangle, const R3OrientedBox& box);
RNBoolean R3Contains(const R3Triangle& triangle, const R3Sphere& sphere);
RNBoolean R3Contains(const R3Triangle& triangle, const R3Shape& shape);

RNBoolean R3Contains(const R3Circle& circle, const R3Point& point);

RNBoolean R3Contains(const R3Ellipse& ellipse, const R3Point& point);

RNBoolean R3Contains(const R3Rectangle& rectangle, const R3Point& point);

RNBoolean R3Contains(const R3Halfspace& halfspace, const R3Point& point);
RNBoolean R3Contains(const R3Halfspace& halfspace, const R3Line& line);
RNBoolean R3Contains(const R3Halfspace& halfspace, const R3Ray& ray);
RNBoolean R3Contains(const R3Halfspace& halfspace, const R3Span& span);
RNBoolean R3Contains(const R3Halfspace& halfspace, const R3Plane& plane);
RNBoolean R3Contains(const R3Halfspace& halfspace, const R3Triangle& triangle);
RNBoolean R3Contains(const R3Halfspace& halfspace, const R3Circle& circle);
RNBoolean R3Contains(const R3Halfspace& halfspace1, const R3Halfspace& halfspace2);
RNBoolean R3Contains(const R3Halfspace& halfspace, const R3Box& box);
RNBoolean R3Contains(const R3Halfspace& halfspace, const R3OrientedBox& box);
RNBoolean R3Contains(const R3Halfspace& halfspace, const R3Sphere& sphere);
RNBoolean R3Contains(const R3Halfspace& halfspace, const R3Cylinder& cylinder);
RNBoolean R3Contains(const R3Halfspace& halfspace, const R3Cone& cone);
RNBoolean R3Contains(const R3Halfspace& halfspace, const R3Shape& shape);

RNBoolean R3Contains(const R3Box& box, const R3Point& point);
RNBoolean R3Contains(const R3Box& box, const R3Line& line);
RNBoolean R3Contains(const R3Box& box, const R3Ray& ray);
RNBoolean R3Contains(const R3Box& box, const R3Span& span);
RNBoolean R3Contains(const R3Box& box, const R3Plane& plane);
RNBoolean R3Contains(const R3Box& box, const R3Halfspace& halfspace);
RNBoolean R3Contains(const R3Box& box1, const R3Box& box2);
RNBoolean R3Contains(const R3Box& box1, const R3OrientedBox& box2);
RNBoolean R3Contains(const R3Box& box, const R3Sphere& sphere);
RNBoolean R3Contains(const R3Box& box, const R3Shape& shape);

RNBoolean R3Contains(const R3OrientedBox& box, const R3Point& point);
RNBoolean R3Contains(const R3OrientedBox& box, const R3Line& point);
RNBoolean R3Contains(const R3OrientedBox& box, const R3Ray& point);
RNBoolean R3Contains(const R3OrientedBox& box, const R3Span& span);
RNBoolean R3Contains(const R3OrientedBox& box, const R3Plane& plane);
RNBoolean R3Contains(const R3OrientedBox& box, const R3Halfspace& halfspace);
RNBoolean R3Contains(const R3OrientedBox& box1, const R3Box& box2);
RNBoolean R3Contains(const R3OrientedBox& box1, const R3OrientedBox& box2);
RNBoolean R3Contains(const R3OrientedBox& box, const R3Sphere& sphere);
RNBoolean R3Contains(const R3OrientedBox& box, const R3Shape& shape);

RNBoolean R3Contains(const R3Sphere& sphere, const R3Point& point);
RNBoolean R3Contains(const R3Sphere& sphere, const R3Line& line);
RNBoolean R3Contains(const R3Sphere& sphere, const R3Ray& ray);
RNBoolean R3Contains(const R3Sphere& sphere, const R3Span& span);
RNBoolean R3Contains(const R3Sphere& sphere, const R3Plane& plane);
RNBoolean R3Contains(const R3Sphere& sphere, const R3Halfspace& halfspace);
RNBoolean R3Contains(const R3Sphere& sphere, const R3Box& box);
RNBoolean R3Contains(const R3Sphere& sphere, const R3OrientedBox& box);
RNBoolean R3Contains(const R3Sphere& sphere1, const R3Sphere& sphere2);
RNBoolean R3Contains(const R3Sphere& sphere, const R3Shape& shape);

RNBoolean R3Contains(const R3Cylinder& cylinder, const R3Point& point);

RNBoolean R3Contains(const R3Cone& cone, const R3Point& point);

RNBoolean R3Contains(const R3Shape& shape, const R3Point& point);
RNBoolean R3Contains(const R3Shape& shape, const R3Line& line);
RNBoolean R3Contains(const R3Shape& shape, const R3Ray& ray);
RNBoolean R3Contains(const R3Shape& shape, const R3Span& span);
RNBoolean R3Contains(const R3Shape& shape, const R3Plane& plane);
RNBoolean R3Contains(const R3Shape& shape, const R3Halfspace& halfspace);
RNBoolean R3Contains(const R3Shape& shape, const R3Box& box);
RNBoolean R3Contains(const R3Shape& shape, const R3OrientedBox& box);
RNBoolean R3Contains(const R3Shape& shape, const R3Sphere& sphere);
RNBoolean R3Contains(const R3Shape& shape1, const R3Shape& shape2);


