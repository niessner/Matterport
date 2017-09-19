/* Source file for the intersection utility */



/* Include files */

#include "R3Shapes/R3Shapes.h"



RNClassID R3Intersects(const R3Point& point1, const R3Point& point2)
{
    // Check if two points are the same within tolerance
    if (R3Contains(point1, point2)) {
	// Points are same
	return R3_POINT_CLASS_ID;
    }
    else {
	// Points are different
	return RN_NULL_CLASS_ID;
    }
}



RNClassID R3Intersects(const R3Point& point, const R3Line& line)
{
    // Check if line contains point
    if (R3Contains(line, point)) {
	// Point lies on line
	return R3_POINT_CLASS_ID;
    }
    else {
	// Point is not on line
	return RN_NULL_CLASS_ID;
    }
}



RNClassID R3Intersects(const R3Point& point, const R3Ray& ray, RNScalar *hit_t)
{
    // Check if ray contains point
    if (R3Contains(ray, point)) {
	// Point lies on ray
	if (hit_t) *hit_t = ray.T(point);
	return R3_POINT_CLASS_ID;
    }
    else {
	// Point is not on ray
	return RN_NULL_CLASS_ID;
    }
}



RNClassID R3Intersects(const R3Point& point, const R3Span& span, RNScalar *hit_t)
{
    // Check if span contains point
    if (R3Contains(span, point)) {
	// Point lies on span
	if (hit_t) *hit_t = span.T(point);
	return R3_POINT_CLASS_ID;
    }
    else {
	// Point is not on span
	return RN_NULL_CLASS_ID;
    }
}



RNClassID R3Intersects(const R3Point& point, const R3Plane& plane)
{
    // Check whether plane contains point
    if (R3Contains(plane, point)) {
	// Point lies in plane
	return R3_POINT_CLASS_ID;
    }
    else {
	// Point does not lie in plane
	return RN_NULL_CLASS_ID;
    }
}



RNClassID R3Intersects(const R3Point& point, const R3Triangle& triangle)
{
    // Check if triangle contains point
    if (R3Contains(triangle, point)) {
	// Point lies on triangle
	return R3_POINT_CLASS_ID;
    }
    else {
	// Point is not on triangle
	return RN_NULL_CLASS_ID;
    }
}



RNClassID R3Intersects(const R3Point& point, const R3Circle& circle)
{
    // Check if circle contains point
    if (R3Contains(circle, point)) {
	// Point lies on circle
	return R3_POINT_CLASS_ID;
    }
    else {
	// Point is not on circle
	return RN_NULL_CLASS_ID;
    }
}



RNClassID R3Intersects(const R3Point& point, const R3Ellipse& ellipse)
{
    // Check if ellipse contains point
    if (R3Contains(ellipse, point)) {
	// Point lies on ellipse
	return R3_POINT_CLASS_ID;
    }
    else {
	// Point is not on ellipse
	return RN_NULL_CLASS_ID;
    }
}



RNClassID R3Intersects(const R3Point& point, const R3Rectangle& rectangle)
{
    // Check if rectangle contains point
    if (R3Contains(rectangle, point)) {
	// Point lies on rectangle
	return R3_POINT_CLASS_ID;
    }
    else {
	// Point is not on rectangle
	return RN_NULL_CLASS_ID;
    }
}



RNClassID R3Intersects(const R3Point& point, const R3Halfspace& halfspace)
{
    // Check whether halfspace contains point
    if (R3Contains(halfspace, point)) {
	// Point lies in halfspace
	return R3_POINT_CLASS_ID;
    }
    else {
	// Point does not lie in halfspace
	return RN_NULL_CLASS_ID;
    }
}



RNClassID R3Intersects(const R3Point& point, const R3Box& box)
{
    // Check whether box contains point
    if (R3Contains(box, point)) {
	// Point lies in box
	return R3_POINT_CLASS_ID;
    }
    else {
	// Point does not lie in box
	return RN_NULL_CLASS_ID;
    }
}



RNClassID R3Intersects(const R3Point& point, const R3OrientedBox& box)
{
    // Check whether box contains point
    if (R3Contains(box, point)) {
	// Point lies in box
	return R3_POINT_CLASS_ID;
    }
    else {
	// Point does not lie in box
	return RN_NULL_CLASS_ID;
    }
}



RNClassID R3Intersects(const R3Point& point, const R3Sphere& sphere)
{
    // Check whether sphere contains point
    if (R3Contains(sphere, point)) {
	// Point lies in sphere
	return R3_POINT_CLASS_ID;
    }
    else {
	// Point does not lie in sphere
	return RN_NULL_CLASS_ID;
    }
}



RNClassID R3Intersects(const R3Point& point, const R3Cylinder& cylinder)
{
    // Check whether cylinder contains point
    if (R3Contains(cylinder, point)) {
	// Point lies in cylinder
	return R3_POINT_CLASS_ID;
    }
    else {
	// Point does not lie in cylinder
	return RN_NULL_CLASS_ID;
    }
}



RNClassID R3Intersects(const R3Point& point, const R3Cone& cone)
{
    // Check whether cone contains point
    if (R3Contains(cone, point)) {
	// Point lies in cone
	return R3_POINT_CLASS_ID;
    }
    else {
	// Point does not lie in cone
	return RN_NULL_CLASS_ID;
    }
}



RNClassID R3Intersects(const R3Point& point, const R3Shape& shape)
{
    // Return whether shape intersects point
    return shape.Intersects(point);
}



RNClassID R3Intersects(const R3Line& line1, const R3Line& line2, 
    R3Point *hit_point)
{
    // There's got to be a better way ???

    // Get vectors in more convenient form
    const R3Vector v1 = line1.Vector();
    const R3Vector v2 = line2.Vector();

    // Compute useful intermediate values
    const RNScalar v1v1 = 1.0;  // v1.Dot(v1);
    const RNScalar v2v2 = 1.0;  // v2.Dot(v2);
    RNScalar v1v2 = v1.Dot(v2);
    RNScalar denom = v1v2*v1v2 - v1v1*v2v2;

    // Check if lines are parallel
    if (RNIsZero(denom)) {
	// Check if lines are coincident
	if (R3Contains(line1, line2.Point())) {
	    // Lines are coincident
	    return R3_LINE_CLASS_ID;
	}
	else {
	    // Lines are parallel and do not intersect
	    return RN_NULL_CLASS_ID;
	}
    }
    else {
	// Find point on this line closest to other line
	const R3Vector p1 = line1.Point().Vector();
	const R3Vector p2 = line2.Point().Vector();
	RNScalar p1v1 = v1.Dot(p1);
	RNScalar p2v2 = v2.Dot(p2);
	RNScalar p1v2 = v2.Dot(p1);
	RNScalar p2v1 = v1.Dot(p2);
	RNScalar t = (v1v2*p2v2 + v2v2*p1v1 - v1v2*p1v2 - v2v2*p2v1) / denom;
	R3Point point = line1.Point() + line1.Vector() * t;

	// Check if point is on other line
	if (R3Contains(line2, point)) {
	    // Lines intersect at point
	    if (hit_point) *hit_point = point;
	    return R3_POINT_CLASS_ID;
	}
	else {
	    // Lines do not intersect
	    return RN_NULL_CLASS_ID;
	}
    }
}



RNClassID R3Intersects(const R3Line& line, const R3Ray& ray, 
    R3Point *hit_point, RNScalar *hit_t)
{
    // There's got to be a better way ???

    // Get vectors in more convenient form
    const R3Vector& v1 = ray.Vector();
    const R3Vector& v2 = line.Vector();

    // Compute useful intermediate values
    const RNScalar v1v1 = 1.0;  // v1.Dot(v1);
    const RNScalar v2v2 = 1.0;  // v2.Dot(v2);
    RNScalar v1v2 = v1.Dot(v2);
    RNScalar denom = v1v2*v1v2 - v1v1*v2v2;

    // Check if line and ray are parallel
    if (RNIsZero(denom)) {
	// Check if line and ray are coincident
	if (R3Contains(line, ray.Start())) {
	    // Line and ray are coincident
	    return R3_RAY_CLASS_ID;
	}
	else {
	    // Lines are parallel and do not intersect
	    return RN_NULL_CLASS_ID;
	}
    }
    else {
	// Find point on this line closest to other line
	const R3Vector p1 = ray.Start().Vector();
	const R3Vector p2 = line.Point().Vector();
	RNScalar p1v1 = v1.Dot(p1);
	RNScalar p2v2 = v2.Dot(p2);
	RNScalar p1v2 = v2.Dot(p1);
	RNScalar p2v1 = v1.Dot(p2);
	RNScalar s = (v1v2*p2v2 + v2v2*p1v1 - v1v2*p1v2 - v2v2*p2v1) / denom;
	R3Point point = ray.Start() + ray.Vector() * s;

	// Check if point is on other line
	if (R3Contains(line, point)) {
	    // Lines intersect at point
	    if (hit_point) *hit_point = point;
	    if (hit_t) *hit_t = s;
	    return R3_POINT_CLASS_ID;
	}
	else {
	    // Lines do not intersect
	    return RN_NULL_CLASS_ID;
	}
    }
}



RNClassID R3Intersects(const R3Line& line, const R3Span& span, 
    R3Point *hit_point, RNScalar *hit_t)
{
    // Intersect line with span's ray
    RNScalar s;
    switch (R3Intersects(line, span.Ray(), hit_point, &s)) {
    case RN_NULL_CLASS_ID:
	// Span and line do not intersect
	return RN_NULL_CLASS_ID;

    case R3_RAY_CLASS_ID:
	// Line contains span
	return R3_SPAN_CLASS_ID;

    case R3_POINT_CLASS_ID: 
	// Check if intersection point is on span
	if (RNIsGreater(s, span.Length())) {
	    // Span's ray intersects line beyond span endpoint
	    return RN_NULL_CLASS_ID;
	}
	else {
	    // Span and line intersect at point
	    if (hit_t) *hit_t = s;
	    return R3_POINT_CLASS_ID;
	}
	break;

    default:
	RNAbort("Illegal intersection");
	break;
    }

    // Should never get here
    RNAbort("Should never get here");
    return RN_NULL_CLASS_ID;
}



RNClassID R3Intersects(const R3Line& line, const R3Plane& plane, 
    R3Point *hit_point)
{
    // Compute dot product of line vector and plane normal
    RNScalar denom = plane.Normal().Dot(line.Vector());

    // Check if line and plane are parallel
    if (RNIsZero(denom)) {
	// Check if line lies in plane	
	if (R3Contains(plane, line.Point())) {
	    // Line lies in plane
	    return R3_LINE_CLASS_ID;
	}
	else {
	    // Line is parallel to plane, but not in it
	    return RN_NULL_CLASS_ID;
	}
    }
    else {
	// Line and plane are not parallel - they must intersect at a point
	if (hit_point) {
	    RNScalar t = -(R3SignedDistance(plane, line.Point())) / denom;
	    *hit_point = line.Point() + line.Vector() * t;
	}
	return R3_POINT_CLASS_ID;
    }
}



RNClassID R3Intersects(const R3Line& line, const R3Triangle& triangle, 
    R3Point *hit_point)
{
    // Find intersection of line and plane
    R3Point p;
    switch (R3Intersects(line, triangle.Plane(), &p)) {
    case RN_NULL_CLASS_ID:
	// Line does not intersect triangle plane
	return RN_NULL_CLASS_ID;
	
    case R3_POINT_CLASS_ID:
	// Line intersects triangle plane in a single point
	// Check whether triangle contains intersection point
	if (!R3Contains(triangle, p)) return RN_NULL_CLASS_ID;
	if (hit_point) *hit_point = p;
	return R3_POINT_CLASS_ID;

    case R3_LINE_CLASS_ID: 
	// Line lies in triangle plane 
	// Intersect line and triangle in 2D
        // RNAbort("Case not implemented");
	return R3_SPAN_CLASS_ID;

    default:
	RNAbort("Invalid intersection result");
	break;
    }

    // Should never get here
    RNAbort("Should never get here");
    return RN_NULL_CLASS_ID;
}



RNClassID R3Intersects(const R3Line& line, const R3Halfspace& halfspace, 
    R3Ray *result)
{
    // Compute dot product of line vector and halfspace plane normal
    RNScalar denom = halfspace.Plane().Normal().Dot(line.Vector());

    // Check if line and halfspace plane are parallel
    if (RNIsZero(denom)) {
	// Check if line lies in halfspace
	if (R3Contains(halfspace, line.Point())) {
	    // Line lies in halfspace
	    return R3_LINE_CLASS_ID;
	}
	else {
	    // Line is parallel to halfspace plane, but not in halfspace
	    return RN_NULL_CLASS_ID;
	}
    }
    else {
	// Line and halfspace plane are not parallel - intersection must be a ray
	if (result) {
	    RNScalar t = -(R3SignedDistance(halfspace.Plane(), line.Point())) / denom;
	    R3Point point = line.Point() + line.Vector() * t;
	    R3Vector vector = (RNIsPositive(denom)) ? line.Vector() : -(line.Vector());
	    *result = R3Ray(point, vector);
	}
	return R3_RAY_CLASS_ID;
    }
}



RNClassID R3Intersects(const R3Line& line, const R3Box& box, 
    R3Point *hit_point)
{
    RNClassID id = RN_NULL_CLASS_ID;
    R3Point point1;
    RNScalar tval1;

    // Check if box contains point on line
    if (R3Contains(box, line.Point())) {
	if (hit_point) *hit_point = line.Point();
	return R3_SPAN_CLASS_ID;
    }
    else {
	// Find parametric distance to front plane in each dimension
	for (RNDimension dim = RN_X; dim <= RN_Z; dim++) {
	    if (RNIsPositive(line.Vector()[dim])) {
		RNScalar delta = box.Min()[dim] - line.Point()[dim];
		tval1 = delta / line.Vector()[dim];
	    }
	    else if (RNIsNegative(line.Vector()[dim])) {
		RNScalar delta = box.Max()[dim] - line.Point()[dim];
		tval1 = delta / line.Vector()[dim];
	    }
	    else {
		continue;
	    }
	    
	    // Compute intersection point on plane
	    RNDimension dim1 = (dim + 1) % 3;
	    RNDimension dim2 = (dim + 2) % 3;
	    point1[dim] = (line.Vector()[dim] > 0.0) ? box.Min()[dim] : box.Max()[dim];
	    assert(RNIsEqual(point1[dim], line.Point()[dim] + line.Vector()[dim] * tval1));
	    point1[dim1] = line.Point()[dim1] + line.Vector()[dim1] * tval1;
	    point1[dim2] = line.Point()[dim2] + line.Vector()[dim2] * tval1;
	    
	    // Check if intersection point on plane is inside box
	    if ((RNIsLessOrEqual(point1[dim1], box.Max()[dim1])) &&
		(RNIsGreaterOrEqual(point1[dim1], box.Min()[dim1])) &&
		(RNIsLessOrEqual(point1[dim2], box.Max()[dim2])) &&
		(RNIsGreaterOrEqual(point1[dim2], box.Min()[dim2]))) {
		// Found first intersection point 
		id = R3_SPAN_CLASS_ID;
		break;
	    }
	}

	// Fill in results
	if (id != RN_NULL_CLASS_ID) {
	    if (hit_point) *hit_point = point1;
	}
	
	// Return intersection 
	return id;
    }
}



RNClassID R3Intersects(const R3Line& line, const R3OrientedBox& box, 
    R3Point *hit_point)
{
    // Check if box is empty
    if (box.IsEmpty()) return RN_NULL_CLASS_ID;

    // Transform line into std coordinate system
    R3Line transformed_line = line;
    R3Affine cs_to_std_xform(box.CoordSystem().InverseMatrix());
    transformed_line.Transform(cs_to_std_xform);

    // Intersect in std coordinate system
    R3Box axis_aligned_box(-box.Radius(0), -box.Radius(1), -box.Radius(2), box.Radius(0), box.Radius(1), box.Radius(2));
    int status = R3Intersects(transformed_line, axis_aligned_box, hit_point);
    if (status == RN_NULL_CLASS_ID) return RN_NULL_CLASS_ID;

    // Compute hit point
    if (hit_point) {
      R3Affine std_to_cs_xform(box.CoordSystem().Matrix());
      hit_point->Transform(std_to_cs_xform);
    }

    // Return whether hit
    return status;
}



RNClassID R3Intersects(const R3Line& line, const R3Sphere& sphere, 
    R3Point *hit_point)
{
    // Check if sphere contains line point
    if (R3Contains(sphere, line.Point())) {
	if (hit_point) *hit_point = line.Point();
	return R3_SPAN_CLASS_ID;
    }
    else {
	// Check if line points towards sphere
	R3Vector EO = sphere.Center() - line.Point();
	RNScalar v = EO.Dot(line.Vector());

	// Compute possible intersection -- Graphics Gems I, p. 388    
	RNScalar disc = sphere.Radius() * sphere.Radius() - (EO.Dot(EO) - v*v);
	if (RNIsNegative(disc)) {
	    // No intersection
	    return RN_NULL_CLASS_ID;
	}
	else {
	    // Line intersects sphere (it grazes if disc is zero)
	    if (hit_point) {
		RNScalar d = sqrt(disc);
		RNScalar t = v - d;
		*hit_point = line.Point() + t * line.Vector();
	    }
	    
	    // Return result
	    if (RNIsZero(disc)) return R3_POINT_CLASS_ID;
	    else return R3_SPAN_CLASS_ID;
	}
    }
}



RNClassID R3Intersects(const R3Line& line, const R3Shape& shape)
{
    // Return whether line intersects shape
    return shape.Intersects(line);
}



RNClassID R3Intersects(const R3Ray& ray1, const R3Ray& ray2, 
    R3Point *hit_point, RNScalar *hit_t1, RNScalar *hit_t2)
{
    // Intersect ray1 with ray2's line
    R3Point point;
    RNScalar s1 = 0, s2 = 0;
    switch (R3Intersects(ray1, ray2.Line(), &point, &s1)) {
    case RN_NULL_CLASS_ID:
	// Rays do not intersect
	return RN_NULL_CLASS_ID;

    case R3_RAY_CLASS_ID:
	// Ray2's line contains ray1 -- check directions
	if (R3Contains(ray1.Vector(), ray2.Vector())) {
	    // Rays share same line in same direction
	    return R3_RAY_CLASS_ID;
	}
	else {
	    // Rays share same line, but in different directions
	    assert(R3Contains(ray1.Vector(), -(ray2.Vector())));

	    // Check intersection point
	    s2 = ray2.T(ray1.Start());
	    if (RNIsNegative(s2)) {
		// Rays share line, but intersection is behind endpoint of ray2
		return RN_NULL_CLASS_ID;
	    }
	    else if (RNIsZero(s2)) {
		// Rays share line and endpoint
		if (hit_t1) *hit_t1 = s1;
		if (hit_t2) *hit_t2 = s2;
		if (hit_point) *hit_point = point;
		return R3_POINT_CLASS_ID;
	    }
	    else {
		// Rays share line and intersect in span 
		return R3_SPAN_CLASS_ID;
	    }
	}

    case R3_POINT_CLASS_ID: 
	// Check if intersection point is on ray2
	s2 = ray2.T(point);
	if (RNIsNegative(s1)) {
	    // Line intersection is behind ray1
	    return RN_NULL_CLASS_ID;
	}
	else {
	    // Rays intersect at point
	    if (hit_t1) *hit_t1 = s1;
	    if (hit_t2) *hit_t2 = s2;
	    if (hit_point) *hit_point = point;
	    return R3_POINT_CLASS_ID;
	}
	break;

    default:
	RNAbort("Illegal intersection");
	break;
    }

    // Should never get here
    RNAbort("Should never get here");
    return RN_NULL_CLASS_ID;
}



RNClassID R3Intersects(const R3Ray& ray, const R3Span& span, 
    R3Point *hit_point, RNScalar *hit_tray, RNScalar *hit_tspan)
{
    // Intersect ray with span's ray
    RNScalar s1 = 0, s2 = 0;
    switch (R3Intersects(ray, span.Ray(), hit_point, &s1, &s2)) {
    case RN_NULL_CLASS_ID:
	// Span and ray do not intersect
	return RN_NULL_CLASS_ID;

    case R3_SPAN_CLASS_ID:
	// Ray and span point at each other on same line
	return R3_SPAN_CLASS_ID;

    case R3_RAY_CLASS_ID:
	// Ray and span point in same direction on same line
	if (RNIsGreater(s2, span.Length())) {
	    // Intersection is beyond span endpoint 
	    return RN_NULL_CLASS_ID;
	}
	else if (RNIsEqual(s2, span.Length())) {
	    // Span endpoint is ray start point
	    if (hit_tray) *hit_tray = s1;
	    if (hit_tspan) *hit_tspan = s2;
	    return R3_POINT_CLASS_ID;
	}
	else {
	    // Span and ray have non-null intersection in same line
	    return R3_SPAN_CLASS_ID;
	}
	break;

    case R3_POINT_CLASS_ID: 
	// Check if intersection point is on span
	if (RNIsGreater(s2, span.Length())) {
	    // Intersection is beyond span endpoint 
	    return RN_NULL_CLASS_ID;
	}
	else {
	    // Span and ray intersect at point
	    if (hit_tray) *hit_tray = s1;
	    if (hit_tspan) *hit_tspan = s2;
	    return R3_POINT_CLASS_ID;
	}
	break;

    default:
	RNAbort("Illegal intersection");
	break;
    }

    // Should never get here
    RNAbort("Should never get here");
    return RN_NULL_CLASS_ID;
}



RNClassID R3Intersects(const R3Ray& ray, const R3Plane& plane, 
    R3Point *hit_point, RNScalar *hit_t)
{
    // Compute dot product of ray vector and plane normal
    RNScalar denom = plane.Normal().Dot(ray.Vector());

    // Check if ray and plane are parallel
    if (RNIsZero(denom)) {
	// Check if ray lies in plane	
	if (R3Contains(plane, ray.Start())) {
	    // Ray lies in plane
	    return R3_RAY_CLASS_ID;
	}
	else {
	    // Ray is parallel to plane, but not in it
	    return RN_NULL_CLASS_ID;
	}
    }
    else {
	// Ray and plane are not parallel - compute parametric value on ray
	RNScalar s = -(R3SignedDistance(plane, ray.Start())) / denom;
	if (RNIsNegative(s)) {
	    // Ray is directed away from plane
	    return RN_NULL_CLASS_ID;
	}
	else {
	    // Ray is not parallel to plane and is pointing at plane - it must intersect
	    if (hit_point) *hit_point = ray.Start() + ray.Vector() * s;
	    if (hit_t) *hit_t = s;
	    return R3_POINT_CLASS_ID;
	}
    }
}



RNClassID R3Intersects(const R3Ray& ray, const R3Halfspace& halfspace, 
    R3Point *hit_point, RNScalar *hit_t)
{
    // Check if halfspace contains point
    RNBoolean pinside = R3Contains(halfspace, ray.Start());

    // Check whether vectors point in same direction
    RNBoolean vinside = RNIsPositive(halfspace.Normal().Dot(ray.Vector()));

    // Classify ray with respect to halfspace
    if (pinside) {
	// Start point of ray is inside halfspace
	if (hit_point) *hit_point = ray.Start();
	if (hit_t) *hit_t = 0.0;
	if (vinside) return R3_RAY_CLASS_ID;
	else return R3_SPAN_CLASS_ID;
    }
    else {
	if (vinside) return R3Intersects(ray, halfspace.Plane(), hit_point, hit_t);
	else return RN_NULL_CLASS_ID;
    }
}



RNClassID R3Intersects(const R3Ray& ray, const R3Triangle& triangle, 
    R3Point *hit_point, R3Vector *hit_normal, RNScalar *hit_t)
{
    // Find intersection of ray and plane
    R3Point p1;
    RNScalar t1;
    switch (R3Intersects(ray, triangle.Plane(), &p1, &t1)) {
    case RN_NULL_CLASS_ID:
	// Ray does not intersect triangle plane
	return RN_NULL_CLASS_ID;
	
    case R3_POINT_CLASS_ID:
	// Ray intersects triangle plane in a single point
	// Check whether triangle contains intersection point
	if (!R3Contains(triangle, p1)) return RN_NULL_CLASS_ID;
	if (hit_point) *hit_point = p1;
	if (hit_normal) *hit_normal = triangle.Normal();
	if (hit_t) *hit_t = t1;
	return R3_POINT_CLASS_ID;

    case R3_RAY_CLASS_ID: 
	// Ray lies in triangle plane 
	// Intersect ray and triangle in 2D
	// RNAbort("Case not implemented");
	// return R3_SPAN_CLASS_ID;
	return RN_NULL_CLASS_ID;

    default:
	RNAbort("Invalid intersection result");
	break;
    }

    // Should never get here
    RNAbort("Should never get here");
    return RN_NULL_CLASS_ID;
}



RNClassID R3Intersects(const R3Ray& ray, const R3TriangleArray& array,
    R3Point *hit_point, R3Vector *hit_normal, RNScalar *hit_t)
{
    // This could be more efficient ???
    // Sort intersections with front-facing planes 
    //   before do any point-in-triangle operations

    // Check bounding volume for intersection 
    if (!R3Intersects(ray, array.Box())) 
	return RN_NULL_CLASS_ID;

    // Check each triangle for intersection
    RNClassID status = RN_NULL_CLASS_ID;
    RNScalar min_t = FLT_MAX;
    for (int i = 0; i < array.NTriangles(); i++) {
        R3Point point;
        R3Vector normal;
	RNScalar t;
        if (R3Intersects(ray, *(array.Triangle(i)), &point, &normal, &t) == R3_POINT_CLASS_ID) {
	    if (t < min_t) {
	        status = R3_POINT_CLASS_ID;
	        if (hit_point) *hit_point = point;
	        if (hit_normal) *hit_normal = normal;
	        min_t = t;
	    }
	}
    }

    // Update hit t
    if (hit_t) *hit_t = min_t;

    // Return whether hit any triangle
    return status;
}



RNClassID R3Intersects(const R3Ray& ray, const R3Circle& circle, 
    R3Point *hit_point, R3Vector *hit_normal, RNScalar *hit_t)
{
    // Check for empty circle
    if (circle.IsEmpty()) return RN_NULL_CLASS_ID;
    
    // Find intersection of ray and plane
    R3Point p1;
    RNScalar t1;
    switch (R3Intersects(ray, circle.Plane(), &p1, &t1)) {
    case RN_NULL_CLASS_ID:
	// Ray does not intersect circle plane
	return RN_NULL_CLASS_ID;
	
    case R3_POINT_CLASS_ID:
    {
	// Ray intersects circle plane in a single point
	// Check whether circle contains intersection point
	RNScalar radius_squared = circle.Radius() * circle.Radius();
	R3Vector v = circle.Center() - p1;
	RNScalar distance_squared = v.X() * v.X() + v.Y() * v.Y() + v.Z() * v.Z();
	if (RNIsGreater(distance_squared, radius_squared)) return RN_NULL_CLASS_ID;
	if (hit_point) *hit_point = p1;
	if (hit_normal) *hit_normal = circle.Normal();
	if (hit_t) *hit_t = t1;
	return R3_POINT_CLASS_ID;
    }

    case R3_RAY_CLASS_ID: 
	// Ray lies in circle plane 
	// Intersect ray and circle in 2D
        // RNAbort("Case not implemented");
	return R3_SPAN_CLASS_ID;

    default:
	RNAbort("Invalid intersection result");
	break;
    }

    // Should never get here
    RNAbort("Should never get here");
    return RN_NULL_CLASS_ID;
}



RNClassID R3Intersects(const R3Ray& ray, const R3Ellipse& ellipse, 
    R3Point *hit_point, R3Vector *hit_normal, RNScalar *hit_t)
{
    // Check for empty ellipse
    if (ellipse.IsEmpty()) return RN_NULL_CLASS_ID;
    
    // Find intersection of ray and plane
    R3Point p1;
    RNScalar t1;
    switch (R3Intersects(ray, ellipse.Plane(), &p1, &t1)) {
    case RN_NULL_CLASS_ID:
	// Ray does not intersect ellipse plane
	return RN_NULL_CLASS_ID;
	
    case R3_POINT_CLASS_ID:
    {
	// Ray intersects ellipse plane in a single point
	// Check whether ellipse contains intersection point
        if (!R3Contains(ellipse, p1)) return RN_NULL_CLASS_ID;
  	if (hit_point) *hit_point = p1;
	if (hit_normal) *hit_normal = ellipse.Normal();
	if (hit_t) *hit_t = t1;
	return R3_POINT_CLASS_ID;
    }

    case R3_RAY_CLASS_ID: 
	// Ray lies in ellipse plane 
	// Intersect ray and ellipse in 2D
        // RNAbort("Case not implemented");
	return R3_SPAN_CLASS_ID;

    default:
	RNAbort("Invalid intersection result");
	break;
    }

    // Should never get here
    RNAbort("Should never get here");
    return RN_NULL_CLASS_ID;
}



RNClassID R3Intersects(const R3Ray& ray, const R3Rectangle& rectangle, 
    R3Point *hit_point, R3Vector *hit_normal, RNScalar *hit_t)
{
    // Check for empty rectangle
    if (rectangle.IsEmpty()) return RN_NULL_CLASS_ID;
    
    // Find intersection of ray and plane
    R3Point p1;
    RNScalar t1;
    switch (R3Intersects(ray, rectangle.Plane(), &p1, &t1)) {
    case RN_NULL_CLASS_ID:
	// Ray does not intersect rectangle plane
	return RN_NULL_CLASS_ID;
	
    case R3_POINT_CLASS_ID:
    {
	// Ray intersects rectangle plane in a single point
	// Check whether rectangle contains intersection point
        if (!R3Contains(rectangle, p1)) return RN_NULL_CLASS_ID;
  	if (hit_point) *hit_point = p1;
	if (hit_normal) *hit_normal = rectangle.Normal();
	if (hit_t) *hit_t = t1;
	return R3_POINT_CLASS_ID;
    }

    case R3_RAY_CLASS_ID: 
	// Ray lies in rectangle plane 
	// Intersect ray and rectangle in 2D
        // RNAbort("Case not implemented");
	return R3_SPAN_CLASS_ID;

    default:
	RNAbort("Invalid intersection result");
	break;
    }

    // Should never get here
    RNAbort("Should never get here");
    return RN_NULL_CLASS_ID;
}



RNClassID R3Intersects(const R3Ray& ray, const R3Box& box, 
    R3Point *hit_point1, R3Vector *hit_normal1, RNScalar *hit_t1)
{
    // Check box 
    if (box.IsEmpty()) {
        // Box is empty
        return RN_NULL_CLASS_ID;
    }
    else {
        R3Point point1;
        RNScalar tval1;

        // Check if ray start is inside box
        RNBoolean start_inside = R3Contains(box, ray.Start());

	// Find parametric distance to front plane in each dimension
	for (RNDimension dim = RN_X; dim <= RN_Z; dim++) {
	    if (RNIsPositive(ray.Vector()[dim])) {
                RNCoord box_coord = (start_inside) ? box.Max()[dim] : box.Min()[dim];
		RNScalar delta = box_coord - ray.Start()[dim];
		if (delta < 0.0) continue;
		tval1 = delta / ray.Vector()[dim];
	    }
	    else if (RNIsNegative(ray.Vector()[dim])) {
                RNCoord box_coord = (start_inside) ? box.Min()[dim] : box.Max()[dim];
		RNScalar delta = box_coord - ray.Start()[dim];
		if (delta > 0.0) continue;
		tval1 = delta / ray.Vector()[dim];
	    }
	    else {
		continue;
	    }
	    
	    // Compute intersection point on plane
	    RNDimension dim1 = (dim + 1) % 3;
	    RNDimension dim2 = (dim + 2) % 3;
	    point1[dim] = ray.Start()[dim] + ray.Vector()[dim] * tval1;
	    point1[dim1] = ray.Start()[dim1] + ray.Vector()[dim1] * tval1;
	    point1[dim2] = ray.Start()[dim2] + ray.Vector()[dim2] * tval1;
	    
	    // Check if intersection point on plane is inside box
	    if ((RNIsLessOrEqual(point1[dim1], box.Max()[dim1])) &&
		(RNIsGreaterOrEqual(point1[dim1], box.Min()[dim1])) &&
		(RNIsLessOrEqual(point1[dim2], box.Max()[dim2])) &&
		(RNIsGreaterOrEqual(point1[dim2], box.Min()[dim2]))) {
		// Found first intersection point 
  	        if (hit_t1) *hit_t1 = tval1;
	        if (hit_point1) *hit_point1 = point1;
	        if (hit_normal1) {
	    	    if (RNIsNegative(ray.Vector()[dim])) *hit_normal1 = R3xyz_triad.Axis(dim);
		    else *hit_normal1 = -(R3xyz_triad.Axis(dim));
	        }
		return R3_SPAN_CLASS_ID;
	    }
	}

	// Return no intersection 
	return RN_NULL_CLASS_ID;
    }
}



RNClassID R3Intersects(const R3Ray& ray, const R3OrientedBox& box, 
    R3Point *hit_point1, R3Vector *hit_normal1, RNScalar *hit_t1)
{
    // Check if box is empty
    if (box.IsEmpty()) return RN_NULL_CLASS_ID;

    // Transform ray into std coordinate system
    R3Ray transformed_ray = ray;
    R3Affine cs_to_std_xform(box.CoordSystem().InverseMatrix());
    transformed_ray.Transform(cs_to_std_xform);

    // Intersect in std coordinate system
    R3Box axis_aligned_box(-box.Radius(0), -box.Radius(1), -box.Radius(2), box.Radius(0), box.Radius(1), box.Radius(2));
    int status = R3Intersects(transformed_ray, axis_aligned_box, hit_point1, hit_normal1, hit_t1);
    if (status == RN_NULL_CLASS_ID) return RN_NULL_CLASS_ID;

    // Compute hit point
    if (hit_point1 || hit_normal1) {
      R3Affine std_to_cs_xform(box.CoordSystem().Matrix());
      if (hit_point1) hit_point1->Transform(std_to_cs_xform);
      if (hit_normal1) hit_normal1->Transform(std_to_cs_xform);
    }

    // Return whether hit
    return status;
}



RNClassID R3Intersects(const R3Ray& ray, const R3Sphere& sphere, 
    R3Point *hit_point1, R3Vector *hit_normal1, RNScalar *hit_t1)
{
    // Check if sphere contains ray start point
    RNBoolean start_inside = R3Contains(sphere, ray.Start());

    // Check if ray points towards sphere
    R3Vector EO = sphere.Center() - ray.Start();
    RNScalar v = EO.Dot(ray.Vector());
    if (!start_inside && RNIsNegativeOrZero(v)) {
        // Ray points away from sphere
        return RN_NULL_CLASS_ID;
    }
    else {
        // Compute possible intersection -- Graphics Gems I, p. 388    
        RNScalar disc = sphere.Radius() * sphere.Radius() - (EO.Dot(EO) - v*v);
        if (RNIsNegative(disc)) {
	    // No intersection
	    return RN_NULL_CLASS_ID;
        }
        else {
            // Ray intersects sphere (it grazes if disc is zero)
            // Compute first intersection
            if (hit_t1 || hit_point1 || hit_normal1) {
                RNScalar d = sqrt(disc);
                RNScalar t = (start_inside) ? v + d : v - d;
                R3Point p = ray.Start() + t * ray.Vector();
                if (hit_t1) *hit_t1 = t;
                if (hit_point1) *hit_point1 = p;
                if (hit_normal1) *hit_normal1 = (p - sphere.Centroid()) / sphere.Radius();
            }

            // Compute second intersection
            // if (hit_t2 || hit_point2 || hit_normal2) {
            //     if (!start_inside) {
            //         t = v + d;
            //         if (hit_point2) *hit_point2 = ray.Start() + t * ray.Vector();
            //         if (hit_t2) *hit_t2 = t;
            //     }
            // }
            
            // Return result
            if (RNIsZero(disc)) return R3_POINT_CLASS_ID;
            else return R3_SPAN_CLASS_ID;
        }
    }
}



RNClassID R3Intersects(const R3Ray& ray, const R3Cylinder& cylinder, 
    R3Point *hit_point1, R3Vector *hit_normal1, RNScalar *hit_t1)
{
    // See Graphics Gems IV, page 356

    // Compute parametric values of intersection with infinite cylinder
    RNScalar cyl_t1, cyl_t2;
    R3Vector R = ray.Vector();
    R3Vector A = cylinder.Axis().Vector();
    R3Vector D = R % A;
    RNScalar a = D.Length();
    if (RNIsZero(a)) {
	// Ray and axis are parallel
	RNScalar d = R3Distance(ray.Start(), cylinder.Axis().Line());
	if (RNIsLessOrEqual(d, cylinder.Radius())) {
	    // Ray lies entirely inside infinite cylinder
	    cyl_t1 = 0.0;
	    cyl_t2 = RN_INFINITY;
	}
	else {
	    // Ray lies entirely outside infinite cylinder
	    return RN_NULL_CLASS_ID;
	}
    }
    else {
	// Ray and axis are not parallel
	// Compute distance from axis to closest point on ray's line
	D /= a;
	R3Vector RC = ray.Start() - cylinder.Axis().Start();
	RNScalar d = RC.Dot(D);
	if (d < 0.0) d = -d;
	if (RNIsGreater(d, cylinder.Radius())) {
	    // Ray does not intersect infinite cylinder
	    return RN_NULL_CLASS_ID;
	}
	else {
	    // Find parametric value (cyl_t) of point on ray closest to axis
	    assert(RNIsNotZero(a));
	    RNScalar t = (RC % A).Dot(D) / -a;
	    
	    // Find parametric value delta (s) to intersection points
	    RNScalar s;
	    R3Vector O = D % A;
	    O.Normalize();
	    RNScalar b = cylinder.Radius() * cylinder.Radius() - d*d;
	    if (RNIsPositive(b)) {
		// Crossing span intersection
		RNScalar e = R.Dot(O);
		assert(RNIsNotZero(e));
		s = sqrt(b) / e;
		if (s < 0.0) s = -s;
	    }
	    else {
		// Grazing point intersection
		s = 0.0;
	    }
	    
	    // Find parametric value (cyl_t2) of second intersection
	    cyl_t2 = t + s;
	    if (RNIsNegative(cyl_t2)) {
		// Intersection point is behind ray
		return RN_NULL_CLASS_ID;
	    }
	    else {
		// Find parametric value (cyl_t1) of first intersection
		cyl_t1 = t - s;
		
		// Check if ray start is in infinite cylinder
		if (cyl_t1 < 0.0) cyl_t1 = 0.0;
	    }
	}
    }
	
    // Compute parametric values of intersection with caps
    RNScalar cap_t1 = 0, cap_t2 = 0;
    RNScalar dot = ray.Vector().Dot(cylinder.Axis().Vector());
    RNScalar d_top = R3SignedDistance(cylinder.Top().Plane(), ray.Start());
    RNScalar d_base = R3SignedDistance(cylinder.Base().Plane(), ray.Start());
    if (RNIsPositive(d_top)) {
	// Ray starts above top cap
	if (RNIsPositiveOrZero(dot)) {
	    // Ray points upward or orthogonally from above cylinder top
	    return RN_NULL_CLASS_ID;
	}
	else {
	    // Ray points downward from above cylinder top
	    if ((R3Intersects(ray, cylinder.Top().Plane(), NULL, &cap_t1)) &&
		(RNIsGreater(cap_t1, cyl_t2))) {
		// Ray misses cylinder
		return RN_NULL_CLASS_ID;
	    }
	    else if ((R3Intersects(ray, -(cylinder.Base().Plane()), NULL, &cap_t2)) &&
		     (RNIsLess(cap_t2, cyl_t1))) {
		// Ray misses cylinder
		return RN_NULL_CLASS_ID;
	    }
	    else if (RNIsGreater(cap_t1, cyl_t1)) {
		// Ray hits top cap 
		assert(cap_t1 >= 0.0);
		if (hit_t1) *hit_t1 = cap_t1;
		if (hit_point1) *hit_point1 = ray.Point(cap_t1);
		if (hit_normal1) *hit_normal1 = cylinder.Top().Normal();
		return R3_POINT_CLASS_ID;
	    }
	    else {
		// Ray hits side of cylinder
	    }
	}
    }
    else if (RNIsPositive(d_base)) {
	// Ray starts below base cap
	if (RNIsNegativeOrZero(dot)) {
	    // Ray points downward or orthogonally
	    return RN_NULL_CLASS_ID;
	}
	else {
	    // Ray points upward from below cylinder base
	    if ((R3Intersects(ray, cylinder.Base().Plane(), NULL, &cap_t1)) &&
		(RNIsGreater(cap_t1, cyl_t2))) {
		// Ray misses cylinder
		return RN_NULL_CLASS_ID;
	    }
	    else if ((R3Intersects(ray, -(cylinder.Top().Plane()), NULL, &cap_t2)) &&
		     (RNIsLess(cap_t2, cyl_t1))) {
		// Ray misses cylinder
		return RN_NULL_CLASS_ID;
	    }
	    else if (RNIsGreater(cap_t1, cyl_t1)) {
		// Ray hits base cap 
		assert(cap_t1 >= 0.0);
		if (hit_t1) *hit_t1 = cap_t1;
		if (hit_point1) *hit_point1 = ray.Point(cap_t1);
		if (hit_normal1) *hit_normal1 = cylinder.Base().Normal();
		return R3_POINT_CLASS_ID;
	    }
	    else {
		// Ray hits side of cylinder
	    }
	}
    }
    else {
	// Ray starts on or between caps
	if (RNIsPositive(dot)) {
	    if (R3Intersects(ray, -(cylinder.Top().Plane()), NULL, &cap_t1)) {
		if (RNIsLess(cap_t1, cyl_t1)) {
		    // Intersection with infinite cylinder is beyond cap
		    return RN_NULL_CLASS_ID;
		}
	    }
	}
	else if (RNIsNegative(dot)) {
	    if (R3Intersects(ray, -(cylinder.Base().Plane()), NULL, &cap_t1)) {
		if (RNIsLess(cap_t1, cyl_t1)) {
		    // Intersection with infinite cylinder is beyond cap
		    return RN_NULL_CLASS_ID;
		}
	    }
	}
    }

    // Ray hit side of cylinder
    assert(cyl_t1 >= 0.0);
    if (hit_t1) *hit_t1 = cyl_t1;
    if (hit_point1) *hit_point1 = ray.Point(cyl_t1);
    if (hit_normal1) {
	R3Point p, *H;
	if (hit_point1) H = hit_point1;
	else { p = ray.Point(cyl_t1); H = &p; }
        R3Vector HB = *H - cylinder.Axis().Start();
	*hit_normal1 = (HB - HB.Dot(A) * A) / cylinder.Radius();
    }
    return R3_POINT_CLASS_ID;
}



RNClassID R3Intersects(const R3Ray& ray, const R3Cone& cone, 
    R3Point *hit_point1, R3Vector *hit_normal1, RNScalar *hit_t1)
{
    // Transform ray into canonical coordinate system with
    // axis end point at origin and axis going along negative z axis
    R3Ray ray1(ray);
    R3Affine xform(R4identity_matrix);
    xform.Translate(- (cone.Axis().End().Vector()) );
    xform.Rotate(cone.Axis().Vector(), R3negz_vector);
    ray1.Transform(xform);

    // Compute coefficients for quadratic formula
    RNScalar r = cone.Radius();
    const R3Point& p = ray1.Start();
    const R3Vector& v = ray1.Vector();
    RNScalar A = v.X()*v.X() + v.Y()*v.Y() - r*r*v.Z()*v.Z();
    RNScalar B = 2.0*p.X()*v.X() + 2.0*p.Y()*v.Y() - 2.0*r*r*p.Z()*v.Z();
    RNScalar C = p.X()*p.X() + p.Y()*p.Y() - r*r*p.Z()*p.Z();

    // Solve for t1 and t2 on sides
    RNScalar ray_t;
    if (RNIsZero(A)) {
        if (RNIsZero(B)) {
	    // No intersection
	    return RN_NULL_CLASS_ID;
	}
	else {
	    // ???
	    ray_t = -C/B;
	}
    }
    else {
        RNScalar sqr = B*B - 4.0*A*C;
	if (RNIsNegative(sqr)) {
	    // No intersection
	    return RN_NULL_CLASS_ID;
	}
	else if (RNIsPositive(sqr)) {
	    // Intersection at two points
	    RNScalar s = sqrt(sqr);
	    RNScalar ray_t1 = 0.5 * (-B - s) / A;
	    RNScalar ray_t2 = 0.5 * (-B + s) / A;
	    ray_t = (ray_t1 < ray_t2) ? ray_t1 : ray_t2;
	}
	else {
	    // Intersection at one point
	    ray_t = -0.5 * B / A;
	}
    }

    // Check whether intersection point is on right side of ray
    if (RNIsNegative(ray_t)) return RN_NULL_CLASS_ID;

    // Check whether intersection point is on right side of double cone or beyond cap
    R3Point p1 = ray1.Point(ray_t);
    if (RNIsNegative(p1.Z())) return RN_NULL_CLASS_ID;
    if (RNIsGreater(p1.Z(), cone.Axis().Length())) return RN_NULL_CLASS_ID;

    // Check for intersection with base
    if (RNIsGreater(ray1.Start().Z(), cone.Axis().Length())) {
        // ???
    }

    // This is not complete ???

    // Transform intersection back to world coordinates
    R3Point intersection = ray.Point(ray_t);
    intersection.InverseTransform(xform);

    // Ray hit cone
    if (hit_t1) *hit_t1 = ray_t;
    if (hit_point1) { *hit_point1 = intersection; }
    if (hit_normal1) { *hit_normal1 = R3zero_vector; }
    return R3_POINT_CLASS_ID;
}



RNClassID R3Intersects(const R3Ray& ray, const R3Shape& shape, 
    R3Point *hit_point, R3Vector *hit_normal, RNScalar *hit_t)
{
    // Return whether ray intersects shape
    return shape.Intersects(ray, hit_point, hit_normal, hit_t);
}



RNClassID R3Intersects(const R3Span& span1, const R3Span& span2, 
    R3Point *hit_point, RNScalar *hit_t1, RNScalar *hit_t2)
{
    // Intersect span1 with span2's ray
    RNScalar s1, s2;
    switch (R3Intersects(span1, span2.Ray(), hit_point, &s1, &s2)) {
    case RN_NULL_CLASS_ID:
	// Spans do not intersect
	return RN_NULL_CLASS_ID;

    case R3_SPAN_CLASS_ID:
	// Span1 lies on span2's ray -- check if they overlap ???
	if ((R3Contains(span2, span1.Start())) || (R3Contains(span2, span1.End()))) {
	    // Spans are on same line, and overlap
	    return R3_SPAN_CLASS_ID;
	}
	else {
	    // Spans are on same line, but do not overlap
	    return RN_NULL_CLASS_ID;
	}
	break;

    case R3_POINT_CLASS_ID: 
	// Check if intersection point is on span
	if (RNIsGreater(s1, span1.Length())) {
	    // Span's ray intersects plane beyond span endpoint
	    return RN_NULL_CLASS_ID;
	}
	else {
	    // Span and plane intersect at point
	    if (hit_t1) *hit_t1 = s1;
	    if (hit_t2) *hit_t2 = s2;
	    return R3_POINT_CLASS_ID;
	}
	break;

    default:
	RNAbort("Illegal intersection");
	break;
    }

    // Should never get here
    RNAbort("Should never get here");
    return RN_NULL_CLASS_ID;
}



RNClassID R3Intersects(const R3Span& span, const R3Plane& plane, 
    R3Point *hit_point, RNScalar *hit_t)
{
#if FALSE
    // Intersect plane with span's ray
    // Doesn't get backfacing planes ???
    RNScalar s;
    switch (R3Intersects(plane, span.Ray(), hit_point, &s)) {
    case RN_NULL_CLASS_ID:
	// Span and plane do not intersect
	return RN_NULL_CLASS_ID;

    case R3_RAY_CLASS_ID:
	// Plane contains span
	return R3_SPAN_CLASS_ID;

    case R3_POINT_CLASS_ID: 
	// Check if intersection point is on span
	if (RNIsGreater(s, span.Length())) {
	    // Span's ray intersects plane beyond span endpoint
	    return RN_NULL_CLASS_ID;
	}
	else {
	    // Span and plane intersect at point
	    if (hit_t) *hit_t = s;
	    return R3_POINT_CLASS_ID;
	}
	break;

    default:
	RNAbort("Illegal intersection");
	break;
    }

    // Should never get here
    RNAbort("Should never get here");
    return RN_NULL_CLASS_ID;
#else
    // Compute signed distances from plane to span endpoints
    RNScalar d1 = R3SignedDistance(plane, span.Start());
    RNScalar d2 = R3SignedDistance(plane, span.End());

    // Compute intersection relationship based on endpoint distances
    if (RNIsPositive(d1)) {
        if (RNIsPositive(d2)) {
	    // Span lies above plane
	    return RN_NULL_CLASS_ID;
	}
	else if (d2 >= 0.0) {
	    // Span abuts plane at end point
	    if (hit_point) *hit_point = span.End();
	    if (hit_t) *hit_t = span.Length();
	    return R3_POINT_CLASS_ID;
	}
	else {
	    // Span crosses plane
	    RNScalar t = span.Length() * d1 / (d1 - d2);
	    if (hit_point) *hit_point = span.Point(t);
	    if (hit_t) *hit_t = t;
	    return R3_POINT_CLASS_ID;
	}
    }
    else if (RNIsNegative(d1)) {
        if (RNIsNegative(d2)) {
	    // Span lies below plane
	    return RN_NULL_CLASS_ID;
	}
	else if (d2 <= 0.0) {
	    // Span abuts plane at end point
	    if (hit_point) *hit_point = span.End();
	    if (hit_t) *hit_t = span.Length();
	    return R3_POINT_CLASS_ID;
	}
	else {
	    // Span crosses plane
	    RNScalar t = span.Length() * d1 / (d1 - d2);
	    if (hit_point) *hit_point = span.Point(t);
	    if (hit_t) *hit_t = t;
	    return R3_POINT_CLASS_ID;
	}
    }
    else {
        if (RNIsPositive(d2)) {	  
	    if (d1 >= 0.0) {
		// Span abuts plane at start point
	        if (hit_point) *hit_point = span.Start();
		if (hit_t) *hit_t = 0.0;
		return R3_POINT_CLASS_ID;
	    }
	    else {
	        // Span crosses plane
	        RNScalar t = span.Length() * d1 / (d1 - d2);
		if (hit_point) *hit_point = span.Point(t);
		if (hit_t) *hit_t = t;
		return R3_POINT_CLASS_ID;
	    }
	}
	else if (RNIsNegative(d2)) {	  
	    if (d1 <= 0.0) {
		// Span abuts plane at start point
	        if (hit_point) *hit_point = span.Start();
		if (hit_t) *hit_t = 0.0;
		return R3_POINT_CLASS_ID;
	    }
	    else {
	        // Span crosses plane
	        RNScalar t = span.Length() * d1 / (d1 - d2);
		if (hit_point) *hit_point = span.Point(t);
		if (hit_t) *hit_t = t;
		return R3_POINT_CLASS_ID;
	    }
	}
	else {
	    // Span lies on plane
	    return R3_SPAN_CLASS_ID;
	}
    }
#endif
}



RNClassID R3Intersects(const R3Span& span, const R3Triangle& triangle, 
    R3Point *hit_point1, RNScalar *hit_t1)
{
    // Find intersection of span and plane
    R3Point p;
    RNScalar t;
    switch (R3Intersects(span, triangle.Plane(), &p, &t)) {
    case RN_NULL_CLASS_ID:
	// Span does not intersect triangle plane
	return RN_NULL_CLASS_ID;
	
    case R3_POINT_CLASS_ID:
	// Span intersects triangle plane in a single point
	// Check whether triangle contains intersection point
	if (!R3Contains(triangle, p)) return RN_NULL_CLASS_ID;
	if (hit_point1) *hit_point1 = p;
	if (hit_t1) *hit_t1 = t;
	return R3_POINT_CLASS_ID;

    case R3_SPAN_CLASS_ID: 
	// Span lies in triangle plane 
	{
	    // Clip span to all edges of triangle
	    R3Span s(span);
	    const R3Point *p1 = &(triangle.Vertex(2)->Position());
	    for (int i = 0; i < 3; i++) {
		// Get kth point
		const R3Point *p2 = &(triangle.Vertex(i)->Position());
		
		// Construct edge plane (containing edge and normal to triangle)
		R3Plane edge_plane(*p2, triangle.Normal(), *p2 - *p1);
		
		// Clip span to edge plane
		RNClassID s_type = s.Clip(edge_plane);
		if (s_type == RN_NULL_CLASS_ID) {
		    return RN_NULL_CLASS_ID;
		}
		else if (s_type == R3_POINT_CLASS_ID) {
		    if (hit_point1) *hit_point1 = s.Start();
		    if (hit_t1) *hit_t1 = span.T(s.Start());
		    return R3_POINT_CLASS_ID;
		}
		
		// Remember previous point
		p1 = p2;
	    }
	    
	    // Span intersects triangle in span
	    if (hit_point1) *hit_point1 = s.Start();
	    if (hit_t1) *hit_t1 = span.T(s.Start());
	    // if (hit_point2) *hit_point2 = s.End();
	    // if (hit_t2) *hit_t2 = span.T(s.End());
	    return R3_SPAN_CLASS_ID;
        }
	break;

    default:
	RNAbort("Invalid intersection result");
	break;
    }

    // Should never get here
    RNAbort("Should never get here");
    return RN_NULL_CLASS_ID;
}



RNClassID R3Intersects(const R3Span& span, const R3Halfspace& halfspace, 
    R3Point *hit_point, RNScalar *hit_t)
{
    // Intersect halfspace with span's ray
    RNScalar s;
    switch (R3Intersects(halfspace, span.Ray(), hit_point, &s)) {
    case RN_NULL_CLASS_ID:
	// Span and halfspace do not intersect
	return RN_NULL_CLASS_ID;

    case R3_RAY_CLASS_ID:
    case R3_SPAN_CLASS_ID:
	// Halfspace contains span
	return R3_SPAN_CLASS_ID;

    case R3_POINT_CLASS_ID: 
	// Check if intersection point is on span
	if (RNIsGreater(s, span.Length())) {
	    // Span's ray intersects halfspace beyond span endpoint
	    return RN_NULL_CLASS_ID;
	}
	else {
	    // Span and halfspace intersect at point
	    if (hit_t) *hit_t = s;
	    return R3_POINT_CLASS_ID;
	}
	break;

    default:
	RNAbort("Illegal intersection");
	break;
    }

    // Should never get here
    RNAbort("Should never get here");
    return RN_NULL_CLASS_ID;
}



RNClassID R3Intersects(const R3Span& span, const R3Shape& shape, 
    R3Point *hit_point, RNScalar *hit_t)
{
    // Intersect shape with span's ray
    RNScalar s;
    switch (R3Intersects(shape, span.Ray(), hit_point, NULL, &s)) {
    case RN_NULL_CLASS_ID:
	// Span and shape do not intersect
	return RN_NULL_CLASS_ID;

    case R3_RAY_CLASS_ID:
    case R3_SPAN_CLASS_ID:
	// Shape contains span
	return R3_SPAN_CLASS_ID;

    case R3_POINT_CLASS_ID: 
	// Check if intersection point is on span
	if (RNIsGreater(s, span.Length())) {
	    // Span's ray intersects shape beyond span endpoint
	    return RN_NULL_CLASS_ID;
	}
	else {
	    // Span and shape intersect at point
	    if (hit_t) *hit_t = s;
	    return R3_POINT_CLASS_ID;
	}
	break;

    default:
	RNAbort("Illegal intersection");
	break;
    }

    // Should never get here
    RNAbort("Should never get here");
    return RN_NULL_CLASS_ID;
}



RNClassID R3Intersects(const R3Plane& plane1, const R3Plane& plane2, 
    R3Line *result)
{
    // Check if planes are parallel
    if (R3Parallel(plane1, plane2)) {
	// Planes are parallel - check if planes are equal
	if (RNIsEqual(plane1.D(), plane2.D())) {
	    // Planes are coincident
	    return R3_PLANE_CLASS_ID;
	}
	else {
	    // Planes are parallel, but not coincident
	    return RN_NULL_CLASS_ID;
	}
    }
    else {
	// Planes are not parallel -- they must intersect
	if (result) {
	    // Compute direction vector of intersection line
	    R3Vector v = plane1.Normal();
	    v.Cross(plane2.Normal());

	    // Find v's largest dimension
	    RNDimension dim = v.MaxDimension();
	    RNDimension dim1 = (dim+1) % 3;
	    RNDimension dim2 = (dim+2) % 3;
	    assert(RNIsNotZero(v[dim]));

	    // Solve for point on zero plane of v's largest dimension
	    R3Point p;
	    p[dim] = 0.0;
	    p[dim1] = (plane1[dim2] * plane2.D() - plane2[dim2] * plane1.D()) / v[dim];
	    p[dim2] = (plane2[dim1] * plane1.D() - plane1[dim1] * plane2.D()) / v[dim];
	    *result = R3Line(p, v);
	}

	// Planes intersect in a line
	return R3_LINE_CLASS_ID;
    }
}



RNClassID R3Intersects(const R3Plane& plane1, const R3Plane& plane2, const R3Plane& plane3, 
    R3Point *hit_point)
{
    // From Graphics Gems I, page 305
    RNScalar denom = R4MatrixDet3(plane1.Normal().X(), plane2.Normal().X(), plane3.Normal().X(), 
				  plane1.Normal().Y(), plane2.Normal().Y(), plane3.Normal().Y(), 
				  plane1.Normal().Z(), plane2.Normal().Z(), plane3.Normal().Z());

    // Check for any two parallel planes 
    if (RNIsZero(denom)) {
	// Check for case of all planes coincident
	if (R3Contains(plane1, plane2) && R3Contains(plane1, plane3)) return R3_PLANE_CLASS_ID;
	else return RN_NULL_CLASS_ID;
    }

    // Compute unique intersection point
    if (hit_point) {
	*hit_point = R3zero_point;
	*hit_point += (plane2.Normal() % plane3.Normal()) * -(plane1.D());
	*hit_point += (plane3.Normal() % plane1.Normal()) * -(plane2.D());
	*hit_point += (plane1.Normal() % plane2.Normal()) * -(plane3.D());
	*hit_point /= denom;
    }

    // Return intersection type
    return R3_POINT_CLASS_ID;
}



RNClassID R3Intersects(const R3Plane& plane, const R3Halfspace& halfspace)
{
    // Check if normals
    if (R3Parallel(plane.Normal(), halfspace.Normal())) {
	// Check if halfspace intersects any point on plane
	return R3Intersects(halfspace, plane.Point());
    }
    else {
	// Normals are not parallel -- halfspace and plane must intersect
	return RN_UNKNOWN_CLASS_ID;
    }
}



RNClassID R3Intersects(const R3Plane& plane, const R3Box& box)
{
    // Return whether sphere intersects plane
    RNOctant octant = plane.Normal().Octant();
    if (R3SignedDistance(plane, box.Corner(octant)) < 0) return RN_NULL_CLASS_ID;
    if (R3SignedDistance(plane, box.Corner(~octant & 0x7)) > 0) return RN_NULL_CLASS_ID;
    return RN_UNKNOWN_CLASS_ID;
}




RNClassID R3Intersects(const R3Plane& plane, const R3OrientedBox& box)
{
    // Check if box is empty
    if (box.IsEmpty()) return RN_NULL_CLASS_ID;

    // Transform plane into std coordinate system and intersect there
    R3Plane transformed_plane = plane;
    R3Affine cs_to_std_xform(box.CoordSystem().InverseMatrix());
    transformed_plane.Transform(cs_to_std_xform);
    R3Box axis_aligned_box(-box.Radius(0), -box.Radius(1), -box.Radius(2), box.Radius(0), box.Radius(1), box.Radius(2));
    return R3Intersects(transformed_plane, axis_aligned_box);
}



RNClassID R3Intersects(const R3Plane& plane, const R3Sphere& sphere)
{
    // Return whether sphere intersects plane
    RNScalar d = R3Distance(plane, sphere.Center());
    if (RNIsGreater(d, sphere.Radius())) return RN_NULL_CLASS_ID;
    else if (RNIsEqual(d, sphere.Radius())) return R3_POINT_CLASS_ID;
    else return RN_UNKNOWN_CLASS_ID;
}



RNClassID R3Intersects(const R3Plane& plane, const R3Shape& shape)
{
    // Return whether plane intersects shape
    return shape.Intersects(plane);
}



RNClassID R3Intersects(const R3Plane& plane, const R3Ellipsoid& ellipsoid,
                       R3Ellipse *result)
{
  // Compute transformation from ellipse to canonical coordinate systems
  R3Affine xform(ellipsoid.CoordSystem().Matrix());

  // Transform plane to canonical coordinate system
  R3Plane transformed_plane(plane);
  transformed_plane.InverseTransform(xform);

  // Compute intersection
  R3Ellipse ellipse(R3xyz_coordinate_system, R2null_vector); 

  // If no intersection
  RNAbort("Not implemented");
  return RN_NULL_CLASS_ID;

  // Transform intersection back to original coordinate system
  if (result) {
    *result = ellipse;
    result->Transform(xform);
  }

  // Return success
  return R3_ELLIPSE_CLASS_ID;
}



RNClassID R3Intersects(const R3Triangle& triangle, const R3Halfspace& halfspace)
{
    // Check if triangle lies above plane
    if (R3Splits(halfspace.Plane(), triangle) != RN_BELOW) return RN_UNKNOWN_CLASS_ID;
    else return RN_NULL_CLASS_ID;
}



RNClassID R3Intersects(const R3Halfspace& halfspace1, const R3Halfspace& halfspace2)
{
    // Check halfspace directions
    if (R3Contains(halfspace1.Normal(), -(halfspace2.Normal()))) {
	// Halfspace plane normals point in opposite directions
	return R3Intersects(halfspace1, halfspace2.Plane().Point());
    }
    else { 
	// Halfspaces are not anti-parallel -- they must intersect
	return RN_UNKNOWN_CLASS_ID;
    }
}



RNClassID R3Intersects(const R3Halfspace& halfspace, const R3Circle& circle)
{
    // Return whether halfspace intersects circle
    RNScalar d = R3SignedDistance(halfspace.Plane(), circle.Center());
    if (RNIsPositiveOrZero(d)) return RN_UNKNOWN_CLASS_ID;
    else if (RNIsLess(d, -(circle.Radius()))) return RN_NULL_CLASS_ID;
    else {
        RNScalar cos_theta = halfspace.Plane().Normal().Dot(circle.Normal());
	if (cos_theta < 0.0) cos_theta = halfspace.Plane().Normal().Dot(-circle.Normal());
	if (RNIsGreaterOrEqual(d, -(cos_theta * circle.Radius()))) return RN_UNKNOWN_CLASS_ID;
	else return RN_NULL_CLASS_ID;
    }
}



RNClassID R3Intersects(const R3Halfspace& halfspace, const R3Box& box)
{
    // Return whether halfspace intersects box
    RNOctant octant = halfspace.Normal().Octant();
    if (!R3Contains(halfspace, box.Corner(octant))) return RN_NULL_CLASS_ID;
    return RN_UNKNOWN_CLASS_ID;
}



RNClassID R3Intersects(const R3Halfspace& halfspace, const R3OrientedBox& box)
{
    // Check if box is empty
    if (box.IsEmpty()) return RN_NULL_CLASS_ID;

    // Transform halfspace into std coordinate system and intersect there
    R3Halfspace transformed_halfspace = halfspace;
    R3Affine cs_to_std_xform(box.CoordSystem().InverseMatrix());
    transformed_halfspace.Transform(cs_to_std_xform);
    R3Box axis_aligned_box(-box.Radius(0), -box.Radius(1), -box.Radius(2), box.Radius(0), box.Radius(1), box.Radius(2));
    return R3Intersects(transformed_halfspace, axis_aligned_box);
}



RNClassID R3Intersects(const R3Halfspace& halfspace, const R3Sphere& sphere)
{
    // Return whether sphere intersects halfspace
    RNScalar d = R3SignedDistance(halfspace.Plane(), sphere.Center());
    if (RNIsLess(d, -(sphere.Radius()))) return RN_NULL_CLASS_ID;
    else if (RNIsEqual(d, -(sphere.Radius()))) return R3_POINT_CLASS_ID;
    else return RN_UNKNOWN_CLASS_ID;
}



RNClassID R3Intersects(const R3Halfspace& halfspace, const R3Cylinder& cylinder)
{
    // Return whether halfspace intersects cylinder
    if (R3Intersects(halfspace, cylinder.Top())) return RN_UNKNOWN_CLASS_ID;
    if (R3Intersects(halfspace, cylinder.Base())) return RN_UNKNOWN_CLASS_ID;
    return RN_NULL_CLASS_ID;
}



RNClassID R3Intersects(const R3Halfspace& halfspace, const R3Cone& cone)
{
    // Return whether halfspace intersects cone
    if (R3Intersects(halfspace, cone.Apex())) return RN_UNKNOWN_CLASS_ID;
    if (R3Intersects(halfspace, cone.Base())) return RN_UNKNOWN_CLASS_ID;
    return RN_NULL_CLASS_ID;
}



RNClassID R3Intersects(const R3Halfspace& halfspace, const R3Shape& shape)
{
    // Return whether halfspace intersects shape
    return shape.Intersects(halfspace);
}


RNClassID R3Intersects(const R3Triangle& triangle, const R3Box& box)
{
  // Check triangle bounding box and plane
  if (!R3Intersects(triangle.BBox(), box)) return RN_NULL_CLASS_ID;
  if (!R3Intersects(triangle.Plane(), box)) return RN_NULL_CLASS_ID;

  // Make polygon
  R3Point points[16];
  points[0] = triangle.Vertex(0)->Position();
  points[1] = triangle.Vertex(1)->Position();
  points[2] = triangle.Vertex(2)->Position();
  int npoints = 3;

  // Clip polygon against each side
  for (RNDirection dir = RN_LO; dir <= RN_HI; dir++) {
    for (RNDimension dim = RN_X; dim <= RN_Z; dim++) {
      points[npoints] = points[0];
      for (int i = 0; i < npoints; i++) {
        R3Point& p1 = points[i];
        R3Point& p2 = points[i+1];
        RNScalar d1 = p1[dim] - box[dir][dim];
        RNScalar d2 = p2[dim] - box[dir][dim];
        if (dir == RN_LO) { d1 = -d1; d2 = -d2; }
        if (d1 < 0) { // Inside
          if (d2 > 0) { // Outside
            // Insert a point at crossing
            RNScalar denom = d2 + -d1;
            R3Point p = (d2/denom)*p1 + (-d1/denom)*p2;
            for (int j = npoints; j > i; j--) points[j+1] = points[j];
            points[i+1] = p;
            npoints++;
            i += 2;
          }
        }
        else if (d1 > 0) { // Outside
          if (d2 < 0) { // Inside
            // Replace p1 with point at crossing
            RNScalar denom = -d2 + d1;
            R3Point p = (-d2/denom)*p1 + (d1/denom)*p2;
            points[i] = p;
          }
          else {
            // Remove p1
            for (int j = i; j < npoints; j++) points[j] = points[j+1];
            npoints--;
            i--;
          }
        }
      }

      // Check number of points
      if (npoints == 0) return RN_NULL_CLASS_ID;
      assert(npoints < 16);
      assert(npoints > 0);
    }
  }

  // Triangle survived all clips
  if (npoints <= 0) return RN_NULL_CLASS_ID;
  else if (npoints == 1) return R3_POINT_CLASS_ID;
  else if (npoints == 2) return R3_SPAN_CLASS_ID;
  else return R3_POLYGON_CLASS_ID;
}



RNClassID R3Intersects(const R3Box& box1, const R3Box& box2, 
    R3Box *result)
{
    // Compute whether box intersects box
    if ((RNIsGreater(box1.XMin(), box2.XMax())) ||
        (RNIsGreater(box1.YMin(), box2.YMax())) ||
        (RNIsGreater(box1.ZMin(), box2.ZMax())) ||
        (RNIsGreater(box2.XMin(), box1.XMax())) ||
        (RNIsGreater(box2.YMin(), box1.YMax())) ||
        (RNIsGreater(box2.ZMin(), box1.ZMax()))) return RN_NULL_CLASS_ID;

    // Compute intersection
    if (result) {
        RNCoord xmin = box1.XMin();
        if (xmin < box2.XMin()) xmin = box2.XMin();
        RNCoord ymin = box1.YMin();
        if (ymin < box2.YMin()) ymin = box2.YMin();
        RNCoord zmin = box1.ZMin();
        if (zmin < box2.ZMin()) zmin = box2.ZMin();
        RNCoord xmax = box1.XMax();
        if (xmax > box2.XMax()) xmax = box2.XMax();
        RNCoord ymax = box1.YMax();
        if (ymax > box2.YMax()) ymax = box2.YMax();
        RNCoord zmax = box1.ZMax();
        if (zmax > box2.ZMax()) zmax = box2.ZMax();
        *result = R3Box(xmin, ymin, zmin, xmax, ymax, zmax);
    }
    
    // Return ID of intersection ???
    return RN_UNKNOWN_CLASS_ID;
}



RNClassID R3Intersects(const R3Box& box1, const R3OrientedBox& box2)
{
    // Check if box is empty
    if (box1.IsEmpty()) return RN_NULL_CLASS_ID;
    if (box2.IsEmpty()) return RN_NULL_CLASS_ID;

    // Check for intersection with all halfspaces
    for (int dir = 0; dir < 2; dir++) {
      for (int dim = 0; dim < 3; dim++) {
        R3Halfspace halfspace(-box2.Plane(dir, dim), 0);
        if (!R3Intersects(halfspace, box1)) {
          return RN_NULL_CLASS_ID;
        }
      }
    }

    // Passed all tests
    return RN_UNKNOWN_CLASS_ID;
}



RNClassID R3Intersects(const R3Box& box, const R3Sphere& sphere)
{
    // Return whether box intersects sphere
    return RNIsZero(R3Distance(box, sphere));
}



RNClassID R3Intersects(const R3Box& box, const R3Shape& shape)
{
    // Return whether box intersects shape
    return shape.Intersects(box);
}



RNClassID R3Intersects(const R3OrientedBox& box1, const R3OrientedBox& box2)
{
    // Check if box is empty
    if (box1.IsEmpty()) return RN_NULL_CLASS_ID;
    if (box2.IsEmpty()) return RN_NULL_CLASS_ID;

    // Check for intersection with all halfspaces
    for (int dir = 0; dir < 2; dir++) {
      for (int dim = 0; dim < 3; dim++) {
        R3Halfspace halfspace(-box2.Plane(dir, dim), 0);
        if (!R3Intersects(halfspace, box1)) {
          return RN_NULL_CLASS_ID;
        }
      }
    }

    // Passed all tests
    return RN_UNKNOWN_CLASS_ID;
}



RNClassID R3Intersects(const R3OrientedBox& box, const R3Sphere& sphere)
{
    // Check if box is empty
    if (box.IsEmpty()) return RN_NULL_CLASS_ID;

    // Check for intersection with all halfspaces
    for (int dir = 0; dir < 2; dir++) {
      for (int dim = 0; dim < 3; dim++) {
        R3Halfspace halfspace(-box.Plane(dir, dim), 0);
        if (!R3Intersects(halfspace, sphere)) {
          return RN_NULL_CLASS_ID;
        }
      }
    }

    // Passed all tests
    return RN_UNKNOWN_CLASS_ID;
}



RNClassID R3Intersects(const R3OrientedBox& box, const R3Shape& shape)
{
    // Check if box is empty
    if (box.IsEmpty()) return RN_NULL_CLASS_ID;

    // Check for intersection with all halfspaces
    for (int dir = 0; dir < 2; dir++) {
      for (int dim = 0; dim < 3; dim++) {
        R3Halfspace halfspace(-box.Plane(dir, dim), 0);
        if (!R3Intersects(halfspace, shape)) {
          return RN_NULL_CLASS_ID;
        }
      }
    }

    // Passed all tests
    return RN_UNKNOWN_CLASS_ID;
}



RNClassID R3Intersects(const R3Sphere& sphere1, const R3Sphere& sphere2)
{
    // Return whether sphere intersects sphere
    R3Vector v = sphere1.Center() - sphere2.Center();
    RNScalar distance_squared = v.X() * v.X() + v.Y() * v.Y() + v.Z() * v.Z();
    RNScalar radius_sum = sphere1.Radius() + sphere2.Radius();
    RNScalar radius_sum_squared = radius_sum * radius_sum;
    if (RNIsGreater(distance_squared, radius_sum_squared)) return RN_NULL_CLASS_ID;
    else if (RNIsEqual(distance_squared, radius_sum_squared)) return R3_POINT_CLASS_ID;
    else return RN_UNKNOWN_CLASS_ID;
}



RNClassID R3Intersects(const R3Sphere& sphere, const R3Shape& shape)
{
    // Return whether sphere intersects shape
    return shape.Intersects(sphere);
}



RNClassID R3Intersects(const R3Shape& shape1, const R3Shape& shape2)
{
    // Return whether shape intersects shape
    return shape1.Intersects(shape2);
}





