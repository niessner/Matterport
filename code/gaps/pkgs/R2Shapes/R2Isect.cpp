/* Source file for the intersection utility */



/* Include files */

#include "R2Shapes/R2Shapes.h"



RNClassID R2Intersects(const R2Point& point1, const R2Point& point2)
{
    // Check if two points are the same within tolerance
    if (R2Contains(point1, point2)) {
	// Points are same
	return R2_POINT_CLASS_ID;
    }
    else {
	// Points are different
	return RN_NULL_CLASS_ID;
    }
}



RNClassID R2Intersects(const R2Point& point, const R2Line& line)
{
    // Check if line contains point
    if (R2Contains(line, point)) {
	// Point lies on line
	return R2_POINT_CLASS_ID;
    }
    else {
	// Point is not on line
	return RN_NULL_CLASS_ID;
    }
}



RNClassID R2Intersects(const R2Point& point, const R2Ray& ray, RNScalar *hit_t)
{
    // Check if ray contains point
    if (R2Contains(ray, point)) {
	// Point lies on ray
	if (hit_t) *hit_t = ray.T(point);
	return R2_POINT_CLASS_ID;
    }
    else {
	// Point is not on ray
	return RN_NULL_CLASS_ID;
    }
}



RNClassID R2Intersects(const R2Point& point, const R2Span& span, RNScalar *hit_t)
{
    // Check if span contains point
    if (R2Contains(span, point)) {
	// Point lies on span
	if (hit_t) *hit_t = span.T(point);
	return R2_POINT_CLASS_ID;
    }
    else {
	// Point is not on span
	return RN_NULL_CLASS_ID;
    }
}



RNClassID R2Intersects(const R2Point& point, const R2Halfspace& halfspace)
{
    // Check whether halfspace contains point
    if (R2Contains(halfspace, point)) {
	// Point lies in halfspace
	return R2_POINT_CLASS_ID;
    }
    else {
	// Point does not lie in halfspace
	return RN_NULL_CLASS_ID;
    }
}



RNClassID R2Intersects(const R2Point& point, const R2Box& box)
{
    // Check whether box contains point
    if (R2Contains(box, point)) {
	// Point lies in box
	return R2_POINT_CLASS_ID;
    }
    else {
	// Point does not lie in box
	return RN_NULL_CLASS_ID;
    }
}



RNClassID R2Intersects(const R2Point& point, const R2Circle& circle)
{
    // Check whether circle contains point
    if (R2Contains(circle, point)) {
	// Point lies in circle
	return R2_POINT_CLASS_ID;
    }
    else {
	// Point does not lie in circle
	return RN_NULL_CLASS_ID;
    }
}



RNClassID R2Intersects(const R2Line& line1, const R2Line& line2, 
    R2Point *hit_point)
{
    // Get variables in more convenient form
    RNScalar a1 = line1.A();
    RNScalar b1 = line1.B();
    RNScalar c1 = line1.C();
    RNScalar a2 = line2.A();
    RNScalar b2 = line2.B();
    RNScalar c2 = line2.C();

    // Check if line and ray are parallel
    RNScalar denom = a2*b1 - a1*b2;
    if (RNIsZero(denom)) {
	// Check if lines are coincident
	if (R2Contains(line1, line2.Point())) {
	    // Lines are coincident
	    return R2_LINE_CLASS_ID;
	}
	else {
	    // Lines are parallel and do not intersect
	    return RN_NULL_CLASS_ID;
	}
    }
    else {
	// Find hit point
        if (hit_point) {
	    (*hit_point)[0] = (b2*c1 - b1*c2) / denom;
	    (*hit_point)[1] = (a1*c2 - a2*c1) / denom;
	}
	return R2_POINT_CLASS_ID;
    }
}



RNClassID R2Intersects(const R2Line& line, const R2Ray& ray, 
    R2Point *hit_point, RNScalar *hit_t)
{
    // Get variables in more convenient form
    const R2Vector& n1 = line.Normal();
    const RNScalar c1 = line.C();
    const R2Vector p2 = ray.Start().Vector();
    const R2Vector& v2 = ray.Vector();

    // Check if line and ray are parallel
    RNScalar denom = n1.Dot(v2);
    if (RNIsZero(denom)) {
	// Check if line and ray are coincident
	if (R2Contains(line, ray.Start())) {
	    // Line and ray are coincident
	    return R2_RAY_CLASS_ID;
	}
	else {
	    // Lines are parallel and do not intersect
	    return RN_NULL_CLASS_ID;
	}
    }
    else {
	// Find parametric value on ray of intersection point between lines
        RNScalar t = (c1 + n1.Dot(p2)) / -denom;

	// Check if intersection point is on ray
	if (RNIsPositiveOrZero(t)) {
	    if (hit_point) *hit_point = ray.Point(t);
	    if (hit_t) *hit_t = t;
	    return R2_POINT_CLASS_ID;
	}
	else {
	    // Line and ray do not intersect
	    return RN_NULL_CLASS_ID;
	}
    }
}



RNClassID R2Intersects(const R2Line& line, const R2Span& span, 
    R2Point *hit_point, RNScalar *hit_t)
{
    // Intersect line with span's ray
    RNScalar s;
    switch (R2Intersects(line, span.Ray(), hit_point, &s)) {
    case RN_NULL_CLASS_ID:
	// Span and line do not intersect
	return RN_NULL_CLASS_ID;

    case R2_RAY_CLASS_ID:
	// Line contains span
	return R2_SPAN_CLASS_ID;

    case R2_POINT_CLASS_ID: 
	// Check if intersection point is on span
	if (RNIsGreater(s, span.Length())) {
	    // Span's ray intersects line beyond span endpoint
	    return RN_NULL_CLASS_ID;
	}
	else {
	    // Span and line intersect at point
	    if (hit_t) *hit_t = s;
	    return R2_POINT_CLASS_ID;
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



RNClassID R2Intersects(const R2Line& line, const R2Halfspace& halfspace, 
    R2Ray *result)
{
    RNAbort("Not implemented");
    return RN_NULL_CLASS_ID;
}



RNClassID R2Intersects(const R2Line& line, const R2Box& box, 
    R2Point *hit_point)
{
    RNAbort("Not implemented");
    return RN_NULL_CLASS_ID;
}



RNClassID R2Intersects(const R2Line& line, const R2Circle& circle, 
    R2Point *hit_point)
{
    RNAbort("Not implemented");
    return RN_NULL_CLASS_ID;
}



RNClassID R2Intersects(const R2Ray& ray1, const R2Ray& ray2, 
    R2Point *hit_point, RNScalar *hit_t1, RNScalar *hit_t2)
{
    // Get variables in more convenient form
    const R2Vector p1 = ray1.Start().Vector();
    const R2Vector p2 = ray2.Start().Vector();
    const R2Vector& v1 = ray1.Vector();
    const R2Vector& v2 = ray2.Vector();
    const R2Vector& n1 = ray1.Normal();
    const R2Vector& n2 = ray2.Normal();
    const RNScalar c1 = ray1.Line().C();
    const RNScalar c2 = ray2.Line().C();

    // Check if ray1 and ray2 are parallel
    RNScalar denom = n1.Dot(v2);
    if (RNIsZero(denom)) {
	// Check if ray1 and ray2 are coincident
	if (R2Contains(ray1.Line(), ray2.Start())) {
	    // Ray1 and ray2 share the same line
	    if (R2Contains(v1, v2)) {
	        // Rays are colinear and in same direction
	        return R2_RAY_CLASS_ID;
	    }
	    else {
	        // Rays are colinear and in opposite direction
	        if (R2Contains(ray1, ray2.Start())) return R2_SPAN_CLASS_ID;
		else return RN_NULL_CLASS_ID;
	    }
	}
	else {
	    // Rays are parallel and do not intersect
	    return RN_NULL_CLASS_ID;
	}
    }
    else {
	// Find parametric values of intersection point 
        RNScalar t1 = (c2 + n2.Dot(p1)) / denom;
        RNScalar t2 = (c1 + n1.Dot(p2)) / -denom;

	// Check if intersection point is on ray2
	if (RNIsPositiveOrZero(t1) && RNIsPositiveOrZero(t2)) {
	    if (hit_point) *hit_point = ray2.Point(t2);
	    if (hit_t1) *hit_t1 = t1;
	    if (hit_t2) *hit_t2 = t2;
	    return R2_POINT_CLASS_ID;
	}
	else {
	    // Ray1 and ray2 do not intersect
	    return RN_NULL_CLASS_ID;
	}
    }
}



RNClassID R2Intersects(const R2Ray& ray, const R2Span& span, 
    R2Point *hit_point, RNScalar *hit_tray, RNScalar *hit_tspan)
{
    // Intersect ray with span's ray
    RNScalar s;
    switch (R2Intersects(ray, span.Ray(), hit_point, hit_tray, &s)) {
    case RN_NULL_CLASS_ID:
	// Span and ray do not intersect
	return RN_NULL_CLASS_ID;

    case R2_RAY_CLASS_ID:
	// Ray contains span
	return R2_SPAN_CLASS_ID;

    case R2_POINT_CLASS_ID: 
	// Check if intersection point is on span
	if (RNIsGreater(s, span.Length())) {
	    // Span's ray intersects ray beyond span endpoint
	    return RN_NULL_CLASS_ID;
	}
	else {
	    // Span and ray intersect at point
	    if (hit_tspan) *hit_tspan = s;
	    return R2_POINT_CLASS_ID;
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



RNClassID R2Intersects(const R2Ray& ray, const R2Halfspace& halfspace, 
    R2Point *hit_point, RNScalar *hit_t)
{
    RNAbort("Not implemented");
    return RN_NULL_CLASS_ID;
}



RNClassID R2Intersects(const R2Ray& ray, const R2Box& box, 
    R2Point *hit_point1, R2Vector *hit_normal1, RNScalar *hit_t1)
{
    RNClassID id = RN_NULL_CLASS_ID;
    R2Point point1;
    RNScalar tval1;

    // Check box 
    if (box.IsEmpty()) {
        // Box is empty
        return RN_NULL_CLASS_ID;
    }
    else if (R2Contains(box, ray.Start())) {
        // Box contains ray start point
	if (hit_point1) *hit_point1 = ray.Start();
	if (hit_normal1) *hit_normal1 = ray.Vector();
	if (hit_t1) *hit_t1 = 0.0;
	return R2_SPAN_CLASS_ID;
    }
    else {
	// Find parametric distance to front plane in each dimension
	for (RNDimension dim = RN_X; dim <= RN_Y; dim++) {
	    if (RNIsPositive(ray.Vector()[dim])) {
		RNScalar delta = box.Min()[dim] - ray.Start()[dim];
		if (delta < 0.0) continue;
		tval1 = delta / ray.Vector()[dim];
	    }
	    else if (RNIsNegative(ray.Vector()[dim])) {
		RNScalar delta = box.Max()[dim] - ray.Start()[dim];
		if (delta > 0.0) continue;
		tval1 = delta / ray.Vector()[dim];
	    }
	    else {
		continue;
	    }
	    
	    // Compute intersection point on plane
	    RNDimension dim1 = 1 - dim;
	    point1[dim] = (ray.Vector()[dim] > 0.0) ? box.Min()[dim] : box.Max()[dim];
	    point1[dim1] = ray.Start()[dim1] + ray.Vector()[dim1] * tval1;
	    // assert(RNIsEqual(point1[dim], ray.Start()[dim] + ray.Vector()[dim] * tval1));
	    
	    // Check if intersection point on plane is inside box
	    if ((RNIsLessOrEqual(point1[dim1], box.Max()[dim1])) &&
		(RNIsGreaterOrEqual(point1[dim1], box.Min()[dim1]))) {
		// Found first intersection point 
		id = R2_SPAN_CLASS_ID;
                if (hit_t1) *hit_t1 = tval1;
                if (hit_point1) *hit_point1 = point1;
                if (hit_normal1) {
	            *hit_normal1 = R2zero_vector;
		    if (RNIsNegative(ray.Vector()[dim])) *hit_normal1 = (*hit_normal1)[dim] = 1.0;
		    else (*hit_normal1)[dim] = -1.0;
	        }
		break;
	    }
	}

	// Return intersection 
	return id;
    }
}



RNClassID R2Intersects(const R2Ray& ray, const R2Circle& circle, 
    R2Point *hit_point1, R2Vector *hit_normal1, RNScalar *hit_t1)
{
    // Check if circle contains ray start point
    if (R2Contains(circle, ray.Start())) {
	if (hit_point1) *hit_point1 = ray.Start();
	if (hit_normal1) *hit_normal1 = ray.Vector();
	if (hit_t1) *hit_t1 = 0.0;
	return R2_SPAN_CLASS_ID;
    }
    else {
	// Check if ray points towards circle
	R2Vector EO = circle.Center() - ray.Start();
	RNScalar v = EO.Dot(ray.Vector());
	if (RNIsNegativeOrZero(v)) {
	    // Ray points away from circle
	    return RN_NULL_CLASS_ID;
	}
	else {
	    // Compute possible intersection -- Graphics Gems I, p. 388    
	    RNScalar disc = circle.Radius() * circle.Radius() - (EO.Dot(EO) - v*v);
	    if (RNIsNegative(disc)) {
		// No intersection
		return RN_NULL_CLASS_ID;
	    }
	    else {
		// Ray intersects circle (it grazes if disc is zero)
		// Compute first intersection
		if (hit_t1 || hit_point1 || hit_normal1) {
		    RNScalar d = sqrt(disc);
		    RNScalar t = v - d;
		    if (hit_t1) *hit_t1 = t;
		    if (hit_point1) *hit_point1 = ray.Start() + t * ray.Vector();
		    if (hit_normal1) {
			R2Point p, *pp;
			if (hit_point1) pp = hit_point1;
			else { p = ray.Start() + t * ray.Vector(); pp = &p; }
			*hit_normal1 = (*pp - circle.Centroid()) / circle.Radius();
		    }
		}

		// Compute second intersection
		// if (hit_t2 || hit_point2 || hit_normal2) {
		//     t = v + d;
		//     if (hit_point2) *hit_point2 = ray.Start() + t * ray.Vector();
		//     if (hit_t2) *hit_t2 = t;
		// }

		// Return result
		if (RNIsZero(disc)) return R2_POINT_CLASS_ID;
	        else return R2_SPAN_CLASS_ID;
	    }
	}
    }
}



RNClassID R2Intersects(const R2Span& span1, const R2Span& span2, 
    R2Point *hit_point, RNScalar *hit_t1, RNScalar *hit_t2)
{
    // Intersect span1 with span2's ray
    RNScalar s1, s2;
    switch (R2Intersects(span1, span2.Ray(), hit_point, &s1, &s2)) {
    case RN_NULL_CLASS_ID:
	// Spans do not intersect
	return RN_NULL_CLASS_ID;

    case R2_SPAN_CLASS_ID:
	// Span1 lies on span2's ray -- check if they overlap ???
	if ((R2Contains(span2, span1.Start())) || (R2Contains(span2, span1.End()))) {
	    // Spans are on same line, and overlap
	    return R2_SPAN_CLASS_ID;
	}
	else {
	    // Spans are on same line, but do not overlap
	    return RN_NULL_CLASS_ID;
	}
	break;

    case R2_POINT_CLASS_ID: 
	// Check if intersection point is on span
	if (RNIsGreater(s1, span1.Length())) {
	    // Span's ray intersects plane beyond span endpoint
	    return RN_NULL_CLASS_ID;
	}
	else {
	    // Spans intersect at point
	    if (hit_t1) *hit_t1 = s1;
	    if (hit_t2) *hit_t2 = s2;
	    return R2_POINT_CLASS_ID;
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



RNClassID R2Intersects(const R2Span& span, const R2Halfspace& halfspace, 
    R2Point *hit_point, RNScalar *hit_t)
{
    // Intersect halfspace with span's ray
    RNScalar s;
    switch (R2Intersects(halfspace, span.Ray(), hit_point, &s)) {
    case RN_NULL_CLASS_ID:
	// Span and halfspace do not intersect
	return RN_NULL_CLASS_ID;

    case R2_RAY_CLASS_ID:
    case R2_SPAN_CLASS_ID:
	// Halfspace contains span
	return R2_SPAN_CLASS_ID;

    case R2_POINT_CLASS_ID: 
	// Check if intersection point is on span
	if (RNIsGreater(s, span.Length())) {
	    // Span's ray intersects halfspace beyond span endpoint
	    return RN_NULL_CLASS_ID;
	}
	else {
	    // Span and halfspace intersect at point
	    if (hit_t) *hit_t = s;
	    return R2_POINT_CLASS_ID;
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



RNClassID R2Intersects(const R2Span& span, const R2Box& box, 
    R2Point *hit_point, RNScalar *hit_t)
{
    // Intersect box with span's ray
    RNScalar s;
    switch (R2Intersects(box, span.Ray(), hit_point, NULL, &s)) {
    case RN_NULL_CLASS_ID:
	// Span and box do not intersect
	return RN_NULL_CLASS_ID;

    case R2_SPAN_CLASS_ID:
	// Box contains span
        if (RNIsGreater(s, span.Length()))
            return RN_NULL_CLASS_ID;
        else
            return R2_SPAN_CLASS_ID;

    case R2_POINT_CLASS_ID: 
	// Check if intersection point is on span
	if (RNIsGreater(s, span.Length())) {
	    // Span's ray intersects box beyond span endpoint
	    return RN_NULL_CLASS_ID;
	}
	else {
	    // Span and box intersect at point
	    if (hit_t) *hit_t = s;
	    return R2_POINT_CLASS_ID;
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



RNClassID R2Intersects(const R2Span& span, const R2Circle& circle, 
    R2Point *hit_point, RNScalar *hit_t)
{
    // Intersect circle with span's ray
    RNScalar s;
    switch (R2Intersects(circle, span.Ray(), hit_point, NULL, &s)) {
    case RN_NULL_CLASS_ID:
	// Span and circle do not intersect
	return RN_NULL_CLASS_ID;

    case R2_SPAN_CLASS_ID:
	// Circle contains span
	return R2_SPAN_CLASS_ID;

    case R2_POINT_CLASS_ID: 
	// Check if intersection point is on span
	if (RNIsGreater(s, span.Length())) {
	    // Span's ray intersects circle beyond span endpoint
	    return RN_NULL_CLASS_ID;
	}
	else {
	    // Span and circle intersect at point
	    if (hit_t) *hit_t = s;
	    return R2_POINT_CLASS_ID;
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



RNClassID R2Intersects(const R2Halfspace& halfspace1, const R2Halfspace& halfspace2)
{
    // Check halfspace directions
    if (R2Contains(halfspace1.Normal(), -(halfspace2.Normal()))) {
	// Halfspace normals point in opposite directions
	return R2Intersects(halfspace1, halfspace2.Line().Point());
    }
    else { 
	// Halfspaces are not anti-parallel -- they must intersect
	return RN_UNKNOWN_CLASS_ID;
    }
}



RNClassID R2Intersects(const R2Halfspace& halfspace, const R2Box& box)
{
    // Return whether halfspace intersects box
    RNQuadrant quadrant = halfspace.Normal().Quadrant();
    if (!R2Contains(halfspace, box.Corner(quadrant))) return RN_NULL_CLASS_ID;
    return RN_UNKNOWN_CLASS_ID;
}



RNClassID R2Intersects(const R2Halfspace& halfspace, const R2Circle& circle)
{
    // Return whether circle intersects halfspace
    RNScalar d = R2SignedDistance(halfspace.Line(), circle.Center());
    if (RNIsLess(d, -(circle.Radius()))) return RN_NULL_CLASS_ID;
    else if (RNIsEqual(d, -(circle.Radius()))) return R2_POINT_CLASS_ID;
    else return RN_UNKNOWN_CLASS_ID;
}



RNClassID R2Intersects(const R2Box& box1, const R2Box& box2, 
    R2Box *result)
{
    // Compute whether box intersects box
    if ((RNIsGreater(box1.XMin(), box2.XMax())) ||
        (RNIsGreater(box1.YMin(), box2.YMax())) ||
        (RNIsGreater(box2.XMin(), box1.XMax())) ||
        (RNIsGreater(box2.YMin(), box1.YMax()))) return RN_NULL_CLASS_ID;

    // Compute intersection
    if (result) {
        RNCoord xmin = box1.XMin();
        if (xmin < box2.XMin()) xmin = box2.XMin();
        RNCoord ymin = box1.YMin();
        if (ymin > box2.YMin()) ymin = box2.YMin();
        RNCoord xmax = box1.XMax();
        if (xmax < box2.XMax()) xmax = box2.XMax();
        RNCoord ymax = box1.YMax();
        if (ymax > box2.YMax()) ymax = box2.YMax();
        *result = R2Box(xmin, ymin, xmax, ymax);
    }
    
    // Return ID of intersection ???
    return RN_UNKNOWN_CLASS_ID;
}



RNClassID R2Intersects(const R2Box& box, const R2Circle& circle)
{
    // Return whether box intersects circle
    return RNIsZero(R2Distance(box, circle));
}



RNClassID R2Intersects(const R2Circle& circle1, const R2Circle& circle2)
{
    // Return whether circle intersects circle
    R2Vector v = circle1.Center() - circle2.Center();
    RNScalar distance_squared = v.X() * v.X() + v.Y() * v.Y();
    RNScalar radius_sum = circle1.Radius() + circle2.Radius();
    RNScalar radius_sum_squared = radius_sum * radius_sum;
    if (RNIsGreater(distance_squared, radius_sum_squared)) return RN_NULL_CLASS_ID;
    else if (RNIsEqual(distance_squared, radius_sum_squared)) return R2_POINT_CLASS_ID;
    else return RN_UNKNOWN_CLASS_ID;
}





