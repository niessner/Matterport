/* Source file for GAPS distance utility */



/* Include files */

#include "R3Shapes/R3Shapes.h"



/* Public variables */

#ifdef R3_USER_DEFINED_SHAPE_TYPES
R3DistanceFunction *R3distance_functions[R3_NUM_SHAPE_TYPES][R3_NUM_SHAPE_TYPES] = { NULL };
#endif



/* Public functions */

int R3InitDistance()
{
#ifdef R3_USER_DEFINED_SHAPE_TYPES
   // Register default distance functions
    RNLength (*d1)(const R3Box&, const R3Box&) = &R3Distance;
    R3RegisterDistanceFunction(R3Box::CLASS_ID(), R3Box::CLASS_ID(), (R3DistanceFunction *) d1);
#endif

    // Return success 
    return TRUE;
}



void R3StopDistance()
{
}



RNLength R3SquaredDistance(const R3Point& point1, const R3Point& point2)
{
    // Return length of vector between points
    R3Vector v = point1 - point2;
    return v.Dot(v);
}



RNLength R3Distance(const R3Point& point1, const R3Point& point2)
{
    // Return length of vector between points
    R3Vector v = point1 - point2;
    return v.Length();
}



RNLength R3Distance(const R3Point& point, const R3Line& line)
{
    // Return distance from point to line (Riddle p. 904)
    R3Vector v = line.Vector();
    v.Cross(line.Point() - point);
    return v.Length();
}



RNLength R3Distance(const R3Point& point, const R3Ray& ray)
{
    // Check if start point is closest
    R3Vector v = point - ray.Start();
    RNScalar dir = v.Dot(ray.Vector());
    if (RNIsNegative(dir)) return v.Length();

    // Return distance from point to ray line
    return R3Distance(point, ray.Line());
}



RNLength R3Distance(const R3Point& point, const R3Span& span)
{
    // Check if span has zero length
    if (RNIsZero(span.Length())) return R3Distance(point, span.Start());

    // Check if start point is closest
    R3Vector v1 = point - span.Start();
    RNScalar dir1 = v1.Dot(span.Vector());
    if (RNIsNegative(dir1)) return v1.Length();

    // Check if end point is closest
    R3Vector v2 = point - span.End();
    RNScalar dir2 = v2.Dot(span.Vector());
    if (RNIsPositive(dir2)) return v2.Length();

    // Return distance from point to span line
    return R3Distance(point, span.Line());
}



RNLength R3Distance(const R3Point& point, const R3Plane& plane)
{
    // Return distance from point to plane
    RNScalar d = R3SignedDistance(plane, point);
    if (RNIsPositive(d)) return d;
    else if (RNIsNegative(d)) return -d;
    else return 0.0;
}



RNLength R3Distance(const R3Point& point, const R3Halfspace& halfspace)
{
    // Return distance from point to halfspace
    RNScalar d = R3SignedDistance(halfspace.Plane(), point);
    if (RNIsNegative(d)) return -d;
    else return 0.0;
}



RNLength R3Distance(const R3Point& point, const R3Box& box)
{
    // Check if box is empty
    if (box.IsEmpty()) return FLT_MAX;

    // Find axial distances from point to box
    RNLength dx, dy, dz;
    if (RNIsGreater(point.X(), box.XMax())) dx = point.X() - box.XMax();
    else if (RNIsLess(point.X(), box.XMin())) dx = box.XMin()- point.X();
    else dx = 0.0;
    if (RNIsGreater(point.Y(), box.YMax())) dy = point.Y() - box.YMax();
    else if (RNIsLess(point.Y(), box.YMin())) dy = box.YMin()- point.Y();
    else dy = 0.0;
    if (RNIsGreater(point.Z(), box.ZMax())) dz = point.Z() - box.ZMax();
    else if (RNIsLess(point.Z(), box.ZMin())) dz = box.ZMin()- point.Z();
    else dz = 0.0;
    
    // Return distance between point and closest point in box 
    if ((dy == 0.0) && (dz == 0.0)) return dx;
    else if ((dx == 0.0) && (dz == 0.0)) return dy;
    else if ((dx == 0.0) && (dy == 0.0)) return dz;
    else return sqrt(dx*dx + dy*dy + dz*dz);
}



RNLength R3Distance(const R3Point& point, const R3OrientedBox& box)
{
    // Check if box is empty
    if (box.IsEmpty()) return FLT_MAX;

    // Transform point into std coordinate system
    R3Point transformed_point = point;
    R3Affine cs_to_std_xform(box.CoordSystem().InverseMatrix());
    transformed_point.Transform(cs_to_std_xform);

    // Return distance in std coordinate system
    R3Box axis_aligned_box(-box.Radius(0), -box.Radius(1), -box.Radius(2), box.Radius(0), box.Radius(1), box.Radius(2));
    return R3Distance(transformed_point, axis_aligned_box);
}



RNLength R3Distance(const R3Point& point, const R3Sphere& sphere)
{
    // Return distance from point to sphere
    RNLength d = R3Distance(point, sphere.Center()) - sphere.Radius();
    return ((d > 0.0) ? d : 0.0);
}



RNLength R3Distance(const R3Point& point, const R3Cone& cone)
{
    // Return distance from point to sphere
    RNScalar t = cone.Axis().T(point);
    if (RNIsNegativeOrZero(t)) {
        RNAbort("Not implemented");
        return 0.0;
    }
    else if (RNIsGreater(t, cone.Axis().Length())) {
        return R3Distance(point, cone.Apex());
    }
    else {
        // Get distance between axis and point
        R3Point p = cone.Axis().Point(t);
	RNLength r = (1.0 - t / cone.Axis().Length()) * cone.Radius();
	RNLength d = R3Distance(p, point) - r;
	if (RNIsNegativeOrZero(d)) return 0.0;
	else return d;
    }

}



RNLength R3Distance(const R3Point& point, const R3Shape& shape)
{
    // Return distance from point to shape
    return shape.Distance(point);
}



RNLength R3Distance(const R3Line& line1, const R3Line& line2)
{
    // Return distance from line to line (Riddle p. 905)
    R3Vector v = line1.Vector();
    v.Cross(line2.Vector());
    RNLength vlength = v.Length();
    if (vlength == 0.0) return R3Distance(line1.Point(), line2);
    else v /= vlength;
    return fabs(v.Dot(line1.Point() - line2.Point()));
}



RNLength R3Distance(const R3Line& line, const R3Ray& ray)
{
    // Return distance from line to ray
    RNAbort("Not implemented");
    return 0.0;
}



RNLength R3Distance(const R3Line& line, const R3Span& span)
{
    RNAbort("Not implemented");
    return 0.0;
}



RNLength R3Distance(const R3Line& line, const R3Plane& plane)
{
    // Return distance from line to plane
    RNScalar d = R3SignedDistance(plane, line);
    if (RNIsPositive(d)) return d;
    else if (RNIsNegative(d)) return -d;
    else return 0.0;
}



RNLength R3Distance(const R3Line& line, const R3Halfspace& halfspace)
{
    // Return distance from line to halfspace
    RNScalar d = R3SignedDistance(halfspace.Plane(), line);
    if (RNIsNegative(d)) return -d;
    else return 0.0;
}



RNLength R3Distance(const R3Line& line, const R3Box& box)
{
    RNAbort("Not implemented");
    return 0.0;
}



RNLength R3Distance(const R3Line& line, const R3OrientedBox& box)
{
    // Check if box is empty
    if (box.IsEmpty()) return FLT_MAX;

    // Transform line into std coordinate system
    R3Line transformed_line = line;
    R3Affine cs_to_std_xform(box.CoordSystem().InverseMatrix());
    transformed_line.Transform(cs_to_std_xform);

    // Return distance in std coordinate system
    R3Box axis_aligned_box(-box.Radius(0), -box.Radius(1), -box.Radius(2), box.Radius(0), box.Radius(1), box.Radius(2));
    return R3Distance(transformed_line, axis_aligned_box);
}



RNLength R3Distance(const R3Line& line, const R3Sphere& sphere)
{
    // Return distance from line to sphere
    RNLength d = R3Distance(line, sphere.Center()) - sphere.Radius();
    return ((d > 0.0) ? d : 0.0);
}



RNLength R3Distance(const R3Line& line, const R3Shape& shape)
{
    // Return distance from line to shape
    return shape.Distance(line);
}



RNLength R3Distance(const R3Ray& ray1, const R3Ray& ray2)
{
    RNAbort("Not implemented");
    return 0.0;
}



RNLength R3Distance(const R3Ray& ray, const R3Span& span)
{
    // There's got to be a better way ???

    // Get vectors in more convenient form
    const R3Vector v1 = ray.Vector();
    const R3Vector v2 = span.Vector();

    // Compute useful intermediate values
    const RNScalar v1v1 = 1.0;  // v1.Dot(v1);
    const RNScalar v2v2 = 1.0;  // v2.Dot(v2);
    RNScalar v1v2 = v1.Dot(v2);
    RNScalar denom = v1v2*v1v2 - v1v1*v2v2;

    // Check if ray and span are parallel
    if (RNIsZero(denom)) {
        // Look at directions of vectors, then check relative starts and stops
        RNScalar dot = v1.Dot(v2);
        if (dot > 0) {
          R3Vector vse = span.End() - ray.Start();
          if (v1.Dot(vse) < 0) return R3Distance(span.End(), ray.Start());
          else return R3Distance(span.End(), ray.Line());
        }
        else {
          R3Vector vss = span.Start() - ray.Start();
          if (v1.Dot(vss) < 0) return R3Distance(span.Start(), ray.Start());
          else return R3Distance(span.Start(), ray.Line());
        }
    }
    else {
	// Find closest points
	const R3Vector p1 = ray.Start().Vector();
	const R3Vector p2 = span.Start().Vector();
	RNScalar p1v1 = v1.Dot(p1);
	RNScalar p2v2 = v2.Dot(p2);
	RNScalar p1v2 = v2.Dot(p1);
	RNScalar p2v1 = v1.Dot(p2);
	RNScalar ray_t = (v1v2*p2v2 + v2v2*p1v1 - v1v2*p1v2 - v2v2*p2v1) / denom;
	RNScalar span_t = (v1v2*p1v1 + v1v1*p2v2 - v1v2*p2v1 - v1v1*p1v2) / denom;
	R3Point ray_point = (ray_t <= 0.0) ? ray.Start() : ray.Point(ray_t);
	R3Point span_point = (span_t <= 0.0) ? span.Start() : 
	    (span_t >= span.Length()) ? span.End() : span.Ray().Point(span_t);
	RNLength distance = R3Distance(ray_point, span_point);
	return distance;
    }
}



RNLength R3Distance(const R3Ray& ray, const R3Plane& plane)
{
    // Return distance from ray to plane
    RNScalar d = R3SignedDistance(plane, ray);
    if (RNIsPositive(d)) return d;
    else if (RNIsNegative(d)) return -d;
    else return 0.0;
}



RNLength R3Distance(const R3Ray& ray, const R3Halfspace& halfspace)
{
    // Return distance from ray to halfspace
    RNScalar d = R3SignedDistance(halfspace.Plane(), ray);
    if (RNIsNegative(d)) return -d;
    else return 0.0;
}



RNLength R3Distance(const R3Ray& ray, const R3Box& box)
{
    RNAbort("Not implemented");
    return 0.0;
}



RNLength R3Distance(const R3Ray& ray, const R3OrientedBox& box)
{
    // Check if box is empty
    if (box.IsEmpty()) return FLT_MAX;

    // Transform ray into std coordinate system
    R3Ray transformed_ray = ray;
    R3Affine cs_to_std_xform(box.CoordSystem().InverseMatrix());
    transformed_ray.Transform(cs_to_std_xform);

    // Return distance in std coordinate system
    R3Box axis_aligned_box(-box.Radius(0), -box.Radius(1), -box.Radius(2), box.Radius(0), box.Radius(1), box.Radius(2));
    return R3Distance(transformed_ray, axis_aligned_box);
}



RNLength R3Distance(const R3Ray& ray, const R3Sphere& sphere)
{
    // Return distance from ray to sphere
    RNLength d = R3Distance(ray, sphere.Center()) - sphere.Radius();
    return ((d > 0.0) ? d : 0.0);
}



RNLength R3Distance(const R3Ray& ray, const R3Shape& shape)
{
    // Return distance from ray to shape
    return shape.Distance(ray);
}



RNLength R3Distance(const R3Span& span1, const R3Span& span2)
{
    // There's got to be a better way ???

    // Get vectors in more convenient form
    const R3Vector v1 = span1.Vector();
    const R3Vector v2 = span2.Vector();

    // Compute useful intermediate values
    const RNScalar v1v1 = 1.0;  // v1.Dot(v1);
    const RNScalar v2v2 = 1.0;  // v2.Dot(v2);
    RNScalar v1v2 = v1.Dot(v2);
    RNScalar denom = v1v2*v1v2 - v1v1*v2v2;

    // Check if span1 and span are parallel
    if (RNIsZero(denom)) {
        // Look at directions of vectors, then check relative starts and stops
        RNScalar dot = v1.Dot(v2);
        if (dot > 0) {
          R3Vector ves = span2.Start() - span1.End();
          R3Vector vse = span2.End() - span1.Start();
          if (v1.Dot(ves) > 0) return R3Distance(span1.End(), span2.Start());
          else if (v1.Dot(vse) < 0) return R3Distance(span1.Start(), span2.End());
          else return R3Distance(span1.Start(), span2.Line());
        }
        else {
          R3Vector vee = span2.End() - span1.End();
          R3Vector vss = span2.Start() - span1.Start();
          if (v1.Dot(vss) < 0) return R3Distance(span1.Start(), span2.Start());
          else if (v1.Dot(vee) > 0) return R3Distance(span1.End(), span2.End());
          else return R3Distance(span1.Start(), span2.Line());
        }
    }
    else {
	// Find closest points
	const R3Vector p1 = span1.Start().Vector();
	const R3Vector p2 = span2.Start().Vector();
	RNScalar p1v1 = v1.Dot(p1);
	RNScalar p2v2 = v2.Dot(p2);
	RNScalar p1v2 = v2.Dot(p1);
	RNScalar p2v1 = v1.Dot(p2);
	RNScalar span1_t = (v1v2*p2v2 + v2v2*p1v1 - v1v2*p1v2 - v2v2*p2v1) / denom;
	RNScalar span2_t = (v1v2*p1v1 + v1v1*p2v2 - v1v2*p2v1 - v1v1*p1v2) / denom;
        if ((span1_t >= 0) && (span1_t <= span1.Length()) &&
            (span2_t >= 0) && (span2_t <= span2.Length())) {
	    return R3Distance(span1.Point(span1_t), span2.Point(span2_t));
        }
        else {
            // There has to be a better way!!!
            RNLength d1A = R3Distance(span1.Start(), span2);
            RNLength d1B = R3Distance(span1.End(), span2);
            RNLength d2A = R3Distance(span2.Start(), span1);
            RNLength d2B = R3Distance(span2.End(), span1);
            RNLength d = FLT_MAX;
            if (d1A < d) d = d1A;
            if (d1B < d) d = d1B;
            if (d2A < d) d = d2A;
            if (d2B < d) d = d2B;
            return d;
        }
    }
}



RNLength R3Distance(const R3Span& span, const R3Plane& plane)
{
    // Return distance from span to plane
    RNScalar d = R3SignedDistance(plane, span);
    if (RNIsPositive(d)) return d;
    else if (RNIsNegative(d)) return -d;
    else return 0.0;
}



RNLength R3Distance(const R3Span& span, const R3Halfspace& halfspace)
{
    // Return distance from span to halfspace
    RNScalar d = R3SignedDistance(halfspace.Plane(), span);
    if (RNIsNegative(d)) return -d;
    else return 0.0;
}



RNLength R3Distance(const R3Span& span, const R3Box& box)
{
    RNAbort("Not implemented");
    return 0.0;
}



RNLength R3Distance(const R3Span& span, const R3OrientedBox& box)
{
    // Check if box is empty
    if (box.IsEmpty()) return FLT_MAX;

    // Transform span into std coordinate system
    R3Span transformed_span = span;
    R3Affine cs_to_std_xform(box.CoordSystem().InverseMatrix());
    transformed_span.Transform(cs_to_std_xform);

    // Return distance in std coordinate system
    R3Box axis_aligned_box(-box.Radius(0), -box.Radius(1), -box.Radius(2), box.Radius(0), box.Radius(1), box.Radius(2));
    return R3Distance(transformed_span, axis_aligned_box);
}



RNLength R3Distance(const R3Span& span, const R3Sphere& sphere)
{
    // Return distance from span to sphere
    RNLength d = R3Distance(span, sphere.Center()) - sphere.Radius();
    return ((d > 0.0) ? d : 0.0);
}



RNLength R3Distance(const R3Span& span, const R3Shape& shape)
{
    // Return distance from span to shape
    return shape.Distance(span);
}



RNLength R3Distance(const R3Plane& plane1, const R3Plane& plane2)
{
    // Return distance from plane to plane
    RNScalar d = R3SignedDistance(plane1, plane2);
    if (RNIsPositive(d)) return d;
    else if (RNIsNegative(d)) return -d;
    else return 0.0;
}



RNLength R3Distance(const R3Plane& plane, const R3Halfspace& halfspace)
{
    // Return distance from plane to halfspace
    RNScalar d = R3SignedDistance(halfspace.Plane(), plane);
    if (RNIsNegative(d)) return -d;
    else return 0.0;
}



RNLength R3Distance(const R3Plane& plane, const R3Box& box)
{
    // Return distance from plane to box
    RNScalar d = R3SignedDistance(plane, box);
    if (RNIsPositive(d)) return d;
    else if (RNIsNegative(d)) return -d;
    else return 0.0;
}



RNLength R3Distance(const R3Plane& plane, const R3OrientedBox& box)
{
    // Return distance from plane to box
    RNScalar d = R3SignedDistance(plane, box);
    if (RNIsPositive(d)) return d;
    else if (RNIsNegative(d)) return -d;
    else return 0.0;
}



RNLength R3Distance(const R3Plane& plane, const R3Sphere& sphere)
{
    // Return distance from plane to sphere
    RNScalar d = R3SignedDistance(plane, sphere);
    if (RNIsPositive(d)) return d;
    else if (RNIsNegative(d)) return -d;
    else return 0.0;
}



RNLength R3Distance(const R3Plane& plane, const R3Shape& shape)
{
    // Return distance from plane to shape
    return shape.Distance(plane);
}



RNLength R3SignedDistance(const R3Plane& plane, const R3Point& point)
{
    // Return signed distance from point to plane (Riddle p. 914)
    return (point.X()*plane.A() + point.Y()*plane.B() + point.Z()*plane.C() + plane.D());
}



RNLength R3SignedDistance(const R3Plane& plane, const R3Line& line)
{
    // Return signed distance from plane to line
    if (R3Parallel(plane, line)) {
	// Plane and line are parallel
	return R3SignedDistance(plane, line.Point());
    }
    else {
	// Plane and line are not parallel
	return 0.0;
    }
}



RNLength R3SignedDistance(const R3Plane& plane, const R3Ray& ray)
{
    // Return signed distance from plane to ray
    RNLength d1 = R3SignedDistance(plane, ray.Start());
    if (RNIsPositive(d1)) {
	// Start point is above plane
	RNScalar dot = ray.Vector().Dot(plane.Normal());
	if (RNIsNegative(dot)) return 0.0;
	else return d1;
    }
    else if (RNIsNegative(d1)) {
	// Start point is below plane
	RNScalar dot = ray.Vector().Dot(plane.Normal());
	if (RNIsPositive(dot)) return 0.0;
	else return d1;
    }
    else {
	// Start point is on plane
	return 0.0;
    }
}



RNLength R3SignedDistance(const R3Plane& plane, const R3Span& span)
{
    // Return signed distance from plane to span
    RNLength d1 = R3SignedDistance(plane, span.Start());
    if (RNIsPositive(d1)) {
	// Start point is above plane
	RNLength d2 = R3SignedDistance(plane, span.End());
	if (RNIsPositive(d2)) return ((d1 > d2) ? d2 : d1);
	else return 0.0;
    }
    else if (RNIsNegative(d1)) {
	// Start point is below plane
	RNLength d2 = R3SignedDistance(plane, span.End());
	if (RNIsNegative(d2)) return ((d1 > d2) ? d1 : d2);
	else return 0.0;
    }
    else {
	// Start point is on plane
	return 0.0;
    }
}



RNLength R3SignedDistance(const R3Plane& plane1, const R3Plane& plane2)
{
    // Return signed distance from plane to plane
    RNScalar dot = plane1.Normal().Dot(plane2.Normal());
    if (RNIsEqual(dot, 1.0)) return (plane1.D() - plane2.D());
    else if (RNIsEqual(dot, -1.0)) return (plane1.D() + plane2.D());
    else return 0.0;
}



RNLength R3SignedDistance(const R3Plane& plane, const R3Triangle& triangle)
{
    // Compute signed distance from plane to triangle
    RNLength max_distance = -FLT_MAX;
    RNLength min_distance = FLT_MAX;
    for (int i = 0; i < 3; i++) {
        RNLength distance = R3SignedDistance(plane, triangle.Vertex(i)->Position());
	if (distance < min_distance) min_distance = distance;
	if (distance > max_distance) max_distance = distance;
	if ((min_distance < 0.0) && (max_distance > 0.0)) return 0.0;
    }
    
    // Return signed distance
    if (min_distance > 0.0) return min_distance;
    else if (max_distance < 0.0) return max_distance;
    else return 0.0;
}



RNLength R3SignedDistance(const R3Plane& plane, const R3Halfspace& halfspace)
{
    // Return signed distance from plane to halfspace
    if (plane.Normal() != -(halfspace.Plane().Normal())) return 0.0;
    return (plane.D() + halfspace.Plane().D());
}



RNLength R3SignedDistance(const R3Plane& plane, const R3Box& box)
{
    // Check if box is empty
    if (box.IsEmpty()) return FLT_MAX;

    // Return signed distance from plane to box
    RNOctant octant = plane.Normal().Octant();
    RNScalar d1 = R3SignedDistance(plane, box.Corner(~octant & 0x7));
    if (RNIsPositiveOrZero(d1)) return d1;
    RNScalar d2 = R3SignedDistance(plane, box.Corner(octant));
    if (RNIsNegative(d2)) return d2;
    else return 0.0;
}



RNLength R3SignedDistance(const R3Plane& plane, const R3OrientedBox& box)
{
    // Check if box is empty
    if (box.IsEmpty()) return FLT_MAX;

    // Transform plane into std coordinate system
    R3Plane transformed_plane = plane;
    R3Affine cs_to_std_xform(box.CoordSystem().InverseMatrix());
    transformed_plane.Transform(cs_to_std_xform);

    // Return distance in std coordinate system
    R3Box axis_aligned_box(-box.Radius(0), -box.Radius(1), -box.Radius(2), box.Radius(0), box.Radius(1), box.Radius(2));
    return R3SignedDistance(transformed_plane, axis_aligned_box);
}



RNLength R3SignedDistance(const R3Plane& plane, const R3Sphere& sphere)
{
    RNLength d = R3SignedDistance(plane, sphere.Center());
    if (d < 0.0) {
	d += sphere.Radius();
	if (d > 0.0) return 0.0;
    }
    else if (d > 0.0) {
	d -= sphere.Radius();
	if (d < 0.0) return 0.0;
    }
    return d;
}



RNLength R3SignedDistance(const R3Plane& plane, const R3Shape& shape)
{
    // Return signed distance from plane to shape
    return shape.SignedDistance(plane);
}



RNLength R3Distance(const R3Halfspace& halfspace1, const R3Halfspace& halfspace2)
{
    // Return signed distance from plane to halfspace
    if (halfspace1.Plane().Normal() != -(halfspace2.Plane().Normal())) return 0.0;
    return (halfspace1.Plane().D() + halfspace2.Plane().D());
}



RNLength R3Distance(const R3Halfspace& halfspace, const R3Box& box)
{
    RNAbort("Not implemented");
    return 0.0;
}



RNLength R3Distance(const R3Halfspace& halfspace, const R3OrientedBox& box)
{
    RNAbort("Not implemented");
    return 0.0;
}



RNLength R3Distance(const R3Halfspace& halfspace, const R3Sphere& sphere)
{
    // Return distance from halfspace to sphere
    RNLength d = R3Distance(halfspace, sphere.Center()) - sphere.Radius();
    return ((d > 0.0) ? d : 0.0);
}



RNLength R3Distance(const R3Halfspace& halfspace, const R3Shape& shape)
{
    // Return distance from halfspace to shape
    return shape.Distance(halfspace);
}



RNLength R3Distance(const R3Box& box1, const R3Box& box2)
{
    // Check if either box is empty
    if (box1.IsEmpty()) return FLT_MAX;
    if (box2.IsEmpty()) return FLT_MAX;

    // Find axial distances from box1 to box2
    RNLength dx, dy, dz;
    if (RNIsGreater(box1.XMin(), box2.XMax())) dx = box1.XMin() - box2.XMax();
    else if (RNIsGreater(box2.XMin(), box1.XMax())) dx = box2.XMin() - box1.XMax();
    else dx = 0.0;
    if (RNIsGreater(box1.YMin(), box2.YMax())) dy = box1.YMin() - box2.YMax();
    else if (RNIsGreater(box2.YMin(), box1.YMax())) dy = box2.YMin() - box1.YMax();
    else dy = 0.0;
    if (RNIsGreater(box1.ZMin(), box2.ZMax())) dz = box1.ZMin() - box2.ZMax();
    else if (RNIsGreater(box2.ZMin(), box1.ZMax())) dz = box2.ZMin() - box1.ZMax();
    else dz = 0.0;
    
    // Return distance between point and closest point in box 
    if ((dy == 0.0) && (dz == 0.0)) return dx;
    else if ((dx == 0.0) && (dz == 0.0)) return dy;
    else if ((dx == 0.0) && (dy == 0.0)) return dz;
    else return sqrt(dx*dx + dy*dy + dz*dz);
}



RNLength R3Distance(const R3Box& box1, const R3OrientedBox& box2)
{
    RNAbort("Not implemented");
    return 0.0;
}



RNLength R3Distance(const R3Box& box, const R3Sphere& sphere)
{
    // Return distance from box to sphere
    RNLength d = R3Distance(box, sphere.Center()) - sphere.Radius();
    return ((d > 0.0) ? d : 0.0);
}



RNLength R3Distance(const R3Box& box, const R3Shape& shape)
{
    // Return distance from box to shape
    return shape.Distance(box);
}



RNLength R3Distance(const R3OrientedBox& box1, const R3OrientedBox& box2)
{
    RNAbort("Not implemented");
    return 0.0;
}



RNLength R3Distance(const R3OrientedBox& box, const R3Sphere& sphere)
{
    // Check if box is empty
    if (box.IsEmpty()) return FLT_MAX;

    // Transform sphere into std coordinate system
    R3Sphere transformed_sphere = sphere;
    R3Affine cs_to_std_xform(box.CoordSystem().InverseMatrix());
    transformed_sphere.Transform(cs_to_std_xform);

    // Return distance in std coordinate system
    R3Box axis_aligned_box(-box.Radius(0), -box.Radius(1), -box.Radius(2), box.Radius(0), box.Radius(1), box.Radius(2));
    return R3Distance(transformed_sphere, axis_aligned_box);
}



RNLength R3Distance(const R3OrientedBox& box, const R3Shape& shape)
{
    // Return distance from box to shape
    return shape.Distance(box);
}



RNLength R3Distance(const R3Sphere& sphere1, const R3Sphere& sphere2)
{
    // Return distance from sphere to sphere
    RNLength d = R3Distance(sphere2, sphere1.Center()) - sphere1.Radius();
    return ((d > 0.0) ? d : 0.0);
}



RNLength R3Distance(const R3Sphere& sphere, const R3Shape& shape)
{
    // Return distance from shape to sphere
    RNLength d = shape.Distance(sphere.Center()) - sphere.Radius();
    return ((d > 0.0) ? d : 0.0);
}



RNLength R3Distance(const R3Shape& shape1, const R3Shape& shape2)
{
    // Return distance between shapes
    return shape1.Distance(shape2);
}



#ifdef R3_USER_DEFINED_SHAPE_TYPES

RNLength R3Distance(const R3Shape& shape1, const R3Shape& shape2)
{
    RNClassID id1 = shape1.ClassID();
    assert(id1 < R3_NUM_SHAPE_TYPES);
    RNClassID id2 = shape2.ClassID();
    assert(id2 < R3_NUM_SHAPE_TYPES);
    return R3distance_functions[id1][id2](shape1, shape2);
}    

void R3RegisterDistanceFunction(RNClassID id1, RNClassID id2, R3DistanceFunction *distance)
{
    // Check ids
    assert(id1 <= R3_NUM_SHAPE_TYPES);
    assert(id2 <= R3_NUM_SHAPE_TYPES);

    // Add distance function to array
    R3distance_functions[id1][id2] = distance;
}

#endif




