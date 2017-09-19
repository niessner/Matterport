/* Source file for GAPS containment utility */



/* Include files */

#include "R3Shapes/R3Shapes.h"



RNBoolean R3Contains(const R3Vector& vector1, const R3Vector& vector2)
{
    // Return whether vector1 and vector2 are equal within tolerance
    return (RNIsEqual(vector1.X(), vector2.X()) &&
            RNIsEqual(vector1.Y(), vector2.Y()) &&
            RNIsEqual(vector1.Z(), vector2.Z()));
}



RNBoolean R3Contains(const R3Point& point1, const R3Point& point2)
{
    // Return whether point1 and point2 are equal to within tolerance
    return (RNIsEqual(point1.X(), point2.X()) &&
            RNIsEqual(point1.Y(), point2.Y()) &&
            RNIsEqual(point1.Z(), point2.Z()));
}



RNBoolean R3Contains(const R3Point& /*point*/, const R3Line& /*line*/)
{
    // Impossible
    return FALSE;
}



RNBoolean R3Contains(const R3Point& /*point*/, const R3Ray& /*ray*/)
{
    // Impossible
    return FALSE;
}



RNBoolean R3Contains(const R3Point& point, const R3Span& span)
{
    // Return whether point contains span's endpoints
    return (R3Contains(point, span.Start()) && 
	    R3Contains(point, span.End()));
}



RNBoolean R3Contains(const R3Point& point, const R3Triangle& triangle)
{
    // Get triangle vertices
    R3TriangleVertex *v0 = triangle.V0();
    R3TriangleVertex *v1 = triangle.V1();
    R3TriangleVertex *v2 = triangle.V2();
    if (!v0 || !v1 || !v2) return FALSE;

    // Return whether point contains triangle's vertices
    return (R3Contains(point, v0->Position()) &&
            R3Contains(point, v1->Position()) &&
            R3Contains(point, v2->Position())); 
}



RNBoolean R3Contains(const R3Point& point, const R3Circle& circle)
{
    // Check circle radius and center
    if (RNIsPositive(circle.Radius())) return FALSE;
    return R3Contains(point, circle.Center());
}



RNBoolean R3Contains(const R3Point& point, const R3Ellipse& ellipse)
{
    // Check ellipse radii and center
    if (RNIsPositive(ellipse.Radius(0)) || RNIsPositive(ellipse.Radius(1))) return FALSE;
    return R3Contains(point, ellipse.Center());
}



RNBoolean R3Contains(const R3Point& point, const R3Rectangle& rectangle)
{
    // Check rectangle radii and center
    if (RNIsPositive(rectangle.Radius(0)) || RNIsPositive(rectangle.Radius(1))) return FALSE;
    return R3Contains(point, rectangle.Center());
}



RNBoolean R3Contains(const R3Point& /*point*/, const R3Plane& /*plane*/)
{
    // Impossible
    return FALSE;
}



RNBoolean R3Contains(const R3Point& /*point*/, const R3Halfspace& /*halfspace*/)
{
    // Impossible
    return FALSE;
}



RNBoolean R3Contains(const R3Point& point, const R3Box& box)
{
    // Return whether point contains box
    if (box.IsEmpty()) return FALSE;
    if (!box.IsPoint()) return FALSE;
    if (!R3Contains(point, box.Min())) return FALSE;
    return TRUE;
}



RNBoolean R3Contains(const R3Point& point, const R3OrientedBox& box)
{
    // Return whether point contains box
    if (box.IsEmpty()) return FALSE;
    if (!box.IsPoint()) return FALSE;
    if (!R3Contains(point, box.Center())) return FALSE;
    return TRUE;
}



RNBoolean R3Contains(const R3Point& point, const R3Sphere& sphere)
{
    // Return whether point contains sphere
    if (!sphere.IsPoint()) return FALSE;
    if (!R3Contains(point, sphere.Center())) return FALSE;
    return TRUE;
}



RNBoolean R3Contains(const R3Point& point, const R3Shape& shape)
{
    // Return whether point contains shape
    return shape.Inside(point);
}



RNBoolean R3Contains(const R3Line& line, const R3Point& point)
{
    // Return whether line contains point ???
    return RNIsZero(R3Distance(line, point));
}



RNBoolean R3Contains(const R3Line& line1, const R3Line& line2)
{
    // Return whether two lines are equal within tolerance
    return (R3Contains(line1.Vector(), line2.Vector()) &&
	    R3Contains(line1, line2.Point()));
}



RNBoolean R3Contains(const R3Line& line, const R3Ray& ray)
{
    // Return whether line contains ray
    return (R3Parallel(line, ray) &&
	    R3Contains(line, ray.Start()));
}



RNBoolean R3Contains(const R3Line& line, const R3Span& span)
{
    // Return whether line contains span ???
    return (R3Contains(line, span.Start()) &&
	    R3Contains(line, span.End()));
}



RNBoolean R3Contains(const R3Line& /*line*/, const R3Plane& /*plane*/)
{
    // Impossible
    return FALSE;
}



RNBoolean R3Contains(const R3Line& /*line*/, const R3Halfspace& /*halfspace*/)
{
    // Impossible
    return FALSE;
}



RNBoolean R3Contains(const R3Line& line, const R3Box& box)
{
    // Return whether line contains box 
    if (box.IsEmpty()) return FALSE;
    if (!box.IsLinear()) return FALSE;
    if (!R3Contains(line, box.Min())) return FALSE;
    if (!R3Contains(line, box.Max())) return FALSE;
    return TRUE;
}



RNBoolean R3Contains(const R3Line& line, const R3OrientedBox& box)
{
    // Return whether line contains box 
    if (box.IsEmpty()) return FALSE;
    if (!box.IsLinear()) return FALSE;
    for (int i = 0; i < RN_NUM_OCTANTS; i++) 
      if (!R3Contains(line, box.Corner(i))) return FALSE;
    return TRUE;
}



RNBoolean R3Contains(const R3Line& line, const R3Sphere& sphere)
{
    // Return whether line contains sphere 
    if (!sphere.IsPoint()) return FALSE;
    if (!R3Contains(line, sphere.Center())) return FALSE;
    return TRUE;
}



RNBoolean R3Contains(const R3Line& line, const R3Shape& shape)
{
    // Return whether line contains shape ???
    return shape.Inside(line);
}



RNBoolean R3Contains(const R3Ray& ray, const R3Point& point)
{
    // Return whether ray contains point ???
    if (ray.IsZero()) return FALSE;
    else return RNIsZero(R3Distance(ray, point));
}



RNBoolean R3Contains(const R3Ray& /*ray*/, const R3Line& /*line*/)
{
    // Impossible
    return FALSE;
}



RNBoolean R3Contains(const R3Ray& ray1, const R3Ray& ray2)
{
    // Return whether ray1 contains ray2
    return (R3Contains(ray1.Vector(), ray2.Vector()) &&
	    R3Contains(ray1, ray2.Start()));
}



RNBoolean R3Contains(const R3Ray& ray, const R3Span& span)
{
    // Return whether ray contains span
    return (R3Contains(ray, span.Start()) &&
	    R3Contains(ray, span.End()));
}



RNBoolean R3Contains(const R3Ray& /*ray*/, const R3Plane& /*plane*/)
{
    // Impossible
    return FALSE;
}



RNBoolean R3Contains(const R3Ray& /*ray*/, const R3Halfspace& /*halfspace*/)
{
    // Impossible
    return FALSE;
}



RNBoolean R3Contains(const R3Ray& ray, const R3Box& box)
{
    // Return whether ray contains box
    if (!R3Contains(ray.Line(), box)) return FALSE;
    RNOctant octant = ray.Vector().Octant();
    if (!R3Contains(ray, box.Corner(~octant & 0x7))) return FALSE;
    return TRUE;
}



RNBoolean R3Contains(const R3Ray& ray, const R3OrientedBox& box)
{
    // Return whether ray contains box
    if (box.IsEmpty()) return FALSE;
    if (!box.IsLinear()) return FALSE;
    for (int i = 0; i < RN_NUM_OCTANTS; i++) 
      if (!R3Contains(ray, box.Corner(i))) return FALSE;
    return TRUE;
}




RNBoolean R3Contains(const R3Ray& ray, const R3Sphere& sphere)
{
    // Return whether ray contains sphere 
    if (!sphere.IsPoint()) return FALSE;
    if (!R3Contains(ray, sphere.Center())) return FALSE;
    return TRUE;
}



RNBoolean R3Contains(const R3Ray& ray, const R3Shape& shape)
{
    // Return whether ray contains shape
    return shape.Inside(ray);
}



RNBoolean R3Contains(const R3Span& span, const R3Point& point)
{
    // Return whether span contains point ???
    if (span.IsPoint()) return R3Contains(span.Start(), point);
    else return RNIsZero(R3Distance(span, point));
}



RNBoolean R3Contains(const R3Span& /*span*/, const R3Line& /*line*/)
{
    // Impossible
    return FALSE;
}



RNBoolean R3Contains(const R3Span& /*span*/, const R3Ray& /*ray*/)
{
    // Impossible
    return FALSE;
}



RNBoolean R3Contains(const R3Span& span1, const R3Span& span2)
{
    // Return whether span1 contains span2
    return (R3Contains(span1, span2.Start()) &&
	    R3Contains(span1, span2.End()));
}



RNBoolean R3Contains(const R3Span& /*span*/, const R3Plane& /*plane*/)
{
    // Impossible
    return FALSE;
}



RNBoolean R3Contains(const R3Span& /*span*/, const R3Halfspace& /*halfspace*/)
{
    // Impossible
    return FALSE;
}



RNBoolean R3Contains(const R3Span& span, const R3Box& box)
{
    // Return whether span contains box
    if (!R3Contains(span.Line(), box)) return FALSE;
    if (!R3Contains(span, box.Min())) return FALSE;
    if (!R3Contains(span, box.Max())) return FALSE;
    return TRUE;
}



RNBoolean R3Contains(const R3Span& span, const R3OrientedBox& box)
{
    // Return whether span contains box
    if (box.IsEmpty()) return FALSE;
    if (!box.IsLinear()) return FALSE;
    for (int i = 0; i < RN_NUM_OCTANTS; i++) 
      if (!R3Contains(span, box.Corner(i))) return FALSE;
    return TRUE;
}



RNBoolean R3Contains(const R3Span& span, const R3Sphere& sphere)
{
    // Return whether span contains sphere 
    if (!sphere.IsPoint()) return FALSE;
    if (!R3Contains(span, sphere.Center())) return FALSE;
    return TRUE;
}



RNBoolean R3Contains(const R3Span& span, const R3Shape& shape)
{
    // Return whether span contains shape
    return shape.Inside(span);
}



RNBoolean R3Contains(const R3Plane& plane, const R3Point& point)
{
    // Return whether plane contains point
    return RNIsZero(R3SignedDistance(plane, point));
}



RNBoolean R3Contains(const R3Plane& plane, const R3Line& line)
{
    // Return whether plane contains line
    return (R3Parallel(plane, line) && 
	    R3Contains(plane, line.Point()));
}



RNBoolean R3Contains(const R3Plane& plane, const R3Ray& ray)
{
    // Return whether plane contains ray
    return (R3Parallel(plane, ray) &&
	    R3Contains(plane, ray.Start()));
}



RNBoolean R3Contains(const R3Plane& plane, const R3Span& span)
{
    // Return whether plane contains span
    return (R3Contains(plane, span.Start()) &&
	    R3Contains(plane, span.End()));
}



RNBoolean R3Contains(const R3Plane& plane, const R3Triangle& triangle)
{
    // Return whether plane contains triangle
    return (R3Contains(plane, triangle.Plane()));
}



RNBoolean R3Contains(const R3Plane& plane1, const R3Plane& plane2)
{
    // Return whether plane1 and plane2 are equal within tolerance
    return (R3Contains(plane1.Normal(), plane2.Normal()) &&
	    RNIsEqual(plane1.D(), plane2.D()));
}



RNBoolean R3Contains(const R3Plane& /*plane*/, const R3Halfspace& /*halfspace*/)
{
    // Impossible
    return FALSE;
}



RNBoolean R3Contains(const R3Plane& plane, const R3Box& box)
{
    // Return whether plane contains box
    if (box.IsEmpty()) return FALSE;
    if (!box.IsPlanar()) return FALSE;
    if (!R3Contains(plane, box.Min())) return FALSE;
    if (!R3Contains(plane, box.Max())) return FALSE;
    return TRUE;
}



RNBoolean R3Contains(const R3Plane& plane, const R3OrientedBox& box)
{
    // Return whether plane contains box
    if (box.IsEmpty()) return FALSE;
    if (!box.IsPlanar()) return FALSE;
    for (int i = 0; i < RN_NUM_OCTANTS; i++) 
      if (!R3Contains(plane, box.Corner(i))) return FALSE;
    return TRUE;
}



RNBoolean R3Contains(const R3Plane& plane, const R3Sphere& sphere)
{
    // Return whether plane contains sphere 
    if (!sphere.IsPoint()) return FALSE;
    if (!R3Contains(plane, sphere.Center())) return FALSE;
    return TRUE;
}



RNBoolean R3Contains(const R3Plane& plane, const R3Shape& shape)
{
    // Return whether plane contains shape
    return shape.Inside(plane);
}



RNBoolean R3Contains(const R3Triangle& triangle, const R3Point& point)
{
    // Check whether triangle bounding shape contains point
    if (!R3Contains(triangle.Box(), point)) return FALSE;
  
    // Check whether triangle plane contains point
    if (!R3Contains(triangle.Plane(), point)) return FALSE;
  
    // Compute whether point is on correct side of each edge
    const R3Point& p0 = triangle.Vertex(0)->Position();
    const R3Point& p1 = triangle.Vertex(1)->Position();
    const R3Point& p2 = triangle.Vertex(2)->Position();
    R3Plane h01(p1, triangle.Normal(), p1 - p0);
    if (RNIsNegative(R3SignedDistance(h01, point))) return FALSE;
    R3Plane plane01(p1, triangle.Normal(), p1 - p0);
    if (RNIsNegative(R3SignedDistance(plane01, point))) return FALSE;
    R3Plane plane12(p2, triangle.Normal(), p2 - p1);
    if (RNIsNegative(R3SignedDistance(plane12, point))) return FALSE;
    R3Plane plane20(p0, triangle.Normal(), p0 - p2);
    if (RNIsNegative(R3SignedDistance(plane20, point))) return FALSE;
  
    // Triangle contains point
    return TRUE;
}



RNBoolean R3Contains(const R3Triangle& /* triangle */, const R3Line& /* line */)
{
    // Impossible 
    return FALSE;
}



RNBoolean R3Contains(const R3Triangle& /* triangle */, const R3Ray& /* ray */)
{
    // Impossible
    return FALSE;
}



RNBoolean R3Contains(const R3Triangle& triangle, const R3Span& span)
{
    // Return whether triangle contains span
    return (R3Contains(triangle, span.Start()) &&
	    R3Contains(triangle, span.End()));
}



RNBoolean R3Contains(const R3Triangle& /* triangle */, const R3Plane& /* plane */)
{
    // Impossible
    return FALSE;
}



RNBoolean R3Contains(const R3Triangle& triangle1, const R3Triangle& triangle2)
{
    // Check whether planes are the same
    if (!R3Contains(triangle1.Plane(), triangle2.Plane())) return FALSE;

    // Check whether triangle1's bounding box contains triangle2
    if (!R3Contains(triangle1.BBox(), triangle2)) return FALSE;

    // Check whether points of triangle2 are inside triangle1
    for (int i = 0; i < 3; i++) {
	if (!R3Contains(triangle1, triangle2.Vertex(i)->Position())) 
            return FALSE;
    }

    // Passed all tests
    return TRUE;
}



RNBoolean R3Contains(const R3Triangle& /*triangle*/, const R3Halfspace& /*halfspace*/)
{
    // Impossible
    return FALSE;
}



RNBoolean R3Contains(const R3Triangle& triangle, const R3Box& box)
{
    // Return whether triangle contains box
    if (box.IsEmpty()) return FALSE;
    if (!box.IsPlanar()) return FALSE;
    if (!R3Contains(triangle, box.Min())) return FALSE;
    if (!R3Contains(triangle, box.Max())) return FALSE;
    return TRUE;
}



RNBoolean R3Contains(const R3Triangle& triangle, const R3OrientedBox& box)
{
    // Return whether triangle contains box
    if (box.IsEmpty()) return FALSE;
    if (!box.IsPlanar()) return FALSE;
    for (int i = 0; i < RN_NUM_OCTANTS; i++) 
      if (!R3Contains(triangle, box.Corner(i))) return FALSE;
    return TRUE;
}



RNBoolean R3Contains(const R3Triangle& triangle, const R3Sphere& sphere)
{
    // Return whether triangle contains sphere 
    if (!sphere.IsPoint()) return FALSE;
    if (!R3Contains(triangle, sphere.Center())) return FALSE;
    return TRUE;
}



RNBoolean R3Contains(const R3Triangle& triangle, const R3Shape& shape)
{
    // Return whether triangle contains shape
    // return shape.Inside(triangle);
    RNAbort("Not implemented");
    return FALSE;
}



RNBoolean R3Contains(const R3Circle& circle, const R3Point& point)
{
    // Check whether ellipse is empty
    if (circle.IsEmpty()) return FALSE;

    // Check whether circle plane contains point
    if (!R3Contains(circle.Plane(), point)) return FALSE;

    // Check whether point is within radius of center
    RNScalar dd = R3SquaredDistance(circle.Center(), point);
    if (dd > circle.Radius() * circle.Radius()) return FALSE;

    // Circle contains point
    return TRUE;
}



RNBoolean R3Contains(const R3Ellipse& ellipse, const R3Point& point)
{
    // Check whether ellipse is empty
    if (ellipse.IsEmpty()) return FALSE;
    if (!RNIsPositive(ellipse.Radius(0))) return FALSE;
    if (!RNIsPositive(ellipse.Radius(1))) return FALSE;

    // Check whether ellipse plane contains point
    if (!R3Contains(ellipse.Plane(), point)) return FALSE;

    // Check whether point is within radius of center
    R3Vector v = point - ellipse.Center();
    RNScalar d0 = v.Dot(ellipse.Axis(0)) / ellipse.Radius(0);
    RNScalar d1 = v.Dot(ellipse.Axis(1)) / ellipse.Radius(1);
    RNScalar dd = d0*d0 + d1*d1;
    if (dd > 1.0) return FALSE;

    // Ellipse contains point
    return TRUE;
}



RNBoolean R3Contains(const R3Rectangle& rectangle, const R3Point& point)
{
    // Check whether rectangle is empty
    if (rectangle.IsEmpty()) return FALSE;
    if (!RNIsPositive(rectangle.Radius(0))) return FALSE;
    if (!RNIsPositive(rectangle.Radius(1))) return FALSE;

    // Check whether rectangle plane contains point
    if (!R3Contains(rectangle.Plane(), point)) return FALSE;

    // Check whether point is within radius of center
    R3Vector v = point - rectangle.Center();
    RNScalar d0 = fabs(v.Dot(rectangle.Axis(0)));
    if (d0 > rectangle.Radius(0)) return FALSE;
    RNScalar d1 = fabs(v.Dot(rectangle.Axis(1)));
    if (d1 > rectangle.Radius(1)) return FALSE;

    // Rectangle contains point
    return TRUE;
}



RNBoolean R3Contains(const R3Halfspace& halfspace, const R3Point& point)
{
    // Return whether halfspace contains point
    return (RNIsPositiveOrZero(R3SignedDistance(halfspace.Plane(), point)));
}



RNBoolean R3Contains(const R3Halfspace& halfspace, const R3Line& line)
{
    // Return whether halfspace contains line
    return (R3Parallel(halfspace.Plane(), line) &&
	    R3Contains(halfspace, line.Point()));
}



RNBoolean R3Contains(const R3Halfspace& halfspace, const R3Ray& ray)
{
    // Return whether halfspace contains ray
    return (RNIsPositiveOrZero(halfspace.Normal().Dot(ray.Vector())) &&
	    R3Contains(halfspace, ray.Start()));
}



RNBoolean R3Contains(const R3Halfspace& halfspace, const R3Span& span)
{
    // Return whether halfspace contains span
    return (R3Contains(halfspace, span.Start()) &&
	    R3Contains(halfspace, span.End()));
}



RNBoolean R3Contains(const R3Halfspace& halfspace, const R3Plane& plane)
{
    // Check whether planes are parallel
    if (!R3Parallel(halfspace.Plane(), plane)) return FALSE;

    // Check whether halfspace contains at least one point on plane
    if (!R3Contains(halfspace, plane.Point())) return FALSE;

    // Passed all tests
    return TRUE;
}



RNBoolean R3Contains(const R3Halfspace& halfspace, const R3Triangle& triangle)
{
    // Check whether halfspace contains triangle bounding shape
    if (R3Contains(halfspace, triangle.Box())) return TRUE;

    // Check whether halfspace contains all vertices of triangle
    if (!R3Contains(halfspace, triangle.Vertex(0)->Position())) return FALSE;
    if (!R3Contains(halfspace, triangle.Vertex(1)->Position())) return FALSE;
    if (!R3Contains(halfspace, triangle.Vertex(2)->Position())) return FALSE;

    // Passed all tests
    return TRUE;
}



RNBoolean R3Contains(const R3Halfspace& halfspace, const R3Circle& circle)
{
    // Return whether halfspace contains circle
    RNScalar d = R3SignedDistance(halfspace.Plane(), circle.Center());
    if (RNIsNegative(d)) return FALSE;
    else if (RNIsGreater(d, circle.Radius())) return TRUE;
    else {
        RNScalar cos_theta = halfspace.Plane().Normal().Dot(circle.Normal());
	if (cos_theta < 0.0) cos_theta = halfspace.Plane().Normal().Dot(-circle.Normal());
	return (RNIsGreaterOrEqual(d, cos_theta * circle.Radius()));
    }
}



RNBoolean R3Contains(const R3Halfspace& halfspace1, const R3Halfspace& halfspace2)
{
    // Return whether halfspace1 contains halfspace2
    RNScalar dot = halfspace1.Plane().Normal().Dot(halfspace2.Plane().Normal());
    if (RNIsEqual(dot, 1.0)) {
	// Halfspace and plane are parallel ???
	return RNIsLessOrEqual(halfspace1.Plane().D(), halfspace2.Plane().D());
    }
    else {
	// Halfspace and plane are not parallel
	return FALSE;
    }
}



RNBoolean R3Contains(const R3Halfspace& halfspace, const R3Box& box)
{
    // Return whether halfspace contains box 
    if (box.IsEmpty()) return FALSE;
    RNOctant octant = halfspace.Normal().Octant();
    return (R3Contains(halfspace, box.Corner(~octant & 0x7)));
}



RNBoolean R3Contains(const R3Halfspace& halfspace, const R3OrientedBox& box)
{
    // Return whether halfspace contains box
    if (box.IsEmpty()) return FALSE;
    if (!box.IsPlanar()) return FALSE;
    for (int i = 0; i < RN_NUM_OCTANTS; i++) 
      if (!R3Contains(halfspace, box.Corner(i))) return FALSE;
    return TRUE;
}



RNBoolean R3Contains(const R3Halfspace& halfspace, const R3Sphere& sphere)
{
    // Return whether halfspace contains sphere 
    RNLength d = R3SignedDistance(halfspace.Plane(), sphere.Center());
    return RNIsGreaterOrEqual(d, sphere.Radius());
}



RNBoolean R3Contains(const R3Halfspace& halfspace, const R3Cylinder& cylinder)
{
    // Return whether halfspace contains cylinder 
    return (R3Contains(halfspace, cylinder.Top()) &&
	    R3Contains(halfspace, cylinder.Base()));
}



RNBoolean R3Contains(const R3Halfspace& halfspace, const R3Cone& cone)
{
    // Return whether halfspace contains cone 
    return (R3Contains(halfspace, cone.Apex()) &&
	    R3Contains(halfspace, cone.Base()));
}



RNBoolean R3Contains(const R3Halfspace& halfspace, const R3Shape& shape)
{
    // Return whether halfspace contains shape ???
    return shape.Inside(halfspace);
}



RNBoolean R3Contains(const R3Box& box, const R3Point& point)
{
    // Return whether box contains point
    if (box.IsEmpty()) return FALSE;
    if (RNIsLess(point.X(), box.XMin())) return FALSE;
    if (RNIsLess(point.Y(), box.YMin())) return FALSE;
    if (RNIsLess(point.Z(), box.ZMin())) return FALSE;
    if (RNIsGreater(point.X(), box.XMax())) return FALSE;
    if (RNIsGreater(point.Y(), box.YMax())) return FALSE;
    if (RNIsGreater(point.Z(), box.ZMax())) return FALSE;
    return TRUE;
}



RNBoolean R3Contains(const R3Box& box, const R3Line& line)
{
    // Return whether box contains line
    if (box.IsEmpty()) return FALSE;
    if (box.IsFinite()) return FALSE;
    RNAbort("Not Implemented");
    return FALSE;
}



RNBoolean R3Contains(const R3Box& box, const R3Ray& ray)
{
    // Return whether box contains ray
    if (box.IsEmpty()) return FALSE;
    if (box.IsFinite()) return FALSE;
    RNAbort("Not Implemented");
    return FALSE;
}



RNBoolean R3Contains(const R3Box& box, const R3Span& span)
{
    // Return whether box contains span
    if (!R3Contains(box, span.Start())) return FALSE;
    if (!R3Contains(box, span.End())) return FALSE;
    return TRUE;
}



RNBoolean R3Contains(const R3Box& box, const R3Plane& plane)
{
    // Return whether box contains plane
    if (box.IsEmpty()) return FALSE;
    if (box.IsFinite()) return FALSE;
    RNAbort("Not Implemented");
    return FALSE;
}



RNBoolean R3Contains(const R3Box& box, const R3Halfspace& halfspace)
{
    // Return whether box contains halfspace
    if (box.IsEmpty()) return FALSE;
    RNOctant octant = halfspace.Normal().Octant();
    R3Point corner = box.Corner(octant);
    if (!RNIsFinite(corner.X())) return FALSE;
    if (!RNIsFinite(corner.Y())) return FALSE;
    if (!RNIsFinite(corner.Z())) return FALSE;
    RNAbort("Not Implemented");
    return FALSE;
}



RNBoolean R3Contains(const R3Box& box1, const R3Box& box2) 
{
    // Return whether box1 contains box2
    if (box1.IsEmpty()) return FALSE;
    if (box2.IsEmpty()) return FALSE;
    if (RNIsLess(box2.XMin(), box1.XMin())) return FALSE;
    if (RNIsLess(box2.YMin(), box1.YMin())) return FALSE;
    if (RNIsLess(box2.ZMin(), box1.ZMin())) return FALSE;
    if (RNIsGreater(box2.XMax(), box1.XMax())) return FALSE;
    if (RNIsGreater(box2.YMax(), box1.YMax())) return FALSE;
    if (RNIsGreater(box2.ZMax(), box1.ZMax())) return FALSE;
    return TRUE;
}




RNBoolean R3Contains(const R3Box& box1, const R3OrientedBox& box2)
{
    // Return whether halfspace contains box
    if (box1.IsEmpty()) return FALSE;
    if (box2.IsEmpty()) return FALSE;
    for (int i = 0; i < RN_NUM_OCTANTS; i++) 
      if (!R3Contains(box1, box2.Corner(i))) return FALSE;
    return TRUE;
}



RNBoolean R3Contains(const R3Box& box, const R3Sphere& sphere)
{
    // Return whether box contains sphere 
    if (box.IsEmpty()) return FALSE;
    if (sphere.IsEmpty()) return TRUE;
    if (RNIsLess(sphere.Center().X() - sphere.Radius(), box.XMin())) return FALSE;
    if (RNIsLess(sphere.Center().Y() - sphere.Radius(), box.YMin())) return FALSE;
    if (RNIsLess(sphere.Center().Z() - sphere.Radius(), box.ZMin())) return FALSE;
    if (RNIsGreater(sphere.Center().X() + sphere.Radius(), box.XMax())) return FALSE;
    if (RNIsGreater(sphere.Center().Y() + sphere.Radius(), box.YMax())) return FALSE;
    if (RNIsGreater(sphere.Center().Z() + sphere.Radius(), box.ZMax())) return FALSE;
    return TRUE;
}



RNBoolean R3Contains(const R3Box& box, const R3Shape& shape)
{
    // Return whether box contains shape
    return shape.Inside(box);
}



RNBoolean R3Contains(const R3OrientedBox& box, const R3Point& point)
{
    // Check if box is empty
    if (box.IsEmpty()) return FALSE;

    // Return whether each halfspace contains point
    for (int dir = 0; dir < 2; dir++) {
      for (int dim = 0; dim < 3; dim++) {
        R3Halfspace halfspace(-box.Plane(dir, dim), 0);
        if (!R3Contains(halfspace, point)) {
          return FALSE;
        }
      }
    }

    // Passed all tests
    return TRUE;
}



RNBoolean R3Contains(const R3OrientedBox& box, const R3Line& line)
{
    // Return whether box contains line
    if (box.IsFinite()) return FALSE;
    return TRUE;
}



RNBoolean R3Contains(const R3OrientedBox& box, const R3Ray& ray)
{
    // Return whether box contains ray
    if (box.IsFinite()) return FALSE;
    return TRUE;
}



RNBoolean R3Contains(const R3OrientedBox& box, const R3Span& span)
{
    // Return whether box contains span
    if (!R3Contains(box, span.Start())) return FALSE;
    if (!R3Contains(box, span.End())) return FALSE;
    return TRUE;
}



RNBoolean R3Contains(const R3OrientedBox& box, const R3Plane& plane)
{
    // Return whether box contains plane
    if (box.IsFinite()) return FALSE;
    return TRUE;
}



RNBoolean R3Contains(const R3OrientedBox& box, const R3Halfspace& halfspace)
{
    // Return whether box contains halfspace
    if (box.IsFinite()) return FALSE;
    return TRUE;
}



RNBoolean R3Contains(const R3OrientedBox& box1, const R3Box& box2)
{
    // Check if box is empty
    if (box1.IsEmpty()) return FALSE;
    if (box2.IsEmpty()) return FALSE;

    // Return whether each halfspace contains point
    for (int dir = 0; dir < 2; dir++) {
      for (int dim = 0; dim < 3; dim++) {
        R3Halfspace halfspace(-box1.Plane(dir, dim), 0);
        if (!R3Contains(halfspace, box2)) {
          return FALSE;
        }
      }
    }

    // Passed all tests
    return TRUE;
}



RNBoolean R3Contains(const R3OrientedBox& box1, const R3OrientedBox& box2) 
{
    // Check if box is empty
    if (box1.IsEmpty()) return FALSE;
    if (box2.IsEmpty()) return FALSE;

    // Return whether each halfspace contains point
    for (int dir = 0; dir < 2; dir++) {
      for (int dim = 0; dim < 3; dim++) {
        R3Halfspace halfspace(-box1.Plane(dir, dim), 0);
        if (!R3Contains(halfspace, box2)) {
          return FALSE;
        }
      }
    }

    // Passed all tests
    return TRUE;
}




RNBoolean R3Contains(const R3OrientedBox& box, const R3Sphere& sphere)
{
    // Check if box is empty
    if (box.IsEmpty()) return FALSE;

    // Return whether each halfspace contains sphere
    for (int dir = 0; dir < 2; dir++) {
      for (int dim = 0; dim < 3; dim++) {
        R3Halfspace halfspace(-box.Plane(dir, dim), 0);
        if (!R3Contains(halfspace, sphere)) {
          return FALSE;
        }
      }
    }

    // Passed all tests
    return TRUE;
}



RNBoolean R3Contains(const R3OrientedBox& box, const R3Shape& shape)
{
    // Check if box is empty
    if (box.IsEmpty()) return FALSE;

    // Return whether each halfspace contains shape
    for (int dir = 0; dir < 2; dir++) {
      for (int dim = 0; dim < 3; dim++) {
        R3Halfspace halfspace(-box.Plane(dir, dim), 0);
        if (!R3Contains(halfspace, shape)) {
          return FALSE;
        }
      }
    }

    // Passed all tests
    return TRUE;
}



RNBoolean R3Contains(const R3Sphere& sphere, const R3Point& point)
{
    // Return whether sphere contains point
    RNScalar radius_squared = sphere.Radius() * sphere.Radius();
    R3Vector v = sphere.Center() - point;
    RNScalar distance_squared = v.X() * v.X() + v.Y() * v.Y() + v.Z() * v.Z();
    return RNIsLessOrEqual(distance_squared, radius_squared);
}



RNBoolean R3Contains(const R3Sphere& sphere, const R3Line& line)
{
    // Return whether sphere contains line
    if (sphere.IsFinite()) return FALSE;
    return TRUE;
}



RNBoolean R3Contains(const R3Sphere& sphere, const R3Ray& ray)
{
    // Return whether sphere contains ray
    if (sphere.IsFinite()) return FALSE;
    return TRUE;
}



RNBoolean R3Contains(const R3Sphere& sphere, const R3Span& span)
{
    // Return whether sphere contains span
    if (!R3Contains(sphere, span.Start())) return FALSE;
    if (!R3Contains(sphere, span.End())) return FALSE;
    return TRUE;
}



RNBoolean R3Contains(const R3Sphere& sphere, const R3Plane& plane)
{
    // Return whether sphere contains plane
    if (sphere.IsFinite()) return FALSE;
    return TRUE;
}



RNBoolean R3Contains(const R3Sphere& sphere, const R3Halfspace& halfspace)
{
    // Return whether sphere contains halfspace
    if (sphere.IsFinite()) return FALSE;
    return TRUE;
}



RNBoolean R3Contains(const R3Sphere& sphere, const R3Box& box)
{
    // Return whether sphere contains box 
    if (box.IsEmpty()) return FALSE;
    R3Vector v = box.Centroid() - sphere.Center();
    R3Point corner = box.Corner(v.Octant());
    return R3Contains(sphere, corner);
}



RNBoolean R3Contains(const R3Sphere& sphere, const R3OrientedBox& box)
{
    // Return whether sphere contains box 
    if (box.IsEmpty()) return FALSE;
    for (int i = 0; i < RN_NUM_OCTANTS; i++) 
      if (!R3Contains(sphere, box.Corner(i))) return FALSE;
    return TRUE;
}



RNBoolean R3Contains(const R3Sphere& sphere1, const R3Sphere& sphere2) 
{
    // Return whether sphere1 contains sphere2
    RNLength d = R3Distance(sphere1.Center(), sphere2.Center());
    return RNIsLess(d + sphere2.Radius(), sphere1.Radius());
}



RNBoolean R3Contains(const R3Sphere& sphere, const R3Shape& shape)
{
    // Return whether sphere contains shape
    return shape.Inside(sphere);
}



RNBoolean R3Contains(const R3Cylinder& cylinder, const R3Point& point) 
{
    // Get cylinder axis
    const R3Span& axis = cylinder.Axis();
  
    // Check if outside top or bottom
    RNScalar t = axis.T(point);
    if (RNIsNegative(t) || RNIsGreater(t, axis.Length())) 
        return FALSE;

    // Check if inside radius
    RNLength d = R3Distance(point, axis.Point(t));
    return (RNIsLessOrEqual(d, cylinder.Radius()));
}



RNBoolean R3Contains(const R3Cone& cone, const R3Point& point) 
{
    // Get cone axis
    const R3Span& axis = cone.Axis();
  
    // Check if outside top or bottom
    RNScalar t = axis.T(point);
    if (RNIsNegative(t) || RNIsGreater(t, axis.Length())) 
        return FALSE;

    // Check if inside radius
    RNLength r = (1.0 - t / axis.Length()) * cone.Radius();
    RNLength d = R3Distance(axis.Point(t), point);
    return RNIsGreaterOrEqual(r, d);
}



RNBoolean R3Contains(const R3Shape& shape, const R3Point& point)
{
    // Return whether shape contains point
    return shape.Contains(point);
}



RNBoolean R3Contains(const R3Shape& shape, const R3Line& line)
{
    // Return whether shape contains line
    return shape.Contains(line);
}



RNBoolean R3Contains(const R3Shape& shape, const R3Ray& ray)
{
    // Return whether shape contains ray
    return shape.Contains(ray);
}



RNBoolean R3Contains(const R3Shape& shape, const R3Span& span)
{
    // Return whether shape contains span
    return shape.Contains(span);
}



RNBoolean R3Contains(const R3Shape& shape, const R3Plane& plane)
{
    // Return whether shape contains plane
    return shape.Contains(plane);
}



RNBoolean R3Contains(const R3Shape& shape, const R3Halfspace& halfspace)
{
    // Return whether shape contains halfspace
    return shape.Contains(halfspace);
}



RNBoolean R3Contains(const R3Shape& shape, const R3Box& box)
{
    // Return whether shape contains box
    return shape.Contains(box);
}



RNBoolean R3Contains(const R3Shape& shape, const R3OrientedBox& box)
{
    // Return whether shape contains sphere
    return shape.Contains(box);
}



RNBoolean R3Contains(const R3Shape& shape, const R3Sphere& sphere)
{
    // Return whether shape contains sphere
    return shape.Contains(sphere);
}



RNBoolean R3Contains(const R3Shape& shape1, const R3Shape& shape2)
{
    // Return whether shape1 contains shape2
    return shape1.Contains(shape2);
}





