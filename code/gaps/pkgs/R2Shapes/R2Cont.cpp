/* Source file for GAPS containment utility */



/* Include files */

#include "R2Shapes/R2Shapes.h"



RNBoolean R2Contains(const R2Vector& vector1, const R2Vector& vector2)
{
    // Return whether vector1 and vector2 are equal within tolerance
    if (!RNIsEqual(vector1.X(), vector2.X())) return FALSE;
    if (!RNIsEqual(vector1.Y(), vector2.Y())) return FALSE;
    return TRUE;
}



RNBoolean R2Contains(const R2Point& point1, const R2Point& point2)
{
    // Return whether point1 and point2 are equal to within tolerance
    if (!RNIsEqual(point1.X(), point2.X())) return FALSE;
    if (!RNIsEqual(point1.Y(), point2.Y())) return FALSE;
    return TRUE;
}



RNBoolean R2Contains(const R2Point& /*point*/, const R2Line& /*line*/)
{
    // Impossible
    return FALSE;
}



RNBoolean R2Contains(const R2Point& /*point*/, const R2Ray& /*ray*/)
{
    // Impossible
    return FALSE;
}



RNBoolean R2Contains(const R2Point& point, const R2Span& span)
{
    // Return whether point contains span's endpoints
    return (R2Contains(point, span.Start()) && 
	    R2Contains(point, span.End()));
}



RNBoolean R2Contains(const R2Point& /*point*/, const R2Halfspace& /*halfspace*/)
{
    // Impossible
    return FALSE;
}



RNBoolean R2Contains(const R2Point& point, const R2Box& box)
{
    // Return whether point contains box
    if (!box.IsPoint()) return FALSE;
    if (!R2Contains(point, box.Min())) return FALSE;
    return TRUE;
}



RNBoolean R2Contains(const R2Point& point, const R2Circle& circle)
{
    // Return whether point contains circle
    if (!circle.IsPoint()) return FALSE;
    if (!R2Contains(point, circle.Center())) return FALSE;
    return TRUE;
}



RNBoolean R2Contains(const R2Point& point, const R2Shape& shape)
{
    // Return whether point contains shape
    RNAbort("Not implemented");
    return FALSE;
}



RNBoolean R2Contains(const R2Line& line, const R2Point& point)
{
    // Return whether line contains point
    return RNIsZero(R2SignedDistance(line, point));
}



RNBoolean R2Contains(const R2Line& line1, const R2Line& line2)
{
    // Return whether two lines are equal within tolerance
    return (R2Contains(line1.Vector(), line2.Vector()) &&
	    RNIsEqual(line1.C(), line2.C()));
}



RNBoolean R2Contains(const R2Line& line, const R2Ray& ray)
{
    // Return whether line contains ray
    return (R2Parallel(line, ray) &&
	    R2Contains(line, ray.Start()));
}



RNBoolean R2Contains(const R2Line& line, const R2Span& span)
{
    // Return whether line contains span ???
    return (R2Contains(line, span.Start()) &&
	    R2Contains(line, span.End()));
}



RNBoolean R2Contains(const R2Line& /*line*/, const R2Halfspace& /*halfspace*/)
{
    // Impossible
    return FALSE;
}



RNBoolean R2Contains(const R2Line& line, const R2Box& box)
{
    // Return whether line contains box 
    if (!box.IsLinear()) return FALSE;
    if (!R2Contains(line, box.Min())) return FALSE;
    if (!R2Contains(line, box.Max())) return FALSE;
    return TRUE;
}



RNBoolean R2Contains(const R2Line& line, const R2Circle& circle)
{
    // Return whether line contains circle 
    if (!circle.IsPoint()) return FALSE;
    if (!R2Contains(line, circle.Center())) return FALSE;
    return TRUE;
}



RNBoolean R2Contains(const R2Line& line, const R2Shape& shape)
{
    // Return whether line contains shape ???
    RNAbort("Not implemented");
    return FALSE;
}



RNBoolean R2Contains(const R2Ray& ray, const R2Point& point)
{
    // Return whether ray contains point ???
    if (ray.IsZero()) return FALSE;
    else return RNIsZero(R2Distance(ray, point));
}



RNBoolean R2Contains(const R2Ray& /*ray*/, const R2Line& /*line*/)
{
    // Impossible
    return FALSE;
}



RNBoolean R2Contains(const R2Ray& ray1, const R2Ray& ray2)
{
    // Return whether ray1 contains ray2
    return (R2Contains(ray1.Vector(), ray2.Vector()) &&
	    R2Contains(ray1, ray2.Start()));
}



RNBoolean R2Contains(const R2Ray& ray, const R2Span& span)
{
    // Return whether ray contains span
    return (R2Contains(ray, span.Start()) &&
	    R2Contains(ray, span.End()));
}



RNBoolean R2Contains(const R2Ray& /*ray*/, const R2Halfspace& /*halfspace*/)
{
    // Impossible
    return FALSE;
}



RNBoolean R2Contains(const R2Ray& ray, const R2Box& box)
{
    // Return whether ray contains box
    if (!R2Contains(ray.Line(), box)) return FALSE;
    RNQuadrant quadrant = ray.Vector().Quadrant();
    if (!R2Contains(ray, box.Corner(~quadrant & 0x3))) return FALSE;
    return TRUE;
}



RNBoolean R2Contains(const R2Ray& ray, const R2Circle& circle)
{
    // Return whether ray contains circle 
    if (!circle.IsPoint()) return FALSE;
    if (!R2Contains(ray, circle.Center())) return FALSE;
    return TRUE;
}



RNBoolean R2Contains(const R2Ray& ray, const R2Shape& shape)
{
    // Return whether ray contains shape
    RNAbort("Not implemented");
    return FALSE;
}



RNBoolean R2Contains(const R2Span& span, const R2Point& point)
{
    // Return whether span contains point ???
    if (span.IsPoint()) return R2Contains(span.Start(), point);
    else return RNIsZero(R2Distance(span, point));
}



RNBoolean R2Contains(const R2Span& /*span*/, const R2Line& /*line*/)
{
    // Impossible
    return FALSE;
}



RNBoolean R2Contains(const R2Span& /*span*/, const R2Ray& /*ray*/)
{
    // Impossible
    return FALSE;
}



RNBoolean R2Contains(const R2Span& span1, const R2Span& span2)
{
    // Return whether span1 contains span2
    return (R2Contains(span1, span2.Start()) &&
	    R2Contains(span1, span2.End()));
}



RNBoolean R2Contains(const R2Span& /*span*/, const R2Halfspace& /*halfspace*/)
{
    // Impossible
    return FALSE;
}



RNBoolean R2Contains(const R2Span& span, const R2Box& box)
{
    // Return whether span contains box
    if (!R2Contains(span.Line(), box)) return FALSE;
    if (!R2Contains(span, box.Min())) return FALSE;
    if (!R2Contains(span, box.Max())) return FALSE;
    return TRUE;
}



RNBoolean R2Contains(const R2Span& span, const R2Circle& circle)
{
    // Return whether span contains circle 
    if (!circle.IsPoint()) return FALSE;
    if (!R2Contains(span, circle.Center())) return FALSE;
    return TRUE;
}



RNBoolean R2Contains(const R2Span& span, const R2Shape& shape)
{
    // Return whether span contains shape
    RNAbort("Not implemented");
    return FALSE;
}



RNBoolean R2Contains(const R2Halfspace& halfspace, const R2Point& point)
{
    // Return whether halfspace contains point
    return (!RNIsNegative(R2SignedDistance(halfspace.Line(), point)));
}



RNBoolean R2Contains(const R2Halfspace& halfspace, const R2Line& line)
{
    // Return whether halfspace contains line
    return (R2Parallel(halfspace.Line(), line) &&
	    R2Contains(halfspace, line.Point()));
}



RNBoolean R2Contains(const R2Halfspace& halfspace, const R2Ray& ray)
{
    // Return whether halfspace contains ray
    return (RNIsPositiveOrZero(halfspace.Normal().Dot(ray.Vector())) &&
	    R2Contains(halfspace, ray.Start()));
}



RNBoolean R2Contains(const R2Halfspace& halfspace, const R2Span& span)
{
    // Return whether halfspace contains span
    return (R2Contains(halfspace, span.Start()) &&
	    R2Contains(halfspace, span.End()));
}



RNBoolean R2Contains(const R2Halfspace& halfspace1, const R2Halfspace& halfspace2)
{
    // Return whether halfspace1 contains halfspace2
    RNScalar dot = halfspace1.Normal().Dot(halfspace2.Normal());
    if (RNIsEqual(dot, 1.0)) {
	// Halfspaces are parallel ???
	return RNIsLessOrEqual(halfspace1.Line().C(), halfspace2.Line().C());
    }
    else {
	// Halfspaces are not parallel
	return FALSE;
    }
}



RNBoolean R2Contains(const R2Halfspace& halfspace, const R2Box& box)
{
    // Return whether halfspace contains box 
    RNQuadrant quadrant = halfspace.Normal().Quadrant();
    return (R2Contains(halfspace, box.Corner(~quadrant & 0x3)));
}



RNBoolean R2Contains(const R2Halfspace& halfspace, const R2Circle& circle)
{
    // Return whether halfspace contains circle 
    RNLength d = R2SignedDistance(halfspace.Line(), circle.Center());
    return RNIsGreaterOrEqual(d, circle.Radius());
}



RNBoolean R2Contains(const R2Halfspace& halfspace, const R2Shape& shape)
{
    // Return whether halfspace contains shape ???
    RNAbort("Not implemented");
    return FALSE;
}



RNBoolean R2Contains(const R2Box& box, const R2Point& point)
{
    // Return whether box contains point
    if (box.IsEmpty()) return FALSE;
    if (RNIsLess(point.X(), box.XMin())) return FALSE;
    if (RNIsLess(point.Y(), box.YMin())) return FALSE;
    if (RNIsGreater(point.X(), box.XMax())) return FALSE;
    if (RNIsGreater(point.Y(), box.YMax())) return FALSE;
    return TRUE;
}



RNBoolean R2Contains(const R2Box& box, const R2Line& line)
{
    // Return whether box contains line
    if (box.IsFinite()) return FALSE;
    RNAbort("Not Implemented");
    return FALSE;
}



RNBoolean R2Contains(const R2Box& box, const R2Ray& ray)
{
    // Return whether box contains ray
    if (box.IsFinite()) return FALSE;
    RNAbort("Not Implemented");
    return FALSE;
}



RNBoolean R2Contains(const R2Box& box, const R2Span& span)
{
    // Return whether box contains span
    if (!R2Contains(box, span.Start())) return FALSE;
    if (!R2Contains(box, span.End())) return FALSE;
    return TRUE;
}



RNBoolean R2Contains(const R2Box& box, const R2Halfspace& halfspace)
{
    // Return whether box contains halfspace
    RNQuadrant quadrant = halfspace.Normal().Quadrant();
    R2Point corner = box.Corner(quadrant);
    if (!RNIsFinite(corner.X())) return FALSE;
    if (!RNIsFinite(corner.Y())) return FALSE;
    RNAbort("Not Implemented");
    return FALSE;
}



RNBoolean R2Contains(const R2Box& box1, const R2Box& box2) 
{
    // Return whether box1 contains box2
    if (box1.IsEmpty()) return FALSE;
    if (box2.IsEmpty()) return TRUE;
    if (RNIsLess(box2.XMin(), box1.XMin())) return FALSE;
    if (RNIsLess(box2.YMin(), box1.YMin())) return FALSE;
    if (RNIsGreater(box2.XMax(), box1.XMax())) return FALSE;
    if (RNIsGreater(box2.YMax(), box1.YMax())) return FALSE;
    return TRUE;
}




RNBoolean R2Contains(const R2Box& box, const R2Circle& circle)
{
    // Return whether box contains circle 
    if (box.IsEmpty()) return FALSE;
    if (circle.IsEmpty()) return TRUE;
    if (RNIsLess(circle.Center().X() - circle.Radius(), box.XMin())) return FALSE;
    if (RNIsLess(circle.Center().Y() - circle.Radius(), box.YMin())) return FALSE;
    if (RNIsGreater(circle.Center().X() + circle.Radius(), box.XMax())) return FALSE;
    if (RNIsGreater(circle.Center().Y() + circle.Radius(), box.YMax())) return FALSE;
    return TRUE;
}



RNBoolean R2Contains(const R2Box& box, const R2Shape& shape)
{
    // Return whether box contains shape
    RNAbort("Not implemented");
    return FALSE;
}



RNBoolean R2Contains(const R2Circle& circle, const R2Point& point)
{
    // Return whether circle contains point
    RNScalar radius_squared = circle.Radius() * circle.Radius();
    R2Vector v = circle.Center() - point;
    RNScalar distance_squared = v.X() * v.X() + v.Y() * v.Y();
    return RNIsLessOrEqual(distance_squared, radius_squared);
}



RNBoolean R2Contains(const R2Circle& circle, const R2Line& line)
{
    // Return whether circle contains line
    if (circle.IsFinite()) return FALSE;
    return TRUE;
}



RNBoolean R2Contains(const R2Circle& circle, const R2Ray& ray)
{
    // Return whether circle contains ray
    if (circle.IsFinite()) return FALSE;
    return TRUE;
}



RNBoolean R2Contains(const R2Circle& circle, const R2Span& span)
{
    // Return whether circle contains span
    if (!R2Contains(circle, span.Start())) return FALSE;
    if (!R2Contains(circle, span.End())) return FALSE;
    return TRUE;
}



RNBoolean R2Contains(const R2Circle& circle, const R2Halfspace& halfspace)
{
    // Return whether circle contains halfspace
    if (circle.IsFinite()) return FALSE;
    return TRUE;
}



RNBoolean R2Contains(const R2Circle& circle, const R2Box& box)
{
    // Return whether circle contains box 
    R2Vector v = box.Centroid() - circle.Center();
    R2Point corner = box.Corner(v.Quadrant());
    return R2Contains(circle, corner);
}



RNBoolean R2Contains(const R2Circle& circle1, const R2Circle& circle2) 
{
    // Return whether circle1 contains circle2
    RNLength d = R2Distance(circle1.Center(), circle2.Center());
    return RNIsLess(d + circle2.Radius(), circle1.Radius());
}



RNBoolean R2Contains(const R2Circle& circle, const R2Shape& shape)
{
    // Return whether circle contains shape
    RNAbort("Not implemented");
    return FALSE;
}



RNBoolean R2Contains(const R2Polygon& polygon, const R2Point& point)
{
    // Return whether polygon contains point
    int ncrossings = 0;
    R2Ray ray(point, R2Vector(0.12345, 0.9923508));
    for (int i = 0; i < polygon.NPoints(); i++) {
        const R2Point& p0 = polygon.Point(i);
        const R2Point& p1 = polygon.Point((i+1)%polygon.NPoints());
        if (R2Intersects(ray, R2Span(p0, p1))) ncrossings++;
    }

    // Check odd-crossing rule
    return (ncrossings % 2) ? TRUE : FALSE;
}



RNBoolean R2Contains(const R2Polygon& polygon, const R2Line& line)
{
    // Return whether polygon contains line
    if (polygon.IsFinite()) return FALSE;
    return TRUE;
}



RNBoolean R2Contains(const R2Polygon& polygon, const R2Ray& ray)
{
    // Return whether polygon contains ray
    if (polygon.IsFinite()) return FALSE;
    return TRUE;
}



RNBoolean R2Contains(const R2Polygon& polygon, const R2Span& span)
{
    // Return whether polygon contains span (this is not right, because polygon can be concave)
    if (!R2Contains(polygon, span.Start())) return FALSE;
    if (!R2Contains(polygon, span.End())) return FALSE;
    return TRUE;
}



RNBoolean R2Contains(const R2Polygon& polygon, const R2Halfspace& halfspace)
{
    // Return whether polygon contains halfspace
    if (polygon.IsFinite()) return FALSE;
    return TRUE;
}



RNBoolean R2Contains(const R2Polygon& polygon, const R2Box& box)
{
    // Return whether polygon contains box 
    RNAbort("Not implemented");
    return FALSE;
}



RNBoolean R2Contains(const R2Polygon& polygon, const R2Circle& circle)
{
    // Return whether polygon contains circle 
    RNAbort("Not implemented");
    return FALSE;
}



RNBoolean R2Contains(const R2Polygon& polygon1, const R2Polygon& polygon2) 
{
    // Return whether polygon1 contains polygon2
    RNAbort("Not implemented");
    return FALSE;
}



RNBoolean R2Contains(const R2Polygon& polygon, const R2Shape& shape)
{
    // Return whether polygon contains shape
    RNAbort("Not implemented");
    return FALSE;
}



RNBoolean R2Contains(const R2Shape& shape, const R2Point& point)
{
    // Return whether shape contains point
    RNAbort("Not implemented");
    return FALSE;
}



RNBoolean R2Contains(const R2Shape& shape, const R2Line& line)
{
    // Return whether shape contains line
    RNAbort("Not implemented");
    return FALSE;
}



RNBoolean R2Contains(const R2Shape& shape, const R2Ray& ray)
{
    // Return whether shape contains ray
    RNAbort("Not implemented");
    return FALSE;
}



RNBoolean R2Contains(const R2Shape& shape, const R2Span& span)
{
    // Return whether shape contains span
    RNAbort("Not implemented");
    return FALSE;
}



RNBoolean R2Contains(const R2Shape& shape, const R2Halfspace& halfspace)
{
    // Return whether shape contains halfspace
    RNAbort("Not implemented");
    return FALSE;
}



RNBoolean R2Contains(const R2Shape& shape, const R2Box& box)
{
    // Return whether shape contains box
    RNAbort("Not implemented");
    return FALSE;
}



RNBoolean R2Contains(const R2Shape& shape, const R2Circle& circle)
{
    // Return whether shape contains circle
    RNAbort("Not implemented");
    return FALSE;
}



RNBoolean R2Contains(const R2Shape& shape1, const R2Shape& shape2)
{
    // Return whether shape1 contains shape2
    RNAbort("Not implemented");
    return FALSE;
}





