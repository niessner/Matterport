/* Source file for GAPS distance utility */



/* Include files */

#include "R2Shapes/R2Shapes.h"



/* Public functions */

int R2InitDistance()
{
    // Return success 
    return TRUE;
}



void R2StopDistance()
{
}



RNLength R2Distance(const R2Point& point1, const R2Point& point2)
{
    // Return length of vector between points
    R2Vector v = point1 - point2;
    return v.Length();
}



RNLength R2Distance(const R2Point& point, const R2Line& line)
{
    // Return distance from point to line 
    RNLength d = point.X() * line.A() + point.Y() * line.B() + line.C();
    return (d < 0.0) ? -d : d;
}



RNLength R2Distance(const R2Point& point, const R2Ray& ray)
{
    // Check if start point is closest
    R2Vector v = point - ray.Start();
    RNScalar dir = v.Dot(ray.Vector());
    if (RNIsNegative(dir)) return v.Length();

    // Return distance from point to ray line
    return R2Distance(point, ray.Line());
}



RNLength R2Distance(const R2Point& point, const R2Span& span)
{
    // Check span
    if (RNIsZero(span.Length())) return R2Distance(point, span.Start());

    // Check if start point is closest
    R2Vector v1 = point - span.Start();
    RNScalar dir1 = v1.Dot(span.Vector());
    if (RNIsNegative(dir1)) return v1.Length();

    // Check if end point is closest
    R2Vector v2 = point - span.End();
    RNScalar dir2 = v2.Dot(span.Vector());
    if (RNIsPositive(dir2)) return v2.Length();

    // Return distance from point to span line
    return R2Distance(point, span.Line());
}



RNLength R2Distance(const R2Point& point, const R2Halfspace& halfspace)
{
    // Return distance from point to halfspace
    RNScalar d = R2SignedDistance(halfspace.Line(), point);
    if (RNIsNegative(d)) return -d;
    else return 0.0;
}



RNLength R2Distance(const R2Point& point, const R2Box& box)
{
    // Find axial distances from point to box
    RNLength dx, dy;
    if (RNIsGreater(point.X(), box.XMax())) dx = point.X() - box.XMax();
    else if (RNIsLess(point.X(), box.XMin())) dx = box.XMin()- point.X();
    else dx = 0.0;
    if (RNIsGreater(point.Y(), box.YMax())) dy = point.Y() - box.YMax();
    else if (RNIsLess(point.Y(), box.YMin())) dy = box.YMin()- point.Y();
    else dy = 0.0;
    
    // Return distance between point and closest point in box 
    if (dy == 0.0) return dx;
    else if (dx == 0.0) return dy;
    else return sqrt(dx*dx + dy*dy);
}



RNLength R2Distance(const R2Point& point, const R2Polyline& polyline)
{
    // Compute distance to each segment
    // This could be a lot faster
    RNLength closest_distance = FLT_MAX;
    for (int i = 1; i < polyline.NPoints(); i++) {
        const R2Point& p0 = polyline.Point(i-1);
        const R2Point& p1 = polyline.Point(i);
        R2Span span(p0, p1);
        RNLength d = R2Distance(span, point);
        if (d < closest_distance) closest_distance = d;
    }

    // Return distance
    return closest_distance;
}


RNLength R2Distance(const R2Point& point, const R2Polygon& polygon)
{
    // Check if polygon contains point
    if (R2Contains(polygon, point)) return 0;
  
    // Compute distance to each segment
    // This could be a lot faster
    RNLength closest_distance = FLT_MAX;
    for (int i = 0; i < polygon.NPoints(); i++) {
        const R2Point& p0 = polygon.Point(i);
        const R2Point& p1 = polygon.Point((i+1)%polygon.NPoints());
        R2Span span(p0, p1);
        RNLength d = R2Distance(span, point);
        if (d < closest_distance) closest_distance = d;
    }

    // Return distance
    return closest_distance;
}




RNLength R2Distance(const R2Point& point, const R2Arc& arc)
{
    // Figure sector - positive ds are inside arc sector
    RNLength d1 = R2SignedDistance(R2Line(arc.StartPoint(), arc.Center()), point);
    RNLength d2 = R2SignedDistance(R2Line(arc.Center(), arc.StopPoint()), point);
    if (d1 > 0.0) {
        if (d2 > 0.0) {
  	    // Point is in arc sector
	    RNLength d = R2Distance(arc.Center(), point) - arc.Radius();
	    return (d < 0) ? -d : d;
	}
        else {
	    // Point may be in arc sector 
	    if (arc.SweepAngle() > RN_PI) {
	        // Point is in arc sector
	        RNLength d = R2Distance(arc.Center(), point) - arc.Radius();
		return (d < 0) ? -d : d;
	    }
	    else {
	        // Point is closest to stop point
	        return R2Distance(arc.StopPoint(), point);
	    }
        }
    }
    else {
        if (d2 > 0.0) {
	    // Point may be in arc sector 
	    if (arc.SweepAngle() > RN_PI) {
	        // Point is in arc sector
	        RNLength d = R2Distance(arc.Center(), point) - arc.Radius();
		return (d < 0) ? -d : d;
	    }
	    else {
	        // Point is closest to start point
	        return R2Distance(arc.StartPoint(), point);
	    }
        }
        else {
  	    // Point is outside arc sector
	    RNLength d = R2Distance(arc.Center(), point) - arc.Radius();
	    return (d < 0) ? -d : d;
        }
    }
}



RNLength R2Distance(const R2Point& point, const R2Circle& circle)
{
    // Return distance from point to circle
    RNLength d = R2Distance(point, circle.Center()) - circle.Radius();
    return ((d > 0.0) ? d : 0.0);
}



RNLength R2SquaredDistance(const R2Point& point1, const R2Point& point2)
{
    // Return squared length of vector between points
    R2Vector v = point1 - point2;
    return v.Dot(v);
}



RNLength R2Distance(const R2Line& line1, const R2Line& line2)
{
    // Return distance from line to line 
    RNAbort("Not implemented");
    return 0.0;
}



RNLength R2Distance(const R2Line& line, const R2Ray& ray)
{
    RNAbort("Not implemented");
    return 0.0;
}



RNLength R2Distance(const R2Line& line, const R2Span& span)
{
    // Return distance from span to line
    RNScalar d = R2SignedDistance(line, span);
    if (RNIsPositive(d)) return d;
    else if (RNIsNegative(d)) return -d;
    else return 0.0;
}



RNLength R2Distance(const R2Line& line, const R2Halfspace& halfspace)
{
    // Return distance from line to halfspace
    RNAbort("Not implemented");
    return 0.0;
}



RNLength R2Distance(const R2Line& line, const R2Box& box)
{
    RNAbort("Not implemented");
    return 0.0;
}



RNLength R2Distance(const R2Line& line, const R2Circle& circle)
{
    // Return distance from line to circle
    RNLength d = R2Distance(line, circle.Center()) - circle.Radius();
    return ((d > 0.0) ? d : 0.0);
}



RNLength R2SignedDistance(const R2Line& line, const R2Point& point)
{
    // Return signed distance from point to line 
    return (point.X()*line.A() + point.Y()*line.B() + line.C());
}



RNLength R2SignedDistance(const R2Line& line1, const R2Line& line2)
{
    // Return signed distance from line to line
    RNScalar dot = line1.Vector().Dot(line2.Vector());
    if (RNIsEqual(dot, 1.0)) return (line1.C() - line2.C());
    else if (RNIsEqual(dot, -1.0)) return (line1.C() + line2.C());
    else return 0.0;
}



RNLength R2SignedDistance(const R2Line& line, const R2Ray& ray)
{
    RNAbort("Not implemented");
    return 0.0;
}



RNLength R2SignedDistance(const R2Line& line, const R2Span& span)
{
    // Return signed distance from line to span
    RNLength d1 = R2SignedDistance(line, span.Start());
    if (RNIsPositive(d1)) {
	// Start point is above line
	RNLength d2 = R2SignedDistance(line, span.End());
	if (RNIsPositive(d2)) return ((d1 > d2) ? d2 : d1);
	else return 0.0;
    }
    else if (RNIsNegative(d1)) {
	// Start point is below line
	RNLength d2 = R2SignedDistance(line, span.End());
	if (RNIsNegative(d2)) return ((d1 > d2) ? d1 : d2);
	else return 0.0;
    }
    else {
	// Start point is on line
	return 0.0;
    }
}



RNLength R2SignedDistance(const R2Line& line, const R2Halfspace& halfspace)
{
    // Return signed distance from line to halfspace
    RNAbort("Not implemented");
    return 0.0;
}



RNLength R2SignedDistance(const R2Line& line, const R2Box& box)
{
    // Return signed distance from line to box
    RNQuadrant quadrant = line.Normal().Quadrant();
    RNScalar d1 = R2SignedDistance(line, box.Corner(~quadrant & 0x3));
    if (RNIsPositiveOrZero(d1)) return d1;
    RNScalar d2 = R2SignedDistance(line, box.Corner(quadrant));
    if (RNIsNegative(d2)) return d2;
    else return 0.0;
}



RNLength R2SignedDistance(const R2Line& line, const R2Circle& circle)
{
    RNLength d = R2SignedDistance(line, circle.Center());
    if (d < 0.0) {
	d += circle.Radius();
	if (d > 0.0) return 0.0;
    }
    else if (d > 0.0) {
	d -= circle.Radius();
	if (d < 0.0) return 0.0;
    }
    return d;
}



RNLength R2Distance(const R2Ray& ray1, const R2Ray& ray2)
{
    RNAbort("Not implemented");
    return 0.0;
}



RNLength R2Distance(const R2Ray& ray, const R2Span& span)
{
    RNAbort("Not implemented");
    return 0.0;
}



RNLength R2Distance(const R2Ray& ray, const R2Halfspace& halfspace)
{
    // Return distance from ray to halfspace
    RNScalar d = R2SignedDistance(halfspace.Line(), ray);
    if (RNIsNegative(d)) return -d;
    else return 0.0;
}



RNLength R2Distance(const R2Ray& ray, const R2Box& box)
{
    RNAbort("Not implemented");
    return 0.0;
}



RNLength R2Distance(const R2Ray& ray, const R2Circle& circle)
{
    // Return distance from ray to circle
    RNLength d = R2Distance(ray, circle.Center()) - circle.Radius();
    return ((d > 0.0) ? d : 0.0);
}



RNLength R2Distance(const R2Span& span1, const R2Span& span2)
{
    // Look at this ???
    RNAbort("Not implemented");
    return 0.0;
}



RNLength R2Distance(const R2Span& span, const R2Halfspace& halfspace)
{
    // Return distance from span to halfspace
    RNScalar d = R2SignedDistance(halfspace.Line(), span);
    if (RNIsNegative(d)) return -d;
    else return 0.0;
}



RNLength R2Distance(const R2Span& span, const R2Box& box)
{
    RNAbort("Not implemented");
    return 0.0;
}



RNLength R2Distance(const R2Span& span, const R2Circle& circle)
{
    // Return distance from span to circle
    RNLength d = R2Distance(span, circle.Center()) - circle.Radius();
    return ((d > 0.0) ? d : 0.0);
}



RNLength R2Distance(const R2Halfspace& halfspace1, const R2Halfspace& halfspace2)
{
    RNAbort("Not implemented");
    return 0.0;
}



RNLength R2Distance(const R2Halfspace& halfspace, const R2Box& box)
{
    RNAbort("Not implemented");
    return 0.0;
}



RNLength R2Distance(const R2Halfspace& halfspace, const R2Circle& circle)
{
    // Return distance from halfspace to circle
    RNLength d = R2Distance(halfspace, circle.Center()) - circle.Radius();
    return ((d > 0.0) ? d : 0.0);
}



RNLength R2Distance(const R2Box& box1, const R2Box& box2)
{
    // Find axial distances from box1 to box2
    RNLength dx, dy;
    if (RNIsGreater(box1.XMin(), box2.XMax())) dx = box1.XMin() - box2.XMax();
    else if (RNIsGreater(box2.XMin(), box1.XMax())) dx = box2.XMin() - box1.XMax();
    else dx = 0.0;
    if (RNIsGreater(box1.YMin(), box2.YMax())) dy = box1.YMin() - box2.YMax();
    else if (RNIsGreater(box2.YMin(), box1.YMax())) dy = box2.YMin() - box1.YMax();
    else dy = 0.0;
    
    // Return distance between point and closest point in box 
    if (dy == 0.0) return dx;
    else if (dx == 0.0) return dy;
    else return sqrt(dx*dx + dy*dy);
}



RNLength R2Distance(const R2Box& box, const R2Circle& circle)
{
    // Return distance from box to circle
    RNLength d = R2Distance(box, circle.Center()) - circle.Radius();
    return ((d > 0.0) ? d : 0.0);
}



RNLength R2Distance(const R2Circle& circle1, const R2Circle& circle2)
{
    // Return distance from circle to circle
    RNLength d = R2Distance(circle2, circle1.Center()) - circle1.Radius();
    return ((d > 0.0) ? d : 0.0);
}








