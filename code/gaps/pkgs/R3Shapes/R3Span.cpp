/* Source file for the GAPS span class */



/* Include files */

#include "R3Shapes/R3Shapes.h"



/* Public variables */

R3Span R3null_span(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);



/* Public functions */

int R3InitSpan()
{
    /* Return success */
    return TRUE;
}



void R3StopSpan()
{
}



R3Span::
R3Span(void)
{
}



R3Span::
R3Span(const R3Span& span)
    : ray(span.ray),
      end(span.end),
      length(span.length)
{
}



R3Span::
R3Span(const R3Point& point, const R3Vector& vector)
    : end(point + vector)
{
    R3Vector v(vector);
    length = v.Length();
    if (RNIsNotZero(length)) v /= length;
    ray.Reset(point, v, TRUE);
}



R3Span::
R3Span(const R3Point& point1, const R3Point& point2)
    : end(point2)
{
    R3Vector v(point2 - point1);
    length = v.Length();
    if (RNIsNotZero(length)) v /= length;
    ray.Reset(point1, v, TRUE);
}



R3Span::
R3Span(RNCoord x1, RNCoord y1, RNCoord z1, RNCoord x2, RNCoord y2, RNCoord z2)
    : end(x2, y2, z2)
{
    R3Point p1(x1, y1, z1);
    R3Point p2(x2, y2, z2);
    R3Vector v(p2 - p1);
    length = v.Length();
    if (RNIsNotZero(length)) v /= length;
    ray.Reset(p1, v, TRUE);
}



const R3Point R3Span::
Point(RNScalar t) const
{
    // Return point along span
    return (Start() + Vector() * t);
}



const R3Point R3Span::
Midpoint(void) const
{
    // Return midpoint of span
    return (Start() + End()) * 0.5;
}



const R3Box R3Span::
BBox(void) const
{
    // Return bounding box 
    R3Box bbox(Start(), Start());
    bbox.Union(End());
    return bbox;
}



const R3Sphere R3Span::
BSphere(void) const
{
    // Return bounding sphere 
    return R3Sphere(Midpoint(), 0.5 * Length());
}



const RNScalar R3Span::
T(const R3Point& point) const
{
    // Return parametric value of closest point on span
    if (length == 0.0) return 0.0;
    R3Vector topoint = point - Start();
    return Vector().Dot(topoint);
}



const RNBoolean R3Span::
IsPoint(void) const
{
    // Return whether span covers a single point
    return (length == 0.0);
}



void R3Span::
Reposition(int k, const R3Point& point)
{
    // Set one endpoint of span
    if (k == 0) ray = R3Ray(point, end);
    else { end = point; ray.Align(end - ray.Start()); }
    length = R3Distance(Start(), End());
}



void R3Span::
Align(const R3Vector& vector)
{
    // Set vector of span
    ray.Align(vector);
    end = Start() + Vector() * length;
}



void R3Span::
Transform (const R3Transformation& transformation)
{
    // Transform span
    end.Transform(transformation);
    ray.Transform(transformation);
    length = R3Distance(Start(), End());
}



void R3Span::
InverseTransform (const R3Transformation& transformation)
{
    // Transform span
    end.InverseTransform(transformation);
    ray.InverseTransform(transformation);
    length = R3Distance(Start(), End());
}



void R3Span::
Flip(void)
{
    // Flip direction of span
    R3Point swap = Start();
    ray.Reposition(end);
    end = swap;
    ray.Flip();
}



void R3Span::
Translate(const R3Vector& vector)
{
    // Move endpoints of span
    ray.Translate(vector);
    end += vector;
}



void R3Span::
Reset(const R3Point& point1, const R3Point& point2)
{
    // Reset span
    ray = R3Ray(point1, point2);
    end = point2;
    length = R3Distance(point1, point2);
}


RNClassID R3Span::
Clip(const R3Plane& plane)
{
    // Characterize endpoints with respect to plane
    RNScalar d1 = R3SignedDistance(plane, Start());
    RNScalar d2 = R3SignedDistance(plane, End());

    // Clip span to plane
    if (RNIsNegative(d1)) {
	if (RNIsNegative(d2)) {
	    // Both points are below plane ???
	    Reset(Start(), Start());
	    return RN_NULL_CLASS_ID;
	}
	else if (RNIsPositive(d2)) {
	    // Start is below, end is above -- move start to plane
	    Reset((Start() * d2 + End() * -d1) / (d2 - d1), End());
	    return R3_SPAN_CLASS_ID;
	}
	else {
	    // Start is below, end is on -- move start to end
	    Reset(End(), End());
	    return R3_POINT_CLASS_ID;
	}
    }
    else if (RNIsPositive(d1)) {
	if (RNIsNegative(d2)) {
	    // Start is above, end is below -- move end to plane
	    Reset(Start(), (Start() * -d2 + End() * d1) / (d1 - d2));
	    return R3_SPAN_CLASS_ID;
	}
	else {
	    // Start is above, end is on or above
	    return R3_SPAN_CLASS_ID;
	}
    }
    else {
	if (RNIsNegative(d2)) {
	    // Start is on, end is below -- move end to start
	    Reset(Start(), Start());
	    return R3_POINT_CLASS_ID;
	}
	else if (RNIsPositive(d2)){
	    // Start is on, End is above
	    return R3_SPAN_CLASS_ID;
	}
	else {
	    // Start is on, end is on
	    return R3_POINT_CLASS_ID;
	}
    }
}



const RNBoolean R3Span::
operator==(const R3Span& span) const
{
    // Return whether span is equal
    if (Start() != span.Start()) return FALSE;
    if (End() != span.End()) return FALSE;
    return TRUE;
}




