/* Source file for the GAPS span class */



/* Include files */

#include "R2Shapes/R2Shapes.h"



/* Public functions */

int R2InitSpan()
{
    /* Return success */
    return TRUE;
}



void R2StopSpan()
{
}



R2Span::
R2Span(void)
{
}



R2Span::
R2Span(const R2Span& span)
    : ray(span.ray),
      end(span.end),
      length(span.length)
{
}



R2Span::
R2Span(const R2Point& point, const R2Vector& vector)
    : end(point + vector)
{
    R2Vector v(vector);
    length = v.Length();
    if (RNIsNotZero(length)) v /= length;
    ray.Reset(point, v, TRUE);
}



R2Span::
R2Span(const R2Point& point1, const R2Point& point2)
    : end(point2)
{
    R2Vector v(point2 - point1);
    length = v.Length();
    if (RNIsNotZero(length)) v /= length;
    ray.Reset(point1, v, TRUE);
}



R2Span::
R2Span(RNCoord x1, RNCoord y1, RNCoord x2, RNCoord y2)
    : end(x2, y2)
{
    R2Point p1(x1, y1);
    R2Point p2(x2, y2);
    R2Vector v(p2 - p1);
    length = v.Length();
    if (RNIsNotZero(length)) v /= length;
    ray.Reset(R2Point(x1, y1), v, TRUE);
}



const R2Point R2Span::
Point(RNScalar t) const
{
    // Return point along span
    return (Start() + Vector() * t);
}



const R2Point R2Span::
Midpoint(void) const
{
    // Return midpoint of span
    return (Start() + End()) * 0.5;
}



const R2Box R2Span::
BBox(void) const
{
    // Return bounding box 
    R2Box bbox(Start(), Start());
    bbox.Union(End());
    return bbox;
}



const R2Circle R2Span::
BCircle(void) const
{
    // Return bounding circle
    return R2Circle(Midpoint(), 0.5 * Length());
}



const RNBoolean R2Span::
IsPoint(void) const
{
    // Return whether span covers a single point
    return (length == 0.0);
}



void R2Span::
Reposition(int k, const R2Point& point)
{
    // Set one endpoint of span
    if (k == 0) ray = R2Ray(point, end);
    else { end = point; ray.Align(end - ray.Start()); }
    length = R2Distance(Start(), End());
}



void R2Span::
Project(const R2Line& line)
{
    // Project span onto line
    length *= Vector().Dot(line.Vector());
    if (length < 0.0) length = -length;
    ray.Project(line);
    end.Project(line);
}



void R2Span::
Mirror(const R2Line& line)
{
    // Mirror span over line
    ray.Mirror(line);
    end.Mirror(line);
}



void R2Span::
Align(const R2Vector& vector)
{
    // Set vector of span
    ray.Align(vector);
    end = Start() + Vector() * length;
}



void R2Span::
Transform (const R2Transformation& transformation)
{
    // Transform span
    end.Transform(transformation);
    ray.Transform(transformation);
    length = R2Distance(Start(), End());
}



void R2Span::
InverseTransform (const R2Transformation& transformation)
{
    // Transform span
    end.InverseTransform(transformation);
    ray.InverseTransform(transformation);
    length = R2Distance(Start(), End());
}



void R2Span::
Flip(void)
{
    // Flip direction of span
    R2Point swap = Start();
    ray.Reposition(end);
    end = swap;
    ray.Flip();
}



void R2Span::
Translate(const R2Vector& vector)
{
    // Move endpoints of span
    ray.Translate(vector);
    end += vector;
}



void R2Span::
Reset(const R2Point& point1, const R2Point& point2)
{
    // Reset span
    ray = R2Ray(point1, point2);
    end = point2;
    length = R2Distance(point1, point2);
}



RNClassID R2Span::
Clip(const R2Line& line)
{
    // Characterize endpoints with respect to line
    RNScalar d1 = R2SignedDistance(line, Start());
    RNScalar d2 = R2SignedDistance(line, End());

    // Clip span to line
    if (RNIsNegative(d1)) {
	if (RNIsNegative(d2)) {
	    // Both points are below line ???
	    Reset(Start(), Start());
	    return RN_NULL_CLASS_ID;
	}
	else if (RNIsPositive(d2)) {
	    // Start is below, end is above -- move start to line
	    Reset((Start() * d2 + End() * -d1) / (d2 - d1), End());
	    return R2_SPAN_CLASS_ID;
	}
	else {
	    // Start is below, end is on -- move start to end
	    Reset(End(), End());
	    return R2_POINT_CLASS_ID;
	}
    }
    else if (RNIsPositive(d1)) {
	if (RNIsNegative(d2)) {
	    // Start is above, end is below -- move end to line
	    Reset(Start(), (Start() * -d2 + End() * d1) / (d1 - d2));
	    return R2_SPAN_CLASS_ID;
	}
	else {
	    // Start is above, end is on or above
	    return R2_SPAN_CLASS_ID;
	}
    }
    else {
	if (RNIsNegative(d2)) {
	    // Start is on, end is below -- move end to start
	    Reset(Start(), Start());
	    return R2_POINT_CLASS_ID;
	}
	else {
	    // Start is on, end is on or above
	    return R2_SPAN_CLASS_ID;
	}
    }
}



const RNBoolean R2Span::
operator==(const R2Span& span) const
{
    // Return whether span is equal
    if (Start() != span.Start()) return FALSE;
    if (End() != span.End()) return FALSE;
    return TRUE;
}




