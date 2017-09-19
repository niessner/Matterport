/* Source file for R3 miscellaneous relationship utility */



/* Include files */

#include "R3Shapes/R3Shapes.h"



R3Point 
R3Interpolate(const R3Point& point1, const R3Point& point2, RNScalar t1, RNScalar t2, RNScalar t)
{
  // Check boundary cases
  if (RNIsEqual(t, t1)) return point1;
  if (RNIsEqual(t, t2)) return point2;
  RNScalar denom = t1 - t2;
  if (RNIsZero(denom)) return point1;

  // Interpolate points
  RNScalar mu = (t - t1) / (t2 - t1);
  R3Vector vector = point2 - point1;
  R3Point point = point1 + mu * vector;
  return point;
}



RNAngle 
R3InteriorAngle(const R3Vector& vector1, const R3Vector& vector2)
{
    // Return angle between vectors
    RNScalar d1 = vector1.Length();
    if (RNIsZero(d1)) return 0.0;
    RNScalar d2 = vector2.Length();
    if (RNIsZero(d2)) return 0.0;
    RNScalar cosine = vector1.Dot(vector2) / (d1 * d2);
    if (cosine >= 1.0) return 0.0;
    else if (cosine <= -1.0) return RN_PI;
    else return acos(cosine); 
}



RNAngle 
R3InteriorAngle(const R3Vector& vector, const R3Line& line)
{
    // Return angle between vector and line
    return R3InteriorAngle(vector, line.Vector());
}



RNAngle 
R3InteriorAngle(const R3Vector& vector, const R3Ray& ray)
{
    // Return angle between vector and ray
    return R3InteriorAngle(vector, ray.Vector());
}



RNAngle 
R3InteriorAngle(const R3Vector& vector, const R3Span& span)
{
    // Return angle between vector and span
    return R3InteriorAngle(vector, span.Vector());
}



RNAngle 
R3InteriorAngle(const R3Vector& vector, const R3Plane& plane)
{
    // Return angle between vector and plane
    return (RN_PI_OVER_TWO - R3InteriorAngle(vector, plane.Normal()));
}



RNBoolean 
R3Abuts(const R3Box& box1, const R3Box& box2) 
{
    // Return whether box1 abuts box2
    if (!R3Intersects(box1, box2)) return FALSE;
    if (RNIsEqual(box1.XMin(), box2.XMax())) return TRUE;
    if (RNIsEqual(box1.XMax(), box2.XMin())) return TRUE;
    if (RNIsEqual(box1.YMin(), box2.YMax())) return TRUE;
    if (RNIsEqual(box1.YMax(), box2.YMin())) return TRUE;
    if (RNIsEqual(box1.ZMin(), box2.ZMax())) return TRUE;
    if (RNIsEqual(box1.ZMax(), box2.ZMin())) return TRUE;
    return FALSE;
}



int 
R3Splits(const R3Plane& plane, const R3Span& span)
{
    // Result values:
    // 0 = on, 1 = below, 2 = above, 3 = crossing
    int result = 0;

    // Classify start point with respect to plane
    RNScalar d1 = R3SignedDistance(plane, span.Start());
    if (RNIsNegative(d1)) result |= 0x1;
    if (RNIsPositive(d1)) result |= 0x2;

    // Classify end point with respect to plane
    RNScalar d2 = R3SignedDistance(plane, span.End());
    if (RNIsNegative(d2)) result |= 0x1;
    if (RNIsPositive(d2)) result |= 0x2;

    // Return result
    return result;
}

    

int 
R3Splits(const R3Plane& plane, const R3Span& span, R3Span *below_result, R3Span *above_result)
{
    // Characterize endpoints with respect to plane
    RNScalar d1 = R3SignedDistance(plane, span.Start());
    RNScalar d2 = R3SignedDistance(plane, span.End());

    // Clip span to plane
    if (RNIsNegative(d1)) {
	if (RNIsPositive(d2)) {
	    // Start is below, end is above
	    R3Point intersection_point = (span.Start() * d2 + span.End() * -d1) / (d2 - d1);
	    if (below_result) below_result->Reset(span.Start(), intersection_point);
	    if (above_result) above_result->Reset(intersection_point, span.End());
	    return 0x3;
	}
	else {
	    // Both points are on or below plane 
	    if (below_result) *below_result = span;
	    return 0x1;
	}
    }
    else if (RNIsPositive(d1)) {
	if (RNIsNegative(d2)) {
	    // Start is above, end is below
	    R3Point intersection_point = (span.Start() * -d2 + span.End() * d1) / (d1 - d2);
	    if (below_result) below_result->Reset(intersection_point, span.End());
	    if (above_result) above_result->Reset(span.Start(), intersection_point);
	    return 0x3;
	}
	else {
	    // Start is above, end is on or above
	    if (above_result) *above_result = span;
	    return 0x2;
	}
    }
    else {
	if (RNIsNegative(d2)) {
	    // Start is on, end is below -- move end to start
	    if (below_result) *below_result = span;
	    return 0x1;
	}
	else {
	    // Start is on, end is on or above
	    if (above_result) *above_result = span;
	    return 0x2;
	}
    }
}

    

int 
R3Splits(const R3Plane& plane, const R3Triangle& triangle)
{
    // Result values:
    // 0 = on, 1 = below, 2 = above, 3 = crossing
    int result = 0;

    // Classify triangle with respect to plane
    for (int i = 0; i < 3; i++) {
	RNScalar d = R3SignedDistance(plane, triangle.Vertex(i)->Position());
	if (!(result & 0x1) && RNIsNegative(d)) result |= 0x1;
	if (!(result & 0x2) && RNIsPositive(d)) result |= 0x2;
	if (result == 0x3) break;
    }

    // Return result
    return result;
}

    

