/* Source file for R2 miscellaneous relationship utility */



/* Include files */

#include "R2Shapes/R2Shapes.h"



RNAngle 
R2InteriorAngle(const R2Vector& vector1, const R2Vector& vector2)
{
    // Return angle between vectors
    RNScalar d1 = vector1.Length();
    if (RNIsZero(d1)) return 0.0;
    RNScalar d2 = vector2.Length();
    if (RNIsZero(d2)) return 0.0;
    RNScalar cosine = vector1.Dot(vector2) / (d1 * d2);
    assert(RNIsGreaterOrEqual(cosine, -1.0) && RNIsLessOrEqual(cosine, 1.0));
    if (cosine >= 1.0) return 0.0;
    else if (cosine <= -1.0) return RN_PI;
    else return acos(cosine); 
}



RNAngle 
R2InteriorAngle(const R2Vector& vector, const R2Line& line)
{
    // Return angle between vector and line
    return R2InteriorAngle(vector, line.Vector());
}



RNAngle 
R2InteriorAngle(const R2Vector& vector, const R2Ray& ray)
{
    // Return angle between vector and ray
    return R2InteriorAngle(vector, ray.Vector());
}



RNAngle 
R2InteriorAngle(const R2Vector& vector, const R2Span& span)
{
    // Return angle between vector and span
    return R2InteriorAngle(vector, span.Vector());
}



RNArea 
R2IntersectionArea(const R2Box& box1, const R2Box& box2)
{
  R2Box intersection(box1);
  intersection.Intersect(box2);
  if (intersection.IsEmpty()) return 0;
  return intersection.Area();
}



RNArea 
R2IntersectionArea(const R2Circle& circle1, const R2Circle& circle2)
{
  RNArea area = 0;
  RNLength r1 = circle1.Radius();
  if (RNIsZero(r1)) return 0;
  RNLength r2 = circle2.Radius();
  if (RNIsZero(r2)) return 0;
  RNLength d = R2Distance(circle1.Center(), circle2.Center());
  RNLength rmin = (r1 < r2) ? r1 : r2;
  RNLength rmax = (r1 > r2) ? r1 : r2;
  if (RNIsLessOrEqual(d + rmin, rmax)) {
    // Compute area of smaller circle
    area = RN_PI * rmin * rmin;
  }
  else if (RNIsLessOrEqual(d, r1 + r2)) {
    // Compute area of "asymmetric lens"
    // http://mathworld.wolfram.com/Circle-CircleIntersection.html
    RNScalar a1 = r1*r1*acos((d*d + r1*r1 - r2*r2) / (2*d*r1));
    RNScalar a2 = r2*r2*acos((d*d + r2*r2 - r1*r1) / (2*d*r2));
    RNScalar a3 = (-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2);
    if (a3 > 0) area = a1 + a2 - 0.5 * sqrt(a3);
  }

  // Return area
  return area;
}




RNBoolean 
R2Abuts(const R2Box& box1, const R2Box& box2) 
{
    // Return whether box1 abuts box2
    if (!R2Intersects(box1, box2)) return FALSE;
    if (RNIsEqual(box1.XMin(), box2.XMax())) return TRUE;
    if (RNIsEqual(box1.XMax(), box2.XMin())) return TRUE;
    if (RNIsEqual(box1.YMin(), box2.YMax())) return TRUE;
    if (RNIsEqual(box1.YMax(), box2.YMin())) return TRUE;
    return FALSE;
}



int 
R2Splits(const R2Line& line, const R2Span& span)
{
    // This is just a simple version 

    // Classify span
    // Result values: 0 = on, 1 = below, 2 = above, 3 = crossing
    RNScalar d1 = R2SignedDistance(line, span.Start());
    RNScalar d2 = R2SignedDistance(line, span.End());
    if (RNIsNegative(d1)) {
        if (RNIsPositive(d2)) return 3;
	else return 1;
    }
    else if (RNIsPositive(d1)) {
        if (RNIsNegative(d2)) return 3;
	else return 2;
    }
    else {
        if (RNIsNegative(d2)) return 1;
        else if (RNIsPositive(d2)) return 2;
	else return 0;
    }
}



int 
R2Splits(const R2Line& line, const R2Span& span, R2Span *below_result, R2Span *above_result)
{
    // Classify span
    // Result values: 0 = on, 1 = below, 2 = above, 3 = crossing
    RNScalar d1 = R2SignedDistance(line, span.Start());
    RNScalar d2 = R2SignedDistance(line, span.End());
    if (RNIsNegative(d1)) {
        if (RNIsPositive(d2)) {
  	    // Span crosses line
	    if (below_result || above_result) {
	        // Compute intersection
		R2Point point = span.Point(-d1);
		if (below_result) *below_result = R2Span(span.Start(), point);
		if (above_result) *above_result = R2Span(point, span.End());
	    }
	    return 3;
	}
	else {
  	    // Span is below line
	    if (below_result) *below_result = span; 
	    return 1;
	}
    }
    else if (RNIsPositive(d1)) {
        if (RNIsNegative(d2)) {
  	    // Span crosses line
	    if (below_result || above_result) {
	        // Compute intersection
		R2Point point = span.Point(d1);
		if (below_result) *below_result = R2Span(span.End(), point);
		if (above_result) *above_result = R2Span(point, span.Start());
	    }
	    return 3;
	}
	else {
  	    // Span is above line
	    if (above_result) *above_result = span; 
	    return 2;
	}
    }
    else {
        if (RNIsNegative(d2)) {
  	    // Span is below line
	    if (below_result) *below_result = span; 
	    return 1;
	}
        else if (RNIsPositive(d2)) {
  	    // Span is above line
	    if (above_result) *above_result = span; 
	    return 2;
	}
	else {
  	    // Span is on line
	    if (below_result) *below_result = span; 
	    if (above_result) *above_result = span; 
	    return 0;
	}
    }
}





