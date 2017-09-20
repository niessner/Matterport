/* Source file for the R2 arc class */



/* Include files */

#include "R2Shapes/R2Shapes.h"



/* Public variables */

const R2Arc R2null_arc(R2Point(0.0, 0.0), -1.0, 0.0, 0.0);
const R2Arc R2zero_arc(R2Point(0.0, 0.0), 0.0, 0.0, 0.0);



/* Public functions */

int 
R2InitArc()
{
    /* Return success */
    return TRUE;
}



void 
R2StopArc()
{
}



R2Arc::
R2Arc(void)
{
}



R2Arc::
R2Arc(const R2Arc& arc)
    : circle(arc.Circle()),
      start(arc.start),
      stop(arc.stop)
{
}



R2Arc::
R2Arc(const R2Point& center, RNLength radius, RNAngle theta1, RNAngle theta2)
    : circle(center, radius),
      start(theta1),
      stop(theta2)
{
    // Needed because don't know if RN_TWO_PI is initialized
    RNScalar TWO_PI = 2.0 * 3.14159265358979323846;
    while (start >= TWO_PI) start -= TWO_PI;
    while (start < 0.0) start += TWO_PI;
    while (stop >= (start + TWO_PI)) stop -= TWO_PI;
    while (stop < start) stop += TWO_PI;
    assert((0.0 <= start) && (start <= TWO_PI));
    assert((start <= stop) && (stop <= (start + TWO_PI)));
}



const RNBoolean R2Arc::
IsPoint(void) const
{
    return (circle.IsPoint() || (start == stop));
}



const RNBoolean R2Arc::
IsLinear(void) const
{
    return IsPoint();
}



const RNBoolean R2Arc::
IsConvex(void) const
{
    return IsPoint();
}



#if FALSE

const RNArea R2Arc::
Area(void) const
{
    // Return area of arc
    return circle.Area() * (stop - start) / RN_TWO_PI;
}

#endif



const R2Point R2Arc::
Point(RNAngle angle) const
{
    // Return start point of arc ???
    return Center() + Radius() * R2Vector(cos(angle), sin(angle));
}



const R2Point R2Arc::
Centroid(void) const
{
    // Return centroid of arc ???
    return Center();
}



const R2Shape& R2Arc::
BShape(void) const
{
    // Return circle
    return circle;
}



const R2Box R2Arc::
BBox(void) const
{
    // Return bounding box of arc ???
    R2Box bbox(R2null_box);
    bbox.Union(StartPoint());
    bbox.Union(StopPoint());
    if ((start < 0.0) && (stop > 0.0)) bbox.Union(Center() + R2posx_point * Radius());
    if ((start < RN_PI_OVER_TWO) && (stop > RN_PI_OVER_TWO)) bbox.Union(Center() + R2posy_point * Radius());
    if ((start < RN_PI) && (stop > RN_PI)) bbox.Union(Center() + R2negx_point * Radius());
    if ((start < RN_THREE_PI_OVER_TWO) && (stop > RN_THREE_PI_OVER_TWO)) bbox.Union(Center() + R2negy_point * Radius());
    return bbox;
}



const R2Circle R2Arc::
BCircle(void) const
{
    // Return circle
    return circle;
}



void R2Arc::
Empty(void)
{
    // Empty arc
    *this = R2null_arc;
}



void R2Arc::
Translate(const R2Vector& vector)
{
    // Move arc center
    circle.Translate(vector);
}



void R2Arc::
Reposition(const R2Point& center)
{
    // Set arc center
    circle.Reposition(center);
}



void R2Arc::
Resize(RNLength radius) 
{
    // Set arc radius
    circle.Resize(radius);
}



void R2Arc::
SetStartAngle(RNAngle theta)
{
    // Set arc start angle
    start = theta;
    while (start >= RN_TWO_PI) start -= RN_TWO_PI;
    while (start < 0.0) start += RN_TWO_PI;
    while (stop >= (start + RN_TWO_PI)) { stop -= RN_TWO_PI; }
    while (stop < start) { stop += RN_TWO_PI; }
    assert((0.0 <= start) && (start <= RN_TWO_PI));
    assert((start <= stop) && (stop <= (start + RN_TWO_PI)));
}



void R2Arc::
SetStopAngle(RNAngle theta)
{
    // Set arc stop angle
    stop = theta;
    while (stop >= (start + RN_TWO_PI)) { stop -= RN_TWO_PI; }
    while (stop < start) { stop += RN_TWO_PI; }
    assert((0.0 <= start) && (start <= RN_TWO_PI));
    assert((start <= stop) && (stop <= (start + RN_TWO_PI)));
}



void R2Arc::
Transform (const R2Transformation& transformation)
{
    // Transform circle
    circle.Transform(transformation);

    // Transform angles
    RNAbort("Not implemented");
}













