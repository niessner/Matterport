/* Source file for the R2 Box class */



/* Include files */

#include "R2Shapes/R2Shapes.h"



/* Public variables */

const R2Box R2zero_box(0.0, 0.0, 0.0, 0.0);
const R2Box R2unit_box(-1.0, -1.0, 1.0, 1.0);
const R2Box R2infinite_box(-RN_INFINITY, -RN_INFINITY, RN_INFINITY, RN_INFINITY);
const R2Box R2null_box( RN_INFINITY, RN_INFINITY, -RN_INFINITY, -RN_INFINITY);



/* Public functions */

int 
R2InitBox()
{
    /* Return success */
    return TRUE;
}



void 
R2StopBox()
{
}



R2Box::
R2Box(void)
{
}



R2Box::
R2Box(const R2Box& box)
    : minpt(box.minpt),
      maxpt(box.maxpt)
{
}



R2Box::
R2Box(const R2Point& minpt, const R2Point& maxpt)
    : minpt(minpt), 
      maxpt(maxpt)
{
}



R2Box::
R2Box(const R2Point& center, RNLength xradius, RNLength yradius)
    : minpt(center.X() - xradius, center.Y() - yradius), 
      maxpt(center.X() + xradius, center.Y() + yradius)
{
}



R2Box::
R2Box(RNCoord xmin, RNCoord ymin,
      RNCoord xmax, RNCoord ymax)
    : minpt(xmin, ymin),
      maxpt(xmax, ymax)
{
}



const RNBoolean R2Box::
IsPoint (void) const
{
    // Return whether bounding box contains just one point
    return ((minpt.X() == maxpt.X()) && (minpt.Y() == maxpt.Y()));
}



const RNBoolean R2Box::
IsLinear (void) const
{
    // Return whether bounding box contains just one line segment
    return (this->NDimensions() <= 1);
}



const RNBoolean R2Box::
IsConvex (void) const
{
    // All boxes are convex
    return TRUE;
}



const RNBoolean R2Box::
IsFinite (void) const
{
    // Return whether bounding box contains a finite amount of space
    if (IsEmpty()) return TRUE;
    if (!RNIsFinite(minpt.X())) return FALSE;
    if (!RNIsFinite(minpt.Y())) return FALSE;
    if (!RNIsFinite(maxpt.X())) return FALSE;
    if (!RNIsFinite(maxpt.Y())) return FALSE;
    return TRUE;
}



const R2Point R2Box::
Corner (RNQuadrant quadrant) const
{
    // Return corner point 
    switch (quadrant) {
    case RN_NN_QUADRANT: return minpt;
    case RN_NP_QUADRANT: return R2Point(minpt.X(), maxpt.Y());
    case RN_PN_QUADRANT: return R2Point(maxpt.X(), minpt.Y());
    case RN_PP_QUADRANT: return maxpt;
    default: RNAbort("Invalid quadrant for box corner"); return R2null_point;
    }
}



const R2Point R2Box::
Corner (RNDirection xdir, RNDirection ydir) const
{
    // Return corner point 
    return R2Point(Coord(xdir, RN_X), Coord(ydir, RN_Y));
}



const R2Box R2Box::
Quadrant (RNDirection xdir, RNDirection ydir) const
{
    // Return box in quadrant (xdir, ydir)
    R2Box result;
    R2Point centroid = Centroid();
    if (xdir == RN_HI) { result[RN_LO][RN_X] = centroid[RN_X]; result[RN_HI][RN_X] = maxpt[RN_X]; }
    else { result[RN_LO][RN_X] = minpt[RN_X]; result[RN_HI][RN_X] = centroid[RN_X]; }
    if (ydir == RN_HI) { result[RN_LO][RN_Y] = centroid[RN_Y]; result[RN_HI][RN_Y] = maxpt[RN_Y]; }
    else { result[RN_LO][RN_Y] = minpt[RN_Y]; result[RN_HI][RN_Y] = centroid[RN_Y]; }
    return result;
}



const int R2Box::
NDimensions (void) const
{
    // Return number of non-zero dimensions in box
    int ndimensions = 0;
    if (minpt.X() > maxpt.X()) return -1;
    if (minpt.X() < maxpt.X()) ndimensions++;
    if (minpt.Y() > maxpt.Y()) return -1;
    if (minpt.Y() < maxpt.Y()) ndimensions++;
    return ndimensions;
}



const RNAxis R2Box::
ShortestAxis (void) const
{
    // Return shortest axis
    RNLength dx = this->XLength();
    RNLength dy = this->YLength();
    return (dx < dy) ? RN_XAXIS : RN_YAXIS;
}



const RNAxis R2Box::
LongestAxis (void) const
{
    // Return longest axis
    RNLength dx = this->XLength();
    RNLength dy = this->YLength();
    return (dx > dy) ? RN_XAXIS : RN_YAXIS;
}



const RNLength R2Box::
ShortestAxisLength (void) const
{
    // Return shortest axis length
    RNLength dx = this->XLength();
    RNLength dy = this->YLength();
    return (dx < dy) ? dx : dy;
}



const RNLength R2Box::
LongestAxisLength (void) const
{
    // Return longest axis length
    RNLength dx = this->XLength();
    RNLength dy = this->YLength();
    return (dx > dy) ? dx : dy;
}



const RNLength R2Box::
DiagonalLength (void) const
{
    // Return length of box along diagonal
    return R2Distance(minpt, maxpt);
}



const R2Point R2Box::
ClosestPoint(const R2Point& point) const
{
    // Return closest point in box
    R2Point closest(point);
    if (closest.X() < XMin()) closest[RN_X] = XMin();
    else if (closest.X() > XMax()) closest[RN_X] = XMax();
    if (closest.Y() < YMin()) closest[RN_Y] = YMin();
    else if (closest.Y() > YMax()) closest[RN_Y] = YMax();
    return closest;
}



const RNArea R2Box::
Area (void) const
{
    // Return area of bounding box
    return (XLength() * YLength());
}



const R2Point R2Box::
Centroid (void) const
{
    // Return center point
    return R2Point(XCenter(), YCenter());
}



const R2Shape& R2Box::
BShape (void) const
{
    // Return self
    return *this;
}



const R2Box R2Box::
BBox (void) const
{
    // Return self - shape virtual function
    return *this;
}



const R2Circle R2Box::
BCircle (void) const
{
    // Return bounding circle 
    return R2Circle(Centroid(), DiagonalRadius());
}



void R2Box::
Empty (void) 
{
    // Copy empty box 
    *this = R2null_box;
}



void R2Box::
Inflate (RNScalar fraction) 
{
  // Scale box around centroid by fraction
  if (IsEmpty()) return;
  R2Point centroid = Centroid();
  R2Vector diagonal = Max() - centroid;
  Reset(centroid - fraction * diagonal, centroid + fraction * diagonal);
}



void R2Box::
Translate(const R2Vector& vector)
{
    // Move box by vector
    minpt += vector;
    maxpt += vector;
}



void R2Box::
Union (const R2Point& point) 
{
    // Expand this to include point
    if (minpt.X() > point.X()) minpt[0] = point.X();
    if (minpt.Y() > point.Y()) minpt[1] = point.Y();
    if (maxpt.X() < point.X()) maxpt[0] = point.X();
    if (maxpt.Y() < point.Y()) maxpt[1] = point.Y();
}



void R2Box::
Union (const R2Box& box) 
{
    // Expand this to include box
    if (minpt.X() > box.XMin()) minpt[0] = box.XMin();
    if (minpt.Y() > box.YMin()) minpt[1] = box.YMin();
    if (maxpt.X() < box.XMax()) maxpt[0] = box.XMax();
    if (maxpt.Y() < box.YMax()) maxpt[1] = box.YMax();
}



void R2Box::
Union (const R2Circle& circle) 
{
    // Expand this to include circle
    Union(circle.BBox());
}



void R2Box::
Intersect (const R2Box& box) 
{
    // Intersect with box
    if (minpt.X() < box.XMin()) minpt[0] = box.XMin();
    if (minpt.Y() < box.YMin()) minpt[1] = box.YMin();
    if (maxpt.X() > box.XMax()) maxpt[0] = box.XMax();
    if (maxpt.Y() > box.YMax()) maxpt[1] = box.YMax();
}



void R2Box::
Transform (const R2Transformation& transformation)
{
    // Do not transform empty box
    if (IsEmpty()) return;

    // Transform box ???
    R2Box tmp = R2null_box;
    for (RNQuadrant quadrant = 0; quadrant < RN_NUM_QUADRANTS; quadrant++) {
	R2Point corner(Corner(quadrant));
	corner.Transform(transformation);
	tmp.Union(corner);
    }
    *this = tmp;
}




void R2Box::
Reset(const R2Point& min, const R2Point& max)
{
    // Move box by vector
    minpt = min;
    maxpt = max;
}












