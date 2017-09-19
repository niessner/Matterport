/* Source file for the R3 Box class */



/* Include files */

#include "R3Shapes/R3Shapes.h"



/* Public variables */

const R3Box R3zero_box(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
const R3Box R3unit_box(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0);
const R3Box R3infinite_box(-RN_INFINITY, -RN_INFINITY, -RN_INFINITY,
			    RN_INFINITY, RN_INFINITY, RN_INFINITY);
const R3Box R3null_box( FLT_MAX, FLT_MAX, FLT_MAX, 
		       -FLT_MAX, -FLT_MAX, -FLT_MAX);



/* Class type definitions */

RN_CLASS_TYPE_DEFINITIONS(R3Box);



/* Public functions */

int 
R3InitBox()
{
    /* Return success */
    return TRUE;
}



void 
R3StopBox()
{
}



R3Box::
R3Box(void)
{
}



R3Box::
R3Box(const R3Box& box)
    : minpt(box.minpt),
      maxpt(box.maxpt)
{
}



R3Box::
R3Box(const R3Point& minpt, const R3Point& maxpt)
    : minpt(minpt), 
      maxpt(maxpt)
{
}



R3Box::
R3Box(RNCoord xmin, RNCoord ymin, RNCoord zmin,
      RNCoord xmax, RNCoord ymax, RNCoord zmax)
    : minpt(xmin, ymin, zmin),
      maxpt(xmax, ymax, zmax)
{
}



const RNBoolean R3Box::
IsPoint (void) const
{
    // Return whether bounding box contains just one point
    return ((minpt.X() == maxpt.X()) && (minpt.Y() == maxpt.Y()) && (minpt.Z() == maxpt.Z()));
}



const RNBoolean R3Box::
IsLinear (void) const
{
    // Return whether bounding box contains just one line segment
    return (this->NDimensions() <= 1);
}



const RNBoolean R3Box::
IsPlanar (void) const
{
    // Return whether bounding box contains just one rectangle
    return (this->NDimensions() <= 2);
}



const RNBoolean R3Box::
IsConvex (void) const
{
    // All boxes are convex
    return TRUE;
}



const RNBoolean R3Box::
IsFinite (void) const
{
    // Return whether bounding box contains a finite amount of space
    return ((!IsEmpty()) && 
	    (RNIsFinite(minpt.X()) && RNIsFinite(minpt.Y()) && RNIsFinite(minpt.Z())) && 
	    (RNIsFinite(maxpt.X()) && RNIsFinite(maxpt.Y()) && RNIsFinite(maxpt.Z())));
}



const R3Point R3Box::
Corner (RNOctant octant) const
{
    // Return corner point 
    switch (octant) {
    case RN_NNN_OCTANT: return minpt;
    case RN_NNP_OCTANT: return R3Point(minpt.X(), minpt.Y(), maxpt.Z());
    case RN_NPN_OCTANT: return R3Point(minpt.X(), maxpt.Y(), minpt.Z());
    case RN_NPP_OCTANT: return R3Point(minpt.X(), maxpt.Y(), maxpt.Z());
    case RN_PNN_OCTANT: return R3Point(maxpt.X(), minpt.Y(), minpt.Z());
    case RN_PNP_OCTANT: return R3Point(maxpt.X(), minpt.Y(), maxpt.Z());
    case RN_PPN_OCTANT: return R3Point(maxpt.X(), maxpt.Y(), minpt.Z());
    case RN_PPP_OCTANT: return maxpt;
    default: RNAbort("Invalid octant for box corner"); return R3null_point;
    }
}



const R3Point R3Box::
Corner (RNDirection xdir, RNDirection ydir, RNDirection zdir) const
{
    // Return corner point 
    return R3Point(Coord(xdir, RN_X), Coord(ydir, RN_Y), Coord(zdir, RN_Z));
}



const R3Plane R3Box::
Plane (RNDirection dir, RNDimension dim) const
{
    // Return plane along side
    R3Vector v = R3xyz_triad[dim]; 
    RNCoord d = (*this)[dir][dim];
    if (dir == RN_LO) return R3Plane(-v, d);
    else return R3Plane(v, -d);
}



const R3Box R3Box::
Side (RNDirection dir, RNDimension dim) const
{
    // Return box along side
    R3Box result = *this;
    result[1-dir][dim] = result[dir][dim];
    return result;
}



const R3Box R3Box::
Octant (RNDirection xdir, RNDirection ydir, RNDirection zdir) const
{
    // Return box in octant (xdir, ydir, zdir)
    R3Box result;
    R3Point centroid = Centroid();
    if (xdir == RN_HI) { result[RN_LO][RN_X] = centroid[RN_X]; result[RN_HI][RN_X] = maxpt[RN_X]; }
    else { result[RN_LO][RN_X] = minpt[RN_X]; result[RN_HI][RN_X] = centroid[RN_X]; }
    if (ydir == RN_HI) { result[RN_LO][RN_Y] = centroid[RN_Y]; result[RN_HI][RN_Y] = maxpt[RN_Y]; }
    else { result[RN_LO][RN_Y] = minpt[RN_Y]; result[RN_HI][RN_Y] = centroid[RN_Y]; }
    if (zdir == RN_HI) { result[RN_LO][RN_Z] = centroid[RN_Z]; result[RN_HI][RN_Z] = maxpt[RN_Z]; }
    else { result[RN_LO][RN_Z] = minpt[RN_Z]; result[RN_HI][RN_Z] = centroid[RN_Z]; }
    return result;
}



const int R3Box::
NDimensions (void) const
{
    // Return number of non-zero dimensions in box
    int ndimensions = 0;
    if (minpt.X() > maxpt.X()) return -1;
    if (minpt.X() < maxpt.X()) ndimensions++;
    if (minpt.Y() > maxpt.Y()) return -1;
    if (minpt.Y() < maxpt.Y()) ndimensions++;
    if (minpt.Z() > maxpt.Z()) return -1;
    if (minpt.Z() < maxpt.Z()) ndimensions++;
    return ndimensions;
}



const int R3Box::
NDimensionsAlongSide (RNDirection, RNDimension dim) const
{
    // Return number of non-zero dimensions along side
    int ndimensions = 0;
    int dim1 = (dim+1) % 3;
    int dim2 = (dim+2) % 3;
    if (minpt[dim1] < maxpt[dim1]) ndimensions++;
    if (minpt[dim2] < maxpt[dim2]) ndimensions++;
    return ndimensions;
}



const int R3Box::
NDimensionsAlongAxis (const RNAxis axis) const
{
    // Return number of non-zero dimensions along axis
    if (minpt[axis] <= maxpt[axis]) return 0;
    else return 1;
}



const RNAxis R3Box::
ShortestAxis (void) const
{
    // Compute length of each axis
    RNLength dx = this->XLength();
    RNLength dy = this->YLength();
    RNLength dz = this->ZLength();

    // Return shortest axis
    if (dx < dy) {
	if (dx < dz) return RN_XAXIS;
	else return RN_ZAXIS;
    }
    else {
	if (dy < dz) return RN_YAXIS;
	else return RN_ZAXIS;
    }    
}



const RNAxis R3Box::
LongestAxis (void) const
{
    // Compute length of each axis
    RNLength dx = this->XLength();
    RNLength dy = this->YLength();
    RNLength dz = this->ZLength();

    // Return longest axis
    if (dx > dy) {
	if (dx > dz) return RN_XAXIS;
	else return RN_ZAXIS;
    }
    else {
	if (dy > dz) return RN_YAXIS;
	else return RN_ZAXIS;
    }    
}



const RNLength R3Box::
ShortestAxisLength (void) const
{
    // Compute length of each axis
    RNLength dx = this->XLength();
    RNLength dy = this->YLength();
    RNLength dz = this->ZLength();

    // Return length of shortest axis
    if (dx < dy) {
	if (dx < dz) return dx;
	else return dz;
    }
    else {
	if (dy < dz) return dy;
	else return dz;
    }    
}



const RNLength R3Box::
LongestAxisLength (void) const
{
    // Compute length of each axis
    RNLength dx = this->XLength();
    RNLength dy = this->YLength();
    RNLength dz = this->ZLength();

    // Return longest axis
    if (dx > dy) {
	if (dx > dz) return dx;
	else return dz;
    }
    else {
	if (dy > dz) return dy;
	else return dz;
    }    
}



const RNLength R3Box::
DiagonalLength (void) const
{
    // Return length of Box along diagonal
    return R3Distance(minpt, maxpt);
}



const R3Point R3Box::
ClosestPoint(const R3Point& point) const
{
    // Return closest point in box
    R3Point closest(point);
    if (closest.X() < XMin()) closest[RN_X] = XMin();
    else if (closest.X() > XMax()) closest[RN_X] = XMax();
    if (closest.Y() < YMin()) closest[RN_Y] = YMin();
    else if (closest.Y() > YMax()) closest[RN_Y] = YMax();
    if (closest.Z() < ZMin()) closest[RN_Z] = ZMin();
    else if (closest.Z() > ZMax()) closest[RN_Z] = ZMax();
    return closest;
}



const R3Point R3Box::
FurthestPoint(const R3Point& point) const
{
    // Return furthest point in box
    R3Point furthest(point);
    RNScalar x = (point.X() < XCenter()) ? XMax() : XMin();
    RNScalar y = (point.Y() < YCenter()) ? YMax() : YMin();
    RNScalar z = (point.Z() < ZCenter()) ? ZMax() : ZMin();
    return R3Point(x, y, z);
}



const RNInterval R3Box::
NFacets(void) const
{
    // Return number of facets (polygons)
    return RNInterval(6.0, 6.0);
}



const RNArea R3Box::
Area (void) const
{
    // Return surface area of bounding box
    RNLength dx = XLength();
    RNLength dy = YLength();
    RNLength dz = ZLength();
    return (2.0*dx*dy + 2.0*dx*dz + 2.0*dy*dz);
}



const RNVolume R3Box::
Volume (void) const
{
    // Return volume of bounding box
    return (XLength() * YLength() * ZLength());
}



const R3Point R3Box::
Centroid (void) const
{
    // Return center point
    return R3Point(XCenter(), YCenter(), ZCenter());
}



const R3Shape& R3Box::
BShape (void) const
{
    // Return self
    return *this;
}



const R3Box R3Box::
BBox (void) const
{
    // Return self - shape virtual function
    return *this;
}



const R3Sphere R3Box::
BSphere (void) const
{
    // Return bounding sphere 
    return R3Sphere(Centroid(), DiagonalRadius());
}



void R3Box::
Empty (void) 
{
    // Copy empty box 
    *this = R3null_box;
}



void R3Box::
Inflate (RNScalar fraction) 
{
  // Scale box around centroid by fraction
  if (IsEmpty()) return;
  R3Point centroid = Centroid();
  R3Vector diagonal = Max() - centroid;
  Reset(centroid - fraction * diagonal, centroid + fraction * diagonal);
}



void R3Box::
Translate(const R3Vector& vector)
{
    // Move box by vector
    minpt += vector;
    maxpt += vector;
}



void R3Box::
Union (const R3Point& point) 
{
    // Expand this to include point
    if (IsEmpty()) {
        minpt = point;
        maxpt = point;
    }
    else {
        if (minpt.X() > point.X()) minpt[0] = point.X();
        if (minpt.Y() > point.Y()) minpt[1] = point.Y();
        if (minpt.Z() > point.Z()) minpt[2] = point.Z();
        if (maxpt.X() < point.X()) maxpt[0] = point.X();
        if (maxpt.Y() < point.Y()) maxpt[1] = point.Y();
        if (maxpt.Z() < point.Z()) maxpt[2] = point.Z();
    }
}



void R3Box::
Union (const R3Box& box) 
{
    // Expand this to include box
    if (IsEmpty()) {
        *this = box;
    }
    else if (!box.IsEmpty()) {
        if (minpt.X() > box.XMin()) minpt[0] = box.XMin();
        if (minpt.Y() > box.YMin()) minpt[1] = box.YMin();
        if (minpt.Z() > box.ZMin()) minpt[2] = box.ZMin();
        if (maxpt.X() < box.XMax()) maxpt[0] = box.XMax();
        if (maxpt.Y() < box.YMax()) maxpt[1] = box.YMax();
        if (maxpt.Z() < box.ZMax()) maxpt[2] = box.ZMax();
    }
}



void R3Box::
Union (const R3Sphere& sphere) 
{
    // Expand this to include sphere
    Union(sphere.BBox());
}



void R3Box::
Intersect (const R3Box& box) 
{
    // Intersect with box
    if (!IsEmpty()) {
        if (box.IsEmpty()) {
            *this = box;
        }
        else {
            if (minpt.X() < box.XMin()) minpt[0] = box.XMin();
            if (minpt.Y() < box.YMin()) minpt[1] = box.YMin();
            if (minpt.Z() < box.ZMin()) minpt[2] = box.ZMin();
            if (maxpt.X() > box.XMax()) maxpt[0] = box.XMax();
            if (maxpt.Y() > box.YMax()) maxpt[1] = box.YMax();
            if (maxpt.Z() > box.ZMax()) maxpt[2] = box.ZMax();
        }
    }
}



void R3Box::
Transform (const R3Transformation& transformation)
{
    // Do not transform empty box
    if (IsEmpty()) return;

    // Transform box ???
    R3Box tmp = R3null_box;
    for (RNOctant octant = 0; octant < RN_NUM_OCTANTS; octant++) {
	R3Point corner(Corner(octant));
	corner.Transform(transformation);
	tmp.Union(corner);
    }
    *this = tmp;
}



void R3Box::
Reset(const R3Point& min, const R3Point& max)
{
    // Move box by vector
    minpt = min;
    maxpt = max;
}





