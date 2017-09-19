/* Source file for the R3 ellipsoid class */



/* Include files */

#include "R3Shapes/R3Shapes.h"



/* Class type definitions */

RN_CLASS_TYPE_DEFINITIONS(R3Ellipsoid);



/* Public functions */

int 
R3InitEllipsoid()
{
    /* Return success */
    return TRUE;
}



void 
R3StopEllipsoid()
{
}



R3Ellipsoid::
R3Ellipsoid(void)
{
}



R3Ellipsoid::
R3Ellipsoid(const R3Ellipsoid& ellipsoid)
    : cs(ellipsoid.cs),
      radii(ellipsoid.radii),
      bbox(ellipsoid.bbox)
{
}



R3Ellipsoid::
R3Ellipsoid(const R3CoordSystem& cs, const R3Vector& radii)
    : cs(cs),
      radii(radii)
{
  // Compute bounding box
  UpdateBBox();
}



const RNBoolean R3Ellipsoid::
IsPoint(void) const
{
    // A ellipsoid only lies on a single point if it has radii zero
    return ((radii.X() == 0.0) && (radii.Y() == 0.0) && (radii.Z() == 0.0));
}



const RNBoolean R3Ellipsoid::
IsLinear(void) const
{
    // A ellipsoid only lies on a single line if two radii are zero
    int nzeros = 0;
    if (radii.X() == 0.0) nzeros++;
    if (radii.Y() == 0.0) nzeros++;
    if (radii.Z() == 0.0) nzeros++;
    return (nzeros >= 2);
}



const RNBoolean R3Ellipsoid::
IsPlanar(void) const
{
    // A ellipsoid only lies on a plane if at least one radius is zero
    return ((radii.X() == 0.0) || (radii.Y() == 0.0) || (radii.Z() == 0.0));
}



const RNBoolean R3Ellipsoid::
IsConvex(void) const
{
    // All ellipsoids are convex
    return TRUE;
}



const RNInterval R3Ellipsoid::
NFacets(void) const
{
    // Return number of facets (polygons)
    // Min is stretched icosahedron, max is 32 stacks X 32 slices
    return RNInterval(20.0, 1024.0);
}



const RNArea R3Ellipsoid::
Area(void) const
{
    // Return surface area of ellipsoid ???
    return 0.0;
}



const RNVolume R3Ellipsoid::
Volume(void) const
{
    // Return volume of ellipsoid ???
    return 0.0;
}



const R3Point R3Ellipsoid::
Centroid(void) const
{
    // Return centroid of ellipsoid
    return Center();
}



const R3Shape& R3Ellipsoid::
BShape(void) const
{
    // Return bounding box;
    return bbox;
}



const R3Box R3Ellipsoid::
BBox(void) const
{
    // Return bounding box 
    return bbox;
}



const R3Sphere R3Ellipsoid::
BSphere(void) const
{
    // Return bounding sphere
    switch(radii.Sextant()) {
    case RN_NX_SEXTANT: return R3Sphere(Center(), radii.X());
    case RN_PX_SEXTANT: return R3Sphere(Center(), radii.X());
    case RN_NY_SEXTANT: return R3Sphere(Center(), radii.Y());
    case RN_PY_SEXTANT: return R3Sphere(Center(), radii.Y());
    case RN_NZ_SEXTANT: return R3Sphere(Center(), radii.Z());
    case RN_PZ_SEXTANT: return R3Sphere(Center(), radii.Z());
    default:
	RNAbort("Illegal sextant found");
	return R3null_sphere;
    }
}



void R3Ellipsoid::
Empty(void)
{
    // Empty ellipsoid
    cs = R3xyz_coordinate_system;
    radii = R3null_vector;
    bbox = R3null_box;
}



void R3Ellipsoid::
Translate(const R3Vector& vector)
{
    // Move ellipsoid coordinate system
    cs.Translate(vector);

    // Move bounding box
    bbox.Translate(vector);
}



void R3Ellipsoid::
Reposition(const R3Point& center)
{
    // Set ellipsoid coordinate system
    Translate(center - this->cs.Origin());
}



void R3Ellipsoid::
Resize(const R3Vector& radii) 
{
    // Set ellipsoid radii
    this->radii = radii;

    // Recompute bounding box
    UpdateBBox();
}


void R3Ellipsoid::
Transform (const R3Transformation& transformation)
{
    // Transform center and radii
    cs.Transform(transformation);
    radii *= transformation.ScaleFactor();

    // Recompute bounding box
    UpdateBBox();
}



void R3Ellipsoid::
UpdateBBox(void) 
{
    // Recompute bounding box
    bbox = BSphere().BBox();
}





