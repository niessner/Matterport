/* Source file for the R3 ellipse class */



/* Include files */

#include "R3Shapes/R3Shapes.h"



/* Class type definitions */

RN_CLASS_TYPE_DEFINITIONS(R3Ellipse);



/* Public functions */

int 
R3InitEllipse()
{
    /* Return success */
    return TRUE;
}



void 
R3StopEllipse()
{
}



R3Ellipse::
R3Ellipse(void)
{
}



R3Ellipse::
R3Ellipse(const R3Ellipse& ellipse)
    : cs(ellipse.cs),
      radii(ellipse.radii),
      bbox(ellipse.bbox)
{
}



R3Ellipse::
R3Ellipse(const R3CoordSystem& cs, const R2Vector& radii)
    : cs(cs),
      radii(radii)
{
  // Compute bounding box
  UpdateBBox();
}



const RNBoolean R3Ellipse::
IsPoint(void) const
{
    // A ellipse only lies on a single point if it has radii zero
    return ((radii.X() == 0.0) && (radii.Y() == 0.0));
}



const RNBoolean R3Ellipse::
IsLinear(void) const
{
    // A ellipse only lies on a single line if two radii are zero
    int nzeros = 0;
    if (radii.X() == 0.0) nzeros++;
    if (radii.Y() == 0.0) nzeros++;
    return (nzeros >= 2);
}



const RNBoolean R3Ellipse::
IsPlanar(void) const
{
    // A ellipse only lies on a plane if at least one radius is zero
    return ((radii.X() == 0.0) || (radii.Y() == 0.0));
}



const RNBoolean R3Ellipse::
IsConvex(void) const
{
    // All ellipses are convex
    return TRUE;
}



const RNInterval R3Ellipse::
NFacets(void) const
{
    // Return number of facets (polygons)
    return RNInterval(1.0, 1.0);
}



const RNArea R3Ellipse::
Area(void) const
{
    // Return surface area of ellipse 
    return RN_PI * radii.X() * radii.Y();
}



const R3Point R3Ellipse::
Centroid(void) const
{
    // Return centroid of ellipse
    return Center();
}



const R3Shape& R3Ellipse::
BShape(void) const
{
    // Return bounding box;
    return bbox;
}



const R3Box R3Ellipse::
BBox(void) const
{
    // Return bounding box 
    return bbox;
}



const R3Sphere R3Ellipse::
BSphere(void) const
{
    // Return bounding sphere
    if (radii.X() > radii.Y()) return R3Sphere(Center(), radii.X());
    else return R3Sphere(Center(), radii.Y());
}



void R3Ellipse::
Empty(void)
{
    // Empty ellipse
    cs = R3xyz_coordinate_system;
    radii = R2null_vector;
    bbox = R3null_box;
}



void R3Ellipse::
Flip(void)
{
    // Flip ellipse
    R3Triad triad(cs.Axes().Axis(RN_X), cs.Axes().Axis(RN_Y), -(cs.Axes().Axis(RN_Z)));
    cs.SetAxes(triad);
}



void R3Ellipse::
Translate(const R3Vector& vector)
{
    // Move ellipse coordinate system
    cs.Translate(vector);

    // Move bounding box
    bbox.Translate(vector);
}



void R3Ellipse::
Reposition(const R3Point& center)
{
    // Set ellipse coordinate system
    Translate(center - this->cs.Origin());
}



void R3Ellipse::
Resize(const R2Vector& radii) 
{
    // Set ellipse radii
    this->radii = radii;

    // Recompute bounding box
    UpdateBBox();
}


void R3Ellipse::
Transform (const R3Transformation& transformation)
{
    // Transform center and radii
    cs.Transform(transformation);
    radii *= transformation.ScaleFactor();

    // Recompute bounding box
    UpdateBBox();
}



void R3Ellipse::
UpdateBBox(void) 
{
    // Recompute bounding box
    bbox = BSphere().BBox();
}









