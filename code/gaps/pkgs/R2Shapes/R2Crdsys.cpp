/* Source file for the GAPS coordinatesystem class */



/* Include files */

#include "R2Shapes/R2Shapes.h"



/* Public variables */

const R2CoordSystem R2xy_coordinate_system(R2Point(0.0, 0.0), R2Diad(R2Vector(1.0, 0.0), R2Vector(0.0, 1.0)));



/* Public functions */

int R2InitCoordinateSystem()
{
    /* Return success */
    return TRUE;
}



void R2StopCoordinateSystem()
{
}



R2CoordSystem::
R2CoordSystem(void)
{
}



R2CoordSystem::
R2CoordSystem(const R2CoordSystem& cs)
    : origin(cs.origin),
      axes(cs.axes)
{
}



R2CoordSystem::
R2CoordSystem(const R2Point& origin, const R2Diad& axes)
    : origin(origin),
      axes(axes)
{
}



const R3Matrix R2CoordSystem::
Matrix(void) const
{
    // Return matrix (std -> cs)
    R3Matrix m(R3identity_matrix);
    m.Translate(origin.Vector());
    m.Transform(axes.Matrix());
    return m;
}



const R3Matrix R2CoordSystem::
InverseMatrix(void) const
{
    // Return matrix (cs -> std)
    R3Matrix m(axes.InverseMatrix());
    m.Translate(-origin.Vector());
    return m;
}



void R2CoordSystem::
SetOrigin(const R2Point& origin)
{
    // Set origin
    this->origin = origin;
}



void R2CoordSystem::
SetAxes(const R2Diad& axes)
{
    // Set axes
    this->axes = axes;
}



void R2CoordSystem::
Translate(const R2Vector& offset)
{
    // Translate origin
    origin.Translate(offset);
}



void R2CoordSystem::
Rotate(RNAngle radians)
{
    // Rotate axes only
    axes.Rotate(radians);
    // Normalize to avoid accumulation of error ???
    axes.Normalize();
}



void R2CoordSystem::
Mirror (const R2Line& line)
{
    // Mirror origin and axes
    origin.Mirror(line);
    axes.Mirror(line);
}



void R2CoordSystem::
Transform (const R2Transformation& transformation)
{
    // Transform origin and axes
    origin.Transform(transformation);
    axes.Transform(transformation);
    // Normalize to avoid accumulation of error ???
    axes.Normalize();
}



void R2CoordSystem::
InverseTransform (const R2Transformation& transformation)
{
    // Transform origin and axes
    origin.InverseTransform(transformation);
    axes.InverseTransform(transformation);
    // Normalize to avoid accumulation of error ???
    axes.Normalize();
}




