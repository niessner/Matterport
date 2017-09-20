/* Source file for the GAPS coordinatesystem class */



/* Include files */

#include "R3Shapes/R3Shapes.h"



/* Public variables */

const R3CoordSystem R3xyz_coordinate_system(R3Point(0.0, 0.0, 0.0), 
    R3Triad(R3Vector(1.0, 0.0, 0.0), R3Vector(0.0, 1.0, 0.0), R3Vector(0.0, 0.0, 1.0)));



/* Public functions */

int R3InitCoordinateSystem()
{
    /* Return success */
    return TRUE;
}



void R3StopCoordinateSystem()
{
}



R3CoordSystem::
R3CoordSystem(void)
{
}



R3CoordSystem::
R3CoordSystem(const R3CoordSystem& cs)
    : origin(cs.origin),
      axes(cs.axes)
{
}



R3CoordSystem::
R3CoordSystem(const R3Point& origin, const R3Triad& axes)
    : origin(origin),
      axes(axes)
{
}



const R4Matrix R3CoordSystem::
Matrix(void) const
{
    // Return matrix (std -> cs)
    R4Matrix m(R4identity_matrix);
    m.Translate(origin.Vector());
    m.Transform(axes.Matrix());
    return m;
}



const R4Matrix R3CoordSystem::
InverseMatrix(void) const
{
    // Return matrix (cs -> std)
    R4Matrix m(axes.InverseMatrix());
    m.Translate(-origin.Vector());
    return m;
}



void R3CoordSystem::
SetOrigin(const R3Point& origin)
{
    // Set origin
    this->origin = origin;
}



void R3CoordSystem::
SetAxes(const R3Triad& axes)
{
    // Set axes
    this->axes = axes;
}



void R3CoordSystem::
Reset(const R3Point& origin, const R3Triad& axes)
{
    // Set origin and axes
    this->origin = origin;
    this->axes = axes;
}



void R3CoordSystem::
Translate(const R3Vector& offset)
{
    // Translate origin
    origin.Translate(offset);
}



void R3CoordSystem::
Rotate(RNAxis axis, RNAngle radians)
{
    // Rotate axes only
    axes.Rotate(axis, radians);
    // Normalize to avoid accumulation of error ???
    axes.Normalize();
}



void R3CoordSystem::
Rotate(const R3Vector& axis, RNAngle radians)
{
    // Rotate axes only
    axes.Rotate(axis, radians);
    // Normalize to avoid accumulation of error ???
    axes.Normalize();
}



void R3CoordSystem::
Rotate(const R3Vector& from, const R3Vector& to)
{
    // Rotate axes only
    axes.Rotate(from, to);
    // Normalize to avoid accumulation of error ???
    axes.Normalize();
}



void R3CoordSystem::
Mirror (const R3Plane& plane)
{
    // Mirror origin and axes
    origin.Mirror(plane);
    axes.Mirror(plane);
}



void R3CoordSystem::
Transform (const R3Transformation& transformation)
{
    // Transform origin and axes
    origin.Transform(transformation);
    axes.Transform(transformation);
    // Normalize to avoid accumulation of error ???
    axes.Normalize();
}



void R3CoordSystem::
InverseTransform (const R3Transformation& transformation)
{
    // Transform origin and axes
    origin.InverseTransform(transformation);
    axes.InverseTransform(transformation);
    // Normalize to avoid accumulation of error ???
    axes.Normalize();
}




