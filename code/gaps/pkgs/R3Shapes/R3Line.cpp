/* Source file for the GAPS line class */



/* Include files */

#include "R3Shapes/R3Shapes.h"



/* Public variables */

const R3Line R3null_line(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
const R3Line R3posx_line(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
const R3Line R3posy_line(0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
const R3Line R3posz_line(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
const R3Line R3negx_line(0.0, 0.0, 0.0, -1.0, 0.0, 0.0);
const R3Line R3negy_line(0.0, 0.0, 0.0, 0.0, -1.0, 0.0);
const R3Line R3negz_line(0.0, 0.0, 0.0, 0.0, 0.0, -1.0);



/* Public functions */

int R3InitLine()
{
    /* Return success */
    return TRUE;
}



void R3StopLine()
{
}



R3Line::
R3Line(void)
{
}



R3Line::
R3Line(const R3Line& line)
    : point(line.point),
      vector(line.vector)
{
}



R3Line::
R3Line(const R3Point& point, const R3Vector& vector, RNBoolean normalized)
    : point(point),
      vector(vector)
{
    // Normalize vector
    if (!normalized) this->vector.Normalize();
}



R3Line::
R3Line(const R3Point& point1, const R3Point& point2)
    : point(point1),
      vector(point2 - point1)
{
    // Normalize vector
    vector.Normalize();
}



R3Line::
R3Line(RNCoord x1, RNCoord y1, RNCoord z1, RNCoord x2, RNCoord y2, RNCoord z2)
    : point(x1, y1, z1),
      vector(x2-x1, y2-y1, z2-z1)
{
    // Normalize vector
    vector.Normalize();
}



R3Line::
R3Line(const RNArray<R3Point *>& points)
    : point(0, 0, 0),
      vector(0, 0, 0)
{
    // Check number of points
    if (points.NEntries() < 2) return;

    // Compute centroid
    point = R3Centroid(points);

    // Compute principle axes
    R3Triad triad = R3PrincipleAxes(point, points);

    // Select direction of most variation
    vector = triad[0];
}



R3Line::
R3Line(R3Point *points, int npoints)
    : point(0, 0, 0),
      vector(0, 0, 0)
{
    // Check number of points
    if (npoints < 2) return;

    // Compute centroid
    point = R3Centroid(npoints, points);

    // Compute principle axes
    R3Triad triad = R3PrincipleAxes(point, npoints, points);

    // Select direction of most variation
    vector = triad[0];
}



void R3Line::
Transform (const R3Transformation& transformation)
{
    // Transform point and vector
    point.Transform(transformation);
    vector.Transform(transformation);
    vector.Normalize();
}



void R3Line::
InverseTransform (const R3Transformation& transformation)
{
    // Transform point and vector
    point.InverseTransform(transformation);
    vector.InverseTransform(transformation);
    vector.Normalize();
}



const RNBoolean R3Line::
operator==(const R3Line& line) const
{
    // Check if vectors are equal
    if (vector != line.vector) return FALSE;

    // Return whether point lies on this line
    R3Vector v = vector;
    v.Cross(line.Point() - point);
    return v.IsZero();
}




