/* Source file for the GAPS plane class */



/* Include files */

#include "R3Shapes/R3Shapes.h"



/* Public variables */

const R3Plane R3null_plane(0.0, 0.0, 0.0, 0.0);
const R3Plane R3posyz_plane(1.0, 0.0, 0.0, 0.0);
const R3Plane R3posxz_plane(0.0, 1.0, 0.0, 0.0);
const R3Plane R3posxy_plane(0.0, 0.0, 1.0, 0.0);
const R3Plane R3negyz_plane(-1.0, 0.0, 0.0, 0.0);
const R3Plane R3negxz_plane(0.0, -1.0, 0.0, 0.0);
const R3Plane R3negxy_plane(0.0, 0.0, -1.0, 0.0);



/* Public functions */

int 
R3InitPlane()
{
    /* Return success */
    return TRUE;
}



void 
R3StopPlane()
{
}



R3Plane::
R3Plane(void)
{
}



R3Plane::
R3Plane(const R3Plane& plane)
    : v(plane.v), 
      d(plane.d)
{
}



R3Plane::
R3Plane(RNScalar a, RNScalar b, RNScalar c, RNScalar d)
    : v(a, b, c), 
      d(d)
{
}



R3Plane::
R3Plane(const RNScalar a[4])
    : v(&a[0]), 
      d(a[3])
{
}



R3Plane::
R3Plane(const R3Vector& normal, RNScalar d)
    : v(normal),
      d(d)
{
}



R3Plane::
R3Plane(const R3Point& point, const R3Vector& normal)
{
    // Construct plane from point and normal vector
    v = normal;
    d = -(normal[0]*point[0] + normal[1]*point[1] + normal[2]*point[2]);
}



R3Plane::
R3Plane(const R3Point& point, const R3Line& line)
{
    // Construct plane through point and line
    v = point - line.Point();
    v.Cross(line.Vector());
    v.Normalize();
    d = -(v[0]*point[0] + v[1]*point[1] + v[2]*point[2]);
}



R3Plane::
R3Plane(const R3Point& point, const R3Vector& vector1, const R3Vector& vector2)
{
    // Construct plane through point and two vectors
    v = vector1 % vector2;
    v.Normalize();
    d = -(v[0]*point[0] + v[1]*point[1] + v[2]*point[2]);
}



R3Plane::
R3Plane(const R3Point& point1, const R3Point& point2, const R3Point& point3)
{
    // Construct plane through three points
    v = point2 - point1;
    R3Vector v3 = point3 - point1;
    v.Cross(v3);
    v.Normalize();
    d = -(v[0]*point1[0] + v[1]*point1[1] + v[2]*point1[2]);
}



R3Plane::
R3Plane(const RNArray<R3Point *>& points, RNBoolean polygon_vertices)
{
    // Initialize plane
    v = R3null_vector;
    d = 0;

    // Check number of points
    int npoints = points.NEntries();
    if (npoints < 3) {
      *this = R3null_plane;
      return;
    }

    // Compute centroid
    R3Point c = R3Centroid(points);

    // Check if points form (counter-clockwise) boundary of polygon
    if (polygon_vertices) {
        // Compute best normal for counter-clockwise array of points using newell's method 
        const R3Point *p1 = points[npoints-1];
        for (int i = 0; i < npoints; i++) {
            const R3Point *p2 = points[i];
            v[0] += (p1->Y() - p2->Y()) * (p1->Z() + p2->Z());
            v[1] += (p1->Z() - p2->Z()) * (p1->X() + p2->X());
            v[2] += (p1->X() - p2->X()) * (p1->Y() + p2->Y());
            p1 = p2;
        }

        // Normalize 
        v.Normalize();
    }
    else {
        // Compute principle axes
        R3Triad triad = R3PrincipleAxes(c, points);

        // Select direction of least variation
        v = triad[2];
    }

    // Compute d from centroid and normal
    d = -(v[0]*c[0] + v[1]*c[1] + v[2]*c[2]);
}



R3Plane::
R3Plane(R3Point *points, int npoints, RNBoolean polygon_vertices)
{
    // Initialize plane
    v = R3null_vector;
    d = 0;

    // Check number of points
    if (npoints < 3) {
      *this = R3null_plane;
      return;
    }

    // Compute centroid
    R3Point c = R3Centroid(npoints, points);

    // Check if points form (counter-clockwise) boundary of polygon
    if (polygon_vertices) {
        // Compute best normal for counter-clockwise array of points using newell's method 
        const R3Point *p1 = &points[npoints-1];
        for (int i = 0; i < npoints; i++) {
            const R3Point *p2 = &points[i];
            v[0] += (p1->Y() - p2->Y()) * (p1->Z() + p2->Z());
            v[1] += (p1->Z() - p2->Z()) * (p1->X() + p2->X());
            v[2] += (p1->X() - p2->X()) * (p1->Y() + p2->Y());
            p1 = p2;
        }

        // Normalize 
        v.Normalize();
    }
    else {
        // Compute principle axes
        R3Triad triad = R3PrincipleAxes(c, npoints, points);

        // Select direction of least variation
        v = triad[2];
    }

    // Compute d from centroid and normal
    d = -(v[0]*c[0] + v[1]*c[1] + v[2]*c[2]);
}



const R3Point R3Plane::
Point (void) const
{
    // Return point on plane
    return R3zero_point + v * -d;
}



void R3Plane::
Mirror(const R3Plane& plane)
{
    // Mirror plane ???
    R3Point p = Point();
    p.Mirror(plane);
    v.Mirror(plane);
    Reposition(p);
}



void R3Plane::
Reposition(const R3Point& point)
{
    // Move plane
    d = -(v[0]*point[0] + v[1]*point[1] + v[2]*point[2]);
}



void R3Plane::
Translate(const R3Vector& vector) 
{
    // Move plane by vector - there's got to be a better way ???
    Reposition(Point() + vector);
}



void R3Plane::
Transform (const R3Transformation& transformation)
{
    // Transform plane ???
    R3Point p = Point();
    p.Transform(transformation);
    transformation.ApplyInverseTranspose(v);
    v.Normalize();
    Reposition(p);
}



void R3Plane::
InverseTransform (const R3Transformation& transformation)
{
    // Transform plane ???
    R3Point p = Point();
    p.InverseTransform(transformation);
    transformation.ApplyTranspose(v);
    v.Normalize();
    Reposition(p);
}



void R3Plane::
Reset(const R3Point& point, const R3Vector& normal) 
{
    // Reset plane
    v = normal;
    d = -(normal[0]*point[0] + normal[1]*point[1] + normal[2]*point[2]);
}





