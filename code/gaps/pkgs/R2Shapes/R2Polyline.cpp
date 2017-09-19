/* Source file for the R2 polyline class */



/* Include files */

#include "R2Shapes/R2Shapes.h"



/* Public functions */

int 
R2InitPolyline()
{
  return TRUE;
}



void 
R2StopPolyline()
{
}



R2Polyline::
R2Polyline(void) 
  : points(NULL),
    npoints(0),
    bbox(R2null_box)
{
}



R2Polyline::
R2Polyline(const R2Polyline& polyline) 
  : points(NULL),
    npoints(polyline.npoints),
    bbox(polyline.bbox)
{
  // Copy points
  if (polyline.npoints > 0) {
    points = new R2Point [ npoints ];
    for (int i = 0; i < npoints; i++) {
      points[i] = polyline.points[i];
    }
  }
}



R2Polyline::
R2Polyline(const RNArray<R2Point *>& p)
  : points(NULL),
    npoints(p.NEntries()),
    bbox(R2null_box)
{
  // Copy points
  if (npoints > 0) {
    points = new R2Point [ npoints ];
    for (int i = 0; i < npoints; i++) {
      points[i] = *(p[i]);
      bbox.Union(*(p[i]));
    }
  }
}



R2Polyline::
R2Polyline(const R2Point *p, int np)
  : points(NULL),
    npoints(np),
    bbox(R2null_box)
{
  // Copy points
  if (npoints > 0) {
    points = new R2Point [ npoints ];
    for (int i = 0; i < npoints; i++) {
      points[i] = p[i];
      bbox.Union(p[i]);
    }
  }
}



R2Polyline::
~R2Polyline(void) 
{
  // Delete points
  if (points) delete [] points;
}



const R2Point R2Polyline::
ClosestPoint(const R2Point& point) const
{
  // Return closest point on polyline
  R2Point closest_point(0,0);
  RNLength closest_squared_distance = FLT_MAX;
  for (int i = 0; i < npoints; i++) {
    const R2Point& position = points[i];
    RNLength squared_distance = R2SquaredDistance(point, position);
    if (squared_distance < closest_squared_distance) {
      closest_point = position;
      closest_squared_distance = squared_distance;
    }
  }

  // Return closest point
  return closest_point;
}



const RNBoolean R2Polyline::
IsPoint(void) const
{
    // A polyline only lies on a single point if it has one point
    return (npoints == 1);
}



const RNBoolean R2Polyline::
IsLinear(void) const
{
    // A polyline only lies within a line if it has two points
    return (npoints == 2);
}



const RNLength R2Polyline::
Length(void) const
{
  // Check number of points
  if (npoints < 2) return 0;

  // Compute perimeter
  RNLength sum = 0;
  R2Point *p1 = &points[0];
  for (int i = 1; i < npoints; i++) {
    R2Point *p2 = &points[i];
    sum += R2Distance(*p1, *p2);
    p1 = p2;
  }

  // Return perimeter
  return sum;
}



const R2Point R2Polyline::
Centroid(void) const
{
  // Check number of points
  if (npoints == 0) return R2zero_point;

  // Return centroid
  R2Point sum = R2zero_point;
  for (int i = 0; i < npoints; i++) sum += points[i];
  return sum / npoints;
}



const R2Shape& R2Polyline::
BShape(void) const
{
    // Return bounding box
    return bbox;
}



const R2Box R2Polyline::
BBox(void) const
{
    // Return bounding box of polyline
    return bbox;
}



const R2Circle R2Polyline::
BCircle(void) const
{
    // Return bounding circle
    return bbox.BCircle();
}



R2Vector R2Polyline::
Normal(int k, RNLength radius) const
{
  // Just checking
  if (npoints < 3) return R2zero_vector;

  // Find first point
  R2Point *pointA = &points[(k - 1 + npoints) % npoints];
  if (radius > 0) {
    R2Point *prev = pointA;
    RNLength distance = R2Distance(*pointA, points[k]);
    if (distance < radius) {
      for (int i = 2;  i < npoints/3; i++) {
        pointA = &points[(k - i + npoints) % npoints];
        distance += R2Distance(*pointA, *prev);
        if (distance >= radius) break;
        prev = pointA;
      }
    }
  }

  // Find second point
  R2Point *pointB = &points[(k + 1) % npoints];
  if (radius > 0) {
    R2Point *prev = pointB;
    RNLength distance = R2Distance(*pointB, points[k]);
    if (distance < radius) {
      for (int i = 2;  i < npoints/3; i++) {
        pointB = &points[(k + i) % npoints];
        distance += R2Distance(*pointB, *prev);
        if (distance >= radius) break;
        prev = pointB;
      }
    }
  }

  // Return normal
  R2Vector vector = *pointB - *pointA;
  R2Vector n(vector[1], -vector[0]);
  n.Normalize();

  // Return normal
  return n;
}



R2Line R2Polyline::
Tangent(int k, RNLength radius) const
{
  // Just checking
  if (npoints < 3) return R2null_line;

  // Create array with points of patch within radius
  RNArray<R2Point *> patch;
  patch.Insert(&points[k]);

  // Add points searching one way
  RNLength distance = 0;
  R2Point *prev = &points[k];
  for (int i = 1; i < npoints/3; i++) {
    R2Point *current = &points[(k - i + npoints) % npoints];
    patch.Insert(current);
    distance += R2Distance(*current, *prev);
    if (distance >= radius) break;
    prev = current;
  }

  // Add points searching other way
  distance = 0;
  prev = &points[k];
  for (int i = 1; i < npoints/3; i++) {
    R2Point *current = &points[(k + i) % npoints];
    patch.Insert(current);
    distance += R2Distance(*current, *prev);
    if (distance >= radius) break;
    prev = current;
  }

  // Compute best fitting line to points in patch
  R2Point centroid = R2Centroid(patch);
  R2Diad diad = R2PrincipleAxes(centroid, patch);
  R2Vector vector = diad.Axis(0);

  // Flip vector if facing wrong direction
  R2Point& p1 = points[(k - 1 + npoints) % npoints];
  R2Point& p2 = points[(k + 1) % npoints];
  if (vector.Dot(p2 - p1) < 0) vector.Flip();

  // Return tangent line
  return R2Line(points[k], vector);
}



RNAngle R2Polyline::
InteriorAngle(int k, RNLength radius) const
{
  // Just checking
  if (npoints < 3) return 0;

  // Find first point
  R2Point *pointA = &points[(k - 1 + npoints) % npoints];
  if (radius > 0) {
    R2Point *prev = pointA;
    RNLength distance = R2Distance(*pointA, points[k]);
    if (distance < radius) {
      for (int i = 2;  i < npoints/3; i++) {
        pointA = &points[(k - i + npoints) % npoints];
        distance += R2Distance(*pointA, *prev);
        if (distance >= radius) break;
        prev = pointA;
      }
    }
  }

  // Find second point
  R2Point *pointB = &points[(k + 1) % npoints];
  if (radius > 0) {
    R2Point *prev = pointB;
    RNLength distance = R2Distance(*pointB, points[k]);
    if (distance < radius) {
      for (int i = 2;  i < npoints/3; i++) {
        pointB = &points[(k + i) % npoints];
        distance += R2Distance(*pointB, *prev);
        if (distance >= radius) break;
        prev = pointB;
      }
    }
  }

  // Compute interior angle
  R2Vector va = *pointA - points[k];
  R2Vector vb = *pointB - points[k];
  RNLength lena = va.Length();
  RNLength lenb = vb.Length();
  if (RNIsZero(lena) || RNIsZero(lenb)) return RN_UNKNOWN;
  RNAngle angle = R2InteriorAngle(va, vb);
  return angle;
}



RNScalar R2Polyline::
Curvature(int k, RNLength radius) const
{
  // Just checking
  if (npoints < 3) return 0;

  // Find first point
  R2Point *pointA = &points[(k - 1 + npoints) % npoints];
  if (radius > 0) {
    R2Point *prev = pointA;
    RNLength distance = R2Distance(*pointA, points[k]);
    if (distance < radius) {
      for (int i = 2;  i < npoints/3; i++) {
        pointA = &points[(k - i + npoints) % npoints];
        distance += R2Distance(*pointA, *prev);
        if (distance >= radius) break; 
        prev = pointA;
      }
    }
  }

  // Find second point
  R2Point *pointB = &points[(k + 1) % npoints];
  if (radius > 0) {
    R2Point *prev = pointB;
    RNLength distance = R2Distance(*pointB, points[k]);
    if (distance < radius) {
      for (int i = 2;  i < npoints/3; i++) {
        pointB = &points[(k + i) % npoints];
        distance += R2Distance(*pointB, *prev);
        if (distance >= radius) break;
        prev = pointB;
      }
    }
  }

  // Compute curvature
  R2Vector va = *pointA - points[k];
  R2Vector vb = *pointB - points[k];
  RNLength lena = va.Length();
  RNLength lenb = vb.Length();
  if (RNIsZero(lena) || RNIsZero(lenb)) return RN_UNKNOWN;
  RNScalar angle = R2InteriorAngle(va, vb);
  RNScalar curvature = fabs((RN_PI - angle) / (lena + lenb));

  // Return curvature
  return curvature;
}



void R2Polyline::
Empty(void)
{
    // Empty polyline
    if (points) delete [] points;
    points = NULL;
    npoints = 0;
    bbox = R2null_box;
}



void R2Polyline::
Transform(const R2Transformation& transformation) 
{
  // Transform points
  bbox = R2null_box;
  for (int i = 0; i < npoints; i++) {
    points[i].Transform(transformation);
    bbox.Union(points[i]);
  }
}



void R2Polyline::
Print(FILE *fp) const
{
  // Print points
  fprintf(fp, "%d\n", npoints);
  for (int i = 0; i < npoints; i++) {
    fprintf(fp, "%12.6f %12.6f\n", points[i].X(), points[i].Y());
  }
}



