/* Source file for the R2 polygon class */



/* Include files */

#include "R2Shapes/R2Shapes.h"



/* Public functions */

int 
R2InitPolygon()
{
  return TRUE;
}



void 
R2StopPolygon()
{
}



R2Polygon::
R2Polygon(void) 
  : points(NULL),
    npoints(0),
    bbox(R2null_box),
    clockwise(FALSE)
{
}



R2Polygon::
R2Polygon(const R2Polygon& polygon) 
  : points(NULL),
    npoints(polygon.npoints),
    bbox(polygon.bbox),
    clockwise(polygon.clockwise)
{
  // Copy points
  if (polygon.npoints > 0) {
    points = new R2Point [ npoints ];
    for (int i = 0; i < npoints; i++) {
      points[i] = polygon.points[i];
    }
  }
}



R2Polygon::
R2Polygon(const RNArray<R2Point *>& p, RNBoolean clockwise)
  : points(NULL),
    npoints(p.NEntries()),
    bbox(R2null_box),
    clockwise(clockwise)
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



R2Polygon::
R2Polygon(const R2Point *p, int np, RNBoolean clockwise)
  : points(NULL),
    npoints(np),
    bbox(R2null_box),
    clockwise(clockwise)
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



R2Polygon::
~R2Polygon(void) 
{
  // Delete points
  if (points) delete [] points;
}



const R2Point R2Polygon::
ClosestPoint(const R2Point& point) const
{
  // Return closest point on polygon
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



const RNBoolean R2Polygon::
IsPoint(void) const
{
    // A polygon only lies on a single point if it has one point
    return (npoints == 1);
}



const RNBoolean R2Polygon::
IsLinear(void) const
{
    // A polygon only lies within a line if it has two points
    return (npoints == 2);
}



const RNBoolean R2Polygon::
IsConvex(void) const
{
    RNAbort("Not implemented");
    return FALSE;
}



const RNArea R2Polygon::
Area(void) const
{
  // Check number of points
  if (npoints < 3) return 0;

  // Compute twicearea by sum of cross products
  RNLength sum = 0;
  R2Point *p1 = &points[npoints-1];
  for (int i = 0; i < npoints; i++) {
    R2Point *p2 = &points[i];
    sum += p1->X()*p2->Y() - p2->X()*p1->Y();
    p1 = p2;
  }

  // Compute area
  RNArea area = 0.5 * sum;

  // Flip if clockwise
  if (clockwise) area = -area;

  // Return area
  return area;
}



const RNArea R2Polygon::
Convexity(void) const
{
  // Check number of points
  if (npoints < 3) return 0;

  // Compute area
  RNArea area = Area();

  // Compute area of convex hull
  R2Polygon convex_hull = ConvexHull();
  RNArea convex_hull_area = convex_hull.Area();
  if (convex_hull_area == 0) return 0;

  // Return ratio of areas
  return area / convex_hull_area;
}



const RNLength R2Polygon::
Perimeter(void) const
{
  // Check number of points
  if (npoints < 2) return 0;

  // Compute perimeter
  RNLength sum = 0;
  R2Point *p1 = &points[npoints-1];
  for (int i = 0; i < npoints; i++) {
    R2Point *p2 = &points[i];
    sum += R2Distance(*p1, *p2);
    p1 = p2;
  }

  // Return perimeter
  return sum;
}



const R2Point R2Polygon::
Centroid(void) const
{
  // Check number of points
  if (npoints == 0) return R2zero_point;

  // Return centroid
  R2Point sum = R2zero_point;
  for (int i = 0; i < npoints; i++) sum += points[i];
  return sum / npoints;
}



const R2Shape& R2Polygon::
BShape(void) const
{
    // Return bounding box
    return bbox;
}



const R2Box R2Polygon::
BBox(void) const
{
    // Return bounding box of polygon
    return bbox;
}



const R2Circle R2Polygon::
BCircle(void) const
{
    // Return bounding circle
    return bbox.BCircle();
}



const R2Polygon R2Polygon::
ConvexHull(void) const
{
  // Check number of points
  if (npoints <= 3) return *this;

  // Create convex hull polygon
  R2Polygon convex_hull;
  convex_hull.bbox = R2null_box;
  convex_hull.clockwise = FALSE;
  convex_hull.points = new R2Point [ npoints ];
  convex_hull.npoints = 0;

  // Find index of leftmost point
  int leftmost_index = 0;
  for (int i = 1; i < npoints; i++) {
    if (points[i][0] < points[leftmost_index][0]) {
      leftmost_index = i;
    }
  }

  // Iteratively find next vertex for convex hull via gift-wrapping
  int index = leftmost_index;
  for (int i = 0; i < npoints; i++) {
    // Insert convex hull point
    convex_hull.points[convex_hull.npoints++] = points[index];
    convex_hull.bbox.Union(points[index]);
    R2Point last_point = points[index];
    int last_index = index;

    // Find next giftwrapping point index
    index = (last_index+1)%npoints;
    R2Line line(last_point, points[index]);
    for (int j = 0; j < npoints; j++) {
      if (j == last_index) continue;
      if (j == leftmost_index) {
        RNScalar d = R2SignedDistance(line, points[j]) + 1000 * RN_EPSILON;
        if (d > 0) {
          line = R2Line(last_point, points[j]);
          index = j;
        }
      }
      else {
        if (R2Contains(last_point, points[j])) continue;
        RNScalar d = R2SignedDistance(line, points[j]);
        if (d > 0) {
          line = R2Line(last_point, points[j]);
          index = j;
        }
      }
    }

    // Check if completed loop
    if (index == leftmost_index) break;
  }

  // Return the convex polygon
  return convex_hull;
}



R2Vector R2Polygon::
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

  // Flip if polygon is clockwise
  if (clockwise) n = -n;

  // Return normal
  return n;
}



R2Line R2Polygon::
Tangent(int k, RNLength radius) const
{
  // Just checking
  if (npoints < 3) return R2null_line;

#if 0
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

  // Compute vector along line
  R2Vector vector = *pointB - *pointA;
  vector.Normalize();
#else
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
#endif

  // Flip if polygon is clockwise
  if (clockwise) vector = -vector;

  // Return tangent line
  return R2Line(points[k], vector);
}



RNAngle R2Polygon::
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
  RNScalar cross = (va[0]*vb[1] - va[1]*vb[0]);
  if (clockwise) { if (cross < 0) angle = RN_TWO_PI - angle; }
  else { if (cross > 0) angle = RN_TWO_PI - angle; }
  return angle;
}



RNScalar R2Polygon::
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
  RNScalar cross = (va[0]*vb[1] - va[1]*vb[0]);
  RNScalar sign = (cross > 0) ? -1 : 1;
  RNScalar curvature = sign * (RN_PI - angle) / (lena + lenb);

  // Flip if polygon is clockwise
  if (clockwise) curvature = -curvature;

  // Return curvature
  return curvature;
}



void R2Polygon::
Empty(void)
{
    // Empty polygon
    if (points) delete [] points;
    points = NULL;
    npoints = 0;
    bbox = R2null_box;
}



void R2Polygon::
Clip(const R2Line& line) 
{
  // Check number of points
  if (npoints == 0) {
    return;
  }
  else if (npoints == 1) {
    if (R2SignedDistance(line, points[0]) < 0) {
      bbox = R2null_box;
      delete [] points;
      points = NULL;
      npoints = 0;
    }
  }
  else {
    // Check bounding box
    R2Halfspace halfspace(line, 0);

    // Check if bounding box is entirely on positive side of line
    if (R2Contains(halfspace, bbox)) {
      return;
    }

    // Check if bounding box is entirely on negative side of line
    if (R2Contains(-halfspace, bbox)) {
      bbox = R2null_box;
      delete [] points;
      points = NULL;
      npoints = 0;
      return;
    }

    // Create new array for points
    int nbuffer = 0;
    R2Point *buffer = new R2Point [ 4 * npoints ];
    if (!buffer) RNAbort("Unable to allocate points during clip");

    // Build buffer with clipped points
    const R2Point *p1 = &points[npoints-1];
    RNScalar d1 = R2SignedDistance(line, *p1);
    for (int i = 0; i < npoints; i++) {
      const R2Point *p2 = &points[i];
      RNScalar d2 = R2SignedDistance(line, *p2);
      if (d2 >= 0) {
        // Insert crossing from negative to positive
        if (d1 < 0) {
          RNScalar t = d2 / (d2 - d1);
          buffer[nbuffer++] = *p2 + t * (*p1 - *p2);
        }

        // Insert point on positive side
        buffer[nbuffer++] = *p2;
      }
      else {
        // Insert crossing from positive to negative
        if (d1 >= 0) {
          RNScalar t = d1 / (d1 - d2);
          buffer[nbuffer++] = *p1 + t * (*p2 - *p1);
        }
      }

      // Remember previous point
      p1 = p2;
      d1 = d2;
    }

    // Copy points
    bbox = R2null_box;
    npoints = nbuffer;
    delete [] points;
    points = new R2Point [ npoints ];
    for (int i = 0; i < npoints; i++) {
      points[i] = buffer[i];
      bbox.Union(points[i]);
    }

    // Delete the buffer of points
    delete [] buffer;
  }
}



void R2Polygon::
Clip(const R2Box& box) 
{
  // Clip to each side of box
  if (npoints == 0) return;
  Clip(R2Line(1, 0, -(box.XMin())));
  if (npoints == 0) return;
  Clip(R2Line(-1, 0, box.XMax()));
  if (npoints == 0) return;
  Clip(R2Line(0, 1, -(box.YMin())));
  if (npoints == 0) return;
  Clip(R2Line(0, -1, box.YMax()));
}



void R2Polygon::
Transform(const R2Transformation& transformation) 
{
  // Transform points
  bbox = R2null_box;
  for (int i = 0; i < npoints; i++) {
    points[i].Transform(transformation);
    bbox.Union(points[i]);
  }
}



void R2Polygon::
Print(FILE *fp) const
{
  // Print points
  fprintf(fp, "%d\n", npoints);
  for (int i = 0; i < npoints; i++) {
    fprintf(fp, "%12.6f %12.6f\n", points[i].X(), points[i].Y());
  }
}



int R2Polygon::
ReadTheraFile(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open polygon file: %s\n", filename);
    return 0;
  }

  // Read header line
  int ndimensions;
  double polygon_depth, sample_spacing;
  if (fscanf(fp, "%d %d %lf %lf", &ndimensions, &npoints, &polygon_depth, &sample_spacing) != 4) {
    fprintf(stderr, "Unable to read polygon file header: %s\n", filename);
    return 0;
  }

  // Allocate points
  points = new R2Point [ npoints ];
  if (!points) {
    fprintf(stderr, "Unable to allocate points for polygon: %s\n", filename);
    return 0;
  }

  // Read points
  for (int i = 0; i < npoints; i++) {
    double x, y;
    if (fscanf(fp, "%lf %lf", &x, &y) != 2) {
      fprintf(stderr, "Unable to read point %d from polygon file %s\n", npoints, filename);
      return 0;
    }

    // Add point
    points[i].Reset(x,y);
    bbox.Union(points[i]);
  }

  // Close file
  fclose(fp);

  // Return success
  return npoints;
}











