////////////////////////////////////////////////////////////////////////
// Source file for RGBDSegmentation class
////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "RGBD.h"



////////////////////////////////////////////////////////////////////////
// Parameters
////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////
// Point member functions
////////////////////////////////////////////////////////////////////////

RGBDPoint::
RGBDPoint(void)
  : depth(-1),
    position(0,0,0),
    normal(0,0,0),
    radius(0),
    boundary(0),
    neighbors(),
    segment(NULL),
    segment_affinity(0),
    segment_index(-1),
    grid_index(-1),
    mark(0)
{
}



////////////////////////////////////////////////////////////////////////
// Primitive member functions
////////////////////////////////////////////////////////////////////////

RGBDPrimitive::
RGBDPrimitive(int primitive_type)
  : primitive_type(primitive_type),
    bbox(R3null_box),
    centroid(R3zero_point),
    line(R3null_line),
    plane(R3null_plane)
{
}



RGBDPrimitive::
RGBDPrimitive(const RGBDPrimitive& primitive)
  : primitive_type(primitive.primitive_type),
    bbox(primitive.bbox),
    centroid(primitive.centroid),
    line(primitive.line),
    plane(primitive.plane)
{
}



RGBDPrimitive::
RGBDPrimitive(RGBDPoint *seed_point, const RNArray<RGBDPoint *> *points)
  : primitive_type(RGBD_NULL_PRIMITIVE_TYPE),
    bbox(R3null_box),
    centroid(R3zero_point),
    line(R3null_line),
    plane(R3null_plane)
{
  // Initialize primitive based on points
  Update(seed_point, points); 
}



RNLength RGBDPrimitive::
Distance(const R3Point& position) const
{
  // Return distance from primitive to point
  if (primitive_type == RGBD_POINT_PRIMITIVE_TYPE) return R3Distance(centroid, position);
  else if (primitive_type == RGBD_LINE_PRIMITIVE_TYPE) return R3Distance(line, position);
  else if (primitive_type == RGBD_PLANE_PRIMITIVE_TYPE) return R3Distance(plane, position);
#if 0
  else if (primitive_type == RGBD_PLANAR_GRID_PRIMITIVE_TYPE) {
    R2Point grid_position = planar_grid.GridPosition(position);
    int ix = (int) (grid_position.X() + 0.5);
    if ((ix < 0) || (ix >= planar_grid.XResolution())) return RN_INFINITY;
    int iy = (int) (grid_position.Y() + 0.5);
    if ((iy < 0) || (iy >= planar_grid.YResolution())) return RN_INFINITY;
    RNScalar value = planar_grid.GridValue(ix, iy);
    if (value <= 0) return RN_INFINITY;
    return R3Distance(planar_grid.Plane(), position);
  }
#endif
  else {
    RNAbort("Unrecognized primitive type");
    return RN_INFINITY;
  }
}



void RGBDPrimitive::
Update(const R3Point& point)
{
  // Set everything
  primitive_type = RGBD_POINT_PRIMITIVE_TYPE;
  this->centroid = point;
  line = R3null_line;
  plane = R3null_plane;
}



void RGBDPrimitive::
Update(const R3Line& line)
{
  // Set everything
  primitive_type = RGBD_LINE_PRIMITIVE_TYPE;
  centroid = R3zero_point;
  centroid.Project(line);
  this->line = line;
  plane = R3null_plane;
}



void RGBDPrimitive::
Update(const R3Plane& plane)
{
  // Set everything
  primitive_type = RGBD_PLANE_PRIMITIVE_TYPE;
  centroid = R3zero_point;
  centroid.Project(plane);
  line = R3null_line;
  this->plane = plane;
}



void RGBDPrimitive::
Update(RGBDPoint *seed_point, const RNArray<RGBDPoint *> *points)
{
  // Remember stuff about primitive (so can set same orientation)
  R3Vector previous_vector(0,0,0);
  if (primitive_type == RGBD_LINE_PRIMITIVE_TYPE) previous_vector = line.Vector();
  if (primitive_type == RGBD_PLANE_PRIMITIVE_TYPE) previous_vector = plane.Normal();

  // Update bounding box
  bbox = R3null_box;
  if (points) {
    for (int i = 0; i < points->NEntries(); i++) {
      RGBDPoint *point = points->Kth(i);
      bbox.Union(point->position);
    }
  }

  // Initialize everything
  if (seed_point) {
    R3Point seed_position = seed_point->position;
    centroid = seed_position;
    line.Reset(seed_position, line.Vector()); 
    plane.Reset(seed_position, seed_point->normal);
    bbox.Union(seed_position);
  }
  else {
    // Temporary
    RNAbort("Need seed point");
  }

  // Update based on points
  if (points && (points->NEntries() > 0)) {
    // Allocate arrays of point positions and weights
    const int max_positions = 1024;
    R3Point *positions = new R3Point [ max_positions ];
    RNScalar *weights = new RNScalar [ max_positions ];

    // Fill arrays of point positions and weights
    int npositions = 0;
    int skip = points->NEntries() / max_positions + 1;
    for (int i = 0; i < points->NEntries(); i += skip) {
      RGBDPoint *point = points->Kth(i);
      if (npositions >= max_positions-1) break;
      positions[npositions] = point->position;
      if (primitive_type == RGBD_PLANE_PRIMITIVE_TYPE) weights[npositions] = fabs(plane.Normal().Dot(point->normal));
      else if (primitive_type == RGBD_LINE_PRIMITIVE_TYPE) weights[npositions] = 1.0 - fabs(plane.Normal().Dot(point->normal));
      else weights[npositions] = 1.0;
      npositions++;
    }

    // Add seed point with 20% of the total weight
    if (seed_point) {
      positions[npositions] = seed_point->position;
      weights[npositions] = 0.2 * points->NEntries();
      npositions++;
    }

    // Compute centroid
    centroid = R3Centroid(npositions, positions, weights);

    // Update primitive parameters
    if ((primitive_type == RGBD_NULL_PRIMITIVE_TYPE) && (npositions >= 2)) {
      RNScalar variances[3];
      R3Triad axes = R3PrincipleAxes(centroid, npositions, positions, weights, variances);
      if (variances[0] > RN_EPSILON) {
        if (variances[1] > RN_EPSILON) {
          RNScalar ratio10 = variances[1] / variances[0];
          RNScalar ratio21 = variances[2] / variances[1];
          if (ratio10 < ratio21) {
            primitive_type = RGBD_LINE_PRIMITIVE_TYPE;
            line.Reset(centroid, axes[0]);
          }
          else {
            primitive_type = RGBD_PLANE_PRIMITIVE_TYPE;
            plane.Reset(centroid, axes[2]);
          }
        }
      }
    }
    else if ((primitive_type == RGBD_LINE_PRIMITIVE_TYPE) && (npositions >= 2)) {
      // Compute principle directions
      RNScalar variances[3];
      R3Triad axes = R3PrincipleAxes(centroid, npositions, positions, weights, variances);
      if (variances[0] > RN_EPSILON) {
        // Update line
        R3Vector direction = axes[0];
        line.Reset(centroid, direction);

        // Check if should flip line
        RNScalar dot = direction.Dot(previous_vector);
        if (dot < 0) line.Flip();
      }
    }
    else if ((primitive_type == RGBD_PLANE_PRIMITIVE_TYPE) && (npositions >= 3)) {
      // Compute principle directions
      RNScalar variances[3];
      R3Triad axes = R3PrincipleAxes(centroid, npositions, positions, weights, variances);
      if (variances[1] > RN_EPSILON) {
        // Update plane
        R3Vector normal = axes[2];
        plane.Reset(centroid, normal);

        // Check if should flip plane
        if (seed_point) {
          RNScalar dot = normal.Dot(seed_point->normal);
          if (dot < 0) plane.Flip();
        }
        else {
          RNScalar dot = normal.Dot(previous_vector);
          if (dot < 0) plane.Flip();
        }
      }
    }

#if 0
    // Rasterize planar grid
    if (primitive_type == RGBD_PLANAR_GRID_PRIMITIVE_TYPE) {
      // Rasterize points into planar grid
      R3PlanarGrid density(plane, bbox, min_segment_spacing);
      for (int i = 0; i < points->NEntries(); i++) {
        RGBDPoint *point = points->Kth(i);
        R3Point position = point->position;
        density.RasterizeWorldPoint(position.X(), position.Y(), position.Z(), 1.0);
      }

      // Reset planar grid
      planar_grid.Reset(plane, bbox, min_segment_spacing);
      R3Point seed_position = seed_point->position;
      R2Point grid_position = planar_grid.GridPosition(seed_position);
      int ix = (int) (grid_position.X() + 0.5);
      int iy = (int) (grid_position.Y() + 0.5);
      FloodCopy(density.grid, planar_grid.grid, ix, iy);
    }
#endif

    // Delete array of point positions and weights
    delete [] positions;
    delete [] weights;
  }
}



void RGBDPrimitive::
Update(RGBDPrimitive primitive1, RGBDPrimitive primitive2, RNScalar weight1, RNScalar weight2)
{
  // Just checking
  if (weight1 == 0) {
    primitive_type = primitive2.primitive_type;
    bbox = primitive2.bbox;
    centroid = primitive2.centroid;
    line = primitive2.line;
    plane = primitive2.plane;
  }
  else if (weight2 == 0) {
    primitive_type = primitive1.primitive_type;
    bbox = primitive1.bbox;
    centroid = primitive1.centroid;
    line = primitive1.line;
    plane = primitive1.plane;
  }
  else {
    // Update primitive type
    if (primitive1.primitive_type > primitive2.primitive_type) {
      primitive_type = primitive1.primitive_type;
      weight2 = 0;
    }
    else if (primitive2.primitive_type > primitive1.primitive_type) {
      primitive_type = primitive2.primitive_type;
      weight1 = 0;
    }
    else {
      primitive_type = primitive1.primitive_type;
    }

    // Update centroid
    centroid = R3zero_point;
    centroid += weight1 * primitive1.centroid;
    centroid += weight2 * primitive2.centroid;
    centroid /= weight1 + weight2;

    // Update bbox
    bbox = R3null_box;
    bbox.Union(primitive1.bbox);
    bbox.Union(primitive2.bbox);

    // Update other stuff
    line = R3null_line;
    plane = R3null_plane;
    if (primitive_type == RGBD_LINE_PRIMITIVE_TYPE) {
      // Compute line
      R3Vector vector1 = primitive1.line.Vector();
      R3Vector vector2 = primitive2.line.Vector();
      if (vector1.Dot(vector2) < 0) vector2.Flip();
      R3Vector vector = R3zero_vector;
      vector += weight1 * vector1;
      vector += weight2 * vector2;
      vector /= weight1 + weight2;
      vector.Normalize();
      line.Reset(centroid, vector);
    }
    else if (primitive_type == RGBD_PLANE_PRIMITIVE_TYPE) {
      // Compute plane
      R3Vector normal1 = primitive1.plane.Normal();
      R3Vector normal2 = primitive2.plane.Normal();
      if (normal1.Dot(normal2) < 0) normal2.Flip();
      R3Vector normal = R3zero_vector;
      normal += weight1 * normal1;
      normal += weight2 * normal2;
      normal /= weight1 + weight2;
      normal.Normalize();
      plane.Reset(centroid, normal);
    }
  }
}



////////////////////////////////////////////////////////////////////////
// Segment member functions
////////////////////////////////////////////////////////////////////////

RGBDSegment::
RGBDSegment(RGBDSegmentation *segmentation, RGBDPoint *seed_point, int primitive_type)
  : segmentation(NULL),
    segmentation_index(-1),
    seed_point(seed_point),
    points(),
    parent(NULL),
    children(),
    pairs(),
    primitive(primitive_type),
    possible_affinity(0),
    total_affinity(0)
{
  // Update primitive
  if (seed_point) primitive.Update(seed_point);

  // Insert into segmentation
  if (segmentation) {
    this->segmentation_index = segmentation->segments.NEntries();
    this->segmentation = segmentation;
    segmentation->segments.Insert(this);
  }
}



RGBDSegment::
RGBDSegment(RGBDSegmentation *segmentation, RGBDPoint *seed_point, const RGBDPrimitive& primitive)
  : segmentation(NULL),
    segmentation_index(-1),
    seed_point(seed_point),
    points(),
    parent(NULL),
    children(),
    pairs(),
    primitive(primitive),
    possible_affinity(0),
    total_affinity(0)
{
  // Insert into segmentation
  if (segmentation) {
    this->segmentation_index = segmentation->segments.NEntries();
    this->segmentation = segmentation;
    segmentation->segments.Insert(this);
  }
}



RGBDSegment::
RGBDSegment(RGBDSegmentation *segmentation, RGBDSegment *child1, RGBDSegment *child2)
  : segmentation(NULL),
    segmentation_index(-1),
    seed_point(NULL),
    points(),
    parent(NULL),
    children(),
    pairs(),
    primitive(),
    possible_affinity(0),
    total_affinity(0)
{
  // Assign seed point
  seed_point = child1->seed_point;

  // Update primitive
  primitive.Update(child1->primitive, child2->primitive, child1->points.NEntries(), child2->points.NEntries());

  // Insert points from child1
  while (!child1->points.IsEmpty()) {
    RGBDPoint *point = child1->points.Tail();
    child1->RemovePoint(point);
    RNScalar affinity = Affinity(point);
    // RNScalar affinity = point->segment_affinity; // THIS IS WRONG, USING AFFINITY TO OLD SEGMENT FOR SPEED
    if (affinity < 0) affinity = 0;
    possible_affinity += affinity;
    InsertPoint(point, affinity);
  }

  // Insert points from child2
  while (!child2->points.IsEmpty()) {
    RGBDPoint *point = child2->points.Tail();
    child2->RemovePoint(point);
    RNScalar affinity = Affinity(point);
    // RNScalar affinity = point->segment_affinity; // THIS IS WRONG, USING AFFINITY TO OLD SEGMENT FOR SPEED
    if (affinity < 0) affinity = 0;
    possible_affinity += affinity;
    InsertPoint(point, affinity);
  }

  // Update hierarchy
  child1->parent = this;
  child2->parent = this;
  children.Insert(child1);
  children.Insert(child2);

  // Insert into segmentation
  if (segmentation) {
    this->segmentation_index = segmentation->segments.NEntries();
    this->segmentation = segmentation;
    segmentation->segments.Insert(this);
  }
}



RGBDSegment::
~RGBDSegment(void)
{
  // Delete children
  // for (int i = 0; i < children.NEntries(); i++) {
  //   delete children.Kth(i);
  // }

  // Remove from parent
  if (parent) parent->children.Remove(this);

  // Empty points
  EmptyPoints();
}



RNScalar RGBDSegment::
Coverage(void)
{
  // Return metric of how well segment covers points
  if (possible_affinity == 0) return 0;
  return total_affinity / possible_affinity;
}



void RGBDSegment::
EmptyPoints(void)
{
  // Update points
  for (int i = 0; i < points.NEntries(); i++) {
    RGBDPoint *point = points.Kth(i);
    point->segment = NULL;
    point->segment_affinity = 0;
    point->segment_index = -1;
  }

  // Empty points
  points.Empty();

  // Update affinity
  total_affinity = 0;
}



void RGBDSegment::
InsertPoint(RGBDPoint *point, RNScalar affinity)
{
  // Remove from previous segment
  if (point->segment ) {
    if (point->segment == this) return;
    else point->segment->RemovePoint(point);
  }

  // Update point
  point->segment = this;
  point->segment_index = points.NEntries();
  point->segment_affinity = affinity;

  // Insert point
  points.Insert(point);

  // Update segment
  total_affinity += point->segment_affinity;
}



void RGBDSegment::
RemovePoint(RGBDPoint *point)
{
  // Just checking
  assert(point->segment == this);
  assert(point->segment_index >= 0);

  // Update segment
  total_affinity -= point->segment_affinity;

  // Remove point
  RNArrayEntry *entry = points.KthEntry(point->segment_index);
  RGBDPoint *tail = points.Tail();
  tail->segment_index = point->segment_index;
  points.EntryContents(entry) = tail;
  points.RemoveTail();

  // Update point
  point->segment = NULL;
  point->segment_index = -1;
  point->segment_affinity = 0;
}



void RGBDSegment::
InsertChild(RGBDSegment *child)
{
  // Update primitive
  primitive.Update(this->primitive, child->primitive, this->points.NEntries(), child->points.NEntries());

  // Update affinities for current points
  if (points.NEntries() < 4 * child->points.NEntries()) {
    for (int i = 0; i < points.NEntries(); i++) {
      RGBDPoint *point = points.Kth(i);
      RNScalar affinity = Affinity(point);
      if (affinity < 0) affinity = 0;
      possible_affinity += affinity - point->segment_affinity;
      point->segment_affinity = affinity;
    }
  }

  // Insert points from child
  while (!child->points.IsEmpty()) {
    RGBDPoint *point = child->points.Tail();
    child->RemovePoint(point);
    RNScalar affinity = Affinity(point);
    if (affinity < 0) affinity = 0;
    possible_affinity += affinity;
    InsertPoint(point, affinity);
  }

  // Update hierarchy
  child->parent = this;
  children.Insert(child);
}



void RGBDSegment::
RemoveChild(RGBDSegment *child)
{
  // Remove child
  this->children.Remove(child);
  child->parent = NULL;
}



int RGBDSegment::
UpdatePoints(const R3Kdtree<RGBDPoint *> *kdtree)
{
  // Empty points
  // If do this, some points may end up as part of no segment
  // If don't do this, some segments may end up with outlier points if primitive changes a lot
  // EmptyPoints();  

  // Find points near primitive
  RNArray<RGBDPoint *> points1;
  if (seed_point) {
    // Find connected set of points near primitive
    static int mark = 1;
    RNArray<RGBDPoint *> stack;
    stack.Insert(seed_point);
    seed_point->mark = ++mark;
    while (!stack.IsEmpty()) {
      RGBDPoint *point = stack.Tail();
      stack.RemoveTail();
      points1.Insert(point);
      for (int i = 0; i < point->neighbors.NEntries(); i++) {
        RGBDPoint *neighbor = point->neighbors.Kth(i);
        if (neighbor->mark == mark) continue;
        neighbor->mark = mark;
        RNLength d = primitive.Distance(neighbor->position);
        RNLength max_segment_primitive_distance = (segmentation) ? segmentation->max_segment_primitive_distance : 0;
        if (d > max_segment_primitive_distance) continue;
        stack.Insert(neighbor);
      }
    }
  }
  else if (kdtree) {
    // Find all points near primitive
    RNLength max_segment_primitive_distance = (segmentation) ? segmentation->max_segment_primitive_distance : 0;
    if (primitive.primitive_type == RGBD_POINT_PRIMITIVE_TYPE) kdtree->FindAll(primitive.centroid, 0, max_segment_primitive_distance, points1);
    else if (primitive.primitive_type == RGBD_LINE_PRIMITIVE_TYPE) kdtree->FindAll(primitive.line, 0, max_segment_primitive_distance, points1);
    else if (primitive.primitive_type == RGBD_PLANE_PRIMITIVE_TYPE) kdtree->FindAll(primitive.plane, 0, max_segment_primitive_distance, points1);
    else RNAbort("Unrecognized primitive type");
  }

  // Check points
  int min_segment_points = (segmentation) ? segmentation->min_segment_points : 0;
  if ((min_segment_points > 0) && (points1.NEntries() < min_segment_points)) {
    return 0;
  }
  
  // Allocate affinities
  RNScalar *affinities = new RNScalar [ points1.NEntries() ];
  if (!affinities) return 0;

  // Compute affinities
  possible_affinity = 0;
  RNArray<RGBDPoint *> points2;
  for (int i = 0; i < points1.NEntries(); i++) {
    RGBDPoint *point = points1.Kth(i);
    RNScalar affinity = Affinity(point);
    if (affinity <= 0) continue;
    affinities[points2.NEntries()] = affinity;
    points2.Insert(point);
    possible_affinity += affinity;
  }

  // Check points
  if ((min_segment_points > 0) && (points2.NEntries() < min_segment_points)) {
    delete [] affinities;
    return 0;
  }

  // Compute assignments
  int npoints = 0;
  for (int i = 0; i < points2.NEntries(); i++) {
    RGBDPoint *point = points2.Kth(i);
    if ((point->segment) && (point->segment != this)) {
      if (point->segment->possible_affinity > 0) {
        RNScalar factor = possible_affinity / point->segment->possible_affinity;
        if (point->segment_affinity > factor * factor * affinities[i]) continue;
      }
    }
    npoints++;
  }

  // Check assignments
  if ((min_segment_points > 0) && (npoints < min_segment_points)) {
    delete [] affinities;
    return 0;
  }

  // Insert points (should match previous loop)
  for (int i = 0; i < points2.NEntries(); i++) {
    RGBDPoint *point = points2.Kth(i);
    if ((point->segment) && (point->segment != this)) {
      if (point->segment->possible_affinity > 0) {
        RNScalar factor = possible_affinity / point->segment->possible_affinity;
        if (point->segment_affinity > factor * factor * affinities[i]) continue;
      }
    }
    InsertPoint(point, affinities[i]);
  }

  // Delete affinities
  delete [] affinities;

  // Return success
  return 1;
}



int RGBDSegment::
UpdatePrimitive(void)
{
  // Update primitive
  primitive.Update(seed_point, &points);
  if (primitive.primitive_type == RGBD_NULL_PRIMITIVE_TYPE) return 0;
  else return 1;
}



RNScalar RGBDSegment::
Affinity(RGBDPoint *point) const
{
  // Initialize affinity
  RNScalar affinity = 1.0;

  // Get useful variables
  R3Point position = point->position;

  // Check primitive distance 
  RNLength max_segment_primitive_distance = (segmentation) ? segmentation->max_segment_primitive_distance : 0;
  if (max_segment_primitive_distance > 0) {
    RNLength primitive_distance = primitive.Distance(position);
    if (primitive_distance > max_segment_primitive_distance) return 0;
    RNScalar primitive_distance_affinity = 1.0 - primitive_distance / max_segment_primitive_distance;
    affinity *= primitive_distance_affinity;
  }

  // Check centroid distance
  RNLength max_segment_diameter = (segmentation) ? segmentation->max_segment_diameter : 0;
  if (max_segment_diameter > 0) {
    RNLength centroid_distance = R3Distance(primitive.centroid, position);
    if (centroid_distance > max_segment_diameter) return 0;
    RNScalar centroid_distance_affinity = 1.0 - centroid_distance / max_segment_diameter;
    affinity *= centroid_distance_affinity;
  }

  // Check normal angle
  RNAngle max_segment_normal_angle = (segmentation) ? segmentation->max_segment_normal_angle : 0;
  if (max_segment_normal_angle > 0) {
    if (primitive.primitive_type == RGBD_LINE_PRIMITIVE_TYPE) {
      RNScalar dot = fabs(primitive.line.Vector().Dot(point->normal));
      RNAngle normal_angle = (dot < 1) ? RN_PI_OVER_TWO - acos(dot) : RN_PI_OVER_TWO;
      if (normal_angle > max_segment_normal_angle) return 0;
      RNScalar normal_angle_affinity = 1.0 - normal_angle / max_segment_normal_angle;
      affinity *= normal_angle_affinity;
    }
    else if (primitive.primitive_type == RGBD_PLANE_PRIMITIVE_TYPE) {
      RNScalar dot = fabs(primitive.plane.Normal().Dot(point->normal));
      RNAngle normal_angle = (dot < 1) ? acos(dot) : 0;
      if (normal_angle > max_segment_normal_angle) return 0;
      RNScalar normal_angle_affinity = 1.0 - normal_angle / max_segment_normal_angle;
      affinity *= normal_angle_affinity;
    }
  }

  // Return affinity
  return affinity;
}



RNScalar RGBDSegment::
Affinity(RGBDSegment *segment) const
{
  // Initialize affinity
  RNScalar affinity = 1;

  // Compute centroid distance
  RNLength max_pair_centroid_distance = (segmentation) ? segmentation->max_pair_centroid_distance : 0;
  if (max_pair_centroid_distance > 0) {
    RNLength segment_distance = R3Distance(primitive.centroid, segment->primitive.centroid);
    if (segment_distance > max_pair_centroid_distance) return 0;
    RNScalar segment_distance_factor = 1.0 - segment_distance / max_pair_centroid_distance;
    affinity *= segment_distance_factor;
  }

  // Compute primitive distances
  RNLength max_pair_primitive_distance = (segmentation) ? segmentation->max_pair_primitive_distance : 0;
  if (max_pair_primitive_distance > 0) {
    // Compute point0-primitive1 distance
    RNLength primitive0_distance = primitive.Distance(segment->primitive.centroid);
    if (primitive0_distance > max_pair_primitive_distance) return 0;
    RNScalar primitive0_distance_factor = 1.0 - primitive0_distance / max_pair_primitive_distance;
    affinity *= primitive0_distance_factor;

    // Compute point1-primitive0 distance
    RNLength primitive1_distance = segment->primitive.Distance(primitive.centroid);
    if (primitive1_distance > max_pair_primitive_distance) return 0;
    RNScalar primitive1_distance_factor = 1.0 - primitive1_distance / max_pair_primitive_distance;
    affinity *= primitive1_distance_factor;
  }

  // Compute normal angle
  RNLength max_pair_normal_angle = (segmentation) ? segmentation->max_pair_normal_angle : 0;
  if (max_pair_normal_angle > 0) {
    if ((primitive.primitive_type == RGBD_LINE_PRIMITIVE_TYPE) && (segment->primitive.primitive_type == RGBD_LINE_PRIMITIVE_TYPE)) {
      RNScalar dot = fabs(primitive.line.Vector().Dot(segment->primitive.line.Vector()));
      RNAngle normal_angle = (dot < 1) ? acos(dot) : 0;
      RNAngle max_segment_normal_angle = (segmentation) ? segmentation->max_segment_normal_angle : 0;
      if ((max_segment_normal_angle > 0) && (normal_angle > max_segment_normal_angle)) return 0;
      RNScalar normal_angle_affinity = 1.0 - normal_angle / max_segment_normal_angle;
      affinity *= normal_angle_affinity;
    }
    else if ((primitive.primitive_type == RGBD_PLANE_PRIMITIVE_TYPE) && (segment->primitive.primitive_type == RGBD_LINE_PRIMITIVE_TYPE)) {
      RNScalar dot = fabs(primitive.plane.Normal().Dot(segment->primitive.line.Vector()));
      RNAngle normal_angle = (dot < 1) ? RN_PI_OVER_TWO - acos(dot) : RN_PI_OVER_TWO;
      if (normal_angle > max_pair_normal_angle) return 0;
      RNScalar normal_angle_factor = 1.0 - normal_angle / max_pair_normal_angle;
      affinity *= normal_angle_factor;
    }
    else if ((primitive.primitive_type == RGBD_LINE_PRIMITIVE_TYPE) && (segment->primitive.primitive_type == RGBD_PLANE_PRIMITIVE_TYPE)) {
      RNScalar dot = fabs(primitive.line.Vector().Dot(segment->primitive.plane.Normal()));
      RNAngle normal_angle = (dot < 1) ? RN_PI_OVER_TWO - acos(dot) : RN_PI_OVER_TWO;
      if (normal_angle > max_pair_normal_angle) return 0;
      RNScalar normal_angle_factor = 1.0 - normal_angle / max_pair_normal_angle;
      affinity *= normal_angle_factor;
    }
    else if ((primitive.primitive_type == RGBD_PLANE_PRIMITIVE_TYPE) && (segment->primitive.primitive_type == RGBD_PLANE_PRIMITIVE_TYPE)) {
      RNScalar dot = fabs(primitive.plane.Normal().Dot(segment->primitive.plane.Normal()));
      RNAngle normal_angle = (dot < 1) ? acos(dot) : 0;
      if (normal_angle > max_pair_normal_angle) return 0;
      RNScalar normal_angle_factor = 1.0 - normal_angle / max_pair_normal_angle;
      affinity *= normal_angle_factor;
    }
  }

#if 0
  // Compute imbalance
  double min_imbalance = 0.01;
  double max_imbalance = 0;
  if ((min_imbalance > 0) || (max_imbalance > 0)) {
    int npoints1 = (points.NEntries() < segment->points.NEntries()) ? points.NEntries() : segment->points.NEntries();
    int npoints2 = (points.NEntries() > segment->points.NEntries()) ? points.NEntries() : segment->points.NEntries();
    if (npoints2 == 0) return 0;
    RNScalar imbalance = (double) npoints1 / (double) npoints2;
    if ((max_imbalance > 0) && (imbalance > max_imbalance)) return 0;
    if (imbalance < min_imbalance) imbalance = min_imbalance;
    affinity *= imbalance;
  }
#endif
  
  // Return affinity
  return affinity;
}



static int
RGBDCompareSegments(const void *data1, const void *data2)
{
  RGBDSegment *segment1 = *((RGBDSegment **) data1);
  RGBDSegment *segment2 = *((RGBDSegment **) data2);
  if (segment2->total_affinity > segment1->total_affinity) return 1;
  else if (segment1->total_affinity > segment2->total_affinity) return -1;
  else return 0;
}



////////////////////////////////////////////////////////////////////////
// RGBDSegmentPair member functions
////////////////////////////////////////////////////////////////////////

RGBDSegmentPair::
RGBDSegmentPair(RGBDSegment *segment1, RGBDSegment *segment2, RNScalar affinity)
  : affinity(affinity),
    heapentry(NULL)
{
  // Insert pair into segments
  if (segment1 && segment2) {
    // Remember segments
    segments[0] = segment1;
    segments[1] = segment2;

    // Remember position of pair in segments
    segment_index[0] = segment1->pairs.NEntries();
    segment_index[1] = segment2->pairs.NEntries();

    // Update segments
    segment1->pairs.Insert(this);
    segment2->pairs.Insert(this);
  }
  else {
    // Initialize segments
    segments[0] = NULL;
    segments[1] = NULL;

    // Initialize segment index
    segment_index[0] = -1;
    segment_index[1] = -1;
  }
}



RGBDSegmentPair::
~RGBDSegmentPair(void)
{
  // Remove this pair from first segment
  if (segments[0]) {
    assert(segment_index[0] >= 0);
    RNArrayEntry *entry = segments[0]->pairs.KthEntry(segment_index[0]);
    RGBDSegmentPair *tail = segments[0]->pairs.Tail();
    if (tail->segments[0] == segments[0]) tail->segment_index[0] = segment_index[0];
    else if (tail->segments[1] == segments[0]) tail->segment_index[1] = segment_index[0];
    segments[0]->pairs.EntryContents(entry) = tail;
    segments[0]->pairs.RemoveTail();
  }

  // Remove this pair from second segment
  if (segments[1]) {
    assert(segment_index[1] >= 0);
    RNArrayEntry *entry = segments[1]->pairs.KthEntry(segment_index[1]);
    RGBDSegmentPair *tail = segments[1]->pairs.Tail();
    if (tail->segments[0] == segments[1]) tail->segment_index[0] = segment_index[1];
    else if (tail->segments[1] == segments[1]) tail->segment_index[1] = segment_index[1];
    segments[1]->pairs.EntryContents(entry) = tail;
    segments[1]->pairs.RemoveTail();
  }
}



static RGBDSegmentPair *
FindPair(RGBDSegment *segment1, RGBDSegment *segment2) 
{
  // Swap segments so that segment1 has fewer pairs
  if (segment1->pairs.NEntries() > segment2->pairs.NEntries()) {
    RGBDSegment *swap = segment1; 
    segment1 = segment2; 
    segment2 = swap;
  }

  // Search for pair
  for (int i = 0; i < segment1->pairs.NEntries(); i++) {
    RGBDSegmentPair *pair = segment1->pairs.Kth(i);
    if (pair->segments[0] == segment2) return pair;
    if (pair->segments[1] == segment2) return pair;
  }

  // RGBDSegmentPair not found
  return NULL;
}


#if 0
static int
FindPairIndex(RNArray<RGBDSegmentPair *>& pairs, RGBDSegment *segment0, RGBDSegment *segment1)
{
  // Search for primitive marker
  for (int i = 0; i < pairs.NEntries(); i++) {
    RGBDSegmentPair *pair = pairs.Kth(i);
    if ((segment0 == pair->segments[0]) && (segment1 == pair->segments[1])) return i;
    else if ((segment0 == pair->segments[1]) && (segment1 == pair->segments[0])) return i;
  }

  // Not found
  return -1;
}
#endif



////////////////////////////////////////////////////////////////////////
// Segmentation functions
////////////////////////////////////////////////////////////////////////

RGBDSegmentation::
RGBDSegmentation(void)
  : points(),
    kdtree(NULL),
    segments(),
    point_buffer(NULL)
{
  // Set default parameters
  SetDefaultParameters();
}



RGBDSegmentation::
RGBDSegmentation(const R2Grid& px_image, const R2Grid& py_image, const R2Grid& pz_image,
  const R2Grid& nx_image, const R2Grid& ny_image, const R2Grid& nz_image,
  const R2Grid& depth_image, const R2Grid& radius_image, 
  const R2Grid& boundary_image, const R2Image& color_image,
  const R3Point& viewpoint, const R3Vector& towards, const R3Vector& up,
  int primitive_type)
  : points(),
    kdtree(NULL),
    segments(),
    point_buffer(NULL)
{
  // Set parameters 
  SetDefaultParameters();
  min_segment_points = 10 * depth_image.NEntries() / (640 * 480);
  
  // Create points
  CreatePoints(px_image, py_image, pz_image, nx_image, ny_image, nz_image,
    depth_image, radius_image, boundary_image, color_image);

  // Create segments
  CreateSegments(primitive_type);
}



RGBDSegmentation::
~RGBDSegmentation(void)
{
  // Delete segments
  for (int i = 0; i < segments.NEntries(); i++) delete segments[i];
  
  // Delete kdtree
  if (kdtree) delete kdtree;

  // Delete points
  if (point_buffer) delete [] point_buffer;
  else { for (int i = 0; i < points.NEntries(); i++) delete points[i]; }
}



RNScalar RGBDSegmentation::
Affinity(void) const
{
  RNScalar sum = 0;
  for (int i = 0; i < segments.NEntries(); i++) {
    sum += segments[i]->total_affinity;
  }
  return sum;
}



int RGBDSegmentation::
NUnsegmentedPoints(void) const
{
  // Count unsegmented points
  int count = 0;
  for (int i = 0; i < points.NEntries(); i++) {
    RGBDPoint *point = points.Kth(i);
    if (!point->segment) count++;
  }

  // Return number of unsegmented points
  return count;
}



int RGBDSegmentation::
CreatePoints(const R2Grid& px_image, const R2Grid& py_image, const R2Grid& pz_image, 
  const R2Grid& nx_image, const R2Grid& ny_image, const R2Grid& nz_image,
  const R2Grid& depth_image, const R2Grid& radius_image,
  const R2Grid& boundary_image, const R2Image& color_image)
{
  // Allocate points 
  point_buffer = new RGBDPoint [ depth_image.NEntries() ];
  if (!point_buffer) {
    fprintf(stderr, "Unable to allocate points\n");
    return 0;
  }

  // Fill points
  for (int ix = 0; ix < depth_image.XResolution(); ix++) {
    for (int iy = 0; iy < depth_image.YResolution(); iy++) {
      int i;
      depth_image.IndicesToIndex(ix, iy, i);
      RGBDPoint *point = &point_buffer[i];

      // Check depth
      RNScalar depth = depth_image.GridValue(i);
      if (RNIsNegativeOrZero(depth)) continue;
      point->depth = depth;

      // Get position
      RNScalar px = px_image.GridValue(i);
      RNScalar py = py_image.GridValue(i);
      RNScalar pz = pz_image.GridValue(i);
      point->position.Reset(px, py, pz);

      // Get normal
      RNScalar nx = nx_image.GridValue(i);
      RNScalar ny = ny_image.GridValue(i);
      RNScalar nz = nz_image.GridValue(i);
      point->normal.Reset(nx, ny, nz);

      // Get radius
      RNScalar radius = radius_image.GridValue(i);
      point->radius = radius;

      // Get color
      point->color = color_image.PixelRGB(ix, iy);
    
      // Get flags
      point->boundary = (unsigned int) (boundary_image.GridValue(i) + 0.5);

      // Set grid index
      point->grid_index = i;

      // Insert point
      points.Insert(point);
    }
  }

  // Create kdtree of points
  RGBDPoint tmp; int position_offset = (unsigned char *) &(tmp.position) - (unsigned char *) &tmp;
  kdtree = new R3Kdtree<RGBDPoint *>(points, position_offset);
  if (!kdtree) {
    fprintf(stderr, "Unable to create kdtree\n");
    return 0;
  }
  
  // Create arrays of neighbor points
  for (int i = 0; i < points.NEntries(); i++) {
    RGBDPoint *point = points.Kth(i);
    int ix, iy, neighbor_index;
    depth_image.IndexToIndices(point->grid_index, ix, iy);
    for (int s = -1; s <= 1; s++) {
      if ((ix+s < 0) || (ix+s >= depth_image.XResolution())) continue;
      for (int t = -1; t <= 1; t++) {
        if ((s == 0) && (t == 0)) continue;
        if ((iy+t < 0) || (iy+t >= depth_image.YResolution())) continue;
        depth_image.IndicesToIndex(ix+s, iy+t, neighbor_index);
        RGBDPoint *neighbor = &point_buffer[neighbor_index];
        if ((point->boundary & RGBD_SHADOW_BOUNDARY) && (neighbor->boundary & RGBD_SILHOUETTE_BOUNDARY)) continue;
        if ((point->boundary & RGBD_SILHOUETTE_BOUNDARY) && (neighbor->boundary & RGBD_SHADOW_BOUNDARY)) continue;
        point->neighbors.Insert(neighbor);
      }
    }
  }

  // Return success
  return 1;
}
               


int RGBDSegmentation::
CreateSingletonSegments(int primitive_type)
{
  // Create segment for every point
  for (int i = 0; i < points.NEntries(); i++) {
    RGBDPoint *point = points.Kth(i);

    // Create primitive
    RGBDPrimitive primitive(primitive_type);
    primitive.Update(point);
    
    // Create segment
    RGBDSegment *segment = new RGBDSegment(this, point, primitive);

    // Insert point
    segment->InsertPoint(point, 1.0);

    // Insert segment
    segment->segmentation_index = segments.NEntries();    
    segments.Insert(segment);
  }

  // Return success
  return 1;
}



int RGBDSegmentation::
CreateRansacSegments(int primitive_type)
{
  // Check number of ransac iterations
  if (max_ransac_iterations == 0) return 1;
  
  // Determine how many seed points to skip each iteration
  int skip = 1;
  if ((max_segments > 0) && (points.NEntries()/(4*max_segments) > skip))
    skip = points.NEntries()/(4*max_segments);
  if ((min_segment_points > 0) & ((min_segment_points/4) > skip))
    skip = min_segment_points/4;
  if ((min_segments > 0) && (skip > points.NEntries()/min_segments))
    skip = points.NEntries()/min_segments;

  // Search seed points
  int seed_index = 0;
  while (seed_index < points.NEntries()) {
    // Find next seed point
    RGBDPoint *seed_point = NULL;
    while ((seed_index < points.NEntries()) && !seed_point) {
      RGBDPoint *point = points.Kth(seed_index);
      if (!point->segment || (point->segment_affinity < 0.1)) seed_point = point;
      seed_index += skip; 
    }

    // Check seed point
    if (!seed_point) break;

    // Create segment
    RGBDPrimitive primitive(primitive_type);
    primitive.Update(seed_point);
    RGBDSegment *segment = new RGBDSegment(this, seed_point, primitive);
    if (!segment->UpdatePoints(kdtree)) { delete segment; continue; }

    // Iteratively update everything
    RNBoolean error = FALSE;
    for (int iter = 0; iter < max_ransac_iterations; iter++) {
      if (!segment->UpdatePrimitive()) { error = TRUE; break; }
      if (!segment->UpdatePoints(kdtree)) { error = TRUE; break; }
    }

    // Check for error
    if (error) { 
      delete segment; 
      continue; 
    }

    // Insert segment
    segment->segmentation_index = segments.NEntries();    
    segments.Insert(segment);
  } 

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Segmenting manipulation functions
////////////////////////////////////////////////////////////////////////

int RGBDSegmentation::
RefineSegments(void)
{
  // Iteratively update everything
  int max_iterations = 1;
  for (int iter = 0; iter < max_iterations; iter++) {
    // Copy list of segments
    RNArray<RGBDSegment *> tmp = segments;
    tmp.Sort(RGBDCompareSegments);

    // Will rebuild list of segments
    segments.Empty();

    // Refine each segment
    RNBoolean converged = TRUE;
    for (int i = 0; i < tmp.NEntries(); i++) {
      RGBDSegment *segment = tmp.Kth(i);
      int prev_npoints = segment->points.NEntries();

      // Refine segment
      RNBoolean error = FALSE;
      if (!error && !segment->UpdatePrimitive()) error = TRUE;
      if (!error && !segment->UpdatePoints(kdtree)) error = TRUE; 

      // Insert segment
      segment->segmentation_index = segments.NEntries();    
      if (!error) segments.Insert(segment);
      else delete segment;

      // Check for convergence
      if (error || (prev_npoints != segment->points.NEntries())) {
        converged = FALSE;
      }
    }

    // Check if converged
    if (converged) break;
  }

  // Return success
  return 1;
}



int RGBDSegmentation::
DeleteSegments(void)
{
  // Sort segments
  segments.Sort(RGBDCompareSegments);

  // Separate viable from nonviable ones
  RNArray<RGBDSegment *> viable_segments;
  RNArray<RGBDSegment *> nonviable_segments;
  for (int i = 0; i < segments.NEntries(); i++) {
    RGBDSegment *segment = segments.Kth(i);

    // Check min_segments
    if ((min_segments <= 0) || (i >= min_segments)) {
      // Check segment points
      if (min_segment_points > 0) {
        if (segment->points.NEntries() < min_segment_points) {
          nonviable_segments.Insert(segment);
          continue;
        }
      }

      // Check segment coverage
      if (min_segment_coverage > 0) {
        if (segment->Coverage() < min_segment_coverage) {
          nonviable_segments.Insert(segment);
          continue;
        }
      }

      // Check max_segments
      if (max_segments > 0) {
        if (viable_segments.NEntries() > max_segments) {
          nonviable_segments.Insert(segment);
          continue;
        }
      }
    }
    
    // Segment is viable
    segment->segmentation_index = viable_segments.NEntries();    
    viable_segments.Insert(segment);
  }

  // Delete nonviable segments
  for (int i = 0; i < nonviable_segments.NEntries(); i++) {
    RGBDSegment *segment = nonviable_segments.Kth(i);
    delete segment;
  }

  // Replace segments with viable ones
  segments = viable_segments;

  // Return success
  return 1;
}



int RGBDSegmentation::
MergeSegments(void)
{
  // Initialize statistics
  int merge_count = 0;
  int push_count = 0;

  //////////

  // Create pairs between segments with nearby points
  RNArray<RGBDSegmentPair *> pairs;
  for (int i = 0; i < segments.NEntries(); i++) {
    RGBDSegment *segment0 = segments.Kth(i);
    if (segment0->points.IsEmpty()) continue;
    RGBDPoint *point0 = segment0->points.Head();

    // Create pairs
    for (int j = 0; j < point0->neighbors.NEntries(); j++) {
      RGBDPoint *point1 = point0->neighbors.Kth(j);
      if (point0 == point1) continue;
      RGBDSegment *segment1 = point1->segment;
      if (!segment1) continue;
      if (segment0 == segment1) continue;

      // Check if within max neighbor distance
      if (max_neighbor_distance_factor > 0) {
        RNScalar radius = (point0->radius > point1->radius) ? point1->radius : point0->radius;
        RNScalar max_neighbor_distance = (radius > 0) ? max_neighbor_distance_factor * radius : 0.25;
        RNScalar dd = R3SquaredDistance(point0->position, point1->position);
        if (dd > max_neighbor_distance * max_neighbor_distance) continue;
      }

      // Check if already have pair
      if (FindPair(segment0, segment1)) continue;

      // Compute affinity
      RNScalar affinity = segment0->Affinity(segment1);
      if (affinity < min_pair_affinity) continue;

      // Create pair
      RGBDSegmentPair *pair = new RGBDSegmentPair(segment0, segment1, affinity);
      if (!pair) continue;

      // Insert pair
      pairs.Insert(pair);
    }
  }

  // Check if there are any pairs
  if (pairs.IsEmpty()) return 1;

  //////////

  // Initialize heap
  RGBDSegmentPair tmp;
  RNHeap<RGBDSegmentPair *> heap(&tmp, &tmp.affinity, &tmp.heapentry, FALSE);
  for (int i = 0; i < pairs.NEntries(); i++) {
    RGBDSegmentPair *pair = pairs.Kth(i);
    heap.Push(pair);
  }

  // Merge segments hierarchically
  while (!heap.IsEmpty()) {
    // Get pair
    RGBDSegmentPair *pair = heap.Pop();

    // Check if we are done
    if (pair->affinity < min_pair_affinity) break;

    // Get segments
    RGBDSegment *segment0 = pair->segments[0];
    RGBDSegment *segment1 = pair->segments[1];

    // Check if either segment has already been merged
    if (segment0->parent || segment1->parent) {
      // Find ancestors
      RGBDSegment *ancestor0 = segment0;
      RGBDSegment *ancestor1 = segment1;
      while (ancestor0->parent) ancestor0 = ancestor0->parent;
      while (ancestor1->parent) ancestor1 = ancestor1->parent;
      if (ancestor0 != ancestor1) {
        if (!FindPair(ancestor0, ancestor1)) {
          RNScalar affinity = ancestor0->Affinity(ancestor1);
          if (affinity > min_pair_affinity) {
            // Create a pair between the ancestors
            RGBDSegmentPair *pair = new RGBDSegmentPair(ancestor0, ancestor1, affinity);
            heap.Push(pair);
            push_count++;
          }
        }
      }
    }
    else {
      if (0 && print_progress) {
        static unsigned long count = 0;
        if ((count++ % 1000) == 0) {
          printf("        %15.12f : %9d %9d : %15d %15d %15d\n", pair->affinity, 
                 segment0->points.NEntries(), segment1->points.NEntries(), 
                 heap.NEntries(), merge_count, push_count);
        }
      }

#if 0
      // Create merged segment
      RGBDSegment *segment = new RGBDSegment(this, segment0, segment1);
      segments.Insert(segment);
      merge_count++;
#else
      // Merge smaller segment into bigger one
      RGBDSegment *parent = (segment0->points.NEntries() > segment1->points.NEntries()) ? segment0 : segment1;
      RGBDSegment *child = (segment0->points.NEntries() > segment1->points.NEntries()) ? segment1 : segment0;
      parent->InsertChild(child);
      merge_count++;
#endif
    }

    // Delete pair
    delete pair;
  }

  // Remove merged segments
  RNArray<RGBDSegment *> merged_segments;
  RNArray<RGBDSegment *> all_segments = segments;
  segments.Empty();
  for (int i = 0; i < all_segments.NEntries(); i++) {
    RGBDSegment *segment = all_segments.Kth(i);
    segment->segmentation_index = segments.NEntries();    
    if (!segment->parent) { segments.Insert(segment); continue; }
    segment->parent->RemoveChild(segment);
    merged_segments.Insert(segment);
  }

  // Delete merged segments
  for (int i = 0; i < merged_segments.NEntries(); i++) {
    RGBDSegment *segment = merged_segments.Kth(i);
    delete segment;
  }

  // Return success
  return 1;
}



int RGBDSegmentation::
SplitSegments(void)
{
#if 0
  // Check min segments spacing
  if (min_segment_spacing <= 0) return 1;

  // Split connected components
  RNArray<RGBDSegment *> tmp = segments;
  segments.Empty();
  for (int i = 0; i < tmp.NEntries(); i++) {
    RGBDSegment *segment = tmp.Kth(i);

    // Check segment
    if (segment->points.NEntries() < min_segment_points) continue;

    // Rasterize points into planar grid
    R3PlanarGrid grid(segment->primitive.plane, segment->primitive.bbox, min_segment_spacing);
    for (int j = 0; j < segment->points.NEntries(); j++) {
      RGBDPoint *point = segment->points.Kth(j);
      R3Point position = point->position;
      grid.RasterizeWorldPoint(position.X(), position.Y(), position.Z(), 1.0);
    }

    // Compute connected components
    int max_components = grid.NEntries();
    int *components = new int [ max_components ];
    int ncomponents = grid.ConnectedComponents(RN_EPSILON, max_components, NULL, NULL, components);

    // Check connected components
    if (ncomponents == 1) {
      // One connected component - simply insert segment
      segment->segmentation_index = segments.NEntries();    
      segments.Insert(segment);
    }
    else {
      // Create segment for each connnected component
      for (int j = 0; j < ncomponents; j++) {
        // Make array of points in component
        RNArray<RGBDPoint *> component_points;
        for (int k = 0; k < segment->points.NEntries(); k++) {
          RGBDPoint *point = segment->points.Kth(k);
          R3Point world_position = point->position;
          R2Point grid_position = grid.GridPosition(world_position);
          int ix = (int) (grid_position.X() + 0.5);
          int iy = (int) (grid_position.Y() + 0.5);
          int index; grid.Grid().IndicesToIndex(ix, iy, index);
          if (components[index] != j) continue;
          component_points.Insert(point);
        }

        // Check number of points
        if (component_points.NEntries() > min_segment_points) {

          // Find centroid
          R3Point centroid = R3zero_point;
          for (int k = 0; k < component_points.NEntries(); k++) {
            RGBDPoint *point = component_points.Kth(k);
            R3Point world_position = point->position;
            centroid += world_position;
          }
          centroid /= component_points.NEntries();

          // Find seed point
          RGBDPoint *seed_point = NULL;
          RNLength min_dd = FLT_MAX;
          for (int k = 0; k < component_points.NEntries(); k++) {
            RGBDPoint *point = component_points.Kth(k);
            R3Point world_position = point->position;
            RNLength dd = R3SquaredDistance(centroid, world_position);
            if (dd < min_dd) { seed_point = point; min_dd = dd; }
          }

          // Check seed point
          if (seed_point) {
            // Create segment
            RGBDSegment *c = new RGBDSegment(this, seed_point, RGBD_PLANE_PRIMITIVE_TYPE);
            c->possible_affinity = segment->possible_affinity;

            // Insert points into segment
            for (int k = 0; k < component_points.NEntries(); k++) {
              RGBDPoint *point = component_points.Kth(k);
              c->InsertPoint(point);
            }

            // Update primitive
            c->UpdatePrimitive();

            // Update planar grid
            // c->UpdatePlanarGrid();

            // Insert segment
            c->segmentation_index = segments.NEntries();    
            segments.Insert(c);
          }
        }
      }

      // Delete the original segment
      delete segment;
    }

    // Delete components
    delete [] components;
  }
#endif

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Top-level segmentation functions
////////////////////////////////////////////////////////////////////////

int RGBDSegmentation::
CreateSegments(int primitive_type)
{
  // Check points
  if (points.NEntries() == 0) return 0;
  
  // Print debug message
  RNTime step_time;
  step_time.Read();
  if (print_progress) {
    printf("      SA %.3f\n", step_time.Elapsed());
    step_time.Read();
  }

  // Create segments
  if (initialize_hierarchically) {
    if (!CreateSingletonSegments(primitive_type)) return 0;
    if (!MergeSegments()) return 0;
  }
  else {
    if (!CreateRansacSegments(primitive_type)) return 0;
  }

  // Check segments
  if (segments.IsEmpty()) return 0;

  // Print debug message
  if (print_progress) {
    printf("      SB %.3f %d %d %g\n", step_time.Elapsed(), segments.NEntries(), NUnsegmentedPoints(), Affinity());
    step_time.Read();
  }

  // Iteratively update segments
  for (int i = 0; i < max_refinement_iterations; i++) {
    // Refine segments
    if (!RefineSegments()) return 0;

    // Print debug message
    if (print_progress) {
      printf("      SC %d : %.3f %d %d %g\n", i, step_time.Elapsed(), segments.NEntries(), NUnsegmentedPoints(), Affinity());
      step_time.Read();
    }

    // Create segments
    if (!CreateRansacSegments(primitive_type)) return 0;

    // Print debug message
    if (print_progress) {
      printf("      SD %d : %.3f %d %d %g\n", i, step_time.Elapsed(), segments.NEntries(), NUnsegmentedPoints(), Affinity());
      step_time.Read();
    }

    // Merge segments
    if (!MergeSegments()) return 0;

    // Print debug message
    if (print_progress) {
      printf("      SE %d : %.3f %d %d %g\n", i, step_time.Elapsed(), segments.NEntries(), NUnsegmentedPoints(), Affinity());
      step_time.Read();
    }

    // Delete segments
    if (!DeleteSegments()) return 0;

    // Print debug message
    if (print_progress) {
      printf("      SF %d : %.3f %d %d %g\n", i, step_time.Elapsed(), segments.NEntries(), NUnsegmentedPoints(), Affinity());
      step_time.Read();
    }
  }

  // Split segments
  // if (!SplitSegments()) return 0;

  // Print debug message
  if (print_progress) {
    printf("      SG %.3f %d %d %g\n", step_time.Elapsed(), segments.NEntries(), NUnsegmentedPoints(), Affinity());
    step_time.Read();
  }

  // Delete segments
  if (!DeleteSegments()) return 0;

  // Print debug message
  if (print_progress) {
    printf("      SH %.3f %d %d %g\n", step_time.Elapsed(), segments.NEntries(), NUnsegmentedPoints(), Affinity());
    step_time.Read();
  }

  // Sort segments
  segments.Sort(RGBDCompareSegments);

  // Return success
  return 1;
}




int RGBDSegmentation::
ReadSegmentImage(const char *filename)
{
  // Read file
  R2Grid segmentation_grid;
  if (!segmentation_grid.ReadFile(filename)) return 0;

  // Process segmentation grid
  segmentation_grid.Substitute(0, R2_GRID_UNKNOWN_VALUE);
  segmentation_grid.Subtract(1.0);

  // Create segments
  for (int i = 0; i < segmentation_grid.NEntries(); i++) {
    RGBDPoint *point = &point_buffer[i];
    if (point->depth < 0) continue;
    RNScalar grid_value = segmentation_grid.GridValue(i);
    if (grid_value == R2_GRID_UNKNOWN_VALUE) continue;
    int segment_index = (int) (grid_value + 0.5);
    while (segments.NEntries() <= segment_index) {
      // Create segment
      RGBDSegment *segment = new RGBDSegment(this, NULL, RGBD_PLANE_PRIMITIVE_TYPE);
      segment->segmentation_index = segments.NEntries();    
      segments.Insert(segment);
    }
    RGBDSegment *segment = segments.Kth(segment_index);
    if (!segment->seed_point) segment->seed_point = point;
    segment->InsertPoint(point);
  }

  // Update segment primitives
  for (int i = 0; i < segments.NEntries(); i++) {
    RGBDSegment *segment = segments.Kth(i);
    segment->UpdatePrimitive();
  }
  
  // Sort segments
  segments.Sort(RGBDCompareSegments);

  // Return success
  return 1;
}



int RGBDSegmentation::
WriteSegmentImage(int xres, int yres, const char *filename) const
{
  // Fill image
  R2Grid image(xres, yres);
  for (int i = 0; i < segments.NEntries(); i++) {
    RGBDSegment *segment = segments.Kth(i);
    for (int j = 0; j < segment->points.NEntries(); j++) {
      RGBDPoint *point = segment->points.Kth(j);
      if (point->grid_index < 0) continue;
      image.SetGridValue(point->grid_index, i+1);
    }
  }

  // Write image
  if (!image.WriteFile(filename)) return 0;
  
  // Return success
  return 1;   
}



void RGBDSegmentation::
SetDefaultParameters(void)
{
  min_segment_points = 10;
  min_segments = 0;
  max_segments = 0;
  min_segment_coverage = 0;
  max_segment_diameter = 16;
  max_segment_primitive_distance = 0.1;
  max_segment_normal_angle = RN_PI / 4.0;
  max_neighbor_distance_factor = 16;
  max_pair_centroid_distance = 16;
  max_pair_primitive_distance = 0.1;
  max_pair_normal_angle = RN_PI / 4.0;
  min_pair_affinity = 1.0E-6;
  initialize_hierarchically = TRUE;
  max_refinement_iterations = 0;
  max_ransac_iterations = 0;
  print_progress = FALSE;
}
