// Source file for R3Kdtree class

#ifndef __R3KDTREE__C__
#define __R3KDTREE__C__




////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "R3Shapes/R3Shapes.h"





////////////////////////////////////////////////////////////////////////
// Constant definitions
////////////////////////////////////////////////////////////////////////

static const int R3kdtree_max_points_per_node = 32;





////////////////////////////////////////////////////////////////////////
// Node class definition
////////////////////////////////////////////////////////////////////////

// Node declaration

template <class PtrType>
class R3KdtreeNode {
public:
  R3KdtreeNode(R3KdtreeNode<PtrType> *parent);

public:
  class R3KdtreeNode<PtrType> *parent;
  class R3KdtreeNode<PtrType> *children[2];
  RNScalar split_coordinate;
  RNDimension split_dimension;
  PtrType points[R3kdtree_max_points_per_node];
  int npoints;
};



// Node constructor

template <class PtrType>
R3KdtreeNode<PtrType>::
R3KdtreeNode(R3KdtreeNode<PtrType> *parent)
  : parent(parent),
    npoints(0)
{
  // Initialize everything
  children[0] = NULL;
  children[1] = NULL;
  split_coordinate = 0;
  split_dimension = 0;
}




////////////////////////////////////////////////////////////////////////
// Public tree-level functions
////////////////////////////////////////////////////////////////////////

template <class PtrType>
R3Kdtree<PtrType>::
R3Kdtree(const R3Box& bbox, int position_offset)
  : bbox(bbox),
    position_offset(position_offset),
    position_callback(NULL),
    position_callback_data(NULL),
    npoints(0),
    nnodes(1)
{
  // Create root node
  root = new R3KdtreeNode<PtrType>(NULL);
  assert(root);
}



template <class PtrType>
R3Kdtree<PtrType>::
R3Kdtree(const R3Box& bbox, R3Point (*position_callback)(PtrType, void *), void *position_callback_data)
  : bbox(bbox),
    position_offset(-1),
    position_callback(position_callback),
    position_callback_data(position_callback_data),
    npoints(0),
    nnodes(1)
{
  // Create root node
  root = new R3KdtreeNode<PtrType>(NULL);
  assert(root);
}



template <class PtrType>
R3Kdtree<PtrType>::
R3Kdtree(const RNArray<PtrType>& points, int position_offset)
  : position_offset(position_offset),
    position_callback(NULL),
    position_callback_data(NULL),
    npoints(points.NEntries()),
    nnodes(1)
{
  // Create root node
  root = new R3KdtreeNode<PtrType>(NULL);
  assert(root);

  // Determine bounding box 
  bbox = R3null_box;
  for (int i = 0; i < points.NEntries(); i++) {
    bbox.Union(Position(points[i]));
  }

  // Allocate copy of points array (so that it can be sorted)
  PtrType *copy = new PtrType [ points.NEntries() ];
  assert(copy);

  // Copy points
  for (int i = 0; i < points.NEntries(); i++) 
    copy[i] = points[i];

  // Insert points into root
  InsertPoints(root, bbox, copy, points.NEntries());

  // Delete copy of points array
  delete [] copy;
}



template <class PtrType>
R3Kdtree<PtrType>::
R3Kdtree(const RNArray<PtrType>& points, R3Point (*position_callback)(PtrType, void *), void *position_callback_data)
  : position_offset(-1),
    position_callback(position_callback),
    position_callback_data(position_callback_data),
    npoints(points.NEntries()),
    nnodes(1)
{
  // Create root node
  root = new R3KdtreeNode<PtrType>(NULL);
  assert(root);

  // Determine bounding box 
  bbox = R3null_box;
  for (int i = 0; i < points.NEntries(); i++) {
    bbox.Union(Position(points[i]));
  }

  // Allocate copy of points array (so that it can be sorted)
  PtrType *copy = new PtrType [ points.NEntries() ];
  assert(copy);

  // Copy points
  for (int i = 0; i < points.NEntries(); i++) 
    copy[i] = points[i];

  // Insert points into root
  InsertPoints(root, bbox, copy, points.NEntries());

  // Delete copy of points array
  delete [] copy;
}



template <class PtrType>
R3Kdtree<PtrType>::
R3Kdtree(const R3Kdtree<PtrType>& kdtree)
  : bbox(kdtree.bbox),
    position_offset(kdtree.position_offset),
    position_callback(kdtree.position_callback),
    position_callback_data(kdtree.position_callback_data),
    npoints(kdtree.npoints),
    nnodes(1)
{
  // Create root node
  root = new R3KdtreeNode<PtrType>(NULL);
  assert(root);

  // Copy nodes
  RNArray<R3KdtreeNode<PtrType> *> this_stack;
  RNArray<R3KdtreeNode<PtrType> *> that_stack;
  this_stack.Insert(root);
  that_stack.Insert(kdtree.root);
  while (!this_stack.IsEmpty()) {
    assert(!that_stack.IsEmpty());
    
    // Get nodes from stack
    R3KdtreeNode<PtrType> *this_node = this_stack.Tail(); this_stack.RemoveTail();
    R3KdtreeNode<PtrType> *that_node = that_stack.Tail(); that_stack.RemoveTail();

    // Copy properties
    this_node->split_coordinate = that_node->split_coordinate;
    this_node->split_dimension = that_node->split_dimension;
    this_node->npoints = that_node->npoints;

    // Copy contents
    if (that_node->children[0]) {
      // Create children
      this_node->children[0] = new R3KdtreeNode<PtrType>(this_node);
      this_node->children[1] = new R3KdtreeNode<PtrType>(this_node);
      nnodes += 2;

      // Push children onto stack
      assert(that_node->children[1]);
      this_stack.Insert(this_node->children[0]);
      that_stack.Insert(that_node->children[0]);
      this_stack.Insert(this_node->children[1]);
      that_stack.Insert(that_node->children[1]);
    }
    else {
      // Copy points
      assert(!that_node->children[1]);
      for (int i = 0; i < that_node->npoints; i++) {
        this_node->points[i] = that_node->points[i];
      }
    }
  }
}



template <class PtrType>
R3Kdtree<PtrType>::
~R3Kdtree(void)
{
  // Check root
  if (!root) return;

  // Traverse tree deleting nodes
  RNArray<R3KdtreeNode<PtrType> *> stack;
  stack.InsertTail(root);
  while (!stack.IsEmpty()) {
    R3KdtreeNode<PtrType> *node = stack.Tail();
    stack.RemoveTail();
    if (node->children[0]) stack.Insert(node->children[0]);
    if (node->children[1]) stack.Insert(node->children[1]);
    delete node;
  }
}



template <class PtrType>
const R3Box& R3Kdtree<PtrType>::
BBox(void) const
{
  // Return bounding box of the whole KD tree
  return bbox;
}



template <class PtrType>
int R3Kdtree<PtrType>::
NPoints(void) const
{
  // Return number of points
  return npoints;
}



template <class PtrType>
int R3Kdtree<PtrType>::
NNodes(void) const
{
  // Return number of nodes
  return nnodes;
}



#if 0

template <class PtrType>
int R3Kdtree<PtrType>::
NEntries(void) const
{
  // Check root
  if (!root) return 0;

  // Traverse tree to count number of points
  int npoints = 0;
  RNArray<R3KdtreeNode<PtrType> *> stack;
  stack.InsertTail(root);
  while (!stack.IsEmpty()) {
    R3KdtreeNode<PtrType> *node = stack.Tail();
    stack.RemoveTail();
    npoints += node->npoints;
    if (node->children[0]) stack.Insert(node->children[0]);
    if (node->children[1]) stack.Insert(node->children[1]);
  }

  // Return total number of points
  return npoints;
}

#endif



////////////////////////////////////////////////////////////////////////
// Finding the closest one point to a query point
////////////////////////////////////////////////////////////////////////

template <class PtrType>
void R3Kdtree<PtrType>::
FindClosest(R3KdtreeNode<PtrType> *node, const R3Box& node_box, 
  PtrType query_point, const R3Point& query_position, 
  RNScalar min_distance_squared, RNScalar max_distance_squared,
  int (*IsCompatible)(PtrType, PtrType, void *), void *compatible_data, 
  PtrType& closest_point, RNScalar& closest_distance_squared) const
{
  // Check if node is interior
  if (node->children[0]) {
    assert(node->children[1]);

    // Find and check axial distances from point to node box
    RNLength dx, dy, dz;
    if (RNIsGreater(query_position.X(), node_box.XMax())) dx = query_position.X() - node_box.XMax();
    else if (RNIsLess(query_position.X(), node_box.XMin())) dx = node_box.XMin()- query_position.X();
    else dx = 0.0;
    RNLength dx_squared = dx * dx;
    if (dx_squared >= closest_distance_squared) return;
    if (RNIsGreater(query_position.Y(), node_box.YMax())) dy = query_position.Y() - node_box.YMax();
    else if (RNIsLess(query_position.Y(), node_box.YMin())) dy = node_box.YMin()- query_position.Y();
    else dy = 0.0;
    RNLength dy_squared = dy * dy;
    if (dy_squared >= closest_distance_squared) return;
    if (RNIsGreater(query_position.Z(), node_box.ZMax())) dz = query_position.Z() - node_box.ZMax();
    else if (RNIsLess(query_position.Z(), node_box.ZMin())) dz = node_box.ZMin()- query_position.Z();
    else dz = 0.0;
    RNLength dz_squared = dz * dz;
    if (dz_squared >= closest_distance_squared) return;
    
    // Find and check actual distance from point to node box
    RNLength distance_squared = 0;
    if ((dy == 0.0) && (dz == 0.0)) distance_squared = dx_squared;
    else if ((dx == 0.0) && (dz == 0.0)) distance_squared = dy_squared;
    else if ((dx == 0.0) && (dy == 0.0)) distance_squared = dz_squared;
    else distance_squared = dx_squared + dy_squared + dz_squared;
    if (distance_squared >= closest_distance_squared) return;

    // Compute distance from point to split plane
    RNLength side = query_position[node->split_dimension] - node->split_coordinate;

    // Search children nodes
    if (side <= 0) {
      // Search negative side first
      R3Box child_box(node_box);
      child_box[RN_HI][node->split_dimension] = node->split_coordinate;
      FindClosest(node->children[0], child_box, 
        query_point, query_position, 
        min_distance_squared, max_distance_squared, 
        IsCompatible, compatible_data,
        closest_point, closest_distance_squared);
      if (side*side < closest_distance_squared) {
        R3Box child_box(node_box);
        child_box[RN_LO][node->split_dimension] = node->split_coordinate;
        FindClosest(node->children[1], child_box, 
          query_point, query_position, 
          min_distance_squared, max_distance_squared, 
          IsCompatible, compatible_data,
          closest_point, closest_distance_squared);
      }
    }
    else {
      // Search positive side first
      R3Box child_box(node_box);
      child_box[RN_LO][node->split_dimension] = node->split_coordinate;
      FindClosest(node->children[1], child_box, 
        query_point, query_position, 
        min_distance_squared, max_distance_squared, 
        IsCompatible, compatible_data,
        closest_point, closest_distance_squared);
      if (side*side < closest_distance_squared) {
        R3Box child_box(node_box);
        child_box[RN_HI][node->split_dimension] = node->split_coordinate;
        FindClosest(node->children[0], child_box, 
          query_point, query_position, 
          min_distance_squared, max_distance_squared, 
          IsCompatible, compatible_data,
          closest_point, closest_distance_squared);
      }
    }
  }
  else {
    for (int i = 0; i < node->npoints; i++) {
      PtrType point = node->points[i];
      RNLength distance_squared = R3SquaredDistance(query_position, Position(point));
      if ((distance_squared >= min_distance_squared) && 
         (distance_squared <= closest_distance_squared)) {
        if (!IsCompatible || !query_point || IsCompatible(query_point, point, compatible_data)) {
          closest_distance_squared = distance_squared;
          closest_point = point;
        }
      }
    }
  }
}



template <class PtrType>
PtrType R3Kdtree<PtrType>::
FindClosest(PtrType query_point, 
  RNScalar min_distance, RNScalar max_distance, 
  int (*IsCompatible)(PtrType, PtrType, void *), void *compatible_data, 
  RNScalar *closest_distance) const
{
  // Check root
  if (!root) return NULL;

  // Use squared distances for efficiency
  if (max_distance < 0) return NULL;
  if (min_distance < 0) min_distance = 0;
  RNLength min_distance_squared = min_distance * min_distance;
  RNLength max_distance_squared = max_distance * max_distance;

  // Initialize nearest point 
  PtrType nearest_point = NULL;
  RNLength nearest_distance_squared = max_distance_squared;

  // Search nodes recursively
  FindClosest(root, bbox, 
    query_point, Position(query_point),
    min_distance_squared, max_distance_squared, 
    IsCompatible, compatible_data, 
    nearest_point, nearest_distance_squared);

  // Return closest distance
  if (closest_distance) *closest_distance = sqrt(nearest_distance_squared);

  // Return closest point
  return nearest_point;
}



template <class PtrType>
PtrType R3Kdtree<PtrType>::
FindClosest(PtrType query_point, 
  RNScalar min_distance, RNScalar max_distance, 
  RNScalar *closest_distance) const
{
  // Find the closest point
  return FindClosest(query_point, min_distance, max_distance, NULL, NULL, closest_distance);
}



template <class PtrType>
PtrType R3Kdtree<PtrType>::
FindClosest(const R3Point& query_position, 
  RNScalar min_distance, RNScalar max_distance, 
  RNScalar *closest_distance) const
{
  // Check root
  if (!root) return NULL;

  // Use squared distances for efficiency
  if (max_distance < 0) return NULL;
  if (min_distance < 0) min_distance = 0;
  RNLength min_distance_squared = min_distance * min_distance;
  RNLength max_distance_squared = max_distance * max_distance;

  // Initialize nearest point 
  PtrType nearest_point = NULL;
  RNLength nearest_distance_squared = max_distance_squared;

  // Search nodes recursively
  FindClosest(root, bbox, 
    NULL, query_position, 
    min_distance_squared, max_distance_squared, 
    NULL, NULL, 
    nearest_point, nearest_distance_squared);

  // Return closest distance
  if (closest_distance) *closest_distance = sqrt(nearest_distance_squared);

  // Return closest point
  return nearest_point;
}



////////////////////////////////////////////////////////////////////////
// Finding the closest K points to a query point
////////////////////////////////////////////////////////////////////////

template <class PtrType>
void R3Kdtree<PtrType>::
FindClosest(R3KdtreeNode<PtrType> *node, const R3Box& node_box, 
  PtrType query_point, const R3Point& query_position, 
  RNScalar min_distance_squared, RNScalar max_distance_squared, int max_points,
  int (*IsCompatible)(PtrType, PtrType, void *), void *compatible_data, 
  RNArray<PtrType>& points, RNLength *distances_squared) const
{
  // Update max distance squared
  if (points.NEntries() == max_points) {
    max_distance_squared = distances_squared[max_points-1];
  }

  // Check if node is interior
  if (node->children[0]) {
    assert(node->children[1]);

    // Find and check axial distances from point to node box
    RNLength dx, dy, dz;
    if (RNIsGreater(query_position.X(), node_box.XMax())) dx = query_position.X() - node_box.XMax();
    else if (RNIsLess(query_position.X(), node_box.XMin())) dx = node_box.XMin()- query_position.X();
    else dx = 0.0;
    RNLength dx_squared = dx * dx;
    if (dx_squared > max_distance_squared) return;
    if (RNIsGreater(query_position.Y(), node_box.YMax())) dy = query_position.Y() - node_box.YMax();
    else if (RNIsLess(query_position.Y(), node_box.YMin())) dy = node_box.YMin()- query_position.Y();
    else dy = 0.0;
    RNLength dy_squared = dy * dy;
    if (dy_squared > max_distance_squared) return;
    if (RNIsGreater(query_position.Z(), node_box.ZMax())) dz = query_position.Z() - node_box.ZMax();
    else if (RNIsLess(query_position.Z(), node_box.ZMin())) dz = node_box.ZMin()- query_position.Z();
    else dz = 0.0;
    RNLength dz_squared = dz * dz;
    if (dz_squared > max_distance_squared) return;
    
    // Find and check actual distance from point to node box
    RNLength distance_squared = 0;
    if ((dy == 0.0) && (dz == 0.0)) distance_squared = dx_squared;
    else if ((dx == 0.0) && (dz == 0.0)) distance_squared = dy_squared;
    else if ((dx == 0.0) && (dy == 0.0)) distance_squared = dz_squared;
    else distance_squared = dx_squared + dy_squared + dz_squared;
    if (distance_squared > max_distance_squared) return;

    // Compute distance from point to split plane
    RNLength side = query_position[node->split_dimension] - node->split_coordinate;

    // Search children nodes
    if ((side <= 0) || (side*side <= max_distance_squared)) {
      // Search negative side 
      R3Box child_box(node_box);
      child_box[RN_HI][node->split_dimension] = node->split_coordinate;
      FindClosest(node->children[0], child_box, query_point, query_position, 
        min_distance_squared, max_distance_squared, max_points, IsCompatible, compatible_data,
        points, distances_squared);
    }
    if ((side >= 0) || (side*side <= max_distance_squared)) {
      R3Box child_box(node_box);
      child_box[RN_LO][node->split_dimension] = node->split_coordinate;
      FindClosest(node->children[1], child_box, query_point, query_position, 
        min_distance_squared, max_distance_squared, max_points, IsCompatible, compatible_data,
        points, distances_squared);
    }
  }
  else {
    // Search points
    for (int i = 0; i < node->npoints; i++) {
      PtrType point = node->points[i];
      RNLength distance_squared = R3SquaredDistance(query_position, Position(point));
      if ((distance_squared >= min_distance_squared) && 
          (distance_squared <= max_distance_squared)) {

        // Check if point is compatible
        if (!IsCompatible || !query_point || IsCompatible(query_point, point, compatible_data)) {

          // Find slot for point (points are sorted by distance)
          int slot = 0;
          while (slot < points.NEntries()) {
            if (distance_squared < distances_squared[slot]) break;
            slot++;
          }
          
          // Insert point and distance into sorted arrays
          if (slot < max_points) {
            int first = points.NEntries();
            if (first >= max_points) first = max_points-1;
            for (int j = first; j > slot; j--) distances_squared[j] = distances_squared[j-1];
            distances_squared[slot] = distance_squared;
            points.InsertKth(point, slot);
            points.Truncate(max_points);
          }
        }
      }
    }
  }
}



template <class PtrType>
int R3Kdtree<PtrType>::
FindClosest(PtrType query_point, RNScalar min_distance, RNScalar max_distance, int max_points, 
  RNArray<PtrType>& points, RNLength *distances) const
{
  // Find closest within some distance
  return FindClosest(query_point, min_distance, max_distance, max_points, NULL, NULL, points, distances);
}




template <class PtrType>
int R3Kdtree<PtrType>::
FindClosest(PtrType query_point, 
  RNScalar min_distance, RNScalar max_distance, int max_points, 
  int (*IsCompatible)(PtrType, PtrType, void *), void *compatible_data, 
  RNArray<PtrType>& points, RNLength *distances) const
{
  // Check root
  if (!root) return 0;

  // Use squared distances for efficiency
  if (max_distance < 0) return 0;
  if (min_distance < 0) min_distance = 0;
  RNLength min_distance_squared = min_distance * min_distance;
  RNLength max_distance_squared = max_distance * max_distance;

  // Allocate temporary array of squared distances to max_points closest points
  RNLength *distances_squared = new RNLength [ max_points ];

  // Search nodes recursively
  FindClosest(root, bbox, 
    query_point, Position(query_point),
    min_distance_squared, max_distance_squared, max_points, 
    IsCompatible, compatible_data,
    points, distances_squared);

  // Update return distances
  if (distances) {
    for (int i = 0; i < points.NEntries(); i++) {
      distances[i] = sqrt(distances_squared[i]);
    }
  }

  // Delete temporary array of squared distances
  delete [] distances_squared;

  // Return number of points
  return points.NEntries();
}



template <class PtrType>
int R3Kdtree<PtrType>::
FindClosest(const R3Point& query_position, RNScalar min_distance, RNScalar max_distance, int max_points, 
  RNArray<PtrType>& points, RNLength *distances) const
{
  // Check root
  if (!root) return 0;

  // Use squared distances for efficiency
  if (max_distance < 0) return 0;
  if (min_distance < 0) min_distance = 0;
  RNLength min_distance_squared = min_distance * min_distance;
  RNLength max_distance_squared = max_distance * max_distance;

  // Allocate temporary array of squared distances to max_points closest points
  RNLength *distances_squared = new RNLength [ max_points ];

  // Search nodes recursively
  FindClosest(root, bbox, 
    NULL, query_position, 
    min_distance_squared, max_distance_squared, max_points, 
    NULL, NULL, 
    points, distances_squared);

  // Update return distances
  if (distances) {
    for (int i = 0; i < points.NEntries(); i++) {
      distances[i] = sqrt(distances_squared[i]);
    }
  }

  // Delete temporary array of squared distances
  delete [] distances_squared;

  // Return number of points
  return points.NEntries();
}



////////////////////////////////////////////////////////////////////////
// Finding all points within some distance to a query point
////////////////////////////////////////////////////////////////////////

template <class PtrType>
void R3Kdtree<PtrType>::
FindAll(R3KdtreeNode<PtrType> *node, const R3Box& node_box, 
  PtrType query_point, const R3Point& query_position, 
  RNScalar min_distance_squared, RNScalar max_distance_squared, 
  int (*IsCompatible)(PtrType, PtrType, void *), void *compatible_data, 
  RNArray<PtrType>& points) const
{
  // Check if node is interior
  if (node->children[0]) {
    assert(node->children[1]);

    // Find and check axial distances from point to node box
    RNLength dx, dy, dz;
    if (RNIsGreater(query_position.X(), node_box.XMax())) dx = query_position.X() - node_box.XMax();
    else if (RNIsLess(query_position.X(), node_box.XMin())) dx = node_box.XMin()- query_position.X();
    else dx = 0.0;
    RNLength dx_squared = dx * dx;
    if (dx_squared > max_distance_squared) return;
    if (RNIsGreater(query_position.Y(), node_box.YMax())) dy = query_position.Y() - node_box.YMax();
    else if (RNIsLess(query_position.Y(), node_box.YMin())) dy = node_box.YMin()- query_position.Y();
    else dy = 0.0;
    RNLength dy_squared = dy * dy;
    if (dy_squared > max_distance_squared) return;
    if (RNIsGreater(query_position.Z(), node_box.ZMax())) dz = query_position.Z() - node_box.ZMax();
    else if (RNIsLess(query_position.Z(), node_box.ZMin())) dz = node_box.ZMin()- query_position.Z();
    else dz = 0.0;
    RNLength dz_squared = dz * dz;
    if (dz_squared > max_distance_squared) return;
    
    // Find and check actual distance from point to node box
    RNLength distance_squared = 0;
    if ((dy == 0.0) && (dz == 0.0)) distance_squared = dx_squared;
    else if ((dx == 0.0) && (dz == 0.0)) distance_squared = dy_squared;
    else if ((dx == 0.0) && (dy == 0.0)) distance_squared = dz_squared;
    else distance_squared = dx_squared + dy_squared + dz_squared;
    if (distance_squared > max_distance_squared) return;

    // Compute distance from point to split plane
    RNLength side = query_position[node->split_dimension] - node->split_coordinate;

    // Search children nodes
    if ((side <= 0) || (side*side <= max_distance_squared)) {
      // Search negative side 
      R3Box child_box(node_box);
      child_box[RN_HI][node->split_dimension] = node->split_coordinate;
      FindAll(node->children[0], child_box, query_point, query_position, 
        min_distance_squared, max_distance_squared, 
        IsCompatible, compatible_data, points);
    }
    if ((side >= 0) || (side*side <= max_distance_squared)) {
      R3Box child_box(node_box);
      child_box[RN_LO][node->split_dimension] = node->split_coordinate;
      FindAll(node->children[1], child_box, query_point, query_position, 
        min_distance_squared, max_distance_squared, 
        IsCompatible, compatible_data, points);
    }
  }
  else {
    // Search points
    for (int i = 0; i < node->npoints; i++) {
      PtrType point = node->points[i];
      RNLength distance_squared = R3SquaredDistance(query_position, Position(point));
      if ((distance_squared >= min_distance_squared) && 
          (distance_squared <= max_distance_squared)) {
        if (!IsCompatible || !query_point || IsCompatible(query_point, point, compatible_data)) {
          points.Insert(point);
        }
      }
    }
  }
}



template <class PtrType>
int R3Kdtree<PtrType>::
FindAll(PtrType query_point,
  RNScalar min_distance, RNScalar max_distance, 
  int (*IsCompatible)(PtrType, PtrType, void *), void *compatible_data, 
  RNArray<PtrType>& points) const
{
  // Check root
  if (!root) return 0;

  // Use squared distances for efficiency
  if (max_distance < 0) return 0;
  if (min_distance < 0) min_distance = 0;
  RNLength min_distance_squared = min_distance * min_distance;
  RNLength max_distance_squared = max_distance * max_distance;

  // Search nodes recursively
  FindAll(root, bbox, 
    query_point, Position(query_point), 
    min_distance_squared, max_distance_squared, 
    IsCompatible, compatible_data, 
    points);

  // Return number of points
  return points.NEntries();
}



template <class PtrType>
int R3Kdtree<PtrType>::
FindAll(PtrType query_point, RNScalar min_distance, RNScalar max_distance, RNArray<PtrType>& points) const
{
  // Find all within some distance
  return FindAll(query_point, min_distance, max_distance, NULL, NULL, points);
}




template <class PtrType>
int R3Kdtree<PtrType>::
FindAll(const R3Point& query_position, RNScalar min_distance, RNScalar max_distance, RNArray<PtrType>& points) const
{
  // Check root
  if (!root) return 0;

  // Use squared distances for efficiency
  if (max_distance < 0) return 0;
  if (min_distance < 0) min_distance = 0;
  RNLength min_distance_squared = min_distance * min_distance;
  RNLength max_distance_squared = max_distance * max_distance;

  // Search nodes recursively
  FindAll(root, bbox, 
    NULL, query_position, 
    min_distance_squared, max_distance_squared, 
    NULL, NULL, 
    points);

  // Return number of points
  return points.NEntries();
}



////////////////////////////////////////////////////////////////////////
// Finding the any one point to a query point
////////////////////////////////////////////////////////////////////////

template <class PtrType>
PtrType R3Kdtree<PtrType>::
FindAny(R3KdtreeNode<PtrType> *node, const R3Box& node_box, 
  PtrType query_point, const R3Point& query_position, 
  RNScalar min_distance_squared, RNScalar max_distance_squared,
  int (*IsCompatible)(PtrType, PtrType, void *), void *compatible_data) const
{
  // Check if node is interior
  if (node->children[0]) {
    assert(node->children[1]);

    // Find and check axial distances from point to node box
    RNLength dx, dy, dz;
    if (RNIsGreater(query_position.X(), node_box.XMax())) dx = query_position.X() - node_box.XMax();
    else if (RNIsLess(query_position.X(), node_box.XMin())) dx = node_box.XMin()- query_position.X();
    else dx = 0.0;
    RNLength dx_squared = dx * dx;
    if (dx_squared >= max_distance_squared) return NULL;
    if (RNIsGreater(query_position.Y(), node_box.YMax())) dy = query_position.Y() - node_box.YMax();
    else if (RNIsLess(query_position.Y(), node_box.YMin())) dy = node_box.YMin()- query_position.Y();
    else dy = 0.0;
    RNLength dy_squared = dy * dy;
    if (dy_squared >= max_distance_squared) return NULL;
    if (RNIsGreater(query_position.Z(), node_box.ZMax())) dz = query_position.Z() - node_box.ZMax();
    else if (RNIsLess(query_position.Z(), node_box.ZMin())) dz = node_box.ZMin()- query_position.Z();
    else dz = 0.0;
    RNLength dz_squared = dz * dz;
    if (dz_squared >= max_distance_squared) return NULL;
    
    // Find and check actual distance from point to node box
    RNLength distance_squared = 0;
    if ((dy == 0.0) && (dz == 0.0)) distance_squared = dx_squared;
    else if ((dx == 0.0) && (dz == 0.0)) distance_squared = dy_squared;
    else if ((dx == 0.0) && (dy == 0.0)) distance_squared = dz_squared;
    else distance_squared = dx_squared + dy_squared + dz_squared;
    if (distance_squared >= max_distance_squared) return NULL;

    // Compute distance from point to split plane
    RNLength side = query_position[node->split_dimension] - node->split_coordinate;

    // Search children nodes
    if (side <= 0) {
      // Search negative side first
      R3Box child_box(node_box);
      child_box[RN_HI][node->split_dimension] = node->split_coordinate;
      PtrType any_point = FindAny(node->children[0], child_box, 
        query_point, query_position, 
        min_distance_squared, max_distance_squared, 
        IsCompatible, compatible_data);
      if (any_point) return any_point;
      if (side*side < max_distance_squared) {
        R3Box child_box(node_box);
        child_box[RN_LO][node->split_dimension] = node->split_coordinate;
        return FindAny(node->children[1], child_box, 
          query_point, query_position, 
          min_distance_squared, max_distance_squared, 
          IsCompatible, compatible_data);
      }
    }
    else {
      // Search positive side first
      R3Box child_box(node_box);
      child_box[RN_LO][node->split_dimension] = node->split_coordinate;
      PtrType any_point = FindAny(node->children[1], child_box, 
        query_point, query_position, 
        min_distance_squared, max_distance_squared, 
        IsCompatible, compatible_data);
      if (any_point) return any_point;
      if (side*side < max_distance_squared) {
        R3Box child_box(node_box);
        child_box[RN_HI][node->split_dimension] = node->split_coordinate;
        return FindAny(node->children[0], child_box, 
          query_point, query_position, 
          min_distance_squared, max_distance_squared, 
          IsCompatible, compatible_data);

      }
    }
  }
  else {
    for (int i = 0; i < node->npoints; i++) {
      PtrType point = node->points[i];
      RNLength distance_squared = R3SquaredDistance(query_position, Position(point));
      if ((distance_squared >= min_distance_squared) && 
         (distance_squared <= max_distance_squared)) {
        if (!IsCompatible || !query_point || IsCompatible(query_point, point, compatible_data)) {
          return point;
        }
      }
    }
  }

  // No point found
  return NULL;
}



template <class PtrType>
PtrType R3Kdtree<PtrType>::
FindAny(PtrType query_point, 
  RNScalar min_distance, RNScalar max_distance, 
  int (*IsCompatible)(PtrType, PtrType, void *), void *compatible_data) const
{
  // Check root
  if (!root) return NULL;

  // Use squared distances for efficiency
  if (max_distance < 0) return NULL;
  if (min_distance < 0) min_distance = 0;
  RNLength min_distance_squared = min_distance * min_distance;
  RNLength max_distance_squared = max_distance * max_distance;

  // Search nodes recursively
  return FindAny(root, bbox, 
    query_point, Position(query_point),
    min_distance_squared, max_distance_squared, 
    IsCompatible, compatible_data);
}



template <class PtrType>
PtrType R3Kdtree<PtrType>::
FindAny(PtrType query_point, 
  RNScalar min_distance, RNScalar max_distance) const
{
  // Find the any point
  return FindAny(query_point, min_distance, max_distance, NULL, NULL);
}



template <class PtrType>
PtrType R3Kdtree<PtrType>::
FindAny(const R3Point& query_position, 
  RNScalar min_distance, RNScalar max_distance) const
{
  // Check root
  if (!root) return NULL;

  // Use squared distances for efficiency
  if (max_distance < 0) return NULL;
  if (min_distance < 0) min_distance = 0;
  RNLength min_distance_squared = min_distance * min_distance;
  RNLength max_distance_squared = max_distance * max_distance;

  // Search nodes recursively
  return FindAny(root, bbox, 
    NULL, query_position, 
    min_distance_squared, max_distance_squared,
    NULL, NULL);
}



////////////////////////////////////////////////////////////////////////
// Finding the closest one point to a query shape
////////////////////////////////////////////////////////////////////////

template <class PtrType>
template <class Shape>
void R3Kdtree<PtrType>::
FindClosest(R3KdtreeNode<PtrType> *node, const R3Box& node_box, 
  const Shape& query_shape, 
  RNScalar min_distance, RNScalar max_distance,
  PtrType& closest_point, RNScalar& closest_distance) const
{
  // Check if node is interior
  if (node->children[0]) {
    assert(node->children[1]);

    // Check distance from shape to node box
    RNLength distance = R3Distance(query_shape, node_box);
    if (distance >= closest_distance) return;

    // Search negative side 
    R3Box child_box1(node_box);
    child_box1[RN_HI][node->split_dimension] = node->split_coordinate;
    FindClosest(node->children[0], child_box1, 
      query_shape, 
      min_distance, max_distance, 
      closest_point, closest_distance);

    // Search positive side 
    R3Box child_box2(node_box);
    child_box2[RN_LO][node->split_dimension] = node->split_coordinate;
    FindClosest(node->children[1], child_box2, 
      query_shape, 
      min_distance, max_distance, 
      closest_point, closest_distance);
  }
  else {
    for (int i = 0; i < node->npoints; i++) {
      PtrType point = node->points[i];
      RNLength distance = R3Distance(query_shape, Position(point));
      if ((distance >= min_distance) && 
         (distance <= closest_distance)) {
        closest_distance = distance;
        closest_point = point;
      }
    }
  }
}



template <class PtrType>
PtrType R3Kdtree<PtrType>::
FindClosest(const R3Line& query_line, 
  RNScalar min_distance, RNScalar max_distance, 
  RNScalar *closest_distance) const
{
  // Check root
  if (!root) return NULL;

  // Initialize nearest point 
  PtrType nearest_point = NULL;
  RNLength nearest_distance = max_distance;

  // Search nodes recursively
  FindClosest<R3Line>(root, bbox, 
    query_line, 
    min_distance, max_distance, 
    NULL, NULL, 
    nearest_point, nearest_distance);

  // Return closest distance
  if (closest_distance) *closest_distance = nearest_distance;

  // Return closest point
  return nearest_point;
}



template <class PtrType>
PtrType R3Kdtree<PtrType>::
FindClosest(const R3Plane& query_plane, 
  RNScalar min_distance, RNScalar max_distance, 
  RNScalar *closest_distance) const
{
  // Check root
  if (!root) return NULL;

  // Initialize nearest point 
  PtrType nearest_point = NULL;
  RNLength nearest_distance = max_distance;

  // Search nodes recursively
  FindClosest<R3Plane>(root, bbox, 
    query_plane, 
    min_distance, max_distance, 
    NULL, NULL, 
    nearest_point, nearest_distance);

  // Return closest distance
  if (closest_distance) *closest_distance = nearest_distance;

  // Return closest point
  return nearest_point;
}



template <class PtrType>
PtrType R3Kdtree<PtrType>::
FindClosest(const R3Shape& query_shape, 
  RNScalar min_distance, RNScalar max_distance, 
  RNScalar *closest_distance) const
{
  // Check root
  if (!root) return NULL;

  // Initialize nearest point 
  PtrType nearest_point = NULL;
  RNLength nearest_distance = max_distance;

  // Search nodes recursively
  FindClosest<R3Shape>(root, bbox, 
    query_shape, 
    min_distance, max_distance, 
    NULL, NULL, 
    nearest_point, nearest_distance);

  // Return closest distance
  if (closest_distance) *closest_distance = nearest_distance;

  // Return closest point
  return nearest_point;
}



////////////////////////////////////////////////////////////////////////
// Finding the closest K points to a query shape
////////////////////////////////////////////////////////////////////////

template <class PtrType>
template <class Shape>
void R3Kdtree<PtrType>::
FindClosest(R3KdtreeNode<PtrType> *node, const R3Box& node_box, 
  const Shape& query_shape, 
  RNScalar min_distance, RNScalar max_distance, int max_points,
  RNArray<PtrType>& points, RNLength *distances) const
{
  // Update max distance 
  if (points.NEntries() == max_points) {
    max_distance = distances[max_points-1];
  }

  // Check if node is interior
  if (node->children[0]) {
    assert(node->children[1]);

    // Check distance from shape to node box
    RNLength distance = R3Distance(query_shape, node_box);
    if (distance >= max_distance) return;

    // Search negative side 
    R3Box child_box1(node_box);
    child_box1[RN_HI][node->split_dimension] = node->split_coordinate;
    FindClosest(node->children[0], child_box1, 
      query_shape,
      min_distance, max_distance, max_points, 
      points, distances);

    // Search positive side
    R3Box child_box2(node_box);
    child_box2[RN_LO][node->split_dimension] = node->split_coordinate;
    FindClosest(node->children[1], child_box2, 
      query_shape, 
      min_distance, max_distance, max_points, 
      points, distances);
  }
  else {
    // Search points
    for (int i = 0; i < node->npoints; i++) {
      PtrType point = node->points[i];
      RNLength distance = R3Distance(query_shape, Position(point));
      if ((distance >= min_distance) && 
          (distance <= max_distance)) {

        // Find slot for point (points are sorted by distance)
        int slot = 0;
        while (slot < points.NEntries()) {
          if (distance < distances[slot]) break;
          slot++;
        }
        
        // Insert point and distance into sorted arrays
        if (slot < max_points) {
          int first = points.NEntries();
          if (first >= max_points) first = max_points-1;
          for (int j = first; j > slot; j--) distances[j] = distances[j-1];
          distances[slot] = distance;
          points.InsertKth(point, slot);
          points.Truncate(max_points);
        }
      }
    }
  }
}



template <class PtrType>
int R3Kdtree<PtrType>::
FindClosest(const R3Line& query_line, 
  RNScalar min_distance, RNScalar max_distance, int max_points, 
  RNArray<PtrType>& points, RNLength *distances) const
{
  // Check root
  if (!root) return 0;

  // Allocate temporary array of distances to max_points closest points
  RNLength *tmp_distances = (distances) ? distances : new RNLength [ max_points ];

  // Search nodes recursively
  FindClosest<R3Line>(root, bbox, 
    query_line, 
    min_distance, max_distance, max_points, 
    points, tmp_distances);

  // Delete temporary array of squared distances
  if (!distances) delete [] tmp_distances;

  // Return number of points
  return points.NEntries();
}



template <class PtrType>
int R3Kdtree<PtrType>::
FindClosest(const R3Plane& query_plane, 
  RNScalar min_distance, RNScalar max_distance, int max_points, 
  RNArray<PtrType>& points, RNLength *distances) const
{
  // Check root
  if (!root) return 0;

  // Allocate temporary array of distances to max_points closest points
  RNLength *tmp_distances = (distances) ? distances : new RNLength [ max_points ];

  // Search nodes recursively
  FindClosest<R3Plane>(root, bbox, 
    query_plane, 
    min_distance, max_distance, max_points, 
    points, tmp_distances);

  // Delete temporary array of squared distances
  if (!distances) delete [] tmp_distances;

  // Return number of points
  return points.NEntries();
}



template <class PtrType>
int R3Kdtree<PtrType>::
FindClosest(const R3Shape& query_shape, 
  RNScalar min_distance, RNScalar max_distance, int max_points, 
  RNArray<PtrType>& points, RNLength *distances) const
{
  // Check root
  if (!root) return 0;

  // Allocate temporary array of distances to max_points closest points
  RNLength *tmp_distances = (distances) ? distances : new RNLength [ max_points ];

  // Search nodes recursively
  FindClosest<R3Shape>(root, bbox, 
    query_shape, 
    min_distance, max_distance, max_points, 
    points, tmp_distances);

  // Delete temporary array of squared distances
  if (!distances) delete [] tmp_distances;

  // Return number of points
  return points.NEntries();
}



////////////////////////////////////////////////////////////////////////
// Finding all points within some distance to a query shape
////////////////////////////////////////////////////////////////////////

template <class PtrType>
template <class Shape>
void R3Kdtree<PtrType>::
FindAll(R3KdtreeNode<PtrType> *node, const R3Box& node_box, 
  const Shape& query_shape, 
  RNScalar min_distance, RNScalar max_distance, 
  RNArray<PtrType>& points) const
{
  // Check if node is interior
  if (node->children[0]) {
    assert(node->children[1]);

    // Check distance from shape to node box
    RNLength distance = R3Distance(query_shape, node_box);
    if (distance >= max_distance) return;

    // Search negative side 
    R3Box child_box1(node_box);
    child_box1[RN_HI][node->split_dimension] = node->split_coordinate;
    FindAll(node->children[0], child_box1, query_shape,
      min_distance, max_distance, 
      points);

    // Search positive side
    R3Box child_box2(node_box);
    child_box2[RN_LO][node->split_dimension] = node->split_coordinate;
    FindAll(node->children[1], child_box2, query_shape,
      min_distance, max_distance, 
      points);
  }
  else {
    // Search points
    for (int i = 0; i < node->npoints; i++) {
      PtrType point = node->points[i];
      RNLength distance = R3Distance(query_shape, Position(point));
      if ((distance >= min_distance) && 
          (distance <= max_distance)) {
        points.Insert(point);
      }
    }
  }
}



template <class PtrType>
int R3Kdtree<PtrType>::
FindAll(const R3Line& query_line, 
  RNScalar min_distance, RNScalar max_distance, 
  RNArray<PtrType>& points) const
{
  // Check root
  if (!root) return 0;

  // Search nodes recursively
  FindAll<R3Line>(root, bbox, 
    query_line,
    min_distance, max_distance, 
    points);

  // Return number of points
  return points.NEntries();
}



template <class PtrType>
int R3Kdtree<PtrType>::
FindAll(const R3Plane& query_plane, 
  RNScalar min_distance, RNScalar max_distance, 
  RNArray<PtrType>& points) const
{
  // Check root
  if (!root) return 0;

  // Search nodes recursively
  FindAll<R3Plane>(root, bbox, 
    query_plane,
    min_distance, max_distance, 
    points);

  // Return number of points
  return points.NEntries();
}



template <class PtrType>
int R3Kdtree<PtrType>::
FindAll(const R3Shape& query_shape, 
  RNScalar min_distance, RNScalar max_distance, 
  RNArray<PtrType>& points) const
{
  // Check root
  if (!root) return 0;

  // Search nodes recursively
  FindAll<R3Shape>(root, bbox, 
    query_shape,
    min_distance, max_distance, 
    points);

  // Return number of points
  return points.NEntries();
}



////////////////////////////////////////////////////////////////////////
// Internal tree creation functions
////////////////////////////////////////////////////////////////////////

template <class PtrType>
int R3Kdtree<PtrType>::
PartitionPoints(PtrType *points, int npoints, RNDimension dim, int imin, int imax)
{
  // Check range
  assert(imin <= imax);
  assert(imin >= 0);
  assert(imin < npoints);
  assert(imax >= 0);
  assert(imax < npoints);
  if (imin == imax) return imin;

  // Choose a coordinate at random to split upon
  int irand = (int) (imin + RNRandomScalar() * (imax - imin + 1));
  if (irand < imin) irand = imin;
  if (irand > imax) irand = imax;
  RNCoord split_coord = Position(points[irand])[dim];

  // Swap values at irand and imax
  PtrType swap = points[irand];
  points[irand] = points[imax];
  points[imax] = swap;

  // Partition points according to coordinate
  int split_index = imin;
  int middle_index = (imin + imax) / 2;
  for (int i = imin; i < imax; i++) {
    assert(split_index <= i);
    RNCoord coord = Position(points[i])[dim];
    if (coord < split_coord) {
      PtrType swap = points[split_index];
      points[split_index] = points[i];
      points[i] = swap;
      split_index++;
    }
    else if (coord == split_coord) {
      if (split_index < middle_index) {
        PtrType swap = points[split_index];
        points[split_index] = points[i];
        points[i] = swap;
        split_index++;
      }
    }
  }

  // Swap values at split_index and imax
  swap = points[split_index];
  points[split_index] = points[imax];
  points[imax] = swap;

  // Now split_index has value split_coord
  // All values to the left of split_index have values < split_coord
  // All values to the right of split_index have values >= split_coord

  // Recurse until we find the median
  if (split_index == imin) {
    if (imin >= npoints/2) return imin;
    else return PartitionPoints(points, npoints, dim, imin+1, imax);
  }
  else if (split_index == imax) {
    if (imax <= npoints/2) return imax;
    else return PartitionPoints(points, npoints, dim, imin, imax-1);
  }
  else {
    if (split_index == npoints/2) return split_index;
    else if (split_index > npoints/2) return PartitionPoints(points, npoints, dim, imin, split_index-1);
    else return PartitionPoints(points, npoints, dim, split_index+1, imax);
  }
}



template <class PtrType>
void R3Kdtree<PtrType>::
InsertPoints(R3KdtreeNode<PtrType> *node, const R3Box& node_box, PtrType *points, int npoints) 
{
  // Make sure node is an empty leaf
  assert(node);
  assert(node->children[0] == NULL);
  assert(node->children[1] == NULL);
  assert(node->npoints == 0);

  // Check number of points
  if (npoints < R3kdtree_max_points_per_node) {
    // Insert new points into leaf node and return
    for (int i = 0; i < npoints; i++) {
      node->points[node->npoints++] = points[i];
    }
  }
  else {
    // Find dimension to split along
    node->split_dimension = node_box.LongestAxis();
    
    // Partition points according to coordinates in split_dimension
    int split_index = PartitionPoints(points, npoints, node->split_dimension, 0, npoints-1);
    assert((split_index >= 0) && (split_index < npoints));

    // Determine split coordinate
    node->split_coordinate = Position(points[split_index])[node->split_dimension];
    assert(node->split_coordinate >= node_box[RN_LO][node->split_dimension]);
    assert(node->split_coordinate <= node_box[RN_HI][node->split_dimension]);

    // Construct children node boxes
    R3Box node0_box(node_box);
    R3Box node1_box(node_box);
    node0_box[RN_HI][node->split_dimension] = node->split_coordinate;
    node1_box[RN_LO][node->split_dimension] = node->split_coordinate;

    // Create children
    node->children[0] = new R3KdtreeNode<PtrType>(node);
    node->children[1] = new R3KdtreeNode<PtrType>(node);

    // Insert points into children
    InsertPoints(node->children[0], node0_box, points, split_index);
    InsertPoints(node->children[1], node1_box, &points[split_index], npoints - split_index);

    // Increment number of nodes
    nnodes += 2;
  }
}



template <class PtrType>
void R3Kdtree<PtrType>::
InsertPoint(R3KdtreeNode<PtrType> *node, const R3Box& node_box, PtrType point) 
{
  // Check if interior node
  if (node->children[0]) {
    // Inserting point into an interior node
    assert(node->children[1]);

    // Get point position
    R3Point position = Position(point);

    // Determine side of split
    RNLength side = position[node->split_dimension] - node->split_coordinate;
    if (side <= 0) {
      // Insert into negative side
      R3Box child_box(node_box);
      child_box[RN_HI][node->split_dimension] = node->split_coordinate;
      InsertPoint(node->children[0], child_box, point);
    }
    else {
      // Insert into positive side
      R3Box child_box(node_box);
      child_box[RN_LO][node->split_dimension] = node->split_coordinate;
      InsertPoint(node->children[1], child_box, point);
    }
  }
  else {
    // Inserting point into a leaf
    assert(!node->children[1]);
    
    // Check number of points
    if (node->npoints < R3kdtree_max_points_per_node) {
      // Insert new point into leaf node and return
      node->points[node->npoints++] = point;
    }
    else {
      // Find dimension to split along
      node->split_dimension = node_box.LongestAxis();
    
      // Partition points according to coordinates in split_dimension
      int split_index = PartitionPoints(node->points, node->npoints, node->split_dimension, 0, node->npoints-1);
      assert((split_index >= 0) && (split_index < node->npoints));

      // Determine split coordinate
      node->split_coordinate = Position(node->points[split_index])[node->split_dimension];
      assert(node->split_coordinate >= node_box[RN_LO][node->split_dimension]);
      assert(node->split_coordinate <= node_box[RN_HI][node->split_dimension]);

      // Construct children node boxes
      R3Box node0_box(node_box);
      R3Box node1_box(node_box);
      node0_box[RN_HI][node->split_dimension] = node->split_coordinate;
      node1_box[RN_LO][node->split_dimension] = node->split_coordinate;

      // Create children
      node->children[0] = new R3KdtreeNode<PtrType>(node);
      node->children[1] = new R3KdtreeNode<PtrType>(node);

      // Move points into children
      InsertPoints(node->children[0], node0_box, node->points, split_index);
      InsertPoints(node->children[1], node1_box, &node->points[split_index], node->npoints - split_index);
      node->npoints = 0;
      
      // Increment number of nodes
      nnodes += 2;

      // Insert point into this node
      InsertPoint(node, node_box, point);
    }
  }
}



template <class PtrType>
void R3Kdtree<PtrType>::
InsertPoint(PtrType point)
{
  // Insert point
  InsertPoint(root, bbox, point);

  // Update number of points
  npoints++;
}



template <class PtrType>
void R3Kdtree<PtrType>::
Outline(R3KdtreeNode<PtrType> *node, const R3Box& node_box) const
{
  // Draw kdtree nodes recursively
  if (node->children[0]) {
    assert(node->children[1]);
    assert(node->split_coordinate >= node_box[RN_LO][node->split_dimension]);
    assert(node->split_coordinate <= node_box[RN_HI][node->split_dimension]);
    R3Box child0_box(node_box);
    R3Box child1_box(node_box);
    child0_box[RN_HI][node->split_dimension] = node->split_coordinate;
    child1_box[RN_LO][node->split_dimension] = node->split_coordinate;
    Outline(node->children[0], child0_box);
    Outline(node->children[1], child1_box);
  }
  else {
    node_box.Outline();
  }
}



template <class PtrType>
void R3Kdtree<PtrType>::
Outline(void) const
{
  // Draw kdtree nodes recursively
  if (!root) return;
  Outline(root, bbox);
}



template <class PtrType>
int R3Kdtree<PtrType>::
PrintBalance(R3KdtreeNode<PtrType> *node, int depth) const
{
  // Check node
  if (!node) return 0;

  // Initialize number of decendents
  int ndecendents0 = 0;
  int ndecendents1 = 0;

  // Process interior node
  if (node->children[0] && node->children[1]) {
    // Print balance of children
    ndecendents0 = PrintBalance(node->children[0], depth+1);
    ndecendents1 = PrintBalance(node->children[1], depth+1);

    // Print balance of this node
    printf("%d", depth);
    for (int i = 0; i <= depth; i++) printf("  ");
    printf("I %d %d %g\n", ndecendents0, ndecendents1, (double) ndecendents0 / (double) ndecendents1);
  }
  else {
    printf("%d", depth);
    for (int i = 0; i <= depth; i++) printf("  ");
    printf("L %d\n", node->npoints);
  }

  // Return number of nodes rooted in this subtree
  return 1 + ndecendents0 + ndecendents1;
}



template <class PtrType>
void R3Kdtree<PtrType>::
PrintDebugInfo(void) const
{
  // Check root
  if (!root) return;
  PrintBalance(root, 0);
}



////////////////////////////////////////////////////////////////////////
// Finding the furthest one point
////////////////////////////////////////////////////////////////////////

#if 0  // This should work but is not tested

template <class PtrType>
void R3Kdtree<PtrType>::
FindFurthest(R3KdtreeNode<PtrType> *node, const R3Box& node_box, 
  PtrType query_point, const R3Point& query_position, 
  RNScalar min_distance_squared, RNScalar max_distance_squared,
  int (*IsCompatible)(PtrType, PtrType, void *), void *compatible_data, 
  PtrType& furthest_point, RNScalar& furthest_distance_squared) const
{
  // Check if node is interior
  if (node->children[0]) {
    assert(node->children[1]);

    // Find and check axial distances from point to node box
    RNLength dx, dy, dz;
    if (RNIsGreater(query_position.X(), node_box.XMax())) dx = query_position.X() - node_box.XMin();
    else if (RNIsLess(query_position.X(), node_box.XMin())) dx = node_box.XMax()- query_position.X();
    else dx = 0.0;
    RNLength dx_squared = dx * dx;
    if (RNIsGreater(query_position.Y(), node_box.YMax())) dy = query_position.Y() - node_box.YMin();
    else if (RNIsLess(query_position.Y(), node_box.YMin())) dy = node_box.YMax()- query_position.Y();
    else dy = 0.0;
    RNLength dy_squared = dy * dy;
    if (RNIsGreater(query_position.Z(), node_box.ZMax())) dz = query_position.Z() - node_box.ZMin();
    else if (RNIsLess(query_position.Z(), node_box.ZMin())) dz = node_box.ZMax()- query_position.Z();
    else dz = 0.0;
    RNLength dz_squared = dz * dz;
    
    // Find and check actual distance from point to node box
    RNLength distance_squared = 0;
    if ((dy == 0.0) && (dz == 0.0)) distance_squared = dx_squared;
    else if ((dx == 0.0) && (dz == 0.0)) distance_squared = dy_squared;
    else if ((dx == 0.0) && (dy == 0.0)) distance_squared = dz_squared;
    else distance_squared = dx_squared + dy_squared + dz_squared;
    if (distance_squared <= furthest_distance_squared) return;

    // Compute distance from point to split plane
    RNLength side = query_position[node->split_dimension] - node->split_coordinate;

    // Search children nodes
    if (side <= 0) {
      // Search negative side first
      R3Box child_box(node_box);
      child_box[RN_LO][node->split_dimension] = node->split_coordinate;
      FindFurthest(node->children[1], child_box, 
        query_point, query_position, 
        min_distance_squared, max_distance_squared, 
        IsCompatible, compatible_data,
        furthest_point, furthest_distance_squared);
      if (side*side > furthest_distance_squared) {
        R3Box child_box(node_box);
        child_box[RN_HI][node->split_dimension] = node->split_coordinate;
        FindFurthest(node->children[0], child_box, 
          query_point, query_position, 
          min_distance_squared, max_distance_squared, 
          IsCompatible, compatible_data,
          furthest_point, furthest_distance_squared);
      }
    }
    else {
      // Search positive side first
      R3Box child_box(node_box);
      child_box[RN_HI][node->split_dimension] = node->split_coordinate;
      FindFurthest(node->children[0], child_box, 
        query_point, query_position, 
        min_distance_squared, max_distance_squared, 
        IsCompatible, compatible_data,
        furthest_point, furthest_distance_squared);
      if (side*side > furthest_distance_squared) {
        R3Box child_box(node_box);
        child_box[RN_LO][node->split_dimension] = node->split_coordinate;
        FindFurthest(node->children[1], child_box, 
          query_point, query_position, 
          min_distance_squared, max_distance_squared, 
          IsCompatible, compatible_data,
          furthest_point, furthest_distance_squared);
      }
    }
  }
  else {
    for (int i = 0; i < node->npoints; i++) {
      PtrType point = node->points[i];
      RNLength distance_squared = R3SquaredDistance(query_position, Position(point));
      if ((distance_squared >= min_distance_squared) && 
          (distance_squared <= max_distance_squared) && 
          (distance_squared >= furthest_distance_squared)) {
        if (!IsCompatible || !query_point || IsCompatible(query_point, point, compatible_data)) {
          furthest_distance_squared = distance_squared;
          furthest_point = point;
        }
      }
    }
  }
}



template <class PtrType>
PtrType R3Kdtree<PtrType>::
FindFurthest(PtrType query_point, 
  RNScalar min_distance, RNScalar max_distance, 
  int (*IsCompatible)(PtrType, PtrType, void *), void *compatible_data, 
  RNScalar *furthest_distance) const
{
  // Check root
  if (!root) return NULL;

  // Use squared distances for efficiency
  if (max_distance < 0) return NULL;
  if (min_distance < 0) min_distance = 0;
  RNLength min_distance_squared = min_distance * min_distance;
  RNLength max_distance_squared = max_distance * max_distance;

  // Initialize furthest point 
  PtrType furthest_point = NULL;
  RNLength furthest_distance_squared = max_distance_squared;

  // Search nodes recursively
  FindFurthest(root, bbox, 
    query_point, Position(query_point),
    min_distance_squared, max_distance_squared, 
    IsCompatible, compatible_data, 
    furthest_point, furthest_distance_squared);

  // Return furthest distance
  if (furthest_distance) *furthest_distance = sqrt(furthest_distance_squared);

  // Return furthest point
  return furthest_point;
}



template <class PtrType>
PtrType R3Kdtree<PtrType>::
FindFurthest(PtrType query_point, 
  RNScalar min_distance, RNScalar max_distance, 
  RNScalar *furthest_distance) const
{
  // Find the furthest point
  return FindFurthest(query_point, min_distance, max_distance, furthest_distance);
}



template <class PtrType>
PtrType R3Kdtree<PtrType>::
FindFurthest(const R3Point& query_position, 
  RNScalar min_distance, RNScalar max_distance, 
  RNScalar *furthest_distance) const
{
  // Check root
  if (!root) return NULL;

  // Use squared distances for efficiency
  if (max_distance < 0) return NULL;
  if (min_distance < 0) min_distance = 0;
  RNLength min_distance_squared = min_distance * min_distance;
  RNLength max_distance_squared = max_distance * max_distance;

  // Initialize furthest point 
  PtrType furthest_point = NULL;
  RNLength furthest_distance_squared = max_distance_squared;

  // Search nodes recursively
  FindFurthest(root, bbox, 
    NULL, query_position, 
    min_distance_squared, max_distance_squared, 
    NULL, NULL, 
    furthest_point, furthest_distance_squared);

  // Return furthest distance
  if (furthest_distance) *furthest_distance = sqrt(furthest_distance_squared);

  // Return furthest point
  return furthest_point;
}

#endif



#endif


