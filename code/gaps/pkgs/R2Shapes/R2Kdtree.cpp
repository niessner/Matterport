// Source file for R2Kdtree class

#ifndef __R2KDTREE__C__
#define __R2KDTREE__C__




////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "R2Shapes/R2Shapes.h"





////////////////////////////////////////////////////////////////////////
// Constant definitions
////////////////////////////////////////////////////////////////////////

static const int R2kdtree_max_points_per_node = 32;





////////////////////////////////////////////////////////////////////////
// Node class definition
////////////////////////////////////////////////////////////////////////

// Node declaration

template <class PtrType>
class R2KdtreeNode {
public:
  R2KdtreeNode(R2KdtreeNode<PtrType> *parent);

public:
  class R2KdtreeNode<PtrType> *parent;
  class R2KdtreeNode<PtrType> *children[2];
  RNScalar split_coordinate;
  RNDimension split_dimension;
  PtrType points[R2kdtree_max_points_per_node];
  int npoints;
};



// Node constructor

template <class PtrType>
R2KdtreeNode<PtrType>::
R2KdtreeNode(R2KdtreeNode<PtrType> *parent)
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
R2Kdtree<PtrType>::
R2Kdtree(const R2Box& bbox, int position_offset)
  : bbox(bbox),
    position_offset(position_offset),
    position_callback(NULL),
    position_callback_data(NULL),
    nnodes(0)
{
  // Create root node
  root = new R2KdtreeNode<PtrType>(NULL);
  assert(root);
}



template <class PtrType>
R2Kdtree<PtrType>::
R2Kdtree(const RNArray<PtrType>& points, int position_offset)
  : position_offset(position_offset),
    position_callback(NULL),
    position_callback_data(NULL),
    nnodes(1)
{
  // Create root node
  root = new R2KdtreeNode<PtrType>(NULL);
  assert(root);

  // Determine bounding box 
  bbox = R2null_box;
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
R2Kdtree<PtrType>::
R2Kdtree(const RNArray<PtrType>& points, R2Point (*position_callback)(PtrType, void *), void *position_callback_data)
  : position_offset(-1),
    position_callback(position_callback),
    position_callback_data(position_callback_data),
    nnodes(1)
{
  // Create root node
  root = new R2KdtreeNode<PtrType>(NULL);
  assert(root);

  // Determine bounding box 
  bbox = R2null_box;
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
R2Kdtree<PtrType>::
~R2Kdtree(void)
{
  // Check root
  if (!root) return;

  // Traverse tree deleting nodes
  RNArray<R2KdtreeNode<PtrType> *> stack;
  stack.InsertTail(root);
  while (!stack.IsEmpty()) {
    R2KdtreeNode<PtrType> *node = stack.Tail();
    stack.RemoveTail();
    if (node->children[0]) stack.Insert(node->children[0]);
    if (node->children[1]) stack.Insert(node->children[1]);
    delete node;
  }
}



template <class PtrType>
const R2Box& R2Kdtree<PtrType>::
BBox(void) const
{
  // Return bounding box of the whole KD tree
  return bbox;
}



template <class PtrType>
int R2Kdtree<PtrType>::
NNodes(void) const
{
  // Return number of nodes
  return nnodes;
}



#if 0

template <class PtrType>
int R2Kdtree<PtrType>::
NEntries(void) const
{
  // Check root
  if (!root) return 0;

  // Traverse tree to count number of points
  int npoints = 0;
  RNArray<R2KdtreeNode<PtrType> *> stack;
  stack.InsertTail(root);
  while (!stack.IsEmpty()) {
    R2KdtreeNode<PtrType> *node = stack.Tail();
    stack.RemoveTail();
    npoints += node->npoints;
    if (node->children[0]) stack.Insert(node->children[0]);
    if (node->children[1]) stack.Insert(node->children[1]);
  }

  // Return total number of points
  return npoints;
}

#endif



template <class PtrType>
void R2Kdtree<PtrType>::
FindClosest(R2KdtreeNode<PtrType> *node, const R2Box& node_box, const R2Point& position, 
  RNScalar min_distance_squared, RNScalar max_distance_squared,
  PtrType& closest_point, RNScalar& closest_distance_squared) const
{
  // Check if node is interior
  if (node->children[0]) {
    assert(node->children[1]);

    // Find and check axial distances from point to node box
    RNLength dx, dy;
    if (RNIsGreater(position.X(), node_box.XMax())) dx = position.X() - node_box.XMax();
    else if (RNIsLess(position.X(), node_box.XMin())) dx = node_box.XMin()- position.X();
    else dx = 0.0;
    RNLength dx_squared = dx * dx;
    if (dx_squared >= closest_distance_squared) return;
    if (RNIsGreater(position.Y(), node_box.YMax())) dy = position.Y() - node_box.YMax();
    else if (RNIsLess(position.Y(), node_box.YMin())) dy = node_box.YMin()- position.Y();
    else dy = 0.0;
    RNLength dy_squared = dy * dy;
    if (dy_squared >= closest_distance_squared) return;
    
    // Find and check actual distance from point to node box
    RNLength distance_squared = dx_squared + dy_squared;
    if (distance_squared >= closest_distance_squared) return;

    // Compute distance from point to split plane
    RNLength side = position[node->split_dimension] - node->split_coordinate;

    // Search children nodes
    if (side <= 0) {
      // Search negative side first
      R2Box child_box(node_box);
      child_box[RN_HI][node->split_dimension] = node->split_coordinate;
      FindClosest(node->children[0], child_box, position, 
        min_distance_squared, max_distance_squared, 
        closest_point, closest_distance_squared);
      if (side*side < closest_distance_squared) {
        R2Box child_box(node_box);
        child_box[RN_LO][node->split_dimension] = node->split_coordinate;
        FindClosest(node->children[1], child_box, position, 
          min_distance_squared, max_distance_squared, 
          closest_point, closest_distance_squared);
      }
    }
    else {
      // Search positive side first
      R2Box child_box(node_box);
      child_box[RN_LO][node->split_dimension] = node->split_coordinate;
      FindClosest(node->children[1], child_box, position, 
        min_distance_squared, max_distance_squared, 
        closest_point, closest_distance_squared);
      if (side*side < closest_distance_squared) {
        R2Box child_box(node_box);
        child_box[RN_HI][node->split_dimension] = node->split_coordinate;
        FindClosest(node->children[0], child_box, position, 
          min_distance_squared, max_distance_squared, 
          closest_point, closest_distance_squared);
      }
    }
  }
  else {
    for (int i = 0; i < node->npoints; i++) {
      PtrType point = node->points[i];
      R2Vector v = position - Position(point);
      RNLength distance_squared = v.Dot(v);
      if ((distance_squared >= min_distance_squared) && 
         (distance_squared <= closest_distance_squared)) {
        closest_distance_squared = distance_squared;
        closest_point = point;
      }
    }
  }
}



template <class PtrType>
PtrType R2Kdtree<PtrType>::
FindClosest(const R2Point& position, 
  RNScalar min_distance, RNScalar max_distance, 
  RNScalar *closest_distance) const
{
  // Check root
  if (!root) return NULL;

  // Use squared distances for efficiency
  RNLength min_distance_squared = min_distance * min_distance;
  RNLength max_distance_squared = max_distance * max_distance;

  // Initialize nearest point 
  PtrType nearest_point = NULL;
  RNLength nearest_distance_squared = max_distance_squared;

  // Search nodes recursively
  FindClosest(root, bbox, position, 
    min_distance_squared, max_distance_squared, 
    nearest_point, nearest_distance_squared);

  // Return closest distance
  if (closest_distance) *closest_distance = sqrt(nearest_distance_squared);

  // Return closest point
  return nearest_point;
}



template <class PtrType>
PtrType R2Kdtree<PtrType>::
FindClosest(PtrType point, 
  RNScalar min_distance, RNScalar max_distance, 
  RNScalar *closest_distance) const
{
  // Find the closest point
  return FindClosest(Position(point), min_distance, max_distance, closest_distance);
}



template <class PtrType>
void R2Kdtree<PtrType>::
FindAll(R2KdtreeNode<PtrType> *node, const R2Box& node_box, const R2Point& position, 
  RNScalar min_distance_squared, RNScalar max_distance_squared, RNArray<PtrType>& points) const
{
  // Check if node is interior
  if (node->children[0]) {
    assert(node->children[1]);

    // Find and check axial distances from point to node box
    RNLength dx, dy;
    if (RNIsGreater(position.X(), node_box.XMax())) dx = position.X() - node_box.XMax();
    else if (RNIsLess(position.X(), node_box.XMin())) dx = node_box.XMin()- position.X();
    else dx = 0.0;
    RNLength dx_squared = dx * dx;
    if (dx_squared > max_distance_squared) return;
    if (RNIsGreater(position.Y(), node_box.YMax())) dy = position.Y() - node_box.YMax();
    else if (RNIsLess(position.Y(), node_box.YMin())) dy = node_box.YMin()- position.Y();
    else dy = 0.0;
    RNLength dy_squared = dy * dy;
    if (dy_squared > max_distance_squared) return;
    
    // Find and check actual distance from point to node box
    RNLength distance_squared = dx_squared + dy_squared;
    if (distance_squared > max_distance_squared) return;

    // Compute distance from point to split plane
    RNLength side = position[node->split_dimension] - node->split_coordinate;

    // Search children nodes
    if ((side <= 0) || (side*side <= max_distance_squared)) {
      // Search negative side 
      R2Box child_box(node_box);
      child_box[RN_HI][node->split_dimension] = node->split_coordinate;
      FindAll(node->children[0], child_box, position, 
        min_distance_squared, max_distance_squared, points);
    }
    if ((side >= 0) || (side*side <= max_distance_squared)) {
      R2Box child_box(node_box);
      child_box[RN_LO][node->split_dimension] = node->split_coordinate;
      FindAll(node->children[1], child_box, position, 
        min_distance_squared, max_distance_squared, points);
    }
  }
  else {
    // Search points
    for (int i = 0; i < node->npoints; i++) {
      PtrType point = node->points[i];
      R2Vector v = position - Position(point);
      RNLength distance_squared = v.Dot(v);
      if ((distance_squared >= min_distance_squared) && 
          (distance_squared <= max_distance_squared)) {
        points.Insert(point);
      }
    }
  }
}



template <class PtrType>
int R2Kdtree<PtrType>::
FindAll(const R2Point& position, RNScalar min_distance, RNScalar max_distance, RNArray<PtrType>& points) const
{
  // Check root
  if (!root) return 0;

  // Use squared distances for efficiency
  RNLength min_distance_squared = min_distance * min_distance;
  RNLength max_distance_squared = max_distance * max_distance;

  // Search nodes recursively
  FindAll(root, bbox, position, min_distance_squared, max_distance_squared, points);

  // Return number of points
  return points.NEntries();
}



template <class PtrType>
int R2Kdtree<PtrType>::
FindAll(PtrType point, RNScalar min_distance, RNScalar max_distance, RNArray<PtrType>& points) const
{
  // Find all within some distance
  return FindAll(Position(point), min_distance, max_distance, points);
}




////////////////////////////////////////////////////////////////////////
// Internal tree creation functions
////////////////////////////////////////////////////////////////////////

template <class PtrType>
int R2Kdtree<PtrType>::
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

  // Partition points according to coordinate
  int left_index = imin;
  int right_index = imax;
  int different_coord = 0;
  while (left_index <= right_index) {
    while (left_index <= right_index) {
      RNCoord left_coord = Position(points[left_index])[dim];
      if (left_coord != split_coord) different_coord = 1;
      if (left_coord > split_coord) break; 
      left_index++;
    }
    while (left_index <= right_index) {
      RNCoord right_coord = Position(points[right_index])[dim];
      if (right_coord != split_coord) different_coord = 1;
      if (right_coord <= split_coord) break;
      right_index--;
    }
    if (left_index < right_index) {
      PtrType swap = points[right_index];
      points[right_index] = points[left_index];
      points[left_index] = swap;
      left_index++;
      right_index--;
    }
  }

  // Check for sanity
  assert(left_index == right_index+1);
  assert(left_index > imin);
  assert(right_index <= imax);

  // Recurse until we find the median
  if (!different_coord) return (imin + imax) / 2;
  else if (left_index > npoints/2) return PartitionPoints(points, npoints, dim, imin, left_index-1);
  else return PartitionPoints(points, npoints, dim, left_index, imax);
}



template <class PtrType>
void R2Kdtree<PtrType>::
InsertPoints(R2KdtreeNode<PtrType> *node, const R2Box& node_box, PtrType *points, int npoints) 
{
  // Make sure node is an empty leaf
  assert(node);
  assert(node->children[0] == NULL);
  assert(node->children[1] == NULL);
  assert(node->npoints == 0);

  // Check number of points
  if (npoints <= R2kdtree_max_points_per_node) {
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
    R2Box node0_box(node_box);
    R2Box node1_box(node_box);
    node0_box[RN_HI][node->split_dimension] = node->split_coordinate;
    node1_box[RN_LO][node->split_dimension] = node->split_coordinate;

    // Create children
    node->children[0] = new R2KdtreeNode<PtrType>(node);
    node->children[1] = new R2KdtreeNode<PtrType>(node);

    // Insert points into children
    InsertPoints(node->children[0], node0_box, points, split_index);
    InsertPoints(node->children[1], node1_box, &points[split_index], npoints - split_index);

    // Increment number of nodes
    nnodes += 2;
  }
}



template <class PtrType>
void R2Kdtree<PtrType>::
Outline(R2KdtreeNode<PtrType> *node, const R2Box& node_box) const
{
  // Draw kdtree nodes recursively
  if (node->children[0]) {
    assert(node->children[1]);
    assert(node->split_coordinate >= node_box[RN_LO][node->split_dimension]);
    assert(node->split_coordinate <= node_box[RN_HI][node->split_dimension]);
    R2Box child0_box(node_box);
    R2Box child1_box(node_box);
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
void R2Kdtree<PtrType>::
Outline(void) const
{
  // Draw kdtree nodes recursively
  if (!root) return;
  Outline(root, bbox);
}



template <class PtrType>
int R2Kdtree<PtrType>::
PrintBalance(R2KdtreeNode<PtrType> *node, int depth) const
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
void R2Kdtree<PtrType>::
PrintDebugInfo(void) const
{
  // Check root
  if (!root) return;
  PrintBalance(root, 0);
}


#endif


