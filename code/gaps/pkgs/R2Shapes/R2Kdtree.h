// Include file for KDTree class

#ifndef __R2KDTREE__H__
#define __R2KDTREE__H__



// Node declaration

template <class PtrType>
class R2KdtreeNode;



// Class declaration

template <class PtrType>
class R2Kdtree {
public:
  // Constructor/destructors
  R2Kdtree(const R2Box& bbox, int position_offset = 0);
  R2Kdtree(const RNArray<PtrType>& points, int position_offset = 0);
  R2Kdtree(const RNArray<PtrType>& points, R2Point (*position_callback)(PtrType, void *), void *data);
  ~R2Kdtree(void);

  // Property functions
  const R2Box& BBox(void) const;
  int NNodes(void) const;

  // Search for closest to a point
  PtrType FindClosest(PtrType point, RNLength min_distance = 0, RNLength max_distance = RN_INFINITY, RNLength *closest_distance = NULL) const;
  PtrType FindClosest(const R2Point& position, RNLength min_distance = 0, RNLength max_distance = RN_INFINITY, RNLength *closest_distance = NULL) const;

  // Search for all within some distance to a point
  int FindAll(PtrType point, RNLength min_distance, RNLength max_distance, RNArray<PtrType>& points) const;
  int FindAll(const R2Point& position, RNLength min_distance, RNLength max_distance, RNArray<PtrType>& points) const;

public:
  // Internal search functions
  void FindClosest(R2KdtreeNode<PtrType> *node, const R2Box& node_box, const R2Point& position, 
    RNLength min_distance_squared, RNLength max_distance_squared, 
    PtrType& closest_point, RNLength& closest_distance_squared) const;
  void FindAll(R2KdtreeNode<PtrType> *node, const R2Box& node_box, const R2Point& position, 
    RNLength min_distance_squared, RNLength max_distance_squared, RNArray<PtrType>& points) const;

  // Internal manipulation functions
  void InsertPoints(R2KdtreeNode<PtrType> *node, const R2Box& node_box, PtrType *points, int npoints);
  int PartitionPoints(PtrType *points, int npoints, RNDimension dim, int imin, int imax);

  // Internal visualization functions
  void Outline(R2KdtreeNode<PtrType> *node, const R2Box& bbox) const;
  void Outline(void) const;

  // Internal debugging functions
  void PrintDebugInfo(void) const;
  int PrintBalance(R2KdtreeNode<PtrType> *node, int depth) const;

  // Internal point extraction function
  const R2Point Position(PtrType point) const { 
    if (position_offset >= 0) return *((R2Point *) ((unsigned char *) point + position_offset)); 
    else if (position_callback) return (*position_callback)(point, position_callback_data);
    else { fprintf(stderr, "Invalid position callback\n"); abort(); return R2null_point; }
  } ;

public:
  // Internal data
  R2Box bbox;
  int position_offset;
  R2Point (*position_callback)(PtrType, void *);
  void *position_callback_data;
  R2KdtreeNode<PtrType> *root;
  int nnodes;
};



// Include templated definitions

#include "R2Kdtree.cpp"


#endif
