// Source file for mesh search tree class



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "R3Shapes/R3Shapes.h"





////////////////////////////////////////////////////////////////////////
// Constant definitions
////////////////////////////////////////////////////////////////////////

static const int max_faces_per_node = 128;
static const RNScalar max_area_ratio = 0.1;
static const int max_depth = 64;



////////////////////////////////////////////////////////////////////////
// Face class definition
////////////////////////////////////////////////////////////////////////

class R3MeshSearchTreeFace {
public:
  R3MeshSearchTreeFace(R3Mesh *mesh, R3MeshFace *face) 
  : face(face), area(mesh->FaceArea(face)), reference_count(0), mark(0) {};

public:
  R3MeshFace *face;
  RNArea area;
  int reference_count;
  RNMark mark;
};



////////////////////////////////////////////////////////////////////////
// Node class definition
////////////////////////////////////////////////////////////////////////

// Node declaration

class R3MeshSearchTreeNode {
public:
  R3MeshSearchTreeNode(R3MeshSearchTreeNode *parent = NULL)
    : parent(parent), split_coordinate(0), split_dimension(0), big_faces(), small_faces()
  { children[0] = NULL; children[1] = NULL; };

public:
  class R3MeshSearchTreeNode *parent;
  class R3MeshSearchTreeNode *children[2];
  RNScalar split_coordinate;
  RNDimension split_dimension;
  RNArray<R3MeshSearchTreeFace *> big_faces;
  RNArray<R3MeshSearchTreeFace *> small_faces;
};




////////////////////////////////////////////////////////////////////////
// Constructor/destructor functions
////////////////////////////////////////////////////////////////////////

R3MeshSearchTree::
R3MeshSearchTree(R3Mesh *mesh)
  : mesh(mesh),
    nnodes(1),
    mark(1)
{
  // Create root 
  root = new R3MeshSearchTreeNode(NULL);
  assert(root);

  // Insert faces into tree
  for (int i = 0; i < mesh->NFaces(); i++) {
    R3MeshFace *face = mesh->Face(i);
    InsertFace(face);
  }
}



R3MeshSearchTree::
~R3MeshSearchTree(void)
{
  // Empty tree
  Empty();
}



////////////////////////////////////////////////////////////////////////
// Insert/delete functions
////////////////////////////////////////////////////////////////////////

static RNBoolean
R3Intersects(R3Mesh *mesh, R3MeshFace *face, const R3Box& box)
{
  // Check triangle bounding box and plane
  if (R3Contains(box, mesh->FaceBBox(face))) return TRUE;
  if (!R3Intersects(mesh->FaceBBox(face), box)) return FALSE;
  if (!R3Intersects(mesh->FacePlane(face), box)) return FALSE;

  // Make polygon
  R3Point points[16];
  points[0] = mesh->VertexPosition(mesh->VertexOnFace(face, 0));
  points[1] = mesh->VertexPosition(mesh->VertexOnFace(face, 1));
  points[2] = mesh->VertexPosition(mesh->VertexOnFace(face, 2));
  int npoints = 3;

  // Clip polygon against each side
  for (RNDirection dir = RN_LO; dir <= RN_HI; dir++) {
    for (RNDimension dim = RN_X; dim <= RN_Z; dim++) {
      assert(npoints < 15);
      points[npoints] = points[0];
      for (int i = 0; i < npoints; i++) {
        R3Point& p1 = points[i];
        R3Point& p2 = points[i+1];
        RNScalar d1 = p1[dim] - box[dir][dim];
        RNScalar d2 = p2[dim] - box[dir][dim];
        if (dir == RN_LO) { d1 = -d1; d2 = -d2; }
        if (d1 < 0) { // Inside
          if (d2 > 0) { // Outside
            // Insert a point at crossing
            RNScalar denom = d2 + -d1;
            R3Point p = (d2/denom)*p1 + (-d1/denom)*p2;
            for (int j = npoints; j > i; j--) points[j+1] = points[j];
            points[i+1] = p;
            npoints++;
            i += 2;
          }
        }
        else if (d1 > 0) { // Outside
          if (d2 < 0) { // Inside
            // Replace p1 with point at crossing
            RNScalar denom = -d2 + d1;
            R3Point p = (-d2/denom)*p1 + (d1/denom)*p2;
            points[i] = p;
          }
          else {
            // Remove p1
            for (int j = i; j < npoints; j++) points[j] = points[j+1];
            npoints--;
            i--;
          }
        }
      }

      // Check number of points 
      if (npoints == 0) return FALSE; 
      assert(npoints < 16);
      assert(npoints > 0);
    }
  }

  // Triangle survived all clips
  return TRUE;
}



void R3MeshSearchTree::
InsertFace(R3MeshFace *face)
{
  // Check if face intersects box
  if (!R3Intersects(mesh, face, BBox())) return;

  // Create container
  R3MeshSearchTreeFace *face_container = new R3MeshSearchTreeFace(mesh, face);
  assert(face_container);

  // Insert face into root
  InsertFace(face_container, root, BBox(), 0);
}



void R3MeshSearchTree::
InsertFace(R3MeshSearchTreeFace *face, R3MeshSearchTreeNode *node, const R3Box& node_box, int depth) 
{
  // Check if face intersects box
  if (!R3Intersects(mesh, face->face, node_box)) return;

  // Check if interior node
  if (node->children[0]) {
    // Interior node -- Insert into children
    assert(node->children[1]);
    const R3Box& face_box = mesh->FaceBBox(face->face);
    if (face_box[RN_LO][node->split_dimension] <= node->split_coordinate) {
      R3Box node0_box(node_box);
      node0_box[RN_HI][node->split_dimension] = node->split_coordinate;
      InsertFace(face, node->children[0], node0_box, depth + 1);
    }
    if (face_box[RN_HI][node->split_dimension] >= node->split_coordinate) {
      R3Box node1_box(node_box);
      node1_box[RN_LO][node->split_dimension] = node->split_coordinate;
      InsertFace(face, node->children[1], node1_box, depth + 1);
    }
  }
  else {
    // Check face area 
    RNScalar node_diagonal = node_box.DiagonalLength();
    if (node_diagonal == 0) return;
    RNScalar node_area = node_diagonal * node_diagonal;
    RNScalar area_ratio = face->area / node_area;
    if ((area_ratio >= max_area_ratio) || 
        (depth >= max_depth)) {
      // Face is too big/deep to be sorted into children, insert into big faces list
      node->big_faces.Insert(face);
      assert(face->reference_count >= 0);
      face->reference_count++;
    }
    else {
      // Leaf node -- Check if there is room for this face
      if (node->small_faces.NEntries() < max_faces_per_node) {
        // Simply insert face into list
        node->small_faces.Insert(face);
        face->reference_count++;
      }
      else {
        // Create two children 
        node->children[0] = new R3MeshSearchTreeNode(node);
        node->children[1] = new R3MeshSearchTreeNode(node);
        node->split_dimension = node_box.LongestAxis();
        node->split_coordinate = node_box.AxisCenter(node->split_dimension);
        nnodes += 2;
      
        // Re-insert faces into subtree
        InsertFace(face, node, node_box, depth+1);
        for (int i = 0; i < node->small_faces.NEntries(); i++) {
          InsertFace(node->small_faces[i], node, node_box, depth+1);
        }
      
        // Clear out faces from node that is now interior
        for (int i = 0; i < node->small_faces.NEntries(); i++) {
          assert(node->small_faces[i]->reference_count > 0);
          node->small_faces[i]->reference_count--;
        }
        node->small_faces.Empty();
      }
    }
  }
}



void R3MeshSearchTree::
Empty(R3MeshSearchTreeNode *node)
{
  // Delete children
  if (node->children[0]) {
    Empty(node->children[0]);
    delete node->children[0];
    node->children[0] = NULL;
  }
  if (node->children[1]) {
    Empty(node->children[1]);
    delete node->children[1];
    node->children[1] = NULL;
  }

  // Delete small face containers
  for (int i = 0; i < node->small_faces.NEntries(); i++) {
    R3MeshSearchTreeFace *face = node->small_faces[i];
    assert(face->reference_count > 0);
    if (--(face->reference_count) == 0) delete face;
  }

  // Delete big face containers
  for (int i = 0; i < node->big_faces.NEntries(); i++) {
    R3MeshSearchTreeFace *face = node->big_faces[i];
    assert(face->reference_count > 0);
    if (--(face->reference_count) == 0) delete face;
  }

  // Empty face containers
  node->small_faces.Empty();
  node->big_faces.Empty();
}


void R3MeshSearchTree::
Empty(void)
{
  // Check root
  if (root) Empty(root);
}



////////////////////////////////////////////////////////////////////////
// Closest point search functions
////////////////////////////////////////////////////////////////////////

void R3MeshSearchTree::
FindClosest(const R3Point& query_position, const R3Vector& query_normal, R3MeshIntersection& closest, 
  RNScalar min_distance_squared, RNScalar& max_distance_squared, 
  int (*IsCompatible)(const R3Point&, const R3Vector&, R3Mesh *, R3MeshFace *, void *), void *compatible_data,
  R3MeshFace *face) const
{
  // Check distance to plane
  const R3Plane& plane = mesh->FacePlane(face);
  RNScalar plane_signed_distance = R3SignedDistance(plane, query_position);
  RNScalar plane_distance_squared = plane_signed_distance * plane_signed_distance;
  if (plane_distance_squared >= max_distance_squared) return;

  // Check distance to bounding box
  RNScalar bbox_distance_squared = DistanceSquared(query_position, mesh->FaceBBox(face), max_distance_squared);
  if (bbox_distance_squared >= max_distance_squared) return;

  // Check compatibility 
  if (IsCompatible) {
    if (!(*IsCompatible)(query_position, query_normal, mesh, face, compatible_data)) return;
  }

  // Get face vertices
  R3MeshVertex *v0 = mesh->VertexOnFace(face, 0);
  R3MeshVertex *v1 = mesh->VertexOnFace(face, 1);
  R3MeshVertex *v2 = mesh->VertexOnFace(face, 2);

  // Get vertex positions
  const R3Point& p0 = mesh->VertexPosition(v0);
  const R3Point& p1 = mesh->VertexPosition(v1);
  const R3Point& p2 = mesh->VertexPosition(v2);

  // Project query point onto face plane
  const R3Vector& face_normal = mesh->FaceNormal(face);
  R3Point plane_point = query_position - plane_signed_distance * face_normal;

  // Check sides of edges
  R3Vector e0 = p1 - p0;
  e0.Normalize();
  R3Vector n0 = mesh->FaceNormal(face) % e0;
  R3Plane s0(p0, n0);
  RNScalar b0 = R3SignedDistance(s0, plane_point);
  R3Vector e1 = p2 - p1;
  e1.Normalize();
  R3Vector n1 = mesh->FaceNormal(face) % e1;
  R3Plane s1(p1, n1);
  RNScalar b1 = R3SignedDistance(s1, plane_point);
  R3Vector e2 = p0 - p2;
  e2.Normalize();
  R3Vector n2 = mesh->FaceNormal(face) % e2;
  R3Plane s2(p2, n2);
  RNScalar b2 = R3SignedDistance(s2, plane_point);

  // Consider plane_point's position in relation to edges of the triangle
  if ((b0 >= 0) && (b1 >= 0) && (b2 >= 0)) {
    // Point is inside face
    if (plane_distance_squared >= min_distance_squared) {
      closest.type = R3_MESH_FACE_TYPE;
      closest.face = face;
      closest.point = plane_point;
      max_distance_squared = plane_distance_squared;
    }
  }
  else {
    // Point is outside face -- check each edge
    if (b0 < 0) {
      // Outside edge0
      R3Vector edge_vector = p1 - p0;
      RNScalar edge_length = edge_vector.Length();
      if (edge_length > 0) {
        edge_vector /= edge_length;
        R3Vector point_vector = plane_point - p0;
        RNScalar t = edge_vector.Dot(point_vector);
        if (t <= 0) {
          RNScalar distance_squared = DistanceSquared(query_position, p0);
          if ((distance_squared >= min_distance_squared) && (distance_squared < max_distance_squared)) {
            closest.type = R3_MESH_VERTEX_TYPE;
            closest.vertex = v0;
            closest.point = p0;
            max_distance_squared = distance_squared;
          }
        }
        else if (t >= edge_length) {
          RNScalar distance_squared = DistanceSquared(query_position, p1);
          if ((distance_squared >= min_distance_squared) && (distance_squared < max_distance_squared)) {
            closest.type = R3_MESH_VERTEX_TYPE;
            closest.vertex = v1;
            closest.point = p1;
            max_distance_squared = distance_squared;
          }
        }
        else {
          R3Point point = p0 + t * edge_vector;
          RNScalar distance_squared = DistanceSquared(query_position, point);
          if ((distance_squared >= min_distance_squared) && (distance_squared < max_distance_squared)) {
            closest.type = R3_MESH_EDGE_TYPE;
            closest.edge = mesh->EdgeOnFace(face, 0);
            closest.point = point;
            max_distance_squared = distance_squared;
          }
        }
      }
    }
    if (b1 < 0) {
      // Outside edge1
      R3Vector edge_vector = p2 - p1;
      RNScalar edge_length = edge_vector.Length();
      if (edge_length > 0) {
        edge_vector /= edge_length;
        R3Vector point_vector = plane_point - p1;
        RNScalar t = edge_vector.Dot(point_vector);
        if (t <= 0) {
          RNScalar distance_squared = DistanceSquared(query_position, p1);
          if ((distance_squared >= min_distance_squared) && (distance_squared < max_distance_squared)) {
            closest.type = R3_MESH_VERTEX_TYPE;
            closest.vertex = v1;
            closest.point = p1;
            max_distance_squared = distance_squared;
          }
        }
        else if (t >= edge_length) {
          RNScalar distance_squared = DistanceSquared(query_position, p2);
          if ((distance_squared >= min_distance_squared) && (distance_squared < max_distance_squared)) {
            closest.type = R3_MESH_VERTEX_TYPE;
            closest.vertex = v2;
            closest.point = p2;
            max_distance_squared = distance_squared;
          }
        }
        else {
          R3Point point = p1 + t * edge_vector;
          RNScalar distance_squared = DistanceSquared(query_position, point);
          if ((distance_squared >= min_distance_squared) && (distance_squared < max_distance_squared)) {
            closest.type = R3_MESH_EDGE_TYPE;
            closest.edge = mesh->EdgeOnFace(face, 1);
            closest.point = point;
            max_distance_squared = distance_squared;
          }
        }
      }
    }
    if (b2 < 0) {
      // Outside edge2
      R3Vector edge_vector = p0 - p2;
      RNScalar edge_length = edge_vector.Length();
      if (edge_length > 0) {
        edge_vector /= edge_length;
        R3Vector point_vector = plane_point - p2;
        RNScalar t = edge_vector.Dot(point_vector);
        if (t <= 0) {
          RNScalar distance_squared = DistanceSquared(query_position, p2);
          if ((distance_squared >= min_distance_squared) && (distance_squared < max_distance_squared)) {
            closest.type = R3_MESH_VERTEX_TYPE;
            closest.vertex = v2;
            closest.point = p2;
            max_distance_squared = distance_squared;
          }
        }
        else if (t >= edge_length) {
          RNScalar distance_squared = DistanceSquared(query_position, p0);
          if ((distance_squared >= min_distance_squared) && (distance_squared < max_distance_squared)) {
            closest.type = R3_MESH_VERTEX_TYPE;
            closest.vertex = v0;
            closest.point = p0;
            max_distance_squared = distance_squared;
          }
        }
        else {
          R3Point point = p2 + t * edge_vector;
          RNScalar distance_squared = DistanceSquared(query_position, point);
          if ((distance_squared >= min_distance_squared) && (distance_squared < max_distance_squared)) {
            closest.type = R3_MESH_EDGE_TYPE;
            closest.edge = mesh->EdgeOnFace(face, 2);
            closest.point = point;
            max_distance_squared = distance_squared;
          }
        }
      }
    }
  }
}



void R3MeshSearchTree::
FindClosest(const R3Point& query_position, const R3Vector& query_normal, R3MeshIntersection& closest, 
  RNScalar min_distance_squared, RNScalar& max_distance_squared, 
  int (*IsCompatible)(const R3Point&, const R3Vector&, R3Mesh *, R3MeshFace *, void *), void *compatible_data,
  R3MeshSearchTreeNode *node, const R3Box& node_box) const
{
  // Compute distance (squared) from query point to node bbox
  RNScalar distance_squared = DistanceSquared(query_position, node_box, max_distance_squared);
  if (distance_squared >= max_distance_squared) return;

  // Update based on distance to each big face
  for (int i = 0; i < node->big_faces.NEntries(); i++) {
    // Get face container and check mark
    R3MeshSearchTreeFace *face_container = node->big_faces[i];
    if (face_container->mark == mark) continue;
    face_container->mark = mark;
  
    // Find closest point in mesh face
    FindClosest(query_position, query_normal, closest, 
      min_distance_squared, max_distance_squared, 
      IsCompatible, compatible_data, face_container->face);
  }

  // Check if node is interior
  if (node->children[0]) {
    assert(node->children[1]);
    assert(node->small_faces.IsEmpty());

    // Compute distance from query point to split plane
    RNScalar side = query_position[node->split_dimension] - node->split_coordinate;

    // Search children nodes
    if (side <= 0) {
      // Search negative side first
      R3Box child_box(node_box);
      child_box[RN_HI][node->split_dimension] = node->split_coordinate;
      FindClosest(query_position, query_normal, closest, 
        min_distance_squared, max_distance_squared, IsCompatible, compatible_data,
        node->children[0], child_box);
      if (side*side < max_distance_squared) {
        R3Box child_box(node_box);
        child_box[RN_LO][node->split_dimension] = node->split_coordinate;
        FindClosest(query_position, query_normal, closest, 
          min_distance_squared, max_distance_squared, IsCompatible, compatible_data,
          node->children[1], child_box);
      }
    }
    else {
      // Search positive side first
      R3Box child_box(node_box);
      child_box[RN_LO][node->split_dimension] = node->split_coordinate;
      FindClosest(query_position, query_normal, closest, 
        min_distance_squared, max_distance_squared, IsCompatible, compatible_data,
        node->children[1], child_box);
      if (side*side < max_distance_squared) {
        R3Box child_box(node_box);
        child_box[RN_HI][node->split_dimension] = node->split_coordinate;
        FindClosest(query_position, query_normal, closest, 
          min_distance_squared, max_distance_squared, IsCompatible, compatible_data,
          node->children[0], child_box);
      }
    }
  }
  else {
    // Update based on distance to each small face
    for (int i = 0; i < node->small_faces.NEntries(); i++) {
      // Get face container and check mark
      R3MeshSearchTreeFace *face_container = node->small_faces[i];
      if (face_container->mark == mark) continue;
      face_container->mark = mark;

      // Find closest point in mesh face
      FindClosest(query_position, query_normal, closest, 
        min_distance_squared, max_distance_squared, 
        IsCompatible, compatible_data, face_container->face);
    }
  }
}



void R3MeshSearchTree::
FindClosest(const R3Point& query_position, const R3Vector& query_normal, R3MeshIntersection& closest,
  RNScalar min_distance, RNScalar max_distance, 
  int (*IsCompatible)(const R3Point&, const R3Vector&, R3Mesh *, R3MeshFace *, void *), void *compatible_data)
{
  // Initialize result
  closest.type = R3_MESH_NULL_TYPE;
  closest.vertex = NULL;
  closest.edge = NULL;
  closest.face = NULL;
  closest.point = R3zero_point;
  closest.t = 0;

  // Check root
  if (!root) return;

  // Update mark (used to avoid checking same face twice)
  mark++;

  // Use squared distances for efficiency
  RNScalar min_distance_squared = min_distance * min_distance;
  RNScalar closest_distance_squared = max_distance * max_distance;

  // Search nodes recursively
  FindClosest(query_position, query_normal, closest, 
    min_distance_squared, closest_distance_squared, 
    IsCompatible, compatible_data, 
    root, BBox());

  // Update result
  closest.t = sqrt(closest_distance_squared);
  if (closest.type == R3_MESH_VERTEX_TYPE) { 
    closest.face = mesh->FaceOnVertex(closest.vertex);
    closest.edge = mesh->EdgeOnVertex(closest.vertex, closest.face, RN_CCW); 
  }
  else if (closest.type == R3_MESH_EDGE_TYPE) { 
    closest.vertex = NULL; 
    closest.face = mesh->FaceOnEdge(closest.edge);
  }
  else if (closest.type == R3_MESH_FACE_TYPE) { 
    closest.vertex = NULL; 
    closest.edge = NULL; 
  }

  // Just checking
  assert((closest.type == R3_MESH_NULL_TYPE) || (RNIsZero(R3SquaredDistance(closest.point, mesh->ClosestPointOnFace(closest.face, closest.point)), RN_SMALL_EPSILON)));
}



void R3MeshSearchTree::
FindClosest(const R3Point& query_position, R3MeshIntersection& closest,
  RNScalar min_distance, RNScalar max_distance,
  int (*IsCompatible)(const R3Point&, const R3Vector&, R3Mesh *, R3MeshFace *, void *), void *compatible_data)
{
  // Find closest point, ignoring normal
  FindClosest(query_position, R3zero_vector, closest, min_distance, max_distance, IsCompatible, compatible_data);
}



////////////////////////////////////////////////////////////////////////
// Find all search functions (up to distance cutoff) 
////////////////////////////////////////////////////////////////////////

void R3MeshSearchTree::
FindAll(const R3Point& query_position, const R3Vector& query_normal, RNArray<R3MeshIntersection *>& hits, 
  RNScalar min_distance_squared, RNScalar max_distance_squared, 
  int (*IsCompatible)(const R3Point&, const R3Vector&, R3Mesh *, R3MeshFace *, void *), void *compatible_data,
  R3MeshFace *face) const
{
  // Check distance to plane
  const R3Plane& plane = mesh->FacePlane(face);
  RNScalar plane_signed_distance = R3SignedDistance(plane, query_position);
  RNScalar plane_distance_squared = plane_signed_distance * plane_signed_distance;
  if (plane_distance_squared >= max_distance_squared) return;

  // Check distance to bounding box
  RNScalar bbox_distance_squared = DistanceSquared(query_position, mesh->FaceBBox(face), max_distance_squared);
  if (bbox_distance_squared >= max_distance_squared) return;

  // Check compatibility 
  if (IsCompatible) {
    if (!(*IsCompatible)(query_position, query_normal, mesh, face, compatible_data)) return;
  }

  // Get face vertices
  R3MeshVertex *v0 = mesh->VertexOnFace(face, 0);
  R3MeshVertex *v1 = mesh->VertexOnFace(face, 1);
  R3MeshVertex *v2 = mesh->VertexOnFace(face, 2);

  // Get vertex positions
  const R3Point& p0 = mesh->VertexPosition(v0);
  const R3Point& p1 = mesh->VertexPosition(v1);
  const R3Point& p2 = mesh->VertexPosition(v2);

  // Project query point onto face plane
  const R3Vector& face_normal = mesh->FaceNormal(face);
  R3Point plane_point = query_position - plane_signed_distance * face_normal;

  // Check sides of edges
  R3Vector e0 = p1 - p0;
  e0.Normalize();
  R3Vector n0 = mesh->FaceNormal(face) % e0;
  R3Plane s0(p0, n0);
  RNScalar b0 = R3SignedDistance(s0, plane_point);
  R3Vector e1 = p2 - p1;
  e1.Normalize();
  R3Vector n1 = mesh->FaceNormal(face) % e1;
  R3Plane s1(p1, n1);
  RNScalar b1 = R3SignedDistance(s1, plane_point);
  R3Vector e2 = p0 - p2;
  e2.Normalize();
  R3Vector n2 = mesh->FaceNormal(face) % e2;
  R3Plane s2(p2, n2);
  RNScalar b2 = R3SignedDistance(s2, plane_point);

  // Initialize hit info
  R3MeshIntersection hit;
  hit.type = R3_MESH_NULL_TYPE;

  // Consider plane_point's position in relation to edges of the triangle
  if ((b0 >= 0) && (b1 >= 0) && (b2 >= 0)) {
    // Point is inside face
    if (plane_distance_squared >= min_distance_squared) {
      hit.type = R3_MESH_FACE_TYPE;
      hit.vertex = NULL;
      hit.edge = NULL;
      hit.face = face;
      hit.point = plane_point;
      hit.t = sqrt(plane_distance_squared);
    }
  }
  else {
    // Point is outside face -- check each edge
    if (b0 < 0) {
      // Outside edge0
      R3Vector edge_vector = p1 - p0;
      RNScalar edge_length = edge_vector.Length();
      if (edge_length > 0) {
        edge_vector /= edge_length;
        R3Vector point_vector = plane_point - p0;
        RNScalar t = edge_vector.Dot(point_vector);
        if (t <= 0) {
          RNScalar distance_squared = DistanceSquared(query_position, p0);
          if ((distance_squared >= min_distance_squared) && (distance_squared < max_distance_squared)) {
            hit.type = R3_MESH_VERTEX_TYPE;
            hit.vertex = v0;
            hit.edge = mesh->EdgeOnVertex(hit.vertex, face);
            hit.face = face;
            hit.point = p0;
            hit.t = sqrt(distance_squared);
            max_distance_squared = distance_squared;
          }
        }
        else if (t >= edge_length) {
          RNScalar distance_squared = DistanceSquared(query_position, p1);
          if ((distance_squared >= min_distance_squared) && (distance_squared < max_distance_squared)) {
            hit.type = R3_MESH_VERTEX_TYPE;
            hit.vertex = v1;
            hit.edge = mesh->EdgeOnVertex(hit.vertex, face);
            hit.face = face;
            hit.point = p1;
            hit.t = sqrt(distance_squared);
            max_distance_squared = distance_squared;
          }
        }
        else {
          R3Point point = p0 + t * edge_vector;
          RNScalar distance_squared = DistanceSquared(query_position, point);
          if ((distance_squared >= min_distance_squared) && (distance_squared < max_distance_squared)) {
            hit.type = R3_MESH_EDGE_TYPE;
            hit.vertex = NULL;
            hit.edge = mesh->EdgeOnFace(face, 0);
            hit.face = face;
            hit.point = point;
            hit.t = sqrt(distance_squared);
            max_distance_squared = distance_squared;
          }
        }
      }
    }
    if (b1 < 0) {
      // Outside edge1
      R3Vector edge_vector = p2 - p1;
      RNScalar edge_length = edge_vector.Length();
      if (edge_length > 0) {
        edge_vector /= edge_length;
        R3Vector point_vector = plane_point - p1;
        RNScalar t = edge_vector.Dot(point_vector);
        if (t <= 0) {
          RNScalar distance_squared = DistanceSquared(query_position, p1);
          if ((distance_squared >= min_distance_squared) && (distance_squared < max_distance_squared)) {
            hit.type = R3_MESH_VERTEX_TYPE;
            hit.vertex = v1;
            hit.edge = mesh->EdgeOnVertex(hit.vertex, face);
            hit.face = face;
            hit.point = p1;
            hit.t = sqrt(distance_squared);
            max_distance_squared = distance_squared;
          }
        }
        else if (t >= edge_length) {
          RNScalar distance_squared = DistanceSquared(query_position, p2);
          if ((distance_squared >= min_distance_squared) && (distance_squared < max_distance_squared)) {
            hit.type = R3_MESH_VERTEX_TYPE;
            hit.vertex = v2;
            hit.edge = mesh->EdgeOnVertex(hit.vertex, face);
            hit.face = face;
            hit.point = p2;
            hit.t = sqrt(distance_squared);
            max_distance_squared = distance_squared;
          }
        }
        else {
          R3Point point = p1 + t * edge_vector;
          RNScalar distance_squared = DistanceSquared(query_position, point);
          if ((distance_squared >= min_distance_squared) && (distance_squared < max_distance_squared)) {
            hit.type = R3_MESH_EDGE_TYPE;
            hit.vertex = NULL;
            hit.edge = mesh->EdgeOnFace(face, 1);
            hit.face = face;
            hit.point = point;
            hit.t = sqrt(distance_squared);
            max_distance_squared = distance_squared;
          }
        }
      }
    }
    if (b2 < 0) {
      // Outside edge2
      R3Vector edge_vector = p0 - p2;
      RNScalar edge_length = edge_vector.Length();
      if (edge_length > 0) {
        edge_vector /= edge_length;
        R3Vector point_vector = plane_point - p2;
        RNScalar t = edge_vector.Dot(point_vector);
        if (t <= 0) {
          RNScalar distance_squared = DistanceSquared(query_position, p2);
          if ((distance_squared >= min_distance_squared) && (distance_squared < max_distance_squared)) {
            hit.type = R3_MESH_VERTEX_TYPE;
            hit.vertex = v2;
            hit.edge = mesh->EdgeOnVertex(hit.vertex, face);
            hit.face = face;
            hit.point = p2;
            hit.t = sqrt(distance_squared);
            max_distance_squared = distance_squared;
          }
        }
        else if (t >= edge_length) {
          RNScalar distance_squared = DistanceSquared(query_position, p0);
          if ((distance_squared >= min_distance_squared) && (distance_squared < max_distance_squared)) {
            hit.type = R3_MESH_VERTEX_TYPE;
            hit.vertex = v0;
            hit.edge = mesh->EdgeOnVertex(hit.vertex, face);
            hit.face = face;
            hit.point = p0;
            hit.t = sqrt(distance_squared);
            max_distance_squared = distance_squared;
          }
        }
        else {
          R3Point point = p2 + t * edge_vector;
          RNScalar distance_squared = DistanceSquared(query_position, point);
          if ((distance_squared >= min_distance_squared) && (distance_squared < max_distance_squared)) {
            hit.type = R3_MESH_EDGE_TYPE;
            hit.vertex = NULL;
            hit.edge = mesh->EdgeOnFace(face, 2);
            hit.face = face;
            hit.point = point;
            hit.t = sqrt(distance_squared);
            max_distance_squared = distance_squared;
          }
        }
      }
    }
  }

  // Insert hit
  if (hit.type != R3_MESH_NULL_TYPE) {
    hits.Insert(new R3MeshIntersection(hit));
  }
}



void R3MeshSearchTree::
FindAll(const R3Point& query_position, const R3Vector& query_normal, RNArray<R3MeshIntersection *>& hits, 
  RNScalar min_distance_squared, RNScalar max_distance_squared, 
  int (*IsCompatible)(const R3Point&, const R3Vector&, R3Mesh *, R3MeshFace *, void *), void *compatible_data,
  R3MeshSearchTreeNode *node, const R3Box& node_box) const
{
  // Compute distance (squared) from query point to node bbox
  RNScalar distance_squared = DistanceSquared(query_position, node_box, max_distance_squared);
  if (distance_squared >= max_distance_squared) return;

  // Check each big face
  for (int i = 0; i < node->big_faces.NEntries(); i++) {
    // Get face container and check mark
    R3MeshSearchTreeFace *face_container = node->big_faces[i];
    if (face_container->mark == mark) continue;
    face_container->mark = mark;
  
    // Find point in mesh face
    FindAll(query_position, query_normal, hits, 
      min_distance_squared, max_distance_squared, 
      IsCompatible, compatible_data, face_container->face);
  }

  // Check if node is interior
  if (node->children[0]) {
    assert(node->children[1]);
    assert(node->small_faces.IsEmpty());

    // Compute distance from query point to split plane
    RNScalar side = query_position[node->split_dimension] - node->split_coordinate;

    // Search children nodes
    if (side <= 0) {
      // Search negative side first
      R3Box child_box(node_box);
      child_box[RN_HI][node->split_dimension] = node->split_coordinate;
      FindAll(query_position, query_normal, hits, 
        min_distance_squared, max_distance_squared, IsCompatible, compatible_data,
        node->children[0], child_box);
      if (side*side < max_distance_squared) {
        R3Box child_box(node_box);
        child_box[RN_LO][node->split_dimension] = node->split_coordinate;
        FindAll(query_position, query_normal, hits, 
          min_distance_squared, max_distance_squared, IsCompatible, compatible_data,
          node->children[1], child_box);
      }
    }
    else {
      // Search positive side first
      R3Box child_box(node_box);
      child_box[RN_LO][node->split_dimension] = node->split_coordinate;
      FindAll(query_position, query_normal, hits, 
        min_distance_squared, max_distance_squared, IsCompatible, compatible_data,
        node->children[1], child_box);
      if (side*side < max_distance_squared) {
        R3Box child_box(node_box);
        child_box[RN_HI][node->split_dimension] = node->split_coordinate;
        FindAll(query_position, query_normal, hits, 
          min_distance_squared, max_distance_squared, IsCompatible, compatible_data,
          node->children[0], child_box);
      }
    }
  }
  else {
    // Check each small face
    for (int i = 0; i < node->small_faces.NEntries(); i++) {
      // Get face container and check mark
      R3MeshSearchTreeFace *face_container = node->small_faces[i];
      if (face_container->mark == mark) continue;
      face_container->mark = mark;

      // Find point in mesh face
      FindAll(query_position, query_normal, hits, 
        min_distance_squared, max_distance_squared, 
        IsCompatible, compatible_data, face_container->face);
    }
  }
}



void R3MeshSearchTree::
FindAll(const R3Point& query_position, const R3Vector& query_normal, RNArray<R3MeshIntersection *>& hits, 
  RNScalar min_distance, RNScalar max_distance, 
  int (*IsCompatible)(const R3Point&, const R3Vector&, R3Mesh *, R3MeshFace *, void *), void *compatible_data)
{
  // Check root
  if (!root) return;

  // Update mark (used to avoid checking same face twice)
  mark++;

  // Use squared distances for efficiency
  RNScalar min_distance_squared = min_distance * min_distance;
  RNScalar max_distance_squared = max_distance * max_distance;

  // Search nodes recursively
  FindAll(query_position, query_normal, hits,
    min_distance_squared, max_distance_squared, 
    IsCompatible, compatible_data, 
    root, BBox());
}



void R3MeshSearchTree::
FindAll(const R3Point& query_position, RNArray<R3MeshIntersection *>& hits, 
  RNScalar min_distance, RNScalar max_distance,
  int (*IsCompatible)(const R3Point&, const R3Vector&, R3Mesh *, R3MeshFace *, void *), void *compatible_data)
{
  // Find closest point, ignoring normal
  FindAll(query_position, R3zero_vector, hits, min_distance, max_distance, IsCompatible, compatible_data);
}



////////////////////////////////////////////////////////////////////////
// Ray intersection search functions
////////////////////////////////////////////////////////////////////////

void R3MeshSearchTree::
FindIntersection(const R3Ray& ray, R3MeshIntersection& closest, 
  RNScalar min_t, RNScalar& max_t, 
  int (*IsCompatible)(const R3Point&, const R3Vector&, R3Mesh *, R3MeshFace *, void *), void *compatible_data,
  R3MeshFace *face) const
{
  // Check compatibility 
  if (IsCompatible) {
    if (!(*IsCompatible)(ray.Start(), ray.Vector(), mesh, face, compatible_data)) return;
  }

  // Check intersection with plane (this is redundant, but allows checking min_t and max_t)
  RNScalar plane_t;
  if (!R3Intersects(ray, mesh->FacePlane(face), NULL, &plane_t)) return;
  if (plane_t >= max_t) return;
  if (plane_t < min_t) return;

  // Check intersection with face
  R3MeshIntersection face_intersection;
  if (!mesh->Intersection(ray, face, &face_intersection)) return;
  if (face_intersection.t >= max_t) return;
  if (face_intersection.t < min_t) return;

  // Update closest intersection
  closest.type = R3_MESH_FACE_TYPE;
  closest.face = face;
  closest.point = face_intersection.point;
  closest.t = face_intersection.t;
  max_t = face_intersection.t;
}



void R3MeshSearchTree::
FindIntersection(const R3Ray& ray, R3MeshIntersection& closest, 
  RNScalar min_t, RNScalar& max_t, 
  int (*IsCompatible)(const R3Point&, const R3Vector&, R3Mesh *, R3MeshFace *, void *), void *compatible_data,
  R3MeshSearchTreeNode *node, const R3Box& node_box) const
{
  // Find intersection with bounding box
  RNScalar node_box_t;
  if (!R3Intersects(ray, node_box, NULL, NULL, &node_box_t)) return;
  if (node_box_t > max_t) return;

  // Update based on closest intersection to each big face
  for (int i = 0; i < node->big_faces.NEntries(); i++) {
    // Get face container and check mark
    R3MeshSearchTreeFace *face_container = node->big_faces[i];
    if (face_container->mark == mark) continue;
    face_container->mark = mark;
  
    // Find closest point in mesh face
    FindIntersection(ray, closest, min_t, max_t, 
      IsCompatible, compatible_data, face_container->face);
  }

  // Check if node is interior
  if (node->children[0]) {
    assert(node->children[1]);
    assert(node->small_faces.IsEmpty());

    // Compute distance from query point to split plane
    RNScalar side = ray.Start()[node->split_dimension] - node->split_coordinate;
    RNScalar vec = ray.Vector()[node->split_dimension];
    RNScalar plane_t = (side*vec < 0) ? -side/vec : RN_INFINITY;
    
    // Search children nodes
    if (side <= 0) {
      // Search negative side first
      if (plane_t >= min_t) {
        R3Box child_box(node_box);
        child_box[RN_HI][node->split_dimension] = node->split_coordinate;
        FindIntersection(ray, closest, min_t, max_t,
          IsCompatible, compatible_data, node->children[0], child_box);
      }
      if (plane_t < max_t) {
        R3Box child_box(node_box);
        child_box[RN_LO][node->split_dimension] = node->split_coordinate;
        FindIntersection(ray, closest, min_t, max_t, 
          IsCompatible, compatible_data, node->children[1], child_box);
      }
    }
    else {
      // Search positive side first
      if (plane_t >= min_t) {
        R3Box child_box(node_box);
        child_box[RN_LO][node->split_dimension] = node->split_coordinate;
        FindIntersection(ray, closest, min_t, max_t, 
          IsCompatible, compatible_data, node->children[1], child_box);
      }
      if (plane_t < max_t) {
        R3Box child_box(node_box);
        child_box[RN_HI][node->split_dimension] = node->split_coordinate;
        FindIntersection(ray, closest, min_t, max_t,
          IsCompatible, compatible_data, node->children[0], child_box);
      }
    }
  }
  else {
    // Update based on distance to each small face
    for (int i = 0; i < node->small_faces.NEntries(); i++) {
      // Get face container and check mark
      R3MeshSearchTreeFace *face_container = node->small_faces[i];
      if (face_container->mark == mark) continue;
      face_container->mark = mark;

      // Find closest point in mesh face
      FindIntersection(ray, closest, min_t, max_t,
        IsCompatible, compatible_data, face_container->face);
    }
  }
}



void R3MeshSearchTree::
FindIntersection(const R3Ray& ray, R3MeshIntersection& closest,
  RNScalar min_t, RNScalar max_t, 
  int (*IsCompatible)(const R3Point&, const R3Vector&, R3Mesh *, R3MeshFace *, void *), void *compatible_data)
{
  // Initialize result
  closest.type = R3_MESH_NULL_TYPE;
  closest.vertex = NULL;
  closest.edge = NULL;
  closest.face = NULL;
  closest.point = R3zero_point;
  closest.t = 0;

  // Check root
  if (!root) return;

  // Update mark (used to avoid checking same face twice)
  mark++;

  // Search nodes recursively
  FindIntersection(ray, closest,
    min_t, max_t,
    IsCompatible, compatible_data, 
    root, BBox());
}



////////////////////////////////////////////////////////////////////////
// Distance functions
////////////////////////////////////////////////////////////////////////

#if 0

static RNLength 
R3Distance(R3MeshSearchTree *mesh1, R3MeshFace *face1, R3MeshSearchTree *tree2, R3MeshFace *face2,
  R3MeshIntersection *closest_point1, R3MeshIntersection *closest_point2)
{
  RNLength closest_distance = FLT_MAX;
  for (int i = 0; i < 3; i++) {
    R3MeshIntersection intersection1;
    R3MeshVertex *vertex1 = mesh1->VertexOnFace(face1, i);
    const R3Point& position1 = mesh1->VertexPosition(vertex1);
    FindClosest(vertex1, R3zero_vector, intersection1, 0, FLT_MAX, NULL, NULL, face2)
  RNScalar min_distance_squared, RNScalar& max_distance_squared, 
  int (*IsCompatible)(const R3Point&, const R3Vector&, R3Mesh *, R3MeshFace *, void *), void *compatible_data,
  R3MeshFace *face) const
}



static RNLength 
R3Distance(const R3MeshSearchTree& tree1, const R3MeshSearchTree& tree2, 
  RNScalar min_distance, RNScalar max_distance,
  int (*IsCompatible)(R3Mesh *, R3MeshFace *, R3Mesh *, R3MeshFace *, void *), void *compatible_data,
  R3MeshIntersection *closest_point1, R3MeshIntersection *closest_point2,
  R3MeshSearchTreeNode *node1, const R3Box& node_box1, R3MeshSearchTreeNode *node2, const R3Box& node_box2)

{
}



RNLength 
R3Distance(const R3MeshSearchTree& tree1, const R3MeshSearchTree& tree2, 
  RNScalar min_distance, RNScalar max_distance,
  int (*IsCompatible)(R3Mesh *, R3MeshFace *, R3Mesh *, R3MeshFace *, void *), void *compatible_data,
  R3MeshIntersection *closest_point1, R3MeshIntersection *closest_point2)
{
}

#endif



////////////////////////////////////////////////////////////////////////
// Visualization and debugging functions
////////////////////////////////////////////////////////////////////////

void R3MeshSearchTree::
Outline(R3MeshSearchTreeNode *node, const R3Box& node_box) const
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



void R3MeshSearchTree::
Outline(void) const
{
  // Draw kdtree nodes recursively
  if (!root) return;
  Outline(root, BBox());
}



int R3MeshSearchTree::
Print(R3MeshSearchTreeNode *node, int depth) const
{
  // Check node
  if (!node) return 0;

  // Initialize number of decendents
  int ndecendents0 = 0;
  int ndecendents1 = 0;

  // Process interior node
  if (node->children[0] && node->children[1]) {
    // Print balance of children
    ndecendents0 = Print(node->children[0], depth+1);
    ndecendents1 = Print(node->children[1], depth+1);

    // Print balance of this node
    printf("%d", depth);
    for (int i = 0; i <= depth; i++) printf("  ");
    printf("I %d : %d %d %g\n", node->big_faces.NEntries(), ndecendents0, ndecendents1, (double) ndecendents0 / (double) ndecendents1);
  }
  else {
    printf("%d", depth);
    for (int i = 0; i <= depth; i++) printf("  ");
    printf("L %d \n", node->small_faces.NEntries());
  }

  // Return number of nodes rooted in this subtree
  return 1 + ndecendents0 + ndecendents1;
}



void R3MeshSearchTree::
Print(void) const
{
  // Print recursively
  Print(root, 0);
}


int R3MeshSearchTree::
NNodes(void) const
{
  // Return number of nodes
  return nnodes;
}



////////////////////////////////////////////////////////////////////////
// Utility functions
////////////////////////////////////////////////////////////////////////

RNScalar R3MeshSearchTree::
DistanceSquared(const R3Point& query_position, const R3Point& point) const
{
  // Compute squared distance from query to point
  RNScalar dx = query_position[0] - point[0];
  RNScalar dy = query_position[1] - point[1];
  RNScalar dz = query_position[2] - point[2];
  return dx*dx + dy*dy + dz*dz;
}



RNScalar R3MeshSearchTree::
DistanceSquared(const R3Point& query_position, const R3Box& box, RNScalar max_distance_squared) const
{
  // Find and check axial distances from face to node box
  RNScalar dx, dy, dz;
  if (query_position.X() > box.XMax()) dx = query_position.X() - box.XMax();
  else if (query_position.X() < box.XMin()) dx = box.XMin()- query_position.X();
  else dx = 0.0;
  RNScalar dx_squared = dx * dx;
  if (dx_squared >= max_distance_squared) return dx_squared;
  if (query_position.Y() > box.YMax()) dy = query_position.Y() - box.YMax();
  else if (query_position.Y() < box.YMin()) dy = box.YMin()- query_position.Y();
  else dy = 0.0;
  RNScalar dy_squared = dy * dy;
  if (dy_squared >= max_distance_squared) return dy_squared;
  if (query_position.Z() > box.ZMax()) dz = query_position.Z() - box.ZMax();
  else if (query_position.Z() < box.ZMin()) dz = box.ZMin()- query_position.Z();
  else dz = 0.0;
  RNScalar dz_squared = dz * dz;
  if (dz_squared >= max_distance_squared) return dz_squared;
    
  // Find and check actual distance from face to node box
  RNScalar distance_squared = 0;
  if ((dy == 0.0) && (dz == 0.0)) distance_squared = dx_squared;
  else if ((dx == 0.0) && (dz == 0.0)) distance_squared = dy_squared;
  else if ((dx == 0.0) && (dy == 0.0)) distance_squared = dz_squared;
  else distance_squared = dx_squared + dy_squared + dz_squared;

  // Return distance squared
  return distance_squared;
}



#if 0

XXX THIS IS VERY INCOMPLETE >>>
XXX THE MARKS CANNOT BE USED LIKE THEY ARE NOW XXX
XXX THE AFFINE TRANFORMATIONS NEED TO BE APPLIED TO THE BOXES AND FACES DURING CHECKS XXX
XXX THE CLOSEST POINTS NEED TO BE INVERSE TRANSFORMED BEFORE RETURNING RESULTS XXX

RNLength
R3Distance(R3Mesh *mesh1, R3MeshFace *face1, 
  R3Mesh *mesh2, R3MeshFace *face2, 
  R3MeshIntersection* closest_point1, 
  R3MeshIntersection* closest_point2)
{
  // XXX This is not right when closest points are in middles of both edges XXX
  // XXX but that's probably OK for now XXX 

  // Initialize result
  RNScalar closest_distance = FLT_MAX;

  // Find closest point on face2 for each vertex of face1
  for (int i1 = 0; i1 < 3; i1++) {
    R3MeshIntersection intersection2;
    R3MeshVertex *vertex1 = mesh1->VertexOnFace(face1, i1);
    const R3Point& position1 = mesh1->VertexPosition(vertex1);
    RNLength distance = R3Distance(mesh2, face2, position1, (closest_point2) ? &intersection2 : NULL);
    if (distance < closest_distance) {
      closest_distance = distance;
      if (closest_point1) {
        closest_point1->type = R3_MESH_VERTEX_TYPE;
        closest_point1->vertex = vertex1;
        closest_point1->edge = NULL;
        closest_point1->face = face1;
        closest_point1->point = position1;
        closest_point1->t = intersection2.t;
      }
      if (closest_point2) {
        *closest_point2 = intersection2;
      }
    }
  }

  // Find closest point on face1 for each vertex of face2
  for (int i2 = 0; i2 < 3; i2++) {
    R3MeshIntersection intersection1;
    R3MeshVertex *vertex2 = mesh2->VertexOnFace(face2, i2);
    const R3Point& position2 = mesh2->VertexPosition(vertex2);
    RNLength distance = R3Distance(mesh1, face1, position2, &intersection1);
    if (distance < closest_distance) {
      closest_distance = distance;
      if (closest_point2) {
        closest_point2->type = R3_MESH_VERTEX_TYPE;
        closest_point2->vertex = vertex2;
        closest_point2->edge = NULL;
        closest_point2->face = face2;
        closest_point2->point = position2;
        closest_point2->t = intersection1.t;
      }
      if (closest_point1) {
        *closest_point1 = intersection1;
      }
    }
  }

  // Return closest distance
  return closest_distance;
}



static void 
FindClosest(const R3MeshSearchTree& tree1, const R3MeshSearchTree& tree2, const R3Affine& xform2, 
  R3MeshIntersection& closest1, R3MeshIntersection& closest2, 
  RNScalar min_distance, RNScalar& max_distance, 
  int (*IsCompatible)(R3Mesh *, R3MeshFace *, R3Mesh *, R3MeshFace *, void *), void *compatible_data,
  R3MeshFace *face1, R3MeshFace *face2) const
{
}



static void 
FindClosest(const R3MeshSearchTree& tree1, const R3MeshSearchTree& tree2, const R3Affine& xform2, 
  R3MeshIntersection& closest1, R3MeshIntersection& closest2, 
  RNScalar min_distance, RNScalar& max_distance, 
  int (*IsCompatible)(R3Mesh *, R3MeshFace *, R3Mesh *, R3MeshFace *, void *), void *compatible_data,
  R3MeshSearchTreeNode *node1, const R3Box& box1, R3MeshFace *face2) const
{
}



static void 
FindClosest(const R3MeshSearchTree& tree1, const R3MeshSearchTree& tree2, 
  const R3Affine& xform12, const R3Affine& xform21, R3MeshIntersection& closest, 
  RNScalar min_distance, RNScalar& max_distance, 
  int (*IsCompatible)(R3Mesh *, R3MeshFace *, R3Mesh *, R3MeshFace *, void *), void *compatible_data,
  R3MeshSearchTreeNode *node1, const R3Box& box1, R3MeshSearchTreeNode *node2, const R3Box& box2) const
{
  // Compute distance between boxes
  RNScalar distance = R3Distance(box1, box2);
  if (distance > max_distance) return;
  if (distance < min_distance) return;

  // Update based on distance to each big face in node1
  for (int i = 0; i < node1->big_faces.NEntries(); i++) {
    // Get face container and check mark
    R3MeshSearchTreeFace *face_container = node1->big_faces[i];
    if (face_container->mark == mark) continue;
    face_container->mark = mark;
  
    // Find closest point for mesh face
    FindClosest(tree2, tree1, affine12, closest2, closest1, 
      min_distance, max_distance, 
      IsCompatible, compatible_data, 
      node2, box2, face_container->face);
  }

  // Update based on distance to each big face in node2
  for (int i = 0; i < node2->big_faces.NEntries(); i++) {
    // Get face container and check mark
    R3MeshSearchTreeFace *face_container = node2->big_faces[i];
    if (face_container->mark == mark) continue;
    face_container->mark = mark;
  
    // Find closest point for mesh face
    FindClosest(tree1, tree2, affine21, closest1, closest2, 
      min_distance, max_distance, 
      IsCompatible, compatible_data, 
      node1, box1, face_container->face);
  }

  // Check if node1 and node2 are interior
  if (node1->children[0]) {
    // node1 is interior
    assert(node1->children[1]);
    assert(node1->small_faces.IsEmpty());

    // Check if should visit child1 on negative side of split plane
    if (node1->split_coordinate > box2[0][node1->split_dimension] - max_distance) {
      if (node1->split_coordinate < box2[1][node1->split_dimension] + max_distance) {
        R3Box child_box1(box1);
        child_box1[RN_HI][node1->split_dimension] = node1->split_coordinate;
        FindClosest(tree1, tree2, xform12, xform21, closest1, closest2,
          min_distance, max_distance, IsCompatible, compatible_data,
          node1->children[0], child_box1, node2, box2);
      }
    }

    // Check if should visit child1 on positive side of split plane
    if (node1->split_coordinate < box2[1][node1->split_dimension] + max_distance) {
      if (node1->split_coordinate > box2[0][node1->split_dimension] - max_distance) {
        R3Box child_box1(box1);
        child_box1[RN_LO][node1->split_dimension] = node1->split_coordinate;
        FindClosest(tree1, tree2, xform12, xform21, closest1, closest2,
          min_distance, max_distance, IsCompatible, compatible_data,
          node1->children[1], child_box1, node2, box2);
      }
    }
  }
  else if (node2->children[0]) {
    // node2 is interior
    assert(node2->children[1]);
    assert(node2->small_faces.IsEmpty());

    // Check if should visit child2 on negative side of split plane
    if (node2->split_coordinate > box1[0][node2->split_dimension] - max_distance) {
      if (node2->split_coordinate < box1[1][node2->split_dimension] + max_distance) {
        R3Box child_box2(box2);
        child_box2[RN_HI][node2->split_dimension] = node2->split_coordinate;
        FindClosest(tree1, tree2, xform12, xform21, closest1, closest2,
          min_distance, max_distance, IsCompatible, compatible_data,
          node1, box1, node2->children[0], child_box2);
      }
    }

    // Check if should visit child2 on positive side of split plane
    if (node2->split_coordinate < box1[1][node2->split_dimension] + max_distance) {
      if (node2->split_coordinate > box1[0][node2->split_dimension] - max_distance) {
        R3Box child_box2(box2);
        child_box2[RN_LO][node2->split_dimension] = node2->split_coordinate;
        FindClosest(tree1, tree2, xform12, xform21, closest1, closest2,
          min_distance, max_distance, IsCompatible, compatible_data,
          node1, box1, node2->children[1], child_box2);
      }
    }
  }
  else {
    // Both node1 and node2 are leaves
    for (int i = 0; i < node2->small_faces.NEntries(); i++) {
      // Get face container and check mark
      R3MeshSearchTreeFace *face_container = node2->small_faces[i];
      if (face_container->mark == mark) continue;
      face_container->mark = mark;

      // Find closest point in node1 to mesh face2
      FindClosest(tree1, tree2, affine21, closest1, closest2, 
        min_distance, max_distance, 
        IsCompatible, compatible_data, 
        node1, box1, face_container->face);
    }
  }
}



static void 
FindClosest(const R3MeshSearchTree& tree1, const R3MeshSearchTree& tree2, 
  R3MeshIntersection& closest1, R3MeshIntersection& closest2,
  RNScalar min_distance, RNScalar max_distance, 
  int (*IsCompatible)(const R3Point&, const R3Vector&, R3Mesh *, R3MeshFace *, void *), void *compatible_data)
{
  // Initialize result
  closest1.type = closest2.type = R3_MESH_NULL_TYPE;
  closest1.vertex = closest2.vertex = NULL;
  closest1.edge = closest2.edge = NULL;
  closest1.face = closest2.face = NULL;
  closest1.point = closest2.point = R3zero_point;
  closest1.t = closest2.t = 0;

  // Check root
  if (!tree1.root || !tree2.root) return;

  // Update mark (used to avoid checking same face twice)
  mark++;

  // Initialize closest distance
  RNLength closest_distance = max_distance;

  // Search nodes recursively
  FindClosest(tree, affine, closest1, closest2,
    min_distance, closest_distance, 
    IsCompatible, compatible_data, 
    tree1.root, tree1.BBox(), tree2.root, tree2.BBox());

  // Update closest1 result
  closest1.t = closest_distance;
  if (closest1.type == R3_MESH_VERTEX_TYPE) { 
    closest1.edge = mesh->EdgeOnVertex(closest.vertex); 
    closest1.face = mesh->FaceOnEdge(closest.edge);
  }
  else if (closest1.type == R3_MESH_EDGE_TYPE) { 
    closest1.vertex = NULL; 
    closest1.face = mesh->FaceOnEdge(closest.edge);
  }
  else if (closest1.type == R3_MESH_FACE_TYPE) { 
    closest1.vertex = NULL; 
    closest1.edge = NULL; 
  }

  // Update closest2 result
  closest2.t = closest_distance;
  if (closest2.type == R3_MESH_VERTEX_TYPE) { 
    closest2.edge = mesh->EdgeOnVertex(closest.vertex); 
    closest2.face = mesh->FaceOnEdge(closest.edge);
  }
  else if (closest2.type == R3_MESH_EDGE_TYPE) { 
    closest2.vertex = NULL; 
    closest2.face = mesh->FaceOnEdge(closest.edge);
  }
  else if (closest2.type == R3_MESH_FACE_TYPE) { 
    closest2.vertex = NULL; 
    closest2.edge = NULL; 
  }
}



#endif



