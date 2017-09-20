/* Source file for the R3 scene node class */



////////////////////////////////////////////////////////////////////////
// INCLUDE FILES
////////////////////////////////////////////////////////////////////////

#include "R3Graphics.h"





////////////////////////////////////////////////////////////////////////
// PKG INITIALIZATION FUNCTIONS
////////////////////////////////////////////////////////////////////////

int 
R3InitSceneNode()
{
  /* Return success */
  return TRUE;
}



void 
R3StopSceneNode()
{
}



////////////////////////////////////////////////////////////////////////
// MEMBER FUNCTIONS
////////////////////////////////////////////////////////////////////////

R3SceneNode::
R3SceneNode(R3Scene *scene) 
  : scene(NULL),
    scene_index(-1),
    parent(NULL),
    parent_index(-1),
    children(),
    elements(),
    references(),
    info(),
    transformation(R3identity_affine),
    bbox(FLT_MAX, FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX),
    name(NULL),
    data(NULL)
{
  // Insert node into scene
  if (scene) {
    scene->InsertNode(this);
  }
}



R3SceneNode::
R3SceneNode(const R3SceneNode& node) 
  : scene(NULL),
    scene_index(-1),
    parent(NULL),
    parent_index(-1),
    children(),
    elements(),
    references(),
    info(node.info),
    transformation(node.transformation),
    bbox(FLT_MAX, FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX),
    name((node.name) ? strdup(node.name) : NULL),
    data(NULL)
{
}



R3SceneNode::
~R3SceneNode(void)
{
  // Remove elements
  for (int i = 0; i < elements.NEntries(); i++) {
    R3SceneElement *element = elements.Kth(i);
    RemoveElement(element);
  }

  // Remove references
  for (int i = 0; i < references.NEntries(); i++) {
    R3SceneReference *reference = references.Kth(i);
    RemoveReference(reference);
  }

  // Remove children
  for (int i = 0; i < children.NEntries(); i++) {
    R3SceneNode *child = children.Kth(i);
    RemoveChild(child);
  }

  // Remove node from parent
  if (parent) {
    parent->RemoveChild(this);
  }

  // Remove node from scene
  if (scene) {
    scene->RemoveNode(this);
  }

  // Delete name
  if (name) free(name);
}



const RNInterval R3SceneNode::
NFacets(void) const
{
  // Initialize nfacets
  RNInterval nfacets(0,0);
  
  // Add facets from elements
  for (int i = 0; i < NElements(); i++) {
    R3SceneElement *element = Element(i);
    nfacets += element->NFacets();
  }
  
  // Add facets from references
  for (int i = 0; i < NReferences(); i++) {
    R3SceneReference *reference = Reference(i);
    R3Scene *referenced_scene = reference->ReferencedScene();
    if (referenced_scene) nfacets += referenced_scene->NFacets();
  }
  
  // Add facets from children
  for (int i = 0; i < NChildren(); i++) {
    R3SceneNode *child = Child(i);
    nfacets += child->NFacets();
  }

  // Return nfacets
  return nfacets;
}



const RNLength R3SceneNode::
Length(void) const
{
  // Initialize length
  RNLength length = 0;
  
  // Add length from elements
  for (int i = 0; i < NElements(); i++) {
    R3SceneElement *element = Element(i);
    length += element->Length();
  }
  
  // Add length from references
  for (int i = 0; i < NReferences(); i++) {
    R3SceneReference *reference = Reference(i);
    R3Scene *referenced_scene = reference->ReferencedScene();
    if (referenced_scene) length += referenced_scene->Length();
  }
  
  // Add length from children
  for (int i = 0; i < NChildren(); i++) {
    R3SceneNode *child = Child(i);
    length += child->Length();
  }

  // Return length
  return length;
}



const RNArea R3SceneNode::
Area(void) const
{
  // Initialize area
  RNArea area = 0;
  
  // Add area from elements
  for (int i = 0; i < NElements(); i++) {
    R3SceneElement *element = Element(i);
    area += element->Area();
  }
  
  // Add area from references
  for (int i = 0; i < NReferences(); i++) {
    R3SceneReference *reference = Reference(i);
    R3Scene *referenced_scene = reference->ReferencedScene();
    if (referenced_scene) area += referenced_scene->Area();
  }
  
  // Add area from children
  for (int i = 0; i < NChildren(); i++) {
    R3SceneNode *child = Child(i);
    area += child->Area();
  }

  // Return area
  return area;
}



const RNVolume R3SceneNode::
Volume(void) const
{
  // Initialize volume
  RNVolume volume = 0;
  
  // Add volume from elements
  for (int i = 0; i < NElements(); i++) {
    R3SceneElement *element = Element(i);
    volume += element->Volume();
  }
  
  // Add volume from references
  for (int i = 0; i < NReferences(); i++) {
    R3SceneReference *reference = Reference(i);
    R3Scene *referenced_scene = reference->ReferencedScene();
    if (referenced_scene) volume += referenced_scene->Volume();
  }
  
  // Add volume from children
  for (int i = 0; i < NChildren(); i++) {
    R3SceneNode *child = Child(i);
    volume += child->Volume();
  }

  // Return volume
  return volume;
}



R3Point R3SceneNode::
ClosestPoint(const R3Point& point) const
{
  // Return closest point
  R3Point result;
  if (!FindClosest(point, NULL, NULL, NULL, &result)) return R3unknown_point;
  else return result;
}



const R3Shape& R3SceneNode::
BShape(void) const
{
  // Return bounding box
  if (bbox[0][0] == FLT_MAX) 
    ((R3SceneNode *) this)->UpdateBBox();
  return bbox;
}



const R3Box& R3SceneNode::
BBox(void) const
{
  // Return bounding box
  if (bbox[0][0] == FLT_MAX) 
    ((R3SceneNode *) this)->UpdateBBox();
  return bbox;
}



R3Box R3SceneNode::
WorldBBox(void) const
{
  // This returns bounding box in world coordinate frame
  // Note that BBox() returns bounding box after xform (in parent's coordinate frame)
  // Note that this bbox will be over-estimating, since it is a transformed AABB

  // Start for bbox in parent's coordinate frame
  R3Box box = BBox();

  // Transform box
  R3SceneNode *ancestor = Parent();
  while (ancestor) {
    box.Transform(ancestor->Transformation());
    ancestor = ancestor->Parent();
  }

  // Return box;
  return bbox;
}



R3Box R3SceneNode::
LocalBBox(void) const
{
  // This returns bounding box in this node's coordinate frame (before transformation)
  // Note that BBox() returns bounding box after xform (in parent's coordinate frame)
  
  // Compute bounding box in local coordinates
  R3Box bbox = R3null_box;

  // Include elements
  for (int i = 0; i < NElements(); i++) {
    R3SceneElement *element = Element(i);
    bbox.Union(element->BBox());
  }

  // Include references
  for (int i = 0; i < NReferences(); i++) {
    R3SceneReference *reference = Reference(i);
    bbox.Union(reference->BBox());
  }

  // Include elements
  for (int i = 0; i < NChildren(); i++) {
    R3SceneNode *child = Child(i);
    bbox.Union(child->BBox());
  }

  // Return box
  return bbox;
}



const R3Point R3SceneNode::
Centroid(void) const
{
  // Return centroid of bounding box
  return BBox().Centroid();
}



RNBoolean R3SceneNode::
IsDecendent(const R3SceneNode *node) const
{
  // Return whether this is an decendent of node
  const R3SceneNode *ancestor = this;
  while (ancestor) {
    if (ancestor == node) return TRUE;
    ancestor = ancestor->Parent();
  }

  // Node not found
  return FALSE;
}



R3Affine R3SceneNode::
CumulativeTransformation(void) const
{
  // Initialize cumulative transformation
  R3Affine cumulative_transformation(R3identity_affine);

  // Concatenate transformations back to the root
  const R3SceneNode *ancestor = this;
  while (ancestor) {
    R3Affine t(R3identity_affine);
    t.Transform(ancestor->Transformation());
    t.Transform(cumulative_transformation);
    cumulative_transformation = t;
    ancestor = ancestor->Parent();
  }

  // Return cumulative transformation from root to node (including node's transformation)
  return cumulative_transformation;
}



R3Affine R3SceneNode::
CumulativeParentTransformation(void) const
{
  // Initialize cumulative transformation
  R3Affine cumulative_transformation(R3identity_affine);

  // Concatenate transformations back to the root
  const R3SceneNode *ancestor = Parent();
  while (ancestor) {
    R3Affine t(R3identity_affine);
    t.Transform(ancestor->Transformation());
    t.Transform(cumulative_transformation);
    cumulative_transformation = t;
    ancestor = ancestor->Parent();
  }

  // Return cumulative transformation from root to parent (not including node's transformation)
  return cumulative_transformation;
}



void R3SceneNode::
InsertChild(R3SceneNode *node) 
{
  // Check scene
  assert(node->scene == scene);
  assert(node->scene_index >= 0);

  // Insert node into children
  assert(!node->parent);
  assert(node->parent_index == -1);
  node->parent = this;
  node->parent_index = children.NEntries();
  children.Insert(node);

  // Invalidate bounding box
  InvalidateBBox();
}



void R3SceneNode::
RemoveChild(R3SceneNode *node) 
{
  // Remove node from children
  assert(node->parent == this);
  RNArrayEntry *entry = children.KthEntry(node->parent_index);
  assert(children.EntryContents(entry) == node);
  R3SceneNode *tail = children.Tail();
  children.EntryContents(entry) = tail;
  tail->parent_index = node->parent_index;
  children.RemoveTail();
  node->parent = NULL;
  node->parent_index = -1;

  // Invalidate bounding box
  InvalidateBBox();
}



void R3SceneNode::
InsertElement(R3SceneElement *element) 
{
  // Update goemetry
  assert(!element->node);
  element->node = this;

  // Insert element
  elements.Insert(element);

  // Invalidate bounding box
  InvalidateBBox();
}



void R3SceneNode::
RemoveElement(R3SceneElement *element) 
{
  // Update goemetry
  assert(element->node == this);
  element->node = NULL;

  // Remove element
  elements.Remove(element);

  // Invalidate bounding box
  InvalidateBBox();
}



void R3SceneNode::
InsertReference(R3SceneReference *reference) 
{
  // Insert reference
  references.Insert(reference);

  // Invalidate bounding box
  InvalidateBBox();
}



void R3SceneNode::
RemoveReference(R3SceneReference *reference) 
{
  // Remove reference
  references.Remove(reference);

  // Invalidate bounding box
  InvalidateBBox();
}



void R3SceneNode::
InsertInfo(const char *key, const char *value) 
{
  // Insert key-value pair
  info.Insert(key, value);
}



void R3SceneNode::
ReplaceInfo(const char *key, const char *value) 
{
  // Replace key-value pair
  info.Replace(key, value);
}



void R3SceneNode::
RemoveInfo(const char *key) 
{
  // Insert key-value pair
  info.Remove(key);
}



void R3SceneNode::
SetTransformation(const R3Affine& transformation)
{
  // Set transformation
  this->transformation = transformation;

  // Invalidate bounding box
  InvalidateBBox();
}



void R3SceneNode::
Transform(const R3Affine& transformation)
{
  // Set transformation
  this->transformation.Transform(transformation);

  // Invalidate bounding box
  InvalidateBBox();
}



void R3SceneNode::
SetName(const char *name)
{
  // Set name
  if (this->name) free(this->name);
  if (name) this->name = strdup(name);
  else this->name = NULL;
}



void R3SceneNode::
SetData(void *data)
{
  // Set data
  this->data = data;
}



RNLength R3SceneNode::
Distance(const R3Point& point) const
{
  // Return distance from point to closest point in node subtree
  RNLength distance = RN_INFINITY;
  if (!FindClosest(point, NULL, NULL, NULL, NULL, NULL, &distance)) return RN_INFINITY;
  return distance;
}



RNBoolean R3SceneNode::
FindClosest(const R3Point& point,
  R3SceneNode **hit_node, R3Material **hit_material, R3Shape **hit_shape,
  R3Point *hit_point, R3Vector *hit_normal, RNLength *hit_d,
  RNLength min_d, RNLength max_d) const
{
  // Temporary variables
  RNBoolean found = FALSE;
  RNScalar d;

  // Check if bounding box is within max_d
  RNScalar bbox_d = R3Distance(point, BBox());
  if (bbox_d > max_d) return FALSE;

  // Apply inverse transformation to point
  R3Point node_point = point;
  node_point.InverseTransform(transformation);

  // Apply inverse transform to min_d and max_d
  RNScalar scale = transformation.ScaleFactor();
  if (RNIsZero(scale)) return FALSE;
  if (min_d < RN_INFINITY) min_d /= scale;
  if (max_d < RN_INFINITY) max_d /= scale;

  // Find closest element 
  for (int i = 0; i < elements.NEntries(); i++) {
    R3SceneElement *element = elements.Kth(i);
    if (element->FindClosest(node_point, hit_shape, hit_point, hit_normal, &d, min_d, max_d)) {
      if ((d >= min_d) && (d <= max_d)) {
        if (hit_node) *hit_node = (R3SceneNode *) this;
        if (hit_material) *hit_material = element->Material();
        if (hit_d) *hit_d = d;
        found = TRUE;
        max_d = d;
      }
    }
  }

  // Find closest reference 
  for (int i = 0; i < references.NEntries(); i++) {
    R3SceneReference *reference = references.Kth(i);
    R3Scene *referenced_scene = reference->ReferencedScene();
    if (referenced_scene) {
      if (referenced_scene->FindClosest(node_point, hit_node, hit_material, hit_shape, hit_point, hit_normal, &d, min_d, max_d)) {
        if ((d >= min_d) && (d <= max_d)) {
          if (hit_d) *hit_d = d;
          found = TRUE;
          max_d = d;
        }
      }
    }
  }

  // Find closest node
  for (int i = 0; i < children.NEntries(); i++) {
    R3SceneNode *child = children.Kth(i);
    if (child->FindClosest(node_point, hit_node, hit_material, hit_shape, hit_point, hit_normal, &d, min_d, max_d)) {
      if ((d >= min_d) && (d <= max_d)) {
        if (hit_d) *hit_d = d;
        found = TRUE;
        max_d = d;
      }
    }
  }

  // Transform hit into parent's coordinate system
  if (found) {
    if (hit_point) hit_point->Transform(transformation);
    if (hit_normal) hit_normal->Transform(transformation);
    if (hit_normal) hit_normal->Normalize();
    if (hit_d) *hit_d *= scale;
  }

  // Return if found a point
  return found;
}



RNBoolean R3SceneNode::
Intersects(const R3Ray& ray,
  R3SceneNode **hit_node, R3Material **hit_material, R3Shape **hit_shape,
  R3Point *hit_point, R3Vector *hit_normal, RNScalar *hit_t,
  RNScalar min_t, RNScalar max_t) const
{
  // Temporary variables
  R3SceneNode *closest_node = NULL;
  R3Point closest_point = R3zero_point;
  RNScalar closest_t = max_t;
  RNScalar bbox_t;
  R3SceneNode *node;
  R3Material *material;
  R3Shape *shape;
  R3Point point;
  R3Vector normal;
  RNScalar t;

  // Check if ray intersects bounding box
  if (!R3Contains(BBox(), ray.Start())) {
    if (!R3Intersects(ray, BBox(), NULL, NULL, &bbox_t)) return FALSE;
    if (RNIsGreater(bbox_t, max_t)) return FALSE;
  }

  // Apply inverse transformation to ray
  R3Ray node_ray = ray;
  node_ray.InverseTransform(transformation);

  // Apply inverse transform to min_t and closest_t
  RNScalar scale = 1.0;
  R3Vector v(ray.Vector());
  v.InverseTransform(transformation);
  RNScalar length = v.Length();
  if (RNIsNegativeOrZero(length)) return FALSE;
  if (RNIsNotEqual(length, 1.0)) {
    scale = length;
    min_t *= scale;
    closest_t *= scale;
  }

  // Find closest element intersection
  for (int i = 0; i < elements.NEntries(); i++) {
    R3SceneElement *element = elements.Kth(i);
    if (element->Intersects(node_ray, &shape, &point, &normal, &t, min_t, closest_t)) {
      if ((t >= min_t) && (t <= closest_t)) {
        if (hit_node) *hit_node = (R3SceneNode *) this;
        if (hit_material) *hit_material = element->Material();
        if (hit_shape) *hit_shape = shape; 
        if (hit_normal) *hit_normal = normal; 
        closest_node = (R3SceneNode *) this;
        closest_point = point;
        closest_t = t;
      }
    }
  }

  // Find closest reference intersection
  for (int i = 0; i < references.NEntries(); i++) {
    R3SceneReference *reference = references.Kth(i);
    R3Scene *referenced_scene = reference->ReferencedScene();
    if (referenced_scene) {
      if (referenced_scene->Intersects(node_ray, &node, &material, &shape, &point, &normal, &t, min_t, closest_t)) {
        if ((t >= min_t) && (t <= closest_t)) {
          if (hit_node) *hit_node = node;
          if (hit_material) *hit_material = material;
          if (hit_shape) *hit_shape = shape; 
          if (hit_normal) *hit_normal = normal; 
          closest_node = (R3SceneNode *) this;
          closest_point = point;
          closest_t = t;
        }
      }
    }
  }

  // Find closest node intersection
  for (int i = 0; i < children.NEntries(); i++) {
    R3SceneNode *child = children.Kth(i);
    if (child->Intersects(node_ray, &node, &material, &shape, &point, &normal, &t, min_t, closest_t)) {
      if ((t >= min_t) && (t <= closest_t)) {
        if (hit_node) *hit_node = node;
        if (hit_material) *hit_material = material;
        if (hit_shape) *hit_shape = shape; 
        if (hit_normal) *hit_normal = normal; 
        closest_node = node;
        closest_point = point;
        closest_t = t;
      }
    }
  }

  // Check if found hit
  if (!closest_node) return FALSE;

  // Transform hit point into parent's coordinate system
  if (hit_t || hit_point) {
    closest_point.Transform(transformation); 
    if (hit_point) *hit_point = closest_point;
    if (hit_t) *hit_t = closest_t / scale; 
  }

  // Transform hit normal into parent's coordinate system
  if (hit_normal) {
    hit_normal->Transform(transformation); 
    hit_normal->Normalize();
  }

  // Return success
  return TRUE;
}



void R3SceneNode::
Draw(const R3DrawFlags draw_flags, const RNArray<R3Material *> *materials) const
{
  // Push transformation
  transformation.Push();

  // Draw elements
  for (int i = 0; i < elements.NEntries(); i++) {
    R3SceneElement *element = elements.Kth(i);
    element->Draw(draw_flags, materials);
  }

  // Draw references
  for (int i = 0; i < references.NEntries(); i++) {
    R3SceneReference *reference = references.Kth(i);
    reference->Draw(draw_flags, materials);
  }

  // Draw children
  for (int i = 0; i < children.NEntries(); i++) {
    R3SceneNode *child = children.Kth(i);
    child->Draw(draw_flags, materials);
  }

  // Pop transformation
  transformation.Pop();
}



void R3SceneNode::
UpdateBBox(void)
{
  // Initialize bounding box
  bbox = R3null_box;

  // Add bounding box of elements
  for (int i = 0; i < elements.NEntries(); i++) {
    R3SceneElement *element = elements.Kth(i);
    R3Box element_bbox = element->BBox();
    element_bbox.Transform(transformation);
    bbox.Union(element_bbox);
  }
    
  // Add bounding box of references
  for (int i = 0; i < references.NEntries(); i++) {
    R3SceneReference *reference = references.Kth(i);
    R3Scene *referenced_scene = reference->ReferencedScene();
    if (!referenced_scene) continue;
    R3Box reference_bbox = referenced_scene->BBox();
    reference_bbox.Transform(transformation);
    bbox.Union(reference_bbox);
  }
    
  // Add bounding box of children nodes
  for (int i = 0; i < children.NEntries(); i++) {
    R3SceneNode *child = children.Kth(i);
    R3Box child_bbox = child->BBox();
    child_bbox.Transform(transformation);
    bbox.Union(child_bbox);
  }
}



void R3SceneNode::
InvalidateBBox(void)
{
  // Invalidate bounding box
  bbox[0][0] = FLT_MAX;

  // Invalidate parent's bounding box
  if (parent) parent->InvalidateBBox();
}



