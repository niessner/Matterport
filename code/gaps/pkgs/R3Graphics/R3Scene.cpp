/* Source file for the R3 scene class */



/* Include files */

#include "R3Graphics.h"



/* Member functions */

int 
R3InitScene()
{
  /* Return success */
  return TRUE;
}



void 
R3StopScene()
{
}



R3Scene::
R3Scene(void)
  : nodes(),
    lights(),
    materials(),
    brdfs(),
    textures(),
    referenced_scenes(),
    info(),
    viewer(),
    ambient(0, 0, 0),
    background(0, 0, 0),
    filename(NULL),
    name(NULL),
    data(NULL)
{
  // Create root node
  root = new R3SceneNode(this);

  // Initialize viewer
  viewer.SetCamera(R3default_camera);
  viewer.SetViewport(R2default_viewport);
}



R3Scene::
~R3Scene(void)
{
  // Delete everything
  // ???

  // Delete filename
  if (filename) free(filename);

  // Delete name
  if (name) free(name);
}



R3SceneNode *R3Scene::
Node(const char *name) const
{
  // Search nodes for matching name
  for (int i = 0; i < NNodes(); i++) {
    R3SceneNode *node = Node(i);
    if (!node->Name()) continue;
    if (!strcmp(node->Name(), name)) return node;
  }

  // Node not found
  return NULL;
}



R3Material *R3Scene::
Material(const char *name) const
{
  // Search materials for matching name
  for (int i = 0; i < NMaterials(); i++) {
    R3Material *material = Material(i);
    if (!material->Name()) continue;
    if (!strcmp(material->Name(), name)) return material;
  }

  // Material not found
  return NULL;
}



R3Light *R3Scene::
Light(const char *name) const
{
  // Search lights for matching name
  for (int i = 0; i < NLights(); i++) {
    R3Light *light = Light(i);
    if (!light->Name()) continue;
    if (!strcmp(light->Name(), name)) return light;
  }

  // Light not found
  return NULL;
}



R3Brdf *R3Scene::
Brdf(const char *name) const
{
  // Search brdfs for matching name
  for (int i = 0; i < NBrdfs(); i++) {
    R3Brdf *brdf = Brdf(i);
    if (!brdf->Name()) continue;
    if (!strcmp(brdf->Name(), name)) return brdf;
  }

  // Brdf not found
  return NULL;
}



R2Texture *R3Scene::
Texture(const char *name) const
{
  // Search textures for matching name
  for (int i = 0; i < NTextures(); i++) {
    R2Texture *texture = Texture(i);
    if (!texture->Name()) continue;
    if (!strcmp(texture->Name(), name)) return texture;
  }

  // Texture not found
  return NULL;
}



R3Scene *R3Scene::
ReferencedScene(const char *name) const
{
  // Search textures for matching name
  for (int i = 0; i < NReferencedScenes(); i++) {
    R3Scene *referenced_scene = ReferencedScene(i);
    if (!referenced_scene->Name()) continue;
    if (!strcmp(referenced_scene->Name(), name)) return referenced_scene;
  }

  // Referenced scene not found
  return NULL;
}



void R3Scene::
InsertNode(R3SceneNode *node) 
{
  // Insert node
  assert(!node->scene);
  assert(node->scene_index == -1);
  node->scene = this;
  node->scene_index = nodes.NEntries();
  nodes.Insert(node);
}



void R3Scene::
RemoveNode(R3SceneNode *node) 
{
  // Remove node
  assert(node->scene == this);
  assert(node->scene_index >= 0);
  RNArrayEntry *entry = nodes.KthEntry(node->scene_index);
  assert(nodes.EntryContents(entry) == node);
  R3SceneNode *tail = nodes.Tail();
  nodes.EntryContents(entry) = tail;
  tail->scene_index = node->scene_index;
  nodes.RemoveTail();
  node->scene_index = -1;
  node->scene = NULL;
}



void R3Scene::
InsertLight(R3Light *light) 
{
  // Insert light
  assert(!light->scene);
  assert(light->scene_index == -1);
  light->scene = this;
  light->scene_index = lights.NEntries();
  lights.Insert(light);
}



void R3Scene::
RemoveLight(R3Light *light) 
{
  // Remove light
  assert(light->scene == this);
  assert(light->scene_index >= 0);
  RNArrayEntry *entry = lights.KthEntry(light->scene_index);
  assert(lights.EntryContents(entry) == light);
  R3Light *tail = lights.Tail();
  lights.EntryContents(entry) = tail;
  tail->scene_index = light->scene_index;
  lights.RemoveTail();
  light->scene_index = -1;
  light->scene = NULL;
}



void R3Scene::
InsertMaterial(R3Material *material) 
{
  // Insert material
  assert(!material->scene);
  assert(material->scene_index == -1);
  material->scene = this;
  material->scene_index = materials.NEntries();
  materials.Insert(material);
}



void R3Scene::
RemoveMaterial(R3Material *material) 
{
  // Remove material
  assert(material->scene == this);
  assert(material->scene_index >= 0);
  RNArrayEntry *entry = materials.KthEntry(material->scene_index);
  assert(materials.EntryContents(entry) == material);
  R3Material *tail = materials.Tail();
  materials.EntryContents(entry) = tail;
  tail->scene_index = material->scene_index;
  materials.RemoveTail();
  material->scene_index = -1;
  material->scene = NULL;
}



void R3Scene::
InsertBrdf(R3Brdf *brdf) 
{
  // Insert brdf
  assert(!brdf->scene);
  assert(brdf->scene_index == -1);
  brdf->scene = this;
  brdf->scene_index = brdfs.NEntries();
  brdfs.Insert(brdf);
}



void R3Scene::
RemoveBrdf(R3Brdf *brdf) 
{
  // Remove brdf
  assert(brdf->scene == this);
  assert(brdf->scene_index >= 0);
  RNArrayEntry *entry = brdfs.KthEntry(brdf->scene_index);
  assert(brdfs.EntryContents(entry) == brdf);
  R3Brdf *tail = brdfs.Tail();
  brdfs.EntryContents(entry) = tail;
  tail->scene_index = brdf->scene_index;
  brdfs.RemoveTail();
  brdf->scene_index = -1;
  brdf->scene = NULL;
}



void R3Scene::
InsertTexture(R2Texture *texture) 
{
  // Insert texture
  assert(!texture->scene);
  assert(texture->scene_index == -1);
  texture->scene = this;
  texture->scene_index = textures.NEntries();
  textures.Insert(texture);
}



void R3Scene::
RemoveTexture(R2Texture *texture) 
{
  // Remove texture
  assert(texture->scene == this);
  assert(texture->scene_index >= 0);
  RNArrayEntry *entry = textures.KthEntry(texture->scene_index);
  assert(textures.EntryContents(entry) == texture);
  R2Texture *tail = textures.Tail();
  textures.EntryContents(entry) = tail;
  tail->scene_index = texture->scene_index;
  textures.RemoveTail();
  texture->scene_index = -1;
  texture->scene = NULL;
}



void R3Scene::
InsertReferencedScene(R3Scene *referenced_scene) 
{
  // Insert referenced scene
  referenced_scenes.Insert(referenced_scene);
}



void R3Scene::
RemoveReferencedScene(R3Scene *referenced_scene) 
{
  // Remove referenced scene
  referenced_scenes.Remove(referenced_scene);
}



void R3Scene::
InsertInfo(const char *key, const char *value) 
{
  // Insert key-value pair
  info.Insert(key, value);
}



void R3Scene::
ReplaceInfo(const char *key, const char *value) 
{
  // Replace key-value pair
  info.Replace(key, value);
}



void R3Scene::
RemoveInfo(const char *key) 
{
  // Insert key-value pair
  info.Remove(key);
}



void R3Scene::
SetCamera(const R3Camera& camera) 
{
  // Remember camera
  viewer.SetCamera(camera);
}



void R3Scene::
SetViewport(const R2Viewport& viewport) 
{
  // Remember viewport
  viewer.SetViewport(viewport);
}



void R3Scene::
SetViewer(const R3Viewer& viewer) 
{
  // Remember viewer
  this->viewer = viewer;
}



void R3Scene::
SetName(const char *name)
{
  // Set name
  if (this->name) free(this->name);
  if (name) this->name = strdup(name);
  else this->name = NULL;
}



void R3Scene::
SetData(void *data)
{
  // Set data
  this->data = data;
}



static void
CopyScene(R3Scene *src_scene, R3Scene *dst_scene,
  R3SceneNode *dst_root_node = NULL, const RNArray<R3Material *> *materials = NULL)
{
  // Get/check inputs
  if (!src_scene || !dst_scene) return;
  if (!dst_root_node) dst_root_node = dst_scene->Root();
  
  // Copy lights
  RNArray<R3Light *> dst_lights;
  for (int i = 0; i < src_scene->NLights(); i++) {
    R3Light *src_light = src_scene->Light(i);
    R3Light *dst_light = src_light->Copy();
    dst_scene->InsertLight(dst_light);
    dst_lights.Insert(dst_light);
  }

  // Copy textures
  RNArray<R2Texture *> dst_textures;
  for (int i = 0; i < src_scene->NTextures(); i++) {
    R2Texture *src_texture = src_scene->Texture(i);
    R2Texture *dst_texture = new R2Texture(*src_texture);
    dst_textures.Insert(dst_texture);
    dst_scene->InsertTexture(dst_texture);
  }

  // Copy brdfs
  RNArray<R3Brdf *> dst_brdfs;
  for (int i = 0; i < src_scene->NBrdfs(); i++) {
    R3Brdf *src_brdf = src_scene->Brdf(i);
    R3Brdf *dst_brdf = new R3Brdf(*src_brdf);
    dst_brdfs.Insert(dst_brdf);
    dst_scene->InsertBrdf(dst_brdf);
  }

  // Copy materials
  RNArray<R3Material *> dst_materials;
  for (int i = 0; i < src_scene->NMaterials(); i++) {
    R3Material *dst_material = (i < materials->NEntries()) ? materials->Kth(i) : NULL;
    if (!dst_material) {
      R3Material *src_material = src_scene->Material(i);
      R3Brdf *dst_brdf = (src_material->Brdf()) ? dst_brdfs[src_material->Brdf()->SceneIndex()] : NULL;
      R2Texture *dst_texture = (src_material->Texture()) ? dst_textures[src_material->Texture()->SceneIndex()] : NULL;
      dst_material = new R3Material(dst_brdf, dst_texture, src_material->Name());
      dst_scene->InsertMaterial(dst_material);
    }
    dst_materials.Insert(dst_material);
  }

  // Copy nodes 
  RNArray<R3SceneNode *> dst_nodes;
  for (int i = 0; i < src_scene->NNodes(); i++) {
    R3SceneNode *src_node = src_scene->Node(i);
    R3SceneNode *dst_node = new R3SceneNode(*src_node);
    dst_scene->InsertNode(dst_node);
    dst_nodes.Insert(dst_node);

    // Copy elements
    for (int i = 0; i < src_node->NElements(); i++) {
      R3SceneElement *src_element = src_node->Element(i);
      R3Material *src_material = src_element->Material();
      R3Material *dst_material = (src_material) ? dst_materials[src_material->SceneIndex()] : NULL;
      R3SceneElement *dst_element = new R3SceneElement(dst_material);
      dst_node->InsertElement(dst_element);
      for (int j = 0; j < src_element->NShapes(); j++) {
        R3Shape *src_shape = src_element->Shape(j);
        R3Shape *dst_shape = src_shape->Copy();
        dst_element->InsertShape(dst_shape);
      }
    }

    // Copy references
    for (int i = 0; i < src_node->NReferences(); i++) {
      R3SceneReference *src_reference = src_node->Reference(i);
      R3Scene *src_referenced_scene = src_reference->ReferencedScene();
      RNArray<R3Material *> dst_reference_materials;
      for (int j = 0; j < src_reference->NMaterials(); j++) {
        R3Material *src_material = src_reference->Material(j);
        R3Material *dst_material = dst_materials[src_material->SceneIndex()];
        dst_scene->InsertMaterial(dst_material);
      }
      CopyScene(src_referenced_scene, dst_scene, dst_node, &dst_reference_materials);
    }
  }
  
  // Copy node hierarchy
  for (int i = 0; i < src_scene->NNodes(); i++) {
    R3SceneNode *src_node = src_scene->Node(i);
    R3SceneNode *dst_node = dst_nodes.Kth(i);
    R3SceneNode *src_parent = src_node->Parent();
    R3SceneNode *dst_parent = (src_parent) ? dst_nodes[src_parent->SceneIndex()] : dst_root_node;
    dst_parent->InsertChild(dst_node);
  }
}

  

static void
R3SceneRemoveReferences(R3Scene *scene, R3SceneNode *node)
{
  // Recurse to children
  for (int i = 0; i < node->NChildren(); i++) {
    R3SceneNode *child = node->Child(i);
    R3SceneRemoveReferences(scene, child);
  }

  // Insert copy of all referenced scenes
  for (int i = 0; i < node->NReferences(); i++) {
    R3SceneReference *reference = node->Reference(i);
    R3Scene *referenced_scene = reference->ReferencedScene();
    CopyScene(referenced_scene, scene, node, &(reference->Materials()));
  }

  // Remove references
  while (node->NReferences() > 0) {
    R3SceneReference *reference = node->Reference(0);
    node->RemoveReference(reference);
    delete reference;
  }
}



void R3Scene::
RemoveReferences(void)
{
  // Replace all references with copies of referenced scenes
  R3SceneRemoveReferences(this, root);

  // Remove/delete referenced scenes
  while (NReferencedScenes() > 0) {
    R3Scene *referenced_scene = ReferencedScene(0);
    RemoveReferencedScene(referenced_scene);
    delete referenced_scene;
  }
}



static void
R3SceneRemoveHierarchy(R3Scene *scene, R3SceneNode *node, const R3Affine& parent_transformation)
{
  // Compute transformation
  R3Affine transformation = R3identity_affine;
  transformation.Transform(parent_transformation);
  transformation.Transform(node->Transformation());

  // Recurse to children
  for (int i = 0; i < node->NChildren(); i++) {
    R3SceneNode *child = node->Child(i);
    R3SceneRemoveHierarchy(scene, child, transformation);
  }

  // Check if should move node to root
  R3SceneNode *root = scene->Root();
  if (node != root) {
    R3SceneNode *parent = node->Parent();
    if (parent && (parent != root)) {
      // Assign transformation
      node->SetTransformation(transformation);

      // Move node to be child of root
      parent->RemoveChild(node);
      root->InsertChild(node);
    }
  }
}



void R3Scene::
RemoveHierarchy(void)
{
  // Maintain topology of scene, but set all node transformations to identity
  R3SceneRemoveHierarchy(this, root, R3identity_affine);

  // Set root node transformation to identity
  root->SetTransformation(R3identity_affine);
}



static void
R3SceneRemoveTransformations(R3Scene *scene, R3SceneNode *node, const R3Affine& parent_transformation)
{
  // Compute transformation
  R3Affine transformation = R3identity_affine;
  transformation.Transform(parent_transformation);
  transformation.Transform(node->Transformation());

  // Check if node has references
  if (node->NReferences() > 0) {
    node->SetTransformation(transformation);
    return;
  }

  // Transform elements
  for (int i = 0; i < node->NElements(); i++) {
    R3SceneElement *element = node->Element(i);
    element->Transform(transformation);
  }

  // Recurse to children
  for (int i = 0; i < node->NChildren(); i++) {
    R3SceneNode *child = node->Child(i);
    R3SceneRemoveTransformations(scene, child, transformation);
  }

  // Remove transformation
  node->SetTransformation(R3identity_affine);
}



void R3Scene::
RemoveTransformations(void)
{
  // Maintain topology of scene, but set all node transformations to identity
  R3SceneRemoveTransformations(this, root, R3identity_affine);
}



static void
R3SceneSubdivideTriangles(R3Scene *scene, R3SceneNode *node, RNLength max_edge_length)
{
  // Check max edge length
  if (RNIsNegativeOrZero(max_edge_length)) return;

  // Transform max_edge_length
  RNScalar scale = node->Transformation().ScaleFactor();
  if (RNIsNotZero(scale)) max_edge_length /= scale;

  // Subdivide triangles in elements
  for (int i = 0; i < node->NElements(); i++) {
    R3SceneElement *element = node->Element(i);
    for (int j = 0; j < element->NShapes(); j++) {
      R3Shape *shape = element->Shape(j);
      if (shape->ClassID() == R3TriangleArray::CLASS_ID()) {
        R3TriangleArray *triangles = (R3TriangleArray *) shape;
        triangles->Subdivide(max_edge_length);
      }
    }
  }

  // Subdivide triangles in references
  for (int i = 0; i < node->NReferences(); i++) {
    R3SceneReference *reference = node->Reference(i);
    R3Scene *referenced_scene = reference->ReferencedScene();
    R3SceneSubdivideTriangles(referenced_scene, referenced_scene->Root(), max_edge_length);
  }

  // Recurse to children
  for (int i = 0; i < node->NChildren(); i++) {
    R3SceneNode *child = node->Child(i);
    R3SceneSubdivideTriangles(scene, child, max_edge_length);
  }
}



void R3Scene::
SubdivideTriangles(RNLength max_edge_length)
{
  // Subdivide triangles until none is longer than max edge length
  R3SceneSubdivideTriangles(this, root, max_edge_length);
}




void R3Scene::
CreateDirectionalLights(void)
{
  // Create directional light pointing diagonally down
  RNRgb color1(1,1,1);
  R3Vector direction1(-3,-4,-5);
  direction1.Normalize();
  R3DirectionalLight *light1 = new R3DirectionalLight(direction1, color1);
  InsertLight(light1);

  // Create second (weaker) directional light pointing diagonally up
  RNRgb color2(0.5, 0.5, 0.5);
  R3Vector direction2(3,2,3);
  direction2.Normalize();
  R3DirectionalLight *light2 = new R3DirectionalLight(direction2, color2);
  InsertLight(light2);
}



RNLength R3Scene::
Distance(const R3Point& point) const
{
  // Find distance to root node
  return root->Distance(point);
}



RNBoolean R3Scene::
FindClosest(const R3Point& point,
  R3SceneNode **hit_node, R3Material **hit_material, R3Shape **hit_shape,
  R3Point *hit_point, R3Vector *hit_normal, RNScalar *hit_d,
  RNScalar min_d, RNScalar max_d) const
{
  // Find closest point in root node
  return root->FindClosest(point, hit_node, hit_material, hit_shape, hit_point, hit_normal, hit_d, min_d, max_d);
}



RNBoolean R3Scene::
Intersects(const R3Ray& ray,
  R3SceneNode **hit_node, R3Material **hit_material, R3Shape **hit_shape,
  R3Point *hit_point, R3Vector *hit_normal, RNScalar *hit_t,
  RNScalar min_t, RNScalar max_t) const
{
  // Intersect with root node
  return root->Intersects(ray, hit_node, hit_material, hit_shape, hit_point, hit_normal, hit_t, min_t, max_t);
}



int R3Scene::
LoadLights(int min_index, int max_index) const
{
  // Set ambient light
  static GLfloat ambient[4];
  ambient[0] = Ambient().R();
  ambient[1] = Ambient().G();
  ambient[2] = Ambient().B();
  ambient[3] = 1;
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);

  // Draw lights 
  int count = 0;
  for (int i = 0; i < NLights(); i++) {
    R3Light *light = Light(i);
    if (min_index + count > max_index) break;
    light->Draw(min_index + count);
    count++;
  }

  // Return number of lights loaded
  return count;
}



int R3Scene::
LoadLights(const R3Box& world_bbox, int min_index, int max_index) const
{
  // Determine the maximum number of lights
  if (NLights() == 0) return 1;
  int max_lights = max_index - min_index + 1;
  if (max_lights == 0) return 1;

  // Allocate scores for lights
  RNScalar *scores = new RNScalar [ NLights() ];
  for (int i = 0; i < NLights(); i++) scores[i] = 0;
  
  // Score lights
  for (int i = 0; i < NLights(); i++) {
    R3Light *light = Light(i);
    R3Sphere sphere = light->SphereOfInfluence(1E-3);
    if (RNIsZero(sphere.Radius())) continue; 
    RNScalar d = R3Distance(sphere.Centroid(), world_bbox);
    if (RNIsZero(d)) d = RN_EPSILON;
    RNScalar score = sphere.Radius() / d;
    if (score < 1.0) continue; 
    scores[i] = score;
  }

  // Sort lights based on score
  RNArray<R3Light *> sorted_lights = lights;
  for (int i = 0; i < sorted_lights.NEntries(); i++) {
    for (int j = i+1; j < sorted_lights.NEntries(); j++) {
      if (scores[j] < scores[i]) {
        sorted_lights.Swap(i, j);
        RNScalar swap = scores[i];
        scores[i] = scores[j];
        scores[j] = swap;
      }    
    }
  }

  // Load top lights
  sorted_lights.Truncate(max_lights);
  for (int i = 0; i < sorted_lights.NEntries(); i++) {
    R3Light *light = sorted_lights.Kth(i);
    light->Draw(min_index + i);
  }

  // Delete scores
  delete [] scores;

  // Return success
  return 1;
}



static void 
DrawNodeWithLights(const R3Scene *scene, R3SceneNode *node,
  const R3Affine& parent_transformation, const R3DrawFlags draw_flags,
  int min_light = 0, int max_light = 7)
{
  // Load lights
  if ((node->NElements() > 0) || (node->NReferences() > 0)) {
    R3Box world_bbox = node->BBox();
    world_bbox.Transform(parent_transformation);
    scene->LoadLights(world_bbox, min_light, max_light);
  }
  
  // Push transformation
  R3Affine transformation = parent_transformation;
  transformation.Transform(node->Transformation());
  transformation.Push();

  // Draw elements
  for (int i = 0; i < node->NElements(); i++) {
    R3SceneElement *element = node->Element(i);
    element->Draw(draw_flags);
  }

  // Draw references
  for (int i = 0; i < node->NReferences(); i++) {
    R3SceneReference *reference = node->Reference(i);
    reference->Draw(draw_flags);
  }

  // Pop transformation
  transformation.Pop();

  // Draw children
  for (int i = 0; i < node->NChildren(); i++) {
    R3SceneNode *child = node->Child(i);
    DrawNodeWithLights(scene, child, transformation, draw_flags, min_light, max_light);
  }
}



void R3Scene::
Draw(const R3DrawFlags draw_flags) const
{
  // Draw null material
  R3null_material.Draw();

  // Draw nodes recursively
  if (NLights() > 7) DrawNodeWithLights(this, root, R3identity_affine, draw_flags, 1, 7);
  else root->Draw(draw_flags);

  // Draw null material
  R3null_material.Draw();
}




////////////////////////////////////////////////////////////////////////
// SCENE I/O FUNCTIONS
////////////////////////////////////////////////////////////////////////

int R3Scene::
ReadFile(const char *filename, R3SceneNode *parent_node)
{
  // Parse input filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    fprintf(stderr, "Filename %s has no extension (e.g., .txt)\n", filename);
    return 0;
  }

  // Read file of appropriate type
  if (!strncmp(extension, ".scn", 4)) {
    if (!ReadPrincetonFile(filename, parent_node)) return 0;
  }
  else if (!strncmp(extension, ".ssc", 4)) {
    if (!ReadParseFile(filename, parent_node)) return 0;
  }
  else if (!strncmp(extension, ".txt", 4)) {
    if (!ReadSupportHierarchyFile(filename, parent_node)) return 0;
  }
  else if (!strncmp(extension, ".hier", 5)) {
    if (!ReadGrammarHierarchyFile(filename, parent_node)) return 0;
  }
  else if (!strncmp(extension, ".obj", 4)) {
    if (!ReadObjFile(filename, parent_node)) return 0;
  }
  else if (!strncmp(extension, ".ply", 4)) {
    if (!ReadPlyFile(filename, parent_node)) return 0;
  }
  else if (!strncmp(extension, ".off", 4)) {
    if (!ReadMeshFile(filename, parent_node)) return 0;
  }
  else if (!strncmp(extension, ".rct", 4)) {
    if (!ReadRectangleFile(filename, parent_node)) return 0;
  }
  else if (!strncmp(extension, ".json", 5)) {
    if (!ReadSUNCGFile(filename, parent_node)) return 0;
  }
  else {
    fprintf(stderr, "Unable to read file %s (unrecognized extension: %s)\n", filename, extension);
    return 0;
  }

  // Provide default camera
  if (Camera() == R3default_camera) {
    double scene_radius = BBox().DiagonalRadius();
    R3Point scene_center = BBox().Centroid();
    R3Vector towards = R3Vector(0, 0, -1);
    R3Vector up = R3Vector(0, 1, 0);
    R3Point eye = scene_center;
    R3Camera camera(eye, towards, up, 0.25, 0.25, 0.01 * scene_radius, 100 * scene_radius);
    SetCamera(camera);
  }

  // Set filename
  SetFilename(filename);

  // Return success
  return 1;
}



int R3Scene::
WriteFile(const char *filename) const
{
  // Parse input filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    fprintf(stderr, "Filename %s has no extension (e.g., .txt)\n", filename);
    return 0;
  }

  // Write file of appropriate type
  if (!strncmp(extension, ".scn", 4)) {
    if (!WritePrincetonFile(filename)) return 0;
  }
  else if (!strncmp(extension, ".obj", 4)) {
    if (!WriteObjFile(filename)) return 0;
  }
  else if (!strncmp(extension, ".json", 4)) {
    if (!WriteSUNCGFile(filename)) return 0;
  }
  else {
    fprintf(stderr, "Unable to write file %s (unrecognized extension: %s)\n", filename, extension);
    return 0;
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// OBJ FILE I/O FUNCTIONS
////////////////////////////////////////////////////////////////////////

static int
InsertSceneElement(R3Scene *scene, R3SceneNode *node, R3Material *material, 
  const RNArray<R3TriangleVertex *>& verts, const RNArray<R3Triangle *>& tris)
{
  // Create material if none
  if (!material) {
    R3Brdf *brdf = new R3Brdf(R3default_brdf, "Default");
    scene->InsertBrdf(brdf);
    material = new R3Material(brdf, "Default");
    scene->InsertMaterial(material);
  }
  
  // Find element
  R3SceneElement *element = NULL;
  for (int i = 0; i < node->NElements(); i++) {
    R3SceneElement *g = node->Element(i);
    if (g->Material() == material) { element = g; break; }
  }

  // Create element
  if (!element) {
    element = new R3SceneElement(material);
  }

  // Set vertex marks
  for (int i = 0; i < verts.NEntries(); i++) {
    R3TriangleVertex *vertex = verts.Kth(i);
    vertex->SetMark(0);
  }

  // Create copies of verts and tris for this element
  RNArray<R3Triangle *> element_tris;
  RNArray<R3TriangleVertex *> element_verts;
  for (int i = 0; i < tris.NEntries(); i++) {
    R3Triangle *triangle = tris.Kth(i);
    R3TriangleVertex *element_vertex[3];
    for (int j = 0; j < 3; j++) {
      R3TriangleVertex *vertex = triangle->Vertex(j);
      if (vertex->Mark() == 0) {
        element_verts.Insert(new R3TriangleVertex(*vertex));
        vertex->SetMark(element_verts.NEntries());
      }
      element_vertex[j] = element_verts.Kth(vertex->Mark()-1); 
    }
    R3Triangle *element_triangle = new R3Triangle(element_vertex[0], element_vertex[1], element_vertex[2]);
    element_tris.Insert(element_triangle);
  }

  // Create shape (triangle array)
  R3TriangleArray *shape = new R3TriangleArray(element_verts, element_tris);
        
  // Insert shape
  element->InsertShape(shape);
        
  // Insert element
  node->InsertElement(element);

  // Return success
  return 1;
}



static int
ReadObjMtlFile(R3Scene *scene, const char *dirname, const char *mtlname, RNArray<R3Material *> *returned_materials)
{
  // Open file
  char filename[1024];
  if (dirname) sprintf(filename, "%s/%s", dirname, mtlname);
  else strncpy(filename, mtlname, 1024);
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open file %s", filename);
    return 0;
  }

  // Initialize returned materials
  if (returned_materials) returned_materials->Empty();

  // Parse file
  char buffer[1024];
  int line_count = 0;
  R3Brdf *brdf = NULL;
  R3Material *material = NULL;
  RNSymbolTable<R2Texture *> texture_symbol_table;
  while (fgets(buffer, 1023, fp)) {
    // Increment line counter
    line_count++;

    // Skip white space
    char *bufferp = buffer;
    while (isspace(*bufferp)) bufferp++;

    // Skip blank lines and comments
    if (*bufferp == '#') continue;
    if (*bufferp == '\0') continue;

    // Get keyword
    char keyword[80];
    if (sscanf(bufferp, "%s", keyword) != 1) {
      fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
      return 0;
    }

    // Check keyword
    if (!strcmp(keyword, "newmtl")) {
      // Parse line
      char name[1024];
      if (sscanf(bufferp, "%s%s", keyword, name) != (unsigned int) 2) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Create new material
      brdf = new R3Brdf(name);
      brdf->SetDiffuse(RNwhite_rgb);
      scene->InsertBrdf(brdf);
      material = new R3Material(brdf, NULL, name);
      scene->InsertMaterial(material);
      if (returned_materials) returned_materials->Insert(material);
    }
    else if (!strcmp(keyword, "Ka")) {
      // Parse line
      double r, g, b;
      if (sscanf(bufferp, "%s%lf%lf%lf", keyword, &r, &g, &b) != (unsigned int) 4) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Set ambient reflectance 
      if (material && brdf) {
        brdf->SetAmbient(RNRgb(r, g, b));
        material->Update();
      }
    }
    else if (!strcmp(keyword, "Kd")) {
      // Parse line
      double r, g, b;
      if (sscanf(bufferp, "%s%lf%lf%lf", keyword, &r, &g, &b) != (unsigned int) 4) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Set diffuse reflectance 
      if (material && brdf) {
        brdf->SetDiffuse(RNRgb(r, g, b));
        material->Update();
      }
    }
    else if (!strcmp(keyword, "Ks")) {
      // Parse line
      double r, g, b;
      if (sscanf(bufferp, "%s%lf%lf%lf", keyword, &r, &g, &b) != (unsigned int) 4) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Set specular reflectance 
      if (material && brdf) {
        brdf->SetSpecular(RNRgb(r, g, b));
        material->Update();
      }
    }
    else if (!strcmp(keyword, "Ke")) {
      // Parse line
      double r, g, b;
      if (sscanf(bufferp, "%s%lf%lf%lf", keyword, &r, &g, &b) != (unsigned int) 4) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Set emission
      if (material && brdf) {
        brdf->SetEmission(RNRgb(r, g, b));
        material->Update();
      }
    }
    else if (!strcmp(keyword, "Tf")) {
      // Parse line
      double r, g, b;
      if (sscanf(bufferp, "%s%lf%lf%lf", keyword, &r, &g, &b) != (unsigned int) 4) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Set transmission
      if (material && brdf) {
        brdf->SetTransmission(RNRgb(r, g, b));
        material->Update();
      }
    }
    else if (!strcmp(keyword, "Tr")) {
      // Parse line
      double transparency;
      if (sscanf(bufferp, "%s%lf", keyword, &transparency) != (unsigned int) 2) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Set transmission
      if (material && brdf) {
        brdf->SetOpacity(1 - transparency);
        material->Update();
      }
    }
    else if (!strcmp(keyword, "d")) {
      // Parse line
      double opacity;
      if (sscanf(bufferp, "%s%lf", keyword, &opacity) != (unsigned int) 2) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Set transmission
      if (material && brdf) {
        brdf->SetOpacity(opacity);
        material->Update();
      }
    }
    else if (!strcmp(keyword, "Ns")) {
      // Parse line
      double ns;
      if (sscanf(bufferp, "%s%lf", keyword, &ns) != (unsigned int) 2) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Set shininess
      if (material && brdf) {
        brdf->SetShininess(ns);
        material->Update();
      }
    }
    else if (!strcmp(keyword, "Ni")) {
      // Parse line
      double index_of_refraction;
      if (sscanf(bufferp, "%s%lf", keyword, &index_of_refraction) != (unsigned int) 2) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Set index of refraction
      if (material && brdf) {
        brdf->SetIndexOfRefraction(index_of_refraction);
        material->Update();
      }
    }
    else if (!strcmp(keyword, "map_Kd")) {
      // Parse line
      char texture_name[1024];
      if (sscanf(bufferp, "%s%s", keyword, texture_name) != (unsigned int) 2) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Set texture
      if (material) {
        char texture_filename[1024];
        if (dirname) sprintf(texture_filename, "%s/%s", dirname, texture_name);
        else strncpy(texture_filename, texture_name, 1024);
        R2Texture *texture = NULL;
        if (!texture_symbol_table.Find(texture_filename, &texture)) {
          R2Image *image = new R2Image();
          if (!image->Read(texture_filename)) return 0;
          texture = new R2Texture(image);
          texture_symbol_table.Insert(texture_filename, texture);
          texture->SetFilename(texture_filename);
          texture->SetName(texture_filename);
          scene->InsertTexture(texture);
        }
        material->SetTexture(texture);
        material->Update();
      }
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



static int
ReadObj(R3Scene *scene, R3SceneNode *node, const char *dirname, FILE *fp, RNArray<R3Material *> *returned_materials = NULL)
{
  // Read body
  char buffer[1024];
  int line_count = 0;
  R3Material *material = NULL;
  RNSymbolTable<R3Material *> material_symbol_table;
  RNArray<R2Point *> texture_coords;
  RNArray<R3Vector *> normals;
  RNArray<R3Triangle *> tris;
  RNArray<R3TriangleVertex *> verts;
  RNArray<R3TriangleVertex *> tmp_verts;
  R3SceneNode *top_node = node;
  while (fgets(buffer, 1023, fp)) {
    // Increment line counter
    line_count++;

    // Skip white space
    char *bufferp = buffer;
    while (isspace(*bufferp)) bufferp++;

    // Skip blank lines and comments
    if (*bufferp == '#') continue;
    if (*bufferp == '\0') continue;

    // Get keyword
    char keyword[80];
    if (sscanf(bufferp, "%s", keyword) != 1) {
      fprintf(stderr, "Syntax error on line %d in OBJ file", line_count);
      return 0;
    }

    // Check keyword
    if (!strcmp(keyword, "v")) {
      // Read vertex coordinates
      double x, y, z;
      if (sscanf(bufferp, "%s%lf%lf%lf", keyword, &x, &y, &z) != 4) {
        fprintf(stderr, "Syntax error on line %d in OBJ file", line_count);
        return 0;
      }

      // Create vertex
      R3TriangleVertex *vertex = new R3TriangleVertex(R3Point(x, y, z));
      vertex->SetSharedFlag();
      verts.Insert(vertex);
    }
    else if (!strcmp(keyword, "vt")) {
      // Read texture coordinates
      double u, v;
      if (sscanf(bufferp, "%s%lf%lf", keyword, &u, &v) != 3) {
        fprintf(stderr, "Syntax error on line %d in OBJ file", line_count);
        return 0;
      }

      // Create texture coordinates
      R2Point *vt = new R2Point(u, v);
      texture_coords.Insert(vt);
    }
    else if (!strcmp(keyword, "vn")) {
      // Read texture coordinates
      double x, y, z;
      if (sscanf(bufferp, "%s%lf%lf%lf", keyword, &x, &y, &z) != 4) {
        fprintf(stderr, "Syntax error on line %d in OBJ file", line_count);
        return 0;
      }

      // Create normal
      R3Vector *vn = new R3Vector(x, y, z);
      normals.Insert(vn);
    }
    else if (!strcmp(keyword, "f")) {
      // Read vertex indices
      int quad = 1;
      char s[4][128] = { { '\0' }, { '\0' }, { '\0' },{ '\0' } }; 
      if (sscanf(bufferp, "%s%s%s%s%s", keyword, s[0], s[1], s[2], s[3]) != 5) {
        quad = 0;;
        if (sscanf(bufferp, "%s%s%s%s", keyword, s[0], s[1], s[2]) != 4) {
          fprintf(stderr, "Syntax error on line %d in OBJ file", line_count);
          return 0;
        }
      }

      // Parse vertex indices
      int n = (quad) ? 4 : 3;
      R3TriangleVertex *v[4] = { NULL, NULL, NULL, NULL };
      for (int i = 0; i < n; i++) {
        char *sv = s[i];
        char *st = strchr(sv, '/');
        if (st) *(st++) = 0;
        char *sn = (st) ? strchr(st, '/') : NULL;
        if (sn) *(sn++) = 0;
        if (sn && (strlen(sn) == 0)) sn = NULL; 
        if (st && (strlen(st) == 0)) st = NULL;
        int vi = (sv) ? atoi(sv) : 0;
        int ti = (st) ? atoi(st) : 0;
        int ni = (sn) ? atoi(sn) : 0;
        v[i] = verts.Kth(vi-1);
        if ((ti > 0) && ((ti-1) < texture_coords.NEntries())) {
          R2Point texcoords = *(texture_coords.Kth(ti-1));
          if (!(v[i]->Flags()[R3_VERTEX_TEXTURE_COORDS_DRAW_FLAG])) v[i]->SetTextureCoords(texcoords);
          else if (!R2Contains(texcoords, v[i]->TextureCoords())) {
            v[i] = new R3TriangleVertex(v[i]->Position(), texcoords); 
            v[i]->SetSharedFlag();
            tmp_verts.Insert(v[i]);
          }
        }
        if ((ni > 0) && ((ni-1) < normals.NEntries())) {
          R3Vector normal = *(normals.Kth(ni-1));
          if (!(v[i]->Flags()[R3_VERTEX_NORMALS_DRAW_FLAG])) v[i]->SetNormal(normal);
          else if (!R3Contains(normal, v[i]->Normal())) {
            v[i] = new R3TriangleVertex(v[i]->Position(), normal, v[i]->TextureCoords());
            v[i]->SetSharedFlag();
            tmp_verts.Insert(v[i]);
          }
        }
      }

      // Check vertices
      if ((v[0] == v[1]) || (v[1] == v[2]) || (v[0] == v[2])) continue;
      if ((quad) && ((v[3] == v[0]) || (v[3] == v[1]) || (v[3] == v[2]))) quad = 0;

      // Create first triangle
      if (RNIsPositive(R3Distance(v[0]->Position(), v[1]->Position())) &&
          RNIsPositive(R3Distance(v[1]->Position(), v[2]->Position())) &&
          RNIsPositive(R3Distance(v[2]->Position(), v[0]->Position()))) {
        R3Triangle *triangle = new R3Triangle(v[0], v[1], v[2]);
        tris.Insert(triangle);
      }

      // Create second triangle
      if (quad) {
        if (RNIsPositive(R3Distance(v[0]->Position(), v[2]->Position())) &&
            RNIsPositive(R3Distance(v[2]->Position(), v[3]->Position())) &&
            RNIsPositive(R3Distance(v[0]->Position(), v[3]->Position()))) {
          R3Triangle *triangle = new R3Triangle(v[0], v[2], v[3]);
          tris.Insert(triangle);
        }
      }
    }
    else if (!strcmp(keyword, "mtllib")) {
      // Read fields
      char mtlname[1024];
      if (sscanf(bufferp, "%s%s", keyword, mtlname) != 2) {
        fprintf(stderr, "Syntax error on line %d in OBJ file", line_count);
        return 0;
      }

      // Read materials
      RNArray<R3Material *> parsed_materials;
      if (!ReadObjMtlFile(scene, dirname, mtlname, &parsed_materials)) return 0;
      if (returned_materials) *returned_materials = parsed_materials;

      // Fill symbol table
      material_symbol_table.Empty();
      for (int i = 0; i < parsed_materials.NEntries(); i++) {
        R3Material *material = parsed_materials.Kth(i);
        if (!material->Name()) continue;
        material_symbol_table.Insert(material->Name(), material);
      }
    }
    else if (!strcmp(keyword, "usemtl")) {
      // Read fields
      char mtlname[1024];
      if (sscanf(bufferp, "%s%s", keyword, mtlname) != 2) {
        fprintf(stderr, "Syntax error on line %d in OBJ file", line_count);
        return 0;
      }

      // Process triangles from previous material
      if ((verts.NEntries() > 0) && (tris.NEntries() > 0)) {
        InsertSceneElement(scene, node, material, verts, tris);
        for (int i = 0; i < tris.NEntries(); i++) delete tris[i];
        tris.Empty();
      }

      // Find material
      if (!material_symbol_table.Find(mtlname, &material)) {
        fprintf(stderr, "Unable to find material %s at on line %d in OBJ file\n", mtlname, line_count);
        return 0;
      }
    }
    else if (!strcmp(keyword, "g") || !strcmp(keyword, "o")) {
      // Read name
      char name[1024];
      if (sscanf(bufferp, "%s%s", keyword, name) != 2) {
        fprintf(stderr, "Syntax error on line %d in OBJ file", line_count);
        return 0;
      }

      // Process triangles from previous object
      if ((verts.NEntries() > 0) && (tris.NEntries() > 0)) {
        InsertSceneElement(scene, node, material, verts, tris);
        for (int i = 0; i < tris.NEntries(); i++) delete tris[i];
        tris.Empty();
      }

      // Create child node
      node = new R3SceneNode(scene);
      node->SetName(name);
      top_node->InsertChild(node);
    }
  }

  // Process triangles from previous material
  if ((verts.NEntries() > 0) && (tris.NEntries() > 0)) {
    InsertSceneElement(scene, node, material, verts, tris);
    for (int i = 0; i < tris.NEntries(); i++) delete tris[i];
    tris.Empty();
  }

  // Delete texture coordinates
  for (int i = 0; i < texture_coords.NEntries(); i++) {
    delete texture_coords.Kth(i);
  }

  // Delete normals
  for (int i = 0; i < normals.NEntries(); i++) {
    delete normals.Kth(i);
  }

  // Delete verts (copied in InsertSceneElement)
  for (int i = 0; i < verts.NEntries(); i++) {
    delete verts.Kth(i);
  }

  // Delete tmp verts (copied in InsertSceneElement)
  for (int i = 0; i < tmp_verts.NEntries(); i++) {
    delete tmp_verts.Kth(i);
  }

  // Return success
  return 1;
}



static int
ReadObj(R3Scene *scene, R3SceneNode *node, const char *filename, RNArray<R3Material *> *returned_materials = NULL)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open file %s", filename);
    return 0;
  }

  // Determine directory name (for texture image files)
  char *dirname = NULL;
  char buffer[1024];
  strncpy(buffer, filename, 1024);
  char *endp = strrchr(buffer, '/');
  if (!endp) endp = strrchr(buffer, '\\');
  if (endp) { *endp = '\0'; dirname = buffer; }

  // Read file
  if (!ReadObj(scene, node, dirname, fp, returned_materials)) {
    fprintf(stderr, "Unable to read OBJ file %s\n", filename);
    fclose(fp);
    return 0;
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R3Scene::
ReadObjFile(const char *filename, R3SceneNode *parent_node)
{
  // Check/set parent node
  if (!parent_node) parent_node = root;
  
  // Read obj file, and put contents in parent node 
  return ReadObj(this, parent_node, filename);
}



////



static int
WriteObjMtlFile(const R3Scene *scene, const char *dirname, const char *mtlname)
{
  // Open file
  char filename[1024];
  if (dirname) sprintf(filename, "%s/%s", dirname, mtlname);
  else strncpy(filename, mtlname, 1024);
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    fprintf(stderr, "Unable to open file %s", filename);
    return 0;
  }

  // Write default material
  fprintf(fp, "newmtl _DEFAULT_\n");
  fprintf(fp, "Kd 0 0 0\n");
  fprintf(fp, "\n");

  // Write materials
  for (int i = 0; i < scene->NMaterials(); i++) {
    R3Material *material = scene->Material(i);  

    // ??? CHANGE NAME OF MATERIAL SO THAT IT IS UNIQUE ???
    char matname[1024];
    if (!material->Name()) sprintf(matname, "m%d", i);
    else sprintf(matname, "m%d_%s", i, material->Name());
    material->SetName(matname);
    
    // Write new material command
    fprintf(fp, "newmtl %s\n", material->Name());

    // Write brdf
    if (material->Brdf()) {
      const R3Brdf *brdf = material->Brdf();
      if (brdf->IsAmbient()) fprintf(fp, "Ka %g %g %g\n", brdf->Ambient().R(), brdf->Ambient().G(), brdf->Ambient().B());
      if (brdf->IsDiffuse()) fprintf(fp, "Kd %g %g %g\n", brdf->Diffuse().R(), brdf->Diffuse().G(), brdf->Diffuse().B());
      if (brdf->IsSpecular()) fprintf(fp, "Ks %g %g %g\n", brdf->Specular().R(), brdf->Specular().G(), brdf->Specular().B());
      if (brdf->IsEmissive()) fprintf(fp, "Ke %g %g %g\n", brdf->Emission().R(), brdf->Emission().G(), brdf->Emission().B());
      if (brdf->IsTransparent()) fprintf(fp, "Tf %g %g %g\n", brdf->Transmission().R(), brdf->Transmission().G(), brdf->Transmission().B());
      if (RNIsNotEqual(brdf->Opacity(), 1.0)) fprintf(fp, "Tr %g\n", 1.0 - brdf->Opacity());
      if (RNIsNotEqual(brdf->Opacity(), 1.0)) fprintf(fp, "d %g\n", brdf->Opacity());
      if (RNIsNotZero(brdf->Shininess())) fprintf(fp, "Ns %g\n", brdf->Shininess());
      if (RNIsNotZero(brdf->IndexOfRefraction())) fprintf(fp, "Ni %g\n", brdf->IndexOfRefraction());
    }

    // Write texture
    if (material->Texture()) {
      const R2Texture *texture = material->Texture();

      // Check if texture has a filename associated with it 
      if (texture->Filename()) {
        // Write texture command to material file (USE SAME FILE AS READ WITHOUT RE-WRITING IT)
        fprintf(fp, "map_Kd %s\n", texture->Filename());
      }
      else {
        // Get texture filename
        char texture_filename[1024];
        const char *texture_extension = "jpg";
        if (dirname) sprintf(texture_filename, "%s/%s.%s", dirname, material->Name(), texture_extension);
        else sprintf(texture_filename, "%s.%s", material->Name(), texture_extension);
        ((R2Texture *) texture)->SetFilename(texture_filename);
        
        // Write texture file
        const R2Image *texture_image = texture->Image();
        texture_image->Write(texture_filename);
        
        // Write texture command to material file
        fprintf(fp, "map_Kd %s.%s\n", material->Name(), texture_extension);
      }
    }

    // Write new line after each material
    fprintf(fp, "\n");
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



static int
WriteObj(const R3Scene *scene, R3SceneNode *node, const R3Affine& transformation, int &ngroups, int& nvertices, int& nnormals, int& ntexture_coords, FILE *fp)
{
  // Write group name
  if (node->Name()) fprintf(fp, "g %s\n", node->Name());
  else fprintf(fp, "g GROUP_%d\n", ++ngroups);
  
  // Write elements
  for (int i = 0; i < node->NElements(); i++) {
    R3SceneElement *element = node->Element(i);

    // Write material
    R3Material *material = element->Material();
    if (material) fprintf(fp, "usemtl %s\n", material->Name());
    else fprintf(fp, "usemtl _DEFAULT_\n");
      
    // Write shapes
    for (int j = 0; j < element->NShapes(); j++) {
      R3Shape *shape = element->Shape(j);
      if (shape->ClassID() == R3TriangleArray::CLASS_ID()) {
        R3TriangleArray *triangles = (R3TriangleArray *) shape;
        if (triangles->NVertices() == 0) continue;
        
        // Allocate indices for positions, texture coordinates, and normals
        int *pi = new int [ triangles->NVertices() ];
        for (int k = 0; k < triangles->NVertices(); k++) pi[k] = -1;
        int *ni = new int [ triangles->NVertices() ];
        for (int k = 0; k < triangles->NVertices(); k++) ni[k] = -1;
        int *ti = new int [ triangles->NVertices() ];
        for (int k = 0; k < triangles->NVertices(); k++) ti[k] = -1;

        // Write vertices
        for (int k = 0; k < triangles->NVertices(); k++) {
          R3TriangleVertex *v = triangles->Vertex(k);
          v->SetMark(k);
          if (TRUE) {
            R3Point p = v->Position();
            p.Transform(transformation);
            fprintf(fp, "v %g %g %g\n", p.X(), p.Y(), p.Z());
            pi[k] = ++nvertices;
          }
          if (v->Flags()[R3_VERTEX_NORMALS_DRAW_FLAG]) {
            R3Vector n = v->Normal();
            n.Transform(transformation);
            fprintf(fp, "vn %g %g %g\n", n.X(), n.Y(), n.Z());
            ni[k] = ++nnormals;
          }
          if (v->Flags()[R3_VERTEX_TEXTURE_COORDS_DRAW_FLAG]) {
            R2Point t = v->TextureCoords();
            fprintf(fp, "vt %g %g\n", t.X(), t.Y());
            ti[k] = ++ntexture_coords;
          }
        }
       
        // Write triangles
        for (int k = 0; k < triangles->NTriangles(); k++) {
          R3Triangle *triangle = triangles->Triangle(k);
          R3TriangleVertex *v0 = triangle->V0();
          R3TriangleVertex *v1 = triangle->V1();
          R3TriangleVertex *v2 = triangle->V2();
          unsigned int i0 = v0->Mark();
          unsigned int i1 = v1->Mark();
          unsigned int i2 = v2->Mark();
          if (triangle->Flags()[R3_VERTEX_TEXTURE_COORDS_DRAW_FLAG] && triangle->Flags()[R3_VERTEX_NORMALS_DRAW_FLAG]) {
            if (transformation.IsMirrored()) fprintf(fp, "f %u/%u/%u %u/%u/%u %u/%u/%u\n", pi[i2], ti[i2], ni[i2], pi[i1], ti[i1], ni[i1], pi[i0], ti[i0], ni[i0]);
            else fprintf(fp, "f %u/%u/%u %u/%u/%u %u/%u/%u\n", pi[i0], ti[i0], ni[i0], pi[i1], ti[i1], ni[i1], pi[i2], ti[i2], ni[i2]);
          }
          else if (triangle->Flags()[R3_VERTEX_TEXTURE_COORDS_DRAW_FLAG]) {
            if (transformation.IsMirrored()) fprintf(fp, "f %u/%u %u/%u %u/%u\n",  pi[i2], ti[i2],  pi[i1], ti[i1],  pi[i0], ti[i0]);
            else fprintf(fp, "f %u/%u %u/%u %u/%u\n",  pi[i0], ti[i0],  pi[i1], ti[i1],  pi[i2], ti[i2]);
          }
          else if (triangle->Flags()[R3_VERTEX_NORMALS_DRAW_FLAG]) {
            if (transformation.IsMirrored()) fprintf(fp, "f %u//%u %u//%u %u//%u\n",  pi[i2], ni[i2],  pi[i1], ni[i1],  pi[i0], ni[i0]);
            else fprintf(fp, "f %u//%u %u//%u %u//%u\n",  pi[i0], ni[i0],  pi[i1], ni[i1],  pi[i2], ni[i2]);
          }
          else {
            if (transformation.IsMirrored()) fprintf(fp, "f %u %u %u\n",  pi[i2], pi[i1], pi[i0]);
            else fprintf(fp, "f %u %u %u\n", pi[i0], pi[i1], pi[i2]);
          }
        }

        // Delete indices for positions, texture coordinates, and normals
        delete [] pi;
        delete [] ni;
        delete [] ti;
      }
    }
  }

  // Write children
  for (int i = 0; i < node->NChildren(); i++) {
    R3SceneNode *child = node->Child(i);

    // Update transformation
    R3Affine child_transformation = R3identity_affine;
    child_transformation.Transform(transformation);
    child_transformation.Transform(child->Transformation());

    // Write child
    if (!WriteObj(scene, child, child_transformation, ngroups, nvertices, nnormals, ntexture_coords, fp)) return 0;
  }

  // Return success
  return 1;
}



static int 
WriteObj(const R3Scene *scene, R3SceneNode *node, const char *filename) 
{
  // Determine directory name (for texture image files)
  char *dirname = NULL;
  char buffer[1024];
  strncpy(buffer, filename, 1024);
  char *endp = strrchr(buffer, '/');
  if (!endp) endp = strrchr(buffer, '\\');
  if (!endp) strcpy(buffer, ".");
  else { *endp = '\0'; dirname = buffer; }

  // Determine material filename
  char mtl_filename[1024];
  const char *startp = strrchr(filename, '/');
  startp = (startp) ? startp+1 : filename;
  strncpy(mtl_filename, startp, 1024);
  int slen = strlen(mtl_filename);
  if (slen > 4) mtl_filename[slen-4] = '\0';
  strncat(mtl_filename, ".mtl", 1024);

  // Create directory
  if (dirname) {
    char mkdir_cmd[1024];
    sprintf(mkdir_cmd, "mkdir -p %s", dirname);
    system(mkdir_cmd);
  }

  // Open file
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    fprintf(stderr, "Unable to open file %s", filename);
    return 0;
  }

  // Write material file
  fprintf(fp, "mtllib %s\n", mtl_filename);
  if (!WriteObjMtlFile(scene, dirname, mtl_filename)) {
    fprintf(stderr, "Unable to write OBJ material file %s\n", mtl_filename);
    fclose(fp);
    return 0;
  }

  // Write nodes 
  int ngroups = 0;
  int nvertices = 0;
  int nnormals = 0;
  int ntexture_coords = 0;
  if (!WriteObj(scene, node, node->Transformation(), ngroups, nvertices, nnormals, ntexture_coords, fp)) {
    fprintf(stderr, "Unable to write OBJ file %s\n", filename);
    fclose(fp);
    return 0;
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R3Scene::
WriteObjFile(const char *filename) const
{
  // Remove references (note that this changes scene structure :( )
  ((R3Scene *) this)->RemoveReferences();

  // Write obj file
  return WriteObj(this, root, filename);
}



////////////////////////////////////////////////////////////////////////
// PLY FILE I/O FUNCTIONS
////////////////////////////////////////////////////////////////////////

int R3Scene::
ReadPlyFile(const char *filename, R3SceneNode *parent_node)
{
  // Check/set parent node
  if (!parent_node) parent_node = root;

  // Read ply file
  R3Mesh mesh;
  if (!mesh.ReadFile(filename)) {
    fprintf(stderr, "Unable to read mesh %s\n", filename);
    return 0;
  }
  
  // Create array of vertices
  RNArray<R3TriangleVertex *> vertices;
  for (int i = 0; i < mesh.NVertices(); i++) {
    R3MeshVertex *mesh_vertex = mesh.Vertex(i);
    const R3Point& position = mesh.VertexPosition(mesh_vertex);
    R3TriangleVertex *triangle_vertex = new R3TriangleVertex(position);
    if (!mesh.VertexColor(mesh_vertex).IsBlack()) triangle_vertex->SetColor(mesh.VertexColor(mesh_vertex));
    triangle_vertex->SetTextureCoords(mesh.VertexTextureCoords(mesh_vertex));
    triangle_vertex->SetSharedFlag();
    vertices.Insert(triangle_vertex);
  }

  // Check face segments and materials
  int max_face_segment = 0;
  int max_face_material = 0;
  for (int i = 0; i < mesh.NFaces(); i++) {
    R3MeshFace *face = mesh.Face(i);
    int face_segment = mesh.FaceSegment(face);
    int face_material = mesh.FaceMaterial(face);
    if (face_segment > max_face_segment) max_face_segment = face_segment;
    if (face_material > max_face_material) max_face_material = face_material;
  }

  // Create arrays of triangles
  char node_name[1024];
  RNArray<R3Triangle *> *tris_0_0 = NULL;
  RNSymbolTable<RNArray<R3Triangle *> *> triangles;
  for (int i = 0; i < mesh.NFaces(); i++) {
    R3MeshFace *mesh_face = mesh.Face(i);
    int face_segment = mesh.FaceSegment(mesh_face);
    if (face_segment < 0) face_segment = 0;
    int face_category = mesh.FaceCategory(mesh_face);
    if (face_category < 0) face_category = 0;
    int i0 = mesh.VertexID(mesh.VertexOnFace(mesh_face, 0));
    int i1 = mesh.VertexID(mesh.VertexOnFace(mesh_face, 1));
    int i2 = mesh.VertexID(mesh.VertexOnFace(mesh_face, 2));
    R3TriangleVertex *v0 = vertices.Kth(i0);
    R3TriangleVertex *v1 = vertices.Kth(i1);
    R3TriangleVertex *v2 = vertices.Kth(i2);
    R3Triangle *triangle = new R3Triangle(v0, v1, v2);
    RNArray<R3Triangle *> *tris = ((face_segment == 0) && (face_category == 0)) ? tris_0_0 : NULL; 
    if (!tris) {
      sprintf(node_name, "%d_%d", face_segment, face_category);
      if (!triangles.Find(node_name, &tris)) {
        tris = new RNArray<R3Triangle *>();
        if ((face_segment == 0) && (face_category == 0)) tris_0_0 = tris;
        triangles.Insert(node_name, tris);
      }
    }
    tris->Insert(triangle);
  }

  // Create node for each set of triangles
  std::map<std::string, RNArray<R3Triangle *> *, RNMapComparator< std::string > >::iterator it;
  for (it = triangles.m->begin(); it != triangles.m->end(); it++) {
    std::string node_name = it->first;
    RNArray<R3Triangle *> *tris = it->second;
    if (tris->IsEmpty()) { delete tris; continue; }
    R3SceneNode *node = new R3SceneNode(this);
    if (!InsertSceneElement(this, node, NULL, vertices, *tris)) return 0;
    node->SetName(node_name.c_str());
    const char *category_separator = strchr(node_name.c_str(), '_');
    if (category_separator) node->InsertInfo("index", strdup(category_separator+1));
    parent_node->InsertChild(node);
    for (int i = 0; i < tris->NEntries(); i++) delete tris->Kth(i);
    delete tris;
  }

  // Delete vertices 
  for (int i = 0; i < vertices.NEntries(); i++) delete vertices[i];
  
  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// MESH FILE I/O FUNCTIONS
////////////////////////////////////////////////////////////////////////

static R3TriangleArray *
ReadMesh(const char *filename)
{
  // Read mesh file
  R3Mesh mesh;
  if (!mesh.ReadFile(filename)) {
    fprintf(stderr, "Unable to read mesh %s\n", filename);
    return NULL;
  }
  
  // Create array of vertices
  RNArray<R3TriangleVertex *> vertices;
  for (int i = 0; i < mesh.NVertices(); i++) {
    R3MeshVertex *mesh_vertex = mesh.Vertex(i);
    const R3Point& position = mesh.VertexPosition(mesh_vertex);
    R3TriangleVertex *triangle_vertex = new R3TriangleVertex(position);
    vertices.Insert(triangle_vertex);

    // Check if should assign vertex normal
    RNBoolean smooth = TRUE;
    const RNAngle max_smooth_angle = RN_PI/5.0;
    for (int j = 0; j < mesh.VertexValence(mesh_vertex); j++) {
      R3MeshEdge *mesh_edge = mesh.EdgeOnVertex(mesh_vertex, j);
      RNAngle angle = mesh.EdgeInteriorAngle(mesh_edge);
      if (fabs(angle - RN_PI) > max_smooth_angle) { smooth = FALSE; break; }
    }

    // Assign vertex normal
    if (smooth) {
      const R3Vector& normal = mesh.VertexNormal(mesh_vertex);
      triangle_vertex->SetNormal(normal);
    }
    
    // Assign vertex color
    if (!mesh.VertexColor(mesh_vertex).IsBlack()) {
      triangle_vertex->SetColor(mesh.VertexColor(mesh_vertex));
    }

    // Assign vertex texture coordinates
    triangle_vertex->SetTextureCoords(mesh.VertexTextureCoords(mesh_vertex));
  }

  // Create array of triangles
  RNArray<R3Triangle *> triangles;
  for (int i = 0; i < mesh.NFaces(); i++) {
    R3MeshFace *mesh_face = mesh.Face(i);
    int i0 = mesh.VertexID(mesh.VertexOnFace(mesh_face, 0));
    int i1 = mesh.VertexID(mesh.VertexOnFace(mesh_face, 1));
    int i2 = mesh.VertexID(mesh.VertexOnFace(mesh_face, 2));
    R3TriangleVertex *v0 = vertices.Kth(i0);
    R3TriangleVertex *v1 = vertices.Kth(i1);
    R3TriangleVertex *v2 = vertices.Kth(i2);
    R3Triangle *triangle = new R3Triangle(v0, v1, v2);
    triangles.Insert(triangle);
  }

  // Return triangle array
  return new R3TriangleArray(vertices, triangles);
}

 

int R3Scene::
ReadMeshFile(const char *filename, R3SceneNode *parent_node)
{
  // Check/set parent node                                                                                                                                   
  if (!parent_node) parent_node = root;

  // Load triangles into scene                                                                                                                               
  R3TriangleArray *shape = ReadMesh(filename);
  if (!shape) return 0;
  R3SceneElement *element = new R3SceneElement();
  element->InsertShape(shape);
  R3SceneNode *node = new R3SceneNode(this);
  node->InsertElement(element);
  parent_node->InsertChild(node);

  // Return success                                                                                                                                          
  return 1;
}



////////////////////////////////////////////////////////////////////////
// PRINCETON SCENE FILE I/O FUNCTIONS
////////////////////////////////////////////////////////////////////////

static int
FindPrincetonMaterialAndElement(R3Scene *scene, R3SceneNode *node,
  const RNArray<R3Material *>& materials, int m, R3Material *&default_material,
  R3Material *& material, R3SceneElement *& element)
{
  // Find material from m
  material = NULL;
  if (m >= 0) {
    if (m < materials.NEntries()) material = materials[m];
    else return 0;
  }
  else {
    // Get material from 
    material = default_material;
    if (!material) {
      R3Brdf *brdf = new R3Brdf(R3default_brdf);
      scene->InsertBrdf(brdf);
      material = new R3Material(brdf);
      scene->InsertMaterial(material);
      default_material = material;
    }
  }

  // Find element with that material
  element = NULL;
  for (int i = 0; i < node->NElements(); i++) {
    R3SceneElement *e = node->Element(i);
    if (e->Material() == material) { element = e; break; }
  }

  // Create element if none found
  if (!element) {
    element = new R3SceneElement(material);
    node->InsertElement(element);
  }

  // Return success
  return 1;
}



int R3Scene::
ReadPrincetonFile(const char *filename, R3SceneNode *parent_node)
{
  // Get/set parent_node
  if (!parent_node) parent_node = root;
  
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "r"))) {
    fprintf(stderr, "Unable to open file %s", filename);
    return 0;
  }

  // Create array of materials
  RNArray<R3Material *> parsed_materials;
  R3Material *material = NULL;
  R3SceneElement *element = NULL;

  // Create stack of group information
  const int max_depth = 1024;
  R3SceneNode *group_nodes[max_depth] = { NULL };
  R3Material *group_materials[max_depth] = { NULL };
  group_nodes[0] = parent_node;
  int depth = 0;

  // Read body
  char cmd[128];
  int command_number = 1;
  while (fscanf(fp, "%s", cmd) == 1) {
    if (cmd[0] == '#') {
      // Comment -- read everything until end of line
      do { cmd[0] = fgetc(fp); } while ((cmd[0] >= 0) && (cmd[0] != '\n'));
    }
    else if (!strcmp(cmd, "tri")) {
      // Read data
      int m;
      R3Point p1, p2, p3;
      if (fscanf(fp, "%d%lf%lf%lf%lf%lf%lf%lf%lf%lf", &m, 
        &p1[0], &p1[1], &p1[2], &p2[0], &p2[1], &p2[2], &p3[0], &p3[1], &p3[2]) != 10) {
        fprintf(stderr, "Unable to read triangle at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Create triangle
      R3TriangleVertex *v1 = new R3TriangleVertex(p1);
      R3TriangleVertex *v2 = new R3TriangleVertex(p2);
      R3TriangleVertex *v3 = new R3TriangleVertex(p3);
      R3Triangle *triangle = new R3Triangle(v1, v2, v3);

      // Get material and element from m
      if (!FindPrincetonMaterialAndElement(this, group_nodes[depth], parsed_materials, m, group_materials[depth], material, element)) {
        fprintf(stderr, "Invalid material id at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Insert triangle into element
      element->InsertShape(triangle);
    }
    else if (!strcmp(cmd, "box")) {
      // Read data
      int m;
      R3Point p1, p2;
      if (fscanf(fp, "%d%lf%lf%lf%lf%lf%lf", &m, &p1[0], &p1[1], &p1[2], &p2[0], &p2[1], &p2[2]) != 7) {
        fprintf(stderr, "Unable to read box at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Sort coordinates
      if (p1[0] > p2[0]) { RNCoord swap = p1[0]; p1[0] = p2[0]; p2[0] = swap; }
      if (p1[1] > p2[1]) { RNCoord swap = p1[1]; p1[1] = p2[1]; p2[1] = swap; }
      if (p1[2] > p2[2]) { RNCoord swap = p1[2]; p1[2] = p2[2]; p2[2] = swap; }

      // Create box
      R3Box *box = new R3Box(p1, p2);

      // Get material and element from m
      if (!FindPrincetonMaterialAndElement(this, group_nodes[depth], parsed_materials, m, group_materials[depth], material, element)) {
        fprintf(stderr, "Invalid material id at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Insert shape into element
      element->InsertShape(box);
    }
    else if (!strcmp(cmd, "sphere")) {
      // Read data
      int m;
      R3Point c;
      double r;
      if (fscanf(fp, "%d%lf%lf%lf%lf", &m, &c[0], &c[1], &c[2], &r) != 5) {
        fprintf(stderr, "Unable to read sphere at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Create sphere
      R3Sphere *sphere = new R3Sphere(c, r);

      // Get material and element from m
      if (!FindPrincetonMaterialAndElement(this, group_nodes[depth], parsed_materials, m, group_materials[depth], material, element)) {
        fprintf(stderr, "Invalid material id at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Insert shape into element
      element->InsertShape(sphere);
    }
    else if (!strcmp(cmd, "cylinder")) {
      // Read data
      int m;
      R3Point c;
      double r, h;
      if (fscanf(fp, "%d%lf%lf%lf%lf%lf", &m, &c[0], &c[1], &c[2], &r, &h) != 6) {
        fprintf(stderr, "Unable to read cylinder at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Create cylinder
      R3Point p1 = c - 0.5 * h * R3posy_vector;
      R3Point p2 = c + 0.5 * h * R3posy_vector;
      R3Cylinder *cylinder = new R3Cylinder(p1, p2, r);

      // Get material and element from m
      if (!FindPrincetonMaterialAndElement(this, group_nodes[depth], parsed_materials, m, group_materials[depth], material, element)) {
        fprintf(stderr, "Invalid material id at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Insert shape into element
      element->InsertShape(cylinder);
    }
    else if (!strcmp(cmd, "cone")) {
      // Read data
      int m;
      R3Point c;
      double r, h;
      if (fscanf(fp, "%d%lf%lf%lf%lf%lf", &m, &c[0], &c[1], &c[2], &r, &h) != 6) {
        fprintf(stderr, "Unable to read cone at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Create cone
      R3Point p1 = c - 0.5 * h * R3posy_vector;
      R3Point p2 = c + 0.5 * h * R3posy_vector;
      R3Cone *cone = new R3Cone(p1, p2, r);

      // Get material and element from m
      if (!FindPrincetonMaterialAndElement(this, group_nodes[depth], parsed_materials, m, group_materials[depth], material, element)) {
        fprintf(stderr, "Invalid material id at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Insert shape into element
      element->InsertShape(cone);
    }
    else if (!strcmp(cmd, "mesh")) {
      // Read data
      int m;
      char meshname[256];
      if (fscanf(fp, "%d%s", &m, meshname) != 2) {
        fprintf(stderr, "Unable to parse mesh command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Get mesh filename
      char buffer[2048];
      strcpy(buffer, filename);
      char *bufferp = strrchr(buffer, '/');
      if (bufferp) *(bufferp+1) = '\0';
      else buffer[0] = '\0';
      strcat(buffer, meshname);

      // Read mesh
      R3TriangleArray *mesh = ReadMesh(buffer);
      if (!mesh) return 0;

      // Get material and element from m
      if (!FindPrincetonMaterialAndElement(this, group_nodes[depth], parsed_materials, m, group_materials[depth], material, element)) {
        fprintf(stderr, "Invalid material id at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Insert shape into element
      element->InsertShape(mesh);
    }
    else if (!strcmp(cmd, "line")) {
      // Read data
      int m;
      double x1, y1, z1, x2, y2, z2;
      if (fscanf(fp, "%d%lf%lf%lf%lf%lf%lf", &m, &x1, &y1, &z1, &x2, &y2, &z2) != 7) {
        fprintf(stderr, "Unable to read line at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Create cylinder representing line
      R3Point p1(x1, y1, z1);
      R3Point p2(x2, y2, z2);
      R3Cylinder *cylinder = new R3Cylinder(p1, p2, RN_BIG_EPSILON);

      // Get material and element from m
      if (!FindPrincetonMaterialAndElement(this, group_nodes[depth], parsed_materials, m, group_materials[depth], material, element)) {
        fprintf(stderr, "Invalid material id at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Insert shape into element
      element->InsertShape(cylinder);
    }
    else if (!strcmp(cmd, "begin") || !strcmp(cmd, "group")) {
      // Read data
      int m;
      double matrix[16];
      char group_name[4096] = { '\0' };
      if (!strcmp(cmd, "group")) fscanf(fp, "%s", group_name);
      if (fscanf(fp, "%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &m, 
        &matrix[0], &matrix[1], &matrix[2], &matrix[3], 
        &matrix[4], &matrix[5], &matrix[6], &matrix[7], 
        &matrix[8], &matrix[9], &matrix[10], &matrix[11], 
        &matrix[12], &matrix[13], &matrix[14], &matrix[15]) != 17) {
        fprintf(stderr, "Unable to read begin at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Get material from m
      if (m >= 0) {
        if (m < parsed_materials.NEntries()) material = parsed_materials[m];
        else material = NULL;
      }

      // Create new node
      R3SceneNode *node = new R3SceneNode(this);
      if (group_name[0]) node->SetName(group_name);
      node->SetTransformation(R3Affine(R4Matrix(matrix)));
      group_nodes[depth]->InsertChild(node);

      // Push node onto stack
      depth++;
      group_nodes[depth] = node;
      group_materials[depth] = material;
    }
    else if (!strcmp(cmd, "end")) {
      // Check stack depth
      if (depth <= 0) {
        fprintf(stderr, "Extra end statement at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Pop node from stack
      depth--;
    }
    else if (!strcmp(cmd, "material")) {
      // Read data
      RNRgb ka, kd, ks, kt, e;
      double n, ir;
      char texture_name[256];
      if (fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%s", 
          &ka[0], &ka[1], &ka[2], &kd[0], &kd[1], &kd[2], &ks[0], &ks[1], &ks[2], &kt[0], &kt[1], &kt[2], 
          &e[0], &e[1], &e[2], &n, &ir, texture_name) != 18) {
        fprintf(stderr, "Unable to read material at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Create brdf
      R3Brdf *brdf = new R3Brdf(ka, kd, ks, kt, e, n, ir);
      InsertBrdf(brdf);

      // Create texture
      R2Texture *texture = NULL;
      if (strcmp(texture_name, "0")) {
        // Get texture filename
        char buffer[2048];
        strcpy(buffer, filename);
        char *bufferp = strrchr(buffer, '/');
        if (bufferp) *(bufferp+1) = '\0';
        else buffer[0] = '\0';
        strcat(buffer, texture_name);

        // Read texture file
        R2Image *image = new R2Image();
        if (!image->Read(buffer)) {
          fprintf(stderr, "Unable to read texture from %s at command %d in file %s\n", buffer, command_number, filename);
          return 0;
        }
        
        // Create texture
        texture = new R2Texture(image);
        InsertTexture(texture);
      }

      // Create material
      R3Material *material = new R3Material(brdf, texture);
      InsertMaterial(material);
      parsed_materials.Insert(material);
    }
    else if (!strcmp(cmd, "dir_light")) {
      // Read data
      RNRgb c;
      R3Vector d;
      if (fscanf(fp, "%lf%lf%lf%lf%lf%lf", 
        &c[0], &c[1], &c[2], &d[0], &d[1], &d[2]) != 6) {
        fprintf(stderr, "Unable to read directional light at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Normalize direction
      d.Normalize();

      // Create directional light
      R3DirectionalLight *light = new R3DirectionalLight(d, c);
      InsertLight(light);
    }
    else if (!strcmp(cmd, "point_light")) {
      // Read data
      RNRgb c;
      R3Point p;
      double ca, la, qa;
      if (fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf%lf%lf", &c[0], &c[1], &c[2], &p[0], &p[1], &p[2], &ca, &la, &qa) != 9) {
        fprintf(stderr, "Unable to read point light at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Create point light
      R3PointLight *light = new R3PointLight(p, c, 1, TRUE, ca, la, qa);
      InsertLight(light);
    }
    else if (!strcmp(cmd, "spot_light")) {
      // Read data
      RNRgb c;
      R3Point p;
      R3Vector d;
      double ca, la, qa, sc, sd;
      if (fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", 
        &c[0], &c[1], &c[2], &p[0], &p[1], &p[2], &d[0], &d[1], &d[2], &ca, &la, &qa, &sc, &sd) != 14) {
        fprintf(stderr, "Unable to read point light at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Normalize direction
      d.Normalize();

      // Create spot light
      R3SpotLight *light = new R3SpotLight(p, d, c, sd, sc, 1, TRUE, ca, la, qa);
      InsertLight(light);
    }
    else if (!strcmp(cmd, "area_light")) {
      // Read data
      RNRgb c;
      R3Point p;
      R3Vector d;
      double radius, ca, la, qa;
      if (fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", 
        &c[0], &c[1], &c[2], &p[0], &p[1], &p[2], &d[0], &d[1], &d[2], &radius, &ca, &la, &qa) != 13) {
        fprintf(stderr, "Unable to read area light at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Normalize direction
      d.Normalize();

      // Create spot light
      R3AreaLight *light = new R3AreaLight(p, radius, d, c, 1, TRUE, ca, la, qa);
      InsertLight(light);
    }
    else if (!strcmp(cmd, "camera")) {
      // Read data
      R3Point e;
      R3Vector t, u;
      RNScalar xfov, neardist, fardist;
      if (fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &e[0], &e[1], &e[2], &t[0], &t[1], &t[2], &u[0], &u[1], &u[2], &xfov, &neardist, &fardist) != 12) {
        fprintf(stderr, "Unable to read camera at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Assign camera
      R3Camera camera(e, t, u, xfov, xfov, neardist, fardist);
      SetCamera(camera);
    }
    else if (!strcmp(cmd, "include")) {
      // Read data
      char scenename[256];
      if (fscanf(fp, "%s", scenename) != 1) {
        fprintf(stderr, "Unable to read include command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Get scene filename
      char buffer[2048];
      strcpy(buffer, filename);
      char *bufferp = strrchr(buffer, '/');
      if (bufferp) *(bufferp+1) = '\0';
      else buffer[0] = '\0';
      strcat(buffer, scenename);

      // Read scene from included file
      if (!ReadFile(buffer, group_nodes[depth])) {
        fprintf(stderr, "Unable to read included scene: %s\n", buffer);
        return 0;
      }
    }
    else if (!strcmp(cmd, "background")) {
      // Read data
      double r, g, b;
      if (fscanf(fp, "%lf%lf%lf", &r, &g, &b) != 3) {
        fprintf(stderr, "Unable to read background at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Assign background color
      SetBackground(RNRgb(r, g, b));
    }
    else if (!strcmp(cmd, "ambient")) {
      // Read data
      double r, g, b;
      if (fscanf(fp, "%lf%lf%lf", &r, &g, &b) != 3) {
        fprintf(stderr, "Unable to read ambient at command %d in file %s\n", command_number, filename);
        return 0;
      }

      // Assign ambient color
      SetAmbient(RNRgb(r, g, b));
    }
    else {
      fprintf(stderr, "Unrecognized command %d in file %s: %s\n", command_number, filename, cmd);
      return 0;
    }
	
    // Increment command number
    command_number++;
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



//////////////

static int
CreateTextureFile(const R3Scene *scene, const R2Texture *texture, const char *output_directory)
{
  // Return success
  return 1;
}



static int
WritePrincetonMaterial(const R3Scene *scene, const R3Material *material, FILE *fp)
{
  // Write material
  const R3Brdf *brdf = material->Brdf();
  const R2Texture *texture = material->Texture();
  RNRgb ambient = (brdf) ? brdf->Ambient() : RNblack_rgb;
  RNRgb diffuse = (brdf) ? brdf->Diffuse() : RNblack_rgb;
  RNRgb specular = (brdf) ? brdf->Specular() : RNblack_rgb;
  RNRgb transmission = (brdf) ? brdf->Transmission() : RNblack_rgb;
  RNRgb emission = (brdf) ? brdf->Emission() : RNblack_rgb;
  RNScalar shininess = brdf->Shininess();
  RNScalar indexofref = brdf->IndexOfRefraction();
  const char *texture_name = (texture && texture->Filename()) ? texture->Filename() : "0";
  fprintf(fp, "material %g %g %g  %g %g %g  %g %g %g  %g %g %g  %g %g %g  %g %g %s\n",
    ambient.R(), ambient.G(), ambient.B(),
    diffuse.R(), diffuse.G(), diffuse.B(),
    specular.R(), specular.G(), specular.B(),
    transmission.R(), transmission.G(), transmission.B(),
    emission.R(), emission.G(), emission.B(),
    shininess, indexofref, texture_name);

  // Return success
  return 1;
}



static int
WritePrincetonLight(const R3Scene *scene, const R3Light *light, FILE *fp)
{
  // Get light color
  RNRgb color = light->Color();
  
  // Write light in format appropriate for type
  if (light->ClassID() == R3DirectionalLight::CLASS_ID()) {
    R3DirectionalLight *directional_light = (R3DirectionalLight *) light;
    const R3Vector& direction = directional_light->Direction();
    fprintf(fp, "dir_light  %g %g %g  %g %g %g\n",
            color.R(), color.G(), color.B(),
            direction.X(), direction.Y(), direction.Z());
  }
  else if (light->ClassID() == R3PointLight::CLASS_ID()) {
    R3PointLight *point_light = (R3PointLight *) light;
    const R3Point& position = point_light->Position();
    RNScalar ca = point_light->ConstantAttenuation();
    RNScalar la = point_light->LinearAttenuation();
    RNScalar qa = point_light->QuadraticAttenuation();
    fprintf(fp, "point_light  %g %g %g  %g %g %g  %g %g %g\n",
            color.R(), color.G(), color.B(),
            position.X(), position.Y(), position.Z(),
            ca, la, qa);
  }
  else if (light->ClassID() == R3SpotLight::CLASS_ID()) {
    R3SpotLight *spot_light = (R3SpotLight *) light;
    const R3Point& position = spot_light->Position();
    const R3Vector& direction = spot_light->Direction();
    RNScalar ca = spot_light->ConstantAttenuation();
    RNScalar la = spot_light->LinearAttenuation();
    RNScalar qa = spot_light->QuadraticAttenuation();
    RNScalar sc = spot_light->CutOffAngle();
    RNScalar sd = spot_light->DropOffRate();
    fprintf(fp, "spot_light  %g %g %g  %g %g %g  %g %g %g  %g %g %g  %g %g\n",
            color.R(), color.G(), color.B(),
            position.X(), position.Y(), position.Z(),
            direction.X(), direction.Y(), direction.Z(),
            ca, la, qa, sc, sd);
  }
  else if (light->ClassID() == R3AreaLight::CLASS_ID()) {
    R3AreaLight *area_light = (R3AreaLight *) light;
    const R3Point& position = area_light->Position();
    const R3Vector& direction = area_light->Direction();
    RNScalar radius = area_light->Radius();
    RNScalar ca = area_light->ConstantAttenuation();
    RNScalar la = area_light->LinearAttenuation();
    RNScalar qa = area_light->QuadraticAttenuation();
    fprintf(fp, "area_light  %g %g %g  %g %g %g  %g %g %g  %g  %g %g %g\n",
            color.R(), color.G(), color.B(),
            position.X(), position.Y(), position.Z(),
            direction.X(), direction.Y(), direction.Z(),
            radius, ca, la, qa);
  }
  else {
    // fprintf("Unrecognized light type\n");
    // return 0;
  }

  // Return success
  return 1;
}



static int
WritePrincetonElement(const R3Scene *scene, R3SceneElement *element,
  const R3Affine& transformation, FILE *fp, const char *indent)
{
  // Get material index
  R3Material *material = element->Material();
  int material_index = (material) ? material->SceneIndex() : -1;
  
  // Write shapes
  for (int i = 0; i < element->NShapes(); i++) {
    R3Shape *shape = element->Shape(i);
    if (shape->ClassID() == R3Box::CLASS_ID()) {
      R3Box *box = (R3Box *) shape;
      fprintf(fp, "%sbox  %d   %g %g %g   %g %g %g\n", indent, material_index,
        box->XMin(), box->YMin(), box->ZMin(), box->XMax(), box->YMax(), box->ZMax());
    }
    else if (shape->ClassID() == R3Sphere::CLASS_ID()) {
      R3Sphere *sphere = (R3Sphere *) shape;
      R3Point center = sphere->Center();
      fprintf(fp, "%ssphere  %d   %g %g %g   %g\n", indent, material_index,
        center.X(), center.Y(), center.Z(), sphere->Radius());
    }
    else if (shape->ClassID() == R3Cone::CLASS_ID()) {
      R3Cone *cone = (R3Cone *) shape;
      R3Point center = cone->Axis().Midpoint();
      if (RNIsNotEqual(cone->Axis().Vector().Dot(R3posz_vector), 1.0)) fprintf(stderr, "Warning: cone not axis aligned\n");
      fprintf(fp, "%scone  %d   %g %g %g   %g %g\n", indent, material_index,
        center.X(), center.Y(), center.Z(), cone->Radius(), cone->Height());
    }
    else if (shape->ClassID() == R3Cylinder::CLASS_ID()) {
      R3Cylinder *cylinder = (R3Cylinder *) shape;
      R3Point center = cylinder->Axis().Midpoint();
      if (RNIsNotEqual(cylinder->Axis().Vector().Dot(R3posz_vector), 1.0)) fprintf(stderr, "Warning: cylinder not axis aligned\n");
      fprintf(fp, "%scylinder  %d   %g %g %g   %g %g\n", indent, material_index,
        center.X(), center.Y(), center.Z(), cylinder->Radius(), cylinder->Height());
    }
    else if (shape->ClassID() == R3Triangle::CLASS_ID()) {
      R3Triangle *triangle = (R3Triangle *) shape;
      R3TriangleVertex *v0 = (transformation.IsMirrored()) ? triangle->V2() : triangle->V0();
      R3TriangleVertex *v1 = triangle->V1();
      R3TriangleVertex *v2 = (transformation.IsMirrored()) ? triangle->V0() : triangle->V2();
      const R3Point& p0 = v0->Position();
      const R3Point& p1 = v1->Position();
      const R3Point& p2 = v2->Position();
      fprintf(fp, "%stri  %d   %g %g %g   %g %g %g   %g %g %g\n", indent, material_index,
        p0.X(), p0.Y(), p0.Z(), p1.X(), p1.Y(), p1.Z(), p2.X(), p2.Y(), p2.Z());
    }
    else if (shape->ClassID() == R3TriangleArray::CLASS_ID()) {
      R3TriangleArray *array = (R3TriangleArray *) shape;
      for (int j = 0; j < array->NTriangles(); j++) {
        R3Triangle *triangle = array->Triangle(j);
        R3TriangleVertex *v0 = (transformation.IsMirrored()) ? triangle->V2() : triangle->V0();
        R3TriangleVertex *v1 = triangle->V1();
        R3TriangleVertex *v2 = (transformation.IsMirrored()) ? triangle->V0() : triangle->V2();
        const R3Point& p0 = v0->Position();
        const R3Point& p1 = v1->Position();
        const R3Point& p2 = v2->Position();
        fprintf(fp, "%stri  %d   %g %g %g   %g %g %g   %g %g %g\n", indent, material_index,
          p0.X(), p0.Y(), p0.Z(), p1.X(), p1.Y(), p1.Z(), p2.X(), p2.Y(), p2.Z());
      }
    }
    else {
      fprintf(stderr, "Warning: unrecognized shape type %d\n", shape->ClassID());
      // return 0;
    }
  }
  
  // Return success
  return 1;
}



static int
WritePrincetonReference(const R3Scene *scene, R3SceneReference *reference,
  const R3Affine& transformation, FILE *fp, const char *indent)
{
  // Get referenced scene
  R3Scene *referenced_scene = reference->ReferencedScene();

  // Check if referenced scene has a file to reference
  if (!referenced_scene->Filename()) {
    // The materials will not be right if simply write referenced scene!!!
    // So WritePrincetonFile should have removed references 
    RNAbort("Should never get here");
    return 0;
  }

  // Write include statement
  fprintf(fp, "%sinclude %s\n", indent, referenced_scene->Filename());

  // Return success
  return 1;
}



static int
WritePrincetonNode(const R3Scene *scene, const R3SceneNode *node,
  const R3Transformation& parent_transformation,
  FILE *fp, const char *indent)
{
  // Compute indent
  int nindent = strlen(indent);
  char *child_indent = new char [nindent + 3];
  strcpy(child_indent, indent);
  child_indent[nindent] = ' ';
  child_indent[nindent+1] = ' ';
  child_indent[nindent+2] = '\0';

  // Compute transformation
  R3Affine transformation = R3identity_affine;
  transformation.Transform(parent_transformation);
  transformation.Transform(node->Transformation());

  // Don't write root separately, unless transformed
  RNBoolean print_node = (node != scene->Root()) || (!node->Transformation().IsIdentity());
  if (print_node) {
    // Write group start
    fprintf(fp, "\n");
    if (!node->Name()) fprintf(fp, "%sbegin -1\n", indent);
    else fprintf(fp, "%sgroup %s -1\n", indent, node->Name());

    // Write node transformation
    const R4Matrix& m = node->Transformation().Matrix();
    fprintf(fp, "%s%g %g %g %g\n", child_indent, m[0][0], m[0][1], m[0][2], m[0][3]);
    fprintf(fp, "%s%g %g %g %g\n", child_indent, m[1][0], m[1][1], m[1][2], m[1][3]);
    fprintf(fp, "%s%g %g %g %g\n", child_indent, m[2][0], m[2][1], m[2][2], m[2][3]);
    fprintf(fp, "%s%g %g %g %g\n", child_indent, m[3][0], m[3][1], m[3][2], m[3][3]);
  }
  
  // Write elements
  for (int i = 0; i < node->NElements(); i++) {
    R3SceneElement *element = node->Element(i);
    if (!WritePrincetonElement(scene, element, transformation, fp, child_indent)) return 0;
  }

  // Write references
  for (int i = 0; i < node->NReferences(); i++) {
    R3SceneReference *reference = node->Reference(i);
    if (!WritePrincetonReference(scene, reference, transformation, fp, child_indent)) return 0;
  }

  // Write children
  for (int i = 0; i < node->NChildren(); i++) {
    R3SceneNode *child = node->Child(i);
    if (!WritePrincetonNode(scene, child, transformation, fp, child_indent)) return 0;
  }

  // Write end
  if (print_node) fprintf(fp, "%send\n", indent);

  // Delete child indent buffer
  delete [] child_indent;

  // Return success
  return 1;
}



int R3Scene::
WritePrincetonFile(const char *filename) const
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "w"))) {
    fprintf(stderr, "Unable to open file %s", filename);
    return 0;
  }

  // Determine output directory
  const char *output_directory = ".";

  // Determine if need to remove references
  RNBoolean remove_references = FALSE;
  for (int i = 0; i < NReferencedScenes(); i++) {
    R3Scene *referenced_scene = ReferencedScene(i);
    if (referenced_scene->Filename()) continue;
    remove_references = TRUE; break;
  }

  // Remove references if necessary (this will change the scene structure :( )
  if (remove_references) ((R3Scene *) this)->RemoveReferences();

  // Write color and camera 
  fprintf(fp, "background %g %g %g\n", Background().R(), Background().G(), Background().B());
  fprintf(fp, "ambient %g %g %g\n", Ambient().R(), Ambient().G(), Ambient().B());
  fprintf(fp, "camera %g %g %g  %g %g %g  %g %g %g  %g   %g %g  \n",
    Camera().Origin().X(), Camera().Origin().Y(), Camera().Origin().Z(),
    Camera().Towards().X(), Camera().Towards().Y(), Camera().Towards().Z(),
    Camera().Up().X(), Camera().Up().Y(), Camera().Up().Z(),
    Camera().XFOV(), Camera().Near(), Camera().Far());

  // Create textures files
  for (int i = 0; i < NTextures(); i++) {
    R2Texture *texture = Texture(i);
    if (!CreateTextureFile(this, texture, output_directory)) return 0;
  }

  // Write lights
  fprintf(fp, "\n");
  for (int i = 0; i < NLights(); i++) {
    R3Light *light = Light(i);
    if (!WritePrincetonLight(this, light, fp)) return 0;
  }

  // Write materials
  fprintf(fp, "\n");
  for (int i = 0; i < NMaterials(); i++) {
    R3Material *material = Material(i);
    if (!WritePrincetonMaterial(this, material, fp)) return 0;
  }

  // Write nodes recursively
  int status = WritePrincetonNode(this, root, R3identity_affine, fp, "");
  
  // Close file
  fclose(fp);

  // Return status
  return status;
}



////////////////////////////////////////////////////////////////////////
// SUPPORT HIERARCHY FILE I/O FUNCTIONS
////////////////////////////////////////////////////////////////////////

int R3Scene::
ReadSupportHierarchyFile(const char *filename, R3SceneNode *parent_node)
{
  // Get/set parent node
  if (!parent_node) parent_node = root;
  
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "r"))) {
    fprintf(stderr, "Unable to open file %s", filename);
    return 0;
  }

  // Read scene
  char buffer[1024];
  int line_count = 0;
  RNArray<R3SceneNode *> nodes;
  while (fgets(buffer, 1023, fp)) {
    // Increment line counter
    line_count++;

    // Skip white space
    char *bufferp = buffer;
    while (isspace(*bufferp)) bufferp++;

    // Skip blank lines and comments
    if (*bufferp == '#') continue;
    if (*bufferp == '\0') continue;

    // Get keyword
    char keyword[256];
    if (sscanf(bufferp, "%s", keyword) != 1) {
      fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
      return 0;
    }

    // Check keyword
    if (!strcmp(keyword, "newModel")) {
      // Read fields
      int model_index;
      char model_name[1024];
      if (sscanf(bufferp, "%s%d%s", keyword, &model_index, model_name) != (unsigned int) 3) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Create node
      assert(nodes.NEntries() == model_index);
      R3SceneNode *node = new R3SceneNode(this);
      node->SetName(model_name);
      parent_node->InsertChild(node);
      nodes.Insert(node);

      // Read obj file
      char model_filename[1024];
      sprintf(model_filename, "models/%s.obj", model_name);
      if (!ReadObj(this, node, model_filename)) return 0;
    }
    else if (!strcmp(keyword, "parentIndex")) {
      // Read fields
      int parent_index;
      if (sscanf(bufferp, "%s%d", keyword, &parent_index) != (unsigned int) 2) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Check parent index
      if (parent_index < 0) {
        if (nodes.NEntries() != 1) {
          fprintf(stderr, "Root node was not first in file %s", filename);
          return 0;
        }
      }
      else {
        // Just checking
        if (parent_index >= nodes.NEntries()) {
          fprintf(stderr, "Invalid parent node index %d on line %d in file %s", parent_index, line_count, filename);
          return 0;
        }

        // Set last node's parent
        R3SceneNode *node = nodes.Tail();
        R3SceneNode *parent = nodes.Kth(parent_index);
        R3SceneNode *previous_parent = node->Parent();
        if (parent != previous_parent) {
          if (previous_parent) previous_parent->RemoveChild(node);
          if (parent) parent->InsertChild(node);
        }
      }
    }
    else if (!strcmp(keyword, "transform")) {
      // Read fields
      double m[16];
      if (sscanf(bufferp, "%s%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", 
        keyword, &m[0], &m[1], &m[2], &m[3], &m[4], &m[5], &m[6], &m[7], 
        &m[8], &m[9], &m[10], &m[11], &m[12], &m[13], &m[14], &m[15]) != (unsigned int) 17) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // NOTE: Our transformation is stored relative to parent's coordinate system
      // So, we must compute transformation from inverse of transform from parent node to world coordinates so that can 
      // convert file's absolute transform (which goes from node's coordinates to world coordinates)
      // to our transform (which goes from node's coordinates to parent node's coordinates)
      R3Affine transformation = R3identity_affine;
      R3SceneNode *node = nodes.Tail();
      R3SceneNode *ancestor = node->Parent();
      while (ancestor) {
        transformation.InverseTransform(ancestor->Transformation());
        ancestor = ancestor->Parent();
      }

      // Set last node's transformation
      // Note that file's matrix is for post-multiplication, while ours is for pre-multiplication, so need flip
      R4Matrix matrix(m); matrix.Flip();
      transformation.Transform(R3Affine(matrix, 0));
      node->SetTransformation(transformation);
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// GRAMMAR HIERARCHY FILE I/O FUNCTIONS
////////////////////////////////////////////////////////////////////////

int R3Scene::
ReadGrammarHierarchyFile(const char *filename, R3SceneNode *parent_node)
{
  // Get/set parent node
  if (!parent_node) parent_node = root;

  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "r"))) {
    fprintf(stderr, "Unable to open file %s", filename);
    return 0;
  }

  // Extract base name
  char basename[4096] = { '\0' };
  const char *startp = strrchr(filename, '/');
  if (!startp) startp = filename;
  else startp++;
  const char *endp = strrchr(filename, '.');
  if (!endp) endp = startp + strlen(startp);
  int basename_length = endp - startp;
  if (basename_length > 4095) basename_length = 4095;
  strncpy(basename, startp, basename_length);

  // Read scene
  char buffer[1024];
  int line_count = 0;
  RNBoolean leaf = TRUE;
  RNArray<R3SceneNode *> parsed_nodes;
  R3SceneNode *current_node = NULL;
  while (fgets(buffer, 1023, fp)) {
    // Increment line counter
    line_count++;

    // Skip white space
    char *bufferp = buffer;
    while (isspace(*bufferp)) bufferp++;

    // Skip blank lines and comments
    if (*bufferp == '#') continue;
    if (*bufferp == '\0') continue;

    // Get keyword
    char keyword[256];
    if (sscanf(bufferp, "%s", keyword) != 1) {
      fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
      return 0;
    }

    // Check keyword
    if (!strcmp(keyword, "newModel")) {
      // Read fields
      int model_index;
      if (sscanf(bufferp, "%s%d", keyword, &model_index) != (unsigned int) 2) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Create node
      while (parsed_nodes.NEntries() <= model_index) {
        char node_name[4096];
        sprintf(node_name, "%d", parsed_nodes.NEntries());
        R3SceneNode *node = new R3SceneNode(this);
        node->SetName(node_name);
        parsed_nodes.Insert(node);
      }

      // Remember current node
      current_node = parsed_nodes.Kth(model_index);
    }
    else if (!strcmp(keyword, "root")) {
      // Read fields
      int root_index;
      if (sscanf(bufferp, "%s%d", keyword, &root_index) != (unsigned int) 2) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Create node
      while (parsed_nodes.NEntries() <= root_index) {
        char node_name[4096];
        sprintf(node_name, "%d", parsed_nodes.NEntries());
        R3SceneNode *node = new R3SceneNode(this);
        node->SetName(node_name);
        parsed_nodes.Insert(node);
      }

      // Insert root of this parse as child of root
      R3SceneNode *node = parsed_nodes.Kth(root_index);
      parent_node->InsertChild(node);
    }
    else if (!strcmp(keyword, "parent")) {
      // Read fields
      int parent_index;
      if (sscanf(bufferp, "%s%d", keyword, &parent_index) != (unsigned int) 2) {
        fprintf(stderr, "Syntax error on line %d in file %s", line_count, filename);
        return 0;
      }

      // Check parent index
      if (parent_index >= 0) {
        // Create parent node
        while (parsed_nodes.NEntries() <= parent_index) {
          char node_name[4096];
          sprintf(node_name, "I%d", parsed_nodes.NEntries());
          R3SceneNode *node = new R3SceneNode(this);
          node->SetName(node_name);
          parsed_nodes.Insert(node);
        }
      
        // Set last node's parent
        R3SceneNode *parent = parsed_nodes.Kth(parent_index);
        parent->InsertChild(current_node);
      }
    }
    else if (!strcmp(keyword, "children")) {
      const char *token = strtok(bufferp, " \n\t");
      assert(token && (!strcmp(token, "children")));
      token = strtok(NULL, " \n\t");
      leaf = (token) ? FALSE : TRUE;
    }
    else if (!strcmp(keyword, "leaf_group")) {
      // Check current node
      if (!current_node) {
        fprintf(stderr, "leaf_group before first newModel at line %d in %s\n", line_count, filename);
        return 0;
      }

      // Read models
      if (leaf) {
        const char *token = strtok(bufferp, " \n\t");
        assert(token && (!strcmp(token, "leaf_group")));
        while (TRUE) {
          token = strtok(NULL, " \n\t");
          if (!token) break;
          int model_index = atoi(token);
          // char node_name[4096];
          // sprintf(node_name, "L%d", model_index);
          // R3SceneNode *node = new R3SceneNode(this);
          // node->SetName(node_name);
          // current_node->InsertChild(node);
          char model_filename[1024];
          sprintf(model_filename, "models/%s/%d.off", basename, model_index);
          R3TriangleArray *shape = ReadMesh(model_filename);
          if (shape) {
            R3SceneElement *element = new R3SceneElement();
            element->SetMaterial(&R3default_material);
            element->InsertShape(shape);
            current_node->InsertElement(element);
          }
        }    
      }    
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// PARSE FILE I/O FUNCTIONS
////////////////////////////////////////////////////////////////////////

int R3Scene::
ReadParseFile(const char *filename, R3SceneNode *parent_node)
{
  // Get/set parent node
  if (!parent_node) parent_node = root;

  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "r"))) {
    fprintf(stderr, "Unable to open file %s", filename);
    return 0;
  }

  // Read header
  char buffer[4096];
  if (!fgets(buffer, 4096, fp)) {
    fprintf(stderr, "Unable to read object parse file %s\n", filename);
    return 0;
  }

  // Check header
  if (strncmp(buffer, "OBJECT PARSE 1.0", 16)) {
    fprintf(stderr, "Error in header of oject parse file %s\n", filename);
    return 0;
  }

  // Read file
  int line_number = 0;
  int assignment_index = 0;
  RNArray<R3Shape *> shapes;
  char mesh_directory[4096] = { '.', '\0' };
  while (fgets(buffer, 4096, fp)) {
    // Check line
    line_number++;
    char *bufferp = buffer;
    while (*bufferp && isspace(*bufferp)) bufferp++;
    if (!bufferp) continue;
    if (*bufferp == '#') continue;

    // Parse line
    char keyword[4096];
    if (sscanf(bufferp, "%s", keyword) == (unsigned int) 1) {
      if (!strcmp(keyword, "A")) {
        // Parse assignment
        double score, m[16];
        int segmentation_index, model_index, dummy;
        if (sscanf(bufferp, "%s%d%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%d%d%d%d", keyword, 
          &segmentation_index, &model_index, 
          &m[0], &m[1], &m[2], &m[3], &m[4], &m[5], &m[6], &m[7], 
          &m[8], &m[9], &m[10], &m[11], &m[12], &m[13], &m[14], &m[15], 
          &score, &dummy, &dummy, &dummy, &dummy) != (unsigned int) 24) {
          fprintf(stderr, "Error parsing assignment at line %d of %s\n", line_number, filename);
          return 0;
        }

        // Create node
        R3SceneNode *node = new R3SceneNode(this);
        parent_node->InsertChild(node);

        // Create shape element
        R3Shape *shape = shapes.Kth(model_index);
        R3Brdf *brdf = new R3Brdf(RNRgb(0, 0.25 + score, 0), 0.0, 0.25 + score);
        InsertBrdf(brdf);
        R3Material *material = new R3Material(brdf);
        InsertMaterial(material);
        R3SceneElement *element = new R3SceneElement();
        element->SetMaterial(material);
        element->InsertShape(shape);
        node->InsertElement(element);

        // Set node name
        char node_name[1024];
        sprintf(node_name, "A%d_M%d_S%03d", assignment_index++, model_index, (int) (1000*score));
        node->SetName(node_name);

        // Set node transformation
        R4Matrix matrix(m); 
        R3Affine affine(matrix, 0);
        node->SetTransformation(affine);
      }
      else if (!strcmp(keyword, "M")) {
        // Parse model
        int dummy;
        double cx, cy, cz, r, h;
        char model_name[4096], mesh_name[4096];
        if (sscanf(bufferp, "%s%d%lf%lf%lf%lf%lf%s%s%d%d%d%d%d", keyword, 
          &dummy, &cx, &cy, &cz, &r, &h, model_name, mesh_name, &dummy, &dummy, &dummy, &dummy, &dummy) != (unsigned int) 14) {
          fprintf(stderr, "Error parsing model at line %d of %s\n", line_number, filename);
          return 0;
        }

        // Read mesh
        char mesh_filename[4096];
        R3Shape *shape = NULL;
        if (strcmp(mesh_name, "None")) {
          sprintf(mesh_filename, "%s/%s", mesh_directory, mesh_name);
          shape = ReadMesh(mesh_filename);
          if (!shape) return 0;
        }
        else {
          shape = new R3Sphere(R3Point(0,0,0), 0);
        }

        // Insert shape into list
        shapes.Insert(shape);
      }
      else if (!strcmp(keyword, "MD")) {
        // Parse model directory
        if (sscanf(bufferp, "%s%s", keyword, mesh_directory) != (unsigned int) 2) {
          fprintf(stderr, "Error parsing model directory at line %d of %s\n", line_number, filename);
          return 0;
        }
      }
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// RECTANGLE FILE I/O FUNCTIONS
////////////////////////////////////////////////////////////////////////

int R3Scene::
ReadRectangleFile(const char *filename, R3SceneNode *parent_node)
{
  // Get/check parent node
  if (!parent_node) parent_node = root;

  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "r"))) {
    fprintf(stderr, "Unable to open file %s", filename);
    return 0;
  }

  // Read file
  char buffer[4096];
  int line_number = 0;
  int assignment_index = 0;
  while (fgets(buffer, 4096, fp)) {
    // Check line
    line_number++;
    char *bufferp = buffer;
    while (*bufferp && isspace(*bufferp)) bufferp++;
    if (!bufferp) continue;
    if (*bufferp == '#') continue;

    // Parse line
    char name[4096];
    double x1, y1, x2, y2, x3, y3, x4, y4, zmin, zmax, score;
    if (sscanf(bufferp, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%s", 
       &x1, &y1, &x2, &y2, &x3, &y3, &x4, &y4, &zmin, &zmax, &score, name) != (unsigned int) 12) {
      fprintf(stderr, "Error parsing line %d of %s\n", line_number, filename);
      return 0;
    }

    // Create node
    R3SceneNode *node = new R3SceneNode(this);
    parent_node->InsertChild(node);

    // Create shape element
    RNArray<R3Triangle *> triangles;
    RNArray<R3TriangleVertex *> vertices;
    R3TriangleVertex *v1 = new R3TriangleVertex(R3Point(x1, y1, zmin)); vertices.Insert(v1);
    R3TriangleVertex *v2 = new R3TriangleVertex(R3Point(x2, y2, zmin)); vertices.Insert(v2);
    R3TriangleVertex *v3 = new R3TriangleVertex(R3Point(x3, y3, zmin)); vertices.Insert(v3);
    R3TriangleVertex *v4 = new R3TriangleVertex(R3Point(x4, y4, zmin)); vertices.Insert(v4);
    R3Triangle *t1 = new R3Triangle(v1, v2, v3); triangles.Insert(t1);
    R3Triangle *t2 = new R3Triangle(v1, v3, v4); triangles.Insert(t2);
    R3TriangleArray *base = new R3TriangleArray(vertices, triangles);
    R3Point base_centroid = base->BBox().Centroid();
    R3Point top_centroid = base_centroid + (zmax - zmin) * R3posz_vector;
    R3Cylinder *marker = new R3Cylinder(base_centroid, top_centroid, 0.01 * base->BBox().DiagonalRadius());
    R3Brdf *brdf = new R3Brdf(RNRgb(0, 0.25 + score, 0));
    InsertBrdf(brdf);
    R3Material *material = new R3Material(brdf);
    InsertMaterial(material);
    R3SceneElement *element = new R3SceneElement();
    element->SetMaterial(material);
    element->InsertShape(base);
    element->InsertShape(marker);
    node->InsertElement(element);
    
    // Set node name
    char node_name[1024];
    sprintf(node_name, "A%d_S%03d", assignment_index++, (int) (1000*score));
    node->SetName(node_name);
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}




////////////////////////////////////////////////////////////////////////
// SUNCG PARSING FUNCTIONS
////////////////////////////////////////////////////////////////////////

static int
GetJsonObjectMember(Json::Value *&result, Json::Value *object, const char *str, int expected_type = 0)
{
  // Check object type
  if (object->type() != Json::objectValue) {
    // fprintf(stderr, "JSON: not an object\n");
    return 0;
  }

  // Check object member
  if (!object->isMember(str)) {
    // fprintf(stderr, "JSON object has no member named %s\n", str);
    return 0;
  }

  // Get object member
  result = &((*object)[str]);
  if (result->type() == Json::nullValue) {
    // fprintf(stderr, "JSON object has null member named %s\n", str);
    return 0;
  }

  // Check member type
  if (expected_type > 0) {
    if (result->type() != expected_type) {
      // fprintf(stderr, "JSON object member %s has unexpected type %d (rather than %d)\n", str, result->type(), expected_type);
      return 0;
    }
  }
  
  // Check for empty strings
  if (result->type() == Json::stringValue) {
    if (result->asString().length() == 0) {
      // fprintf(stderr, "JSON object has zero length string named %s\n", str);
      return 0;
    }
  }

  // Return success
  return 1;
}



static int
GetJsonArrayEntry(Json::Value *&result, Json::Value *array, unsigned int k, int expected_type = -1)
{
  // Check array type
  if (array->type() != Json::arrayValue) {
    fprintf(stderr, "JSON: not an array\n");
    return 0;
  }

  // Check array size
  if (array->size() <= k) {
    // fprintf(stderr, "JSON array has no member %d\n", k);
    return 0;
  }

  // Get entry
  result = &((*array)[k]);
  if (result->type() == Json::nullValue) {
    // fprintf(stderr, "JSON array has null member %d\n", k);
    return 0;
  }

  // Check entry type
  if (expected_type > 0) {
    if (result->type() != expected_type) {
      // fprintf(stderr, "JSON array entry %d has unexpected type %d (rather than %d)\n", k, result->type(), expected_type);
      return 0;
    }
  }
  
  // Return success
  return 1;
}



static int
ParseSUNCGMaterials(R3Scene *scene,
  RNSymbolTable<R2Texture *>& texture_symbol_table,
  const RNArray<R3Material *>& input_materials,
  RNArray<R3Material *>& output_materials,
  Json::Value *json_materials)
{
  // Parse JSON array of materials
  Json::Value *json_material, *json_value;
  for (Json::ArrayIndex index = 0; index < json_materials->size(); index++) {
    char material_name[1024] = { '\0' };
    char texture_name[1024] = { '\0' };
    char diffuse_string[1024] = { '\0' };
    if (!GetJsonArrayEntry(json_material, json_materials, index)) continue; 
    if (json_material->type() != Json::objectValue) continue;
    if (GetJsonObjectMember(json_value, json_material, "name"))
      strncpy(material_name, json_value->asString().c_str(), 1024);
    if (GetJsonObjectMember(json_value, json_material, "texture")) 
      strncpy(texture_name, json_value->asString().c_str(), 1024);
    if (GetJsonObjectMember(json_value, json_material, "diffuse")) 
      strncpy(diffuse_string, json_value->asString().c_str(), 1024);

    // Get input material
    R3Material *input_material = (input_materials.NEntries() > (int) index) ? input_materials.Kth(index) : NULL;
    const R3Brdf *input_brdf = (input_material) ? input_material->Brdf() : NULL;
    const R2Texture *input_texture = (input_material) ? input_material->Texture() : NULL;
    
    // Get/create output material
    R3Material *output_material = NULL;
    if (*diffuse_string || *texture_name) {
      // Create output brdf
      R3Brdf *output_brdf = NULL;
      if (input_brdf) output_brdf = new R3Brdf(*input_brdf);
      else { output_brdf = new R3Brdf(); output_brdf->SetDiffuse(RNwhite_rgb); }
      scene->InsertBrdf(output_brdf);
      if (*diffuse_string) {
        long int b = strtol(&diffuse_string[5], NULL, 16); diffuse_string[5] = '\0';
        long int g = strtol(&diffuse_string[3], NULL, 16); diffuse_string[3] = '\0';
        long int r = strtol(&diffuse_string[1], NULL, 16); diffuse_string[1] = '\0';
        RNRgb diffuse_rgb(r / 255.0, g / 255.0, b / 255.0);
        output_brdf->SetDiffuse(diffuse_rgb);
      }

      // Create output texture
      R2Texture *output_texture = NULL;
      if (*texture_name) {
        char texture_filename[1024];
        const char *texture_directory = "../../texture";
        sprintf(texture_filename, "%s/%s.png", texture_directory, texture_name);
        if (!RNFileExists(texture_filename)) sprintf(texture_filename, "%s/%s.jpg", texture_directory, texture_name);
        if (!texture_symbol_table.Find(texture_filename, &output_texture)) {
          if (input_texture) output_texture = new R2Texture(*input_texture);
          else output_texture = new R2Texture();
          R2Image *image = new R2Image();
          if (!image->Read(texture_filename)) return 0;
          output_texture->SetImage(image);
          output_texture->SetFilename(texture_filename);
          output_texture->SetName(texture_name);
          scene->InsertTexture(output_texture);
          texture_symbol_table.Insert(texture_filename, output_texture);
        }
      }
      else if (input_texture) {
        const char *texture_filename = input_texture->Filename();
        if (!texture_symbol_table.Find(texture_filename, &output_texture)) {
          output_texture = new R2Texture(*input_texture);
          scene->InsertTexture(output_texture);
          texture_symbol_table.Insert(texture_filename, output_texture);
        }
      }

      // Create output material
      output_material = new R3Material(output_brdf, output_texture);
      output_material->SetName(material_name);
      scene->InsertMaterial(output_material);
    }

    // Insert material into array of results
    output_materials.Insert(output_material);
  }

  // Return success
  return 1;
}



static int
CreateBox(R3Scene *scene, R3SceneNode *node,
  RNScalar dimensions[3],
  const RNArray<R3Material *>& materials)
{
  // Create six sides of box
  for (int dir = 0; dir < 2; dir++) {
    for (int dim = 0; dim < 3; dim++) {
      int dim1 = (dim+1)%3;
      int dim2 = (dim+2)%3;

      // Compute radii
      RNScalar d0 = dimensions[dim];
      RNScalar d1 = dimensions[dim1];
      RNScalar d2 = dimensions[dim2];
      RNScalar r0 = (dir == 0) ? -0.5*d0 :  0.5*d0;
      RNScalar r1 = (dir == 0) ? -0.5*d1 : 0.5*d1;
      RNScalar r2 = 0.5*d2;

      // Compute axes
      R3Vector axis0 = r0 * R3xyz_triad.Axis(dim);
      R3Vector axis1 = r1 * R3xyz_triad.Axis(dim1);
      R3Vector axis2 = r2 * R3xyz_triad.Axis(dim2);

      // Compute vertex positions
      R3Point center(0, 0, 0);
      R3Point p00 = center + axis0 - axis1 - axis2;
      R3Point p10 = center + axis0 + axis1 - axis2;
      R3Point p11 = center + axis0 + axis1 + axis2;
      R3Point p01 = center + axis0 - axis1 + axis2;

      // Compute texture coordinates
      R2Point t00(0, 0);
      R2Point t10(d1, 0);
      R2Point t11(d1, d2);
      R2Point t01(0, d2);

      // Create two triangles
      RNArray<R3TriangleVertex *> vertices;
      R3TriangleVertex *v00 = new R3TriangleVertex(p00, t00); vertices.Insert(v00);
      R3TriangleVertex *v10 = new R3TriangleVertex(p10, t10); vertices.Insert(v10);
      R3TriangleVertex *v11 = new R3TriangleVertex(p11, t11); vertices.Insert(v11);
      R3TriangleVertex *v01 = new R3TriangleVertex(p01, t01); vertices.Insert(v01);
      RNArray<R3Triangle *> triangles;
      R3Triangle *tri0 = new R3Triangle(v00, v10, v11); triangles.Insert(tri0);
      R3Triangle *tri1 = new R3Triangle(v00, v11, v01); triangles.Insert(tri1);
      R3TriangleArray *shape = new R3TriangleArray(vertices, triangles);

      // Create scene element
      static const int material_indices[2][3] = { { 2, 1, 5 }, { 3, 0, 4 } };
      int material_index = material_indices[dir][dim];
      R3Material *material = (materials.NEntries() > material_index) ? materials[material_index] : &R3default_material;
      R3SceneElement *element = new R3SceneElement(material);
      element->InsertShape(shape);
      node->InsertElement(element);
    }
  }

  // Return success
  return 1;
}



int R3Scene::
ReadSUNCGFile(const char *filename, R3SceneNode *parent_node)
{
  // Useful variables
  const char *input_data_directory = "../..";
  RNSymbolTable<R2Texture *> texture_symbol_table;
  RNSymbolTable<R3Scene *> model_symbol_table;
  if (!parent_node) parent_node = root;
  
  // Open file
  FILE* fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open SUNCG file %s\n", filename);
    return 0;
  }

  // Read file 
  std::string text;
  fseek(fp, 0, SEEK_END);
  long const size = ftell(fp);
  fseek(fp, 0, SEEK_SET);
  char* buffer = new char[size + 1];
  unsigned long const usize = static_cast<unsigned long const>(size);
  if (fread(buffer, 1, usize, fp) != usize) { fprintf(stderr, "Unable to read %s\n", filename); return 0; }
  else { buffer[size] = 0; text = buffer; }
  delete[] buffer;

  // Close file
  fclose(fp);

  // Digest file
  Json::Value json_root;
  Json::Reader json_reader;
  Json::Value *json_items, *json_item, *json_value;
  if (!json_reader.parse(text, json_root, false)) {
    fprintf(stderr, "Unable to parse %s\n", filename);
    return 0;
  }

  // Get/check version
  char version[1024];
  strncpy(version, "suncg@1.0.0", 1024);
  if (GetJsonObjectMember(json_value, &json_root, "version", Json::stringValue)) {
    strncpy(version, json_value->asString().c_str(), 1024);
    if (strcmp(version, "suncg@1.0.0")) {
      fprintf(stderr, "Unrecognized version %s in SUNCG file %s\n", version, filename);
      return 0;
    }
  }
  
  // Get scene id
  char scene_id[1024];
  strncpy(scene_id, "NoName", 1024);
  if (GetJsonObjectMember(json_value, &json_root, "id", Json::stringValue)) {
    strncpy(scene_id, json_value->asString().c_str(), 1024);
  }
  
  // Get scene up direction
  R3Vector scene_up(0, 1, 0);
  if (GetJsonObjectMember(json_items, &json_root, "up", Json::arrayValue)) {
    if (json_items->size() >= 3) {
      if (GetJsonArrayEntry(json_item, json_items, 0))
        scene_up[0] = json_item->asDouble();
      if (GetJsonArrayEntry(json_item, json_items, 1))
        scene_up[1] = json_item->asDouble();
      if (GetJsonArrayEntry(json_item, json_items, 2))
        scene_up[2] = json_item->asDouble();
      scene_up.Normalize();
    }
  }

  // Get scene front direction
  R3Vector scene_front(0, 0, 1);
  if (GetJsonObjectMember(json_items, &json_root, "front", Json::arrayValue)) {
    if (json_items->size() >= 3) {
      if (GetJsonArrayEntry(json_item, json_items, 0))
        scene_front[0] = json_item->asDouble();
      if (GetJsonArrayEntry(json_item, json_items, 1))
        scene_front[1] = json_item->asDouble();
      if (GetJsonArrayEntry(json_item, json_items, 2))
        scene_front[2] = json_item->asDouble();
      scene_front.Normalize();
    }
  }

  // Get scene scale factor (to convert to meters)
  double scaleToMeters = 1.0;
  if (GetJsonObjectMember(json_value, &json_root, "scaleToMeters")) {
    scaleToMeters = json_value->asDouble();
  }
  
  // Create scene node 
  R3SceneNode *scene_node = parent_node;
  scene_node->SetName(scene_id);

  // Set scene transformation
  R3Affine scene_transformation = R3identity_affine;
  // Set up and front directions
  scene_transformation.Scale(scaleToMeters);
  scene_node->SetTransformation(scene_transformation);

  // Parse levels
  Json::Value *json_levels, *json_level;
  if (!GetJsonObjectMember(json_levels, &json_root, "levels", Json::arrayValue)) return 0;
  for (Json::ArrayIndex index = 0; index < json_levels->size(); index++) {
    if (!GetJsonArrayEntry(json_level, json_levels, index)) return 0;
    if (json_level->type() != Json::objectValue) continue;
           
    // Parse level attributes
    int level_id = index;
    if (GetJsonObjectMember(json_value, json_level, "valid"))
      if (!json_value->asString().compare(std::string("0")))  continue;
    if (GetJsonObjectMember(json_value, json_level, "id"))
      level_id = atoi(json_value->asString().c_str());

    // Create level node
    char level_name[1024];
    sprintf(level_name, "Level#%d", level_id);
    R3SceneNode *level_node = new R3SceneNode(this);
    level_node->SetName(level_name);
    scene_node->InsertChild(level_node);

    // Parse nodes
    Json::Value *json_nodes, *json_node, *json_materials;
    if (GetJsonObjectMember(json_nodes, json_level, "nodes", Json::arrayValue)) {
      if (json_nodes->size() == 0) continue;
      R3SceneNode **created_nodes = new R3SceneNode * [ json_nodes->size() ];
      for (unsigned int i = 0; i < json_nodes->size(); i++) created_nodes[i] = NULL;
      for (Json::ArrayIndex index = 0; index < json_nodes->size(); index++) {
        if (!GetJsonArrayEntry(json_node, json_nodes, index)) continue; 
        if (json_node->type() != Json::objectValue) continue;
 
        // Parse node attributes
        char node_id[1024] = { '\0' };;
        char modelId[1024] = { '\0' };;
        char node_type[1024] = { '\0' };
        int hideCeiling = 0, hideFloor = 0, hideWalls = 0;
        int isMirrored = 0, state = 0;
        if (GetJsonObjectMember(json_value, json_node, "valid"))
          if (!json_value->asString().compare(std::string("0")))  continue;
        if (GetJsonObjectMember(json_value, json_node, "id"))
          strncpy(node_id, json_value->asString().c_str(), 1024);
        if (GetJsonObjectMember(json_value, json_node, "type")) 
          strncpy(node_type, json_value->asString().c_str(), 1024);
        if (GetJsonObjectMember(json_value, json_node, "modelId"))
          strncpy(modelId, json_value->asString().c_str(), 1024);
        if (GetJsonObjectMember(json_value, json_node, "hideCeiling")) 
          if (!json_value->asString().compare(std::string("1"))) hideCeiling = 1;
        if (GetJsonObjectMember(json_value, json_node, "hideFloor")) 
          if (!json_value->asString().compare(std::string("1"))) hideFloor = 1;
        if (GetJsonObjectMember(json_value, json_node, "hideWalls")) 
          if (!json_value->asString().compare(std::string("1"))) hideWalls = 1;
        if (GetJsonObjectMember(json_value, json_node, "isMirrored")) 
          if (!json_value->asString().compare(std::string("1"))) isMirrored = 1;
        if (GetJsonObjectMember(json_value, json_node, "state")) 
          if (!json_value->asString().compare(std::string("1"))) state = 1;

        // Parse node transformation
        R3Affine transformation = R3identity_affine;
        if (GetJsonObjectMember(json_items, json_node, "transform", Json::arrayValue)) {
          if (json_items->size() >= 16) {
            R4Matrix matrix = R4identity_matrix;
            for (Json::ArrayIndex index = 0; index < json_items->size(); index++) {
              if (!GetJsonArrayEntry(json_item, json_items, index)) continue;
              matrix[index%4][index/4] = json_item->asDouble();
            }
            transformation.Reset(matrix, isMirrored);
          }
        }

        // Create scene node(s) based on type
        char obj_name[4096], node_name[4096];
        if (!strcmp(node_type, "Ground")) {
          // Create node for ground
          sprintf(obj_name, "%s/room/%s/%sf.obj", input_data_directory, scene_id, modelId); 
          if (!hideFloor && RNFileExists(obj_name)) {
            R3Scene *model = new R3Scene();
            if (!ReadObj(model, model->Root(), obj_name)) return 0;
            sprintf(node_name, "Ground#%s", node_id);
            model->SetName(node_name);
            model->Root()->SetName(node_name);
            model->SetFilename(obj_name);
            InsertReferencedScene(model);
            R3SceneNode *node = new R3SceneNode(this);
            node->InsertReference(new R3SceneReference(model, materials));
            node->SetName(node_name);
            node->SetTransformation(transformation);
            level_node->InsertChild(node);
            // R3SceneNode *node = new R3SceneNode(this);
            // sprintf(node_name, "Ground#%s", node_id);
            // node->SetName(node_name);
            // if (!ReadObj(this, node, obj_name)) return 0;
            // level_node->InsertChild(node);
            created_nodes[index] = node;
          }
        }
        else if (!strcmp(node_type, "Room")) {
          // Create room node
          R3SceneNode *room_node = new R3SceneNode(this);
          sprintf(node_name, "Room#%s", node_id);
          room_node->SetName(node_name);
          // room_node->SetTransformation(transformation);
          level_node->InsertChild(room_node);
          created_nodes[index] = room_node;
          RNArray<R3Material *> materials;
          
           // Create node for floor
          sprintf(obj_name, "%s/room/%s/%sf.obj", input_data_directory, scene_id, modelId); 
          if (!hideFloor && RNFileExists(obj_name)) {
            R3Scene *model = new R3Scene();
            if (!ReadObj(model, model->Root(), obj_name)) return 0;
            sprintf(node_name, "Floor#%s", node_id);
            model->SetName(node_name);
            model->Root()->SetName(node_name);
            model->SetFilename(obj_name);
            InsertReferencedScene(model);
            R3SceneNode *node = new R3SceneNode(this);
            node->InsertReference(new R3SceneReference(model, materials));
            node->SetName(node_name);
            room_node->InsertChild(node);
            // R3SceneNode *node = new R3SceneNode(this);
            // sprintf(node_name, "Floor#%s", node_id);
            // node->SetName(node_name);
            // if (!ReadObj(this, node, obj_name)) return 0;
            // room_node->InsertChild(node);
          }

          // Create node for ceiling
          sprintf(obj_name, "%s/room/%s/%sc.obj", input_data_directory, scene_id, modelId); 
          if (!hideCeiling && RNFileExists(obj_name)) {
            R3Scene *model = new R3Scene();
            if (!ReadObj(model, model->Root(), obj_name)) return 0;
            sprintf(node_name, "Ceiling#%s", node_id);
            model->SetName(node_name);
            model->Root()->SetName(node_name);
            model->SetFilename(obj_name);
            InsertReferencedScene(model);
            R3SceneNode *node = new R3SceneNode(this);
            node->InsertReference(new R3SceneReference(model, materials));
            node->SetName(node_name);
            room_node->InsertChild(node);
            // R3SceneNode *node = new R3SceneNode(this);
            // sprintf(node_name, "Ceiling#%s", node_id);
            // node->SetName(node_name);
            // if (!ReadObj(this, node, obj_name)) return 0;
            // room_node->InsertChild(node);
          }

          // Create node for walls
          sprintf(obj_name, "%s/room/%s/%sw.obj", input_data_directory, scene_id, modelId); 
          if (!hideWalls && RNFileExists(obj_name)) {
            R3Scene *model = new R3Scene();
            if (!ReadObj(model, model->Root(), obj_name)) return 0;
            sprintf(node_name, "Wall#%s", node_id);
            model->SetName(node_name);
            model->Root()->SetName(node_name);
            model->SetFilename(obj_name);
            InsertReferencedScene(model);
            R3SceneNode *node = new R3SceneNode(this);
            node->InsertReference(new R3SceneReference(model, materials));
            node->SetName(node_name);
            room_node->InsertChild(node);
            // R3SceneNode *node = new R3SceneNode(this);
            // sprintf(node_name, "Wall#%s", node_id);
            // node->SetName(node_name);
            // if (!ReadObj(this, node, obj_name)) return 0;
            // room_node->InsertChild(node);
          }        
        }
        else if (!strcmp(node_type, "Object")) {
          // Read/get model 
          R3Scene *model = NULL;
          if (state) sprintf(obj_name, "%s/object/%s/%s_0.obj", input_data_directory, modelId, modelId); 
          else sprintf(obj_name, "%s/object/%s/%s.obj", input_data_directory, modelId, modelId); 
          if (!model_symbol_table.Find(obj_name, &model)) {
            model = new R3Scene();
            if (!ReadObj(model, model->Root(), obj_name)) return 0;
            sprintf(node_name, "Model#%s", modelId);
            model->SetName(modelId);
            model->Root()->SetName(node_name);
            model->SetFilename(obj_name);
            InsertReferencedScene(model);
            model_symbol_table.Insert(obj_name, model);
          }

          // Read materials
          RNArray<R3Material *> materials;
          if (GetJsonObjectMember(json_materials, json_node, "materials", Json::arrayValue)) {
            if (!ParseSUNCGMaterials(this, texture_symbol_table, model->materials, materials, json_materials)) return 0;
          }

          // Create node with reference to model
          R3SceneNode *node = new R3SceneNode(this);
          sprintf(node_name, "Object#%s", node_id);
          node->InsertReference(new R3SceneReference(model, materials));
          node->SetName(node_name);
          node->SetTransformation(transformation);
          level_node->InsertChild(node);
          created_nodes[index] = node;
        }
        else if (!strcmp(node_type, "Box")) {
          // Parse box dimensions
          RNScalar box_dimensions[3] = { 1, 1, 1 };
          if (GetJsonObjectMember(json_items, json_node, "dimensions", Json::arrayValue)) {
            if (json_items->size() >= 3) {
              if (GetJsonArrayEntry(json_item, json_items, 0))
                box_dimensions[0] = json_item->asDouble();
              if (GetJsonArrayEntry(json_item, json_items, 1))
                box_dimensions[1] = json_item->asDouble();
              if (GetJsonArrayEntry(json_item, json_items, 2))
                box_dimensions[2] = json_item->asDouble();
            }
          }

          // Read materials
          RNArray<R3Material *> materials;
          if (GetJsonObjectMember(json_materials, json_node, "materials", Json::arrayValue)) {
            if (!ParseSUNCGMaterials(this, texture_symbol_table, materials, materials, json_materials)) return 0;
          }

          // Create node for box
          R3SceneNode *node = new R3SceneNode(this);
          sprintf(node_name, "Box#%s", node_id);
          node->SetName(node_name);
          node->SetTransformation(transformation);
          if (!CreateBox(this, node, box_dimensions, materials)) return 0;
          level_node->InsertChild(node);
          created_nodes[index] = node;
        }
      }

      // Move created nodes to be children of room nodes
      for (Json::ArrayIndex index = 0; index < json_nodes->size(); index++) {
        if (!GetJsonArrayEntry(json_node, json_nodes, index)) continue; 
        if (json_node->type() != Json::objectValue) continue;
        if (!GetJsonObjectMember(json_value, json_node, "type")) continue;
        if (strcmp(json_value->asString().c_str(), "Room")) continue;
        R3SceneNode *room_node = created_nodes[index];
        if (!room_node) continue;
        if (GetJsonObjectMember(json_items, json_node, "nodeIndices", Json::arrayValue)) {
          for (Json::ArrayIndex room_index = 0; room_index < json_items->size(); room_index++) {
            GetJsonArrayEntry(json_item, json_items, room_index);
            if (json_item->isNumeric()) {
              int node_index = json_item->asInt();
              if ((node_index >= 0) && ((unsigned int) node_index < json_nodes->size())) {
                R3SceneNode *node = created_nodes[node_index];
                if (node) {
                  level_node->RemoveChild(node);
                  room_node->InsertChild(node); 
                }
              }
            }
          }
        }
      }
      
      // Delete array of created nodes
      delete [] created_nodes;
    }
  }
    
  // Return success
  return 1;
}



//////

static int
WriteSUNCGNode(const R3Scene *scene, const R3SceneNode *node, FILE *fp)
{
  // Check node
  if (!node->Name()) return 0;

  // NOTE: THE SCENE GRAPH MUST FOLLOW THE SUNCG HIERARCHY CONVENTIONS
  // ROOT -> LEVEL* -> OBJECT* | BOX* | GROUND* | (ROOM -> OBJECT*)*
  
  // Check node name
  static int counter = 0;
  if (!strncmp(node->Name(), "Level#", 6)) {
    // Write header
    counter = 0;
    const char *id = &(node->Name()[6]);
    const R3Box& b = node->WorldBBox();
    fprintf(fp, "  {\n");
    fprintf(fp, "    \"id\":\"%s\",\n", id);
    fprintf(fp, "    \"valid\":1,\n");
    fprintf(fp, "    \"bbox\":{\"min\":[%g,%g,%g],\"max\":[%g,%g,%g]},\n", b[0][0], b[0][1], b[0][2], b[1][0], b[1][1], b[1][2]);

    // Write children
    fprintf(fp, "    \"nodes\":[\n");
    for (int i = 0; i < node->NChildren(); i++) {
      R3SceneNode *child = node->Child(i);
      if (!WriteSUNCGNode(scene, child, fp)) return 0;
    }
    fprintf(fp, "    ]\n");
    
    // Write trailer
    fprintf(fp, "  }");
    if (node->ParentIndex() < node->Parent()->NChildren()-1) fprintf(fp, ",");
    fprintf(fp, "\n");
  }
  else if (!strncmp(node->Name(), "Ground#", 7) || !strncmp(node->Name(), "Room#", 5)) {
    // Get id and type
    char nodetype[4096];
    strncpy(nodetype, node->Name(), 4096);
    char *id = strchr(nodetype, '#');
    if (!id) return 0;
    else *(id++) = '\0';

    // Get modelId
    char modelId[4096], buffer[4096];
    strncpy(buffer, id, 4096);
    char *level = strtok(buffer, "_\n");
    if (!level) return 0;
    char *idx = strtok(NULL, "\n");
    if (!idx) return 0;
    sprintf(modelId, "fr_%srm_%s", level, idx);

    // Write json
    const R3Box& b = node->WorldBBox();
    if (counter++ > 0) fprintf(fp, ",\n");
    fprintf(fp, "      {\n");
    fprintf(fp, "        \"id\":\"%s\",\n", id);
    fprintf(fp, "        \"type\":\"%s\",\n", nodetype);
    fprintf(fp, "        \"modelId\":\"%s\",\n", modelId);
    fprintf(fp, "        \"valid\":1,\n");
    fprintf(fp, "        \"bbox\":{\"min\":[%g,%g,%g],\"max\":[%g,%g,%g]}\n",
      b[0][0], b[0][1], b[0][2], b[1][0], b[1][1], b[1][2]);

#if 0
    // Write node indices
    int first = 1;
    fprintf(fp, "        \"nodeIndices\":[");
    for (int i = 0; i < node->NChildren(); i++) {
      R3SceneNode *child = node->Child(i);
      if (!child->Name()) continue;
      if (strncmp(child->Name(), "Object#", 7)) continue;
      const char *s = strrchr(child->Name(), '_');
      if (!s) continue;
      if (!first) fprintf(fp, ","); 
      fprintf(fp, "%s", xxx <- needs to be index of the json entry in the level);
      first = 0;       
    }
    fprintf(fp, "]\n");
#endif
    
    fprintf(fp, "      }");

    // Write children
    for (int i = 0; i < node->NChildren(); i++) {
      R3SceneNode *child = node->Child(i);
      if (!WriteSUNCGNode(scene, child, fp)) return 0;
    }
  }
  else if (!strncmp(node->Name(), "Object#", 7)) {
    // Get id and type
    char nodetype[4096];
    strncpy(nodetype, node->Name(), 4096);
    char *id = strchr(nodetype, '#');
    if (!id) return 0;
    else *(id++) = '\0';

    // Get modelId
    char modelId[4096];
    if (node->NReferences() == 0) return 0; 
    R3SceneReference *reference = node->Reference(0);
    R3Scene *referenced_scene = reference->ReferencedScene();
    if (!referenced_scene->Name()) return 0;
    if (strncmp(referenced_scene->Root()->Name(), "Model#", 6)) return 0;
    strncpy(modelId, &(referenced_scene->Root()->Name()[6]), 4096);

    // Write json
    const R3Box& b = node->WorldBBox();
    R4Matrix matrix = node->Transformation().Matrix();
    RNBoolean is_mirrored = node->Transformation().IsMirrored();
    if (counter++ > 0) fprintf(fp, ",\n");
    fprintf(fp, "      {\n");
    fprintf(fp, "        \"id\":\"%s\",\n", id);
    fprintf(fp, "        \"type\":\"%s\",\n", nodetype);
    fprintf(fp, "        \"modelId\":\"%s\",\n", modelId);
    fprintf(fp, "        \"isMirrored\":%d,\n", is_mirrored);
    fprintf(fp, "        \"bbox\":{\"min\":[%g,%g,%g],\"max\":[%g,%g,%g]},\n",
      b[0][0], b[0][1], b[0][2], b[1][0], b[1][1], b[1][2]);
    if (reference->NMaterials() > 0) {
      fprintf(fp, "        \"materials\":[");
      for (int i = 0; i < reference->NMaterials(); i++) {
        R3Material *material = reference->Material(i);
        const char *material_name = (material) ? material->Name() : "none";
        const R3Brdf *brdf = (material) ? material->Brdf() : NULL;
        const R2Texture *texture = (material) ? material->Texture() : NULL;
        if (i > 0) fprintf(fp, ", ");
        fprintf(fp, "{ \"name\":\"%s\"", material_name);
        if (brdf) {
          const RNRgb& diffuse = brdf->Diffuse();
          int r = 255.0 * diffuse.R() + 0.5;
          int g = 255.0 * diffuse.G() + 0.5;
          int b = 255.0 * diffuse.B() + 0.5;
          r = (r < 255) ? ((r > 0) ? r : 0) : 255;
          g = (g < 255) ? ((g > 0) ? g : 0) : 255;
          b = (b < 255) ? ((b > 0) ? b : 0) : 255;
          fprintf(fp, ", \"diffuse\":\"#%02x%02x%02x\"", r, g, b);
        }
        if (texture && texture->Name()) {
          fprintf(fp, ", \"texture\":\"%s\"", texture->Name());
        }
        fprintf(fp, " }");
      }
      fprintf(fp, "],\n");
    }
    if (!matrix.IsIdentity()) {
      fprintf(fp, "        \"transform\":[");
      for (int j = 0; j < 4; j++) {
        for (int i = 0; i < 4; i++) {
          if ((i > 0) || (j > 0)) fprintf(fp, ", ");
          fprintf(fp, "%g", matrix[i][j]);
        }
      }
      fprintf(fp, "],\n");
    }
    fprintf(fp, "        \"valid\":1\n");
    fprintf(fp, "      }");
  }
  else if (!strncmp(node->Name(), "Box#", 4)) {
    fprintf(stderr, "Warning: skipping %s\n", node->Name());
  }
  
  // Return success
  return 1;
}



int R3Scene::
WriteSUNCGFile(const char *filename) const
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "w"))) {
    fprintf(stderr, "Unable to open file %s", filename);
    return 0;
  }

  // Write header
  fprintf(fp, "{\n");
  fprintf(fp, "  \"version\":\"suncg@1.0.0\",\n");
  if (root->Name()) fprintf(fp, "  \"id\":\"%s\",\n", root->Name());
  fprintf(fp, "  \"up\":[0,1,0],\n");
  fprintf(fp, "  \"front\":[0,0,1],\n");
  fprintf(fp, "  \"scaleToMeters\":1,\n");  
  fprintf(fp, "  \"levels\":[\n");  

  // Write nodes recursively
  for (int i = 0; i < root->NChildren(); i++) {
    R3SceneNode *child = root->Child(i);
    if (!WriteSUNCGNode(this, child, fp)) {
      fclose(fp);
      return 0;
    }
  }

  // Write trailer
  fprintf(fp, "\n  ]\n");
  fprintf(fp, "}\n");
  
  // Close file
  fclose(fp);
  
  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// SUNCG UTILITY FUNCTIONS
////////////////////////////////////////////////////////////////////////

static RNBoolean
MatchReferencedSceneName(R3SceneNode *node, const char *name)
{
  // Check node name
  if (node->Name() && !strcmp(node->Name(), name)) return TRUE;

  // Check referenced node names
  for (int i = 0; i < node->NReferences(); i++) {
    R3SceneReference *reference = node->Reference(i);
    R3Scene *referenced_scene = reference->ReferencedScene();
    if (referenced_scene->Name() && !strcmp(referenced_scene->Name(), name)) {
      return TRUE;
    }
  }

  // No match
  return FALSE;
}



static int
InsertCopiesOfLight(R3Scene *scene, R3Light *original, const char *reference_frame)
{
  // Iniitalize return status
  int status = 0;
  
  // Insert copies of light into scene
  if (!strcmp(reference_frame, "world")) {
    // Insert one copy of light 
    R3Light *light = original->Copy();
    scene->InsertLight(light);
    status++;
  }
  else {
    // Insert copy of light for each matching node with name matching reference frame
    for (int i = 0; i < scene->NNodes(); i++) {
      R3SceneNode *node = scene->Node(i);
      if (!MatchReferencedSceneName(node, reference_frame)) continue;
      R3Light *light = original->Copy();
      R3Affine transformation_to_world = node->CumulativeTransformation();
      light->Transform(node->CumulativeTransformation());
      scene->InsertLight(light);
      status++;
    }
  }

  // Return whether inserted any lights
  return status;
}



int R3Scene::
ReadSUNCGLightsFile(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open lights file %s\n", filename);
    return 0;
  }

  // Read file
  char buffer[4096];
  int line_number = 0;
  while (fgets(buffer, 4096, fp)) {
    line_number++;
    char cmd[1024], reference_frame[1024];
    if (sscanf(buffer, "%s", cmd) != (unsigned int) 1) continue;
    if (cmd[0] == '#') continue;

    // Check cmd
    if (!strcmp(cmd, "directional_light")) {
      // Parse directional light 
      double intensity, r, g, b, dx, dy, dz;
      if (sscanf(buffer, "%s%s%lf%lf%lf%lf%lf%lf%lf", cmd, reference_frame,
        &intensity, &r, &g, &b, &dx, &dy, &dz) != (unsigned int) 9) {
        fprintf(stderr, "Unable to parse directional light from line %d from %s\n", line_number, filename);
        return 0;
      }

      // Create directional light
      RNRgb color(r, g, b);
      R3Vector direction(dx, dy, dz);
      R3DirectionalLight light(direction, color, intensity);
      InsertCopiesOfLight(this, &light, reference_frame);
    }
    else if (!strcmp(cmd, "point_light")) {
      // Parse point light 
      double intensity, r, g, b, px, py, pz, ca, la, qa;
      if (sscanf(buffer, "%s%s%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", cmd, reference_frame,  
        &intensity, &r, &g, &b, &px, &py, &pz, &ca, &la, &qa) != (unsigned int) 12) {
        fprintf(stderr, "Unable to parse point light from line %d from %s\n", line_number, filename);
        return 0;
      }

      // Create point light
      RNRgb color(r, g, b);
      R3Point position(px, py, pz);
      R3PointLight light(position, color, intensity, TRUE, ca, la, qa);
      InsertCopiesOfLight(this, &light, reference_frame);
    }
    else if (!strcmp(cmd, "spot_light")) {
      // Parse spot light 
      double intensity, r, g, b, px, py, pz, dx, dy, dz, sd, sc, ca, la, qa;
      if (sscanf(buffer, "%s%s%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", cmd, reference_frame,  
        &intensity, &r, &g, &b, &px, &py, &pz, &dx, &dy, &dz, &sd, &sc, &ca, &la, &qa) != (unsigned int) 17) {
        fprintf(stderr, "Unable to parse spot light from line %d from %s\n", line_number, filename);
        return 0;
      }

      // Create spot light
      RNRgb color(r, g, b);
      R3Point position(px, py, pz);
      R3Vector direction(dx, dy, dz);
      R3SpotLight light(position, direction, color, sd, sc, intensity, TRUE, ca, la, qa);
      InsertCopiesOfLight(this, &light, reference_frame);
    }
    else if (!strcmp(cmd, "line_light")) {
      // Parse spot light 
      double intensity, r, g, b, px1, py1, pz1, px2, py2, pz2, dx, dy, dz, sd, sc, ca, la, qa;
      if (sscanf(buffer, "%s%s%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", cmd, reference_frame,  
        &intensity, &r, &g, &b, &px1, &py1, &pz1, &px2, &py2, &pz2, &dx, &dy, &dz, &sd, &sc, &ca, &la, &qa) != (unsigned int) 20) {
        fprintf(stderr, "Unable to parse line light from line %d from %s\n", line_number, filename);
        return 0;
      }

      // Create spot light
      RNRgb color(r, g, b);
      R3Point position1(px1, py1, pz1);
      R3Point position2(px2, py2, pz2);
      R3Point position = 0.5 * (position1 + position2);
      R3Vector direction(dx, dy, dz);
      R3SpotLight light(position, direction, color, sd, sc, intensity, TRUE, ca, la, qa);
      InsertCopiesOfLight(this, &light, reference_frame);
    }
    else {
      fprintf(stderr, "Unrecognized light type %s at line %d of %s\n", cmd, line_number, filename);
      return 0;
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}
  


int R3Scene::
ReadSUNCGModelFile(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open lights file %s\n", filename);
    return 0;
  }

  // Read keys from first line
  int line_number = 1;
  char key_buffer[4096];
  RNArray<char *> keys;
  if (fgets(key_buffer, 4096, fp)) {
    char *token = strtok(key_buffer, ",\r\n");
    while (token) {
      keys.Insert(strdup(token));
      token = strtok(NULL, ",\r\n");
    }
  }

  // Extract index of model_id
  int model_id_k = -1;
  for (int i = 0; i < keys.NEntries(); i++) {
    if (!strcmp(keys[i], "model_id")) {
      model_id_k = i;
      break;
    }
  }

  // Check if found model_id
  if (model_id_k < 0) {
    fprintf(stderr, "Did not find \"model_id\" in header on line %d of %s\n", line_number, filename);
    return 0;
  }

  // Read subsequent lines of file
  char value_buffer[4096];
  while (fgets(value_buffer, 4096, fp)) {
    line_number++;

    // Read values
    RNArray<char *> values;
    char *token = strtok(value_buffer, ",\r\n");
    while (token) {
      values.Insert(strdup(token));
      token = strtok(NULL, ",\r\n");
    }

    // Check number of values
    if (values.NEntries() == 0) continue;
    if (values.NEntries() != keys.NEntries()) {
      fprintf(stderr, "Invalid number of entries at line %d in %s\n", line_number, filename);
      return 0;
    }

    // Get model id
    const char *model_id = values[model_id_k];
    if (!model_id) continue;
    int model_id_length = strlen(model_id);
    if (model_id_length == 0) continue;

    // Assign key-value info to nodes matching model_id
    for (int i = 0; i < NNodes(); i++) {
      R3SceneNode *node = Node(i);

      // Determine if node name matches model_id
      if (node->Name()) {
        char node_name[1024];
        strncpy(node_name, node->Name(), 1024);
        char *model_name = strchr(node_name, '#');
        if (model_name) { *model_name = '\0'; model_name++; }
        else model_name = node_name;
        if (!strcmp(node_name, model_id) || !strcmp(model_name, model_id)) {
          for (int j = 0; j < keys.NEntries(); j++) node->InsertInfo(keys[j], values[j]);
        }
      }

      // Determine if referenced scene name matches model_id
      for (int j = 0; j < node->NReferences(); j++) {
        R3SceneReference *reference = node->Reference(j);
        R3Scene *referenced_scene = reference->ReferencedScene();
        R3SceneNode *root = referenced_scene->Root();
        if (!root->Name()) continue;
        char node_name[1024];
        strncpy(node_name, root->Name(), 1024);
        char *model_name = strchr(node_name, '#');
        if (model_name) { *model_name = '\0'; model_name++; }
        else model_name = node_name;
        if (!strcmp(node_name, model_id) || !strcmp(model_name, model_id)) {
          for (int j = 0; j < keys.NEntries(); j++) node->InsertInfo(keys[j], values[j]);
          break;
        }
      }
    }

    // Assign key-value info to root node of referenced scene matching model_id
    R3Scene *model = ReferencedScene(model_id);
    if (model) {
      for (int i = 0; i < keys.NEntries(); i++) {
        model->Root()->InsertInfo(keys[i], values[i]);
      }
    }
  }

  // Close file
  fclose(fp);

  // Delete inactive nodes
  RNArray<R3SceneNode *> tmp;
  for (int i = 0; i < NNodes(); i++)
    tmp.Insert(Node(i));
  for (int i = 0; i < tmp.NEntries(); i++) {
    R3SceneNode *node = tmp.Kth(i);
    const char *active = node->Info("active");
    if (!active) continue;
    if (strcmp(active, "0")) continue;
    delete node;
  }

  
#if 0
  // Print result
  for (int i = 0; i < NNodes(); i++) {
    R3SceneNode *node = Node(i);
    const char *model_index = node->Info("index");
    if (model_index) printf("%s %s\n", node->Name(), model_index);
  }
#endif
  
  // Return success
  return 1;
}

