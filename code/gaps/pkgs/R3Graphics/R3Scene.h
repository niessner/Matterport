/* Include file for the R3 scene class */



/* Initialization functions */

int R3InitScene();
void R3StopScene();



/* Class definition */

class R3Scene {
public:
  // Constructor functions
  R3Scene(void);
  ~R3Scene(void);

  // Property functions
  const R3Box& BBox(void) const;
  const R3Point Centroid(void) const;
  const RNInterval NFacets(void) const;
  const RNLength Length(void) const;
  const RNArea Area(void) const;
  const RNVolume Volume(void) const;
  const R3Point ClosestPoint(const R3Point& point) const;
  const char *Name(void) const;
  const char *Info(const char *key) const;
  void *Data(void) const;

  // Access functions
  int NNodes(void) const;
  R3SceneNode *Node(int k) const;
  R3SceneNode *Node(const char *name) const;
  R3SceneNode *Root(void) const;
  int NLights(void) const;
  R3Light *Light(int k) const;
  R3Light *Light(const char *name) const;
  int NMaterials(void) const;
  R3Material *Material(int k) const;
  R3Material *Material(const char *name) const;
  int NBrdfs(void) const;
  R3Brdf *Brdf(int k) const;
  R3Brdf *Brdf(const char *name) const;
  int NTextures(void) const;
  R2Texture *Texture(int k) const;
  R2Texture *Texture(const char *name) const;
  int NReferencedScenes(void) const;
  R3Scene *ReferencedScene(int k) const;
  R3Scene *ReferencedScene(const char *name) const;
  const R3Camera& Camera(void) const;
  const R2Viewport& Viewport(void) const;
  const R3Viewer& Viewer(void) const;
  const RNRgb& Ambient(void) const;
  const RNRgb& Background(void) const;
  const char *Filename(void) const;

  // Manipulation functions
  void InsertNode(R3SceneNode *node);
  void RemoveNode(R3SceneNode *node);
  void InsertLight(R3Light *light);
  void RemoveLight(R3Light *light);
  void InsertMaterial(R3Material *material);
  void RemoveMaterial(R3Material *material);
  void InsertBrdf(R3Brdf *brdf);
  void RemoveBrdf(R3Brdf *brdf);
  void InsertTexture(R2Texture *texture);
  void RemoveTexture(R2Texture *texture);
  void InsertReferencedScene(R3Scene *referenced_scene);
  void RemoveReferencedScene(R3Scene *referenced_scene);
  void InsertInfo(const char *key, const char *info);
  void ReplaceInfo(const char *key, const char *info);
  void RemoveInfo(const char *key);
  void SetCamera(const R3Camera& viewer);
  void SetViewport(const R2Viewport& viewport);
  void SetViewer(const R3Viewer& viewer);
  void SetAmbient(const RNRgb& ambient);
  void SetBackground(const RNRgb& background);
  void SetFilename(const char *filename);
  void SetName(const char *name);
  void SetData(void *data);
  void RemoveReferences(void);
  void RemoveHierarchy(void);
  void RemoveTransformations(void);
  void SubdivideTriangles(RNLength max_edge_length);
  void CreateDirectionalLights(void);

  // Query functions
  RNLength Distance(const R3Point& point) const;
  RNBoolean FindClosest(const R3Point& point,
    R3SceneNode **hit_node = NULL, R3Material **hit_material = NULL, R3Shape **hit_shape = NULL,
    R3Point *hit_point = NULL, R3Vector *hit_normal = NULL, RNScalar *hit_d = NULL,
    RNScalar min_d = 0.0, RNScalar max_d = RN_INFINITY) const;
  RNBoolean Intersects(const R3Ray& ray,
    R3SceneNode **hit_node = NULL, R3Material **hit_material = NULL, R3Shape **hit_shape = NULL,
    R3Point *hit_point = NULL, R3Vector *hit_normal = NULL, RNScalar *hit_t = NULL,
    RNScalar min_t = 0.0, RNScalar max_t = RN_INFINITY) const;

  // I/O functions
  int ReadFile(const char *filename, R3SceneNode *parent_node = NULL);
  int ReadObjFile(const char *filename, R3SceneNode *parent_node = NULL);
  int ReadPlyFile(const char *filename, R3SceneNode *parent_node = NULL);
  int ReadMeshFile(const char *filename, R3SceneNode *parent_node = NULL);
  int ReadSUNCGFile(const char *filename, R3SceneNode *parent_node = NULL);
  int ReadPlanner5DFile(const char *filename, R3SceneNode *parent_node = NULL);
  int ReadPrincetonFile(const char *filename, R3SceneNode *parent_node = NULL);
  int ReadParseFile(const char *filename, R3SceneNode *parent_node = NULL);
  int ReadSupportHierarchyFile(const char *filename, R3SceneNode *parent_node = NULL);
  int ReadGrammarHierarchyFile(const char *filename, R3SceneNode *parent_node = NULL);
  int ReadRectangleFile(const char *filename, R3SceneNode *parent_node = NULL);
  int WriteFile(const char *filename) const;
  int WriteObjFile(const char *filename) const;
  int WriteSUNCGFile(const char *filename) const;
  int WritePrincetonFile(const char *filename) const;

  // Draw functions
  void Draw(const R3DrawFlags draw_flags = R3_DEFAULT_DRAW_FLAGS) const;
  
  // Lighting functions
  int LoadLights(int min_index = 0, int max_index = 7) const;
  int LoadLights(const R3Box& world_bbox, int min_index = 0, int max_index = 7) const;

  // SUNCG utility functions
  int ReadSUNCGLightsFile(const char *filename);
  int ReadSUNCGModelFile(const char *filename);

private:
  R3SceneNode *root;
  RNArray<R3SceneNode *> nodes;
  RNArray<R3Light *> lights;
  RNArray<R3Material *> materials;
  RNArray<R3Brdf *> brdfs;
  RNArray<R2Texture *> textures;
  RNArray<R3Scene *> referenced_scenes;
  RNSymbolTable<const char *> info;
  R3Viewer viewer;
  RNRgb ambient;
  RNRgb background;
  char *filename;
  char *name;
  void *data;
};



/* Inline functions */

inline const R3Box& R3Scene::
BBox(void) const
{
  // Return bounding box of root node
  return root->BBox();
}



inline const R3Point R3Scene::
Centroid(void) const
{
  // Return centroid of root node
  return root->Centroid();
}



inline const RNInterval R3Scene::
NFacets(void) const
{
  // Return the range of how many facets can be drawn for scene
  return root->NFacets();
}



inline const RNLength R3Scene::
Length(void) const
{
  // Return total perimeter of all facets 
  return root->Length();
}



inline const RNArea R3Scene::
Area(void) const
{
  // Return surface area of scene
  return root->Area();
}



inline const RNVolume R3Scene::
Volume(void) const
{
  // Return volume of scene
  return root->Volume();
}



inline const R3Point R3Scene::
ClosestPoint(const R3Point& point) const
{
  // Return closest point on surface of scene
  return root->ClosestPoint(point);
}



inline const char *R3Scene::
Info(const char *key) const
{
  // Return info associated with key
  const char *value;
  if (!info.Find(key, &value)) return NULL;
  return value;
}



inline const char *R3Scene::
Name(void) const
{
  // Return name
  return name;
}



inline void *R3Scene::
Data(void) const
{
  // Return user-defined data
  return data;
}



inline R3SceneNode *R3Scene::
Root(void) const
{
  // Return root node
  return root;
}



inline int R3Scene::
NNodes(void) const
{
  // Return number of nodes
  return nodes.NEntries();
}



inline R3SceneNode *R3Scene::
Node(int k) const
{
  // Return kth node
  return nodes.Kth(k);
}



inline int R3Scene::
NLights(void) const
{
  // Return number of lights
  return lights.NEntries();
}



inline R3Light *R3Scene::
Light(int k) const
{
  // Return kth light
  return lights.Kth(k);
}



inline int R3Scene::
NMaterials(void) const
{
  // Return number of materials
  return materials.NEntries();
}



inline R3Material *R3Scene::
Material(int k) const
{
  // Return kth material
  return materials.Kth(k);
}



inline int R3Scene::
NBrdfs(void) const
{
  // Return number of brdfs
  return brdfs.NEntries();
}



inline R3Brdf *R3Scene::
Brdf(int k) const
{
  // Return kth brdf
  return brdfs.Kth(k);
}



inline int R3Scene::
NTextures(void) const
{
  // Return number of textures
  return textures.NEntries();
}



inline R2Texture *R3Scene::
Texture(int k) const
{
  // Return kth texture
  return textures.Kth(k);
}



inline int R3Scene::
NReferencedScenes(void) const
{
  // Return number of referenced scenes
  return referenced_scenes.NEntries();
}



inline R3Scene *R3Scene::
ReferencedScene(int k) const
{
  // Return kth referenced scene
  return referenced_scenes.Kth(k);
}



inline const R3Camera& R3Scene::
Camera(void) const
{
  // Return camera
  return viewer.Camera();
}


inline const R2Viewport& R3Scene::
Viewport(void) const
{
  // Return viewport
  return viewer.Viewport();
}


inline const R3Viewer& R3Scene::
Viewer(void) const
{
  // Return viewer
  return viewer;
}


inline const RNRgb& R3Scene::
Ambient(void) const
{
  // Return ambient light color
  return ambient;
}



inline const RNRgb& R3Scene::
Background(void) const
{
  // Return background color
  return background;
}



inline const char *R3Scene::
Filename(void) const
{
  // Return filename
  return filename;
}



inline void R3Scene::
SetAmbient(const RNRgb& ambient) 
{
  // Set ambient light color
  this->ambient = ambient;
}



inline void R3Scene::
SetBackground(const RNRgb& background) 
{
  // Set background color
  this->background = background;
}



inline void R3Scene::
SetFilename(const char *filename)
{
  // Set filename
  if (this->filename) free(this->filename);
  if (filename) this->filename = strdup(filename);
  else this->filename = NULL;
}


  
