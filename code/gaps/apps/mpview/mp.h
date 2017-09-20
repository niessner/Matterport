////////////////////////////////////////////////////////////////////////
// Include file for mp structures
////////////////////////////////////////////////////////////////////////
#ifndef __MP__
#define __MP__



////////////////////////////////////////////////////////////////////////

#define MP_DEFAULT_DRAW_FLAGS 0xFFFFFFFF



////////////////////////////////////////////////////////////////////////

struct MPImage {
  // Constructor/destructor
  MPImage(void);
  ~MPImage(void);

  // Drawing stuff
  void Draw(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawCamera(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawBBox(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawImage(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawPoints(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawQuads(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;

public:
  struct MPHouse *house;
  int house_index;
  struct MPPanorama *panorama;
  int panorama_index;
  char *name;
  int camera_index;
  int yaw_index;
  RGBDImage rgbd;
  R4Matrix extrinsics;
  R3Matrix intrinsics;
  int width, height;
  R3Point position;
};



////////////////////////////////////////////////////////////////////////

struct MPPanorama {
  // Constructor/destructor
  MPPanorama(void);
  ~MPPanorama(void);

  // Manipulation stuff
  void InsertImage(MPImage *image);
  void RemoveImage(MPImage *image);

  // Drawing stuff
  void Draw(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawPosition(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawName(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawImages(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;

public:
  struct MPHouse *house;
  int house_index;
  struct MPRegion *region;
  int region_index;
  char *name;
  RNArray<MPImage *> images;
  R3Point position;
};



////////////////////////////////////////////////////////////////////////

struct MPSegment {
  // Constructor/destructor
  MPSegment(void);
  ~MPSegment(void);

  // Drawing stuff
  void Draw(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawMesh(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawBBox(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;

public:
  struct MPHouse *house;
  int house_index;
  struct MPObject *object;
  int object_index;
  R3Mesh *mesh;
  RNArray<R3MeshFace *> faces;
  RNArea area;
  R3Point position;
  R3Box bbox;
  int id;
};



////////////////////////////////////////////////////////////////////////

struct MPObject {
  // Constructor/destructor
  MPObject(void);
  ~MPObject(void);

  // Manipulation stuff
  void InsertSegment(MPSegment *segment);
  void RemoveSegment(MPSegment *segment);

  // Drawing stuff
  void Draw(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawBBox(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawLabel(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawSegments(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;

public:
  struct MPHouse *house;
  int house_index;
  struct MPRegion *region;
  int region_index;
  struct MPCategory *category;
  int category_index;
  RNArray<MPSegment *> segments;
  R3Point position;
  R3OrientedBox obb;
};




////////////////////////////////////////////////////////////////////////

struct MPCategory {
  // Constructor/destructor
  MPCategory(void);
  ~MPCategory(void);

  // Manipulation stuff
  void InsertObject(MPObject *object);
  void RemoveObject(MPObject *segment);

public:
  struct MPHouse *house;
  int house_index;
  RNArray<MPObject *> objects;
  int label_id;
  char *label_name;
  int mpcat40_id;
  char *mpcat40_name;  
};



////////////////////////////////////////////////////////////////////////

struct MPVertex {
  // Constructor/destructor
  MPVertex(void);
  ~MPVertex(void);

  // Drawing stuff
  void Draw(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawPosition(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;

public:
  struct MPHouse *house;
  int house_index;
  struct MPSurface *surface;
  int surface_index;
  R3Point position;
  R3Vector normal;
  char *label;
};



////////////////////////////////////////////////////////////////////////

struct MPSurface {
  // Constructor/destructor
  MPSurface(void);
  ~MPSurface(void);

  // Manipulation stuff
  void InsertVertex(MPVertex *vertex, RNBoolean search_for_best_index = FALSE);
  void RemoveVertex(MPVertex *vertex);
  void FlipOrientation(void);

  // Drawing stuff
  void Draw(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawPolygon(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawVertices(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;

  // Processing stuff
  void RecomputeBBox(void);

public:
  struct MPHouse *house;
  int house_index;
  struct MPRegion *region;
  int region_index;
  RNArray<MPVertex *> vertices;
  R3Point position;
  R3Vector normal;
  R3Box bbox;
  char *label;
};



////////////////////////////////////////////////////////////////////////

struct MPRegion {
  // Constructor/destructor
  MPRegion(void);
  ~MPRegion(void);

  // Property stuff
  R2Polygon FloorPolygon(void) const;

  // Manipulation stuff
  void InsertPanorama(MPPanorama *panorama);
  void RemovePanorama(MPPanorama *panorama);
  void InsertSurface(MPSurface *surface);
  void RemoveSurface(MPSurface *surface);
  void InsertObject(struct MPObject *object);
  void RemoveObject(struct MPObject *object);
  void InsertPortal(struct MPPortal *portal, int index);
  void RemovePortal(struct MPPortal *portal);

  // Drawing stuff
  void Draw(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawPosition(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawBBox(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawLabel(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawPortals(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawObjects(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawSurfaces(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawPanoramas(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;

  // Processing stuff
  void RecomputeBBox(RNBoolean preserve_zmax = FALSE);

public:
  struct MPHouse *house;
  int house_index;
  struct MPLevel *level;
  int level_index;
  RNArray<MPPanorama *> panoramas;
  RNArray<MPSurface *> surfaces;
  RNArray<struct MPObject *> objects;
  RNArray<struct MPPortal *> portals;
  R3Point position;
  RNLength height;
  R3Box bbox;
  char *label;
};



////////////////////////////////////////////////////////////////////////

struct MPPortal {
  // Constructor/destructor
  MPPortal(void);
  ~MPPortal(void);

  // Drawing stuff
  void Draw(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawSpan(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawLabel(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  

public:
  struct MPHouse *house;
  int house_index;
  MPRegion *regions[2];
  int region_indices[2];
  R3Span span;
  char *label;
};



////////////////////////////////////////////////////////////////////////

struct MPLevel {
  // Constructor/destructor
  MPLevel(void);
  ~MPLevel(void);

  // Manipulation stuff
  void InsertRegion(MPRegion *region);
  void RemoveRegion(MPRegion *region);

  // Drawing stuff
  void Draw(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawPosition(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawBBox(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawRegions(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;

public:
  struct MPHouse *house;
  int house_index;
  RNArray<MPRegion *> regions;
  R3Point position;
  R3Box bbox;
  char *label;
};



////////////////////////////////////////////////////////////////////////

struct MPHouse {
  // Constructor/destructor
  MPHouse(const char *name = NULL, const char *label = NULL);
  ~MPHouse(void);

  // Manipulation stuff
  void InsertImage(MPImage *image);
  void RemoveImage(MPImage *image);
  void InsertPanorama(MPPanorama *panorama);
  void RemovePanorama(MPPanorama *panorama);
  void InsertCategory(MPCategory *category);
  void RemoveCategory(MPCategory *category);
  void InsertSegment(MPSegment *segment);
  void RemoveSegment(MPSegment *segment);
  void InsertObject(MPObject *object);
  void RemoveObject(MPObject *object);
  void InsertVertex(MPVertex *vertex);
  void RemoveVertex(MPVertex *vertex);
  void InsertSurface(MPSurface *surface);
  void RemoveSurface(MPSurface *surface);
  void InsertRegion(MPRegion *region);
  void RemoveRegion(MPRegion *region);
  void InsertPortal(MPPortal *portal);
  void RemovePortal(MPPortal *portal);
  void InsertLevel(MPLevel *level);
  void RemoveLevel(MPLevel *level);

  // Input/output stuff
  int ReadFile(const char *filename);
  int ReadAsciiFile(const char *filename);
  int WriteFile(const char *filename) const;
  int WriteAsciiFile(const char *filename) const;

  // Other input stuff
  int ReadMeshFile(const char *filename);
  int ReadSceneFile(const char *filename);
  int ReadSegmentFile(const char *filename);
  int ReadObjectFile(const char *filename);
  int ReadCategoryFile(const char *filename);
  int ReadConfigurationFile(const char *filename);

  // Drawing stuff
  void Draw(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawLevels(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawPortals(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawRegions(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawSurfaces(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawVertices(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawObjects(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawSegments(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawPanoramas(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawImages(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawScene(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;
  void DrawMesh(RNFlags draw_flags = MP_DEFAULT_DRAW_FLAGS) const;

public:
  // Utility stuff
  MPRegion *FindRegion(const R3Point& position, RNLength max_distance = 1) const;
  MPLevel *FindLevel(const R3Point& position, RNLength max_dz = 1.0) const;
  MPCategory *FindCategory(const char *label_name) const;
  MPCategory *FindCategory(int label_id) const;
  MPSegment *FindSegment(int segment_id) const;
  
  // Geometric query stuff
  MPImage *FindClosestImage(const R3Point& query_position, const R3Vector& normal = R3zero_vector,
    RNLength max_distance = FLT_MAX, RNBoolean check_normal = FALSE, RNBoolean check_visiblity = FALSE) const;
  MPVertex *FindClosestVertex(const R3Point& query_position, const R3Vector& normal = R3zero_vector,
    RNLength max_distance = FLT_MAX, RNBoolean check_normal = FALSE, RNBoolean check_visiblity = FALSE) const;
  MPRegion *FindClosestRegion(const R3Point& query_position, const R3Vector& normal = R3zero_vector, 
    RNLength max_distance = FLT_MAX, RNBoolean check_normal = FALSE, RNBoolean check_visiblity = FALSE) const;

public:
  RNArray<MPImage *> images;
  RNArray<MPPanorama *> panoramas;
  RNArray<MPCategory *> categories;
  RNArray<MPSegment *> segments;
  RNArray<MPObject *> objects;
  RNArray<MPVertex *> vertices;
  RNArray<MPSurface *> surfaces;
  RNArray<MPRegion *> regions;
  RNArray<MPPortal *> portals;
  RNArray<MPLevel *> levels;
  RGBDConfiguration rgbd;
  R3Scene *scene;
  R3Mesh *mesh;
  R3Box bbox;
  char *label;
  char *name;
};



////////////////////////////////////////////////////////////////////////

// Constants for defining drawing modes

#define MP_DRAW_VERTICES         0x00000001
#define MP_DRAW_EDGES            0x00000002
#define MP_DRAW_FACES            0x00000004
#define MP_DRAW_BBOXES           0x00000008
#define MP_DRAW_DEPICTIONS       0x00000010
#define MP_DRAW_LABELS           0x00000020
#define MP_DRAW_IMAGES           0x00000040
#define MP_DRAW_FLAGS            0x000000FF

#define MP_SHOW_IMAGES           0x00000100
#define MP_SHOW_PANORAMAS        0x00000200
#define MP_SHOW_VERTICES         0x00000400
#define MP_SHOW_SURFACES         0x00000800
#define MP_SHOW_SEGMENTS         0x00001000
#define MP_SHOW_OBJECTS          0x00002000
#define MP_SHOW_REGIONS          0x00004000
#define MP_SHOW_PORTALS          0x00008000
#define MP_SHOW_LEVELS           0x00010000
#define MP_SHOW_MESH             0x00020000
#define MP_SHOW_SCENE            0x00040000
#define MP_SHOW_FLAGS            0x000FFF00

#define MP_COLOR_BY_IMAGE        0x00100000
#define MP_COLOR_BY_PANORAMA     0x00200000
#define MP_COLOR_BY_SURFACE      0x00400000
#define MP_COLOR_BY_SEGMENT      0x00800000
#define MP_COLOR_BY_OBJECT       0x01000000
#define MP_COLOR_BY_REGION       0x02000000
#define MP_COLOR_BY_PORTAL       0x04000000
#define MP_COLOR_BY_LEVEL        0x08000000
#define MP_COLOR_BY_LABEL        0x10000000
#define MP_COLOR_BY_INDEX        0x20000000
#define MP_COLOR_BY_RGB          0x40000000
#define MP_COLOR_FOR_PICK        0x80000000
#define MP_COLOR_FLAGS           0xFFF00000



// Constants for picking tags

#define MP_IMAGE_TAG             0xF0
#define MP_PANORAMA_TAG          0xF1
#define MP_OBJECT_TAG            0xF2
#define MP_SEGMENT_TAG           0xF3
#define MP_VERTEX_TAG            0xF4
#define MP_SURFACE_TAG           0xF5
#define MP_REGION_TAG            0xF6
#define MP_PORTAL_TAG            0xF7
#define MP_LEVEL_TAG             0xF8
#define MP_MESH_TAG              0xF9
#define MP_SCENE_TAG             0xFA



#endif
