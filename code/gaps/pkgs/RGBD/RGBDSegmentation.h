////////////////////////////////////////////////////////////////////////
// Include file for RGBDImage class
////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////
// Constant definitions
////////////////////////////////////////////////////////////////////////

// Primitive types

enum {
  RGBD_NULL_PRIMITIVE_TYPE,
  RGBD_POINT_PRIMITIVE_TYPE,
  RGBD_LINE_PRIMITIVE_TYPE,
  RGBD_PLANE_PRIMITIVE_TYPE,
  RGBD_NUM_PRIMITIVE_TYPES
};



////////////////////////////////////////////////////////////////////////
// Class definitions
////////////////////////////////////////////////////////////////////////

struct RGBDPoint {
public:
  RGBDPoint(void);
public:
  RNScalar depth;
  R3Point position;
  R3Vector normal;
  RNScalar radius;
  RNRgb color;
  unsigned int boundary;
  RNArray<RGBDPoint *> neighbors;
  struct RGBDSegment *segment;
  RNScalar segment_affinity;
  int segment_index;
  int grid_index;
  int mark;
};

struct RGBDPrimitive {
  RGBDPrimitive(int primitive_type = 0);
  RGBDPrimitive(const RGBDPrimitive& primitive);
  RGBDPrimitive(RGBDPoint *seed_point, const RNArray<RGBDPoint *> *points = NULL);
  RNLength Distance(const R3Point& position) const;
  void Update(const R3Point& point);
  void Update(const R3Line& line);
  void Update(const R3Plane& plane);
  void Update(RGBDPoint *seed_point = NULL, const RNArray<RGBDPoint *> *points = NULL);
  void Update(RGBDPrimitive primitive1, RGBDPrimitive primitive2, RNScalar weight1 = 1.0, RNScalar weight2 = 1.0);
public:
  int primitive_type;
  R3Box bbox;
  R3Point centroid;
  R3Line line;
  R3Plane plane;
};

struct RGBDSegment {
public:
  RGBDSegment(struct RGBDSegmentation *segmentation, RGBDPoint *seed_point = NULL, int primitive_type = 0);
  RGBDSegment(struct RGBDSegmentation *segmentation, RGBDPoint *seed_point, const RGBDPrimitive& primitive);
  RGBDSegment(struct RGBDSegmentation *segmentation, RGBDSegment *child1, RGBDSegment *child2);
  ~RGBDSegment(void);
  RNScalar Coverage(void);
  void EmptyPoints(void);
  void InsertPoint(RGBDPoint *point, RNScalar affinity = 1.0);
  void RemovePoint(RGBDPoint *point);
  void InsertChild(RGBDSegment *child);
  void RemoveChild(RGBDSegment *child);
  int UpdatePoints(const R3Kdtree<RGBDPoint *> *kdtree);
  int UpdatePrimitive(void);
  RNScalar Affinity(RGBDPoint *point) const;
  RNScalar Affinity(RGBDSegment *segment) const;
public:
  struct RGBDSegmentation *segmentation;
  int segmentation_index;
  RGBDPoint *seed_point;
  RNArray<RGBDPoint *> points;
  RGBDSegment *parent;
  RNArray<RGBDSegment *> children;
  RNArray<struct RGBDSegmentPair *> pairs;
  RGBDPrimitive primitive;
  RNScalar possible_affinity; 
  RNScalar total_affinity;
};

struct RGBDSegmentPair {
public:
  RGBDSegmentPair(RGBDSegment *segment1 = NULL, RGBDSegment *segment2 = NULL, RNScalar affinity = 0);
  ~RGBDSegmentPair(void);
public:
  RGBDSegment *segments[2];
  int segment_index[2];
  RNScalar affinity; 
  RGBDSegmentPair **heapentry;
};

struct RGBDSegmentation {
public:
  RGBDSegmentation(void);
  RGBDSegmentation(const R2Grid& px_image, const R2Grid& py_image, const R2Grid& pz_image,
    const R2Grid& nx_image, const R2Grid& ny_image, const R2Grid& nz_image,
    const R2Grid& depth_image, const R2Grid& radius_image, 
    const R2Grid& boundary_image, const R2Image& color_image,
    const R3Point& viewpoint, const R3Vector& towards, const R3Vector& up,
    int primitive_type = RGBD_PLANE_PRIMITIVE_TYPE);
  ~RGBDSegmentation(void);
  RNScalar Affinity(void) const;
  int NUnsegmentedPoints(void) const;
public:
  int CreatePoints(const R2Grid& px_image, const R2Grid& py_image, const R2Grid& pz_image, 
    const R2Grid& nx_image, const R2Grid& ny_image, const R2Grid& nz_image,
    const R2Grid& depth_image, const R2Grid& radius_image,
    const R2Grid& boundary_image, const R2Image& color_image);
  int CreateSegments(int primitive_type);
  int CreateSingletonSegments(int primitive_type);
  int CreateRansacSegments(int primitive_type);
  int RefineSegments(void);  
  int DeleteSegments(void);  
  int MergeSegments(void);  
  int SplitSegments(void);
public:
  int ReadSegmentImage(const char *filename);
  int WriteSegmentImage(int xres, int yres, const char *filename) const;
public:
  RNArray<RGBDPoint *> points;
  R3Kdtree<RGBDPoint *> *kdtree;
  RNArray<RGBDSegment *> segments;
  RGBDPoint *point_buffer;

  // Parameters
  void SetDefaultParameters(void);
  int min_segment_points;
  int min_segments;
  int max_segments;
  RNScalar min_segment_coverage;
  RNLength max_segment_diameter;
  RNLength max_segment_primitive_distance;
  RNAngle max_segment_normal_angle;
  RNScalar max_neighbor_distance_factor;
  RNLength max_pair_centroid_distance;
  RNLength max_pair_primitive_distance;
  RNAngle max_pair_normal_angle;
  RNScalar min_pair_affinity;
  RNBoolean initialize_hierarchically; 
  int max_refinement_iterations;
  int max_ransac_iterations;
  int print_progress;
};




