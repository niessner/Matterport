////////////////////////////////////////////////////////////////////////
// Include file for RGBDSurface class
////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////
// Class definition
////////////////////////////////////////////////////////////////////////

class RGBDSurface {
public:
  // Constructors/destructors
  RGBDSurface(void);
  RGBDSurface(const char *texture_filename, R3Rectangle *rectangle, RNLength texel_spacing);
  RGBDSurface(const char *texture_filename, R3Mesh *mesh, RNLength texel_spacing);
  virtual ~RGBDSurface(void);

  // Geometric property functions
  R3Point WorldCentroid(void) const;
  R3Box WorldBBox(void) const;
  RNArea WorldArea(void) const;
  RNScalar WorldTexelSpacing(void) const;

  // Configuration access functions
  RGBDConfiguration *Configuration(void) const;
  int ConfigurationIndex(void) const;
  
  // Channel access functions
  int NChannels(void) const;
  R2Grid *Channel(int channel_index) const;
  R2Grid *RedChannel(void) const;
  R2Grid *GreenChannel(void) const;
  R2Grid *BlueChannel(void) const;

  // Texture pixel property access functions
  int NTexels(int dim) const;
  RNRgb TexelColor(int ix, int iy) const;
  RNScalar TexelChannelValue(int ix, int iy, int channel_index) const;
  R3Point TexelWorldPosition(int ix, int iy) const;
  R3Vector TexelWorldNormal(int ix, int iy) const;

  // Texture point property access functions
  RNRgb TexelColor(const R2Point& texture_position) const;
  R3Point TexelWorldPosition(const R2Point& texture_position) const;
  R3Vector TexelWorldNormal(const R2Point& texture_position) const;
  RNScalar TexelChannelValue(const R2Point& texture_position, int channel_index) const;

  // Surface point property access functions
  RNRgb SurfelColor(const R2Point& surface_position) const;
  R3Point SurfelWorldPosition(const R2Point& surface_position) const;
  R3Vector SurfelWorldNormal(const R2Point& surface_position) const;
  RNScalar SurfelChannelValue(const R2Point& surface_position, int channel_index) const;

  // Manipulation functions
  virtual void SetTexelColor(int ix, int iy, const RNRgb& color);
  virtual void SetTexelChannelValue(int ix, int iy, int channel_index, RNScalar value);
  virtual void SetChannel(int channel_index, const R2Grid& image);
  virtual void SetColorChannels(const R2Image& image);
  virtual void Transform(const R3Transformation& transformation);

  // Transformation functions
  virtual int TransformTextureToSurface(const R2Point& texture_position, R2Point& surface_position) const;
  virtual int TransformSurfaceToTexture(const R2Point& surface_position, R2Point& texture_position) const;
  virtual int TransformSurfaceToWorld(const R2Point& surface_position, R3Point& world_position) const;
  virtual int TransformWorldToSurface(const R3Point& world_position, R2Point& surface_position) const;

  // Draw functions
  virtual void Draw(int color_scheme = RGBD_PHOTO_COLOR_SCHEME) const;
  virtual void DrawFaces(int color_scheme = RGBD_INDEX_COLOR_SCHEME) const;
  virtual void DrawEdges(int color_scheme = RGBD_INDEX_COLOR_SCHEME) const;
  virtual void DrawTexture(int color_scheme = RGBD_PHOTO_COLOR_SCHEME) const;
  virtual void LoadColor(int color_scheme = RGBD_PHOTO_COLOR_SCHEME) const;
  
  // Read/write functions
  virtual int ReadChannels(void);
  virtual int ReadColorChannels(void);
  virtual int WriteChannels(void);
  virtual int WriteColorChannels(void);
  virtual int ReleaseChannels(void);
  virtual int ReleaseColorChannels(void);

public:
  // Channel creation functions
  virtual int CreateColorChannels(const R2Image& image);

  // Texture filename functions
  const char *TextureFilename(void) const;
  virtual void SetTextureFilename(const char *filename);
  
  // Mesh filename functions
  const char *MeshFilename(void) const;
  virtual void SetMeshFilename(const char *filename);

  // Texel spacing functions
  virtual void SetWorldTexelSpacing(RNLength spacing);
  
  // Update functions
  virtual void InvalidateOpenGL(void);
  virtual void UpdateOpenGL(void);

  // Mesh search functions
  virtual R3MeshFace *SearchMeshFaceIndex(const R2Point& texture_position, double *barycentric_coordinates = NULL) const;
  virtual void InvalidateMeshFaceIndex(void);
  virtual void UpdateMeshFaceIndex(void);
  
private:
  // Internal variables
  friend class RGBDConfiguration;
  RGBDConfiguration *configuration;
  int configuration_index;
  RNArray<R2Grid *> channels;
  int width, height;
  int opengl_texture_id;
  char *texture_filename;
  char *mesh_filename;
  R2Box surface_bbox;
  int color_resident_count;

public: // temporary
  // Rectangular surface variables
  R3Rectangle *rectangle;

  // Mesh surface variables
  R3Mesh *mesh;
  R2Grid *mesh_face_index;
};



////////////////////////////////////////////////////////////////////////
// Inline functions
////////////////////////////////////////////////////////////////////////

inline RGBDConfiguration *RGBDSurface::
Configuration(void) const
{
  // Return configuration of which this surface is part
  return configuration;
}



inline int RGBDSurface::
ConfigurationIndex(void) const
{
  // Return index of this surface in configuration (i.e., where configuration->Surface(k) == this)
  return configuration_index;
}



inline int RGBDSurface::
NChannels(void) const
{
  // Return number of channels
  return channels.NEntries();
}



inline R2Grid *RGBDSurface::
Channel(int channel_index) const
{
  // Check channel
  if (!channels[channel_index]) {
    fprintf(stderr, "RGBD Channel is not resident in memory -- cannot get it\n");
    return NULL;
  }
  
  // Return channel
  return channels[channel_index];
}



inline R2Grid *RGBDSurface::
RedChannel(void) const
{
  // Return channel
  return Channel(RGBD_RED_CHANNEL);
}



inline R2Grid *RGBDSurface::
GreenChannel(void) const
{
  // Return channel
  return Channel(RGBD_GREEN_CHANNEL);
}



inline R2Grid *RGBDSurface::
BlueChannel(void) const
{
  // Return channel
  return Channel(RGBD_BLUE_CHANNEL);
}



inline int RGBDSurface::
NTexels(int dim) const
{
  // Return number of pixels in dimension
  assert((dim >= 0) && (dim < 2));
  return (dim == 0) ? width : height;
}



inline RNScalar RGBDSurface::
TexelChannelValue(int ix, int iy, int channel_index) const
{
  // Check channel
  if (!channels[channel_index]) {
    fprintf(stderr, "RGBD Channel is not resident in memory -- cannot get value\n");
    return R2_GRID_UNKNOWN_VALUE;
  }
  
  // Return value in channel at pixel position
  return channels[channel_index]->GridValue(ix, iy);
}



inline RNScalar RGBDSurface::
TexelChannelValue(const R2Point& image_position, int channel_index) const
{
  // Check channel
  if (!channels[channel_index]) {
    fprintf(stderr, "RGBD Channel is not resident in memory -- cannot get value\n");
    return R2_GRID_UNKNOWN_VALUE;
  }
  
  // Return value in channel at pixel position
  return channels[channel_index]->GridValue(image_position);
}



inline const char *RGBDSurface::
TextureFilename(void) const
{
  // Return texture filename
  return texture_filename;
}



inline const char *RGBDSurface::
MeshFilename(void) const
{
  // Return mesh filename
  return mesh_filename;
}



inline R3Point RGBDSurface::
WorldCentroid(void) const
{
  // Return centroid of surface
  if (rectangle) return rectangle->Centroid();
  else if (mesh) return mesh->Centroid();
  return R3zero_point;
}



inline R3Box RGBDSurface::
WorldBBox(void) const
{
  // Return bounding box
  if (rectangle) return rectangle->BBox();
  else if (mesh) return mesh->BBox();
  return R3null_box;
}



inline RNArea RGBDSurface::
WorldArea(void) const
{
  // Return area of surface
  if (rectangle) return rectangle->Area();
  else if (mesh) return mesh->Area();
  return 0.0;
}



inline RNLength RGBDSurface::
WorldTexelSpacing(void) const
{
  // Return world texel spacing
  if (width == 0) return 0;
  return surface_bbox.XLength() / width;
}


