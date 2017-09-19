////////////////////////////////////////////////////////////////////////
// Include file for RGBDConfiguration class
////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////
// Class definition
////////////////////////////////////////////////////////////////////////

class RGBDConfiguration {
public:
  // Constructors/destructors
  RGBDConfiguration(void);
  virtual ~RGBDConfiguration(void);

  // Property functions
  R3Point WorldCentroid(void) const;
  R3Box WorldBBox(void) const;

  // Image and surface access functions
  int NImages(void) const;
  RGBDImage *Image(int k) const;
  int NSurfaces(void) const;
  RGBDSurface *Surface(int k) const;

  // Directory name access functions
  const char *ColorDirectory(void) const;
  const char *DepthDirectory(void) const;
  const char *TextureDirectory(void) const;
  const char *DatasetFormat(void) const;
  
  // Image and surface manipulation functions
  virtual void InsertImage(RGBDImage *image);
  virtual void RemoveImage(RGBDImage *image);
  virtual void InsertSurface(RGBDSurface *surface);
  virtual void RemoveSurface(RGBDSurface *surface);
  virtual void Transform(const R3Transformation& transformation);

  // Direction name manipulation functions
  void SetColorDirectory(const char *name);
  void SetDepthDirectory(const char *name);
  void SetTextureDirectory(const char *name);
  void SetDatasetFormat(const char *name);

  // Input/output
  virtual int ReadFile(const char *filename, int read_every_kth_image = 1);
  virtual int WriteFile(const char *filename, int write_every_kth_image = 1) const;
  
public:
  // Read/release functions
  virtual int ReadChannels(void);
  virtual int ReleaseChannels(void);
  virtual int ReadDepthChannels(void);
  virtual int ReleaseDepthChannels(void);
  virtual int ReadColorChannels(void);
  virtual int ReleaseColorChannels(void);

  // Update functions
  virtual void InvalidateWorldBBox(void);
  virtual void UpdateWorldBBox(void);
  
private:
  // Internal variables
  RNArray<RGBDImage *> images;
  RNArray<RGBDSurface *> surfaces;
  char *color_directory;
  char *depth_directory;
  char *texture_directory;
  char *dataset_format;
  R3Box world_bbox;
};



////////////////////////////////////////////////////////////////////////
// Inline functions
////////////////////////////////////////////////////////////////////////

inline R3Box RGBDConfiguration::
WorldBBox(void) const
{
  // Return world bounding box
  if (world_bbox.IsEmpty()) ((RGBDConfiguration *) this)->UpdateWorldBBox();
  return world_bbox;
}



inline int RGBDConfiguration::
NImages(void) const
{
  // Return number of images
  return images.NEntries();
}



inline RGBDImage *RGBDConfiguration::
Image(int k) const
{
  // Return kth image
  return images.Kth(k);
}



inline int RGBDConfiguration::
NSurfaces(void) const
{
  // Return number of surfaces
  return surfaces.NEntries();
}



inline RGBDSurface *RGBDConfiguration::
Surface(int k) const
{
  // Return kth surface
  return surfaces.Kth(k);
}



inline const char *RGBDConfiguration::
ColorDirectory(void) const
{
  // Return name of color directory
  return color_directory;
}



inline const char *RGBDConfiguration::
DepthDirectory(void) const
{
  // Return name of depth directory
  return depth_directory;
}



inline const char *RGBDConfiguration::
TextureDirectory(void) const
{
  // Return name of texture directory
  return texture_directory;
}



inline const char *RGBDConfiguration::
DatasetFormat(void) const
{
  // Return name of dataset format
  return dataset_format;
}


