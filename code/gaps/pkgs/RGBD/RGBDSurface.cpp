////////////////////////////////////////////////////////////////////////
// Source file for RGBDSurface class
////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "RGBD.h"




////////////////////////////////////////////////////////////////////////
// Constructors/destructors
////////////////////////////////////////////////////////////////////////

RGBDSurface::
RGBDSurface(void)
  : configuration(NULL),
    configuration_index(-1),
    channels(),
    width(0), height(0),
    opengl_texture_id(-1),
    texture_filename(NULL),
    mesh_filename(NULL),
    surface_bbox(FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX),
    color_resident_count(0),
    rectangle(NULL),
    mesh(NULL),
    mesh_face_index(NULL)
{
}



RGBDSurface::
RGBDSurface(const char *texture_filename, R3Rectangle *rectangle, RNLength texel_spacing)
  : configuration(NULL),
    configuration_index(-1),
    channels(),
    width(0), height(0),
    opengl_texture_id(-1),
    texture_filename(NULL),
    mesh_filename(NULL),
    surface_bbox(FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX),
    color_resident_count(0),
    rectangle(rectangle),
    mesh(NULL),
    mesh_face_index(NULL)
{
  // Check rectangle
  if (!rectangle) return;
  
  // Set texture filename
  SetTextureFilename(texture_filename);

  // Compute coordinate stuff
  RNLength xlength = 2*rectangle->Radius(0);
  RNLength ylength = 2*rectangle->Radius(1);
  if (texel_spacing <= 0) texel_spacing = 0.01;
  width = (int) (xlength / texel_spacing + 0.5);
  height = (int) (ylength / texel_spacing + 0.5);
  if (width <= 0) width = 1;
  if (height <= 0) height = 1;

  // Compute bounding box of surface coordinates
  surface_bbox.Reset(R2zero_point, R2Point(xlength, ylength));
}



RGBDSurface::
RGBDSurface(const char *texture_filename, R3Mesh *mesh, RNLength texel_spacing)
  : configuration(NULL),
    configuration_index(-1),
    channels(),
    width(0), height(0),
    opengl_texture_id(-1),
    texture_filename(NULL),
    mesh_filename(NULL),
    surface_bbox(FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX),
    color_resident_count(0),
    rectangle(NULL),
    mesh(mesh),
    mesh_face_index(NULL)
{
  // Check mesh
  if (mesh->NFaces() == 0) return;
  RNArea mesh_area = mesh->Area();
  if (RNIsZero(mesh_area)) {
    fprintf(stderr, "Warning: RGBDSurface mesh has zero area\n");
    return;
  }

  // Compute bounding box of surface coordinates
  for (int i = 0; i < mesh->NVertices(); i++) {
    R3MeshVertex *vertex = mesh->Vertex(i);
    const R2Point& texcoords = mesh->VertexTextureCoords(vertex);
    surface_bbox.Union(texcoords);
  }

  // Check surface_bbox
  if (RNIsZero(surface_bbox.Area())) {
    fprintf(stderr, "Warning: RGBDSurface mesh does not have valid texture coordinates\n");
    surface_bbox.Reset(R2zero_point, R2ones_point);
  }

  // Compute width and height from texel_spacing
  if (texel_spacing <= 0) texel_spacing = 0.01;
  width = (int) (surface_bbox.XLength() / texel_spacing + 0.5);
  height = (int) (surface_bbox.YLength() / texel_spacing + 0.5);
  if (width <= 0) width = 1;
  if (height <= 0) height = 1;

  // Set texture filename
  SetTextureFilename(texture_filename);
}



RGBDSurface::
~RGBDSurface(void)
{
  // Delete opengl texture
  InvalidateOpenGL();
  
  // Remove from configuration
  if (configuration) {
    configuration->RemoveSurface(this);
  }

  // Delete channels
  for (int i = 0; i < channels.NEntries(); i++) {
    if (channels[i]) delete channels[i];
  }

  // Delete filenames
  if (texture_filename) free(texture_filename);
  if (mesh_filename) free(mesh_filename);

  // Delete mesh face index
  if (mesh_face_index) delete mesh_face_index;

  // Delete geometry stuff
  if (rectangle) delete rectangle;
  if (mesh) delete mesh;
}



////////////////////////////////////////////////////////////////////////
// Texture pixel property access functions
////////////////////////////////////////////////////////////////////////

RNRgb RGBDSurface::
TexelColor(int ix, int iy) const
{
  // Return color of pixel
  if (NChannels() < 3) return RNblack_rgb;
  RNScalar r = TexelChannelValue(ix, iy, RGBD_RED_CHANNEL);
  RNScalar g = TexelChannelValue(ix, iy, RGBD_GREEN_CHANNEL);
  RNScalar b = TexelChannelValue(ix, iy, RGBD_BLUE_CHANNEL);
  return RNRgb(r, g, b);
}



R3Point RGBDSurface::
TexelWorldPosition(int ix, int iy) const
{
  // Return position in world coordinates
  R2Point texture_position(ix + 0.5, iy + 0.5);
  return TexelWorldPosition(texture_position);
}



R3Vector RGBDSurface::
TexelWorldNormal(int ix, int iy) const
{
  // Return surface normal
  R2Point texture_position(ix + 0.5, iy + 0.5);
  return TexelWorldNormal(texture_position);
}



////////////////////////////////////////////////////////////////////////
// Texture point property access functions
////////////////////////////////////////////////////////////////////////

RNRgb RGBDSurface::
TexelColor(const R2Point& texture_position) const
{
  // Return color of pixel
  if (NChannels() < 3) return RNblack_rgb;
  RNScalar r = TexelChannelValue(texture_position, RGBD_RED_CHANNEL);
  RNScalar g = TexelChannelValue(texture_position, RGBD_GREEN_CHANNEL);
  RNScalar b = TexelChannelValue(texture_position, RGBD_BLUE_CHANNEL);
  return RNRgb(r, g, b);
}



R3Point RGBDSurface::
TexelWorldPosition(const R2Point& texture_position) const
{
  // Return position in world coordinates
  R3Point world_position(0,0,0);
  RGBDTransformTextureToWorld(texture_position, world_position, this);
  return world_position;
}



R3Vector RGBDSurface::
TexelWorldNormal(const R2Point& texture_position) const
{
  // Return surface normal
  if (rectangle) {
    return rectangle->Normal();
  }
  else if (mesh) {
    double barycentrics[3];
    R3MeshFace *face = SearchMeshFaceIndex(texture_position, barycentrics);
    if (!face) return R3zero_vector;
    else return mesh->FaceNormal(face);
  }
  else {
    return R3zero_vector;
  }
}



////////////////////////////////////////////////////////////////////////
// Surface point property access functions
////////////////////////////////////////////////////////////////////////

RNRgb RGBDSurface::
SurfelColor(const R2Point& surface_position) const
{
  // Return color of texture at surface point
  R2Point texture_position;
  if (!TransformSurfaceToTexture(surface_position, texture_position)) return RNblack_rgb;
  return TexelColor(texture_position);
}



R3Point RGBDSurface::
SurfelWorldPosition(const R2Point& surface_position) const
{
  // Return world position at surface point
  R3Point world_position;
  if (!TransformSurfaceToWorld(surface_position, world_position)) return R3zero_point;
  return world_position;
}



R3Vector RGBDSurface::
SurfelWorldNormal(const R2Point& surface_position) const
{
  // Return normal at surface point
  if (rectangle) return rectangle->Normal();
  else if (mesh) return R3negy_vector; 
  return R3zero_vector;
}



RNScalar RGBDSurface::
SurfelChannelValue(const R2Point& surface_position, int channel_index) const
{
  // Return
  R2Point texture_position;
  if (!TransformSurfaceToTexture(surface_position, texture_position)) return 0;
  return TexelChannelValue(texture_position, channel_index);
}



////////////////////////////////////////////////////////////////////////
// Manipulation functions
////////////////////////////////////////////////////////////////////////

void RGBDSurface::
SetTexelColor(int ix, int iy, const RNRgb& color)
{
  // Set pixel red, green, and blue
  SetTexelChannelValue(ix, iy, RGBD_RED_CHANNEL, color.R());
  SetTexelChannelValue(ix, iy, RGBD_GREEN_CHANNEL, color.G());
  SetTexelChannelValue(ix, iy, RGBD_BLUE_CHANNEL, color.B());

  // Invalidate opengl
  InvalidateOpenGL();
}



void RGBDSurface::
SetTexelChannelValue(int ix, int iy, int channel_index, RNScalar value)
{
  // Check channel
  if (!channels[channel_index]) {
    fprintf(stderr, "RGBD Channel is not resident in memory -- cannot get value\n");
    return;
  }
  
  // Set channel value at pixel
  if ((channel_index < 0) || (channel_index >= channels.NEntries())) return;
  if ((ix < 0) || (ix >= channels[channel_index]->XResolution())) return;
  if ((iy < 0) || (ix >= channels[channel_index]->YResolution())) return;
  channels[channel_index]->SetGridValue(ix, iy, value);

  // Invalidate opengl
  if ((channel_index >= RGBD_RED_CHANNEL) && (channel_index <= RGBD_BLUE_CHANNEL)) InvalidateOpenGL();
}



void RGBDSurface::
SetChannel(int channel_index, const R2Grid& image)
{
  // Check if color channels are resident
  if ((channel_index >= RGBD_RED_CHANNEL) && (channel_index <= RGBD_BLUE_CHANNEL)) {
    if (color_resident_count == 0) {
      fprintf(stderr, "Unable to set channel -- it has not been created\n");
      return;
    }
  }

  // Set channel
  *(channels[channel_index]) = image;

  // Set width and height
  this->width = image.XResolution();
  this->height = image.YResolution();

  // Invalidate opengl
  if ((channel_index >= RGBD_RED_CHANNEL) && (channel_index <= RGBD_BLUE_CHANNEL)) InvalidateOpenGL();
}



void RGBDSurface::
SetColorChannels(const R2Image& image)
{
  // Check if color channels are resident
  if (color_resident_count == 0) {
    fprintf(stderr, "Unable to set color channels -- they have not been created\n");
    return;
  }

  // Copy color values
  for (int iy = 0; iy < image.Height(); iy++) {
    for (int ix = 0; ix < image.Width(); ix++) {
      RNRgb color = image.PixelRGB(ix, iy);
      channels[RGBD_RED_CHANNEL]->SetGridValue(ix, iy, color.R());      
      channels[RGBD_GREEN_CHANNEL]->SetGridValue(ix, iy, color.G());      
      channels[RGBD_BLUE_CHANNEL]->SetGridValue(ix, iy, color.B());
    }
  }

  // Set width and height
  this->width = image.Width();
  this->height = image.Height();

  // Invalidate opengl
  InvalidateOpenGL();
  
}



void RGBDSurface::
SetTextureFilename(const char *filename)
{
  // Set filename
  if (texture_filename) free(texture_filename);
  if (filename && strcmp(filename, "-")) texture_filename = strdup(filename);
  else texture_filename = NULL;
}



void RGBDSurface::
SetMeshFilename(const char *filename)
{
  // Set filename
  if (mesh_filename) free(mesh_filename);
  if (filename && strcmp(filename, "-")) mesh_filename = strdup(filename);
  else mesh_filename = NULL;
}



void RGBDSurface::
SetWorldTexelSpacing(RNLength texel_spacing)
{
  // Compute new width and height
  if (texel_spacing > 0) {
    width = surface_bbox.XLength() / texel_spacing;
    height = surface_bbox.YLength() / texel_spacing;
  }
  
  // Resample channels
  for (int i = 0; i < NChannels(); i++) {
    R2Grid *channel = Channel(i);
    if (channel) channel->Resample(width, height);
  }

  // Invalidate mesh face index
  InvalidateMeshFaceIndex();
}



void RGBDSurface::
Transform(const R3Transformation& transformation)
{
  // Update rectangle
  if (rectangle) rectangle->Transform(transformation);
  else if (mesh) mesh->Transform(transformation);

  // Invalidate the configuration's bounding box
  if (configuration) configuration->InvalidateWorldBBox();
}



////////////////////////////////////////////////////////////////////////
// Transformation functions
////////////////////////////////////////////////////////////////////////

int RGBDSurface::
TransformTextureToSurface(const R2Point& texture_position, R2Point& surface_position) const
{
  // Transform from position in texture coordinates to surface coordinates
  if ((width == 0) || (height == 0)) return 0;
  RNScalar xlength = surface_bbox.XLength();
  RNScalar ylength = surface_bbox.YLength();
  surface_position[0] = surface_bbox[0][0] + xlength * texture_position[0] / width;
  surface_position[1] = surface_bbox[0][1] + ylength * texture_position[1] / height;
  return 1;
}



int RGBDSurface::
TransformSurfaceToTexture(const R2Point& surface_position, R2Point& texture_position) const
{
  // Transform from position in surface coordinates to texture coordinates
  RNScalar xlength = surface_bbox.XLength();
  RNScalar ylength = surface_bbox.YLength();
  if ((xlength == 0) || (ylength == 0)) return 0;
  texture_position[0] = width  * (surface_position[0] - surface_bbox[0][0]) / xlength;
  texture_position[1] = height * (surface_position[1] - surface_bbox[0][1]) / ylength;
  return 1;
}



int RGBDSurface::
TransformSurfaceToWorld(const R2Point& surface_position, R3Point& world_position) const
{
  // Transform from position in surface coordinates to world coordinates
  if (rectangle) {
    world_position = rectangle->Corner(0,0);
    world_position += surface_position[0] * rectangle->Axis(0);
    world_position += surface_position[1] * rectangle->Axis(1);
  }
  else if (mesh) {
    R2Point texture_position;
    double barycentrics[3];
    if (!TransformSurfaceToTexture(surface_position, texture_position)) return 0;
    R3MeshFace *face = SearchMeshFaceIndex(texture_position, barycentrics);
    if (face) world_position = mesh->FacePoint(face, barycentrics);
    else return 0;
  }

  // Return success
  return 1;
}



int RGBDSurface::
TransformWorldToSurface(const R3Point& world_position, R2Point& surface_position) const
{
  // Transform from position in world coordinates to surface coordinates
  if (rectangle) {
    R3Point origin = rectangle->Corner(0,0);
    surface_position[0] = (world_position - origin).Dot(rectangle->Axis(0));;
    surface_position[1] = (world_position - origin).Dot(rectangle->Axis(1));;
  }
  else {
    return 0;
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Draw functions
////////////////////////////////////////////////////////////////////////

void RGBDSurface::
Draw(int color_scheme) const
{
  // Draw surface
  DrawTexture(color_scheme);
}



void RGBDSurface::
DrawFaces(int color_scheme) const
{
  // Enable lighting and material
  if (color_scheme == RGBD_RENDER_COLOR_SCHEME) {
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_LIGHTING);
  }

  // Load color
  LoadColor(color_scheme);

  // Draw faces
  if (rectangle) {
    rectangle->Draw();
  }
  else if (mesh) {
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < mesh->NFaces(); i++) {
      R3MeshFace *face = mesh->Face(i);
      R3MeshVertex *v0 = mesh->VertexOnFace(face, 0);
      R3MeshVertex *v1 = mesh->VertexOnFace(face, 1);
      R3MeshVertex *v2 = mesh->VertexOnFace(face, 2);
      if (color_scheme == RGBD_PHOTO_COLOR_SCHEME) {
        RNLoadRgb(mesh->VertexColor(v0));
        R3LoadPoint(mesh->VertexPosition(v0));
        RNLoadRgb(mesh->VertexColor(v1));
        R3LoadPoint(mesh->VertexPosition(v1));
        RNLoadRgb(mesh->VertexColor(v2));
        R3LoadPoint(mesh->VertexPosition(v2));
      }
      else if (color_scheme == RGBD_RENDER_COLOR_SCHEME) {
        R3LoadNormal(mesh->FaceNormal(face));
        R3LoadPoint(mesh->VertexPosition(v0));
        R3LoadPoint(mesh->VertexPosition(v1));
        R3LoadPoint(mesh->VertexPosition(v2));
      }
      else {
        R3LoadPoint(mesh->VertexPosition(v0));
        R3LoadPoint(mesh->VertexPosition(v1));
        R3LoadPoint(mesh->VertexPosition(v2));
      }
    }
    glEnd();
  }

  // Disable lighting and material
  if (color_scheme == RGBD_RENDER_COLOR_SCHEME) {
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_LIGHTING);
  }
}



void RGBDSurface::
DrawEdges(int color_scheme) const
{
  // Load color
  if (color_scheme == RGBD_NO_COLOR_SCHEME) LoadColor(color_scheme);
  else if (color_scheme == RGBD_HIGHLIGHT_COLOR_SCHEME) LoadColor(color_scheme);
  else LoadColor(RGBD_INDEX_COLOR_SCHEME);

  // Draw rectangle
  if (rectangle) rectangle->Outline();
  else if (mesh) mesh->DrawEdges();
}



void RGBDSurface::
DrawTexture(int color_scheme) const
{
  // Check channels
  if (NChannels() == 0) return;
  
  // Update/select opengl texture
  if (opengl_texture_id <= 0) ((RGBDSurface *) this)->UpdateOpenGL();
  if (opengl_texture_id <= 0) return;
  glBindTexture(GL_TEXTURE_2D, opengl_texture_id);
  glEnable(GL_TEXTURE_2D);

  // Load color
  LoadColor(color_scheme);

  // Draw rectangle
  if (rectangle) {
    glBegin(GL_QUADS);
    R3LoadTextureCoords(0.0, 0.0);
    R3LoadPoint(rectangle->Corner(0,0));
    R3LoadTextureCoords(1.0, 0.0);
    R3LoadPoint(rectangle->Corner(1,0));
    R3LoadTextureCoords(1.0, 1.0);
    R3LoadPoint(rectangle->Corner(1,1));
    R3LoadTextureCoords(0.0, 1.0);
    R3LoadPoint(rectangle->Corner(0,1));
    glEnd();
  }
  else if (mesh) {
    if ((width > 0) && (height > 0)) {
      glBegin(GL_TRIANGLES);
      for (int i = 0; i < mesh->NFaces(); i++) {
        R3MeshFace *face = mesh->Face(i);
        for (int j = 0; j < 3; j ++) {
          R3MeshVertex *vertex = mesh->VertexOnFace(face, j);
          const R3Point& world_position = mesh->VertexPosition(vertex);
          R2Point texture_position, surface_position = mesh->VertexTextureCoords(vertex);
          TransformSurfaceToTexture(surface_position, texture_position);
          texture_position[0] = texture_position[0] / width;
          texture_position[1] = texture_position[1] / height;
          R3LoadTextureCoords(texture_position);
          R3LoadPoint(world_position);
        }
      }
      glEnd();
    }
  }

  // Unselect opengl texture
  glDisable(GL_TEXTURE_2D);
  glBindTexture(GL_TEXTURE_2D, 0);
}



void RGBDSurface::
LoadColor(int color_scheme) const
{
  // Check color scheme
  if (color_scheme == RGBD_INDEX_COLOR_SCHEME) {
    // Load color encoding index
    int k = 255.0 * (ConfigurationIndex()+1) / configuration->NSurfaces();
    unsigned char r = k;
    unsigned char g = 0;
    unsigned char b = 0;
    glColor3ub(r, g, b);
  }
  else if (color_scheme == RGBD_PHOTO_COLOR_SCHEME) {
    // Load white 
    glColor3d(1.0, 1.0, 1.0);
  }
  else if (color_scheme == RGBD_RENDER_COLOR_SCHEME) {
    // Load white 
    glColor3d(1.0, 1.0, 1.0);
  }
  else if (color_scheme == RGBD_HIGHLIGHT_COLOR_SCHEME) {
    // Load highlight color
    glColor3d(1.0, 1.0, 0.0);
  }
}



////////////////////////////////////////////////////////////////////////
// Read/write/release functions
////////////////////////////////////////////////////////////////////////

int RGBDSurface::
ReadChannels(void)
{
  // Read all channels
  if (!ReadColorChannels()) return 0;
  return 1;
}



int RGBDSurface::
ReadColorChannels(void)
{
  // Initialize image
  R2Image rgb_image(width, height);

  // Check filename
  if (texture_filename) {
    // Get full filename
    char full_filename[4096];
    const char *dirname = (configuration) ? configuration->TextureDirectory() : NULL;
    if (dirname) sprintf(full_filename, "%s/%s", dirname, texture_filename);
    else sprintf(full_filename, "%s", texture_filename);
    if (!RNFileExists(full_filename)) return 0;
  
    // Read color image
    if (!rgb_image.Read(full_filename)) return 0;

    // Check image dimensions
    if ((rgb_image.Width() != width) || (rgb_image.Height() != height)) {
      fprintf(stderr, "Mismatching image dimensions (%d,%d) vs (%d,%d) in %s\n",
        rgb_image.Width(), rgb_image.Height(), width, height, texture_filename);
      // return 0;
    }
  }

  // Create color channels
  return CreateColorChannels(rgb_image);
}



int RGBDSurface::
WriteChannels(void)
{
  // Write all channels
  if (!WriteColorChannels()) return 0;
  return 1;
}



int RGBDSurface::
WriteColorChannels(void)
{
  // Check filename
  if (NChannels() <= RGBD_BLUE_CHANNEL) return 0;

  // Create texture filename, if there is none
  if (!texture_filename) {
    char buffer[4096];
    sprintf(buffer, "t%06d.jpg", ConfigurationIndex());
    SetTextureFilename(buffer);
  }

  // Get full filename
  char full_filename[4096];
  const char *dirname = (configuration) ? configuration->TextureDirectory() : NULL;
  if (dirname) sprintf(full_filename, "%s/%s", dirname, texture_filename);
  else sprintf(full_filename, "%s", texture_filename);

  // Compute color image
  R2Image rgb_image(width, height, 3);
  for (int iy = 0; iy < height; iy++) {
    for (int ix = 0; ix < width; ix++) {
      RNRgb rgb = TexelColor(ix, iy);
      rgb_image.SetPixelRGB(ix, iy, rgb);
    }
  }

  // Write color image
  if (!rgb_image.Write(full_filename)) return 0;

  // Return success
  return 1;
}



int RGBDSurface::
ReleaseChannels(void)
{
  // Release all channels
  if (!ReleaseColorChannels()) return 0;
  return 1;
}



int RGBDSurface::
ReleaseColorChannels(void)
{
  // Check/update read count
  if (--color_resident_count > 0) return 1;

  // Write color channel before release???
  // WriteColorChannels();
  
  // Delete color channels
  for (int channel_index = RGBD_RED_CHANNEL; channel_index <= RGBD_BLUE_CHANNEL; channel_index++) {
    if (NChannels() <= channel_index) break;
    if (!channels[channel_index]) continue;
    delete channels[channel_index];
    channels[channel_index] = NULL;
  }

  // Return success;
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Channel creation functions
////////////////////////////////////////////////////////////////////////

int RGBDSurface::
CreateColorChannels(const R2Image& image)
{
  // Check/update read count
  if (color_resident_count++ > 0) return 1;

  // Create color channels
  while (channels.NEntries() <= RGBD_BLUE_CHANNEL) channels.Insert(NULL);
  if (!channels[RGBD_RED_CHANNEL]) channels[RGBD_RED_CHANNEL] = new R2Grid(image.Width(), image.Height());
  if (!channels[RGBD_GREEN_CHANNEL]) channels[RGBD_GREEN_CHANNEL] = new R2Grid(image.Width(), image.Height());
  if (!channels[RGBD_BLUE_CHANNEL]) channels[RGBD_BLUE_CHANNEL] = new R2Grid(image.Width(), image.Height());

  // Initialize color channels
  SetColorChannels(image);

  // Set width and height
  this->width = image.Width();
  this->height = image.Height();

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Update functions
////////////////////////////////////////////////////////////////////////

void RGBDSurface::
InvalidateOpenGL(void) 
{
  // Delete opengl texture
  if (opengl_texture_id > 0) {
    GLuint i = opengl_texture_id;
    glDeleteTextures(1, &i);
  }

  // Reset opengl identifier
  opengl_texture_id = -1;
}


  
void RGBDSurface::
UpdateOpenGL(void) 
{
  // Check identifier
  if (opengl_texture_id > 0) return;
  
  // Create identifier
  GLuint identifier;
  glGenTextures(1, &identifier);

  // Allocate pixels
  unsigned char *pixels = new unsigned char [ width * height * 4];

  // Fill pixels
  unsigned char *pixelsp = pixels;
  for (int iy = 0; iy < height; iy++) {
    for (int ix = 0; ix < width; ix++) {
      RNRgb color = TexelColor(ix, iy);
      *pixelsp++ = (color.R() < 1.0) ? color.R() * 255 : 255;
      *pixelsp++ = (color.G() < 1.0) ? color.G() * 255 : 255;
      *pixelsp++ = (color.B() < 1.0) ? color.B() * 255 : 255;
      *pixelsp++ = 255;
    }
  }

  // Define texture
  glBindTexture(GL_TEXTURE_2D, identifier);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height,
    0, GL_RGBA, GL_UNSIGNED_BYTE, pixels);

  // Delete the pixels
  delete [] pixels;

  //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
  // gluBuild2DMipmaps(GL_TEXTURE_2D, 3, rgb_image.Width(), rgb_image.Height(),
  //   GL_RGB, GL_UNSIGNED_BYTE, (const unsigned char *) rgb_image.Pixels());
 
  // Remember identifier
  opengl_texture_id = identifier;
}



////////////////////////////////////////////////////////////////////////
// Mesh search functions
////////////////////////////////////////////////////////////////////////

R3MeshFace *RGBDSurface::
SearchMeshFaceIndex(const R2Point& texture_position, double *barycentrics) const
{
  // Check mesh
  if (!mesh) return NULL;

  // Update mesh face index
  if (!mesh_face_index) ((RGBDSurface *) this)->UpdateMeshFaceIndex();
  if (!mesh_face_index) return NULL;

  // Check texture coordinates
  int texture_ix = (int) (texture_position.X() + 0.5);
  int texture_iy = (int) (texture_position.Y() + 0.5);
  if ((texture_ix < 0) || (texture_ix >= mesh_face_index->XResolution())) return NULL;
  if ((texture_iy < 0) || (texture_iy >= mesh_face_index->YResolution())) return NULL;

  // Look up face index
  RNScalar face_index_value = mesh_face_index->GridValue(texture_ix, texture_iy);
  if (face_index_value == R2_GRID_UNKNOWN_VALUE) return NULL;
  int face_index = (int) (face_index_value + 0.5);

  // Get face
  if ((face_index < 0) || (face_index >= mesh->NFaces())) return NULL;
  R3MeshFace *face = mesh->Face(face_index);

  // Get barycentric coordinates
  if (barycentrics) {
    // Get useful variables
    R3MeshVertex *va = mesh->VertexOnFace(face, 0);
    R3MeshVertex *vb = mesh->VertexOnFace(face, 1);
    R3MeshVertex *vc = mesh->VertexOnFace(face, 2);
    R2Point ta = mesh->VertexTextureCoords(va);
    R2Point tb = mesh->VertexTextureCoords(vb);
    R2Point tc = mesh->VertexTextureCoords(vc);
    R2Vector vab = tb - ta;
    R2Vector vac = tc - ta;
    RNScalar areaABC = vab[0]*vac[1] - vab[1]*vac[0];
    if (RNIsZero(areaABC)) {
      barycentrics[0] = 0;
      barycentrics[1] = 0;
      barycentrics[2] = 0;
      return NULL;
    }

    // Get/check surface position
    R2Point surface_position;
    if (!RGBDTransformTextureToSurface(texture_position, surface_position, this)) return NULL;

    // Compute barycentric coordinates
    R2Vector vpb = tb - surface_position;
    R2Vector vpc = tc - surface_position;
    R2Vector vpa = ta - surface_position;
    RNScalar areaPBC = vpb[0]*vpc[1] - vpb[1]*vpc[0];
    RNScalar areaPCA = vpc[0]*vpa[1] - vpc[1]*vpa[0];
    barycentrics[0] = areaPBC / areaABC;
    barycentrics[1] = areaPCA / areaABC;
    barycentrics[2] = 1 - barycentrics[0] - barycentrics[1];

    // I guess it's OK if the barycentric coordinates are (barely) outside face (since face_index grid is pixelated)
    // if (RNIsNegative(barycentrics[0]) || RNIsGreater(barycentrics[0], 1.0)) return NULL; 
    // if (RNIsNegative(barycentrics[1]) || RNIsGreater(barycentrics[1], 1.0)) return NULL; 
    // if (RNIsNegative(barycentrics[2]) || RNIsGreater(barycentrics[2], 1.0)) return NULL; 
  }

  // Return face
  return face;
}
  


void RGBDSurface::
UpdateMeshFaceIndex(void)
{
  // Check everything
  if (!mesh) return;
  if (width == 0) return;
  if (height == 0) return;
  if (RNIsZero(surface_bbox.Area())) return;
  
  // Create face index
  mesh_face_index = new R2Grid(width, height, surface_bbox);

  // Update mesh face index grid
  mesh_face_index->Clear(R2_GRID_UNKNOWN_VALUE);
  for (int i = 0; i < mesh->NFaces(); i++) {
    R3MeshFace *face = mesh->Face(i);
    R3MeshVertex *v0 = mesh->VertexOnFace(face, 0);
    R3MeshVertex *v1 = mesh->VertexOnFace(face, 1);
    R3MeshVertex *v2 = mesh->VertexOnFace(face, 2);
    R2Point s0 = mesh->VertexTextureCoords(v0);
    R2Point s1 = mesh->VertexTextureCoords(v1);
    R2Point s2 = mesh->VertexTextureCoords(v2);
    mesh_face_index->RasterizeWorldTriangle(s0, s1, s2, i, R2_GRID_REPLACE_OPERATION);
  }

  // For debugging
  // mesh_face_index->WriteFile("mesh_face_index.pfm");
}
  


void RGBDSurface::
InvalidateMeshFaceIndex(void)
{
  // Delete mesh face index
  if (!mesh_face_index) return;
  delete mesh_face_index;
  mesh_face_index = NULL;
}
