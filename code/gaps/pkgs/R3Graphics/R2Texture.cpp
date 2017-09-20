/* Source file for the R2 texture class */



/* Include files */

#include "R3Graphics.h"



/* Public variables */

R2Texture R2null_texture;



/* Public functions */

int 
R2InitTexture()
{
  /* Return success */
  return TRUE;
}



void 
R2StopTexture()
{
}



R2Texture::
R2Texture(void)
  : scene(NULL),
    scene_index(-1),
    name(NULL),
    filename(NULL),
    image(NULL),
    s_wrap(R2_REPEAT_TEXTURE_WRAP),
    t_wrap(R2_REPEAT_TEXTURE_WRAP),
    min_filter(R2_LINEAR_MIPMAP_LINEAR_TEXTURE_FILTER),
    mag_filter(R2_LINEAR_TEXTURE_FILTER),
    blend(R2_MODULATE_TEXTURE_BLEND),    
    flags(0),
    id(0)
{
}



R2Texture::
R2Texture(const R2Texture& texture, const char *name)
  : scene(NULL),
    scene_index(-1),
    name(NULL),
    filename((texture.filename) ? strdup(texture.filename) : NULL),
    image(texture.image),
    s_wrap(texture.s_wrap),
    t_wrap(texture.t_wrap),
    min_filter(texture.min_filter),
    mag_filter(texture.mag_filter),
    blend(texture.blend),
    flags(texture.flags),
    id((texture.image) ? -1 : 0)
{
    // Set name
    if (name) this->name = strdup(name);

    // Update texture                                                                                                   
    Update();
}



R2Texture::
R2Texture(const R2Image *image,
  R2TextureWrap s_wrap, R2TextureWrap t_wrap,
  R2TextureFilter min_filter, R2TextureFilter mag_filter,
  R2TextureBlend blend,
  const char *name)
  : scene(NULL),
    scene_index(-1),
    name(NULL),
    filename(NULL),
    image(image),
    s_wrap(s_wrap),
    t_wrap(t_wrap),
    min_filter(min_filter),
    mag_filter(mag_filter),
    blend(blend),
    flags(RN_NO_FLAGS),
    id(-1)
{
    // Set name
    if (name) this->name = strdup(name);

    // Update texture
    Update();
}




R2Texture::
R2Texture(const char *filename,
  R2TextureWrap s_wrap, R2TextureWrap t_wrap,
  R2TextureFilter min_filter, R2TextureFilter mag_filter,
  R2TextureBlend blend,
  const char *name)
  : scene(NULL),
    scene_index(-1),
    name(NULL),
    filename((filename) ? strdup(filename) : NULL),
    image(NULL),
    s_wrap(s_wrap),
    t_wrap(t_wrap),
    min_filter(min_filter),
    mag_filter(mag_filter),
    blend(blend),
    flags(RN_NO_FLAGS),
    id(-1)
{
    // Set name
    if (name) this->name = strdup(name);

    // Create image
    image = new R2Image(filename);
    assert(image);

    // Update 
    if (!image) id = 0;
    Update();
}



R2Texture::
~R2Texture(void)
{
  // Unload texture
  if (IsLoaded() > 0) Unload();

  // Remove from scene
  if (scene) scene->RemoveTexture(this);

  // Free name
  if (name) free(name);

  // Free filename
  if (filename) free(filename);
}



void R2Texture::
SetName(const char *name)
{
  // Set name
  if (this->name) free(this->name);
  if (name) this->name = strdup(name);
  else this->name = NULL;
}


  
void R2Texture::
SetFilename(const char *filename)
{
  // Set filename
  if (this->filename) free(this->filename);
  if (filename) this->filename = strdup(filename);
  else this->filename = NULL;
}


  
void R2Texture::
SetImage(const R2Image *image)
{
  // Unload previous texture
  if (IsLoaded()) Unload();

  // Set image
  flags.Remove(R2_TEXTURE_TRANSPARENCY_FLAG);
  if ((image) && ((image->NComponents() == 2) || (image->NComponents() == 4)))
    flags.Add(R2_TEXTURE_TRANSPARENCY_FLAG);
  this->image = image;
  this->id = -1;
}



void R2Texture::
Update (void)
{
  // Update flags
  UpdateFlags(RN_NO_FLAGS);
}



void R2Texture::
UpdateFlags (const RNFlags flags)
{
  // Update flags
  this->flags = flags;
  if ((image) && ((image->NComponents() == 2) || (image->NComponents() == 4)))
    this->flags.Add(R2_TEXTURE_TRANSPARENCY_FLAG);
}



void R2Texture::
Load(void) const
{
  // Check if needs to be loaded
  if (id >= 0) return;

  // Create texture id
  GLuint i;
  glGenTextures(1, &i);
  ((R2Texture *) this)->id = i;
  assert(id > 0);

  // Begin texture definition                                                                                     
  glBindTexture(GL_TEXTURE_2D, id);

  // Define texture parameters
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, s_wrap);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, t_wrap);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, min_filter);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, mag_filter);
  glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, blend);

  // Determine texture format
  GLenum format = GL_RGB;
  if (image->Depth() == 1) format = GL_LUMINANCE;
  else if (image->Depth() == 2) format = GL_LUMINANCE_ALPHA;
  else if (image->Depth() == 3) format = GL_RGB;
  else if (image->Depth() == 4) format = GL_RGBA;
  else RNAbort("Illegal texture image");

#ifdef AT_ONE_POINT_THIS_WAS_NEEDED_FOR_MESA
  // For some reason, gluBuild2DMipmaps segfaults in mesa
  glTexImage2D(GL_TEXTURE_2D, 0, format, image->Width(), image->Height(),
    0, format, GL_UNSIGNED_BYTE, (const unsigned char *) image->Pixels());
#else
  // Define texture
  if ((min_filter == R2_NEAREST_MIPMAP_NEAREST_TEXTURE_FILTER) ||
      (min_filter == R2_NEAREST_MIPMAP_LINEAR_TEXTURE_FILTER) ||
      (min_filter == R2_LINEAR_MIPMAP_NEAREST_TEXTURE_FILTER) ||
      (min_filter == R2_LINEAR_MIPMAP_LINEAR_TEXTURE_FILTER)) {
    gluBuild2DMipmaps(GL_TEXTURE_2D, format, image->Width(), image->Height(),
      format, GL_UNSIGNED_BYTE, (const unsigned char *) image->Pixels());
  }
  else {
    glTexImage2D(GL_TEXTURE_2D, 0, format, image->Width(), image->Height(),
      0, format, GL_UNSIGNED_BYTE, (const unsigned char *) image->Pixels());
  }
#endif

  // End texture definition                                                                                       
  glBindTexture(GL_TEXTURE_2D, 0);
}



void R2Texture::
Unload(void) const
{
  // Check id
  if (!IsLoaded()) return;

  // Delete texture
  GLuint i = id;
  glDeleteTextures(1, &i);

  // Reset texture id
  ((R2Texture *) this)->id = -1;
}



void R2Texture::
Draw(RNBoolean force) const
{
  // Check if same texture as last time
  static const R2Texture *R2current_texture = NULL;
  if (!force && (this == R2current_texture)) return;

  // Check if null texture
  if (id == 0) {
    // Disable texture
    glDisable(GL_TEXTURE_2D);
  }
  else {
    // Load texture
    if (!IsLoaded()) Load();

    // Set texture
    glBindTexture(GL_TEXTURE_2D, id);

    // Enable texture
    glEnable(GL_TEXTURE_2D);
  }

  // Remember current texture
  R2current_texture = this;
}














  
