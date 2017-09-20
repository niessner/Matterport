/* Include file for the R2 texture class */



/* Initialization functions */

int R2InitTexture();
void R2StopTexture();



/* Type definitions */

typedef enum {
  R2_NO_TEXTURE_WRAP = GL_REPEAT,
  R2_REPEAT_TEXTURE_WRAP = GL_REPEAT,
  R2_CLAMP_TEXTURE_WRAP = GL_CLAMP
} R2TextureWrap;

typedef enum {
  R2_NO_TEXTURE_FILTER = GL_NEAREST,
  R2_NEAREST_TEXTURE_FILTER = GL_NEAREST,
  R2_NEAREST_MIPMAP_NEAREST_TEXTURE_FILTER = GL_NEAREST_MIPMAP_NEAREST,
  R2_NEAREST_MIPMAP_LINEAR_TEXTURE_FILTER = GL_NEAREST_MIPMAP_LINEAR,
  R2_LINEAR_TEXTURE_FILTER = GL_LINEAR,
  R2_LINEAR_MIPMAP_NEAREST_TEXTURE_FILTER = GL_LINEAR_MIPMAP_NEAREST,
  R2_LINEAR_MIPMAP_LINEAR_TEXTURE_FILTER = GL_LINEAR_MIPMAP_LINEAR
} R2TextureFilter;

typedef enum {
  R2_NO_TEXTURE_BLEND = GL_MODULATE,
  R2_MODULATE_TEXTURE_BLEND = GL_MODULATE,
  R2_DECAL_TEXTURE_BLEND = GL_DECAL,
  R2_BLEND_TEXTURE_BLEND = GL_BLEND
} R2TextureBlend;



/* Class definition */

class R2Texture {
public:
  // Constructor functions
  R2Texture(void);
  R2Texture(const R2Texture& texture, const char *name = NULL);
  R2Texture(const R2Image *image,
            R2TextureWrap s_wrap = R2_REPEAT_TEXTURE_WRAP,
            R2TextureWrap t_wrap = R2_REPEAT_TEXTURE_WRAP,
            R2TextureFilter min_filter = R2_LINEAR_MIPMAP_LINEAR_TEXTURE_FILTER,
            R2TextureFilter mag_filter = R2_LINEAR_TEXTURE_FILTER,
            R2TextureBlend blend = R2_MODULATE_TEXTURE_BLEND,
            const char *name = NULL);
  R2Texture(const char *filename,
            R2TextureWrap s_wrap = R2_REPEAT_TEXTURE_WRAP,
            R2TextureWrap t_wrap = R2_REPEAT_TEXTURE_WRAP,
            R2TextureFilter min_filter = R2_LINEAR_MIPMAP_LINEAR_TEXTURE_FILTER,
            R2TextureFilter mag_filter = R2_LINEAR_TEXTURE_FILTER,
            R2TextureBlend blend = R2_MODULATE_TEXTURE_BLEND,
            const char *name = NULL);
  ~R2Texture(void);

  // Property functions/operators
  R3Scene *Scene(void) const;
  int SceneIndex(void) const;
  const char *Name(void) const;
  const char *Filename(void) const;
  const R2Image *Image(void) const;
  const R2TextureWrap SWrap(void) const;
  const R2TextureWrap TWrap(void) const;
  const R2TextureFilter MinFilter(void) const;
  const R2TextureFilter MagFilter(void) const;
  const R2TextureBlend Blend(void) const;
  const RNBoolean IsTransparent(void) const;
  const RNBoolean IsLoaded(void) const;
  const RNFlags Flags(void) const;
  const int ID(void) const;

  // Manipulation functions/operations
  void SetName(const char *name);
  void SetFilename(const char *filename);
  void SetImage(const R2Image *image);
  
  // Draw functions/operations
  void Load(void) const;
  void Unload(void) const;
  void Draw(RNBoolean force = FALSE) const;

protected:
  // Upkeep functions/operators
  void Update(void);
  void UpdateFlags(const RNFlags flags);

private:
  friend class R3Scene;
  R3Scene *scene;
  int scene_index;
  char *name;
  char *filename;
  const R2Image *image;
  R2TextureWrap s_wrap, t_wrap;
  R2TextureFilter min_filter, mag_filter;
  R2TextureBlend blend;
  RNFlags flags;
  int id; // 0=notexture, <0=unloaded, >0=loaded
};



/* Flag mask definitions */

#define R2_TEXTURE_FLAGS                0x00000700
#define R2_TEXTURE_TRANSPARENCY_FLAG    0x00000100



/* Public variables */

extern R2Texture R2null_texture;
#define R2default_texture R2null_texture



/* Inline functions */

inline R3Scene *R2Texture::
Scene(void) const
{
  // Return scene 
  return scene;
}



inline int R2Texture::
SceneIndex(void) const
{
  // Return index of texture in scene (can be used with scene->Texture(xxx))
  return scene_index;
}



inline const char *R2Texture::
Name(void) const
{
  // Return name
  return name;
}



inline const char *R2Texture::
Filename(void) const
{
  // Return filename
  return filename;
}



inline const R2Image *R2Texture::
Image(void) const
{
  // Return image
  return image;
}



inline const R2TextureWrap R2Texture::
SWrap(void) const
{
  // Return wrap mode for s coordinate                                                                                
  return s_wrap;
}



inline const R2TextureWrap R2Texture::
TWrap(void) const
{
  // Return wrap mode for t coordinate                                                                                
  return t_wrap;
}



inline const R2TextureFilter R2Texture::
MinFilter(void) const
{
  // Return minification filter                                                                                       
  return min_filter;
}



inline const R2TextureFilter R2Texture::
MagFilter(void) const
{
  // Return magnification filter                                                                                      
  return mag_filter;
}



inline const R2TextureBlend R2Texture::
Blend(void) const
{
  // Return blending mode                                                                                             
  return blend;
}



inline const RNBoolean R2Texture::
IsTransparent(void) const
{
  // Return whether has transparency
  return flags[R2_TEXTURE_TRANSPARENCY_FLAG];
}



inline const RNBoolean R2Texture::
IsLoaded(void) const
{
  // Return whether texture has been loaded to graphics card for rendering
  return (id > 0) ? TRUE : FALSE;
}



inline const RNFlags R2Texture::
Flags(void) const
{
  // Return flags
  return flags;
}



inline const int R2Texture::
ID(void) const
{
  // Return id
  return id;
}




