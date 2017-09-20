/* Include file for the R3 material class */



/* Initialization functions */

int R3InitMaterial();
void R3StopMaterial();



/* Class definition */

class R3Material /* : public R3Attribute */ {
    public:
        // Constructor functions
        R3Material(const char *name = NULL);
        R3Material(const R3Material& material);
        R3Material(const R3Brdf *brdf, const char *name = NULL);
        R3Material(const R2Texture *texture, const char *name = NULL);
        R3Material(const R3Brdf *brdf, const R2Texture *texture, const char *name = NULL);
	~R3Material(void);

        // Material property functions/operations
        R3Scene *Scene(void) const;
        int SceneIndex(void) const;
        const char *Name(void) const;
        const R3Brdf *Brdf(void) const;
        const R2Texture *Texture(void) const;
	const RNFlags& Flags(void) const;
	const RNBoolean IsReflective(void) const;
	const RNBoolean IsAmbient(void) const;
	const RNBoolean IsDiffuse(void) const;
	const RNBoolean IsSpecular(void) const;
	const RNBoolean IsEmissive(void) const;
	const RNBoolean IsShiny(void) const;
	const RNBoolean IsTextured(void) const;
	const RNBoolean IsTransparent(void) const;

	// Manipulation functions/operations
	R3Material& operator=(const R3Material& material);
        void SetBrdf(const R3Brdf *brdf);
        void SetTexture(const R2Texture *texture);
        void SetName(const char *name);

	// Draw functions/operations
        void Load(void) const;
        void Unload(void) const;
        void Draw(RNBoolean force = FALSE) const;

        // Update functions
        void Update(void); 

    private:
        friend class R3Scene;
        R3Scene *scene;
        int scene_index;
        char *name;
        const R3Brdf *brdf;
        const R2Texture *texture;
	RNFlags flags;
};



/* Public variables */

extern R3Material R3null_material;
extern R3Material R3default_material;
extern R3Material R3black_material;
extern R3Material R3red_material;
extern R3Material R3green_material;
extern R3Material R3blue_material;
extern R3Material R3yellow_material;
extern R3Material R3cyan_material;
extern R3Material R3magenta_material;
extern R3Material R3white_material;
extern R3Material R3gray_material;
extern R3Material R3shiny_black_material;
extern R3Material R3shiny_red_material;
extern R3Material R3shiny_green_material;
extern R3Material R3shiny_blue_material;
extern R3Material R3shiny_yellow_material;
extern R3Material R3shiny_cyan_material;
extern R3Material R3shiny_magenta_material;
extern R3Material R3shiny_white_material;
extern R3Material R3shiny_gray_material;



/* Flag mask definitions */

#define R3_MATERIAL_FLAGS		0x00000FFF
#define R3_MATERIAL_BRDF_FLAG           0x00000080
#define R3_MATERIAL_TEXTURE_FLAG	0x00000800
#define R3_MATERIAL_TRANSPARENCY_FLAG	(R3_BRDF_TRANSPARENCY_FLAG | R2_TEXTURE_TRANSPARENCY_FLAG)



/* Inline functions */

inline R3Scene *R3Material::
Scene(void) const
{
  // Return scene 
  return scene;
}



inline int R3Material::
SceneIndex(void) const
{
  // Return index of material in scene (can be used with scene->Material(xxx))
  return scene_index;
}



inline const char *R3Material::
Name(void) const
{
    // Return name
    return name;
}



inline const R3Brdf *R3Material::
Brdf(void) const
{
    // Return brdf
    return brdf;
}



inline const R2Texture *R3Material::
Texture(void) const
{
    // Return texture
    return texture;
}



inline const RNFlags& R3Material::
Flags(void) const
{
    // Return flags
    return flags;
}



inline const RNBoolean R3Material::
IsReflective(void) const
{
    // Return whether has any reflectivity
    return flags[R3_MATERIAL_BRDF_FLAG];
}



inline const RNBoolean R3Material::
IsAmbient(void) const
{
    // Return whether has ambient term
    return flags[R3_BRDF_AMBIENT_FLAG];
}



inline const RNBoolean R3Material::
IsDiffuse(void) const
{
    // Return whether has diffuse term
    return flags[R3_BRDF_DIFFUSE_FLAG];
}



inline const RNBoolean R3Material::
IsSpecular(void) const
{
    // Return whether has specular term
    return flags[R3_BRDF_SPECULAR_FLAG];
}



inline const RNBoolean R3Material::
IsEmissive(void) const
{
    // Return whether has emissive term
    return flags[R3_BRDF_EMISSION_FLAG];
}



inline const RNBoolean R3Material::
IsShiny(void) const
{
    // Return whether has shininess
    return flags[R3_BRDF_SHININESS_FLAG];
}



inline const RNBoolean R3Material::
IsTextured(void) const
{
    // Return whether has texture
    return flags[R3_MATERIAL_TEXTURE_FLAG];
}



inline const RNBoolean R3Material::
IsTransparent(void) const
{
    // Return whether has transparency
    return flags[R3_MATERIAL_TRANSPARENCY_FLAG];
}












