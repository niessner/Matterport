/* Include file for the R3 brdf class */



/* Initialization functions */

int R3InitBrdf();
void R3StopBrdf();



/* Class definition */

class R3Brdf {
    public:
        // Constructor functions
	R3Brdf(const char *name = NULL);
        R3Brdf(const R3Brdf& brdf, const char *name = NULL);
        R3Brdf(const RNRgb& rgb, 
               RNScalar shininess = 0.0,
               RNScalar opacity = 1.0,
               RNScalar indexofrefraction = 1.0,
               const char *name = NULL);
        R3Brdf(RNScalar red, RNScalar green, RNScalar blue,
               RNScalar shininess = 0.0,
               RNScalar opacity = 1.0,
               RNScalar indexofrefraction = 1.0,
               const char *name = NULL);
        R3Brdf(const RNRgb& ambient,
               const RNRgb& diffuse,
               const RNRgb& specular,
               const RNRgb& emission,
               RNScalar shininess = 0.0,
               RNScalar opacity = 1.0,
               RNScalar indexofrefraction = 1.0,
               const char *name = NULL);
        R3Brdf(const RNRgb& ambient,
               const RNRgb& diffuse,
               const RNRgb& specular,
               const RNRgb& transmission,
               const RNRgb& emission,
               RNScalar shininess = 0.0,
               RNScalar indexofrefraction = 1.0,
               const char *name = NULL);
        virtual ~R3Brdf(void);
  
	// Property functions/operators
        R3Scene *Scene(void) const;
        int SceneIndex(void) const;
        const char *Name(void) const;
  	const RNRgb& Ambient(void) const;
  	const RNRgb& Diffuse(void) const;
  	const RNRgb& Specular(void) const;
  	const RNRgb& Transmission(void) const;
  	const RNRgb& Emission(void) const;
        const RNScalar Shininess(void) const;
        const RNScalar Opacity(void) const;
        const RNScalar IndexOfRefraction(void) const;
	const int IsAmbient(void) const;
	const int IsDiffuse(void) const;
	const int IsSpecular(void) const;
	const int IsTransparent(void) const;
	const int IsEmissive(void) const;
	const int IsShiny(void) const;
	const RNFlags Flags(void) const;
	const int ID(void) const;

	// Manipulation functions/operations
        void SetName(const char *name);
        void SetAmbient(const RNRgb& rgb);
  	void SetDiffuse(const RNRgb& rgb);
  	void SetSpecular(const RNRgb& rgb);
  	void SetTransmission(const RNRgb& rgb);
  	void SetEmission(const RNRgb& rgb);
        void SetShininess(RNScalar shininess);
        void SetOpacity(RNScalar opacity);
        void SetIndexOfRefraction(RNScalar indexofrefraction);

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
	RNRgb ambient;
	RNRgb diffuse;
	RNRgb specular;
        RNRgb transmission;
	RNRgb emission;
	RNScalar shininess;
        RNScalar indexofrefraction;
	RNFlags flags;
        int id;
};



/* Public variables */

extern R3Brdf R3null_brdf;
extern R3Brdf R3default_brdf;
extern R3Brdf R3black_brdf;
extern R3Brdf R3red_brdf;
extern R3Brdf R3green_brdf;
extern R3Brdf R3blue_brdf;
extern R3Brdf R3yellow_brdf;
extern R3Brdf R3cyan_brdf;
extern R3Brdf R3magenta_brdf;
extern R3Brdf R3white_brdf;
extern R3Brdf R3gray_brdf;
extern R3Brdf R3shiny_black_brdf;
extern R3Brdf R3shiny_red_brdf;
extern R3Brdf R3shiny_green_brdf;
extern R3Brdf R3shiny_blue_brdf;
extern R3Brdf R3shiny_yellow_brdf;
extern R3Brdf R3shiny_cyan_brdf;
extern R3Brdf R3shiny_magenta_brdf;
extern R3Brdf R3shiny_white_brdf;
extern R3Brdf R3shiny_gray_brdf;



/* Flag mask definitions */

#define R3_BRDF_FLAGS     		0x0000007F
#define R3_BRDF_COLOR_FLAGS             0x0000000F
#define R3_BRDF_AMBIENT_FLAG            0x00000001
#define R3_BRDF_DIFFUSE_FLAG            0x00000002
#define R3_BRDF_SPECULAR_FLAG           0x00000004
#define R3_BRDF_TRANSPARENCY_FLAG       0x00000020
#define R3_BRDF_EMISSION_FLAG           0x00000008
#define R3_BRDF_SHININESS_FLAG          0x00000010
#define R3_BRDF_LOADED_FLAG             0x00000040



/* Inline functions */

inline R3Scene *R3Brdf::
Scene(void) const
{
  // Return scene 
  return scene;
}



inline int R3Brdf::
SceneIndex(void) const
{
  // Return index of brdf in scene (can be used with scene->Brdf(xxx))
  return scene_index;
}



inline const char *R3Brdf::
Name(void) const
{
  // Return name
  return name;
}



inline const RNRgb& R3Brdf::
Ambient(void) const
{
    // Return ambient rgb
    return ambient;
}



inline const RNRgb& R3Brdf::
Diffuse(void) const
{
    // Return diffuse rgb
    return diffuse;
}



inline const RNRgb& R3Brdf::
Specular(void) const
{
    // Return specular rgb
    return specular;
}



inline const RNRgb& R3Brdf::
Transmission(void) const
{
    // Return transmission rgb
    return transmission;
}



inline const RNRgb& R3Brdf::
Emission(void) const
{
    // Return emission rgb
    return emission;
}



inline const RNScalar R3Brdf::
Shininess(void) const
{
    // Return shininess
    return shininess;
}



inline const RNScalar R3Brdf::
Opacity(void) const
{
    // Return opacity
    return 1.0 - Transmission().Luminance();
}



inline const RNScalar R3Brdf::
IndexOfRefraction(void) const
{
    // Return index of refraction
    return indexofrefraction;
}



inline const RNFlags R3Brdf::
Flags(void) const
{
    // Return flags
    return flags;
}



inline const int R3Brdf::
IsAmbient(void) const
{
    // Return whether has ambient term
    return flags[R3_BRDF_AMBIENT_FLAG];
}



inline const int R3Brdf::
IsDiffuse(void) const
{
    // Return whether has diffuse term
    return flags[R3_BRDF_DIFFUSE_FLAG];
}



inline const int R3Brdf::
IsSpecular(void) const
{
    // Return whether has specular term
    return flags[R3_BRDF_SPECULAR_FLAG];
}



inline const int R3Brdf::
IsTransparent(void) const
{
    // Return whether has transparency
    return flags[R3_BRDF_TRANSPARENCY_FLAG];
}



inline const int R3Brdf::
IsEmissive(void) const
{
    // Return whether has emissive term
    return flags[R3_BRDF_EMISSION_FLAG];
}



inline const int R3Brdf::
IsShiny(void) const
{
    // Return whether has shininess
    return flags[R3_BRDF_SHININESS_FLAG];
}



inline const int R3Brdf::
ID(void) const
{
    // Return id
    return id;
}



inline void R3Brdf::
Load(void) const
{
}



inline void R3Brdf::
Unload(void) const
{
}



