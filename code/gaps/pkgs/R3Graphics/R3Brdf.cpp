/* Source file for the R3 brdf class */



/* Include files */

#include "R3Graphics.h"



/* Public variables */

R3Brdf R3null_brdf;
R3Brdf R3default_brdf(
    RNRgb(0.2, 0.2, 0.2), RNRgb(0.8, 0.8, 0.8), 
    RNRgb(0.0, 0.0, 0.0), RNRgb(0.0, 0.0, 0.0), 0.2, 1.0, 1.0);
R3Brdf R3black_brdf(0.0, 0.0, 0.0);
R3Brdf R3red_brdf(1.0, 0.0, 0.0);
R3Brdf R3green_brdf(0.0, 1.0, 0.0);
R3Brdf R3blue_brdf(0.0, 0.0, 1.0);
R3Brdf R3yellow_brdf(1.0, 1.0, 0.0);
R3Brdf R3cyan_brdf(0.0, 1.0, 1.0);
R3Brdf R3magenta_brdf(1.0, 0.0, 1.0);
R3Brdf R3white_brdf(1.0, 1.0, 1.0);
R3Brdf R3gray_brdf(0.5, 0.5, 0.5);
R3Brdf R3shiny_black_brdf(0.0, 0.0, 0.0, 0.8);
R3Brdf R3shiny_red_brdf(1.0, 0.0, 0.0, 0.8);
R3Brdf R3shiny_green_brdf(0.0, 1.0, 0.0, 0.8);
R3Brdf R3shiny_blue_brdf(0.0, 0.0, 1.0, 0.8);
R3Brdf R3shiny_yellow_brdf(1.0, 1.0, 0.0, 0.8);
R3Brdf R3shiny_cyan_brdf(0.0, 1.0, 1.0, 0.8);
R3Brdf R3shiny_magenta_brdf(1.0, 0.0, 1.0, 0.8);
R3Brdf R3shiny_white_brdf(1.0, 1.0, 1.0, 0.8);
R3Brdf R3shiny_gray_brdf(0.5, 0.5, 0.5, 0.8);



/* Public functions */

int 
R3InitBrdf()
{
    /* Return success */
    return TRUE;
}



void 
R3StopBrdf()
{
}



R3Brdf::
R3Brdf(const char *name)
  : scene(NULL),
    scene_index(-1),
    name(NULL),
    ambient(0.0, 0.0, 0.0),
    diffuse(0.0, 0.0, 0.0),
    specular(0.0, 0.0, 0.0),
    transmission(0.0, 0.0, 0.0),
    emission(0.0, 0.0, 0.0),
    shininess(0),
    indexofrefraction(0),
    flags(RN_NO_FLAGS), 
    id(0)
{
    // Set name
    if (name) this->name = strdup(name);
}



R3Brdf::
R3Brdf(const R3Brdf& brdf, const char *name)
    : scene(NULL),
      scene_index(-1),
      name(NULL),
      ambient(brdf.ambient),
      diffuse(brdf.diffuse),
      specular(brdf.specular),
      transmission(brdf.transmission),
      emission(brdf.emission),
      shininess(brdf.shininess),
      indexofrefraction(brdf.indexofrefraction),
      flags(brdf.flags),
      id(-1)
{
    // Set name
    if (name) this->name = strdup(name);
}



R3Brdf::
R3Brdf(const RNRgb& rgb, 
       RNScalar shininess, 
       RNScalar opacity,
       RNScalar indexofrefraction,
       const char *name)
    : scene(NULL),
      scene_index(-1),
      name(NULL),
      ambient(rgb),
      diffuse(rgb),
      specular(0.0, 0.0, 0.0),
      transmission(1 - opacity, 1 - opacity, 1 - opacity),
      emission(0.0, 0.0, 0.0),
      shininess(shininess),
      indexofrefraction(indexofrefraction),
      flags(RN_NO_FLAGS),
      id(-1)
{
    // Set name
    if (name) this->name = strdup(name);

    // Update brdf flags
    Update();
}



R3Brdf::
R3Brdf(RNScalar red, RNScalar green, RNScalar blue,
       RNScalar shininess, 
       RNScalar opacity,
       RNScalar indexofrefraction,
       const char *name)
    : scene(NULL),
      scene_index(-1),
      name(NULL),
      ambient(red, green, blue),
      diffuse(red, green, blue),
      specular(0.0, 0.0, 0.0),
      transmission(1 - opacity, 1 - opacity, 1 - opacity),
      emission(0.0, 0.0, 0.0),
      shininess(shininess),
      indexofrefraction(indexofrefraction),
      flags(RN_NO_FLAGS),
      id(-1)
{
    // Set name
    if (name) this->name = strdup(name);

    // Update brdf
    Update();
}



R3Brdf::
R3Brdf(const RNRgb& ambient, 
       const RNRgb& diffuse, 
       const RNRgb& specular,
       const RNRgb& emission, 
       RNScalar shininess, 
       RNScalar opacity,
       RNScalar indexofrefraction,
       const char *name)
    : scene(NULL),
      scene_index(-1),
      name(NULL),
      ambient(ambient),
      diffuse(diffuse),
      specular(specular),
      transmission(1 - opacity, 1 - opacity, 1 - opacity),
      emission(emission),
      shininess(shininess),
      indexofrefraction(indexofrefraction),
      flags(RN_NO_FLAGS),
      id(-1)
{
    // Set name
    if (name) this->name = strdup(name);

    // Update brdf
    Update();
}



R3Brdf::
R3Brdf(const RNRgb& ambient, 
       const RNRgb& diffuse, 
       const RNRgb& specular,
       const RNRgb& transmission, 
       const RNRgb& emission, 
       RNScalar shininess, 
       RNScalar indexofrefraction,
       const char *name)
    : scene(NULL),
      scene_index(-1),
      name(NULL),
      ambient(ambient),
      diffuse(diffuse),
      specular(specular),
      transmission(transmission),
      emission(emission),
      shininess(shininess),
      indexofrefraction(indexofrefraction),
      flags(RN_NO_FLAGS),
      id(-1)
{
    // Set name
    if (name) this->name = strdup(name);

    // Update brdf
    Update();
}



R3Brdf::
~R3Brdf(void)
{
    // Remove brdf from scene
    if (scene) scene->RemoveBrdf(this);

    // Free name
    if (name) free(name);
}



void R3Brdf::
SetName(const char *name)
{
  // Set name
  if (this->name) free(this->name);
  if (name) this->name = strdup(name);
  else this->name = NULL;
}


  
void R3Brdf::
SetAmbient(const RNRgb& rgb) 
{
    // Set ambient rgb
    if (rgb.IsBlack()) flags.Remove(R3_BRDF_AMBIENT_FLAG);
    else flags.Add(R3_BRDF_AMBIENT_FLAG);
    ambient = rgb;
}



void R3Brdf::
SetDiffuse(const RNRgb& rgb) 
{
    // Set diffuse rgb
    if (rgb.IsBlack()) flags.Remove(R3_BRDF_DIFFUSE_FLAG);
    else flags.Add(R3_BRDF_DIFFUSE_FLAG);
    diffuse = rgb;
}



void R3Brdf::
SetSpecular(const RNRgb& rgb) 
{
    // Set specular rgb
    if (rgb.IsBlack()) flags.Remove(R3_BRDF_SPECULAR_FLAG);
    else flags.Add(R3_BRDF_SPECULAR_FLAG);
    specular = rgb;
}



void R3Brdf::
SetTransmission(const RNRgb& rgb) 
{
    // Set emission rgb
    if (rgb.IsBlack()) flags.Remove(R3_BRDF_TRANSPARENCY_FLAG);
    else flags.Add(R3_BRDF_TRANSPARENCY_FLAG);
    transmission = rgb;
}



void R3Brdf::
SetEmission(const RNRgb& rgb) 
{
    // Set emission rgb
    if (rgb.IsBlack()) flags.Remove(R3_BRDF_EMISSION_FLAG);
    else flags.Add(R3_BRDF_EMISSION_FLAG);
    emission = rgb;
}



void R3Brdf::
SetShininess(RNScalar shininess) 
{
    // Set shininess
    if (shininess == 0.0) flags.Remove(R3_BRDF_SHININESS_FLAG);
    else flags.Add(R3_BRDF_SHININESS_FLAG);
    this->shininess = shininess;
}



void R3Brdf::
SetOpacity(RNScalar opacity)
{
    // Set opacity
    RNScalar t = 1 - opacity;
    if (RNIsEqual(t, 0.0)) t = 0.0;
    else if (RNIsEqual(t, 1.0)) t = 1.0;
    if (t == 0.0) flags.Remove(R3_BRDF_TRANSPARENCY_FLAG);    
    else flags.Add(R3_BRDF_TRANSPARENCY_FLAG);
    this->transmission = RNRgb(t, t, t);
}



void R3Brdf::
SetIndexOfRefraction(RNScalar indexofrefraction)
{
    // Set indexofrefraction
    this->indexofrefraction = indexofrefraction;
}



void R3Brdf::
Update (void)
{
    // Update flags
    UpdateFlags(RN_NO_FLAGS);
}



void R3Brdf::
UpdateFlags (const RNFlags flags)
{
    // Update flags
    this->flags = flags;
    if (!ambient.IsBlack()) this->flags.Add(R3_BRDF_AMBIENT_FLAG);
    if (!diffuse.IsBlack()) this->flags.Add(R3_BRDF_DIFFUSE_FLAG);
    if (!specular.IsBlack()) this->flags.Add(R3_BRDF_SPECULAR_FLAG);
    if (!transmission.IsBlack()) this->flags.Add(R3_BRDF_TRANSPARENCY_FLAG);
    if (!emission.IsBlack()) this->flags.Add(R3_BRDF_EMISSION_FLAG);
    if (shininess != 0.0) this->flags.Add(R3_BRDF_SHININESS_FLAG);
}



void R3Brdf::
Draw(RNBoolean force) const
{
    // Check if same brdf
    static const R3Brdf *R3current_brdf = NULL;
    if (!force && (this == R3current_brdf)) return;

    // Set brdf attributes
    GLfloat c[4];
    c[0] = Ambient().R();
    c[1] = Ambient().G();
    c[2] = Ambient().B();
    c[3] = Opacity();
    glMaterialfv(GL_FRONT, GL_AMBIENT, c);
    c[0] = Diffuse().R();
    c[1] = Diffuse().G();
    c[2] = Diffuse().B();
    c[3] = Opacity();
    glMaterialfv(GL_FRONT, GL_DIFFUSE, c);
    c[0] = Specular().R();
    c[1] = Specular().G();
    c[2] = Specular().B();
    c[3] = Opacity();
    glMaterialfv(GL_FRONT, GL_SPECULAR, c);
    c[0] = Emission().R();
    c[1] = Emission().G();
    c[2] = Emission().B();
    c[3] = Opacity();
    glMaterialfv(GL_FRONT, GL_EMISSION, c);
    glMaterialf(GL_FRONT, GL_SHININESS, Shininess());
    if (IsTransparent()) {
      glDepthMask(FALSE);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glEnable(GL_BLEND);
    }
    else {
      glDisable(GL_BLEND);
      glBlendFunc(GL_ONE, GL_ZERO);
      glDepthMask(TRUE);
    }

    // Remember brdf
    R3current_brdf = this;
}















