/* Source file for the R3 material class */



/* Include files */

#include "R3Graphics.h"



/* Public variables */

R3Material R3null_material;
R3Material R3black_material;
R3Material R3red_material;
R3Material R3green_material;
R3Material R3blue_material;
R3Material R3yellow_material;
R3Material R3cyan_material;
R3Material R3magenta_material;
R3Material R3white_material;
R3Material R3gray_material;
R3Material R3shiny_black_material;
R3Material R3shiny_red_material;
R3Material R3shiny_green_material;
R3Material R3shiny_blue_material;
R3Material R3shiny_yellow_material;
R3Material R3shiny_cyan_material;
R3Material R3shiny_magenta_material;
R3Material R3shiny_white_material;
R3Material R3shiny_gray_material;
R3Material R3default_material;



/* Public functions */

int 
R3InitMaterial()
{
    /* Initialize public variables */
    R3null_material = R3Material(&R3null_brdf, "Null");
    R3black_material = R3Material(&R3black_brdf, "Black");
    R3red_material = R3Material(&R3red_brdf, "Red");
    R3green_material = R3Material(&R3green_brdf, "Green");
    R3blue_material = R3Material(&R3blue_brdf, "Blue");
    R3yellow_material = R3Material(&R3yellow_brdf, "Yellow");
    R3cyan_material = R3Material(&R3cyan_brdf, "Cyan");
    R3magenta_material = R3Material(&R3magenta_brdf, "Magenta");
    R3white_material = R3Material(&R3white_brdf, "White");
    R3gray_material = R3Material(&R3gray_brdf, "Gray");
    R3shiny_black_material = R3Material(&R3shiny_black_brdf, "ShinyBlack");
    R3shiny_red_material = R3Material(&R3shiny_red_brdf, "ShinyRed");
    R3shiny_green_material = R3Material(&R3shiny_green_brdf, "ShinyGreen");
    R3shiny_blue_material = R3Material(&R3shiny_blue_brdf, "ShinyBlue");
    R3shiny_yellow_material = R3Material(&R3shiny_yellow_brdf, "ShinyYellow");
    R3shiny_cyan_material = R3Material(&R3shiny_cyan_brdf, "ShinyCyan");
    R3shiny_magenta_material = R3Material(&R3shiny_magenta_brdf, "ShinyMagenta");
    R3shiny_white_material = R3Material(&R3shiny_white_brdf, "ShinyWhite");
    R3shiny_gray_material = R3Material(&R3shiny_gray_brdf, "ShinyGray");
    R3default_material = R3Material(&R3default_brdf, "Default");
    
    /* Return success */
    return TRUE;
}



void 
R3StopMaterial()
{
}



R3Material::
R3Material(const char *name)
  :   scene(NULL),
      scene_index(-1),
      name(NULL),
      brdf(NULL),
      texture(NULL),
      flags(0)
{
    // Set name
    if (name) this->name = strdup(name);

    // Update flags
    Update();
}



R3Material::
R3Material(const R3Material& material)
    : scene(NULL),
      scene_index(-1),
      name(NULL),
      brdf(material.brdf),
      texture(material.texture),
      flags(material.flags)
{
    // Set name
    if (material.name) this->name = strdup(material.name);
}



R3Material::
R3Material(const R3Brdf *brdf, const char *name)
    : scene(NULL),
      scene_index(-1),
      name(NULL),
      brdf(brdf),
      texture(&R2null_texture),
      flags(0)
{
    // Set name
    if (name) this->name = strdup(name);

    // Update flags
    Update();
}



R3Material::
R3Material(const R2Texture *texture, const char *name)
    : scene(NULL),
      scene_index(-1),
      name(NULL),
      brdf(&R3white_brdf),
      texture(texture),
      flags(0)
{
    // Set name
    if (name) this->name = strdup(name);

    // Update flags
    Update();
}



R3Material::
R3Material(const R3Brdf *brdf, const R2Texture *texture, const char *name)
    : scene(NULL),
      scene_index(-1),
      name(NULL),
      brdf(brdf),
      texture(texture),
      flags(0)
{
    // Set name
    if (name) this->name = strdup(name);

    // Update flags
    Update();
}



R3Material::
~R3Material(void)
{
    // Remove from scene
    if (scene) scene->RemoveMaterial(this);
  
    // Free name
    if (name) free(name);
}


 
R3Material& R3Material::
operator=(const R3Material& material)
{
    // Update everything
    SetBrdf(material.brdf);
    SetTexture(material.texture);
    SetName(material.name);
    Update();
    return *this;
}



void R3Material::
SetName(const char *name)
{
    // Set name
    if (this->name) free(this->name);
    if (name) this->name = strdup(name);
    else this->name = NULL;
}



void R3Material::
SetBrdf(const R3Brdf *brdf)
{
    // Set brdf
    if ((this->brdf) && (!brdf)) flags.Remove(R3_MATERIAL_BRDF_FLAG);
    if ((!this->brdf) && (brdf)) flags.Add(R3_MATERIAL_BRDF_FLAG);
    if (this->brdf) flags.Remove(this->brdf->Flags());
    if (brdf) flags.Add(brdf->Flags());
    this->brdf = brdf;
}




void R3Material::
SetTexture(const R2Texture *texture)
{
    // Set texture
    if ((this->texture) && (!texture)) flags.Remove(R3_MATERIAL_TEXTURE_FLAG);
    if ((!this->texture) && (texture)) flags.Add(R3_MATERIAL_TEXTURE_FLAG);
    if (this->texture) flags.Remove(this->texture->Flags());
    if (texture) flags.Add(texture->Flags());
    this->texture = texture;
}



void R3Material::
Update(void)
{
    // Update flags
    this->flags = RN_NO_FLAGS;
    if ((brdf) && (brdf->ID() != 0)) 
	this->flags.Add(R3_MATERIAL_BRDF_FLAG | brdf->Flags());
    if ((texture) && (texture->ID() != 0)) 
	this->flags.Add(R3_MATERIAL_TEXTURE_FLAG | texture->Flags());
}




void R3Material::
Load(void) const
{
    // Load brdf
    if (brdf) brdf->Load();

    // Load texture
    if (texture) texture->Load();
}



void R3Material::
Unload(void) const
{
    // Unload brdf
    if (brdf) brdf->Unload();

    // Unload texture
    if (texture) texture->Unload();
}



void R3Material::
Draw(RNBoolean force) const
{
    // Check if same material as last time
    static const R3Material *R3current_material = NULL;
    if (!force && (this == R3current_material)) return;

    // Draw brdf
    if (brdf) brdf->Draw(force);
    else R3default_brdf.Draw(force);

    // Draw texture
    if (texture) texture->Draw(force);
    else R2null_texture.Draw(force);

    // Remember current material
    R3current_material = this;
}























