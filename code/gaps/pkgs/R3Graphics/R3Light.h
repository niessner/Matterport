/* Include file for the R3 light class */



/* Initialization functions */

int R3InitLight();
void R3StopLight();



/* Class definition */

class R3Light {
    public:
        // Constructor functions
	R3Light(const char *name = NULL);
        R3Light(const R3Light& light, const char *name = NULL);
        R3Light(const RNRgb& color, RNScalar intensity = 1.0, RNBoolean active = TRUE, const char *name = NULL);
        virtual ~R3Light(void);

	// Property functions/operations
        R3Scene *Scene(void) const;
        int SceneIndex(void) const;
        const char *Name(void) const;
	const RNBoolean IsActive(void) const;
  	const RNScalar Intensity(void) const;
  	const RNRgb& Color(void) const;
        R3Light *Copy(void) const;

        // Manipulation functions/operations
        virtual void SetName(const char *name);
	virtual void SetActive(RNBoolean active);
  	virtual void SetIntensity(RNScalar intensity);
  	virtual void SetColor(const RNRgb& color);
        virtual void Transform(const R3Transformation& transformation);

	// Geometry evaluation functions
	virtual RNRgb IrradianceAtPoint(const R3Point& point) const;
	virtual RNScalar IntensityAtPoint(const R3Point& point) const = 0;
	virtual RNScalar RadiusOfInfluence(RNScalar intensity) const = 0;
	virtual R3Sphere SphereOfInfluence(RNScalar intensity) const = 0;

        // Reflection evaluation functions
	virtual RNRgb Reflection(const R3Brdf& brdf, const R3Point& eye, 
	    const R3Point& point, const R3Vector& normal) const = 0;
	virtual RNRgb DiffuseReflection(const R3Brdf& brdf, 
	    const R3Point& point, const R3Vector& normal) const = 0;
	virtual RNRgb SpecularReflection(const R3Brdf& brdf, const R3Point& eye, 
	    const R3Point& point, const R3Vector& normal) const = 0;

	// Draw functions/operations
        virtual void Draw(int i) const = 0;

	// Class type definitions
	RN_CLASS_TYPE_DECLARATIONS(R3Light);

    private:
        friend class R3Scene;
        R3Scene *scene;
        int scene_index;
        char *name;
	RNBoolean active;
	RNScalar intensity;
	RNRgb color;
        int id;
};



/* Public variables */

extern RNScalar R3ambient_light_intensity;
extern RNRgb R3ambient_light_color;



/* Inline functions */

inline R3Scene *R3Light::
Scene(void) const
{
  // Return scene 
  return scene;
}



inline int R3Light::
SceneIndex(void) const
{
  // Return index of light in scene (can be used with scene->Light(xxx))
  return scene_index;
}



inline const char *R3Light::
Name(void) const
{
    // Return name
    return name;
}



inline const RNBoolean R3Light::
IsActive(void) const
{
    // Return status
    return active;
}



inline const RNScalar R3Light::
Intensity(void) const
{
    // Return intensity 
    return intensity;
}



inline const RNRgb& R3Light::
Color(void) const
{
    // Return color 
    return color;
}



