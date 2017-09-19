/* Include file for the R3 light class */



/* Initialization functions */

int R3InitDirectionalLight();
void R3StopDirectionalLight();



/* Class definition */

class R3DirectionalLight : public R3Light {
    public:
        // Constructor functions
	R3DirectionalLight(void);
        R3DirectionalLight(const R3DirectionalLight& light);
        R3DirectionalLight(const R3Vector& direction, const RNRgb& color, 
            RNScalar intensity = 1.0, RNBoolean active = TRUE);

	// Property functions/operators
  	const R3Vector& Direction(void) const;

	// Manipulation functions/operations
  	virtual void SetDirection(const R3Vector& direction);
        virtual void Transform(const R3Transformation& transformation);

	// Geometry evaluation functions
	virtual RNScalar IntensityAtPoint(const R3Point& point) const;
	virtual R3Vector DirectionFromPoint(const R3Point& point) const;
	virtual RNScalar RadiusOfInfluence(RNScalar intensity) const;
	virtual R3Sphere SphereOfInfluence(RNScalar intensity) const;

	// Reflection evaluation functions
	virtual RNRgb Reflection(const R3Brdf& brdf, const R3Point& eye, 
	    const R3Point& point, const R3Vector& normal) const;
	virtual RNRgb DiffuseReflection(const R3Brdf& brdf, 
	    const R3Point& point, const R3Vector& normal) const;
	virtual RNRgb SpecularReflection(const R3Brdf& brdf, const R3Point& eye, 
	    const R3Point& point, const R3Vector& normal) const;

	// Draw functions/operations
        virtual void Draw(int i) const;

	// Class type definitions
	RN_CLASS_TYPE_DECLARATIONS(R3DirectionalLight);

    private:
	R3Vector direction;
};



/* Public variables */

extern R3DirectionalLight R3null_directional_light;
extern R3DirectionalLight R3default_directional_light;



/* Inline functions */

inline const R3Vector& R3DirectionalLight::
Direction(void) const
{
    // Return direction 
    return direction;
}




