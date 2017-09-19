/* Include file for the R3 point light class */



/* Initialization functions */

int R3InitAreaLight();
void R3StopAreaLight();



/* Class definition */

class R3AreaLight : public R3Light {
    public:
        // Constructor functions
	R3AreaLight(void);
        R3AreaLight(const R3AreaLight& light);
        R3AreaLight(const R3Point& position, RNLength radius, const R3Vector& direction, const RNRgb& color, 
            RNScalar intensity = 1.0, RNBoolean active = TRUE,
            RNScalar ca = 0, RNScalar la = 0, RNScalar qa = 1);

	// Property functions/operators
  	const R3Point& Position(void) const;
  	const R3Vector& Direction(void) const;
        const RNLength Radius(void) const;
        const RNScalar ConstantAttenuation(void) const;
        const RNScalar LinearAttenuation(void) const;
        const RNScalar QuadraticAttenuation(void) const;

	// Manipulation functions/operations
  	virtual void SetPosition(const R3Point& position);
  	virtual void SetDirection(const R3Vector& direction);
  	virtual void SetRadius(RNLength radius);
        virtual void SetConstantAttenuation(RNScalar ca); 
        virtual void SetLinearAttenuation(RNScalar la); 
        virtual void SetQuadraticAttenuation(RNScalar qa); 
        virtual void Transform(const R3Transformation& transformation);

	// Geometry evaluation functions
	virtual RNScalar IntensityAtPoint(const R3Point& point) const;
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
	RN_CLASS_TYPE_DECLARATIONS(R3AreaLight);

    private:
	R3Circle circle;
        RNScalar constant_attenuation;
        RNScalar linear_attenuation;
        RNScalar quadratic_attenuation;
};



/* Public variables */

extern R3AreaLight R3null_area_light;



/* Inline functions */

inline const R3Point& R3AreaLight::
Position(void) const
{
    // Return position 
    return circle.Center();
}



inline const R3Vector& R3AreaLight::
Direction(void) const
{
    // Return direction 
    return circle.Normal();
}



inline const RNLength R3AreaLight::
Radius(void) const
{
    // Return radius 
    return circle.Radius();
}



inline const RNScalar R3AreaLight::
ConstantAttenuation(void) const
{
    // Return constant coefficient of attenuation
    return constant_attenuation;
}



inline const RNScalar R3AreaLight::
LinearAttenuation(void) const
{
    // Return linear coefficient of attenuation
    return linear_attenuation;
}



inline const RNScalar R3AreaLight::
QuadraticAttenuation(void) const
{
    // Return quadratic coefficient of attenuation
    return quadratic_attenuation;
}




