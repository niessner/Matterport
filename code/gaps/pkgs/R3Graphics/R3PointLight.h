/* Include file for the R3 point light class */



/* Initialization functions */

int R3InitPointLight();
void R3StopPointLight();



/* Class definition */

class R3PointLight : public R3Light {
    public:
        // Constructor functions
	R3PointLight(void);
        R3PointLight(const R3PointLight& light);
        R3PointLight(const R3Point& position, const RNRgb& color, 
            RNScalar intensity = 1.0, RNBoolean active = TRUE,
            RNScalar ca = 0, RNScalar la = 0, RNScalar qa = 1);

	// Property functions/operators
  	const R3Point& Position(void) const;
        const RNScalar ConstantAttenuation(void) const;
        const RNScalar LinearAttenuation(void) const;
        const RNScalar QuadraticAttenuation(void) const;

	// Manipulation functions/operations
  	virtual void SetPosition(const R3Point& position);
        virtual void SetConstantAttenuation(RNScalar ca); 
        virtual void SetLinearAttenuation(RNScalar la); 
        virtual void SetQuadraticAttenuation(RNScalar qa);
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
	RN_CLASS_TYPE_DECLARATIONS(R3PointLight);

    private:
	R3Point position;
        RNScalar constant_attenuation;
        RNScalar linear_attenuation;
        RNScalar quadratic_attenuation;
};



/* Public variables */

extern R3PointLight R3null_point_light;



/* Inline functions */

inline const R3Point& R3PointLight::
Position(void) const
{
    // Return position 
    return position;
}



inline const RNScalar R3PointLight::
ConstantAttenuation(void) const
{
    // Return constant coefficient of attenuation
    return constant_attenuation;
}



inline const RNScalar R3PointLight::
LinearAttenuation(void) const
{
    // Return linear coefficient of attenuation
    return linear_attenuation;
}



inline const RNScalar R3PointLight::
QuadraticAttenuation(void) const
{
    // Return quadratic coefficient of attenuation
    return quadratic_attenuation;
}




