/* Include file for the R3 light class */



/* Initialization functions */

int R3InitSpotLight();
void R3StopSpotLight();



/* Class definition */

class R3SpotLight : public R3PointLight {
    public:
        // Constructor functions
	R3SpotLight(void);
        R3SpotLight(const R3SpotLight& light);
        R3SpotLight(const R3Point& position, const R3Vector& direction,
	    const RNRgb& color, RNScalar dropoffrate = 0.0, RNAngle cutoffangle = 0.785398,
            RNScalar intensity = 1.0, RNBoolean active = TRUE,
            RNScalar ca = 0, RNScalar la = 0, RNScalar qa = 1);

	// Property functions/operators
  	const R3Vector& Direction(void) const;
  	const RNScalar DropOffRate(void) const;
  	const RNAngle CutOffAngle(void) const;

	// Manipulation functions/operations
  	virtual void SetDirection(const R3Vector& direction);
  	virtual void SetDropOffRate(RNScalar dropoffrate);
  	virtual void SetCutOffAngle(RNAngle cutoffangle);
        virtual void Transform(const R3Transformation& transformation);

	// Evaluation functions
	virtual RNScalar IntensityAtPoint(const R3Point& point) const;

	// Draw functions/operations
        virtual void Draw(int i) const;

	// Class type definitions
	RN_CLASS_TYPE_DECLARATIONS(R3SpotLight);

    private:
	R3Vector direction;
	RNScalar dropoffrate;
	RNAngle cutoffangle;
};



/* Public variables */

extern R3SpotLight R3null_spot_light;



/* Inline functions */

inline const R3Vector& R3SpotLight::
Direction(void) const
{
    // Return direction 
    return direction;
}



inline const RNScalar R3SpotLight::
DropOffRate(void) const
{
    // Return drop off rate 
    return dropoffrate;
}



inline const RNAngle R3SpotLight::
CutOffAngle(void) const
{
    // Return cut off angle 
    return cutoffangle;
}





