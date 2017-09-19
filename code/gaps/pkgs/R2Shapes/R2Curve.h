/* Include file for the R2 curve class */



/* Initialization functions */

int R2InitCurve();
void R2StopCurve();



/* Class definition */

class R2Curve : public R2Shape {
    public:
        // Constructors/destructors ???
	R2Curve(void);
	~R2Curve(void);

        // Shape property functions/operators
	const RNBoolean IsCurve(void) const;

	// Curve property functions/operators
	// virtual const R2Point Start(void) const = 0;
	// virtual const R2Point End(void) const = 0;
};



