/* Include file for the R3 curve class */



/* Initialization functions */

int R3InitCurve();
void R3StopCurve();



/* Class definition */

class R3Curve : public R3Shape {
    public:
        // Constructors/destructors ???
	R3Curve(void);
	virtual ~R3Curve(void);

        // Curve properties
        virtual const RNScalar StartParameter(void) const = 0;
        virtual const RNScalar EndParameter(void) const = 0;
        virtual const R3Point StartPosition(void) const;
        virtual const R3Point EndPosition(void) const;

        // Point access
        virtual R3Point PointPosition(RNScalar u) const = 0;
        virtual R3Vector PointDerivative(RNScalar u) const = 0;
        virtual R3Vector PointDirection(RNScalar u) const;

        // Shape property functions/operators
	const RNBoolean IsCurve(void) const;
};



