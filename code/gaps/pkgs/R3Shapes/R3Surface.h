/* Include file for the R3 surface class */



/* Initialization functions */

int R3InitSurface();
void R3StopSurface();



/* Class definition */

class R3Surface : public R3Shape {
    public:
        // Constructors/destructors ???
	R3Surface(void);
	~R3Surface(void);

        // Shape property functions/operators
	const RNBoolean IsSurface(void) const;
};



