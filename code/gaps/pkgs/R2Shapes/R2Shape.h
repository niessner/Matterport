/* Include file for the R2 shape class */



/* Initialization functions */

int R2InitShape();
void R2StopShape();



/* Draw flags definition */

typedef RNFlags R2DrawFlags;
extern const R2DrawFlags R2_EDGES_DRAW_FLAG;
extern const R2DrawFlags R2_SURFACES_DRAW_FLAG;
extern const R2DrawFlags R2_NULL_DRAW_FLAGS;
extern const R2DrawFlags R2_EVERYTHING_DRAW_FLAGS;
extern const R2DrawFlags R2_DEFAULT_DRAW_FLAGS;



/* Class definition */

class R2Shape /* : public R2Base */ {
    public:
        // Constructors/destructors
        virtual ~R2Shape(void) {};

        // Draw functions/operators
        virtual void Draw(const R2DrawFlags draw_flags = R2_DEFAULT_DRAW_FLAGS) const = 0;
        void Outline(void) const { Draw(R2_EDGES_DRAW_FLAG); };
};



/* Shape type enumeration ??? */

enum {
    R2_POINT_CLASS_ID = 1,
    R2_LINE_CLASS_ID = 2,
    R2_RAY_CLASS_ID = 3,
    R2_SPAN_CLASS_ID = 4,
    R2_POLYGON_CLASS_ID = 6
};






