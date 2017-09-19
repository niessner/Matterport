/* Include file for GAPS basic stuff */



/* Initialization functions */

int R3InitBase();
void R3StopBase();



/* Class definition */

class R3Base {};



/* Draw flags definition */

typedef RNFlags R3DrawFlags;
#define R3_NULL_DRAW_FLAGS                    (0x000)
#define R3_EDGES_DRAW_FLAG                    (0x001)
#define R3_SURFACES_DRAW_FLAG                 (0x002)
#define R3_SURFACE_NORMALS_DRAW_FLAG          (0x010)
#define R3_SURFACE_TEXTURE_DRAW_FLAG          (0x020)
#define R3_SURFACE_MATERIAL_DRAW_FLAG         (0x040)
#define R3_VERTEX_NORMALS_DRAW_FLAG           (0x100)
#define R3_VERTEX_TEXTURE_COORDS_DRAW_FLAG    (0x200)
#define R3_VERTEX_COLORS_DRAW_FLAG            (0x400)
#define R3_EVERYTHING_DRAW_FLAGS              (0x773)
#define R3_DEFAULT_DRAW_FLAGS                 (R3_EVERYTHING_DRAW_FLAGS & ~R3_EDGES_DRAW_FLAG)

#define R3_VERTEX_SHARED_FLAG                (0x1000)
