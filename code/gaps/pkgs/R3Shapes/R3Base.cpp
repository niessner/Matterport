/* Source file for GAPS basics  */



/* Include files */

#include "R3Shapes/R3Shapes.h"




/* Public draw flags */

#if FALSE

const R3DrawFlags R3_EDGES_DRAW_FLAG                    (0x001);
const R3DrawFlags R3_SURFACES_DRAW_FLAG                 (0x002);
const R3DrawFlags R3_SURFACE_NORMALS_DRAW_FLAG          (0x010);
const R3DrawFlags R3_SURFACE_TEXTURE_DRAW_FLAG          (0x020);
const R3DrawFlags R3_VERTEX_NORMALS_DRAW_FLAG           (0x100);
const R3DrawFlags R3_VERTEX_TEXTURE_COORDS_DRAW_FLAG    (0x200);

const R3DrawFlags R3_NULL_DRAW_FLAGS                    (0x000);
const R3DrawFlags R3_EVERYTHING_DRAW_FLAGS              (0xFFF);
const R3DrawFlags R3_DEFAULT_DRAW_FLAGS                 (R3_EVERYTHING_DRAW_FLAGS & 
							 ~R3_VERTEX_TEXTURE_COORDS_DRAW_FLAG &
							 ~R3_EDGES_DRAW_FLAG);
#endif



int R3InitBase() 
{
    /* Return OK status */
    return TRUE;
}



void R3StopBase()
{
}




