/* Source file for the R2 shape class */



/* Include files */

#include "R2Shapes/R2Shapes.h"



/* Public draw flags */

const R2DrawFlags R2_EDGES_DRAW_FLAG                    (0x001);
const R2DrawFlags R2_SURFACES_DRAW_FLAG                 (0x002);
const R2DrawFlags R2_NULL_DRAW_FLAGS                    (0x000);
const R2DrawFlags R2_EVERYTHING_DRAW_FLAGS              (0xFFF);
const R2DrawFlags R2_DEFAULT_DRAW_FLAGS                 (R2_EVERYTHING_DRAW_FLAGS);



/* Public functions */

int 
R2InitShape()
{
    /* Return success */
    return TRUE;
}



void 
R2StopShape()
{
}




