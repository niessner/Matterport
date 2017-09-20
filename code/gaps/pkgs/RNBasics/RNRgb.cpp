/* Source file for the RN RGB color class */



/* Include files */

#include "RNBasics.h"



/* Public variables */

RNRgb RNnull_rgb(0.0, 0.0, 0.0);
RNRgb RNblack_rgb(0.0, 0.0, 0.0);
RNRgb RNgray_rgb(0.5, 0.5, 0.5);
RNRgb RNred_rgb(1.0, 0.0, 0.0);
RNRgb RNgreen_rgb(0.0, 1.0, 0.0);
RNRgb RNblue_rgb(0.0, 0.0, 1.0);
RNRgb RNyellow_rgb(1.0, 1.0, 0.0);
RNRgb RNcyan_rgb(0.0, 1.0, 1.0);
RNRgb RNmagenta_rgb(1.0, 0.0, 1.0);
RNRgb RNwhite_rgb(1.0, 1.0, 1.0);



/* Public functions */

int 
RNInitRgb()
{
    /* Return success */
    return TRUE;
}



void 
RNStopRgb()
{
}



RNRgb::
RNRgb(void)
{
}



RNRgb::
RNRgb(const RNRgb& rgb)
{
    c[0] = rgb.c[0];
    c[1] = rgb.c[1];
    c[2] = rgb.c[2];
}



RNRgb::
RNRgb(RNScalar red, RNScalar green, RNScalar blue)
{
    c[0] = red; 
    c[1] = green; 
    c[2] = blue;
}



RNRgb::
RNRgb(const RNScalar array[3])
{
    c[0] = array[0]; 
    c[1] = array[1]; 
    c[2] = array[2];
}



const RNBoolean RNRgb::
operator==(const RNRgb& rgb) const
{
    // Return whether rgb is equal
    return ((c[0] == rgb.c[0]) && (c[1] == rgb.c[1]) && (c[2] == rgb.c[2]));
}



const RNBoolean RNRgb::
operator!=(const RNRgb& rgb) const
{
    // Return whether rgb is not equal
    return ((c[0] != rgb.c[0]) || (c[1] != rgb.c[1]) || (c[2] != rgb.c[2]));
}



RNRgb& RNRgb::
operator=(const RNRgb& rgb)
{
    c[0] = rgb.c[0];
    c[1] = rgb.c[1];
    c[2] = rgb.c[2];
    return *this;
}



RNRgb& RNRgb::
operator+=(const RNRgb& rgb)
{
    c[0] += rgb.c[0];
    c[1] += rgb.c[1];
    c[2] += rgb.c[2];
    return *this;
}



RNRgb& RNRgb::
operator-=(const RNRgb& rgb)
{
    c[0] -= rgb.c[0];
    c[1] -= rgb.c[1];
    c[2] -= rgb.c[2];
    return *this;
}



RNRgb& RNRgb::
operator*=(const RNRgb& rgb)
{
    c[0] *= rgb.c[0];
    c[1] *= rgb.c[1];
    c[2] *= rgb.c[2];
    return *this;
}



RNRgb& RNRgb::
operator*=(RNScalar a)
{
    c[0] *= a;
    c[1] *= a;
    c[2] *= a;
    return *this;
}



RNRgb& RNRgb::
operator/=(RNScalar a)
{
    //  assert(!zero(a)); 
    c[0] /= a;
    c[1] /= a;
    c[2] /= a;
    return *this;
}





