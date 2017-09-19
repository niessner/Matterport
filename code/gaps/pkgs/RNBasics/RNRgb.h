/* Include file for the RN rgb class */



/* Initialization functions */

int RNInitRgb();
void RNStopRgb();



/* Class definition */

class RNRgb /* : public RNColor */ {
    public:
        // Constructor functions
	RNRgb(void);
	RNRgb(const RNRgb& rgb);
        RNRgb(RNScalar red, RNScalar green, RNScalar blue);
	RNRgb(const RNScalar array[3]);

        // Property functions/operators
	const RNScalar R(void) const;
	const RNScalar G(void) const;
	const RNScalar B(void) const;
	const RNScalar Coord(int i) const;
	const RNScalar operator[](int i) const;
	const RNScalar *Coords(void) const;
        const RNScalar Luminance(void) const;
	const int IsBlack(void) const;
	const int IsWhite(void) const;
	const int operator==(const RNRgb& rgb) const;
	const int operator!=(const RNRgb& rgb) const;

	// Manipulation functions/operations
	void SetRed(RNScalar red);
	void SetGreen(RNScalar green);
	void SetBlue(RNScalar blue);
	void Reset(RNScalar red, RNScalar green, RNScalar blue);

	// Assignment operators
	RNRgb& operator=(const RNRgb& rgb);
	RNRgb& operator+=(const RNRgb& rgb);
	RNRgb& operator-=(const RNRgb& rgb);
	RNRgb& operator*=(const RNRgb& rgb);
	RNRgb& operator*=(RNScalar a);
	RNRgb& operator/=(RNScalar a);
 
	// Arithmetic operations
	friend RNRgb operator+(const RNRgb& rgb1, const RNRgb& rgb2);
	friend RNRgb operator-(const RNRgb& rgb1, const RNRgb& rgb2);
	friend RNRgb operator*(const RNRgb& rgb1, const RNRgb& rgb2);
	friend RNRgb operator*(const RNRgb& rgb, RNScalar a);
	friend RNRgb operator*(RNScalar a, const RNRgb& rgb);
	friend RNRgb operator/(const RNRgb& rgb, RNScalar a);

        // Undocumented functions/operators
  	RNScalar& operator[](int i);

    private:
	RNScalar c[3];
};



/* Public variables */

extern RNRgb RNnull_rgb;
extern RNRgb RNblack_rgb;
extern RNRgb RNgray_rgb;
extern RNRgb RNred_rgb;
extern RNRgb RNgreen_rgb;
extern RNRgb RNblue_rgb;
extern RNRgb RNyellow_rgb;
extern RNRgb RNcyan_rgb;
extern RNRgb RNmagenta_rgb;
extern RNRgb RNwhite_rgb;



/* Inline functions */

inline const RNScalar RNRgb::
R (void) const
{
    return(c[0]);
}



inline const RNScalar RNRgb::
G (void) const
{
    return(c[1]);
}



inline const RNScalar RNRgb::
B (void) const
{
    return(c[2]);
}



inline const RNScalar RNRgb::
Coord(int i) const
{
    assert ((i>=0)&&(i<=2));
    return(c[i]);
}



inline const RNScalar RNRgb::
operator[](int i) const
{
    return Coord(i);
}



inline const RNScalar *RNRgb::
Coords(void) const
{
    // Return rgb array
    return c;
}



inline const int RNRgb::
IsBlack (void) const
{
    // Return whether color is black
    return ((c[0] == 0.0) && (c[1] == 0.0) && (c[2] == 0.0));
}



inline const int RNRgb::
IsWhite (void) const
{
    // Return whether color is white
    return ((c[0] == 1.0) && (c[1] == 1.0) && (c[2] == 1.0));
}



inline const RNScalar RNRgb::
Luminance(void) const
{
    // Return luminance
    return 0.30 * c[0] + 0.59 * c[1] + 0.11 * c[2];
}



inline void RNRgb::
SetRed(RNScalar red)
{
    // Set red component
    c[0] = red;
}



inline void RNRgb::
SetGreen(RNScalar green)
{
    // Set green component
    c[1] = green;
}



inline void RNRgb::
SetBlue(RNScalar blue)
{
    // Set blue component
    c[2] = blue;
}



inline void RNRgb::
Reset (RNScalar red, RNScalar green, RNScalar blue)
{
    // Set all components
    c[0] = red;
    c[1] = green;
    c[2] = blue;
}



inline RNRgb 
operator+(const RNRgb& rgb1, const RNRgb& rgb2)
{
    return RNRgb(rgb1.c[0] + rgb2.c[0], 
		 rgb1.c[1] + rgb2.c[1], 
		 rgb1.c[2] + rgb2.c[2]);
}



inline RNRgb 
operator-(const RNRgb& rgb1, const RNRgb& rgb2)
{
    return RNRgb(rgb1.c[0] - rgb2.c[0], 
		 rgb1.c[1] - rgb2.c[1], 
		 rgb1.c[2] - rgb2.c[2]);
}



inline RNRgb 
operator*(const RNRgb& rgb1, const RNRgb& rgb2)
{
    return RNRgb(rgb1.c[0] * rgb2.c[0], 
		 rgb1.c[1] * rgb2.c[1], 
		 rgb1.c[2] * rgb2.c[2]);
}



inline RNRgb 
operator*(const RNRgb& rgb, RNScalar a)
{
    return RNRgb(rgb.c[0] * a, 
		 rgb.c[1] * a, 
		 rgb.c[2] * a);
}



inline RNRgb 
operator*(RNScalar a, const RNRgb& rgb)
{
    return rgb * a;
}



inline RNRgb 
operator/(const RNRgb& rgb, RNScalar a)
{
    return RNRgb(rgb.c[0] / a, 
		 rgb.c[1] / a, 
		 rgb.c[2] / a);
}



inline RNScalar& RNRgb::
operator[] (int i) 
{
    assert ((i>=0)&&(i<=2));
    return(c[i]);
}



inline void 
RNLoadRgb(const RNRgb& rgb)
{
    // Load rgb
    RNLoadRgb(rgb.R(), rgb.G(), rgb.B());
}




