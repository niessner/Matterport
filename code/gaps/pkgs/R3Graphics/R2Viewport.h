/* Include file for R3 viewport class */



/* Initialization functions */

int R3InitViewport();
void R3StopViewport();



/* Class definition */

class R2Viewport {
    public:
        // Constructor functions
	R2Viewport(void);
        R2Viewport(int xmin, int ymin, int width, int height);

        // Property functions/operators
  	const int XMin(void) const;
  	const int YMin(void) const;
  	const int XMax(void) const;
  	const int YMax(void) const;

        // Property functions/operators
  	const int Width(void) const;
  	const int Height(void) const;
  	const int XCenter(void) const;
  	const int YCenter(void) const;
	const int Area(void) const;
        const R2Box BBox(void) const;
	const RNBoolean operator==(const R2Viewport& viewport) const;
	const RNBoolean operator!=(const R2Viewport& viewport) const;

	// Manipulation functions/operators
	void Move(int xmin, int ymin);
	void Resize(int width, int height);
        void Resize(int xmin, int ymin, int width, int height);

        // Draw functions 
        void Load(void) const;
	
    private:
	int xmin;
	int ymin;
	int width;
	int height;
};



/* Public variables */

extern R2Viewport R2default_viewport;



/* Inline functions */

inline const int R2Viewport::
XMin (void) const
{
    return xmin;
}



inline const int R2Viewport::
YMin (void) const
{
    return ymin;
}



inline const int R2Viewport::
XMax (void) const
{
    return xmin + width - 1;
}



inline const int R2Viewport::
YMax (void) const
{
    return ymin + height - 1;
}



inline const int R2Viewport::
Width (void) const
{
    return width;
}



inline const int R2Viewport::
Height (void) const
{
    return height;
}



inline const int R2Viewport::
XCenter (void) const
{
    return xmin + (width / 2);
}



inline const int R2Viewport::
YCenter (void) const
{
    return ymin + (height / 2);
}



inline const int R2Viewport::
Area (void) const
{
    return width * height;
}



inline const R2Box R2Viewport::
BBox (void) const
{
    return R2Box(xmin, ymin, xmin + width, ymin + height);
}



inline const RNBoolean R2Viewport::
operator==(const R2Viewport& viewport) const
{
    // Return whether viewport is not equal
    return (xmin == viewport.xmin) && 
           (ymin == viewport.ymin) && 
           (width == viewport.width) && 
           (height == viewport.height);
}



inline const RNBoolean R2Viewport::
operator!=(const R2Viewport& viewport) const
{
    // Return whether viewport is not equal
    return !(*this == viewport);
}



inline void R2Viewport::
Load(void) const
{
    glViewport(xmin, ymin, xmin + width, ymin + height);
}











