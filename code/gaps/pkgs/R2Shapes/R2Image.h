// Include file for image class



// Class definition

class R2Image {
 public:
  // Constructors
  R2Image(void);
  R2Image(const char *filename);
  R2Image(int width, int height, int ncomponents = 3);
  R2Image(int width, int height, int ncomponents, unsigned char *data);
  R2Image(const R2Grid& red, const R2Grid& green, const R2Grid& blue);
  R2Image(const R2Image& image);
  ~R2Image(void);

  // Accessors
  const unsigned char *Pixels(void) const;
  const unsigned char *Pixels(int row) const;
  const unsigned char *Pixel(int row, int column) const;
  const RNRgb PixelRGB(int row, int column) const;
  int Width(void) const;
  int Height(void) const;
  int Depth(void) const;
  int NComponents(void) const;
  int RowSize(void) const;
  int Size(void) const;

  // Manipulation
  R2Image& operator=(const R2Image& image);
  void Clear(const RNRgb& rgb = RNblack_rgb);
  void Add(const R2Image& image);
  void Subtract(const R2Image& image);
  void SetPixel(int row, int column, const unsigned char *pixel);
  void SetPixelRGB(int row, int column, const RNRgb& rgb);

  // Reading/writing
  int Read(const char *filename);
  int ReadBMP(const char *filename);
  int ReadPPM(const char *filename);
  int ReadPFM(const char *filename);
  int ReadJPEG(const char *filename);
  int ReadTIFF(const char *filename);
  int ReadPNG(const char *filename);
  int ReadRAW(const char *filename);
  int ReadGRD(const char *filename);
  int Write(const char *filename) const;
  int WriteBMP(const char *filename) const;
  int WritePPM(const char *filename, int ascii = 0) const;
  int WriteRAW(const char *filename) const;
  int WriteJPEG(const char *filename) const;
  int WriteTIFF(const char *filename) const;
  int WritePNG(const char *filename) const;
  void Capture(void);

  // Draw functions
  void Draw(int x = 0, int y = 0) const;

 private:
  int width;
  int height;
  int ncomponents;
  int rowsize;
  unsigned char *pixels;
};



// Inline functions

inline int R2Image::
Width(void) const
{
  // Return width
  return width;
}



inline int R2Image::
Height(void) const
{
  // Return height
  return height;
}



inline int R2Image::
Depth(void) const
{
  // Return number of bytes per pixel
  return ncomponents;
}



inline int R2Image::
NComponents(void) const
{
  // Return number of bytes per pixel
  return Depth();
}



inline int R2Image::
RowSize(void) const
{
  // Return size of row in bytes
  return rowsize;
}



inline int R2Image::
Size(void) const
{
  // Return size of image in bytes
  return rowsize * height;
}



inline const unsigned char *R2Image::
Pixels(void) const
{
  // Return pixels pointer (pixels start at lower-left)
  return pixels;
}



inline const unsigned char *R2Image::
Pixels(int y) const
{
  // Return pixels pointer for row y
  return &pixels[y*rowsize];
}



inline const unsigned char *R2Image::
Pixel(int x, int y) const
{
  // Return pixel value at (x,y)
  return &pixels[y*rowsize + x*ncomponents];
}












