// Source file for image class



// Include files

#include "R2Shapes.h"
#include "ctype.h"



R2Image::
R2Image(void)
  : width(0), 
    height(0),
    ncomponents(0),
    rowsize(0),
    pixels(NULL)
{
}



R2Image::
R2Image(const char *filename)
  : width(0), 
    height(0),
    ncomponents(0),
    rowsize(0),
    pixels(NULL)
{
  // Read image
  Read(filename);
}



R2Image::
R2Image(int width, int height, int ncomponents)
  : width(width), 
    height(height),
    ncomponents(ncomponents),
    rowsize(0),
    pixels(NULL)
{
  // Initialize pixels
  rowsize = ncomponents * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  if (nbytes > 0) {
    pixels = new unsigned char [nbytes];
    assert(pixels);
    unsigned char *p = pixels;
    while (--nbytes) *(p++) = 0;
  }
}



R2Image::
R2Image(int width, int height, int ncomponents, unsigned char *data)
  : width(width), 
    height(height),
    ncomponents(ncomponents),
    rowsize(0),
    pixels(NULL)
{
  // Copy pixels (first pixel is lower-left)
  rowsize = ncomponents * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  if (nbytes > 0) {
    pixels = new unsigned char [nbytes];
    assert(pixels);
    unsigned char *p = pixels;
    while (--nbytes) *(p++) = *(data++);
  }
}



R2Image::
R2Image(const R2Grid& red, const R2Grid& green, const R2Grid& blue)
  : width(red.XResolution()), 
    height(red.YResolution()),
    ncomponents(3),
    rowsize(0),
    pixels(NULL)
{
  // Just checking
  assert(red.XResolution() != green.XResolution());
  assert(red.YResolution() != green.YResolution());
  assert(red.XResolution() != blue.XResolution());
  assert(red.YResolution() != blue.YResolution());

  // Copy pixels (first pixel is lower-left)
  rowsize = ncomponents * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  if (nbytes > 0) {
    pixels = new unsigned char [nbytes];
    assert(pixels);
    for (int i = 0; i < width*height; i++) {
      int r = 255.0 * red.GridValue(i);
      int g = 255.0 * green.GridValue(i);
      int b = 255.0 * blue.GridValue(i);
      r = (r < 0) ? 0 : ((r > 255) ? 255 : r);
      g = (g < 0) ? 0 : ((g > 255) ? 255 : g);
      b = (b < 0) ? 0 : ((b > 255) ? 255 : b);
      pixels[i*3+0] = r;
      pixels[i*3+1] = g;
      pixels[i*3+2] = b;
    }
  }
}



R2Image::
R2Image(const R2Image& image)
  : width(image.width), 
    height(image.height),
    ncomponents(image.ncomponents),
    rowsize(image.rowsize),
    pixels(NULL)
{
  // Copy pixels
  int nbytes = rowsize * height;
  if (nbytes > 0) {
    pixels = new unsigned char [nbytes];
    assert(pixels);
    unsigned char *p = pixels;
    unsigned char *data = image.pixels;
    while (--nbytes) *(p++) = *(data++);
  }
}



R2Image::
~R2Image(void)
{
  // Free image pixels
  if (pixels) delete [] pixels;
}



const RNRgb R2Image::
PixelRGB(int x, int y) const
{
  // Return pixel RGB at (x,y)
  RNScalar r, g, b;
  const unsigned char *pixel = Pixel(x, y);
  switch (ncomponents) {
  case 1: 
  case 2:
    r = *(pixel++) / 255.0;
    g = b = r;
    break;

  case 3:
  case 4:
    r = *(pixel++) / 255.0;
    g = *(pixel++) / 255.0;
    b = *(pixel++) / 255.0;
    break;

  default:
    r = b = g = 0;
    break;
  }
  return RNRgb(r, g, b);
}



R2Image& R2Image::
operator=(const R2Image& image)
{
  // Assign values
  this->width = image.width;
  this->height = image.height;
  this->ncomponents = image.ncomponents;
  this->rowsize = image.rowsize;

  // Copy pixels
  if (this->pixels) delete [] pixels;
  int nbytes = this->rowsize * this->height;
  this->pixels = new unsigned char [nbytes];
  assert(this->pixels);
  unsigned char *p =this-> pixels;
  unsigned char *data = image.pixels;
  while (--nbytes) *(p++) = *(data++);

  // Return this
  return *this;
}



void R2Image::
Clear(const RNRgb& rgb)
{
  for (int i = 0; i < width; i++) {
    for (int j = 0; j < height; j++) {
      SetPixelRGB(i, j, rgb);
    }
  }
}



void R2Image::
Add(const R2Image& image)
{
  int nbytes = rowsize * height;
  unsigned char *p1 = pixels;
  unsigned char *p2 = image.pixels;
  while (--nbytes) {
    int value = *p1;
    value += *(p2++);
    if (value > 255) value = 255;
    *(p1++) = value;
  }
}



void R2Image::
Subtract(const R2Image& image)
{
  int nbytes = rowsize * height;
  unsigned char *p1 = pixels;
  unsigned char *p2 = image.pixels;
  while (--nbytes) {
    int value = *p1 - *p2 + 128;
    if (value < 0) value = 0;
    if (value > 255) value = 255;
    *p1 = value;
    p1++;
    p2++;
  }
}



void R2Image::
SetPixel(int x, int y, const unsigned char *value)
{
  // Set pixel color
  unsigned char *pixel = &(pixels[y*rowsize + x*ncomponents]);
  for (int i = 0; i < ncomponents; i++) *(pixel++) = *(value++);
}



void R2Image::
SetPixelRGB(int x, int y, const RNRgb& rgb)
{
  // Set pixel color
  unsigned char *pixel = &(pixels[y*rowsize + x*ncomponents]);
  switch (ncomponents) {
  case 1: 
    *(pixel++) = (unsigned char) (255 * rgb.Luminance());
    break;

  case 2:
    *(pixel++) = (unsigned char) (255 * rgb.Luminance());
    *(pixel++) = (unsigned char) 255;
    break;

  case 3:
    *(pixel++) = (unsigned char) (255 * rgb.R());
    *(pixel++) = (unsigned char) (255 * rgb.G());
    *(pixel++) = (unsigned char) (255 * rgb.B());
    break;

  case 4:
    *(pixel++) = (unsigned char) (255 * rgb.R());
    *(pixel++) = (unsigned char) (255 * rgb.G());
    *(pixel++) = (unsigned char) (255 * rgb.B());
    *(pixel++) = (unsigned char) 255;
    break;
  }
}



void R2Image::
Capture(void)
{
  // Check image size
  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);
  assert(width >= viewport[2]);
  assert(height >= viewport[3]);
  assert(pixels);

  // Read pixels from frame buffer 
  switch (ncomponents) {
  case 1: glReadPixels(0, 0, width, height, GL_LUMINANCE, GL_UNSIGNED_BYTE, pixels); break;
  case 2: glReadPixels(0, 0, width, height, GL_LUMINANCE_ALPHA, GL_UNSIGNED_BYTE, pixels); break;
  case 3: glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels); break;
  case 4: glReadPixels(0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, pixels); break;
  }
}



void R2Image::
Draw(int x, int y) const
{
  // Set projection matrix
  glMatrixMode(GL_PROJECTION);  
  glPushMatrix();
  glLoadIdentity();
  gluOrtho2D(0, width, 0, height);

  // Set model view matrix
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  // Set position for image
  glRasterPos2i(x, y);

  // Draw pixels
  switch (ncomponents) {
  case 1: glDrawPixels(width, height, GL_LUMINANCE, GL_UNSIGNED_BYTE, pixels); break;
  case 2: glDrawPixels(width, height, GL_LUMINANCE_ALPHA, GL_UNSIGNED_BYTE, pixels); break;
  case 3: glDrawPixels(width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels); break;
  case 4: glDrawPixels(width, height, GL_RGBA, GL_UNSIGNED_BYTE, pixels); break;
  default: fprintf(stderr, "Unrecognized number of components in image: %d\n", ncomponents); break;
  }

  // Reset projection matrix
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();

  // Reset model view matrix
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}


int R2Image::
Read(const char *filename)
{
  // Initialize everything
  if (pixels) { delete [] pixels; pixels = NULL; }
  width = height = 0;

  // Parse input filename extension
  const char *input_extension;
  if (!(input_extension = strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
    return 0;
  }
  
  // Read file of appropriate type
  if (!strncmp(input_extension, ".bmp", 4)) return ReadBMP(filename);
  else if (!strncmp(input_extension, ".ppm", 4)) return ReadPPM(filename);
  else if (!strncmp(input_extension, ".pgm", 4)) return ReadPPM(filename);
  else if (!strncmp(input_extension, ".pfm", 4)) return ReadPFM(filename);
  else if (!strncmp(input_extension, ".jpg", 4)) return ReadJPEG(filename);
  else if (!strncmp(input_extension, ".jpeg", 5)) return ReadJPEG(filename);
  else if (!strncmp(input_extension, ".tif", 4)) return ReadTIFF(filename);
  else if (!strncmp(input_extension, ".tiff", 5)) return ReadTIFF(filename);
  else if (!strncmp(input_extension, ".png", 4)) return ReadPNG(filename);
  else if (!strncmp(input_extension, ".raw", 4)) return ReadRAW(filename);
  else if (!strncmp(input_extension, ".grd", 4)) return ReadGRD(filename);
  
  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



int R2Image::
Write(const char *filename) const
{
  // Parse input filename extension
  const char *input_extension;
  if (!(input_extension = strrchr(filename, '.'))) {
    fprintf(stderr, "Output file has no extension (e.g., .jpg).\n");
    return 0;
  }
  
  // Write file of appropriate type
  if (!strncmp(input_extension, ".bmp", 4)) return WriteBMP(filename);
  else if (!strncmp(input_extension, ".pgm", 4)) return WritePPM(filename, 0);
  else if (!strncmp(input_extension, ".ppm", 4)) return WritePPM(filename, 0);
  else if (!strncmp(input_extension, ".jpg", 4)) return WriteJPEG(filename);
  else if (!strncmp(input_extension, ".jpeg", 5)) return WriteJPEG(filename);
  else if (!strncmp(input_extension, ".tif", 4)) return WriteTIFF(filename);
  else if (!strncmp(input_extension, ".tiff", 5)) return WriteTIFF(filename);
  else if (!strncmp(input_extension, ".png", 4)) return WritePNG(filename);
  else if (!strncmp(input_extension, ".raw", 4)) return WriteRAW(filename);

  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



////////////////////////////////////////////////////////////////////////
// BMP I/O
////////////////////////////////////////////////////////////////////////

typedef struct BMPBITMAPFILEHEADER {
  unsigned short int bfType;
  unsigned int bfSize;
  unsigned short int bfReserved1;
  unsigned short int bfReserved2;
  unsigned int bfOffBits;
} BMPBITMAPFILEHEADER;

typedef struct BMPBITMAPINFOHEADER {
  unsigned int biSize;
  int biWidth;
  int biHeight;
  unsigned short int biPlanes;
  unsigned short int biBitCount;
  unsigned int biCompression;
  unsigned int biSizeImage;
  int biXPelsPerMeter;
  int biYPelsPerMeter;
  unsigned int biClrUsed;
  unsigned int biClrImportant;
} BMPBITMAPINFOHEADER;

#ifndef BI_BITFIELDS
#define BI_RGB        0L
#define BI_RLE8       1L
#define BI_RLE4       2L
#define BI_BITFIELDS  3L
#endif

#define BMP_BF_TYPE 0x4D42 /* word BM */
#define BMP_BF_OFF_BITS 54 /* 14 for file header + 40 for info header (not sizeof(), but packed size) */
#define BMP_BI_SIZE 40 /* packed size of info header */


/* Reads a unsigned short int from a file in little endian format */
static unsigned short int WordReadLE(FILE *fp)
{
  unsigned short int lsb, msb;
  lsb = getc(fp);
  msb = getc(fp);
  return (msb << 8) | lsb;
}



/* Writes a unsigned short int to a file in little endian format */
static void WordWriteLE(unsigned short int x, FILE *fp)
{
  unsigned char lsb = (unsigned char) (x & 0x00FF); putc(lsb, fp); 
  unsigned char msb = (unsigned char) (x >> 8); putc(msb, fp);
}



/* Reads a unsigned int word from a file in little endian format */
static unsigned int DWordReadLE(FILE *fp)
{
  unsigned int b1 = getc(fp);
  unsigned int b2 = getc(fp);
  unsigned int b3 = getc(fp);
  unsigned int b4 = getc(fp);
  return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



/* Writes a unsigned int to a file in little endian format */
static void DWordWriteLE(unsigned int x, FILE *fp)
{
  unsigned char b1 = (x & 0x000000FF); putc(b1, fp);
  unsigned char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
  unsigned char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
  unsigned char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



/* Reads a int word from a file in little endian format */
static int LongReadLE(FILE *fp)
{
  int b1 = getc(fp);
  int b2 = getc(fp);
  int b3 = getc(fp);
  int b4 = getc(fp);
  return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



/* Writes a int to a file in little endian format */
static void LongWriteLE(int x, FILE *fp)
{
  char b1 = (x & 0x000000FF); putc(b1, fp);
  char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
  char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
  char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



int R2Image::
ReadBMP(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  /* Read file header */
  BMPBITMAPFILEHEADER bmfh;
  bmfh.bfType = WordReadLE(fp);
  bmfh.bfSize = DWordReadLE(fp);
  bmfh.bfReserved1 = WordReadLE(fp);
  bmfh.bfReserved2 = WordReadLE(fp);
  bmfh.bfOffBits = DWordReadLE(fp);
  
  /* Check file header */
  assert(bmfh.bfType == BMP_BF_TYPE);
  /* ignore bmfh.bfSize */
  /* ignore bmfh.bfReserved1 */
  /* ignore bmfh.bfReserved2 */
  assert(bmfh.bfOffBits == BMP_BF_OFF_BITS);
  
  /* Read info header */
  BMPBITMAPINFOHEADER bmih;
  bmih.biSize = DWordReadLE(fp);
  bmih.biWidth = LongReadLE(fp);
  bmih.biHeight = LongReadLE(fp);
  bmih.biPlanes = WordReadLE(fp);
  bmih.biBitCount = WordReadLE(fp);
  bmih.biCompression = DWordReadLE(fp);
  bmih.biSizeImage = DWordReadLE(fp);
  bmih.biXPelsPerMeter = LongReadLE(fp);
  bmih.biYPelsPerMeter = LongReadLE(fp);
  bmih.biClrUsed = DWordReadLE(fp);
  bmih.biClrImportant = DWordReadLE(fp);
  
  // Check info header 
  assert(bmih.biSize == BMP_BI_SIZE);
  assert(bmih.biWidth > 0);
  assert(bmih.biHeight > 0);
  assert(bmih.biPlanes == 1);
  assert(bmih.biBitCount == 24);  /* RGB */
  assert(bmih.biCompression == BI_RGB);   /* RGB */
  int lineLength = bmih.biWidth * 3;  /* RGB */
  if ((lineLength % 4) != 0) lineLength = (lineLength / 4 + 1) * 4;
  assert(bmih.biSizeImage == (unsigned int) lineLength * (unsigned int) bmih.biHeight);

  // Assign width, height, and ncomponents
  width = bmih.biWidth;
  height = bmih.biHeight;
  ncomponents = 3;
  rowsize = ncomponents * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;

  // Allocate pixels
  int nbytes = bmih.biSizeImage;
  pixels = new unsigned char [nbytes];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Read pixels 
  fseek(fp, (long) bmfh.bfOffBits, SEEK_SET);
  if (fread(pixels, 1, bmih.biSizeImage, fp) != bmih.biSizeImage) {
    fprintf(stderr, "Error while reading BMP file %s", filename);
    return 0;
  }

  // Swap blue and red in each pixel
  for (int j = 0; j < height; j++) {
    unsigned char *p = &pixels[j * rowsize];
    for (int i = 0; i < width; i++) {
      unsigned char c = *p;
      *(p) = *(p+2);
      *(p+2) = c;
      p += 3;
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R2Image::
WriteBMP(const char *filename) const
{
  // Check ncomponents
  assert(ncomponents == 3);

  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Write file header 
  BMPBITMAPFILEHEADER bmfh;
  bmfh.bfType = BMP_BF_TYPE;
  bmfh.bfSize = BMP_BF_OFF_BITS + rowsize * height;
  bmfh.bfReserved1 = 0;
  bmfh.bfReserved2 = 0;
  bmfh.bfOffBits = BMP_BF_OFF_BITS;
  WordWriteLE(bmfh.bfType, fp);
  DWordWriteLE(bmfh.bfSize, fp);
  WordWriteLE(bmfh.bfReserved1, fp);
  WordWriteLE(bmfh.bfReserved2, fp);
  DWordWriteLE(bmfh.bfOffBits, fp);

  // Write info header 
  BMPBITMAPINFOHEADER bmih;
  bmih.biSize = BMP_BI_SIZE;
  bmih.biWidth = width;
  bmih.biHeight = height;
  bmih.biPlanes = 1;
  bmih.biBitCount = 24;       /* RGB */
  bmih.biCompression = BI_RGB;    /* RGB */
  bmih.biSizeImage = rowsize * (unsigned int) bmih.biHeight;  /* RGB */
  bmih.biXPelsPerMeter = 2925;
  bmih.biYPelsPerMeter = 2925;
  bmih.biClrUsed = 0;
  bmih.biClrImportant = 0;
  DWordWriteLE(bmih.biSize, fp);
  LongWriteLE(bmih.biWidth, fp);
  LongWriteLE(bmih.biHeight, fp);
  WordWriteLE(bmih.biPlanes, fp);
  WordWriteLE(bmih.biBitCount, fp);
  DWordWriteLE(bmih.biCompression, fp);
  DWordWriteLE(bmih.biSizeImage, fp);
  LongWriteLE(bmih.biXPelsPerMeter, fp);
  LongWriteLE(bmih.biYPelsPerMeter, fp);
  DWordWriteLE(bmih.biClrUsed, fp);
  DWordWriteLE(bmih.biClrImportant, fp);

  // Write image, swapping blue and red in each pixel
  int pad = rowsize - width * 3;
  for (int j = 0; j < height; j++) {
    // Write row of pixels BGR
    unsigned char *p = &pixels[j * rowsize];
    for (int i = 0; i < width; i++) {
      fputc(*(p+2), fp);
      fputc(*(p+1), fp);
      fputc(*(p+0), fp);
      p += 3;
    }

    // Pad row
    for (int i = 0; i < pad; i++) fputc(0, fp);
  }
  
  // Close file
  fclose(fp);

  // Return success
  return 1;  
}



////////////////////////////////////////////////////////////////////////
// PPM I/O
////////////////////////////////////////////////////////////////////////

int R2Image::
ReadPPM(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Read magic identifier
  char buffer[128];
  if (!fgets(buffer, 128, fp)) {
    fprintf(stderr, "Unable to read magic id in PPM file");
    fclose(fp);
    return 0;
  }

  // skip comments
  int c = getc(fp);
  while (c == '#') {
    while (c != '\n') c = getc(fp);
    c = getc(fp);
  }
  ungetc(c, fp);

  // Read width and height
  if (fscanf(fp, "%d%d", &width, &height) != 2) {
    fprintf(stderr, "Unable to read width and height in PPM file");
    fclose(fp);
    return 0;
  }
	
  // Read max value
  int max_value;
  if (fscanf(fp, "%d", &max_value) != 1) {
    fprintf(stderr, "Unable to read max_value in PPM file");
    fclose(fp);
    return 0;
  }
	
  // Assign ncomponents and rowsize
  ncomponents = 3;
  rowsize = ncomponents * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;

  // Allocate image pixels
  int nbytes = rowsize * height;
  pixels = new unsigned char [nbytes];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for PPM file");
    fclose(fp);
    return 0;
  }

  // Check file type
  if (!strcmp(buffer, "P2\n")) {
    // Read asci gray image data 
    for (int j = height-1; j >= 0; j--) {
      unsigned char *p = &pixels[j*rowsize];
      for (int i = 0; i < width; i++) {
	// Read pixel values
	int value;
	if (fscanf(fp, "%d", &value) != 1) {
	  fprintf(stderr, "Unable to read data at (%d,%d) in PGM file", i, j);
	  fclose(fp);
	  return 0;
	}

	// Assign pixel values
	*(p++) = ((255 * value) / max_value) & 0xFF;
	*(p++) = ((255 * value) / max_value) & 0xFF;
	*(p++) = ((255 * value) / max_value) & 0xFF;
      }
    }
  }
  else if (!strcmp(buffer, "P3\n")) {
    // Read asci rgb image data 
    for (int j = height-1; j >= 0; j--) {
      unsigned char *p = &pixels[j*rowsize];
      for (int i = 0; i < width; i++) {
	// Read pixel values
	int r, g, b;
	if (fscanf(fp, "%d%d%d", &r, &g, &b) != 3) {
	  fprintf(stderr, "Unable to read data at (%d,%d) in PPM file", i, j);
	  fclose(fp);
	  return 0;
	}

	// Assign pixel values
	*(p++) = ((255 * r) / max_value) & 0xFF;
	*(p++) = ((255 * g) / max_value) & 0xFF;
	*(p++) = ((255 * b) / max_value) & 0xFF;
      }
    }
  }
  else if (!strcmp(buffer, "P5\n")) {
    // Read up to one character of whitespace (\n) after max_value
    int c = getc(fp);
    if (!isspace(c)) putc(c, fp);

    // Read binary gray image data 
    // First pgm pixel is top-left, so read in opposite scan-line order
    for (int j = height-1; j >= 0; j--) {
      unsigned char *p = &pixels[j*rowsize];
      for (int i = 0; i < width; i++) {
	// Read pixel values
        unsigned char value = ((255 * getc(fp)) / max_value) & 0xFF;
	*(p++) = value;
	*(p++) = value;
	*(p++) = value;
      }
    }
  }
  else if (!strcmp(buffer, "P6\n")) {
    // Read up to one character of whitespace (\n) after max_value
    int c = getc(fp);
    if (!isspace(c)) putc(c, fp);

    // Read binary rgb image data 
    // First ppm pixel is top-left, so read in opposite scan-line order
    for (int j = height-1; j >= 0; j--) {
      unsigned char *p = &pixels[j*rowsize];
      for (int i = 0; i < width; i++) {
	// Read pixel values
	*(p++) = ((255 * getc(fp)) / max_value) & 0xFF;
	*(p++) = ((255 * getc(fp)) / max_value) & 0xFF;
	*(p++) = ((255 * getc(fp)) / max_value) & 0xFF;
      }
    }
  }
  else {
    fprintf(stderr, "Unrecognized magic header in pbm file: %s\n", buffer);
    fclose(fp);
    return 0;
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R2Image::
WritePPM(const char *filename, int ascii) const
{
  // Check ncomponents
  assert(ncomponents == 3);

  // Check pgm/ppm
  int pgm = (strstr(filename, ".pgm")) ? 1 : 0;

  // Check type
  if (ascii) {
    // Open file
    FILE *fp = fopen(filename, "w");
    if (!fp) {
      fprintf(stderr, "Unable to open image file: %s", filename);
      return 0;
    }

    // Print PPM image file 
    // First ppm pixel is top-left, so write in opposite scan-line order
    if (pgm) fprintf(fp, "P2\n");
    else fprintf(fp, "P3\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");
    for (int j = height-1; j >= 0 ; j--) {
      for (int i = 0; i < width; i++) {
        const unsigned char *p = Pixel(i, j);
        int r = *(p++) & 0xFF;
        int g = *(p++) & 0xFF;
        int b = *(p++) & 0xFF;
        if (pgm) fprintf(fp, "%-3d  ", (int) (0.30*r + 0.59*g + 0.11*b));
        else fprintf(fp, "%-3d %-3d %-3d  ", r, g, b);
        if (((i+1) % 4) == 0) fprintf(fp, "\n");
      }
      if ((width % 4) != 0) fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    // Close file
    fclose(fp);
  }
  else {
    // Open file
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
      fprintf(stderr, "Unable to open image file: %s", filename);
      return 0;
    }
    
    // Print PPM image file 
    // First ppm pixel is top-left, so write in opposite scan-line order
    if (pgm) fprintf(fp, "P5\n");
    else fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");
    for (int j = height-1; j >= 0 ; j--) {
      for (int i = 0; i < width; i++) {
        const unsigned char *p = Pixel(i, j);
        int r = *(p++) & 0xFF;
        int g = *(p++) & 0xFF;
        int b = *(p++) & 0xFF;
        if (pgm) fprintf(fp, "%c", (int) (0.30*r + 0.59*g + 0.11*b));
        else fprintf(fp, "%c%c%c", r, g, b);
      }
    }
    
    // Close file
    fclose(fp);
  }

  // Return number of bytes written
  return width * height * ncomponents;  
}



////////////////////////////////////////////////////////////////////////
// PPM I/O
////////////////////////////////////////////////////////////////////////

int R2Image::
ReadPFM(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image: %s\n", filename);
    return 0;
  }

  // Read magic header
  int c;
  c = fgetc(fp); if (c != 'P') { fprintf(stderr, "Bad magic keyword in %s\n", filename); return 0; }
  c = fgetc(fp); if (c != 'f') { fprintf(stderr, "Bad magic keyword in %s\n", filename); return 0; }
  c = fgetc(fp); if (c != '\n') { fprintf(stderr, "Bad magic keyword in %s\n", filename); return 0; }

  // Read width
  int width_count = 0;
  char width_string[256];
  for (int i = 0; i < 256; i++) { 
    c = fgetc(fp); 
    if ((c == ' ') && (width_count == 0)) { continue; }
    else if ((c == ' ') || (c == '\n')) { width_string[width_count] = '\0'; break; }
    else if (!isdigit(c)) { fprintf(stderr, "Bad width character %c in %s\n", c, filename); return 0; }
    else width_string[width_count++] = c;
  }

  // Check width
  if ((width_count == 0) || (width_count > 128)) {
    fprintf(stderr, "Error reading width in %s\n", filename); 
    return 0; 
  }

  // Read height
  int height_count = 0;
  char height_string[256];
  for (int i = 0; i < 256; i++) { 
    c = fgetc(fp); 
    if ((c == ' ') && (height_count == 0)) { continue; }
    else if ((c == ' ') || (c == '\n')) { height_string[height_count] = '\0'; break; }
    else if (!isdigit(c)) { fprintf(stderr, "Bad height character %c in %s\n", c, filename); return 0; }
    else height_string[height_count++] = c;
  }

  // Check height
  if ((height_count == 0) || (height_count > 128)) {
    fprintf(stderr, "Error reading height in %s\n", filename); 
    return 0; 
  }

  // Read endian
  int endian_count = 0;
  char endian_string[256];
  for (int i = 0; i < 256; i++) { 
    c = fgetc(fp); 
    if ((c == ' ') && (endian_count == 0)) { continue; }
    else if ((c == ' ') || (c == '\n')) { endian_string[endian_count] = '\0'; break; }
    if (!isdigit(c) && (c != '.') && (c != '-')) { fprintf(stderr, "Bad endian character %c in %s\n", c, filename); return 0; }
    endian_string[endian_count++] = c;
  }

  // Check endian
  if ((endian_count == 0) || (endian_count > 128)) {
    fprintf(stderr, "Error reading endian in %s\n", filename); 
    return 0; 
  }

  // Parse values
  width = atoi(width_string);
  height = atoi(height_string);
  float endian = atof(endian_string);
  if (endian == -999.0F) fprintf(stderr, "Just trying to avoid compiler warning for unused variable\n");

  // Allocate values
  int npixels = width * height;
  float *values = new float [ npixels ];
  if (!values) {
    fprintf(stderr, "Unable to allocate values for %s\n", filename);
    return 0;
  }

  // Read values
  int count = npixels;
  while (count > 0) {
    int n = (int) fread(&values[npixels-count], sizeof(float), count, fp);
    if (n <= 0) { fprintf(stderr, "read error\n"); abort(); }
    count -= n;
  }

  // Close image
  fclose(fp);

  // Determine offset and scale
  float offset = 0;
  float scale = 1;
  if (npixels > 0) {
    float sum = 0;
    float min_value = FLT_MAX;
    float max_value = -FLT_MAX;
    for (int i = 0; i < npixels; i++) {
      if (values[i] == R2_GRID_UNKNOWN_VALUE) continue;
      if (values[i] > max_value) max_value = values[i];
      if (values[i] < min_value) min_value = values[i];
      sum += values[i];
    }
    float mean = sum / npixels;
    float ssd = 0;
    for (int i = 0; i < npixels; i++) {
      if (values[i] == R2_GRID_UNKNOWN_VALUE) continue;
      float delta = values[i] - mean;
      ssd += delta * delta;
    }
    float stddev = sqrt(ssd);
    float min_clamped = mean - 3 * stddev;
    if (min_clamped < min_value) min_clamped = min_value;
    float max_clamped = mean + 3 * stddev;
    if (max_clamped > max_value) max_clamped = max_value;
    float range = max_clamped - min_clamped;
    if (range > 0) {
      offset = min_clamped;
      scale = 1.0F / range;
    }
  }

  // Allocate pixels
  ncomponents = 3;
  rowsize = ncomponents * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  pixels = new unsigned char [nbytes];
  assert(pixels);

  // Fill in pixels
  float *valuesp = values;
  for (int j = 0; j < height; j++) {
    unsigned char *pixelsp = &pixels[j*rowsize];
    for (int i = 0; i < width; i++) {
      if (*valuesp == R2_GRID_UNKNOWN_VALUE) {
        *(pixelsp++) = 64;
        *(pixelsp++) = 32;
        *(pixelsp++) = 0;
      }
      else {
        int value = (int) ((*valuesp - offset) * 255 * scale);
        if (value < 0) value = 0;
        if (value > 255) value = 255;
        unsigned char p = (unsigned char) value;
        *(pixelsp++) = p;
        *(pixelsp++) = p;
        *(pixelsp++) = p;
      }
      valuesp++;
    }
  }

  // Return success
  return nbytes;
}


////////////////////////////////////////////////////////////////////////
// jpeg I/O
////////////////////////////////////////////////////////////////////////

#define RN_USE_JPEG

#ifdef RN_NO_JPEG
#undef RN_USE_JPEG
#endif



#ifdef RN_USE_JPEG
  extern "C" { 
#   define XMD_H // Otherwise, a conflict with INT32
#   if RN_OS == RN_WINDOWS
#     define HAVE_BOOLEAN
#     undef FAR // Otherwise, a conflict with windows.h
#   endif
#   include "jpeg/jpeglib.h"
  };
#endif



int R2Image::
ReadJPEG(const char *filename)
{
#ifdef RN_USE_JPEG
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Initialize decompression info
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_decompress(&cinfo);
  jpeg_stdio_src(&cinfo, fp);
  jpeg_read_header(&cinfo, TRUE);
  jpeg_start_decompress(&cinfo);

  // Remember image attributes
  width = cinfo.output_width;
  height = cinfo.output_height;
  ncomponents = cinfo.output_components;
  rowsize = ncomponents * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;

  // Allocate image pixels
  int nbytes = rowsize * height;
  this->pixels = new unsigned char [nbytes];
  if (!this->pixels) {
    fprintf(stderr, "Unable to allocate memory for JPEG file");
    fclose(fp);
    return 0;
  }

  // Read scan lines 
  // First jpeg pixel is top-left, so read pixels in opposite scan-line order
  while (cinfo.output_scanline < cinfo.output_height) {
    int scanline = cinfo.output_height - cinfo.output_scanline - 1;
    unsigned char *row_pointer = &pixels[scanline * rowsize];
    jpeg_read_scanlines(&cinfo, &row_pointer, 1);
  }

  // Free everything
  jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);

  // Close file
  fclose(fp);

  // Return success
  return 1;
#else
  RNFail("JPEG not supported");
  return 0;
#endif
}


	

int R2Image::
WriteJPEG(const char *filename) const
{
#ifdef RN_USE_JPEG
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Initialize compression info
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
  jpeg_stdio_dest(&cinfo, fp);
  cinfo.image_width = width; 	/* image width and height, in pixels */
  cinfo.image_height = height;
  cinfo.input_components = ncomponents;		/* # of color components per pixel */
  if (ncomponents == 1) cinfo.in_color_space = JCS_GRAYSCALE;
  else if (ncomponents == 3) cinfo.in_color_space = JCS_RGB;
  else if (ncomponents == 4) cinfo.in_color_space = JCS_CMYK;
  else cinfo.in_color_space = JCS_UNKNOWN;
  cinfo.dct_method = JDCT_ISLOW;
  jpeg_set_defaults(&cinfo);
  cinfo.optimize_coding = TRUE;
  jpeg_set_quality(&cinfo, 75, TRUE);
  jpeg_start_compress(&cinfo, TRUE);
	
  // Output scan lines
  // First jpeg pixel is top-left, so write in opposite scan-line order
  while (cinfo.next_scanline < cinfo.image_height) {
    int scanline = cinfo.image_height - cinfo.next_scanline - 1;
    unsigned char *row_pointer = &pixels[scanline * rowsize];
    jpeg_write_scanlines(&cinfo, &row_pointer, 1);
  }

  // Free everything
  jpeg_finish_compress(&cinfo);
  jpeg_destroy_compress(&cinfo);

  // Close file
  fclose(fp);

  // Return number of bytes written
  return 1;
#else
  RNFail("JPEG not supported");
  return 0;
#endif
}



////////////////////////////////////////////////////////////////////////
// TIFF I/O
////////////////////////////////////////////////////////////////////////

#ifdef RN_NO_TIFF
#undef RN_USE_TIFF
#endif



#ifdef RN_USE_TIFF
# include "tiff/tiffio.h"
#endif



int R2Image::
ReadTIFF(const char *filename)
{
#ifdef RN_USE_TIFF
  // Open file
  TIFF* tif = TIFFOpen(filename, "r");
  if (!tif) {
    fprintf(stderr, "Unable to open TIFF file %s\n", filename);
    return 0;
  }

  // Get image dimensions
  uint32 w, h;
  TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
  TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
  size_t npixels = w * h;

  // Allocate buffer for data
  uint32* raster = (uint32*) _TIFFmalloc(npixels * sizeof (uint32));
  if (!raster) {
    fprintf(stderr, "Unable to allocate data for TIFF file %s\n", filename);
    return 0;
  }

  // Read data into buffer
  if (!TIFFReadRGBAImage(tif, w, h, raster, 0)) {
    fprintf(stderr, "Unable to read TIFF file %s\n", filename);
    return 0;
  }

  // Initialize R2Image data
  width = w;
  height = h;
  ncomponents = 3;
  rowsize = ncomponents * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  pixels = new unsigned char [nbytes];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for TIFF file %s", filename);
    return 0;
  }

  // Fill R2Image pixel data 
  uint32 *rasterp = raster;
  for (int j = 0; j < height; j++) {
    unsigned char *p = &pixels[j*rowsize];
    for (int i = 0; i < width; i++) {
      uint32 pixel = *(rasterp++);
      *(p++) = pixel & 0xFF;
      *(p++) = (pixel >> 8) & 0xFF;
      *(p++) = (pixel >> 16) & 0xFF;
    }
  }

  // Free data
  _TIFFfree(raster);

  // Close file
  TIFFClose(tif);

  // Return success
  return 1;
#else
  RNFail("TIFF not supported");
  return 0;
#endif
}

	

int R2Image::
WriteTIFF(const char *filename) const
{
#ifdef RN_USE_TIFF
  // Open TIFF file
  TIFF *out = TIFFOpen(filename, "w");
  if (!out) {
    fprintf(stderr, "Unable to open TIFF file %s\n", filename);
    return 0;
  }

  // Set TIFF parameters
  TIFFSetField(out, TIFFTAG_IMAGEWIDTH, width);
  TIFFSetField(out, TIFFTAG_IMAGELENGTH, height);
  TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_BOTLEFT);
  TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, 3);
  TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 8);
  TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
  TIFFSetField(out, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
  TIFFSetField(out, TIFFTAG_ROWSPERSTRIP, 1);

  // Allocate data for scan lines
  int scanline_size = TIFFScanlineSize(out);
  unsigned char *buf = (unsigned char *)_TIFFmalloc(scanline_size);
  if (!buf) {
    fprintf(stderr, "Unable to allocate memory for TIFF scan lines\n");
    return 0;
  }

  // Write scan lines to TIFF file
  for (int row = 0; row < height; row++) {
    const unsigned char *p = Pixels(row);
    if (TIFFWriteScanline(out, (tdata_t) p, height - row - 1, 0) < 0) {
      fprintf(stderr, "Unable to write scanline to TIFF image %s\n", filename);
      return 0;
    }
  }

  // Free data for scan lines
  _TIFFfree(buf);

  // Close TIFF file
  TIFFClose(out);

  // Return success
  return 1;
#else
  RNFail("TIFF not supported");
  return 0;
#endif
}



////////////////////////////////////////////////////////////////////////
// PNG I/O
////////////////////////////////////////////////////////////////////////

#define RN_USE_PNG
#ifdef RN_NO_PNG
#undef RN_USE_PNG
#endif



#ifdef RN_USE_PNG
# include "png/png.h"
#endif



int R2Image::
ReadPNG(const char *filename)
{
#ifdef RN_USE_PNG
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open PNG file %s\n", filename);
    return 0;
   }

  // Create and initialize the png_struct 
  png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (png_ptr == NULL) {
    fclose(fp);
    return 0;
  }

  // Allocate/initialize the memory for image information. 
  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (info_ptr == NULL) {
    fclose(fp);
    png_destroy_read_struct(&png_ptr, NULL, NULL);
    return 0;
  }

  // Set up the input control if you are using standard C streams 
  png_init_io(png_ptr, fp);

  // Read the png info 
  png_read_info(png_ptr, info_ptr);

  // Extract image info 
  png_byte color_type = png_get_color_type(png_ptr, info_ptr);
  width = png_get_image_width(png_ptr, info_ptr);
  height = png_get_image_height(png_ptr, info_ptr);
  rowsize = png_get_rowbytes(png_ptr, info_ptr);
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;

  // Set ncomponents
  ncomponents = 0;
  if (color_type == PNG_COLOR_TYPE_GRAY) ncomponents = 1;
  else if (color_type == PNG_COLOR_TYPE_GRAY_ALPHA) ncomponents = 2;
  else if (color_type == PNG_COLOR_TYPE_RGB) ncomponents = 3;
  else if (color_type == PNG_COLOR_TYPE_RGB_ALPHA) ncomponents = 4;
  else { 
    fclose(fp);
    png_destroy_read_struct(&png_ptr, NULL, NULL);
    return 0;
  }

  // Allocate the pixels and row pointers 
  pixels = new unsigned char [ height * rowsize ]; 
  png_bytep *row_pointers = (png_bytep *) png_malloc(png_ptr, height * png_sizeof(png_bytep));
  for (int i = 0; i < height; i++) row_pointers[i] = &pixels[ (height - i - 1) * rowsize ];

  // Read the pixels 
  png_read_image(png_ptr, row_pointers);

  // Finish reading 
  png_read_end(png_ptr, info_ptr);

  // Free the row pointers 
  png_free(png_ptr, row_pointers);

  // Clean up after the read, and free any memory allocated  
  png_destroy_read_struct(&png_ptr, &info_ptr, NULL);

  // Close the file 
  fclose(fp);

  // Return success 
  return 1;
#else
  RNFail("PNG not supported");
  return 0;
#endif
}



int R2Image::
WritePNG(const char *filename) const
{
#ifdef RN_USE_PNG
  // Open the file 
  FILE *fp = fopen(filename, "wb");
  if (fp == NULL) {
    fprintf(stderr, "Unable to open PNG file %s\n", filename);
    return 0;
  }
  
  // Create and initialize the png_struct 
  png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (png_ptr == NULL) {
    fclose(fp);
    return 0;
  }
  
  // Allocate/initialize the image information data. 
  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (info_ptr == NULL) {
    fclose(fp);
    png_destroy_write_struct(&png_ptr,  NULL);
    return 0;
  }
  
  // Determine color type  
  png_byte color_type = 0;
  if (ncomponents == 1) color_type = PNG_COLOR_TYPE_GRAY;
  else if (ncomponents == 2) color_type = PNG_COLOR_TYPE_GRAY_ALPHA;
  else if (ncomponents == 3) color_type = PNG_COLOR_TYPE_RGB;
  else if (ncomponents == 4) color_type = PNG_COLOR_TYPE_RGB_ALPHA;
  else { fprintf(stderr, "Invalid number of components for %s\n", filename); return 0; }

  // Fill in the image data 
  png_set_IHDR(png_ptr, info_ptr, width, height,
    8, color_type, PNG_INTERLACE_NONE,
    PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

  // Allocate the row pointers 
  png_bytep *row_pointers = (png_bytep *) png_malloc(png_ptr, height * png_sizeof(png_bytep));
  for (int i = 0; i < height; i++) row_pointers[i] = &pixels[(height - i - 1) * rowsize];
  
  // Set up the output control 
  png_init_io(png_ptr, fp);
  
  // Write the png info 
  png_write_info(png_ptr, info_ptr);
  
  // Write the pixels 
  png_write_image(png_ptr, row_pointers);
  
  // Finish writing 
  png_write_end(png_ptr, info_ptr);
  
  // Free the row pointers 
  png_free(png_ptr, row_pointers);

  // Clean up after the write, and free any memory allocated 
  png_destroy_write_struct(&png_ptr, &info_ptr);

  // Close the file 
  fclose(fp);

  // Return success
  return 1;
#else
  RNFail("PNG not supported");
  return 0;
#endif
}



////////////////////////////////////////////////////////////////////////
// RAW I/O
////////////////////////////////////////////////////////////////////////

int R2Image::
ReadRAW(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open raw image file %s", filename);
    return 0;
  }

  // Read header
  unsigned int header[4];
  if (fread(header, sizeof(unsigned int), 4, fp) != 4) {
    fprintf(stderr, "Unable to read header to raw image file %s", filename);
    return 0;
  }

  // Check magic number
  static const unsigned int magic = 54321;
  if (header[0] != magic) {
    fprintf(stderr, "Invalid header in raw image file %s", filename);
    return 0;
  }

  // Parse header
  width = header[1];
  height = header[2];
  ncomponents = header[3];
  rowsize = ncomponents * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;

  // Allocate image pixels
  int nbytes = rowsize * height;
  pixels = new unsigned char [nbytes];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for RAW file");
    return 0;
  }

  // Allocate buffer
  float *buf = new float [ width ];
  if (!buf) {
    fprintf(stderr, "Unable to allocate buffer for raw image file %s", filename);
    return 0;
  }

  // Read data for each component
  for (int i = 0; i < ncomponents; i++) {
    for (int j = 0; j < height; j++) {
      // Read buffer for one row
      if (fread(buf, sizeof(float), width, fp) != (unsigned int) width) {
        fprintf(stderr, "Unable to read data from raw image file %s", filename);
        return 0;
      }

      // Parse row
      unsigned char *pixelsp = (unsigned char *) Pixels(j);
      for (int k = 0; k < width; k++) {
        for (int c = 0; c < i; c++) pixelsp++;
        *(pixelsp++) = (int) (buf[k] * 255);
        for (int c = i+1; c < ncomponents; c++) pixelsp++;
      }
    }

  }

  // Free buffer
  delete [] buf;

  // Close file
  fclose(fp);

  // Return number of bytes read
  return width * height * ncomponents * sizeof(float);  
}



int R2Image::
WriteRAW(const char *filename) const
{
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open raw image file %s", filename);
    return 0;
  }

  // Write header
  static const unsigned int magic = 54321;
  unsigned int header[4] = { magic, (unsigned int) width, (unsigned int) height, (unsigned int) ncomponents };
  if (fwrite(header, sizeof(unsigned int), 4, fp) != 4) {
    fprintf(stderr, "Unable to write header to raw image file %s", filename);
    return 0;
  }

  // Allocate buffer
  float *buf = new float [ width ];
  if (!buf) {
    fprintf(stderr, "Unable to allocate buffer for raw image file %s", filename);
    return 0;
  }

  // Write data for each component
  for (int i = 0; i < ncomponents; i++) {
    for (int j = 0; j < height; j++) {
      // Load buffer
      const unsigned char *pixelsp = Pixels(j);
      for (int k = 0; k < width; k++) {
        for (int c = 0; c < i; c++) pixelsp++;
        buf[k] = *(pixelsp++) / 255.0F;
        for (int c = i+1; c < ncomponents; c++) pixelsp++;
      }

      // Write buffer
      if (fwrite(buf, sizeof(float), width, fp) != (unsigned int) width) {
        fprintf(stderr, "Unable to write data to raw image file %s", filename);
        return 0;
      }
    }
  }

  // Free buffer
  delete [] buf;

  // Close file
  fclose(fp);

  // Return number of bytes written
  return width * height * ncomponents * sizeof(float);  
}



////////////////////////////////////////////////////////////////////////
// GRD I/O
////////////////////////////////////////////////////////////////////////

int R2Image::
ReadGRD(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    RNFail("Unable to open file %s", filename);
    return 0;
  }

  // Read grid resolution from file
  int grid_resolution[2];
  if (fread(&grid_resolution, sizeof(int), 2, fp) != 2) {
    RNFail("Unable to read resolution from grid file");
    return 0;
  }

  // Read world_to_grid transformation from file
  RNScalar m[9];
  if (fread(m, sizeof(RNScalar), 9, fp) != 9) {
    RNFail("Unable to read transformation matrix from file");
    return 0;
  }

  // Allocate grid values
  int grid_size = grid_resolution[0] * grid_resolution[1];
  RNScalar *grid_values = new RNScalar [ grid_size ];
  assert(grid_values);

  // Read grid values
  if (fread(grid_values, sizeof(RNScalar), grid_size, fp) != (unsigned int) grid_size) {
    RNFail("Unable to read %d grid values from file", grid_size);
    return 0;
  }

  // Close file
  fclose(fp);

  // Determine offset and scale
  RNScalar offset = 0;
  RNScalar scale = 1;
  if (grid_size > 0) {
    RNScalar sum = 0;
    RNScalar min_value = FLT_MAX;
    RNScalar max_value = -FLT_MAX;
    for (int i = 0; i < grid_size; i++) {
      if (grid_values[i] == R2_GRID_UNKNOWN_VALUE) continue;
      if (grid_values[i] > max_value) max_value = grid_values[i];
      if (grid_values[i] < min_value) min_value = grid_values[i];
      sum += grid_values[i];
    }
    RNScalar mean = sum / grid_size;
    RNScalar ssd = 0;
    for (int i = 0; i < grid_size; i++) {
      if (grid_values[i] == R2_GRID_UNKNOWN_VALUE) continue;
      RNScalar delta = grid_values[i] - mean;
      ssd += delta * delta;
    }
    RNScalar stddev = sqrt(ssd);
    RNScalar min_clamped = mean - 3 * stddev;
    if (min_clamped < min_value) min_clamped = min_value;
    RNScalar max_clamped = mean + 3 * stddev;
    if (max_clamped > max_value) max_clamped = max_value;
    RNScalar range = max_clamped - min_clamped;
    if (range > 0) {
      offset = min_clamped;
      scale = 1.0 / range;
    }
  }

  // Allocate pixels
  ncomponents = 3;
  width = grid_resolution[0];
  height = grid_resolution[1];
  rowsize = ncomponents * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  pixels = new unsigned char [nbytes];
  assert(pixels);

  // Fill in pixels
  RNScalar *grid_valuesp = grid_values;
  for (int j = 0; j < height; j++) {
    unsigned char *pixelsp = &pixels[j*rowsize];
    for (int i = 0; i < width; i++) {
      int value = (int) ((*(grid_valuesp++) - offset) * 255.0 * scale);
      if (value < 0) value = 0;
      if (value > 255) value = 255;
      unsigned char p = (unsigned char) value;
      *(pixelsp++) = p;
      *(pixelsp++) = p;
      *(pixelsp++) = p;
    }
  }

  // Delete grid values
  delete [] grid_values;

  // Return number of bytes
  return nbytes;
}















