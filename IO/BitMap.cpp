// This file is not used in the Lua version of Quanty.
// Older C based routines can use functions in this file
// Do not add changes here, this file will be removed in the future



uint8_t Red(uint32_t RGB)
{
  return (uint8_t)(RGB>>16);
}
uint8_t Green(uint32_t RGB)
{
  return (uint8_t)(RGB>>8);
}
uint8_t Blue(uint32_t RGB)
{
  return (uint8_t)(RGB);
}

typedef struct                       /**** BMP file header structure ****/
{
  unsigned short bfType;        /* Magic number for file */
  unsigned int bfSize;          /* Size of file */
  unsigned short bfReserved1;   /* Reserved */
  unsigned short bfReserved2;   /* Reserved */
  unsigned int bfOffBits;       /* Offset to bitmap data */
} BitmapFileHeader;

typedef struct                       /**** BMP file info structure ****/
{
  unsigned int biSize;          /* Size of info header */
  int    biWidth;               /* Width of image */
  int    biHeight;              /* Height of image */
  unsigned short biPlanes;      /* Number of color planes */
  unsigned short biBitCount;    /* Number of bits per pixel */
  unsigned int biCompression;   /* Type of compression to use */
  unsigned int biSizeImage;     /* Size of image data */
  int    biXPelsPerMeter;       /* X pixels per meter */
  int    biYPelsPerMeter;       /* Y pixels per meter */
  unsigned int biClrUsed;       /* Number of colors used */
  unsigned int biClrImportant;  /* Number of important colors */
} BitmapInfoHeader;


typedef struct                       /**** Colormap entry structure ****/
{
  unsigned char rgbBlue;        /* Blue value */
  unsigned char rgbGreen;       /* Green value */
  unsigned char rgbRed;         /* Red value */
  unsigned char rgbReserved;    /* Reserved */
} RGBQuad;

typedef struct                       /**** Bitmap information structure ****/
{
  BitmapFileHeader bmfHeader;   /* File header */
  BitmapInfoHeader bmiHeader;   /* Image header */
  RGBQuad bmiColors[256];       /* Image colormap */
} BitmapInfo;

/*
 * 'read_word()' - Read a 16-bit unsigned integer.
 */
static unsigned short read_word(FILE * fp)
{
  unsigned char b0, b1;         /* Bytes from file */
  b0 = getc(fp);
  b1 = getc(fp);
  return ((b1 << 8) | b0);
}


/*
 * 'read_dword()' - Read a 32-bit unsigned integer.
 */
static unsigned int read_dword(FILE * fp)
{
  unsigned char b0, b1, b2, b3; /* Bytes from file */
  b0 = getc(fp);
  b1 = getc(fp);
  b2 = getc(fp);
  b3 = getc(fp);
  return ((((((b3 << 8) | b2) << 8) | b1) << 8) | b0);
}


/*
 * 'read_long()' - Read a 32-bit signed integer.
 */
static int read_long(FILE * fp)
{
  unsigned char b0, b1, b2, b3; /* Bytes from file */
  b0 = getc(fp);
  b1 = getc(fp);
  b2 = getc(fp);
  b3 = getc(fp);
  return ((int)(((((b3 << 8) | b2) << 8) | b1) << 8) | b0);
}


/*
 * 'write_word()' - Write a 16-bit unsigned integer.
 */
static int write_word(FILE * fp, unsigned short w)
{
  putc(w, fp);
  return (putc(w >> 8, fp));
}


/*
 * 'write_dword()' - Write a 32-bit unsigned integer.
 */
static int write_dword(FILE * fp, unsigned int dw)
{
  putc(dw, fp);
  putc(dw >> 8, fp);
  putc(dw >> 16, fp);
  return (putc(dw >> 24, fp));
}


/*
 * 'write_long()' - Write a 32-bit signed integer.
 */
static int write_long(FILE * fp, int l)
{
  putc(l, fp);
  putc(l >> 8, fp);
  putc(l >> 16, fp);
  return (putc(l >> 24, fp));
}

/*  A function that transforms a real value to an RGB color value
 *  This function maps values between 0 and 1 to a grey scale with 0 black and 1 white
 */
void RealToRGBColorBW(double V, unsigned char *RGB)
{
  if (V > 1)
  {
    RGB[0] = 255;
    RGB[1] = 255;
    RGB[2] = 255;
    return;
  }
  if (V < 0)
  {
    RGB[0] = 0;
    RGB[1] = 0;
    RGB[2] = 0;
    return;
  }
  RGB[0] = (uint8_t) (255 * V);
  RGB[1] = (uint8_t) (255 * V);
  RGB[2] = (uint8_t) (255 * V);
  return;
}

/*  A function that transforms a real value to an RGB color value
 *  This function maps values between 0 and 1 to a grey scale with 1 black and 0 white
 */
void RealToRGBColorBWInverse(double V, unsigned char *RGB)
{
  if (V < 0)
  {
    RGB[0] = 255;
    RGB[1] = 255;
    RGB[2] = 255;
    return;
  }
  if (V > 1)
  {
    RGB[0] = 0;
    RGB[1] = 0;
    RGB[2] = 0;
    return;
  }
  RGB[0] = (uint8_t) (255 * (1 - V));
  RGB[1] = (uint8_t) (255 * (1 - V));
  RGB[2] = (uint8_t) (255 * (1 - V));
  return;
}

/*  A function that transforms three real values to an RGB color value
 *  This function maps V[0] to the Red channel, V[1] to the Green channel, V[2] to the Blue channel
 *  whith 1 meaning black and 0 meaning the collor
 */
void Real3ToRGBColorSimpleInverse(double V[3], unsigned char *RGB)
{
  if (V[0] < 0)
  {
    RGB[0] = 255;
  }
  else if (V[0] > 1)
  {
    RGB[0] = 0;
  }
  else
  {
    RGB[0] = (uint8_t) (255 * (1 - V[0]));
  }
  if (V[1] < 0)
  {
    RGB[1] = 255;
  }
  else if (V[1] > 1)
  {
    RGB[1] = 0;
  }
  else
  {
    RGB[1] = (uint8_t) (255 * (1 - V[1]));
  }
  if (V[2] < 0)
  {
    RGB[2] = 255;
  }
  else if (V[2] > 1)
  {
    RGB[2] = 0;
  }
  else
  {
    RGB[2] = (uint8_t) (255 * (1 - V[2]));
  }
  return;
}

/*  Print a Matrix to a stream as a bitmap with collor coding acording to the function RealToRGBColor
 *
 */
int StreamPrintRealMatrixToBitmap(MatrixType M, FILE * stream, void (*RealToRGBColor) (double, unsigned char *))
{
  BitmapInfo info;              /* I - Bitmap information */
  uint32_t i, j;
  uint8_t RGB[3];
  unsigned int biSizeImage;

  biSizeImage = (unsigned int)((M.Ni * 3 + ((4 - ((M.Ni * 3) & 3)) & 3)) * M.Nj);

  info.bmfHeader.bfType = 0x4D42; /* unsigned short bfType;           Magic number for file */
  info.bmfHeader.bfSize = 54 + biSizeImage; /*  unsigned int   bfSize;          Size of file */
  info.bmfHeader.bfReserved1 = 0; /*  unsigned short bfReserved1;     Reserved */
  info.bmfHeader.bfReserved2 = 0; /*  unsigned short bfReserved2;     ... */
  info.bmfHeader.bfOffBits = 54;  /*  unsigned int   bfOffBits;       Offset to bitmap data */

  info.bmiHeader.biSize = 40;   /* unsigned int   biSize             Size of info header */
  info.bmiHeader.biWidth = (int)M.Ni; /* int            biWidth;           Width of image */
  info.bmiHeader.biHeight = (int)M.Nj;  /* int            biHeight;          Height of image */
  info.bmiHeader.biPlanes = 1;  /* unsigned short biPlanes;          Number of color planes */
  info.bmiHeader.biBitCount = 24; /* unsigned short biBitCount;        Number of bits per pixel */
  info.bmiHeader.biCompression = 0; /* unsigned int   biCompression;     Type of compression to use */
  info.bmiHeader.biSizeImage = biSizeImage; /* unsigned int   biSizeImage;       Size of image data */
  info.bmiHeader.biXPelsPerMeter = 2835;  /* int            biXPelsPerMeter;   X pixels per meter */
  info.bmiHeader.biYPelsPerMeter = 2835;  /* int            biYPelsPerMeter;   Y pixels per meter */
  info.bmiHeader.biClrUsed = 0; /* unsigned int   biClrUsed;         Number of colors used */
  info.bmiHeader.biClrImportant = 0;  /* unsigned int   biClrImportant;    Number of important colors */


  /* Write the file header, bitmap information, and bitmap pixel data... */
  write_word(stream, info.bmfHeader.bfType);
  write_dword(stream, info.bmfHeader.bfSize);
  write_word(stream, info.bmfHeader.bfReserved1);
  write_word(stream, info.bmfHeader.bfReserved2);
  write_dword(stream, info.bmfHeader.bfOffBits);
  write_dword(stream, info.bmiHeader.biSize);
  write_long(stream, info.bmiHeader.biWidth);
  write_long(stream, info.bmiHeader.biHeight);
  write_word(stream, info.bmiHeader.biPlanes);
  write_word(stream, info.bmiHeader.biBitCount);
  write_dword(stream, info.bmiHeader.biCompression);
  write_dword(stream, info.bmiHeader.biSizeImage);
  write_long(stream, info.bmiHeader.biXPelsPerMeter);
  write_long(stream, info.bmiHeader.biYPelsPerMeter);
  write_dword(stream, info.bmiHeader.biClrUsed);
  write_dword(stream, info.bmiHeader.biClrImportant);
  for (j = 0; j < M.Nj; j++)
  {
    for (i = 0; i < M.Ni; i++)
    {
      RealToRGBColor(M.VRe[i][j], RGB);
      fwrite(RGB, 1, 3, stream);
    }
    RGB[0] = 0;
    for (i = 0; i < ((4 - ((M.Ni * 3) & 3)) & 3); i++)
    {
      fwrite(RGB, 1, 1, stream);
    }
  }
  return 0;
}

/*  Print a Matrix to a file as a bitmap with collor coding acording to the function RealToRGBColor
 *
 */
int FilePrintRealMatrixToBitmap(MatrixType M, const char *const flname, void (*RealToRGBColor) (double, unsigned char *))
{
  FILE  *stream;
  uint32_t retval;
  if ((stream = fopen(flname, "wb")) == NULL)
  {
    printf("ERROR\nCan't open %s for writing with mode %s.\n", flname, "wb");
    return 1;
  }
  retval = StreamPrintRealMatrixToBitmap(M, stream, RealToRGBColor);
  fclose(stream);
  return retval;
}

/*  Print a Matrix to a stdout as a bitmap with collor coding acording to the function RealToRGBColor
 *
 */
uint32_t PrintRealMatrixToBitmap(MatrixType M, void (*RealToRGBColor) (double, unsigned char *))
{
  return StreamPrintRealMatrixToBitmap(M, stdout, RealToRGBColor);
}




/*  Print a list of 3 Matrices (RGB) to a stream as a bitmap with collor coding acording to the function RealToRGBColor
 *
 */
int StreamPrintReal3MatrixToBitmap(MatrixType * M, FILE * stream, void (*Real3ToRGBColor) (double *, unsigned char *))
{
  BitmapInfo info;              /* I - Bitmap information */
  uint32_t i, j;
  uint8_t RGB[3];
  double doubleRGB[3];
  unsigned int biSizeImage;

  biSizeImage = (unsigned int)((M[0].Ni * 3 + ((4 - ((M[0].Ni * 3) & 3)) & 3)) * M[0].Nj);

  info.bmfHeader.bfType = 0x4D42; /* unsigned short bfType;           Magic number for file */
  info.bmfHeader.bfSize = 54 + biSizeImage; /*  unsigned int   bfSize;          Size of file */
  info.bmfHeader.bfReserved1 = 0; /*  unsigned short bfReserved1;     Reserved */
  info.bmfHeader.bfReserved2 = 0; /*  unsigned short bfReserved2;     ... */
  info.bmfHeader.bfOffBits = 54;  /*  unsigned int   bfOffBits;       Offset to bitmap data */

  info.bmiHeader.biSize = 40;   /* unsigned int   biSize             Size of info header */
  info.bmiHeader.biWidth = (int)M[0].Ni;  /* int            biWidth;           Width of image */
  info.bmiHeader.biHeight = (int)M[0].Nj; /* int            biHeight;          Height of image */
  info.bmiHeader.biPlanes = 1;  /* unsigned short biPlanes;          Number of color planes */
  info.bmiHeader.biBitCount = 24; /* unsigned short biBitCount;        Number of bits per pixel */
  info.bmiHeader.biCompression = 0; /* unsigned int   biCompression;     Type of compression to use */
  info.bmiHeader.biSizeImage = biSizeImage; /* unsigned int   biSizeImage;       Size of image data */
  info.bmiHeader.biXPelsPerMeter = 2835;  /* int            biXPelsPerMeter;   X pixels per meter */
  info.bmiHeader.biYPelsPerMeter = 2835;  /* int            biYPelsPerMeter;   Y pixels per meter */
  info.bmiHeader.biClrUsed = 0; /* unsigned int   biClrUsed;         Number of colors used */
  info.bmiHeader.biClrImportant = 0;  /* unsigned int   biClrImportant;    Number of important colors */


  /* Write the file header, bitmap information, and bitmap pixel data... */
  write_word(stream, info.bmfHeader.bfType);
  write_dword(stream, info.bmfHeader.bfSize);
  write_word(stream, info.bmfHeader.bfReserved1);
  write_word(stream, info.bmfHeader.bfReserved2);
  write_dword(stream, info.bmfHeader.bfOffBits);
  write_dword(stream, info.bmiHeader.biSize);
  write_long(stream, info.bmiHeader.biWidth);
  write_long(stream, info.bmiHeader.biHeight);
  write_word(stream, info.bmiHeader.biPlanes);
  write_word(stream, info.bmiHeader.biBitCount);
  write_dword(stream, info.bmiHeader.biCompression);
  write_dword(stream, info.bmiHeader.biSizeImage);
  write_long(stream, info.bmiHeader.biXPelsPerMeter);
  write_long(stream, info.bmiHeader.biYPelsPerMeter);
  write_dword(stream, info.bmiHeader.biClrUsed);
  write_dword(stream, info.bmiHeader.biClrImportant);

  for (j = 0; j < M[0].Nj; j++)
  {
    for (i = 0; i < M[0].Ni; i++)
    {
      doubleRGB[0] = M[0].VRe[i][j];
      doubleRGB[1] = M[1].VRe[i][j];
      doubleRGB[2] = M[2].VRe[i][j];
      Real3ToRGBColor(doubleRGB, RGB);
      fwrite(RGB, 1, 3, stream);
    }
    RGB[0] = 0;
    for (i = 0; i < ((4 - ((M[0].Ni * 3) & 3)) & 3); i++)
    {
      fwrite(RGB, 1, 1, stream);
    }
  }

  return 0;
}

/*  Print a Matrix to a file as a bitmap with collor coding acording to the function RealToRGBColor
 *
 */
int FilePrintReal3MatrixToBitmap(MatrixType * M, const char *const flname, void (*Real3ToRGBColor) (double *, unsigned char *))
{
  FILE  *stream;
  uint32_t retval;
  if ((stream = fopen(flname, "wb")) == NULL)
  {
    printf("ERROR\nCan't open %s for writing with mode %s.\n", flname, "wb");
    return 1;
  }
  retval = StreamPrintReal3MatrixToBitmap(M, stream, Real3ToRGBColor);
  fclose(stream);
  return retval;
}

/*  Print a Matrix to a stdout as a bitmap with collor coding acording to the function RealToRGBColor
 *
 */
uint32_t PrintReal3MatrixToBitmap(MatrixType * M, void (*Real3ToRGBColor) (double *, unsigned char *))
{
  return StreamPrintReal3MatrixToBitmap(M, stdout, Real3ToRGBColor);
}


/*  Print a Matrix to a stream as a bitmap with collor coding acording to the function RealToRGBColor
 *  Rescale the matrix to get a bitmap of size NW by NH
 *  if NW (NH) larger than the size of the matrix use a weighted average of the point left and right in the matrix
 *  if NW (NH) smaller than the size of the matrix use a weighted average over several points in the original matrix
 */
int StreamPrintRealMatrixToBitmapRescale(MatrixType M, FILE * stream, void (*RealToRGBColor) (double, unsigned char *), uint32_t const NW, uint32_t const NH)
{
  //TODO:
  return 0;
}
int FilePrintRealMatrixToBitmapRescale(MatrixType M, const char *const flname, void (*RealToRGBColor) (double, unsigned char *), uint32_t const NW, uint32_t const NH)
{
  FILE  *stream;
  uint32_t retval;
  if ((stream = fopen(flname, "wb")) == NULL)
  {
    printf("ERROR\nCan't open %s for writing with mode %s.\n", flname, "wb");
    return 1;
  }
  retval = StreamPrintRealMatrixToBitmapRescale(M, stream, RealToRGBColor,NW,NH);
  fclose(stream);
  return retval;
}
uint32_t PrintRealMatrixToBitmapRescale(MatrixType M, void (*RealToRGBColor) (double, unsigned char *), uint32_t const NW, uint32_t const NH)
{
  return StreamPrintRealMatrixToBitmapRescale(M, stdout, RealToRGBColor,NW,NH);
}


/*  Print a list of 3 Matrices to a stream as a bitmap with collor coding acording to the function RealToRGBColor
 *  Rescale the matrix to get a bitmap of size NW by NH
 *  if NW (NH) larger than the size of the matrix use a weighted average of the point left and right in the matrix
 *  if NW (NH) smaller than the size of the matrix use a weighted average over several points in the original matrix
 */
int StreamPrintReal3MatrixToBitmapRescale(MatrixType * M, FILE * stream, void (*RealToRGBColor) (double, unsigned char *), uint32_t const NW, uint32_t const NH)
{
  //TODO:
  return 0;
}
int FilePrintReal3MatrixToBitmapRescale(MatrixType * M, const char *const flname, void (*RealToRGBColor) (double, unsigned char *), uint32_t const NW, uint32_t const NH)
{
  FILE  *stream;
  uint32_t retval;
  if ((stream = fopen(flname, "wb")) == NULL)
  {
    printf("ERROR\nCan't open %s for writing with mode %s.\n", flname, "wb");
    return 1;
  }
  retval = StreamPrintReal3MatrixToBitmapRescale(M, stream, RealToRGBColor,NW,NH);
  fclose(stream);
  return retval;
}
uint32_t PrintReal3MatrixToBitmapRescale(MatrixType * M, void (*RealToRGBColor) (double, unsigned char *), uint32_t const NW, uint32_t const NH)
{
  return StreamPrintReal3MatrixToBitmapRescale(M, stdout, RealToRGBColor,NW,NH);
}
