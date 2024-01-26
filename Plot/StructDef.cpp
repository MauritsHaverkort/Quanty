typedef struct PointType
{
  double X;
  double Y;
} PointType;


typedef struct LineType
{
  uint8_t Color[4];     // RGB alpha
  double Width;         // Width of the line in fractions of CanvasSize sqrt(x*y)
  int NPoints;          // Allocation size of Point is NPoints
  PointType* Point;     // positions of the points in fractions of the Canvas
  bool CloseAndFill;    // filling
  uint8_t FillColor[4]; // RGB alpha. If alpha = 0 we do close but add no filling.
} LineType;

typedef struct CurveType
{
  uint8_t Color[4];     // RGB alpha
  double Width;         // Width of the line in percentage of CanvasSize sqrt(x*y)
  int NPoints;          // Allocation size of Point is NPoints
  PointType* Point;     // positions of the points in fractions of the Canvas
  PointType* VL;        // Left  velocity (Bezier) of the points in fractions of the Canvas
  PointType* VR;        // Right velocity (Bezier) of the points in fractions of the Canvas
  bool CloseAndFill;    // filling
  uint8_t FillColor[4]; // RGB alpha. If alpha = 0 we do close but add no filling.
} CurveType;

typedef struct ArcType
{
  uint8_t Color[4];     // RGB alpha
  double Width;         // Width of the line in percentage of CanvasSize
  double R;             // Radius in fractions of the Canvas sqrt(x*y)
  PointType Point;      // positions of the center in fractions of the Canvas
  double ThetaL;        // Starting theta -- arc is drawn clockwise x is horizontal, y vertical theta 0 is x running anti clockwise units radians
  double ThetaR;        // End theta of arc
  bool CloseAndFill;    // filling, via origin
  uint8_t FillColor[4]; // RGB alpha.  If alpha = 0 we do close but add no filling.
} ArcType;

// for fonts we depend on the fonts installed on the users computer
// setting the font to a user defined string is always possible, but might fail
// here a list of fonts generally present on modern machines that I would sugest to choose from
// boaring but safe
const char* AvailableFontTypes[] = { "ArialMT",
                                     "Arial-Black",
                                     "CourierNewPSMT",
                                     "Georgia",
                                     "Impact",
                                     "Tahoma",
                                     "TimesNewRomanPSMT",
                                     "TrebuchetMS",
                                     "Verdana",
                                     "Symbol",
                                     "MicrosoftSansSerif"};
#define ARIAL AvailableFontTypes[0]
#define ARIAL_BLACK AvailableFontTypes[1]
#define COURIER_NEW AvailableFontTypes[2]
#define GEORGIA AvailableFontTypes[3]
#define IMPACT AvailableFontTypes[4]
#define TAHOMA AvailableFontTypes[5]
#define TIMES_NEW_ROMAN AvailableFontTypes[6]
#define TREBUDHET_MS AvailableFontTypes[7]
#define VERDANA AvailableFontTypes[8]
#define SYMBOL AvailableFontTypes[9]
#define MS_SANS_SERIF AvailableFontTypes[10]

#define NUMBER_OF_GLOBAL_FONTS 11
// Postscript standardly has 72 pixels per inch or 2.83464566929133858267717 per mm
// This is not enough, espeically if one zooms.
// The post-script format does not care about the units, only the interpretators do
// and these can normaly scale a page. We thus can freely redifine the number of pixels per mm
// To keep numbers recognizable we can scale by a factor of 10.
#define POSTSCIPT_PIXELS_PER_MM 28.3464566929133858267717

typedef struct TextType
{
  uint8_t Color[4]; //RGB alpha
  PointType Point;  // Position of the string in percentage of CanvasSize
  double Angle;
  char Alignment;   // c (center) l (left) r (right)
  double Size;      // hight of the letter in percentage of CanvasSize.Y
  char* String;     // /0 terminated string with allocation size larger then its length
  const char* Font; // Pointer to list of possible fonts
} TextType;

typedef struct BitMapType
{
  PointType LB;              // left top position in fraction of the canvas
  PointType RT;              // right bottom position in fraction of the canvas
  int NpointsX;              // number of pixels in the x direction
  int NpointsY;              // number of pixels in the y direction
  uint32_t *Pixels;          // pixels in format RGBalpha of size NpointsX * NpointsY
  uint32_t TransparantColor; // if (TransparantColor & 1) == 1 the first 24 bit contain the RGB color that is transparant
} BitMapType;


typedef struct GraphicsType
{
  PointType CanvasSize;    // Size of the Canvas in mm.
  PointType CanvasOrigin;  // Origin of the Canvas in mm -- mostly used for embedding
  int NLayer;              //Allocation size of Layer is NLayer
  GraphicsType** Layer;
  int NLine;               //Allocation size of Line   is 2^n with n the smallest integer such that the number of actual layers present fit
  LineType* Line;
  int NCurve;              //Allocation size of Curve  is 2^n with n the smallest integer such that the number of actual layers present fit
  CurveType* Curve;
  int NArc;                //Allocation size of Arc     is 2^n with n the smallest integer such that the number of actual layers present fit
  ArcType* Arc;
  int NText;               //Allocation size of Text   is 2^n with n the smallest integer such that the number of actual layers present fit
  TextType* Text;
  int NBitMap;             //Allocation size of BitMap is 2^n with n the smallest integer such that the number of actual layers present fit
  BitMapType* BitMap;
} GraphicsType;

typedef struct LinePlotOptions
{
  double   LineWidth;
  uint8_t  LineColor[4];
  double   PointSize;
  uint8_t  PointColor[4];
} LinePlotOptions;

typedef struct FramePlotOptions
{
  PointType PlotRange[2];
  PointType Margins[2];
  double   LineWidth;
  uint8_t  LineColor[4];
  double   FontSize;
  const char* Font;
  char*    XLabel;
  char*    YLabel;
  double   FirstXTick;
  double   dXTick;
  int      NMinorXTick;
  double   FirstYTick;
  double   dYTick;
  int      NMinorYTick;
  char*    PlotXScaleFormat;
  char*    PlotYScaleFormat;
} FramePlotOptions;
