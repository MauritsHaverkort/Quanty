typedef struct PointType
{
  double X;
  double Y;
} PointType;

// A colour is stored as 0xAARRGGBB: red, green and blue, and AA the transparency, 00 fully
// opaque and FF fully transparent. A colour written as 0xRRGGBB is therefore opaque.
// Lua gives colours in the same form, so they are stored as they are given.
static inline uint32_t ColorRed(uint32_t c)          { return (c >> 16) & 0xFF; }
static inline uint32_t ColorGreen(uint32_t c)        { return (c >>  8) & 0xFF; }
static inline uint32_t ColorBlue(uint32_t c)         { return  c        & 0xFF; }
static inline uint32_t ColorTransparency(uint32_t c) { return (c >> 24) & 0xFF; }
static inline double   ColorOpacity(uint32_t c)      { return (double)(255 - ColorTransparency(c)) / 255; }
static inline bool     ColorVisible(uint32_t c)      { return ColorTransparency(c) != 0xFF; }
// the colour c with the transparency T (0 to 255)
static inline uint32_t ColorWithTransparency(uint32_t c, uint32_t T) { return (c & 0x00FFFFFF) | ((T & 0xFF) << 24); }
// the colour c seen through a transparency T and with an opacity Alpha (both 0 to 255) on top of
// its own transparency: the opacities multiply
static inline uint32_t ColorCombine(uint32_t c, uint32_t T, uint32_t Alpha)
{
  double const Opacity = (255.0 - ColorTransparency(c)) * (255.0 - T) * Alpha / (255.0 * 255.0); // 0 to 255
  return ColorWithTransparency(c, 255 - (uint32_t)round(Opacity));
}


typedef struct LineType
{
  uint32_t Color;       // 0xAARRGGBB, AA the transparency
  double Width;         // Width of the line in fractions of CanvasSize sqrt(x*y)
  int NPoints;          // Allocation size of Point is NPoints
  PointType* Point;     // positions of the points in fractions of the Canvas
  bool CloseAndFill;    // filling
  uint32_t FillColor;   // 0xAARRGGBB. If fully transparent we do close but add no filling.
} LineType;

typedef struct CurveType
{
  uint32_t Color;       // 0xAARRGGBB, AA the transparency
  double Width;         // Width of the line in percentage of CanvasSize sqrt(x*y)
  int NPoints;          // Allocation size of Point is NPoints
  PointType* Point;     // positions of the points in fractions of the Canvas
  PointType* VL;        // Left  velocity (Bezier) of the points in fractions of the Canvas
  PointType* VR;        // Right velocity (Bezier) of the points in fractions of the Canvas
  bool CloseAndFill;    // filling
  uint32_t FillColor;   // 0xAARRGGBB. If fully transparent we do close but add no filling.
} CurveType;

typedef struct ArcType
{
  uint32_t Color;       // 0xAARRGGBB, AA the transparency
  double Width;         // Width of the line in percentage of CanvasSize
  double R;             // Radius in fractions of the Canvas sqrt(x*y)
  PointType Point;      // positions of the center in fractions of the Canvas
  double ThetaL;        // Starting theta -- arc is drawn clockwise x is horizontal, y vertical theta 0 is x running anti clockwise units radians
  double ThetaR;        // End theta of arc
  bool CloseAndFill;    // filling, via origin
  uint32_t FillColor;   // 0xAARRGGBB. If fully transparent we do close but add no filling.
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
#define SVG_PIXELS_PER_MM       28.3464566929133858267717

typedef struct TextType
{
  uint32_t Color;   // 0xAARRGGBB, AA the transparency
  PointType Point;  // Position of the string in percentage of CanvasSize
  double Angle;
  char Alignment;   // c (center) l (left) r (right)
  double Size;      // hight of the letter in percentage of CanvasSize.Y
  char* String;     // /0 terminated string with allocation size larger then its length
  const char* Font; // Pointer to list of possible fonts
} TextType;

typedef struct BitMapType
{
  PointType LB;              // left bottom position in fraction of the canvas
  PointType RT;              // right top position in fraction of the canvas
  int NpointsX;              // number of pixels in the x direction
  int NpointsY;              // number of pixels in the y direction
  uint32_t *Pixels;          // pixels 0xAARRGGBB, AA the transparency, NpointsX * NpointsY of them, row 0 at the top
  bool Pixelated;            // draw each pixel as a square, without smoothing between pixels where the format allows it
} BitMapType;



// The left botom coordinates of the canvas are ( CanvasOrigin.X, CanvasOrigin.Y )
// The right top coordinates of the canvas are ( CanvasOrigin.X + CanvasSize.X , CanvasOrigin.Y + CanvasSize.Y )
// X moves from left to right
// Y moves from bottom to top
// all units for the canvas are in mm
// when converted to PS or SVG the units are in 0.1 points a.k.a. 1/720 inch using integers.
// all layers use the same absolute coordinate system
// All elements on a layer are in relative coordinates with values between 0 and 1
// i.e. point (x,y) is found at absolute position ( CanvasOrigin.X + x CanvasSize.X, CanvasOrigin.Y + y CanvasSize.Y )
typedef struct GraphicsType
{
  PointType CanvasSize;    // Size of the Canvas in mm.
  PointType CanvasOrigin;  // Origin of the Canvas in mm -- mostly used for embedding
  int NLayer;              //Allocation size of Layer is NLayer
  GraphicsType** Layer;    // a list of NLayer pointers  -- the first pointer points to the Graphics self.
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

// if you make changes here you also need to make changes in LuaCheckOptions
typedef struct LinePlotOptions
{
  double   LineWidth;
  uint32_t LineColor;   // 0xAARRGGBB, AA the transparency
  double   PointSize;
  uint32_t PointColor;  // 0xAARRGGBB, AA the transparency
} LinePlotOptions;

typedef struct LinePlotOptionsArray
{
  LinePlotOptions *Line;
  int NAllocLine;
} LinePlotOptionsArray;

// if you make changes here you also need to make changes in LuaCheckOptions
typedef struct FramePlotOptions
{
  PointType PlotRange[2];
  PointType Margins[2];
  double   LineWidth;
  uint32_t LineColor;   // 0xAARRGGBB, AA the transparency
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
  double   XAxis;
  double   YAxis;
  bool     XTicks;   // draw the ticks and tick labels of the x axis
  bool     YTicks;   // draw the ticks and tick labels of the y axis
} FramePlotOptions;
