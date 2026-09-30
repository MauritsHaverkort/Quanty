
#include "BitMap.cpp"
#include "GraphicsToPS.cpp"
#include "GraphicsToSVG.cpp"
#include "ScientificPlots.cpp"



// Functions to initialize Graphics and Graphics elements (Line, Curve, Arc, Text, BitMap)

void InitLine(LineType* Line)
{
  Line->Color=0x00000000; // opaque black
  Line->Width=0.001;
  Line->NPoints=0;
  Line->Point=NULL;
  Line->CloseAndFill=false;
  Line->FillColor=0x00000000; // opaque black
}

int InitLine(LineType* Line, int NPoints)
{
  InitLine(Line);
  Line->NPoints=NPoints;
  Line->Point=(PointType*)malloc(NPoints*sizeof(PointType));
  if(Line->Point==NULL)
  {
    printf("ERROR: malloc failed in InitLine\n");
    return 1;
  }
  return 0;
}

LineType* NewLine()
{
  LineType *ret = (LineType*)malloc(sizeof(*ret));
  InitLine(ret);
  return ret;
}

LineType* NewLine(int NPoints)
{
  LineType *ret = (LineType*)malloc(sizeof(*ret));
  InitLine(ret, NPoints);
  return ret;
}

void FreeLine(LineType* Line)
{
  Line->NPoints=0;
  free(Line->Point);
}

void InitCurve(CurveType* Curve)
{
  Curve->Color=0x00000000; // opaque black
  Curve->Width=0.001;
  Curve->NPoints=0;
  Curve->Point=NULL;
  Curve->VL=NULL;
  Curve->VR=NULL;
  Curve->CloseAndFill=false;
  Curve->FillColor=0x00000000; // opaque black
}

int InitCurve(CurveType* Curve, int NPoints)
{
  InitCurve(Curve);
  Curve->NPoints=NPoints;
  Curve->Point=(PointType*)malloc(NPoints*sizeof(PointType));
  Curve->VL=(PointType*)malloc(NPoints*sizeof(PointType));
  Curve->VR=(PointType*)malloc(NPoints*sizeof(PointType));
  if((Curve->Point==NULL)||(Curve->VL==NULL)||(Curve->VR==NULL))
  {
    printf("ERROR: malloc failed in InitCurve\n");
    return 1;
  }
  return 0;
}

CurveType* NewCurve()
{
  CurveType *ret = (CurveType*)malloc(sizeof(*ret));
  InitCurve(ret);
  return ret;
}

CurveType* NewCurve(int NPoints)
{
  CurveType *ret = (CurveType*)malloc(sizeof(*ret));
  InitCurve(ret, NPoints);
  return ret;
}

void FreeCurve(CurveType* Curve)
{
  Curve->NPoints=0;
  free(Curve->Point);
  free(Curve->VL);
  free(Curve->VR);
}

void InitArc(ArcType* Arc)
{
  Arc->Color=0x00000000; // opaque black
  Arc->Width=0.001;
  Arc->R=1;
  Arc->Point.X=0;
  Arc->Point.Y=0;
  Arc->ThetaR=0;
  Arc->ThetaL=2*Pi;
  Arc->CloseAndFill=false;
  Arc->FillColor=0x00000000; // opaque black
}

ArcType* NewArc()
{
  ArcType *ret = (ArcType*)malloc(sizeof(*ret));
  InitArc(ret);
  return ret;
}

void FreeArc(ArcType* Arc)
{
  // nothing to do here
}

void InitText(TextType* Text)
{
  Text->Color=0x00000000; // opaque black
  Text->Point.X=0;
  Text->Point.Y=0;
  Text->Angle=0;
  Text->Alignment='l';
  Text->Size=0.05;      // the default of Graphics.Add
  Text->String=NULL;
  Text->Font=TIMES_NEW_ROMAN;
}

/*  A text keeps a pointer to the name of its font and never frees it, and copies of a text share
 *  the pointer. The name therefore has to live as long as the program. A name found in
 *  AvailableFontTypes is returned as the pointer into that list; any other name is copied once
 *  into FontNames, which is never freed, and the same copy is returned for the same name. Names
 *  given in lua (Graphics.Add, the Font option of the plot routines) point into lua strings that
 *  the garbage collector frees, and must go through here before a text keeps them.
 */
const char* InternFont(const char* Font)
{
  static char** FontNames=NULL;
  static int NFontNames=0;
  char** Temp;
  int i;
  if(Font==NULL) return NULL;
  for(i=0;i<NUMBER_OF_GLOBAL_FONTS;i++)
  {
    if((Font==AvailableFontTypes[i])||(strcmp(Font, AvailableFontTypes[i])==0)) return AvailableFontTypes[i];
  }
  for(i=0;i<NFontNames;i++)
  {
    if(strcmp(Font, FontNames[i])==0) return FontNames[i];
  }
  Temp = (char**)realloc(FontNames, (NFontNames+1)*sizeof(char*));
  if(Temp==NULL)
  {
    printf("Error in InternFont, malloc failed\n");
    return TIMES_NEW_ROMAN;
  }
  FontNames = Temp;
  FontNames[NFontNames] = (char*)malloc((strlen(Font)+1)*sizeof(char));
  if(FontNames[NFontNames]==NULL)
  {
    printf("Error in InternFont, malloc failed\n");
    return TIMES_NEW_ROMAN;
  }
  memcpy(FontNames[NFontNames], Font, (strlen(Font)+1)*sizeof(char));
  NFontNames++;
  return FontNames[NFontNames-1];
}

TextType* NewText()
{
  TextType *ret = (TextType*)malloc(sizeof(*ret));
  InitText(ret);
  return ret;
}

void FreeText(TextType* Text)
{
  free(Text->String);
}

void InitBitMap(BitMapType* BitMap)
{
  BitMap->LB.X=0;
  BitMap->LB.Y=0;
  BitMap->RT.X=1;
  BitMap->RT.Y=1;
  BitMap->NpointsX=0;
  BitMap->NpointsY=0;
  BitMap->Pixels=NULL;
  BitMap->Pixelated=false;
  BitMap->Angle=0;
}

BitMapType* NewBitMap()
{
  BitMapType *ret = (BitMapType*)malloc(sizeof(*ret));
  InitBitMap(ret);
  return ret;
}

void FreeBitMap(BitMapType* BitMap)
{
  free(BitMap->Pixels);
}

void InitGraphics(GraphicsType* G)
{
  G->CanvasOrigin.X=0;
  G->CanvasOrigin.Y=0;
  G->CanvasSize.X=297.30177875068026667937499264012; //A4 landscape
  G->CanvasSize.Y=210.22410381342863575778136905830;
  G->NLayer  = 1;
  G->Layer = (GraphicsType**)malloc(sizeof(GraphicsType*));
  G->Layer[0]=G;
  G->NLine = 0;
  G->Line=NULL;
  G->NCurve = 0;
  G->Curve=NULL;
  G->NArc = 0;
  G->Arc=NULL;
  G->NText=0;
  G->Text=NULL;
  G->NBitMap=0;
  G->BitMap=NULL;
}

GraphicsType* NewGraphics()
{
  GraphicsType *ret = (GraphicsType*)malloc(sizeof(*ret));
  InitGraphics(ret);
  return ret;
}

void FreeGraphics(GraphicsType* Graphics)
{
  int i;
  for(i=1;i<Graphics->NLayer;i++)
  {
    FreeGraphics( Graphics->Layer[i] );
    free(Graphics->Layer[i]);
  }
  free(Graphics->Layer);
  for(i=0;i<Graphics->NLine;i++) FreeLine(&(Graphics->Line[i]));
  free(Graphics->Line);
  for(i=0;i<Graphics->NCurve;i++) FreeCurve(&(Graphics->Curve[i]));
  free(Graphics->Curve);
  for(i=0;i<Graphics->NArc;i++) FreeArc(&(Graphics->Arc[i]));
  free(Graphics->Arc);
  for(i=0;i<Graphics->NText;i++) FreeText(&(Graphics->Text[i]));
  free(Graphics->Text);
  for(i=0;i<Graphics->NBitMap;i++) FreeBitMap(&(Graphics->BitMap[i]));
  free(Graphics->BitMap);
}

// Functions to add Graphics elememts to a Graph

int AddLineInherentPointMemory(GraphicsType* Graphics, uint8_t R, uint8_t G, uint8_t B, uint8_t alpha, double Width, int NPoints, PointType * Points, bool CloseAndFill, uint8_t RF, uint8_t GF, uint8_t BF, uint8_t alphaF )
{
  LineType* Temp = Graphics->Line;
  // test if NLine == 2^n
  if( Graphics->NLine==0 )
  {
    Temp = (LineType*)malloc(sizeof(LineType) );
  }
  else if( ( Graphics->NLine & ( Graphics->NLine-1 ) ) == 0 )
  {
    Temp = (LineType*)realloc(Graphics->Line, 2*Graphics->NLine *sizeof(LineType) );
  }
  if(Temp==NULL)
  {
    printf("Error in AddLine, malloc failed\n");
    return 1;
  }
  Graphics->Line = Temp;
  // alpha and alphaF are opacities, 0xFF full colour and 0x0 fully transparent
  Graphics->Line[Graphics->NLine].Color        = ((uint32_t)(255-alpha)<<24) | ((uint32_t)R<<16) | ((uint32_t)G<<8) | (uint32_t)B;
  Graphics->Line[Graphics->NLine].Width        = Width;
  Graphics->Line[Graphics->NLine].NPoints      = NPoints;
  Graphics->Line[Graphics->NLine].Point        = Points;
  Graphics->Line[Graphics->NLine].CloseAndFill = CloseAndFill;
  Graphics->Line[Graphics->NLine].FillColor    = ((uint32_t)(255-alphaF)<<24) | ((uint32_t)RF<<16) | ((uint32_t)GF<<8) | (uint32_t)BF;
  
  Graphics->NLine = Graphics->NLine + 1;
  return 0;
}

int AddLine(GraphicsType* Graphics, uint8_t R, uint8_t G, uint8_t B, uint8_t alpha, double Width, int NPoints, PointType * Points, bool CloseAndFill, uint8_t RF, uint8_t GF, uint8_t BF, uint8_t alphaF )
{
  if(AddLineInherentPointMemory(Graphics, R, G, B, alpha, Width, NPoints, Points, CloseAndFill, RF, GF, BF, alphaF)) return 1;
  Graphics->Line[Graphics->NLine-1].Point    = (PointType*) malloc(NPoints * sizeof(PointType));
  if(Graphics->Line[Graphics->NLine-1].Point==NULL)
  {
    printf("malloc failed in AddLine");
    return 1;
  }
  memcpy(Graphics->Line[Graphics->NLine-1].Point, Points, NPoints * sizeof(PointType));
  return 0;
}

int AddLineInherentLineMemory(GraphicsType* Graphics, LineType * Line)
{
  LineType* Temp = Graphics->Line;
  // test if NLine == 2^n
  if( Graphics->NLine==0 )
  {
    Temp = (LineType*)malloc(sizeof(LineType) );
  }
  else if( ( Graphics->NLine & ( Graphics->NLine-1 ) ) == 0 )
  {
    Temp = (LineType*)realloc(Graphics->Line, 2*Graphics->NLine *sizeof(LineType) );
  }
  if(Temp==NULL)
  {
    printf("Error in AddLine, malloc failed\n");
    return 1;
  }
  Graphics->Line = Temp;
  Graphics->Line[Graphics->NLine] = *Line;
  Graphics->NLine = Graphics->NLine + 1;
  return 0;
}

int AddLine(GraphicsType* Graphics, LineType * Line)
{
  if(AddLineInherentLineMemory(Graphics, Line)) return 1;
  Graphics->Line[Graphics->NLine-1].Point    = (PointType*) malloc(Line->NPoints * sizeof(PointType));
  if(Graphics->Line[Graphics->NLine-1].Point==NULL)
  {
    printf("malloc failed in AddLine");
    return 1;
  }
  memcpy(Graphics->Line[Graphics->NLine-1].Point, Line->Point, Line->NPoints * sizeof(PointType));  return 0;
}

int AddCurveInherentMemory(GraphicsType *Graphics, CurveType * Curve)
{
  CurveType* Temp = Graphics->Curve;
  // test if NLine == 2^n
  if( Graphics->NCurve==0 )
  {
    Temp = (CurveType*)malloc(sizeof(CurveType) );
  }
  else if( ( Graphics->NCurve & ( Graphics->NCurve-1 ) ) == 0 )
  {
    Temp = (CurveType*)realloc(Graphics->Curve, 2*Graphics->NCurve *sizeof(CurveType) );
  }
  if(Temp==NULL)
  {
    printf("Error in AddCurve, malloc failed\n");
    return 1;
  }
  Graphics->Curve = Temp;
  Graphics->Curve[Graphics->NCurve] = *Curve;
  Graphics->NCurve = Graphics->NCurve + 1;
  return 0;
}

int AddCurve(GraphicsType *Graphics, CurveType * Curve)
{
  if(AddCurveInherentMemory(Graphics, Curve)) return 1;
  Graphics->Curve[Graphics->NCurve-1].Point    = (PointType*) malloc(Curve->NPoints * sizeof(PointType));
  Graphics->Curve[Graphics->NCurve-1].VL       = (PointType*) malloc(Curve->NPoints * sizeof(PointType));
  Graphics->Curve[Graphics->NCurve-1].VR       = (PointType*) malloc(Curve->NPoints * sizeof(PointType));
  if((Graphics->Curve[Graphics->NCurve-1].Point==NULL) || (Graphics->Curve[Graphics->NCurve-1].VL==NULL) || (Graphics->Curve[Graphics->NCurve-1].VR==NULL) )
  {
    printf("malloc failed in AddCurve");
    return 1;
  }
  memcpy(Graphics->Curve[Graphics->NCurve-1].VR, Curve->VR, Curve->NPoints * sizeof(PointType));
  memcpy(Graphics->Curve[Graphics->NCurve-1].VL, Curve->VL, Curve->NPoints * sizeof(PointType));
  memcpy(Graphics->Curve[Graphics->NCurve-1].Point, Curve->Point, Curve->NPoints * sizeof(PointType));
  return 0;
}

int AddArcInherentMemory(GraphicsType *Graphics, ArcType * Arc)
{
  ArcType* Temp = Graphics->Arc;
  // test if NLine == 2^n
  if( Graphics->NArc==0 )
  {
    Temp = (ArcType*)malloc(sizeof(ArcType) );
  }
  else if( ( Graphics->NArc & ( Graphics->NArc-1 ) ) == 0 )
  {
    Temp = (ArcType*)realloc(Graphics->Arc, 2*Graphics->NArc *sizeof(ArcType) );
  }
  if(Temp==NULL)
  {
    printf("Error in AddArc, malloc failed\n");
    return 1;
  }
  Graphics->Arc = Temp;
  Graphics->Arc[Graphics->NArc] = *Arc;
  Graphics->NArc = Graphics->NArc + 1;
  return 0;
}

int AddArc(GraphicsType *Graphics, ArcType * Arc)
{
  return AddArcInherentMemory(Graphics, Arc);
}

int AddTextInherentMemory(GraphicsType* Graphics, TextType* Text)
{
  TextType* Temp=Graphics->Text;
  // test if NLine == 2^n
  if( Graphics->NText==0 )
  {
    Temp = (TextType*)malloc(sizeof(TextType) );
  }
  else if( ( Graphics->NText & ( Graphics->NText-1 ) ) == 0 )
  {
    Temp = (TextType*)realloc(Graphics->Text, 2*Graphics->NText *sizeof(TextType) );
  }
  if(Temp==NULL)
  {
    printf("Error in AddText, malloc failed\n");
    return 1;
  }
  Graphics->Text = Temp;
  Graphics->Text[Graphics->NText] = *Text;
  // the font name may be a lua string that is freed later
  Graphics->Text[Graphics->NText].Font = InternFont(Text->Font);
  Graphics->NText = Graphics->NText + 1;
  return 0;
}

int AddText(GraphicsType* Graphics, TextType* Text)
{
  if (AddTextInherentMemory(Graphics,Text)!=0) return 1;
  Graphics->Text[Graphics->NText-1].String = (char*)malloc(( strlen(Text->String)+1 ) *sizeof(char));
  memcpy(Graphics->Text[Graphics->NText-1].String, Text->String, (( strlen(Text->String)+1 ) *sizeof(char)));
  return 0;
}

int AddBitMapInherentMemory(GraphicsType *Graphics, BitMapType * BitMap)
{
  BitMapType* Temp = Graphics->BitMap;
  // test if NLine == 2^n
  if( Graphics->NBitMap==0 )
  {
    Temp = (BitMapType*)malloc(sizeof(BitMapType) );
  }
  else if( ( Graphics->NBitMap & ( Graphics->NBitMap-1 ) ) == 0 )
  {
    Temp = (BitMapType*)realloc(Graphics->BitMap, 2*Graphics->NBitMap *sizeof(BitMapType) );
  }
  if(Temp==NULL)
  {
    printf("Error in AddBitMap, malloc failed\n");
    return 1;
  }
  Graphics->BitMap = Temp;
  Graphics->BitMap[Graphics->NBitMap] = *BitMap;
  Graphics->NBitMap = Graphics->NBitMap + 1;
  return 0;
}

int AddBitMap(GraphicsType *Graphics, BitMapType * BitMap)
{
  if (AddBitMapInherentMemory(Graphics,BitMap)!=0) return 1;
  Graphics->BitMap[Graphics->NBitMap-1].Pixels = (uint32_t*)malloc( BitMap->NpointsX * BitMap->NpointsY * sizeof(uint32_t) );
  memcpy(Graphics->BitMap[Graphics->NBitMap-1].Pixels, BitMap->Pixels,  BitMap->NpointsX * BitMap->NpointsY * sizeof(uint32_t) );
  return 0;
}


int CopyGraphics(GraphicsType *Graphics, GraphicsType *Copy)
{
  int i;
  InitGraphics(Copy);
  Copy->CanvasOrigin.X=Graphics->CanvasOrigin.X;
  Copy->CanvasOrigin.Y=Graphics->CanvasOrigin.Y;
  Copy->CanvasSize.X=Graphics->CanvasSize.X;
  Copy->CanvasSize.Y=Graphics->CanvasSize.Y;
  
  for(i=0;i<Graphics->NLine;i++)
  {
    if( AddLine(Copy, &(Graphics->Line[i])) ) return 1;
  }
  for(i=0;i<Graphics->NCurve;i++)
  {
    if( AddCurve(Copy, &(Graphics->Curve[i])) ) return 1;
  }
  for(i=0;i<Graphics->NArc;i++)
  {
    if( AddArc(Copy, &(Graphics->Arc[i])) ) return 1;
  }
  for(i=0;i<Graphics->NText;i++)
  {
    if( AddText(Copy, &(Graphics->Text[i])) ) return 1;
  }
  for(i=0;i<Graphics->NBitMap;i++)
  {
    if( AddBitMap(Copy, &(Graphics->BitMap[i])) ) return 1;
  }
  for(i=1;i<Graphics->NLayer;i++)
  {
    if( AddLayer(Copy, Graphics->Layer[i]) ) return 1;
  }
  return 0;
}


int AddLayerInherentMemory(GraphicsType *Graphics, GraphicsType *Layer)
{
  // Graphics->Layer is an array of pointers to type GraphicsType.
  // The first element points to itselve.
  // InitGraphics allocated the array
  GraphicsType** Temp = Graphics->Layer;
  // test if NLayer == 2^n
  if( ( Graphics->NLayer & ( Graphics->NLayer-1 ) ) == 0 )
  {
    Temp = (GraphicsType**)realloc(Graphics->Layer, 2*Graphics->NLayer *sizeof(GraphicsType*) );
  }
  if(Temp==NULL)
  {
    printf("Error in AddLayerInherentMemory, malloc failed\n");
    return 1;
  }
  Graphics->Layer = Temp;
  Graphics->Layer[Graphics->NLayer] = Layer;
  Graphics->NLayer = Graphics->NLayer + 1;
  return 0;
}


int AddLayer(GraphicsType *Graphics, GraphicsType *Layer)
{
  GraphicsType *Copy;
  Copy = (GraphicsType*)malloc(sizeof(GraphicsType));
  if( CopyGraphics(Layer, Copy) ) return 1;
  if( AddLayerInherentMemory(Graphics, Copy) ) return 1;
  return 0;
}


/*  Rotation of a graphics with all its layers by Angle (radians, counterclockwise) about the point
 *  Center (mm, in the absolute coordinates of the page). Everything keeps its size and shape in mm.
 *  The canvas of each layer becomes the bounding box of its rotated canvas, and all positions are
 *  given anew as fractions of it; widths, sizes of text and radii, fractions of sqrt(width*height),
 *  are scaled such that they keep their size in mm.
 *  A bitmap is turned by a multiple of 90 degrees by moving its pixels; by any other angle it keeps
 *  its pixels and gets the angle in BitMap.Angle, about the centre of its rectangle.
 */
typedef struct RotationType
{
  double Angle;      // radians, counterclockwise
  int    Quarter;    // Angle is Quarter times 90 degrees, 0 to 3; -1 for any other angle
  double c, s;       // cos and sin of Angle, exact for multiples of 90 degrees
  PointType Center;  // mm
} RotationType;

// an angle reduced to the interval -pi to pi
static double ReduceAngle(double a)
{
  return a - 2*Pi*floor((a+Pi)/(2*Pi));
}

// the point P (mm) rotated about the centre
static PointType RotatePointMM(RotationType const *R, PointType P)
{
  PointType Q;
  double x = P.X - R->Center.X, y = P.Y - R->Center.Y;
  Q.X = R->Center.X + R->c * x - R->s * y;
  Q.Y = R->Center.Y + R->s * x + R->c * y;
  return Q;
}

// a point in fractions of the canvas O, W, H rotated, in fractions of the canvas On, Wn, Hn
static void RotateFraction(RotationType const *R, PointType *P, PointType O, PointType S, PointType On, PointType Sn)
{
  PointType mm;
  mm.X = O.X + P->X * S.X;
  mm.Y = O.Y + P->Y * S.Y;
  mm = RotatePointMM(R, mm);
  P->X = (mm.X - On.X) / Sn.X;
  P->Y = (mm.Y - On.Y) / Sn.Y;
}

// the pixels of a bitmap turned by Quarter times 90 degrees counterclockwise; row 0 is at the top
static int RotateBitMapPixels(BitMapType *B, int Quarter)
{
  int NX = B->NpointsX, NY = B->NpointsY;
  int NXn = (Quarter%2) ? NY : NX, NYn = (Quarter%2) ? NX : NY;
  if(Quarter==0) return 0;
  uint32_t *P = (uint32_t*)malloc((size_t)NX * (size_t)NY * sizeof(uint32_t));
  if((P==NULL)&&(NX*NY!=0)) { printf("Error in RotateGraphics, malloc failed\n"); return 1; }
  for(int r=0;r<NYn;r++)
  {
    for(int c=0;c<NXn;c++)
    {
      int co, ro; // the pixel of the old bitmap that goes to column c and row r
      if(Quarter==1)      { co = NX-1-r; ro = c;      }
      else if(Quarter==2) { co = NX-1-c; ro = NY-1-r; }
      else                { co = r;      ro = NY-1-c; }
      P[r*NXn+c] = B->Pixels[ro*NX+co];
    }
  }
  free(B->Pixels);
  B->Pixels = P;
  B->NpointsX = NXn;
  B->NpointsY = NYn;
  return 0;
}

// are all canvases in G of a size larger than 0, such that G can be rotated
static bool Rotatable(GraphicsType const *G)
{
  if(!((G->CanvasSize.X>0)&&(G->CanvasSize.Y>0))) return false;
  for(int i=1;i<G->NLayer;i++) if(!Rotatable(G->Layer[i])) return false;
  return true;
}

static int RotateLayer(GraphicsType *G, RotationType const *R)
{
  int i, j;
  PointType O = G->CanvasOrigin, S = G->CanvasSize, On, Sn;
  // the bounding box of the rotated canvas; for a multiple of 90 degrees its size is exact
  PointType Corner[4] = {{O.X, O.Y}, {O.X+S.X, O.Y}, {O.X, O.Y+S.Y}, {O.X+S.X, O.Y+S.Y}};
  PointType Upper;
  for(i=0;i<4;i++)
  {
    PointType Q = RotatePointMM(R, Corner[i]);
    if((i==0)||(Q.X<On.X))    On.X    = Q.X;
    if((i==0)||(Q.Y<On.Y))    On.Y    = Q.Y;
    if((i==0)||(Q.X>Upper.X)) Upper.X = Q.X;
    if((i==0)||(Q.Y>Upper.Y)) Upper.Y = Q.Y;
  }
  if(R->Quarter>=0)
  {
    Sn.X = (R->Quarter%2) ? S.Y : S.X;
    Sn.Y = (R->Quarter%2) ? S.X : S.Y;
  }
  else
  {
    Sn.X = Upper.X - On.X;
    Sn.Y = Upper.Y - On.Y;
  }
  double const f = (R->Quarter>=0) ? 1.0 : sqrt((S.X*S.Y)/(Sn.X*Sn.Y));
  for(i=0;i<G->NLine;i++)
  {
    for(j=0;j<G->Line[i].NPoints;j++) RotateFraction(R, &(G->Line[i].Point[j]), O, S, On, Sn);
    G->Line[i].Width *= f;
  }
  for(i=0;i<G->NCurve;i++)
  {
    for(j=0;j<G->Curve[i].NPoints;j++)
    {
      RotateFraction(R, &(G->Curve[i].Point[j]), O, S, On, Sn);
      RotateFraction(R, &(G->Curve[i].VL[j]),    O, S, On, Sn);
      RotateFraction(R, &(G->Curve[i].VR[j]),    O, S, On, Sn);
    }
    G->Curve[i].Width *= f;
  }
  for(i=0;i<G->NArc;i++)
  {
    ArcType *A = &(G->Arc[i]);
    RotateFraction(R, &(A->Point), O, S, On, Sn);
    A->R     *= f;
    A->Width *= f;
    // both ends turn by Angle and are reduced by the same whole turns, which keeps the arc
    A->ThetaR += R->Angle;
    A->ThetaL += R->Angle;
    double const Turns = 2*Pi*floor((A->ThetaR+Pi)/(2*Pi));
    A->ThetaR -= Turns;
    A->ThetaL -= Turns;
  }
  for(i=0;i<G->NText;i++)
  {
    RotateFraction(R, &(G->Text[i].Point), O, S, On, Sn);
    G->Text[i].Size *= f;
    G->Text[i].Angle = ReduceAngle(G->Text[i].Angle + R->Angle);
  }
  for(i=0;i<G->NBitMap;i++)
  {
    BitMapType *B = &(G->BitMap[i]);
    if((R->Quarter>=0)&&(B->Angle==0))
    {
      // the rotated rectangle is again a rectangle along the axes
      PointType LB = B->LB, RT = B->RT;
      RotateFraction(R, &LB, O, S, On, Sn);
      RotateFraction(R, &RT, O, S, On, Sn);
      B->LB.X = Min(LB.X, RT.X);  B->RT.X = Max(LB.X, RT.X);
      B->LB.Y = Min(LB.Y, RT.Y);  B->RT.Y = Max(LB.Y, RT.Y);
      if(RotateBitMapPixels(B, R->Quarter)) return 1;
    }
    else
    {
      // the centre moves, the rectangle keeps its size in mm and turns about the centre
      PointType M = {(B->LB.X+B->RT.X)/2, (B->LB.Y+B->RT.Y)/2};
      double const hx = (B->RT.X-B->LB.X)*S.X/2, hy = (B->RT.Y-B->LB.Y)*S.Y/2;
      RotateFraction(R, &M, O, S, On, Sn);
      B->LB.X = M.X - hx/Sn.X;  B->RT.X = M.X + hx/Sn.X;
      B->LB.Y = M.Y - hy/Sn.Y;  B->RT.Y = M.Y + hy/Sn.Y;
      B->Angle = ReduceAngle(B->Angle + R->Angle);
    }
  }
  G->CanvasOrigin = On;
  G->CanvasSize   = Sn;
  for(i=1;i<G->NLayer;i++)
  {
    if(RotateLayer(G->Layer[i], R)) return 1;
  }
  return 0;
}

int RotateGraphics(GraphicsType *G, double Angle, PointType Center)
{
  RotationType R;
  if(!Rotatable(G)) return 1; // checked before anything changes
  R.Angle  = Angle;
  R.Center = Center;
  // a multiple of 90 degrees gets exact values of cos and sin, so that four turns give the same graphics
  double const q = Angle / (Pi/2);
  if(Abs(q-round(q)) < 1E-12)
  {
    double m = fmod(round(q), 4.0); // exact, also for large q, where a cast to int is undefined
    if(m<0) m += 4;
    R.Quarter = (int)m;
    double const cq[4] = {1, 0, -1, 0}, sq[4] = {0, 1, 0, -1};
    R.c = cq[R.Quarter];
    R.s = sq[R.Quarter];
  }
  else
  {
    R.Quarter = -1;
    R.c = cos(Angle);
    R.s = sin(Angle);
  }
  return RotateLayer(G, &R);
}
