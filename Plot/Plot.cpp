
#include "BitMap.cpp"
#include "GraphicsToPS.cpp"
#include "GraphicsToSVG.cpp"
#include "ScientificPlots.cpp"



// Functions to initialize Graphics and Graphics elements (Line, Curve, Arc, Text, BitMap)

void InitLine(LineType* Line)
{
  Line->Color[0]=0;
  Line->Color[1]=0;
  Line->Color[2]=0;
  Line->Color[3]=0xFF;
  Line->Width=0.001;
  Line->NPoints=0;
  Line->Point=NULL;
  Line->CloseAndFill=false;
  Line->FillColor[0]=0;
  Line->FillColor[1]=0;
  Line->FillColor[2]=0;
  Line->FillColor[3]=0xFF;
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
  Curve->Color[0]=0;
  Curve->Color[1]=0;
  Curve->Color[2]=0;
  Curve->Color[3]=0xFF;
  Curve->Width=0.001;
  Curve->NPoints=0;
  Curve->Point=NULL;
  Curve->VL=NULL;
  Curve->VR=NULL;
  Curve->CloseAndFill=false;
  Curve->FillColor[0]=0;
  Curve->FillColor[1]=0;
  Curve->FillColor[2]=0;
  Curve->FillColor[3]=0xFF;
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
  Arc->Color[0]=0;
  Arc->Color[1]=0;
  Arc->Color[2]=0;
  Arc->Color[3]=0xFF;
  Arc->Width=0.001;
  Arc->R=1;
  Arc->Point.X=0;
  Arc->Point.Y=0;
  Arc->ThetaR=0;
  Arc->ThetaL=2*Pi;
  Arc->CloseAndFill=false;
  Arc->FillColor[0]=0;
  Arc->FillColor[1]=0;
  Arc->FillColor[2]=0;
  Arc->FillColor[3]=0xFF;
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
  Text->Color[0]=0;
  Text->Color[1]=0;
  Text->Color[2]=0;
  Text->Color[3]=0xFF;
  Text->Point.X=0;
  Text->Point.Y=0;
  Text->Angle=0;
  Text->Alignment='l';
  Text->String=NULL;
  Text->Font=TIMES_NEW_ROMAN;
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
  G->CanvasSize.X=297.3017787506803; //A4 landscape
  G->CanvasSize.Y=210.2241038134286;
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
  }
  if(Graphics->NLayer!=1) free(Graphics->Layer);
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
  Graphics->Line[Graphics->NLine].Color[0]     = R;
  Graphics->Line[Graphics->NLine].Color[1]     = G;
  Graphics->Line[Graphics->NLine].Color[2]     = B;
  Graphics->Line[Graphics->NLine].Color[3]     = alpha; //0xFF full collor 0x0 full transparant
  Graphics->Line[Graphics->NLine].Width        = Width;
  Graphics->Line[Graphics->NLine].NPoints      = NPoints;
  Graphics->Line[Graphics->NLine].Point        = Points;
  Graphics->Line[Graphics->NLine].CloseAndFill = CloseAndFill;
  Graphics->Line[Graphics->NLine].FillColor[0] = RF;
  Graphics->Line[Graphics->NLine].FillColor[1] = GF;
  Graphics->Line[Graphics->NLine].FillColor[2] = BF;
  Graphics->Line[Graphics->NLine].FillColor[3] = alphaF; //0xFF full collor 0x0 full transparant
  
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

