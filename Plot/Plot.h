#include "StructDef.cpp"
#include "GraphicsToPS.h"
#include "GraphicsToSVG.h"
#include "ScientificPlots.h"

void InitLine(LineType* Line);
int InitLine(LineType* Line, int NPoints);
LineType* NewLine();
void FreeLine(LineType* Line);
void InitCurve(CurveType* Curve);
CurveType* NewCurve();
void FreeCurve(CurveType* Curve);
void InitArc(ArcType* Arc);
ArcType* NewArc();
void FreeArc(ArcType* Arc);
void InitText(TextType* Text);
TextType* NewText();
void FreeText(TextType* Text);
void InitBitMap(BitMapType* BitMap);
BitMapType* NewBitMap();
void FreeBitMap(BitMapType* BitMap);
void InitGraphics(GraphicsType* G);
GraphicsType* NewGraphics();
void FreeGraphics(GraphicsType* Graphics);
int AddLineInherentPointMemory(GraphicsType* Graphics, uint8_t R, uint8_t G, uint8_t B, uint8_t alpha, double Width, int NPoints, PointType * Points, bool CloseAndFill, uint8_t RF, uint8_t GF, uint8_t BF, uint8_t alphaF );
int AddLine(GraphicsType* Graphics, uint8_t R, uint8_t G, uint8_t B, uint8_t alpha, double Width, int NPoints, PointType * Points, bool CloseAndFill, uint8_t RF, uint8_t GF, uint8_t BF, uint8_t alphaF );
int AddLineInherentLineMemory(GraphicsType* Graphics, LineType * Line);
int AddLine(GraphicsType* Graphics, LineType * Line);
int AddCurveInherentMemory(GraphicsType *Graphics, CurveType * Curve);
int AddCurve(GraphicsType *Graphics, CurveType * Curve);
int AddArcInherentMemory(GraphicsType *Graphics, ArcType * Arc);
int AddArc(GraphicsType *Graphics, ArcType * Arc);
int AddTextInherentMemory(GraphicsType* Graphics, TextType* Text);
int AddText(GraphicsType* Graphics, TextType* Text);
int AddBitMapInherentMemory(GraphicsType *Graphics, BitMapType * BitMap);
int AddBitMap(GraphicsType *Graphics, BitMapType * BitMap);





