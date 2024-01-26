
// Include the header files for all functions including the definitions of constants, global variables and structures
#include "../Headers/header.h"

// Include the main libraries
#include "../Core/Core.cpp"
#include "../IO/IO.cpp"
#include "../Tight_Binding/Tight_Binding.cpp"
#include "../AtomicOperators/AtomicAux.cpp"

// Include files that should be reoredered into the main libraries
#include "../Model/StandardOperators.cpp"
#include "../Hubbard_Model/CellFunctions.cpp"
#include "../DMFT/SingleSiteFunctions.cpp"
#include "../Yi/functions_yi.cpp"
#include "../Martin/functions_martin.cpp"
#include "../Simon/Simon.cpp"

// Include Plot functions
#include "../Plot/Plot.cpp"
#include "../Plot/LuaBinding.cpp"

// Include the Lua interface
#include "../LuaInterface/lua/includeLua.c"
#include "../LuaInterface/lua_extensions/tools.cpp"
#include "../LuaInterface/lua_extensions/luaInOut.cpp"
#include "../LuaInterface/lua_extensions/TableToArray.cpp"
#include "../LuaInterface/lua_extensions/get_put.cpp"
#include "../LuaInterface/units.cpp"
#include "../LuaInterface/Complex.cpp"
#include "../LuaInterface/Interpolate.cpp" // Interpolation using Polymorphism
#include "../LuaInterface/Operator.cpp"
#include "../LuaInterface/Wavefunction.cpp"
#include "../LuaInterface/Spectra.cpp"
#include "../LuaInterface/Matrix.cpp"
#include "../LuaInterface/TightBinding.cpp"
#include "../LuaInterface/ResponseFunction.cpp"
#include "../LuaInterface/Plot.cpp"
//#include "../LuaInterface/PostScript.cpp"
#include "../LuaInterface/BraKet.cpp"
#include "../LuaInterface/eval.cpp"
#include "../LuaInterface/ssp.cpp"






int main(int argc, char *argv[])
{
  LineType     Line;
  CurveType    Curve;
  TextType     Text;
  char         Label[1024];
  ArcType      Arc;
  BitMapType   BitMap;
  GraphicsType G;
  PointType    Margins[2];
  PointType    PlotRange[2];
  uint8_t      Color[4];
  PointType    Data[20];
  FILE *fptr;
  int i;
  char** PS;
  char* VPS=NULL;
  
  InitGraphics(&G);
  Color[0]=0;
  Color[1]=0;
  Color[2]=0;
  Color[3]=0xFF;
  Margins[0].X=0.08;
  Margins[0].Y=0.08;
  Margins[1].X=0.01;
  Margins[1].Y=0.01*sqrt(2.0);
  PlotRange[0].X=-1.1;
  PlotRange[0].Y=-1.2;
  PlotRange[1].X=5.1;
  PlotRange[1].Y=0.9;
  for(i=0;i<20;i++)
  {
    Data[i].X = i*0.5-3;
    Data[i].Y = sin(Data[i].X * 0.5);
  }
  
  AddPlotAxisToGraphics(&G,PlotRange,Margins,Color,0.001);
  AddXTicksToGraphics(&G, PlotRange, Margins, Color, 0.001, Color, 0.5, -1, 1, 0.01 * sqrt(2), 9);
  AddYTicksToGraphics(&G, PlotRange, Margins, Color, 0.001, Color, 0.5, -2, 0.5, 0.01, 4);
  AddDataPointsToGraphics(&G, PlotRange, Margins, Data, 20, 0.001, Color, 0.002, Color);
  AddXTicksLabelToGraphics(&G, PlotRange, Margins, Color, 0.03, TIMES_NEW_ROMAN, -2, 1, "%.0f");
  AddYTicksLabelToGraphics(&G, PlotRange, Margins, Color, 0.03, TIMES_NEW_ROMAN, -2, 0.5, "%.1f");
  InitText(&Text);
  Text.String = Label;
  strcpy(Label, "Whatever is on the X axis [units]");
  Text.Point.X=0;
  Text.Point.Y=0;
  Text.Angle=0;
  Text.Alignment='c';
  Text.Font=TIMES_NEW_ROMAN;
  Text.Size=0.03;
  AddXLabelToGraphics(&G, Margins, &Text);
  strcpy(Label, "Whatever is on the Y axis [units]");
  Text.Angle=0.5*Pi;
  AddYLabelToGraphics(&G, Margins, &Text);
  
  
  PS = &VPS;
  GraphicsToPS(&G, PS);
  fptr = fopen("test.ps", "w");
  fprintf(fptr,"%s\n",*PS);
  fclose(fptr);
  
  
  
  
  return 0;
  
  
  
  
  
  
  
  InitLine(&Line, 4);
  Line.Point[0].X=0.1;
  Line.Point[0].Y=0.1;
  Line.Point[1].X=0.6;
  Line.Point[1].Y=0.1;
  Line.Point[2].X=0.6;
  Line.Point[2].Y=0.5;
  Line.Point[3].X=0.1;
  Line.Point[3].Y=0.5;
  AddLine(&G, &Line);
  
  Line.Color[0]=255;   Line.Color[1]=0;   Line.Color[2]=0;
  Line.FillColor[0]=255;   Line.FillColor[1]=255;   Line.FillColor[2]=0;   Line.FillColor[3]=0;
  Line.Width=0.01;
  Line.NPoints=3;
  Line.CloseAndFill=true;
  AddLine(&G, &Line);
  
  Line.Point[0].X=0.1;
  Line.Point[0].Y=0.1;
  Line.Point[1].X=0.1;
  Line.Point[1].Y=0.2;
  Line.Point[2].X=0.2;
  Line.Point[2].Y=0.5;
  Line.Point[3].X=0.3;
  Line.Point[3].Y=0.5;
  Line.Width=0.005;
  Line.NPoints=4;
  Line.Color[0]=0;   Line.Color[1]=0;   Line.Color[2]=255;   Line.FillColor[3]=255;
  Line.FillColor[0]=0;   Line.FillColor[1]=128;   Line.FillColor[2]=64;
  AddLine(&G, &Line);

  
  
  InitCurve(&Curve, 4);
  Curve.Point[0].X=0.9;
  Curve.Point[0].Y=0.9;
  Curve.Point[1].X=0.6;
  Curve.Point[1].Y=0.9;
  Curve.Point[2].X=0.6;
  Curve.Point[2].Y=0.5;
  Curve.Point[3].X=0.9;
  Curve.Point[3].Y=0.5;
  
  
  Curve.VL[0].X=0.7;
  Curve.VL[0].Y=0.9;
  Curve.VL[1].X=0.5;
  Curve.VL[1].Y=0.8;
  Curve.VL[2].X=0.6;
  Curve.VL[2].Y=0.6;
  Curve.VL[3].X=0.8;
  Curve.VL[3].Y=0.5;
  
  
  Curve.VR[0].X=1.1;
  Curve.VR[0].Y=0.9;
  Curve.VR[1].X=0.7;
  Curve.VR[1].Y=1.0;
  Curve.VR[2].X=0.6;
  Curve.VR[2].Y=0.7;
  Curve.VR[3].X=1.0;
  Curve.VR[3].Y=0.5;
  
  Curve.CloseAndFill=true;
  AddCurve(&G, &Curve);
  
  InitArc(&Arc);
  Arc.Color[0]=255;
  Arc.Color[1]=128;
  Arc.Color[2]=128;
  Arc.Point.X=0.5;
  Arc.Point.Y=0.5;
  Arc.R=0.3;
  Arc.Width=0.01;
  Arc.CloseAndFill=false;
  Arc.ThetaR=0;
  Arc.ThetaL=2.0*Pi;
  AddArc(&G, &Arc);
  
  InitBitMap(&BitMap);
  BitMap.NpointsX=3;
  BitMap.NpointsY=2;
  BitMap.LB.X=0.1;
  BitMap.LB.Y=0.1;
  BitMap.RT.X=0.6;
  BitMap.RT.Y=0.5;
  
  BitMap.Pixels=(uint32_t*)malloc(3*2*sizeof(uint32_t));
  
  BitMap.Pixels[0] = 0xFFFF00;
  BitMap.Pixels[1] = 0xFF00FF;
  BitMap.Pixels[2] = 0x00FFFF;
  BitMap.Pixels[3] = 0xFF0000;
  BitMap.Pixels[4] = 0x00FF00;
  BitMap.Pixels[5] = 0x0000FF;

  AddBitMap(&G,&BitMap);
  
  
  Text.Color[0]=128;
  Text.Color[1]=255;
  Text.Color[2]=64;
  Text.Alignment='c';
  Text.Point.X = 0.5;
  Text.Point.Y = 0.5;
  Text.Angle=0.25*Pi;
  Text.Font=TIMES_NEW_ROMAN;
  Text.Size=0.05;
  Text.String=(char*)malloc(1024*sizeof(char));
  strcpy(Text.String, "Hello world");
  
  AddText(&G, &Text);
  
  
  
  
  
  PS = &VPS;
  GraphicsToPS(&G, PS);
  
  printf("%s\n",*PS);
  return 0;
}
