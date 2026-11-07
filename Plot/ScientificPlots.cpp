

int AddXTicksToGraphics(GraphicsType* Graphics, PointType *PlotRange, PointType *Margins, uint32_t Color, double LineWidth, uint32_t MinorColor, double MinorScaling, double TickStart, double dTick, double TickSize, int NMinorTicks)
{
  if((LineWidth==0) || (TickSize==0) || ( (!ColorVisible(Color)) && (!ColorVisible(MinorColor)) ) ) return 0;
  double Mval,Xval;
  int i;
  double const Scale = (1-Margins[0].X-Margins[1].X) / (PlotRange[1].X - PlotRange[0].X);
  LineType Line, MLine;
  InitLine(&Line,2);
  InitLine(&MLine,2);
  Line.Color = Color;
  Line.Width = LineWidth;
  MLine.Color = MinorColor;
  MLine.Width = LineWidth * MinorScaling;
  
  Line.Point[0].Y = Margins[0].Y;
  Line.Point[1].Y = Margins[0].Y + TickSize;
  MLine.Point[0].Y = Margins[0].Y;
  MLine.Point[1].Y = Margins[0].Y + TickSize * MinorScaling;
  for(Xval=TickStart-dTick;Xval<PlotRange[1].X;Xval+=dTick)
  {
    if( Xval > PlotRange[0].X + 1E-8 * (PlotRange[1].X-PlotRange[0].X) )
    {
      Line.Point[0].X = (Xval-PlotRange[0].X) * Scale + Margins[0].X;
      Line.Point[1].X = (Xval-PlotRange[0].X) * Scale + Margins[0].X;
      AddLine(Graphics, &Line);
    }
    Mval = Xval;
    for(i=0;i<NMinorTicks;i++)
    {
      Mval += dTick / (NMinorTicks+1.);
      if( (Mval > PlotRange[0].X + 1E-8 * (PlotRange[1].X-PlotRange[0].X) ) && (Mval < PlotRange[1].X - 1E-8 * (PlotRange[1].X-PlotRange[0].X) ) )
      {
        MLine.Point[0].X = (Mval-PlotRange[0].X) * Scale + Margins[0].X;
        MLine.Point[1].X = (Mval-PlotRange[0].X) * Scale + Margins[0].X;
        AddLine(Graphics, &MLine);
      }
    }
  }
  Line.Point[0].Y = 1-Margins[1].Y;
  Line.Point[1].Y = 1-Margins[1].Y - TickSize;
  MLine.Point[0].Y = 1-Margins[1].Y;
  MLine.Point[1].Y = 1-Margins[1].Y - TickSize * MinorScaling;
  for(Xval=TickStart-dTick;Xval<PlotRange[1].X - 1E-8 * (PlotRange[1].X-PlotRange[0].X);Xval+=dTick)
  {
    if( Xval > PlotRange[0].X + 1E-8 * (PlotRange[1].X-PlotRange[0].X) )
    {
      Line.Point[0].X = (Xval-PlotRange[0].X) * Scale + Margins[0].X;
      Line.Point[1].X = (Xval-PlotRange[0].X) * Scale + Margins[0].X;
      AddLine(Graphics, &Line);
    }
    Mval = Xval;
    for(i=0;i<NMinorTicks;i++)
    {
      Mval += dTick / (NMinorTicks+1.);
      if( (Mval > PlotRange[0].X + 1E-8 * (PlotRange[1].X-PlotRange[0].X) ) && (Mval < PlotRange[1].X - 1E-8 * (PlotRange[1].X-PlotRange[0].X) ) )
      {
        MLine.Point[0].X = (Mval-PlotRange[0].X) * Scale + Margins[0].X;
        MLine.Point[1].X = (Mval-PlotRange[0].X) * Scale + Margins[0].X;
        AddLine(Graphics, &MLine);
      }
    }
  }
  FreeLine(&Line);
  FreeLine(&MLine);
  return 0;
}


int AddYTicksToGraphics(GraphicsType* Graphics, PointType *PlotRange, PointType *Margins, uint32_t Color, double LineWidth, uint32_t MinorColor, double MinorScaling, double TickStart, double dTick, double TickSize, int NMinorTicks)
{
  if((LineWidth==0) || (TickSize==0) || ( (!ColorVisible(Color)) && (!ColorVisible(MinorColor)) ) ) return 0;
  double Mval,Yval;
  int i;
  double const Scale = (1-Margins[0].Y-Margins[1].Y) / (PlotRange[1].Y - PlotRange[0].Y);
  LineType Line, MLine;
  InitLine(&Line,2);
  InitLine(&MLine,2);
  Line.Color = Color;
  Line.Width = LineWidth;
  MLine.Color = MinorColor;
  MLine.Width = LineWidth * MinorScaling;
  
  Line.Point[0].X = Margins[0].X;
  Line.Point[1].X = Margins[0].X + TickSize;
  MLine.Point[0].X = Margins[0].X;
  MLine.Point[1].X = Margins[0].X + TickSize * MinorScaling;
  for(Yval=TickStart-dTick;Yval<PlotRange[1].Y-1E-8 * (PlotRange[1].Y-PlotRange[0].Y);Yval+=dTick)
  {
    if( Yval > PlotRange[0].Y + 1E-8 * (PlotRange[1].Y-PlotRange[0].Y) )
    {
      Line.Point[0].Y = (Yval-PlotRange[0].Y) * Scale + Margins[0].Y;
      Line.Point[1].Y = (Yval-PlotRange[0].Y) * Scale + Margins[0].Y;
      AddLine(Graphics, &Line);
    }
    Mval = Yval;
    for(i=0;i<NMinorTicks;i++)
    {
      Mval += dTick / (NMinorTicks+1.);
      if( (Mval > PlotRange[0].Y + 1E-8 * (PlotRange[1].Y-PlotRange[0].Y)) && (Mval < PlotRange[1].Y - 1E-8 * (PlotRange[1].Y-PlotRange[0].Y)) )
      {
        MLine.Point[0].Y = (Mval-PlotRange[0].Y) * Scale + Margins[0].Y;
        MLine.Point[1].Y = (Mval-PlotRange[0].Y) * Scale + Margins[0].Y;
        AddLine(Graphics, &MLine);
      }
    }
  }
  Line.Point[0].X = 1-Margins[1].X;
  Line.Point[1].X = 1-Margins[1].X - TickSize;
  MLine.Point[0].X = 1-Margins[1].X;
  MLine.Point[1].X = 1-Margins[1].X - TickSize * MinorScaling;
  for(Yval=TickStart-dTick;Yval<PlotRange[1].Y-1E-8 * (PlotRange[1].Y-PlotRange[0].Y);Yval+=dTick)
  {
    if( Yval > PlotRange[0].Y + 1E-8 * (PlotRange[1].Y-PlotRange[0].Y))
    {
      Line.Point[0].Y = (Yval-PlotRange[0].Y) * Scale + Margins[0].Y;
      Line.Point[1].Y = (Yval-PlotRange[0].Y) * Scale + Margins[0].Y;
      AddLine(Graphics, &Line);
    }
    Mval = Yval;
    for(i=0;i<NMinorTicks;i++)
    {
      Mval += dTick / (NMinorTicks+1.);
      if( (Mval > PlotRange[0].Y + 1E-8 * (PlotRange[1].Y-PlotRange[0].Y)) && (Mval < PlotRange[1].Y - 1E-8 * (PlotRange[1].Y-PlotRange[0].Y)) )
      {
        MLine.Point[0].Y = (Mval-PlotRange[0].Y) * Scale + Margins[0].Y;
        MLine.Point[1].Y = (Mval-PlotRange[0].Y) * Scale + Margins[0].Y;
        AddLine(Graphics, &MLine);
      }
    }
  }
  FreeLine(&Line);
  FreeLine(&MLine);
  return 0;
}

int AddXTicksLabelToGraphics(GraphicsType* Graphics, PointType *PlotRange, PointType *Margins, uint32_t Color, double TextSize, const char* TextFont, double TickStart, double dTick, const char* format)
{
  if((TextSize==0) || (!ColorVisible(Color)) || (TextFont==NULL) ) return 0;
  double Xval;
  double const Scale = (1-Margins[0].X-Margins[1].X) / (PlotRange[1].X - PlotRange[0].X);
  char String[1024];
  TextType Text;
  InitText(&Text);
  Text.Color = Color;
  Text.Size=TextSize;
  Text.String = String;
  Text.Point.Y = Margins[0].Y-TextSize;
  Text.Alignment='c';
  Text.Angle=0;
  Text.Font = TextFont;
  for(Xval=TickStart-dTick;Xval<PlotRange[1].X-1E-8 * (PlotRange[1].X-PlotRange[0].X) ;Xval+=dTick)
  {
    if( Xval > PlotRange[0].X + 1E-8 * (PlotRange[1].X-PlotRange[0].X) )
    {
      Text.Point.X = (Xval-PlotRange[0].X) * Scale + Margins[0].X;
      sprintf(String, format, Xval);
      AddText(Graphics, &Text);
    }
  }
  return 0;
}

int AddYTicksLabelToGraphics(GraphicsType* Graphics, PointType *PlotRange, PointType *Margins, uint32_t Color, double TextSize, const char* TextFont, double TickStart, double dTick, const char* format)
{
  if((TextSize==0) || (!ColorVisible(Color)) || (TextFont==NULL) ) return 0;
  double Yval;
  double const Scale = (1-Margins[0].Y-Margins[1].Y) / (PlotRange[1].Y - PlotRange[0].Y);
  char String[1024];
  TextType Text;
  InitText(&Text);
  Text.Color = Color;
  Text.Size=TextSize;
  Text.String = String;
  Text.Point.X = Margins[0].X-0.1*TextSize;
  Text.Alignment='r';
  Text.Angle=0;
  Text.Font = TextFont;
  for(Yval=TickStart-dTick;Yval<PlotRange[1].Y-1E-8 * (PlotRange[1].Y-PlotRange[0].Y) ;Yval+=dTick)
  {
    if( Yval > PlotRange[0].Y + 1E-8 * (PlotRange[1].Y-PlotRange[0].Y) )
    {
      Text.Point.Y = (Yval-PlotRange[0].Y) * Scale + Margins[0].Y - 0.3*TextSize;
      sprintf(String, format, Yval);
      AddText(Graphics, &Text);
    }
  }
  return 0;
}

int AddXLabelToGraphics(GraphicsType* Graphics, PointType *Margins, TextType *Text)
{
  if((Text->Size==0) || (!ColorVisible(Text->Color)) || (Text->Font==NULL) ) return 0;
  double shiftX, shiftY;
  shiftY = 0.4*Text->Size;
  shiftX = (1+Margins[0].X - Margins[1].X)/2;
  Text->Point.X += shiftX;
  Text->Point.Y += shiftY;
  AddText(Graphics, Text);
  Text->Point.X -= shiftX;
  Text->Point.Y -= shiftY;
  return 0;
}

int AddYLabelToGraphics(GraphicsType* Graphics, PointType *Margins, TextType *Text)
{
  if((Text->Size==0) || (!ColorVisible(Text->Color)) || (Text->Font==NULL) ) return 0;
  double shiftX, shiftY;
  shiftX = 0.8*Text->Size;
  shiftY = (1+Margins[0].Y - Margins[1].Y)/2;
  Text->Point.X += shiftX;
  Text->Point.Y += shiftY;
  AddText(Graphics, Text);
  Text->Point.X -= shiftX;
  Text->Point.Y -= shiftY;
  return 0;
}

// Input:
// Graphics an initiallized Grahics object
// Plotrange, an array of length 2 containign the range of the data.
//   Point [0] contains XMin and YMin, point 1 contains XMax and YMax.
// Margins, an array of length 2 containign the plot margins.
//   Point [0] contains the left bottom point, in fractions of canvas size, Point [1] contains the right top point in fractions of the canvas size
// LineWidth width in fraction of canvas size (sqrt(Xsize*Ysize)
// LineColor RBGalpha color def
//
int AddPlotAxisToGraphics(GraphicsType* Graphics, PointType *PlotRange, PointType *Margins, uint32_t Color, double Width, double XAxis, double YAxis )
{
  if((Width==0) || (!ColorVisible(Color)) ) return 0;
  LineType Line;
  if( InitLine(&Line, 4) ) return 1;
  Line.Color = Color;
  Line.FillColor = 0xFF000000; // no filling while filling fully transparent.
  Line.CloseAndFill = true;
  Line.Width = Width;
  Line.Point[0].X = Margins[0].X;
  Line.Point[0].Y = Margins[0].Y;
  Line.Point[1].X = Margins[0].X;
  Line.Point[1].Y = 1-Margins[1].Y;
  Line.Point[2].X = 1-Margins[1].X;
  Line.Point[2].Y = 1-Margins[1].Y;
  Line.Point[3].X = 1-Margins[1].X;
  Line.Point[3].Y = Margins[0].Y;
  if( AddLineInherentLineMemory(Graphics, &Line) ) return 1;
  if( (XAxis > PlotRange[0].Y) && (XAxis < PlotRange[1].Y) )
  {
    if( InitLine(&Line, 2) ) return 1;
    Line.Color = Color;
    Line.FillColor = 0xFF000000; // no filling while filling fully transparent.
    Line.CloseAndFill = false;
    Line.Width = Width;
    Line.Point[0].X = Margins[0].X;
    Line.Point[0].Y = Margins[0].Y + (1-Margins[1].Y-Margins[0].Y) * ( XAxis - PlotRange[0].Y) / (PlotRange[1].Y - PlotRange[0].Y);
    Line.Point[1].X = 1-Margins[1].X;
    Line.Point[1].Y = Margins[0].Y + (1-Margins[1].Y-Margins[0].Y) * ( XAxis - PlotRange[0].Y) / (PlotRange[1].Y - PlotRange[0].Y);
    if( AddLineInherentLineMemory(Graphics, &Line) ) return 1;
  }
  if( (YAxis > PlotRange[0].X) && (YAxis < PlotRange[1].X) )
  {
    if( InitLine(&Line, 2) ) return 1;
    Line.Color = Color;
    Line.FillColor = 0xFF000000; // no filling while filling fully transparent.
    Line.CloseAndFill = false;
    Line.Width = Width;
    Line.Point[0].X = Margins[0].X + (1-Margins[1].X-Margins[0].X) * ( YAxis - PlotRange[0].X) / (PlotRange[1].X - PlotRange[0].X);
    Line.Point[0].Y = Margins[0].X;
    Line.Point[1].X = Margins[0].X + (1-Margins[1].X-Margins[0].X) * ( YAxis - PlotRange[0].X) / (PlotRange[1].X - PlotRange[0].X);
    Line.Point[1].Y = 1-Margins[1].X;
    if( AddLineInherentLineMemory(Graphics, &Line) ) return 1;
  }
  return 0;
}

void FindBoundaryPoint(PointType *Boundary, PointType *In, PointType *Out, PointType *OnBoundary)
{
  double alpha,alphax,alphay;
  // a segment that does not move along x (or y) crosses no boundary in that
  // direction; dividing would give 0/0 = NaN for In on that boundary
  alphax = INFINITY;
  if ( Out->X != In->X )
  {
    alphax = ( Boundary[0].X - In->X ) / (Out->X - In->X);
    if ( alphax < 0 ) alphax = ( Boundary[1].X - In->X ) / (Out->X - In->X);
  }
  alphay = INFINITY;
  if ( Out->Y != In->Y )
  {
    alphay = ( Boundary[0].Y - In->Y ) / (Out->Y - In->Y);
    if ( alphay < 0 ) alphay = ( Boundary[1].Y - In->Y ) / (Out->Y - In->Y);
  }
  alpha = Min(alphax, alphay);
  OnBoundary->X = alpha * Out->X + (1-alpha) * In->X;
  OnBoundary->Y = alpha * Out->Y + (1-alpha) * In->Y;
}

int AddDataPointsToGraphics(GraphicsType* Graphics, PointType *PlotRange, PointType *Margins, PointType *Data, int NPoints, double LineWidth, uint32_t LineColor, double PointSize, uint32_t PointColor)
{
  if( ( (LineWidth==0) || (!ColorVisible(LineColor)) ) && ( (PointSize==0) || (!ColorVisible(PointColor)) ) ) return 0;
  LineType Line;
  ArcType Arc;
  int i;
  bool render;
  PointType Ofset, Scale, OnBoundary;
  Scale.X = (1-Margins[0].X-Margins[1].X) / (PlotRange[1].X - PlotRange[0].X);
  Scale.Y = (1-Margins[0].Y-Margins[1].Y) / (PlotRange[1].Y - PlotRange[0].Y);
  Ofset.X = Margins[0].X - PlotRange[0].X * Scale.X;
  Ofset.Y = Margins[0].Y - PlotRange[0].Y * Scale.Y;
  if((LineWidth!=0) && (ColorVisible(LineColor)))
  {
    Line.Color = LineColor;
    Line.Width = LineWidth;
    Line.NPoints = NPoints;
    Line.Point = (PointType*)malloc(NPoints * sizeof(PointType));
    Line.CloseAndFill = false;
    if(Line.Point==NULL)
    {
      printf("Malloc failed in AddDataPointsToGraphics\n");
      return 1;
    }
    render=false;
    for (i=0;i<NPoints;i++)
    {
      if( (Data[i].X >= PlotRange[0].X) && (Data[i].Y >= PlotRange[0].Y) && (Data[i].X <= PlotRange[1].X) && ((Data[i].Y <= PlotRange[1].Y)) )
      { // inside plot range
        if(render==false)
        { // new plot segment starts
          Line.NPoints=0;
          if(i!=0)
          {
            FindBoundaryPoint(PlotRange, &(Data[i]), &(Data[i-1]), &OnBoundary);
            Line.Point[Line.NPoints].X = Ofset.X + Scale.X * OnBoundary.X;
            Line.Point[Line.NPoints].Y = Ofset.Y + Scale.Y * OnBoundary.Y;
            Line.NPoints++;
          }
          Line.Point[Line.NPoints].X = Ofset.X + Scale.X * Data[i].X;
          Line.Point[Line.NPoints].Y = Ofset.Y + Scale.Y * Data[i].Y;
          Line.NPoints++;
        }
        else
        {
          Line.Point[Line.NPoints].X = Ofset.X + Scale.X * Data[i].X;
          Line.Point[Line.NPoints].Y = Ofset.Y + Scale.Y * Data[i].Y;
          Line.NPoints++;
        }
        render=true;
      }
      else
      { // outside plot range
        if(render==true)
        { // find the crossing between the line connecting point i and i-1 and the boundary.
          FindBoundaryPoint(PlotRange, &(Data[i-1]), &(Data[i]), &OnBoundary);
          Line.Point[Line.NPoints].X = Ofset.X + Scale.X * OnBoundary.X;
          Line.Point[Line.NPoints].Y = Ofset.Y + Scale.Y * OnBoundary.Y;
          Line.NPoints++;
          render=false;
          AddLine(Graphics, &Line);
          Line.NPoints=0;
        }
      }
    }
    if(render==true)
    {
      AddLine(Graphics, &Line);
    }
    FreeLine(&Line);
  }
  if((PointSize!=0) && (ColorVisible(PointColor)))
  {
    Arc.Color = PointColor;
    Arc.FillColor = PointColor;
    Arc.CloseAndFill = true;
    Arc.ThetaR = 0;
    Arc.ThetaL = 2*Pi;
    Arc.R = PointSize;
    Arc.Width = 0;
    for(i=0;i<NPoints;i++)
    {
      if( (Data[i].X >= PlotRange[0].X) && (Data[i].Y >= PlotRange[0].Y) && (Data[i].X <= PlotRange[1].X) && (Data[i].Y <= PlotRange[1].Y) )
      {
        Arc.Point.X = Ofset.X + Scale.X * Data[i].X;
        Arc.Point.Y = Ofset.Y + Scale.Y * Data[i].Y;
        AddArc(Graphics, &Arc);
      }
    }
  }
  return 0;
}

// Data->N+1 values at the energies X, or on the even mesh from Data->Emin to Data->Emax if X is NULL
int AddSpectraTypeToGraphics(GraphicsType* Graphics, PointType *PlotRange, PointType *Margins, SpectraType *Data, double LineWidthRe, uint32_t LineColorRe, double PointSizeRe, uint32_t PointColorRe, double LineWidthIm, uint32_t LineColorIm, double PointSizeIm, uint32_t PointColorIm, double YScale, double YShift, double const *X)
{
  PointType *DataPoints;
  int NPoints=Data->N+1;
  int i;
  DataPoints = (PointType*)malloc(NPoints*sizeof(PointType));
  if(DataPoints==NULL)
  {
    printf("Malloc failed in AddSpectraTypeToGraphics\n");
    return 1;
  }
  for(i=0;i<NPoints;i++)
  {
    if(X!=NULL)         DataPoints[i].X = X[i];
    else if(NPoints==1) DataPoints[i].X = Data->Emin;
    else                DataPoints[i].X = i * Data->Emax / (NPoints-1) + (NPoints-i-1) * Data->Emin / (NPoints-1);
    DataPoints[i].Y = YScale * Data->Re[i] + YShift;
  }
  AddDataPointsToGraphics(Graphics, PlotRange, Margins, DataPoints, NPoints, LineWidthRe, LineColorRe, PointSizeRe, PointColorRe);
  for(i=0;i<NPoints;i++)
  {
    DataPoints[i].Y = YScale * Data->Im[i] + YShift;
  }
  AddDataPointsToGraphics(Graphics, PlotRange, Margins, DataPoints, NPoints, LineWidthIm, LineColorIm, PointSizeIm, PointColorIm);
  free(DataPoints);
  return 0;
}

int PlotAllButDataToGraphics(GraphicsType* Graphics, FramePlotOptions *FPO)
{
  TextType Text;
  AddPlotAxisToGraphics(Graphics,FPO->PlotRange,FPO->Margins,FPO->LineColor,FPO->LineWidth, FPO->XAxis, FPO->YAxis);
  if(FPO->XTicks)
  {
    AddXTicksToGraphics(Graphics, FPO->PlotRange, FPO->Margins, FPO->LineColor, FPO->LineWidth, FPO->LineColor, 0.5, FPO->FirstXTick, FPO->dXTick, 0.01 * sqrt(2), FPO->NMinorXTick);
    AddXTicksLabelToGraphics(Graphics, FPO->PlotRange, FPO->Margins, FPO->LineColor, FPO->FontSize, FPO->Font, FPO->FirstXTick, FPO->dXTick, FPO->PlotXScaleFormat);
  }
  if(FPO->YTicks)
  {
    AddYTicksToGraphics(Graphics, FPO->PlotRange, FPO->Margins, FPO->LineColor, FPO->LineWidth, FPO->LineColor, 0.5, FPO->FirstYTick, FPO->dYTick, 0.01          , FPO->NMinorYTick);
    AddYTicksLabelToGraphics(Graphics, FPO->PlotRange, FPO->Margins, FPO->LineColor, FPO->FontSize, FPO->Font, FPO->FirstYTick, FPO->dYTick, FPO->PlotYScaleFormat);
  }
  InitText(&Text);
  Text.String = FPO->XLabel;
  Text.Point.X=0;
  Text.Point.Y=0;
  Text.Angle=0;
  Text.Alignment='c';
  Text.Font=FPO->Font;
  Text.Size=FPO->FontSize;
  Text.Color = FPO->LineColor;
  AddXLabelToGraphics(Graphics, FPO->Margins, &Text);
  Text.String = FPO->YLabel;
  Text.Angle=0.5*Pi;
  AddYLabelToGraphics(Graphics, FPO->Margins, &Text);
  return 0;
}

int AddVerticleLinesForSpecialkPoints(GraphicsType* Graphics, FramePlotOptions *FPO, BandstructurePathType *kPath, double EmptyPathWidth)
{
  char String[2];
  String[1]='\0';
  TextType Text;
  InitText(&Text);
  Text.Color = FPO->LineColor;
  Text.Size=FPO->FontSize;
  Text.String = String;
  Text.Point.Y = FPO->Margins[0].Y-FPO->FontSize;
  Text.Alignment='c';
  Text.Angle=0;
  Text.Font = FPO->Font;
  
  LineType Line;
  if( InitLine(&Line, 2) ) return 1;
  Line.Color = FPO->LineColor;
  Line.FillColor = 0xFF000000; // no filling while filling fully transparent.
  Line.CloseAndFill = false;
  Line.Width = FPO->LineWidth;
  Line.Point[0].Y = FPO->Margins[0].Y;
  Line.Point[1].Y = 1-FPO->Margins[1].Y;
  
  
  double LengthKPath = 0;
  int    NZeroDistance = 0;
  for(int i=0;i<(kPath->NQ)-1;i++)
  {
    if(kPath->Q[i].Nq==0)
    {
      NZeroDistance++;
    }
    else
    {
      LengthKPath += sqrt( Sqr(kPath->Q[i].q[0] - kPath->Q[i+1].q[0]) + Sqr(kPath->Q[i].q[1] - kPath->Q[i+1].q[1]) + Sqr(kPath->Q[i].q[2] - kPath->Q[i+1].q[2]));
    }
  }
  double LengthForKPath = (1 - (NZeroDistance * EmptyPathWidth)) * (1 - FPO->Margins[0].X - FPO->Margins[1].X);
  double XPos = FPO->Margins[0].X;
  Text.String[0] = kPath->Q[0].label;
  Text.Point.X = XPos;
  AddText(Graphics, &Text);
  for(int i=0;i<(kPath->NQ)-1;i++)
  {
    if(kPath->Q[i].Nq==0)
    {
      XPos += EmptyPathWidth * (1 - FPO->Margins[0].X - FPO->Margins[1].X);
    }
    else
    {
      XPos += LengthForKPath * ( sqrt( Sqr(kPath->Q[i].q[0] - kPath->Q[i+1].q[0]) + Sqr(kPath->Q[i].q[1] - kPath->Q[i+1].q[1]) + Sqr(kPath->Q[i].q[2] - kPath->Q[i+1].q[2])) / LengthKPath );
    }
    Text.String[0] = kPath->Q[i+1].label;
    Text.Point.X = XPos;
    AddText(Graphics, &Text);
    Line.Point[0].X = XPos;
    Line.Point[1].X = XPos;
    AddLine(Graphics, &Line);
  }
  return 0;
}




int PlotAllButDataToBandGraphics(GraphicsType* Graphics, FramePlotOptions *FPO, BandstructurePathType *kPath, double EmptyPathWidth)
{
  TextType Text;
  AddPlotAxisToGraphics(Graphics,FPO->PlotRange,FPO->Margins,FPO->LineColor,FPO->LineWidth, FPO->XAxis, FPO->YAxis);
  if(FPO->YTicks)
  {
    AddYTicksToGraphics(Graphics, FPO->PlotRange, FPO->Margins, FPO->LineColor, FPO->LineWidth, FPO->LineColor, 0.5, FPO->FirstYTick, FPO->dYTick, 0.01          , FPO->NMinorYTick);
    AddYTicksLabelToGraphics(Graphics, FPO->PlotRange, FPO->Margins, FPO->LineColor, FPO->FontSize, FPO->Font, FPO->FirstYTick, FPO->dYTick, FPO->PlotYScaleFormat);
  }
  InitText(&Text);
  Text.String = FPO->XLabel;
  Text.Point.X=0;
  Text.Point.Y=0;
  Text.Angle=0;
  Text.Alignment='c';
  Text.Font=FPO->Font;
  Text.Size=FPO->FontSize;
  Text.Color = FPO->LineColor;
  AddXLabelToGraphics(Graphics, FPO->Margins, &Text);
  Text.String = FPO->YLabel;
  Text.Angle=0.5*Pi;
  AddYLabelToGraphics(Graphics, FPO->Margins, &Text);
  AddVerticleLinesForSpecialkPoints(Graphics, FPO, kPath, EmptyPathWidth);
  return 0;
}

int PlotPointsToGraphics(GraphicsType* Graphics, PointType *Data, int NPoints, LinePlotOptions *LPO, FramePlotOptions *FPO)
{
  PlotAllButDataToGraphics(Graphics, FPO);
  AddDataPointsToGraphics(Graphics, FPO->PlotRange, FPO->Margins, Data, NPoints, LPO->LineWidth, LPO->LineColor, LPO->PointSize, LPO->PointColor);
  return 0;
}

int PlotSpectraTypeToGraphics(GraphicsType* Graphics, SpectraType *Data, LinePlotOptions *LPORe, LinePlotOptions *LPOIm, FramePlotOptions *FPO, double YScale)
{
  PlotAllButDataToGraphics(Graphics, FPO);
  AddSpectraTypeToGraphics(Graphics, FPO->PlotRange, FPO->Margins, Data, LPORe->LineWidth, LPORe->LineColor, LPORe->PointSize, LPORe->PointColor, LPOIm->LineWidth, LPOIm->LineColor, LPOIm->PointSize, LPOIm->PointColor, YScale, 0, NULL);
  return 0;
}

// The height of spectrum i in a waterfall plot or a colour map: WaterfallOffset plus
// i*WaterfallShift, or plus its resonant energy
void WaterfallPositions(SpectraTypeArray const *Data, double WaterfallOffset, double WaterfallShift, bool UseResonantEnergy, double *YPos)
{
  for(int i=0;i<(int)Data->NSpectra;i++)
  {
    if(UseResonantEnergy) YPos[i] = WaterfallOffset + Data->Spec[i].ResonantEnergy;
    else                  YPos[i] = WaterfallOffset + i*WaterfallShift;
  }
}

// spectrum i is drawn shifted up by YPos[i], or not shifted if YPos is NULL; its values are at the
// energies X[i], or on its even mesh if X or X[i] is NULL
int PlotSpectraTypeArrayToGraphics(GraphicsType* Graphics, SpectraTypeArray *Data, LinePlotOptions *LPORe, LinePlotOptions *LPOIm, FramePlotOptions *FPO, double YScale, double const *YPos, double * const *X)
{
  int i;
  PlotAllButDataToGraphics(Graphics, FPO);
  for(i=0;i<Data->NSpectra;i++)
  {
    AddSpectraTypeToGraphics(Graphics, FPO->PlotRange, FPO->Margins, &(Data->Spec[i]), LPORe[i].LineWidth, LPORe[i].LineColor, LPORe[i].PointSize, LPORe[i].PointColor, LPOIm[i].LineWidth, LPOIm[i].LineColor, LPOIm[i].PointSize, LPOIm[i].PointColor, YScale, (YPos==NULL)?0:YPos[i], (X==NULL)?NULL:X[i]);
  }
  return 0;
}

// gives the red, green and blue, each between 0 and 1, of the value Re + I Im
typedef void (*ColorFunctionType)(double Re, double Im, void *Context, double RGB[3]);

// an opaque pixel 0x00RRGGBB from red, green and blue between 0 and 1, clipped to that range
uint32_t RGBToPixel(double const RGB[3])
{
  uint32_t Pixel = 0x00; // the transparency, opaque
  for(int c=0;c<3;c++)
  {
    double v = RGB[c];
    if(!(v>0)) v=0;
    if(v>1) v=1;
    Pixel = (Pixel<<8) + (uint32_t)round(255*v);
  }
  return Pixel;
}

/*  Fills the pixels of a colour map of the spectra in Data
 *    Map          - NpointsX and NpointsY set and Pixels allocated; the map covers PlotRange
 *    Data, YPos   - spectrum i lies at height YPos[i]; at least two spectra, YPos strictly increasing or decreasing
 *    YScale       - the spectra are multiplied by YScale before they are coloured
 *    Interpolate  - linear interpolation in energy and between neighbouring spectra, otherwise the nearest
 *                   point, with the boundary between two spectra half way between them
 *    Color        - the colour of a value, Context is handed on to it
 *  Pixels without data are transparent white (0xFFFFFFFF). Row 0 of the pixels is the top of the map.
 */
int SpectraTypeArrayToColorMap(BitMapType *Map, SpectraTypeArray const *Data, double const *YPos, PointType const *PlotRange, double YScale, bool Interpolate, ColorFunctionType Color, void *Context)
{
  int const NX = Map->NpointsX;
  int const NY = Map->NpointsY;
  int const NS = (int)Data->NSpectra;
  if(NS<2) { printf("ERROR in SpectraTypeArrayToColorMap, a colour map needs at least two spectra\n"); return 1;}
  // the spectra in the order of increasing height
  bool const Increasing = YPos[NS-1] > YPos[0];
  int    *Order = (int*)malloc(NS*sizeof(int));
  double *Y     = (double*)malloc(NS*sizeof(double));
  double *Re    = (double*)malloc((size_t)NS*NX*sizeof(double));
  double *Im    = (double*)malloc((size_t)NS*NX*sizeof(double));
  if((Order==NULL)||(Y==NULL)||(Re==NULL)||(Im==NULL))
  {
    printf("ERROR in SpectraTypeArrayToColorMap, malloc failed\n");
    free(Order); free(Y); free(Re); free(Im);
    return 1;
  }
  for(int s=0;s<NS;s++)
  {
    Order[s] = Increasing ? s : NS-1-s;
    Y[s] = YPos[Order[s]];
  }
  // each spectrum on the centres of the pixel columns, NAN outside its energy range
  for(int s=0;s<NS;s++)
  {
    SpectraType const * const S = &(Data->Spec[Order[s]]);
    int const N = (int)S->N;
    for(int ix=0;ix<NX;ix++)
    {
      double const x = PlotRange[0].X + (ix+0.5) * (PlotRange[1].X-PlotRange[0].X) / NX;
      double const t = (N==0) ? 0 : (x - S->Emin) * N / (S->Emax - S->Emin);
      double re=NAN, im=NAN;
      if(Interpolate)
      {
        if((t>=0)&&(t<=N))
        {
          int k = (int)floor(t);
          if(k>=N) k=N-1;
          double const w = t-k;
          if(N==0)      { re = S->Re[0]; im = S->Im[0]; }
          else if(w==0) { re = S->Re[k]; im = S->Im[k]; }
          else          { re = (1-w)*S->Re[k] + w*S->Re[k+1]; im = (1-w)*S->Im[k] + w*S->Im[k+1]; }
        }
      }
      else
      {
        int const k = (int)floor(t+0.5);
        if((k>=0)&&(k<=N)) { re = S->Re[k]; im = S->Im[k]; }
      }
      Re[s*NX+ix] = YScale * re;
      Im[s*NX+ix] = YScale * im;
    }
  }
  // the rows, from the top down
  for(int iy=0;iy<NY;iy++)
  {
    double const y = PlotRange[1].Y - (iy+0.5) * (PlotRange[1].Y-PlotRange[0].Y) / NY;
    int s0=-1, s1=-1;
    double w=0;   // weight of spectrum s1
    if(Interpolate)
    {
      if((y>=Y[0])&&(y<=Y[NS-1]))
      {
        s0=0;
        while((s0<NS-2)&&(y>Y[s0+1])) s0++;
        s1 = s0+1;
        w = (y-Y[s0]) / (Y[s1]-Y[s0]);
      }
    }
    else
    {
      // spectrum s covers the heights from half way to its lower neighbour to half way to its upper neighbour
      double const Bottom = Y[0]    - (Y[1]-Y[0]) / 2;
      double const Top    = Y[NS-1] + (Y[NS-1]-Y[NS-2]) / 2;
      if((y>=Bottom)&&(y<=Top))
      {
        s0=0;
        while((s0<NS-1)&&(y>(Y[s0]+Y[s0+1])/2)) s0++;
        s1=s0;
      }
    }
    for(int ix=0;ix<NX;ix++)
    {
      uint32_t Pixel = 0xFFFFFFFF; // white, fully transparent
      if(s0>=0)
      {
        double re, im;
        if(w==0)      { re = Re[s0*NX+ix]; im = Im[s0*NX+ix]; }
        else if(w==1) { re = Re[s1*NX+ix]; im = Im[s1*NX+ix]; }
        else          { re = (1-w)*Re[s0*NX+ix] + w*Re[s1*NX+ix]; im = (1-w)*Im[s0*NX+ix] + w*Im[s1*NX+ix]; }
        if(!(isnan(re)||isnan(im)))
        {
          double RGB[3] = {1,1,1};
          Color(re, im, Context, RGB);
          Pixel = RGBToPixel(RGB);
        }
      }
      Map->Pixels[iy*NX+ix] = Pixel;
    }
  }
  free(Order); free(Y); free(Re); free(Im);
  return 0;
}

// A colour map of the spectra in Data inside the frame FPO; see SpectraTypeArrayToColorMap
// The graphics takes the pixels, so they are freed with it.
int PlotSpectraTypeArrayToColorMap(GraphicsType* Graphics, SpectraTypeArray const *Data, double const *YPos, FramePlotOptions *FPO, double YScale, int NPixelX, int NPixelY, bool Interpolate, ColorFunctionType Color, void *Context)
{
  BitMapType Map;
  InitBitMap(&Map);
  Map.LB.X = FPO->Margins[0].X;
  Map.LB.Y = FPO->Margins[0].Y;
  Map.RT.X = 1 - FPO->Margins[1].X;
  Map.RT.Y = 1 - FPO->Margins[1].Y;
  Map.NpointsX = NPixelX;
  Map.NpointsY = NPixelY;
  Map.Pixelated = !Interpolate;
  Map.Pixels = (uint32_t*)malloc((size_t)NPixelX*NPixelY*sizeof(uint32_t));
  if(Map.Pixels==NULL) { printf("ERROR in PlotSpectraTypeArrayToColorMap, malloc failed\n"); return 1;}
  if(SpectraTypeArrayToColorMap(&Map, Data, YPos, FPO->PlotRange, YScale, Interpolate, Color, Context)) { free(Map.Pixels); return 1;}
  if(AddBitMapInherentMemory(Graphics, &Map)) { free(Map.Pixels); return 1;}
  PlotAllButDataToGraphics(Graphics, FPO);
  return 0;
}

int PlotSpectraTypeArrayToGraphicsArray(GraphicsType* Graphics, SpectraTypeArray *Data, int NXPlots, int NYPlots, LinePlotOptions *LPORe, LinePlotOptions *LPOIm, FramePlotOptions *FPO, double YScale, double ScaleOffDiagonal, double ShiftOffDiagonal)
{
  if((NXPlots==1)&&(NYPlots==1)) return PlotSpectraTypeArrayToGraphics(Graphics,Data,LPORe,LPOIm,FPO,YScale,NULL,NULL);
  int i;
  FramePlotOptions FPOInset;
  // plot graphs to layer 1 to NXPlots * NYPlots
  // scale and shift the off diagonal plots
  // margin 0 is used for the left and botom margins, as one would for a single plot
  // margin 1 is used for the left and right or bottom and top margin
  FPOInset = *FPO;
  FPOInset.Margins[0] = FPOInset.Margins[1];
  FPOInset.FontSize = 0;
  FPOInset.Font = NULL;
  FPOInset.PlotRange[0].Y -= ShiftOffDiagonal;
  FPOInset.PlotRange[0].Y *= ScaleOffDiagonal;
  FPOInset.PlotRange[1].Y -= ShiftOffDiagonal;
  FPOInset.PlotRange[1].Y *= ScaleOffDiagonal;
  FPOInset.FirstYTick -= ShiftOffDiagonal;
  FPOInset.FirstYTick *= ScaleOffDiagonal;
  FPOInset.dYTick *= ScaleOffDiagonal;
  Graphics->NLayer = 1 + NXPlots * NYPlots;
  Graphics->Layer = (GraphicsType**)malloc(pow2roundup(Graphics->NLayer) * sizeof(GraphicsType*));
  Graphics->Layer[0] = Graphics;
  for(i=1;i<Graphics->NLayer;i++)
  {
    Graphics->Layer[i] = (GraphicsType*)malloc(sizeof(GraphicsType));
    InitGraphics(Graphics->Layer[i]);
    Graphics->Layer[i]->CanvasOrigin.X = Graphics->CanvasSize.X * ( FPO->Margins[0].X + ( (i-1) % NXPlots ) * (1-FPO->Margins[0].X) / NXPlots );
    Graphics->Layer[i]->CanvasOrigin.Y = Graphics->CanvasSize.Y * ( FPO->Margins[0].Y + ( NYPlots - ((i-1) / NXPlots ) -1 ) * (1-FPO->Margins[0].Y) / NYPlots );
    Graphics->Layer[i]->CanvasSize.X = Graphics->CanvasSize.X * (1-FPO->Margins[0].X) / NXPlots;
    Graphics->Layer[i]->CanvasSize.Y = Graphics->CanvasSize.Y * (1-FPO->Margins[0].Y) / NYPlots;
    if(( (i-1) % NXPlots ) == ( (i-1) / NXPlots ))
    {
      FPOInset.FirstYTick /= ScaleOffDiagonal;
      FPOInset.FirstYTick += ShiftOffDiagonal;
      FPOInset.dYTick /= ScaleOffDiagonal;
      FPOInset.PlotRange[0].Y /= ScaleOffDiagonal;
      FPOInset.PlotRange[0].Y += ShiftOffDiagonal;
      FPOInset.PlotRange[1].Y /= ScaleOffDiagonal;
      FPOInset.PlotRange[1].Y += ShiftOffDiagonal;
    }
    PlotSpectraTypeToGraphics(Graphics->Layer[i], &(Data->Spec[i-1]), LPORe, LPOIm, &FPOInset, YScale);
    if(( (i-1) % NXPlots ) == ( (i-1) / NXPlots ))
    {
      FPOInset.FirstYTick -= ShiftOffDiagonal;
      FPOInset.FirstYTick *= ScaleOffDiagonal;
      FPOInset.dYTick *= ScaleOffDiagonal;
      FPOInset.PlotRange[0].Y -= ShiftOffDiagonal;
      FPOInset.PlotRange[0].Y *= ScaleOffDiagonal;
      FPOInset.PlotRange[1].Y -= ShiftOffDiagonal;
      FPOInset.PlotRange[1].Y *= ScaleOffDiagonal;
    }
  }
  // plot tick and axis labels to layer 0
  FPOInset = *FPO;
  const double InsetXWidth   = (1 - FPO->Margins[0].X) / NXPlots;
  const double InsetXMargins = FPO->Margins[1].X * InsetXWidth;
  const double InsetYWidth   = (1 - FPO->Margins[0].Y) / NYPlots;
  const double InsetYMargins = FPO->Margins[1].Y * InsetYWidth;
  for (i=0;(i<NXPlots)&&(FPO->XTicks);i++)
  {
    FPOInset.Margins[0].X = FPO->Margins[0].X + InsetXMargins + i * InsetXWidth;
    FPOInset.Margins[1].X = 1 - FPO->Margins[0].X + InsetXMargins - (i+1) * InsetXWidth;
    AddXTicksLabelToGraphics(Graphics, FPOInset.PlotRange, FPOInset.Margins, FPOInset.LineColor, FPOInset.FontSize, FPOInset.Font, FPOInset.FirstXTick, FPOInset.dXTick, FPOInset.PlotXScaleFormat);
  }
  FPOInset.Margins[0].X = FPO->Margins[0].X;
  for (i=0;(i<NYPlots)&&(FPO->YTicks);i++)
  {
    FPOInset.Margins[0].Y = FPO->Margins[0].Y + InsetYMargins + i * InsetYWidth;
    FPOInset.Margins[1].Y = 1 - FPO->Margins[0].Y + InsetYMargins - (i+1) * InsetYWidth;
    AddYTicksLabelToGraphics(Graphics, FPOInset.PlotRange, FPOInset.Margins, FPOInset.LineColor, FPOInset.FontSize, FPOInset.Font, FPOInset.FirstYTick, FPOInset.dYTick, FPOInset.PlotYScaleFormat);
  }
  TextType Text;
  InitText(&Text);
  Text.String = FPO->XLabel;
  Text.Point.X=0;
  Text.Point.Y=0;
  Text.Angle=0;
  Text.Alignment='c';
  Text.Font=FPO->Font;
  Text.Size=FPO->FontSize;
  Text.Color = FPO->LineColor;
  AddXLabelToGraphics(Graphics, FPO->Margins, &Text);
  Text.String = FPO->YLabel;
  Text.Angle=0.5*Pi;
  AddYLabelToGraphics(Graphics, FPO->Margins, &Text);
  return 0;
}

/*  Given a matrix with the band Energies generates a plot of the band-structure as bitmap
 *  Input
 *    BandEnergies     - a matrix with Nkpoint rows(Ni) and NFermion columns(Nj)
 *    kPath            - a path in the Brillouine zone with Nkpoint points
 *    BandPlot         - an allocated bitmap with the dimensions set
 *      PointType LB;              // left bottom position in fraction of the canvas
 *      PointType RT;              // right top position in fraction of the canvas
 *      int NpointsX;              // number of pixels in the x direction
 *      int NpointsY;              // number of pixels in the y direction
 *      uint32_t *Pixels;          // pixels in format alphaRGB of size NpointsX * NpointsY
 *    Emin             - double for the lower bound of the plot
 *    Emax             - doulbe for the upper bound of the plot
 *    LineWidth        - double line width in ratio of the bitmap height, rounded to integer pixels
 *    LineColor        - RGB alpha color used for the lines
 *    ZeroDistance     - If there are zero k-points between two adjacent points in the kPath we add an aditional distance in ratio of the bitmap width, rounded to integer pixels
 */
int ConvertBandEnergiesToBitMap(CompactMatrixType* BandEnergies, BandstructurePathType *kPath, BitMapType *BandPlot, double Emin, double Emax, double LineWidth, uint32_t LineColor, double ZeroDistance)
{
  if(Emax < Emin) { printf("ERROR in ConvertBandEnergiesToBitMap, Emax < Emin\n"); return 1;}
  if(LineWidth<0) { printf("ERROR in ConvertBandEnergiesToBitMap, LineWidth < 0\n"); return 1;}
  if(ZeroDistance<0) { printf("ERROR in ConvertBandEnergiesToBitMap, ZeroDistance < 0\n"); return 1;}
  MakeBitMatWhiteTransparent(BandPlot);
  uint32_t const PixelColor = LineColor; // a pixel has the form of a colour, 0xAARRGGBB
  double LengthKPath = 0;
  int    NKPointsInPath=kPath->NQ;
  int    NZeroDistance = 0;
  for(int i=0;i<(kPath->NQ)-1;i++)
  {
    NKPointsInPath += kPath->Q[i].Nq;
    if(kPath->Q[i].Nq==0)
    {
      NZeroDistance++;
    }
    else
    {
      LengthKPath += sqrt( Sqr(kPath->Q[i].q[0] - kPath->Q[i+1].q[0]) + Sqr(kPath->Q[i].q[1] - kPath->Q[i+1].q[1]) + Sqr(kPath->Q[i].q[2] - kPath->Q[i+1].q[2]));
    }
  }
  if(NKPointsInPath!=BandEnergies->Ni)
  {
    printf("ERROR in ConvertBandEnergiesToBitMap, Number of points in kPath(%d) does not agree with number of energies in BandEnergies (%d)\n",NKPointsInPath,BandEnergies->Ni);
    return 1;
  }
  double PixelPerE = BandPlot->NpointsY / (Emax - Emin);
  double PixelPerZeroDistance = ZeroDistance * BandPlot->NpointsX;
  double PixelPerk = (BandPlot->NpointsX - NZeroDistance * PixelPerZeroDistance) / LengthKPath;
  double HalfPixelLineWidth = BandPlot->NpointsY * LineWidth / 2;
  double pixelkPosMin=0,pixelkPosMax=0;
  double dk;
  int ik=0;
  for(int i=0;i<(kPath->NQ)-1;i++)
  {
    if(kPath->Q[i].Nq==0)
    {
      pixelkPosMin+=PixelPerZeroDistance;
      ik++;
      continue;
    }
    dk = sqrt( Sqr(kPath->Q[i].q[0] - kPath->Q[i+1].q[0]) + Sqr(kPath->Q[i].q[1] - kPath->Q[i+1].q[1]) + Sqr(kPath->Q[i].q[2] - kPath->Q[i+1].q[2])) / (kPath->Q[i].Nq+1);
    for(int j=-1;j<=kPath->Q[i].Nq;j++)
    {
      if((j==-1)||(j==kPath->Q[i].Nq))
      {
        pixelkPosMax = pixelkPosMin + dk * PixelPerk / 2;
      }
      else
      {
        pixelkPosMax = pixelkPosMin + dk * PixelPerk;
      }
      for(int iE=0;iE<BandEnergies->Nj;iE++)
      {
        for(int ip=Max(0,floor((Emax-BandEnergies->VRe[ik*BandEnergies->Nj + iE])*PixelPerE-HalfPixelLineWidth)); ip<=Min(BandPlot->NpointsY-1,ceil((Emax-BandEnergies->VRe[ik*BandEnergies->Nj + iE])*PixelPerE+HalfPixelLineWidth)); ip++)
        {
          for(int jp=Max(0,floor(pixelkPosMin)); jp<=Min(BandPlot->NpointsX-1,ceil(pixelkPosMax)); jp++)
          {
            BandPlot->Pixels[ip*BandPlot->NpointsX + jp]=PixelColor;
          }
        }
      }
      if(j!=kPath->Q[i].Nq)
      {
        ik++;
      }
      pixelkPosMin = pixelkPosMax;
    }
  }
  return 0;
}
int ConvertBandEnergiesToBitMapInterpolate(CompactMatrixType* BandEnergies, BandstructurePathType *kPath, BitMapType *BandPlot, double Emin, double Emax, double LineWidth, uint32_t LineColor, double ZeroDistance)
{
 if(Emax < Emin) { printf("ERROR in ConvertBandEnergiesToBitMap, Emax < Emin\n"); return 1;}
 if(LineWidth<0) { printf("ERROR in ConvertBandEnergiesToBitMap, LineWidth < 0\n"); return 1;}
 if(ZeroDistance<0) { printf("ERROR in ConvertBandEnergiesToBitMap, ZeroDistance < 0\n"); return 1;}
 MakeBitMatWhiteTransparent(BandPlot);
 uint32_t const PixelColor = LineColor; // a pixel has the form of a colour, 0xAARRGGBB
 double LengthKPath = 0;
 int    NKPointsInPath=kPath->NQ;
 int    NZeroDistance = 0;
 for(int i=0;i<(kPath->NQ)-1;i++)
 {
   NKPointsInPath += kPath->Q[i].Nq;
   if(kPath->Q[i].Nq==0)
   {
     NZeroDistance++;
   }
   else
   {
     LengthKPath += sqrt( Sqr(kPath->Q[i].q[0] - kPath->Q[i+1].q[0]) + Sqr(kPath->Q[i].q[1] - kPath->Q[i+1].q[1]) + Sqr(kPath->Q[i].q[2] - kPath->Q[i+1].q[2]));
   }
 }
 if(NKPointsInPath!=BandEnergies->Ni)
 {
   printf("ERROR in ConvertBandEnergiesToBitMap, Number of points in kPath(%d) does not agree with number of energies in BandEnergies (%d)\n",NKPointsInPath,BandEnergies->Ni);
   return 1;
 }
 double PixelPerE = BandPlot->NpointsY / (Emax - Emin);
 double PixelPerZeroDistance = ZeroDistance * BandPlot->NpointsX;
 double PixelPerk = (BandPlot->NpointsX - NZeroDistance * PixelPerZeroDistance) / LengthKPath;
 double HalfPixelLineWidth = BandPlot->NpointsY * LineWidth / 2;
 double pixelkPosMin=0,pixelkPosMax=0,pixelPosCentre;
 double dk, dEleft, dEright;
 int ik=0;
 for(int i=0;i<(kPath->NQ)-1;i++)
 {
   if(kPath->Q[i].Nq==0)
   {
     pixelkPosMin+=PixelPerZeroDistance;
     ik++;
     continue;
   }
   dk = sqrt( Sqr(kPath->Q[i].q[0] - kPath->Q[i+1].q[0]) + Sqr(kPath->Q[i].q[1] - kPath->Q[i+1].q[1]) + Sqr(kPath->Q[i].q[2] - kPath->Q[i+1].q[2])) / (kPath->Q[i].Nq+1);
   for(int j=-1;j<=kPath->Q[i].Nq;j++)
   {
     if(j==-1)
     {
       pixelPosCentre = pixelkPosMin;
       pixelkPosMax = pixelkPosMin + dk * PixelPerk / 2;
     }
     else if (j==kPath->Q[i].Nq)
     {
       pixelkPosMax = pixelkPosMin + dk * PixelPerk / 2;
       pixelPosCentre = pixelkPosMax;
     }
     else
     {
       pixelPosCentre = pixelkPosMin + dk * PixelPerk / 2;
       pixelkPosMax   = pixelkPosMin + dk * PixelPerk;
     }
     for(int iE=0;iE<BandEnergies->Nj;iE++)
     {
       if(j==-1)
       {
         dEleft = 0;
       }
       else
       {
         dEleft = (BandEnergies->VRe[ik*BandEnergies->Nj + iE] - BandEnergies->VRe[(ik-1)*BandEnergies->Nj + iE])/(dk*PixelPerk);
       }
       if(j==kPath->Q[i].Nq)
       {
         dEright = 0;
       }
       else
       {
         dEright = (BandEnergies->VRe[ik*BandEnergies->Nj + iE] - BandEnergies->VRe[(ik+1)*BandEnergies->Nj + iE])/(dk*PixelPerk);
       }
       for(int jp=Max(0,floor(pixelkPosMin)); jp<=Min(BandPlot->NpointsX-1,ceil(pixelkPosMax)); jp++)
       {
         if(jp<pixelPosCentre)
         {
           for(int ip=Max(0,floor((Emax-BandEnergies->VRe[ik*BandEnergies->Nj + iE] + (pixelPosCentre-jp)*dEleft)*PixelPerE-HalfPixelLineWidth)); ip<=Min(BandPlot->NpointsY-1,ceil((Emax-BandEnergies->VRe[ik*BandEnergies->Nj + iE] + (pixelPosCentre-jp)*dEleft)*PixelPerE+HalfPixelLineWidth)); ip++)
           {
             BandPlot->Pixels[ip*BandPlot->NpointsX + jp]=PixelColor;
           }
         }
         else
         {
           for(int ip=Max(0,floor((Emax-BandEnergies->VRe[ik*BandEnergies->Nj + iE] - (pixelPosCentre-jp)*dEright)*PixelPerE-HalfPixelLineWidth)); ip<=Min(BandPlot->NpointsY-1,ceil((Emax-BandEnergies->VRe[ik*BandEnergies->Nj + iE] - (pixelPosCentre-jp)*dEright)*PixelPerE+HalfPixelLineWidth)); ip++)
           {
             BandPlot->Pixels[ip*BandPlot->NpointsX + jp]=PixelColor;
           }
         }
       }
     }
     if(j!=kPath->Q[i].Nq)
     {
       ik++;
     }
     pixelkPosMin = pixelkPosMax;
   }
 }
 return 0;
}

/*  Given a matrix with the band Energies generates a plot of the band-structure as bitmap
 *  Input
 *    BandEnergies     - a matrix with Nkpoint rows(Ni) and NFermion columns(Nj)
 *    BandCharacter    - a matrix with Nkpoint rows(Ni) and NFermion columns(Nj)
 *    kPath            - a path in the Brillouine zone with Nkpoint points
 *    BandPlot         - an allocated bitmap with the dimensions set
 *      PointType LB;              // left bottom position in fraction of the canvas
 *      PointType RT;              // right top position in fraction of the canvas
 *      int NpointsX;              // number of pixels in the x direction
 *      int NpointsY;              // number of pixels in the y direction
 *      uint32_t *Pixels;          // pixels in format alphaRGB of size NpointsX * NpointsY
 *    Emin             - double for the lower bound of the plot
 *    Emax             - doulbe for the upper bound of the plot
 *    LineWidth        - double line width in ratio of the bitmap height, rounded to integer pixels
 *    LineColor        - RGB alpha color used for the lines
 *    SecondLineColor  - RGB alpha color used for the lines
 *         color is character*LineColor + (1-character)*SecondLineColor
 *    ZeroDistance     - If there are zero k-points between two adjacent points in the kPath we add an aditional distance in ratio of the bitmap width, rounded to integer pixels
 */
int ConvertBandEnergiesToBitMapBW(CompactMatrixType* BandEnergies, CompactMatrixType* BandCharacter, BandstructurePathType *kPath, BitMapType *BandPlot, double Emin, double Emax, double LineWidth, uint32_t LineColor, uint32_t SecondLineColor, double ZeroDistance)
{
  if(Emax < Emin) { printf("ERROR in ConvertBandEnergiesToBitMap, Emax < Emin\n"); return 1;}
  if(LineWidth<0) { printf("ERROR in ConvertBandEnergiesToBitMap, LineWidth < 0\n"); return 1;}
  if(ZeroDistance<0) { printf("ERROR in ConvertBandEnergiesToBitMap, ZeroDistance < 0\n"); return 1;}
  if(BandEnergies->Ni != BandCharacter->Ni) { printf("ERROR in ConvertBandEnergiesToBitMap, Character and Energy matrix not the same dimension\n"); return 1;}
  if(BandEnergies->Nj != BandCharacter->Nj) { printf("ERROR in ConvertBandEnergiesToBitMap, Character and Energy matrix not the same dimension\n"); return 1;}
  MakeBitMatWhiteTransparent(BandPlot);
  uint32_t PixelColor;
  double col;
  double LengthKPath = 0;
  int    NKPointsInPath=kPath->NQ;
  int    NZeroDistance = 0;
  for(int i=0;i<(kPath->NQ)-1;i++)
  {
    NKPointsInPath += kPath->Q[i].Nq;
    if(kPath->Q[i].Nq==0)
    {
      NZeroDistance++;
    }
    else
    {
      LengthKPath += sqrt( Sqr(kPath->Q[i].q[0] - kPath->Q[i+1].q[0]) + Sqr(kPath->Q[i].q[1] - kPath->Q[i+1].q[1]) + Sqr(kPath->Q[i].q[2] - kPath->Q[i+1].q[2]));
    }
  }
  if(NKPointsInPath!=BandEnergies->Ni)
  {
    printf("ERROR in ConvertBandEnergiesToBitMap, Number of points in kPath(%d) does not agree with number of energies in BandEnergies (%d)\n",NKPointsInPath,BandEnergies->Ni);
    return 1;
  }
  double PixelPerE = BandPlot->NpointsY / (Emax - Emin);
  double PixelPerZeroDistance = ZeroDistance * BandPlot->NpointsX;
  double PixelPerk = (BandPlot->NpointsX - NZeroDistance * PixelPerZeroDistance) / LengthKPath;
  double HalfPixelLineWidth = BandPlot->NpointsY * LineWidth / 2;
  double pixelkPosMin=0,pixelkPosMax=0;
  double dk;
  int ik=0;
  for(int i=0;i<(kPath->NQ)-1;i++)
  {
    if(kPath->Q[i].Nq==0)
    {
      pixelkPosMin+=PixelPerZeroDistance;
      ik++;
      continue;
    }
    dk = sqrt( Sqr(kPath->Q[i].q[0] - kPath->Q[i+1].q[0]) + Sqr(kPath->Q[i].q[1] - kPath->Q[i+1].q[1]) + Sqr(kPath->Q[i].q[2] - kPath->Q[i+1].q[2])) / (kPath->Q[i].Nq+1);
    for(int j=-1;j<=kPath->Q[i].Nq;j++)
    {
      if((j==-1)||(j==kPath->Q[i].Nq))
      {
        pixelkPosMax = pixelkPosMin + dk * PixelPerk / 2;
      }
      else
      {
        pixelkPosMax = pixelkPosMin + dk * PixelPerk;
      }
      for(int iE=0;iE<BandEnergies->Nj;iE++)
      {
        // alpha channel
        col = BandCharacter->VRe[ik*BandCharacter->Nj + iE] * (255-ColorTransparency(LineColor)) + (1-BandCharacter->VRe[ik*BandCharacter->Nj + iE]) * (255-ColorTransparency(SecondLineColor)); // the opacity
        if(col<0) col=0;
        if(col>255) col=255;
        PixelColor=255-(uint32_t)col; // the transparency
        //Red
        col = BandCharacter->VRe[ik*BandCharacter->Nj + iE] * ColorRed(LineColor) + (1-BandCharacter->VRe[ik*BandCharacter->Nj + iE]) * ColorRed(SecondLineColor);
        if(col<0) col=0;
        if(col>255) col=255;
        PixelColor = (PixelColor<<8) + (uint32_t)col;
        //Green
        col = BandCharacter->VRe[ik*BandCharacter->Nj + iE] * ColorGreen(LineColor) + (1-BandCharacter->VRe[ik*BandCharacter->Nj + iE]) * ColorGreen(SecondLineColor);
        if(col<0) col=0;
        if(col>255) col=255;
        PixelColor = (PixelColor<<8) + (uint32_t)col;
        // Blue
        col = BandCharacter->VRe[ik*BandCharacter->Nj + iE] * ColorBlue(LineColor) + (1-BandCharacter->VRe[ik*BandCharacter->Nj + iE]) * ColorBlue(SecondLineColor);
        if(col<0) col=0;
        if(col>255) col=255;
        PixelColor = (PixelColor<<8) + (uint32_t)col;
        for(int ip=Max(0,floor((Emax-BandEnergies->VRe[ik*BandEnergies->Nj + iE])*PixelPerE-HalfPixelLineWidth)); ip<=Min(BandPlot->NpointsY-1,ceil((Emax-BandEnergies->VRe[ik*BandEnergies->Nj + iE])*PixelPerE+HalfPixelLineWidth)); ip++)
        {
          for(int jp=Max(0,floor(pixelkPosMin)); jp<=Min(BandPlot->NpointsX-1,ceil(pixelkPosMax)); jp++)
          {
            BandPlot->Pixels[ip*BandPlot->NpointsX + jp]=PixelColor;
          }
        }
      }
      if(j!=kPath->Q[i].Nq)
      {
        ik++;
      }
      pixelkPosMin = pixelkPosMax;
    }
  }
  return 0;
}
int ConvertBandEnergiesToBitMapBWInterpolate(CompactMatrixType* BandEnergies, CompactMatrixType* BandCharacter, BandstructurePathType *kPath, BitMapType *BandPlot, double Emin, double Emax, double LineWidth, uint32_t LineColor, uint32_t SecondLineColor, double ZeroDistance)
{
  if(Emax < Emin) { printf("ERROR in ConvertBandEnergiesToBitMap, Emax < Emin\n"); return 1;}
  if(LineWidth<0) { printf("ERROR in ConvertBandEnergiesToBitMap, LineWidth < 0\n"); return 1;}
  if(ZeroDistance<0) { printf("ERROR in ConvertBandEnergiesToBitMap, ZeroDistance < 0\n"); return 1;}
  if(BandEnergies->Ni != BandCharacter->Ni) { printf("ERROR in ConvertBandEnergiesToBitMap, Character and Energy matrix not the same dimension\n"); return 1;}
  if(BandEnergies->Nj != BandCharacter->Nj) { printf("ERROR in ConvertBandEnergiesToBitMap, Character and Energy matrix not the same dimension\n"); return 1;}
  MakeBitMatWhiteTransparent(BandPlot);
  uint32_t PixelColor;
  double col;
  double LengthKPath = 0;
  int    NKPointsInPath=kPath->NQ;
  int    NZeroDistance = 0;
  for(int i=0;i<(kPath->NQ)-1;i++)
  {
    NKPointsInPath += kPath->Q[i].Nq;
    if(kPath->Q[i].Nq==0)
    {
      NZeroDistance++;
    }
    else
    {
      LengthKPath += sqrt( Sqr(kPath->Q[i].q[0] - kPath->Q[i+1].q[0]) + Sqr(kPath->Q[i].q[1] - kPath->Q[i+1].q[1]) + Sqr(kPath->Q[i].q[2] - kPath->Q[i+1].q[2]));
    }
  }
  if(NKPointsInPath!=BandEnergies->Ni)
  {
    printf("ERROR in ConvertBandEnergiesToBitMap, Number of points in kPath(%d) does not agree with number of energies in BandEnergies (%d)\n",NKPointsInPath,BandEnergies->Ni);
    return 1;
  }
  double PixelPerE = BandPlot->NpointsY / (Emax - Emin);
  double PixelPerZeroDistance = ZeroDistance * BandPlot->NpointsX;
  double PixelPerk = (BandPlot->NpointsX - NZeroDistance * PixelPerZeroDistance) / LengthKPath;
  double HalfPixelLineWidth = BandPlot->NpointsY * LineWidth / 2;
  double pixelkPosMin=0,pixelkPosMax=0,pixelPosCentre;
  double dk, dEleft, dEright;
  int ik=0;
  for(int i=0;i<(kPath->NQ)-1;i++)
  {
    if(kPath->Q[i].Nq==0)
    {
      pixelkPosMin+=PixelPerZeroDistance;
      ik++;
      continue;
    }
    dk = sqrt( Sqr(kPath->Q[i].q[0] - kPath->Q[i+1].q[0]) + Sqr(kPath->Q[i].q[1] - kPath->Q[i+1].q[1]) + Sqr(kPath->Q[i].q[2] - kPath->Q[i+1].q[2])) / (kPath->Q[i].Nq+1);
    for(int j=-1;j<=kPath->Q[i].Nq;j++)
    {
      if(j==-1)
      {
        pixelPosCentre = pixelkPosMin;
        pixelkPosMax = pixelkPosMin + dk * PixelPerk / 2;
      }
      else if (j==kPath->Q[i].Nq)
      {
        pixelkPosMax = pixelkPosMin + dk * PixelPerk / 2;
        pixelPosCentre = pixelkPosMax;
      }
      else
      {
        pixelPosCentre = pixelkPosMin + dk * PixelPerk / 2;
        pixelkPosMax   = pixelkPosMin + dk * PixelPerk;
      }
      for(int iE=0;iE<BandEnergies->Nj;iE++)
      {
        // alpha channel
        col = BandCharacter->VRe[ik*BandCharacter->Nj + iE] * (255-ColorTransparency(LineColor)) + (1-BandCharacter->VRe[ik*BandCharacter->Nj + iE]) * (255-ColorTransparency(SecondLineColor)); // the opacity
        if(col<0) col=0;
        if(col>255) col=255;
        PixelColor=255-(uint32_t)col; // the transparency
        //Red
        col = BandCharacter->VRe[ik*BandCharacter->Nj + iE] * ColorRed(LineColor) + (1-BandCharacter->VRe[ik*BandCharacter->Nj + iE]) * ColorRed(SecondLineColor);
        if(col<0) col=0;
        if(col>255) col=255;
        PixelColor = (PixelColor<<8) + (uint32_t)col;
        //Green
        col = BandCharacter->VRe[ik*BandCharacter->Nj + iE] * ColorGreen(LineColor) + (1-BandCharacter->VRe[ik*BandCharacter->Nj + iE]) * ColorGreen(SecondLineColor);
        if(col<0) col=0;
        if(col>255) col=255;
        PixelColor = (PixelColor<<8) + (uint32_t)col;
        // Blue
        col = BandCharacter->VRe[ik*BandCharacter->Nj + iE] * ColorBlue(LineColor) + (1-BandCharacter->VRe[ik*BandCharacter->Nj + iE]) * ColorBlue(SecondLineColor);
        if(col<0) col=0;
        if(col>255) col=255;
        PixelColor = (PixelColor<<8) + (uint32_t)col;
        if(j==-1)
        {
          dEleft = 0;
        }
        else
        {
          dEleft = (BandEnergies->VRe[ik*BandEnergies->Nj + iE] - BandEnergies->VRe[(ik-1)*BandEnergies->Nj + iE])/(dk*PixelPerk);
        }
        if(j==kPath->Q[i].Nq)
        {
          dEright = 0;
        }
        else
        {
          dEright = (BandEnergies->VRe[ik*BandEnergies->Nj + iE] - BandEnergies->VRe[(ik+1)*BandEnergies->Nj + iE])/(dk*PixelPerk);
        }
        for(int jp=Max(0,floor(pixelkPosMin)); jp<=Min(BandPlot->NpointsX-1,ceil(pixelkPosMax)); jp++)
        {
          if(jp<pixelPosCentre)
          {
            for(int ip=Max(0,floor((Emax-BandEnergies->VRe[ik*BandEnergies->Nj + iE] + (pixelPosCentre-jp)*dEleft)*PixelPerE-HalfPixelLineWidth)); ip<=Min(BandPlot->NpointsY-1,ceil((Emax-BandEnergies->VRe[ik*BandEnergies->Nj + iE] + (pixelPosCentre-jp)*dEleft)*PixelPerE+HalfPixelLineWidth)); ip++)
            {
              BandPlot->Pixels[ip*BandPlot->NpointsX + jp]=PixelColor;
            }
          }
          else
          {
            for(int ip=Max(0,floor((Emax-BandEnergies->VRe[ik*BandEnergies->Nj + iE] - (pixelPosCentre-jp)*dEright)*PixelPerE-HalfPixelLineWidth)); ip<=Min(BandPlot->NpointsY-1,ceil((Emax-BandEnergies->VRe[ik*BandEnergies->Nj + iE] - (pixelPosCentre-jp)*dEright)*PixelPerE+HalfPixelLineWidth)); ip++)
            {
              BandPlot->Pixels[ip*BandPlot->NpointsX + jp]=PixelColor;
            }
          }
        }
      }
      if(j!=kPath->Q[i].Nq)
      {
        ik++;
      }
      pixelkPosMin = pixelkPosMax;
    }
  }
  return 0;
}

/*  Given a matrix with the band Energies generates a plot of the band-structure as bitmap
 *  Input
 *    BandEnergies     - a matrix with Nkpoint rows(Ni) and NFermion columns(Nj)
 *    BandCharacterR   - a matrix with Nkpoint rows(Ni) and NFermion columns(Nj)
 *    BandCharacterG   - a matrix with Nkpoint rows(Ni) and NFermion columns(Nj)
 *    BandCharacterB   - a matrix with Nkpoint rows(Ni) and NFermion columns(Nj)
 *    kPath            - a path in the Brillouine zone with Nkpoint points
 *    BandPlot         - an allocated bitmap with the dimensions set
 *      PointType LB;              // left bottom position in fraction of the canvas
 *      PointType RT;              // right top position in fraction of the canvas
 *      int NpointsX;              // number of pixels in the x direction
 *      int NpointsY;              // number of pixels in the y direction
 *      uint32_t *Pixels;          // pixels in format alphaRGB of size NpointsX * NpointsY
 *    Emin             - double for the lower bound of the plot
 *    Emax             - doulbe for the upper bound of the plot
 *    LineWidth        - double line width in ratio of the bitmap height, rounded to integer pixels
 *    ZeroDistance     - If there are zero k-points between two adjacent points in the kPath we add an aditional distance in ratio of the bitmap width, rounded to integer pixels
 */
int ConvertBandEnergiesToBitMapRGB(CompactMatrixType* BandEnergies, CompactMatrixType* BandCharacterR, CompactMatrixType* BandCharacterG, CompactMatrixType* BandCharacterB, BandstructurePathType *kPath, BitMapType *BandPlot, double Emin, double Emax, double LineWidth, double ZeroDistance)
{
  if(Emax < Emin) { printf("ERROR in ConvertBandEnergiesToBitMap, Emax < Emin\n"); return 1;}
  if(LineWidth<0) { printf("ERROR in ConvertBandEnergiesToBitMap, LineWidth < 0\n"); return 1;}
  if(ZeroDistance<0) { printf("ERROR in ConvertBandEnergiesToBitMap, ZeroDistance < 0\n"); return 1;}
  if(BandEnergies->Ni != BandCharacterR->Ni) { printf("ERROR in ConvertBandEnergiesToBitMap, Character and Energy matrix not the same dimension\n"); return 1;}
  if(BandEnergies->Nj != BandCharacterR->Nj) { printf("ERROR in ConvertBandEnergiesToBitMap, Character and Energy matrix not the same dimension\n"); return 1;}
  if(BandEnergies->Ni != BandCharacterG->Ni) { printf("ERROR in ConvertBandEnergiesToBitMap, Character and Energy matrix not the same dimension\n"); return 1;}
  if(BandEnergies->Nj != BandCharacterG->Nj) { printf("ERROR in ConvertBandEnergiesToBitMap, Character and Energy matrix not the same dimension\n"); return 1;}
  if(BandEnergies->Ni != BandCharacterB->Ni) { printf("ERROR in ConvertBandEnergiesToBitMap, Character and Energy matrix not the same dimension\n"); return 1;}
  if(BandEnergies->Nj != BandCharacterB->Nj) { printf("ERROR in ConvertBandEnergiesToBitMap, Character and Energy matrix not the same dimension\n"); return 1;}
  MakeBitMatWhiteTransparent(BandPlot);
  uint32_t PixelColor;
  double col;
  double LengthKPath = 0;
  int    NKPointsInPath=kPath->NQ;
  int    NZeroDistance = 0;
  for(int i=0;i<(kPath->NQ)-1;i++)
  {
    NKPointsInPath += kPath->Q[i].Nq;
    if(kPath->Q[i].Nq==0)
    {
      NZeroDistance++;
    }
    else
    {
      LengthKPath += sqrt( Sqr(kPath->Q[i].q[0] - kPath->Q[i+1].q[0]) + Sqr(kPath->Q[i].q[1] - kPath->Q[i+1].q[1]) + Sqr(kPath->Q[i].q[2] - kPath->Q[i+1].q[2]));
    }
  }
  if(NKPointsInPath!=BandEnergies->Ni)
  {
    printf("ERROR in ConvertBandEnergiesToBitMap, Number of points in kPath(%d) does not agree with number of energies in BandEnergies (%d)\n",NKPointsInPath,BandEnergies->Ni);
    return 1;
  }
  double PixelPerE = BandPlot->NpointsY / (Emax - Emin);
  double PixelPerZeroDistance = ZeroDistance * BandPlot->NpointsX;
  double PixelPerk = (BandPlot->NpointsX - NZeroDistance * PixelPerZeroDistance) / LengthKPath;
  double HalfPixelLineWidth = BandPlot->NpointsY * LineWidth / 2;
  double pixelkPosMin=0,pixelkPosMax=0;
  double dk;
  int ik=0;
  for(int i=0;i<(kPath->NQ)-1;i++)
  {
    if(kPath->Q[i].Nq==0)
    {
      pixelkPosMin+=PixelPerZeroDistance;
      ik++;
      continue;
    }
    dk = sqrt( Sqr(kPath->Q[i].q[0] - kPath->Q[i+1].q[0]) + Sqr(kPath->Q[i].q[1] - kPath->Q[i+1].q[1]) + Sqr(kPath->Q[i].q[2] - kPath->Q[i+1].q[2])) / (kPath->Q[i].Nq+1);
    for(int j=-1;j<=kPath->Q[i].Nq;j++)
    {
      if((j==-1)||(j==kPath->Q[i].Nq))
      {
        pixelkPosMax = pixelkPosMin + dk * PixelPerk / 2;
      }
      else
      {
        pixelkPosMax = pixelkPosMin + dk * PixelPerk;
      }
      for(int iE=0;iE<BandEnergies->Nj;iE++)
      {
        // the transparency, opaque
        PixelColor=0x00;
        //Red
        col = BandCharacterR->VRe[ik*BandCharacterR->Nj + iE] * 256;
        if(col<0) col=0;
        if(col>255) col=255;
        PixelColor = (PixelColor<<8) + (uint32_t)col;
        //Green
        col = BandCharacterG->VRe[ik*BandCharacterG->Nj + iE] * 256;
        if(col<0) col=0;
        if(col>255) col=255;
        PixelColor = (PixelColor<<8) + (uint32_t)col;
        // Blue
        col = BandCharacterB->VRe[ik*BandCharacterB->Nj + iE] * 256;
        if(col<0) col=0;
        if(col>255) col=255;
        PixelColor = (PixelColor<<8) + (uint32_t)col;
        for(int ip=Max(0,floor((Emax-BandEnergies->VRe[ik*BandEnergies->Nj + iE])*PixelPerE-HalfPixelLineWidth)); ip<=Min(BandPlot->NpointsY-1,ceil((Emax-BandEnergies->VRe[ik*BandEnergies->Nj + iE])*PixelPerE+HalfPixelLineWidth)); ip++)
        {
          for(int jp=Max(0,floor(pixelkPosMin)); jp<=Min(BandPlot->NpointsX-1,ceil(pixelkPosMax)); jp++)
          {
            BandPlot->Pixels[ip*BandPlot->NpointsX + jp]=PixelColor;
          }
        }
      }
      if(j!=kPath->Q[i].Nq)
      {
        ik++;
      }
      pixelkPosMin = pixelkPosMax;
    }
  }
  return 0;
}
int ConvertBandEnergiesToBitMapRGBInterpolate(CompactMatrixType* BandEnergies, CompactMatrixType* BandCharacterR, CompactMatrixType* BandCharacterG, CompactMatrixType* BandCharacterB, BandstructurePathType *kPath, BitMapType *BandPlot, double Emin, double Emax, double LineWidth, double ZeroDistance)
{
  if(Emax < Emin) { printf("ERROR in ConvertBandEnergiesToBitMap, Emax < Emin\n"); return 1;}
  if(LineWidth<0) { printf("ERROR in ConvertBandEnergiesToBitMap, LineWidth < 0\n"); return 1;}
  if(ZeroDistance<0) { printf("ERROR in ConvertBandEnergiesToBitMap, ZeroDistance < 0\n"); return 1;}
  if(BandEnergies->Ni != BandCharacterR->Ni) { printf("ERROR in ConvertBandEnergiesToBitMap, Character and Energy matrix not the same dimension\n"); return 1;}
  if(BandEnergies->Nj != BandCharacterR->Nj) { printf("ERROR in ConvertBandEnergiesToBitMap, Character and Energy matrix not the same dimension\n"); return 1;}
  if(BandEnergies->Ni != BandCharacterG->Ni) { printf("ERROR in ConvertBandEnergiesToBitMap, Character and Energy matrix not the same dimension\n"); return 1;}
  if(BandEnergies->Nj != BandCharacterG->Nj) { printf("ERROR in ConvertBandEnergiesToBitMap, Character and Energy matrix not the same dimension\n"); return 1;}
  if(BandEnergies->Ni != BandCharacterB->Ni) { printf("ERROR in ConvertBandEnergiesToBitMap, Character and Energy matrix not the same dimension\n"); return 1;}
  if(BandEnergies->Nj != BandCharacterB->Nj) { printf("ERROR in ConvertBandEnergiesToBitMap, Character and Energy matrix not the same dimension\n"); return 1;}
  MakeBitMatWhiteTransparent(BandPlot);
  uint32_t PixelColor;
  double col;
  double LengthKPath = 0;
  int    NKPointsInPath=kPath->NQ;
  int    NZeroDistance = 0;
  for(int i=0;i<(kPath->NQ)-1;i++)
  {
    NKPointsInPath += kPath->Q[i].Nq;
    if(kPath->Q[i].Nq==0)
    {
      NZeroDistance++;
    }
    else
    {
      LengthKPath += sqrt( Sqr(kPath->Q[i].q[0] - kPath->Q[i+1].q[0]) + Sqr(kPath->Q[i].q[1] - kPath->Q[i+1].q[1]) + Sqr(kPath->Q[i].q[2] - kPath->Q[i+1].q[2]));
    }
  }
  if(NKPointsInPath!=BandEnergies->Ni)
  {
    printf("ERROR in ConvertBandEnergiesToBitMap, Number of points in kPath(%d) does not agree with number of energies in BandEnergies (%d)\n",NKPointsInPath,BandEnergies->Ni);
    return 1;
  }
  double PixelPerE = BandPlot->NpointsY / (Emax - Emin);
  double PixelPerZeroDistance = ZeroDistance * BandPlot->NpointsX;
  double PixelPerk = (BandPlot->NpointsX - NZeroDistance * PixelPerZeroDistance) / LengthKPath;
  double HalfPixelLineWidth = BandPlot->NpointsY * LineWidth / 2;
  double pixelkPosMin=0,pixelkPosMax=0,pixelPosCentre;
  double dk, dEleft, dEright;
  int ik=0;
  for(int i=0;i<(kPath->NQ)-1;i++)
  {
    if(kPath->Q[i].Nq==0)
    {
      pixelkPosMin+=PixelPerZeroDistance;
      ik++;
      continue;
    }
    dk = sqrt( Sqr(kPath->Q[i].q[0] - kPath->Q[i+1].q[0]) + Sqr(kPath->Q[i].q[1] - kPath->Q[i+1].q[1]) + Sqr(kPath->Q[i].q[2] - kPath->Q[i+1].q[2])) / (kPath->Q[i].Nq+1);
    for(int j=-1;j<=kPath->Q[i].Nq;j++)
    {
      if(j==-1)
      {
        pixelPosCentre = pixelkPosMin;
        pixelkPosMax = pixelkPosMin + dk * PixelPerk / 2;
      }
      else if (j==kPath->Q[i].Nq)
      {
        pixelkPosMax = pixelkPosMin + dk * PixelPerk / 2;
        pixelPosCentre = pixelkPosMax;
      }
      else
      {
        pixelPosCentre = pixelkPosMin + dk * PixelPerk / 2;
        pixelkPosMax   = pixelkPosMin + dk * PixelPerk;
      }
      for(int iE=0;iE<BandEnergies->Nj;iE++)
      {
        // the transparency, opaque
        PixelColor=0x00;
        //Red
        col = BandCharacterR->VRe[ik*BandCharacterR->Nj + iE] * 256;
        if(col<0) col=0;
        if(col>255) col=255;
        PixelColor = (PixelColor<<8) + (uint32_t)col;
        //Green
        col = BandCharacterG->VRe[ik*BandCharacterG->Nj + iE] * 256;
        if(col<0) col=0;
        if(col>255) col=255;
        PixelColor = (PixelColor<<8) + (uint32_t)col;
        // Blue
        col = BandCharacterB->VRe[ik*BandCharacterB->Nj + iE] * 256;
        if(col<0) col=0;
        if(col>255) col=255;
        PixelColor = (PixelColor<<8) + (uint32_t)col;
        if(j==-1)
        {
          dEleft = 0;
        }
        else
        {
          dEleft = (BandEnergies->VRe[ik*BandEnergies->Nj + iE] - BandEnergies->VRe[(ik-1)*BandEnergies->Nj + iE])/(dk*PixelPerk);
        }
        if(j==kPath->Q[i].Nq)
        {
          dEright = 0;
        }
        else
        {
          dEright = (BandEnergies->VRe[ik*BandEnergies->Nj + iE] - BandEnergies->VRe[(ik+1)*BandEnergies->Nj + iE])/(dk*PixelPerk);
        }
        for(int jp=Max(0,floor(pixelkPosMin)); jp<=Min(BandPlot->NpointsX-1,ceil(pixelkPosMax)); jp++)
        {
          if(jp<pixelPosCentre)
          {
            for(int ip=Max(0,floor((Emax-BandEnergies->VRe[ik*BandEnergies->Nj + iE] + (pixelPosCentre-jp)*dEleft)*PixelPerE-HalfPixelLineWidth)); ip<=Min(BandPlot->NpointsY-1,ceil((Emax-BandEnergies->VRe[ik*BandEnergies->Nj + iE] + (pixelPosCentre-jp)*dEleft)*PixelPerE+HalfPixelLineWidth)); ip++)
            {
              BandPlot->Pixels[ip*BandPlot->NpointsX + jp]=PixelColor;
            }
          }
          else
          {
            for(int ip=Max(0,floor((Emax-BandEnergies->VRe[ik*BandEnergies->Nj + iE] - (pixelPosCentre-jp)*dEright)*PixelPerE-HalfPixelLineWidth)); ip<=Min(BandPlot->NpointsY-1,ceil((Emax-BandEnergies->VRe[ik*BandEnergies->Nj + iE] - (pixelPosCentre-jp)*dEright)*PixelPerE+HalfPixelLineWidth)); ip++)
            {
              BandPlot->Pixels[ip*BandPlot->NpointsX + jp]=PixelColor;
            }
          }
        }
      }
      if(j!=kPath->Q[i].Nq)
      {
        ik++;
      }
      pixelkPosMin = pixelkPosMax;
    }
  }
  return 0;
}

/*  Given band spectral functions on the pixel rows generates a plot of the band-structure as bitmap
 *  Input
 *    A                - NSet*Nk*BandPlot->NpointsY values, A[(s*Nk+ik)*NpointsY+ip] is set s at k-point ik
 *                       on pixel row ip (row 0 is the top of the plot, the highest energy).
 *                       Set 0 is the total spectral function, sets 1..NSet-1 the characters.
 *    NSet             - 1 (total only), 2 (total and one character) or 4 (total and three characters)
 *    Nk               - number of k-points
 *    kPath            - a path in the Brillouine zone with Nk points
 *    BandPlot         - an allocated bitmap with the dimensions set
 *    LineColor        - RGB alpha color used for NSet 1 and for character 1 if NSet is 2
 *    SecondLineColor  - RGB alpha color used for character 0 if NSet is 2
 *    Scale            - the opacity of a pixel is min(1, Scale * total)
 *    ZeroDistance     - If there are zero k-points between two adjacent points in the kPath we add an aditional distance in ratio of the bitmap width, rounded to integer pixels
 *  The color follows from the character c = A_c / A_total as for the band lines, character*LineColor +
 *  (1-character)*SecondLineColor for NSet 2 and (c1, c2, c3) as RGB for NSet 4. It is blended with white
 *  by the opacity and written opaque, pixels without weight stay white and transparent. Every k-point
 *  gets a strip of pixel columns, as in ConvertBandEnergiesToBitMap.
 */
int ConvertBandSpectraToBitMap(double const *A, int NSet, int Nk, BandstructurePathType *kPath, BitMapType *BandPlot, uint32_t LineColor, uint32_t SecondLineColor, double Scale, double ZeroDistance)
{
  if(ZeroDistance<0) { printf("ERROR in ConvertBandSpectraToBitMap, ZeroDistance < 0\n"); return 1;}
  if(!(NSet==1 || NSet==2 || NSet==4)) { printf("ERROR in ConvertBandSpectraToBitMap, NSet (%d) should be 1, 2 or 4\n", NSet); return 1;}
  MakeBitMatWhiteTransparent(BandPlot);
  double LengthKPath = 0;
  int    NKPointsInPath=kPath->NQ;
  int    NZeroDistance = 0;
  for(int i=0;i<(kPath->NQ)-1;i++)
  {
    NKPointsInPath += kPath->Q[i].Nq;
    if(kPath->Q[i].Nq==0)
    {
      NZeroDistance++;
    }
    else
    {
      LengthKPath += sqrt( Sqr(kPath->Q[i].q[0] - kPath->Q[i+1].q[0]) + Sqr(kPath->Q[i].q[1] - kPath->Q[i+1].q[1]) + Sqr(kPath->Q[i].q[2] - kPath->Q[i+1].q[2]));
    }
  }
  if(NKPointsInPath!=Nk)
  {
    printf("ERROR in ConvertBandSpectraToBitMap, Number of points in kPath(%d) does not agree with the number of spectra (%d)\n",NKPointsInPath,Nk);
    return 1;
  }
  int    const NY = BandPlot->NpointsY;
  double PixelPerZeroDistance = ZeroDistance * BandPlot->NpointsX;
  double PixelPerk = (BandPlot->NpointsX - NZeroDistance * PixelPerZeroDistance) / LengthKPath;
  double pixelkPosMin=0,pixelkPosMax=0;
  double dk, tot, alpha, c, col[4];
  uint32_t *Column = (uint32_t*)malloc(NY*sizeof(uint32_t));
  int ik=0;
  for(int i=0;i<(kPath->NQ)-1;i++)
  {
    if(kPath->Q[i].Nq==0)
    {
      pixelkPosMin+=PixelPerZeroDistance;
      ik++;
      continue;
    }
    dk = sqrt( Sqr(kPath->Q[i].q[0] - kPath->Q[i+1].q[0]) + Sqr(kPath->Q[i].q[1] - kPath->Q[i+1].q[1]) + Sqr(kPath->Q[i].q[2] - kPath->Q[i+1].q[2])) / (kPath->Q[i].Nq+1);
    for(int j=-1;j<=kPath->Q[i].Nq;j++)
    {
      if((j==-1)||(j==kPath->Q[i].Nq))
      {
        pixelkPosMax = pixelkPosMin + dk * PixelPerk / 2;
      }
      else
      {
        pixelkPosMax = pixelkPosMin + dk * PixelPerk;
      }
      // the colors of the pixel column of k-point ik
      for(int ip=0;ip<NY;ip++)
      {
        tot = A[ik*NY+ip];
        alpha = Scale * tot;
        if(!(alpha>0))
        {
          Column[ip] = 0xFFFFFFFF; // white, fully transparent
          continue;
        }
        if(alpha>1) alpha=1;
        if(NSet==1)
        {
          col[0] = ColorRed(LineColor);
          col[1] = ColorGreen(LineColor);
          col[2] = ColorBlue(LineColor);
          col[3] = 255 - ColorTransparency(LineColor); // the opacity
        }
        else if(NSet==2)
        {
          c = A[(Nk+ik)*NY+ip] / tot;
          if(c<0) c=0;
          if(c>1) c=1;
          col[0] = c * ColorRed(LineColor)   + (1-c) * ColorRed(SecondLineColor);
          col[1] = c * ColorGreen(LineColor) + (1-c) * ColorGreen(SecondLineColor);
          col[2] = c * ColorBlue(LineColor)  + (1-c) * ColorBlue(SecondLineColor);
          col[3] = c * (255 - ColorTransparency(LineColor)) + (1-c) * (255 - ColorTransparency(SecondLineColor)); // the opacity
        }
        else
        {
          for(int ch=0;ch<3;ch++)
          {
            col[ch] = 256 * A[((ch+1)*Nk+ik)*NY+ip] / tot;
            if(col[ch]<0) col[ch]=0;
            if(col[ch]>255) col[ch]=255;
          }
          col[3] = 255;
        }
        // blend with white by the opacity of the pixel times the opacity of the color; the pixel is opaque
        alpha *= col[3] / 255;
        Column[ip] = 0x00;
        for(int ch=0;ch<3;ch++)
        {
          c = 255 - alpha * (255 - col[ch]);
          if(c<0) c=0;
          if(c>255) c=255;
          Column[ip] = (Column[ip]<<8) + (uint32_t)c;
        }
      }
      for(int ip=0;ip<NY;ip++)
      {
        for(int jp=Max(0,floor(pixelkPosMin)); jp<=Min(BandPlot->NpointsX-1,ceil(pixelkPosMax)); jp++)
        {
          BandPlot->Pixels[ip*BandPlot->NpointsX + jp]=Column[ip];
        }
      }
      if(j!=kPath->Q[i].Nq)
      {
        ik++;
      }
      pixelkPosMin = pixelkPosMax;
    }
  }
  free(Column);
  return 0;
}
