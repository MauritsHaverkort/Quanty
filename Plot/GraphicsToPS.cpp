int AddPrologToPS(DynamicString *Buff, int X0, int Y0, int XMAX, int YMAX)
{
  sprintfAdd(Buff, "%%!PS-Adobe-2.0 EPSF-2.0\n");
  sprintfAdd(Buff, "%%%%BoundingBox: %d %d %d %d\n",X0,Y0,XMAX,YMAX);
  sprintfAdd(Buff, "%%%%Creator: Quanty (www.Quanty.org)\n");
  sprintfAdd(Buff, "%%%%EndComments\n");
  sprintfAdd(Buff, "%%%%BeginProlog\nsave\n");
  sprintfAdd(Buff, "countdictstack\n");
  sprintfAdd(Buff, "mark\n");
  sprintfAdd(Buff, "newpath\n");
  sprintfAdd(Buff, "/setpagedevice {pop} def\n");
  sprintfAdd(Buff, "/L {newpath setlinewidth setrgbcolor moveto} def\n");
  sprintfAdd(Buff, "/a {newpath setlinewidth setrgbcolor arc} def\n");
  sprintfAdd(Buff, "/l {lineto} def\n");
  sprintfAdd(Buff, "/c {curveto} def\n");
  sprintfAdd(Buff, "/rightshow { dup stringwidth pop neg 0 rmoveto show} def\n");
  sprintfAdd(Buff, "/centershow { dup stringwidth pop -2 div 0 rmoveto show} def\n");
  sprintfAdd(Buff, "/f {newpath moveto rotate findfont exch scalefont setfont setrgbcolor} def\n");
  sprintfAdd(Buff, "%%%%EndProlog\n%%%%Page 1 1\n");
  return 0;
}

int AddCloseToPS(DynamicString *Buff)
{
  sprintfAdd(Buff,"showpage\n");
  sprintfAdd(Buff,"%%%%Trailer\n");
  sprintfAdd(Buff,"cleartomark\n");
  sprintfAdd(Buff,"countdictstack\n");
  sprintfAdd(Buff,"exch sub { end } repeat\n");
  sprintfAdd(Buff,"restore\n");
  sprintfAdd(Buff,"%%EOF\n");
  return 0;
}

int AddLineToPS(DynamicString *Buff, LineType *Line, int X0, int Y0, double scaleX, double scaleY, double scaleV)
{
  int i;
  if(Line->NPoints<2) return 0;
  if((Line->CloseAndFill) && (Line->FillColor[3] !=0) )
  {
    sprintfAdd(Buff, "%d %d %.4f %.4f %.4f %d L\n",(int)round(X0 + scaleX * Line->Point[0].X),(int)round(Y0 + scaleY * Line->Point[0].Y),
                                                   (double)Line->FillColor[0]/255., (double)Line->FillColor[1]/255., (double)Line->FillColor[2]/255.,0);
    for(i=1;i<Line->NPoints;i++)
    {
      sprintfAdd(Buff, "%d %d l\n",(int)round(X0 + scaleX * Line->Point[i].X),(int)round(Y0 + scaleY * Line->Point[i].Y));
    }
    sprintfAdd(Buff, "closepath\n");
    sprintfAdd(Buff, "gsave\n");
    sprintfAdd(Buff, "fill\n");
    sprintfAdd(Buff, "grestore\n");
    sprintfAdd(Buff, "%.4f %.4f %.4f setrgbcolor\n", (double)Line->Color[0]/255., (double)Line->Color[1]/255., (double)Line->Color[2]/255. );
    sprintfAdd(Buff, "%d setlinewidth\n",Max(1,(int)round(Line->Width * scaleV) ));
    sprintfAdd(Buff, "stroke\n");
  }
  else
  {
    sprintfAdd(Buff, "%d %d %.4f %.4f %.4f %d L\n",(int)round(X0 + scaleX * Line->Point[0].X),(int)round(Y0 + scaleY * Line->Point[0].Y),
               (double)Line->Color[0]/255., (double)Line->Color[1]/255., (double)Line->Color[2]/255., Max(1,(int)round(Line->Width * scaleV) ) );
    for(i=1;i<Line->NPoints;i++)
    {
      sprintfAdd(Buff, "%d %d l\n",(int)round(X0 + scaleX * Line->Point[i].X),(int)round(Y0 + scaleY * Line->Point[i].Y));
    }
    if(Line->CloseAndFill)
    {
      sprintfAdd(Buff, "closepath\n");
    }
    sprintfAdd(Buff, "stroke\n");
  }

  return 0;
}


int AddCurveToPS(DynamicString *Buff, CurveType *Curve, int X0, int Y0, double scaleX, double scaleY, double scaleV)
{
  int i;
  if(Curve->NPoints<2) return 0;
  if((Curve->CloseAndFill) && (Curve->FillColor[3] !=0) )
  {
    sprintfAdd(Buff, "%d %d %.4f %.4f %.4f %d L\n",(int)round(X0 + scaleX * Curve->Point[0].X),(int)round(Y0 + scaleY * Curve->Point[0].Y),
                                             (double)Curve->FillColor[0]/255., (double)Curve->FillColor[1]/255., (double)Curve->FillColor[2]/255., 0 );
    for(i=1;i<Curve->NPoints;i++)
    {
      sprintfAdd(Buff, "%d %d %d %d %d %d c\n",(int)round(X0 + scaleX * Curve->VR[i-1].X),(int)round(Y0 + scaleY * Curve->VR[i-1].Y)
                                              ,(int)round(X0 + scaleX * Curve->VL[i  ].X),(int)round(Y0 + scaleY * Curve->VL[i  ].Y)
                                              ,(int)round(X0 + scaleX * Curve->Point[i].X),(int)round(Y0 + scaleY * Curve->Point[i].Y));
    }
    sprintfAdd(Buff, "%d %d %d %d %d %d c\n",(int)round(X0 + scaleX * Curve->VR[i-1].X),(int)round(Y0 + scaleY * Curve->VR[i-1].Y)
                                            ,(int)round(X0 + scaleX * Curve->VL[0].X),(int)round(Y0 + scaleY * Curve->VL[0].Y)
                                            ,(int)round(X0 + scaleX * Curve->Point[0].X),(int)round(Y0 + scaleY * Curve->Point[0].Y));
    sprintfAdd(Buff, "gsave\n");
    sprintfAdd(Buff, "closepath\n");
    sprintfAdd(Buff, "fill\n");
    sprintfAdd(Buff, "grestore\n");
    sprintfAdd(Buff, "%.4f %.4f %.4f setrgbcolor\n", (double)Curve->Color[0]/255., (double)Curve->Color[1]/255., (double)Curve->Color[2]/255. );
    sprintfAdd(Buff, "%d setlinewidth\n",Max(1,(int)round(Curve->Width * scaleV) ));
    sprintfAdd(Buff, "stroke\n");
  }
  else
  {
    sprintfAdd(Buff, "%d %d %.4f %.4f %.4f %d L\n",(int)round(X0 + scaleX * Curve->Point[0].X),(int)round(Y0 + scaleY * Curve->Point[0].Y),
               (double)Curve->Color[0]/255., (double)Curve->Color[1]/255., (double)Curve->Color[2]/255., Max(1,(int)round(Curve->Width * scaleV) ) );
    for(i=1;i<Curve->NPoints;i++)
    {
      sprintfAdd(Buff, "%d %d %d %d %d %d c\n",(int)round(X0 + scaleX * Curve->VR[i-1].X),(int)round(Y0 + scaleY * Curve->VR[i-1].Y)
                                              ,(int)round(X0 + scaleX * Curve->VL[i  ].X),(int)round(Y0 + scaleY * Curve->VL[i  ].Y)
                                              ,(int)round(X0 + scaleX * Curve->Point[i].X),(int)round(Y0 + scaleY * Curve->Point[i].Y));
    }
    if(Curve->CloseAndFill)
    {
      sprintfAdd(Buff, "%d %d %d %d %d %d c\n",(int)round(X0 + scaleX * Curve->VR[i-1].X),(int)round(Y0 + scaleY * Curve->VR[i-1].Y)
                                              ,(int)round(X0 + scaleX * Curve->VL[0].X),(int)round(Y0 + scaleY * Curve->VL[0].Y)
                                              ,(int)round(X0 + scaleX * Curve->Point[0].X),(int)round(Y0 + scaleY * Curve->Point[0].Y));
    }
    sprintfAdd(Buff, "stroke\n");
  }

  return 0;
}

int AddArcToPS(DynamicString *Buff, ArcType *Arc, int X0, int Y0, double scaleX, double scaleY, double scaleV)
{
  if((Arc->CloseAndFill)&&(Arc->FillColor[3] !=0))
  {
    sprintfAdd(Buff, "%d %d %d %d %d %.4f %.4f %.4f %d a\n", (int)round(X0 + scaleX * Arc->Point.X),(int)round(Y0 + scaleY * Arc->Point.Y)
                                                           , (int)round(scaleV * Arc->R), (int)round(Arc->ThetaR * 360/(2*Pi)), (int)round(Arc->ThetaL * 360/(2*Pi))
                                                           , ((double)Arc->FillColor[0])/255., ((double)Arc->FillColor[1])/255.,((double)Arc->FillColor[2])/255., 0 );
    if(Arc->ThetaL-Arc->ThetaR < 1.99 * Pi)
    {
      sprintfAdd(Buff, "%d %d l\n",(int)round(X0 + scaleX * Arc->Point.X),(int)round(Y0 + scaleY * Arc->Point.Y));
    }
    sprintfAdd(Buff, "closepath\n");
    sprintfAdd(Buff, "gsave\n");
    sprintfAdd(Buff, "fill\n");
    sprintfAdd(Buff, "grestore\n");
    sprintfAdd(Buff, "%.4f %.4f %.4f setrgbcolor\n", ((double)Arc->Color[0])/255., ((double)Arc->Color[1])/255., ((double)Arc->Color[2])/255. );
    sprintfAdd(Buff, "%d setlinewidth\n",Max(1,(int)round(Arc->Width * scaleV) ));
    sprintfAdd(Buff, "stroke\n");
    
  }
  else
  {
    sprintfAdd(Buff, "%d %d %d %d %d %.4f %.4f %.4f %d a\n", (int)round(X0 + scaleX * Arc->Point.X),(int)round(Y0 + scaleY * Arc->Point.Y)
                                                           , (int)round(scaleV * Arc->R), (int)round(Arc->ThetaR * 360/(2*Pi)), (int)round(Arc->ThetaL * 360/(2*Pi))
                                                           , ((double)Arc->Color[0])/255., ((double)Arc->Color[1])/255.,((double)Arc->Color[2])/255., Max(1,(int)round(Arc->Width * scaleV) ) );
    if(Arc->CloseAndFill)
    {
      if(Arc->ThetaL-Arc->ThetaR < 1.99 * Pi)
      {
        sprintfAdd(Buff, "%d %d l\n",(int)round(X0 + scaleX * Arc->Point.X),(int)round(Y0 + scaleY * Arc->Point.Y));
      }
      sprintfAdd(Buff, "closepath\n");
    }
    sprintfAdd(Buff, "stroke\n");
  }
  return 0;
}

int AddBitMapToPS(DynamicString *Buff, BitMapType *BitMap, int X0, int Y0, double scaleX, double scaleY, double scaleV)
{
  int i,j;
  double BscaleX, BscaleY, BX0, BY0;
  BscaleX = BitMap->NpointsX / ((BitMap->RT.X-BitMap->LB.X)*scaleX);
  BscaleY = BitMap->NpointsY / ((BitMap->RT.Y-BitMap->LB.Y)*scaleY);
  BX0 = ( -X0 - BitMap->LB.X * scaleX ) * BscaleX;
  BY0 = ( -Y0 - BitMap->LB.Y * scaleY ) * BscaleY;
  sprintfAdd(Buff,"gsave\n");
  sprintfAdd(Buff,"%d %d 8\n",BitMap->NpointsX,BitMap->NpointsY);
  sprintfAdd(Buff,"[%.6f 0 0 %.6f %.6f %.6f]\n",BscaleX, BscaleY, BX0, BY0);
  sprintfAdd(Buff,"currentfile /ASCIIHexDecode filter\n");
  sprintfAdd(Buff,"false 3 colorimage\n");
  for(j=BitMap->NpointsY-1;j>=0;j--)
  {
    for(i=0;i<BitMap->NpointsX;i++)
    {
      // We store from left to right the alpha channel, the red channel, the green channel and the blue channel.
      // the alpha channel is removed as postscript does not handle transparancy.
      sprintfAdd(Buff,"%06X ",BitMap->Pixels[i+j*BitMap->NpointsX]%0x1000000);
    }
    sprintfAdd(Buff,"\n");
  }
  sprintfAdd(Buff,">\n");
  sprintfAdd(Buff,"grestore\n");
  return 0;
}

int AddTextToPS(DynamicString *Buff, TextType *Text, int X0, int Y0, double scaleX, double scaleY, double scaleV)
{
  if(Text->Color[3]==0) return 0;
  sprintfAdd(Buff, "gsave\n");
  sprintfAdd(Buff, "%.4f %.4f %.4f %d /%s %d %d %d f\n", ((double)Text->Color[0])/255., ((double)Text->Color[1])/255., ((double)Text->Color[2])/255.,
                                                         (int)round(scaleV * Text->Size), Text->Font, (int)round(Text->Angle*180/Pi),
                                                         (int)round(X0 + scaleX * Text->Point.X),(int)round(Y0 + scaleY * Text->Point.Y) );
  sprintfAdd(Buff,"(%s) ",Text->String);
  switch(Text->Alignment)
  {
    case 'l':
      sprintfAdd(Buff,"show\n");
    break;
    case 'r':
      sprintfAdd(Buff,"rightshow\n");
    break;
    case 'c':
      sprintfAdd(Buff,"centershow\n");
    break;
    default:
      printf("ERROR: unrecognized allignment token for text in graph to postscript\n");
  }
  sprintfAdd(Buff, "grestore\n");
  return 0;
}


int AddLayersToPS(DynamicString *Buff, GraphicsType *G)
{
  int i;
  
  int XMax,YMax,X0,Y0;
  double const scaleX =  POSTSCIPT_PIXELS_PER_MM * G->CanvasSize.X ;
  double const scaleY =  POSTSCIPT_PIXELS_PER_MM * G->CanvasSize.Y ;
  double const scaleV = sqrt(abs(scaleX * scaleY));
  
  X0 = round(POSTSCIPT_PIXELS_PER_MM * G->CanvasOrigin.X);
  Y0 = round(POSTSCIPT_PIXELS_PER_MM * G->CanvasOrigin.Y);
  XMax = round(POSTSCIPT_PIXELS_PER_MM * ( G->CanvasSize.X + G->CanvasOrigin.X ) );
  YMax = round(POSTSCIPT_PIXELS_PER_MM * ( G->CanvasSize.Y + G->CanvasOrigin.Y ) );

  
  for(i=0;i<G->NBitMap;i++)
  {
    AddBitMapToPS(Buff, &(G->BitMap[i]), X0, Y0, scaleX, scaleY, scaleV);
  }
  for(i=0;i<G->NLine;i++)
  {
    AddLineToPS(Buff, &(G->Line[i]), X0, Y0, scaleX, scaleY, scaleV);
  }
  for(i=0;i<G->NCurve;i++)
  {
    AddCurveToPS(Buff, &(G->Curve[i]), X0, Y0, scaleX, scaleY, scaleV);
  }
  for(i=0;i<G->NArc;i++)
  {
    AddArcToPS(Buff, &(G->Arc[i]), X0, Y0, scaleX, scaleY, scaleV);
  }
  for(i=0;i<G->NText;i++)
  {
    AddTextToPS(Buff, &(G->Text[i]), X0, Y0, scaleX, scaleY, scaleV);
  }
  for(i=1;i<G->NLayer;i++)
  {
    AddLayersToPS(Buff, G->Layer[i]);
  }
  return 0;
}




// G is a graphics type
// *PS is a string that will contain the ps text

int GraphicsToPS(GraphicsType *G, DynamicString *PS)
{
  int XMax,YMax,X0,Y0;
  
  X0 = round(POSTSCIPT_PIXELS_PER_MM * G->CanvasOrigin.X);
  Y0 = round(POSTSCIPT_PIXELS_PER_MM * G->CanvasOrigin.Y);
  XMax = round(POSTSCIPT_PIXELS_PER_MM * ( G->CanvasSize.X + G->CanvasOrigin.X ) );
  YMax = round(POSTSCIPT_PIXELS_PER_MM * ( G->CanvasSize.Y + G->CanvasOrigin.Y ) );
  
  *PS = InitDynamicString();
  
  AddPrologToPS(PS, X0, Y0, XMax, YMax);

  AddLayersToPS(PS, G);

  AddCloseToPS(PS);
  return 0;
}




int GraphicsToPS(GraphicsType *G, char** PS)
{
  DynamicString Buff;
  GraphicsToPS(G, &Buff);
  *PS = Buff.str;
  return 0;
}







