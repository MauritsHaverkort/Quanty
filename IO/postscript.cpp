/*  allocates a new string and adds the Header of a post-script file with bounding box from 0 0 to XSize YSize
 *
 */
int CreatePS(char ** res, int XSize, int YSize)
{
  *res = (char*)malloc(2024*sizeof(char));
  if(res==NULL) return 1;
  res[0][0]='\0';
  sprintfAdd(res, "%%!PS-Adobe-2.0 EPSF-2.0\n");
  sprintfAdd(res, "%%%%BoundingBox: 0 0 %5d %5d\n",XSize,YSize);
  sprintfAdd(res, "%%%%Creator: Quanty (www.Quanty.org)\n");
  sprintfAdd(res, "%%%%EndComments\n");
  sprintfAdd(res, "%%%%BeginProlog\nsave\n");
  sprintfAdd(res, "countdictstack\n");
  sprintfAdd(res, "mark\n");
  sprintfAdd(res, "newpath\n");
  sprintfAdd(res, "/setpagedevice {pop} def\n");
  sprintfAdd(res, "/p {newpath moveto} def\n");
  sprintfAdd(res, "/m {moveto} def\n");
  sprintfAdd(res, "/L {newpath setlinewidth setrgbcolor moveto} def\n");
  sprintfAdd(res, "/l {lineto} def\n");
  sprintfAdd(res, "/c {closepath} def\n");
  sprintfAdd(res, "/s {stroke} def\n");
  sprintfAdd(res, "/rightshow { dup stringwidth pop neg 0 rmoveto show} def\n");
  sprintfAdd(res, "/centershow { dup stringwidth pop -2 div 0 rmoveto show} def\n");
  sprintfAdd(res, "/lp {newpath moveto /Times-Roman findfont exch scalefont setfont 0 0 0 setrgbcolor show} def\n");        // (text) FontSize PositionX PositionY
  sprintfAdd(res, "/rp {newpath moveto /Times-Roman findfont exch scalefont setfont 0 0 0 setrgbcolor rightshow} def\n");
  sprintfAdd(res, "/cp {newpath moveto /Times-Roman findfont exch scalefont setfont 0 0 0 setrgbcolor centershow} def\n");
  sprintfAdd(res, "/lgp {newpath moveto /Symbol findfont exch scalefont setfont 0 0 0 setrgbcolor show} def\n");        // (text) FontSize PositionX PositionY
  sprintfAdd(res, "/rgp {newpath moveto /Symbol findfont exch scalefont setfont 0 0 0 setrgbcolor rightshow} def\n");
  sprintfAdd(res, "/cgp {newpath moveto /Symbol findfont exch scalefont setfont 0 0 0 setrgbcolor centershow} def\n");
  sprintfAdd(res, "%%%%EndProlog\n%%%%Page 1 1\n");
  return 0;
}

int ClosePS(char ** res)
{
  sprintfAdd(res,"showpage\n");
  sprintfAdd(res,"%%%%Trailer\n");
  sprintfAdd(res,"cleartomark\n");
  sprintfAdd(res,"countdictstack\n");
  sprintfAdd(res,"exch sub { end } repeat\n");
  sprintfAdd(res,"restore\n");
  sprintfAdd(res,"%%EOF\n");
  return 0;
}

/*
 *
 */
int AddTextPS(char ** res, const char *text, double x, double y, double size, char align)
{
  switch (align) {
    case 'c':
      sprintfAdd(res,"(%s) %f %f %f cp\n",text,size,x,y);
      break;
    case 'l':
      sprintfAdd(res,"(%s) %f %f %f lp\n",text,size,x,y);
      break;
    case 'r':
      sprintfAdd(res,"(%s) %f %f %f rp\n",text,size,x,y);
      break;
    default:
      sprintfAdd(res,"(Undefined align character \"%c\" for text \"%s\") %f %f %f lp\n",align,text,size,x,y);
      break;
  }
  return 0;
}

/*
 *
 */
int AddGreekTextPS(char ** res, const char *text, double x, double y, double size, char align)
{
  switch (align) {
    case 'c':
      sprintfAdd(res,"(%s) %f %f %f cgp\n",text,size,x,y);
      break;
    case 'l':
      sprintfAdd(res,"(%s) %f %f %f lgp\n",text,size,x,y);
      break;
    case 'r':
      sprintfAdd(res,"(%s) %f %f %f rgp\n",text,size,x,y);
      break;
    default:
      sprintfAdd(res,"(Undefined align character \"%c\" for text \"%s\") %f %f %f lp\n",align,text,size,x,y);
      break;
  }
  return 0;
}

/*
 *
 */
int AddLinePathPS(char ** res, double *x, double *y, int N, uint32_t RGBColor, double thickness)
{
  int i;
  if(N<2) return 0;
  sprintfAdd(res, "%f %f %f %f %f %f L\n",x[0],y[0],((double)Red(RGBColor))/255,((double)Green(RGBColor))/255,((double)Blue(RGBColor))/255,thickness);
  for(i=1;i<N;i++)
  {
    sprintfAdd(res, "%f %f l\n",x[i],y[i]);
  }
  sprintfAdd(res, "stroke\n");
  return 0;
}

/*
 *
 */
int AddPolygonPS(char ** res, double *x, double *y, int N, uint32_t RGBColor, double thickness)
{
  int i;
  if(N<2) return 0;
  sprintfAdd(res, "%f %f %f %f %f %f L\n",x[0],y[0],((double)(510+Red(RGBColor)))/765,((double)(510+Green(RGBColor)))/765,((double)(510+Blue(RGBColor)))/765,0);
  for(i=1;i<N;i++)
  {
    sprintfAdd(res, "%f %f l\n",x[i],y[i]);
  }
  sprintfAdd(res, "gsave\n");
  sprintfAdd(res, "closepath\n");
  sprintfAdd(res, "fill\n");
  sprintfAdd(res, "grestore\n");
  sprintfAdd(res, "%f %f %f setrgbcolor\n",((double)Red(RGBColor))/255,((double)Green(RGBColor))/255,((double)Blue(RGBColor))/255);
  sprintfAdd(res, "%f setlinewidth\n",thickness);
  sprintfAdd(res, "stroke\n");
  return 0;
}


/*  Place a bitmap with Nwidth and Nhight pixels on the rectangle from x to x+width and y to y+hight
 *  Use the array RGBColor with Nwidth*Nhight elements for the color
 */
int AddBitmapPS(char ** res, double x, double y, double width, double heigth, uint32_t *RGBColor, uint32_t Nwidth, uint32_t Nheigth)
{
  int i;
  sprintfAdd(res, "%f %f translate\n",x,y);
  sprintfAdd(res, "gsave\n");
  sprintfAdd(res, "%f %f scale\n",x+width,y+heigth);
  sprintfAdd(res, "%d %d 8\n",Nwidth,Nheigth);
  sprintfAdd(res, "[%d 0 0 -%d 0 %d]\n",Nwidth,Nheigth,Nheigth);
  sprintfAdd(res, "{ currentfile %d string readhexstring pop }\n",3*Nwidth*Nheigth);
  sprintfAdd(res, "false\n");
  sprintfAdd(res, "3 colorimage\n");
  for(i=0;i<Nwidth*Nheigth;i++)
  {
    sprintfAdd(res, "%06X",RGBColor[i]);
    if(i&32)sprintfAdd(res, "\n");
  }
  sprintfAdd(res, "\ngrestore\n");
  return 0;
}
