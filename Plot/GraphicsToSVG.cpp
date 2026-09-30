int AddPrologToSVG(DynamicString *Buff, int X0, int Y0, int XMAX, int YMAX)
{
  sprintfAdd(Buff, "<svg width=\"%d\" height=\"%d\" viewBox=\"%d %d %d %d\" preserveAspectRatio=\"none\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\"  >\n",XMAX-X0,YMAX-Y0,X0,0,XMAX-X0,YMAX-Y0);
  return 0;
}

int AddRelativeSizePrologToSVG(DynamicString *Buff, int X0, int Y0, int XMAX, int YMAX)
{
  sprintfAdd(Buff, "<svg width=\"100%%\" height=\"100%%\" viewBox=\"%d %d %d %d\" preserveAspectRatio=\"none\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\"  >\n",X0,0,XMAX-X0,YMAX-Y0);
  return 0;
}

int AddCloseToSVG(DynamicString *Buff)
{
  sprintfAdd(Buff,"</svg>\n");
  return 0;
}

int AddLineToSVG(DynamicString *Buff, LineType *Line, int X0, int Y0, double scaleX, double scaleY, double scaleV)
{
  int i;
  if(Line->NPoints<2) return 0;
  sprintfAdd(Buff, "<path ");
  if((Line->CloseAndFill) && (ColorVisible(Line->FillColor)) )
  {
    if( ColorTransparency(Line->FillColor)==ColorTransparency(Line->Color) )
    {
      sprintfAdd(Buff, "opacity=\"%.5f\" ",ColorOpacity(Line->FillColor));
    }
    else
    {
      sprintfAdd(Buff, "fill-opacity=\"%.5f\" ",ColorOpacity(Line->FillColor));
      sprintfAdd(Buff, "stroke-opacity=\"%.5f\" ",ColorOpacity(Line->Color));
    }
    sprintfAdd(Buff, "fill=\"#%02x%02x%02x\" ",ColorRed(Line->FillColor),ColorGreen(Line->FillColor),ColorBlue(Line->FillColor));
  }
  else
  {
    sprintfAdd(Buff, "fill=\"none\" ");
    sprintfAdd(Buff, "opacity=\"%.5f\" ",ColorOpacity(Line->Color));
  }
  sprintfAdd(Buff, "stroke=\"#%02x%02x%02x\" ",ColorRed(Line->Color),ColorGreen(Line->Color),ColorBlue(Line->Color));
  sprintfAdd(Buff, "stroke-width=\"%d\" ",Max(1,(int)round(Line->Width * scaleV) ));
  sprintfAdd(Buff, "d=\"M %d %d ",(int)round(X0 + scaleX * Line->Point[0].X),(int)round(Y0 + scaleY * Line->Point[0].Y) );
  for(i=1;i<Line->NPoints;i++)
  {
    sprintfAdd(Buff, "L %d %d ",(int)round(X0 + scaleX * Line->Point[i].X),(int)round(Y0 + scaleY * Line->Point[i].Y) );
  }
  if(Line->CloseAndFill)
  {
    sprintfAdd(Buff, "Z ");
  }
  sprintfAdd(Buff, "\" />\n");
  return 0;
}


int AddCurveToSVG(DynamicString *Buff, CurveType *Curve, int X0, int Y0, double scaleX, double scaleY, double scaleV)
{
  int i;
  if(Curve->NPoints<2) return 0;
  sprintfAdd(Buff, "<path ");
  if((Curve->CloseAndFill) && (ColorVisible(Curve->FillColor)) )
  {
    if( ColorTransparency(Curve->FillColor)==ColorTransparency(Curve->Color) )
    {
      sprintfAdd(Buff, "opacity=\"%.5f\" ",ColorOpacity(Curve->FillColor));
    }
    else
    {
      sprintfAdd(Buff, "fill-opacity=\"%.5f\" ",ColorOpacity(Curve->FillColor));
      sprintfAdd(Buff, "stroke-opacity=\"%.5f\" ",ColorOpacity(Curve->Color));
    }
    sprintfAdd(Buff, "fill=\"#%02x%02x%02x\" ",ColorRed(Curve->FillColor),ColorGreen(Curve->FillColor),ColorBlue(Curve->FillColor));
  }
  else
  {
    sprintfAdd(Buff, "fill=\"none\" ");
    sprintfAdd(Buff, "opacity=\"%.5f\" ",ColorOpacity(Curve->Color));
  }
  sprintfAdd(Buff, "stroke=\"#%02x%02x%02x\" ",ColorRed(Curve->Color),ColorGreen(Curve->Color),ColorBlue(Curve->Color));
  sprintfAdd(Buff, "stroke-width=\"%d\" ",Max(1,(int)round(Curve->Width * scaleV) ));
  sprintfAdd(Buff, "d=\"M %d %d ",(int)round(X0 + scaleX * Curve->Point[0].X),(int)round(Y0 + scaleY * Curve->Point[0].Y) );
  for(i=1;i<Curve->NPoints;i++)
  {
    sprintfAdd(Buff, "C %d %d %d %d %d %d",(int)round(X0 + scaleX * Curve->VR[i-1].X),(int)round(Y0 + scaleY * Curve->VR[i-1].Y)
                                          ,(int)round(X0 + scaleX * Curve->VL[i  ].X),(int)round(Y0 + scaleY * Curve->VL[i  ].Y)
                                          ,(int)round(X0 + scaleX * Curve->Point[i].X),(int)round(Y0 + scaleY * Curve->Point[i].Y));
  }
  if(Curve->CloseAndFill)
  {
    sprintfAdd(Buff, "C %d %d %d %d %d %d",(int)round(X0 + scaleX * Curve->VR[i-1].X),(int)round(Y0 + scaleY * Curve->VR[i-1].Y)
                                          ,(int)round(X0 + scaleX * Curve->VL[0].X),(int)round(Y0 + scaleY * Curve->VL[0].Y)
                                          ,(int)round(X0 + scaleX * Curve->Point[0].X),(int)round(Y0 + scaleY * Curve->Point[0].Y));
  }
  sprintfAdd(Buff, "\" />\n");
  return 0;
}

int AddArcToSVG(DynamicString *Buff, ArcType *Arc, int X0, int Y0, double scaleX, double scaleY, double scaleV)
{
  if(Arc->ThetaL-Arc->ThetaR > 1.99 * Pi)
  {
    sprintfAdd(Buff, "<circle cx=\"%d\" cy=\"%d\" r=\"%d\" ",(int)round(X0 + scaleX * Arc->Point.X),(int)round(Y0 + scaleY * Arc->Point.Y), (int)round(scaleV * Arc->R));
    sprintfAdd(Buff, "stroke-width=\"%d\" stroke=\"#%02x%02x%02x\" ",Max(1,(int)round(Arc->Width * scaleV) ), ColorRed(Arc->Color), ColorGreen(Arc->Color),ColorBlue(Arc->Color) );
    if((Arc->CloseAndFill)&&(ColorVisible(Arc->FillColor)))
    {
      if( ColorTransparency(Arc->FillColor)==ColorTransparency(Arc->Color) )
      {
        sprintfAdd(Buff, "opacity=\"%.5f\" ",ColorOpacity(Arc->FillColor));
      }
      else
      {
        sprintfAdd(Buff, "fill-opacity=\"%.5f\" ",ColorOpacity(Arc->FillColor));
        sprintfAdd(Buff, "stroke-opacity=\"%.5f\" ",ColorOpacity(Arc->Color));
      }
      sprintfAdd(Buff, "fill=\"#%02x%02x%02x\" ",ColorRed(Arc->FillColor),ColorGreen(Arc->FillColor),ColorBlue(Arc->FillColor));
    }
    else
    {
      sprintfAdd(Buff, "fill=\"none\" ");
      sprintfAdd(Buff, "opacity=\"%.5f\" ",ColorOpacity(Arc->Color));
    }
    sprintfAdd(Buff,"/>\n");
  }
  else
  {
    sprintfAdd(Buff, "<path d=\"M %.2f %.2f ",round(X0 + scaleX * Arc->Point.X) + round(scaleV * Arc->R) * cos(Arc->ThetaR) ,round(Y0 + scaleY * Arc->Point.Y) - round(scaleV * Arc->R) * sin(Arc->ThetaR) );
    sprintfAdd(Buff, "A %d %d ",(int)round(scaleV * Arc->R), (int)round(scaleV * Arc->R));
    if(Arc->ThetaL-Arc->ThetaR <= Pi)
    {
      sprintfAdd(Buff, "0 0 0 ");
    }
    else
    {
      sprintfAdd(Buff, "0 1 0 ");
    }
    sprintfAdd(Buff, "%.2f %.2f ",round(X0 + scaleX * Arc->Point.X) + round(scaleV * Arc->R) * cos(Arc->ThetaL) ,round(Y0 + scaleY * Arc->Point.Y) - round(scaleV * Arc->R) * sin(Arc->ThetaL) );
    if(Arc->CloseAndFill)
    {
      sprintfAdd(Buff, "L %d %d Z\" ",(int)round(X0 + scaleX * Arc->Point.X),(int)round(Y0 + scaleY * Arc->Point.Y));
    }
    else
    {
      sprintfAdd(Buff, "\" ",(int)round(X0 + scaleX * Arc->Point.X),(int)round(Y0 + scaleY * Arc->Point.Y));
    }
    if((Arc->CloseAndFill)&&(ColorVisible(Arc->FillColor)))
    {
      if( ColorTransparency(Arc->FillColor)==ColorTransparency(Arc->Color) )
      {
        sprintfAdd(Buff, "opacity=\"%.5f\" ",ColorOpacity(Arc->FillColor));
      }
      else
      {
        sprintfAdd(Buff, "fill-opacity=\"%.5f\" ",ColorOpacity(Arc->FillColor));
        sprintfAdd(Buff, "stroke-opacity=\"%.5f\" ",ColorOpacity(Arc->Color));
      }
      sprintfAdd(Buff, "fill=\"#%02x%02x%02x\" ",ColorRed(Arc->FillColor),ColorGreen(Arc->FillColor),ColorBlue(Arc->FillColor));
    }
    else
    {
      sprintfAdd(Buff, "fill=\"none\" ");
      sprintfAdd(Buff, "opacity=\"%.5f\" ",ColorOpacity(Arc->Color));
    }
    sprintfAdd(Buff, "stroke=\"#%02x%02x%02x\" ",ColorRed(Arc->Color),ColorGreen(Arc->Color),ColorBlue(Arc->Color));
    sprintfAdd(Buff, "stroke-width=\"%d\" ",Max(1,(int)round(Arc->Width * scaleV) ));
    sprintfAdd(Buff,"/>\n");
  }
  return 0;
}

int AddBitMapToSVG(DynamicString *Buff, BitMapType *BitMap, int X0, int Y0, double scaleX, double scaleY, double scaleV)
{
  
  DynamicString A64;
  char Rotate[128] = "";
  BitMapToBase64String(BitMap, &A64);
  if(BitMap->Angle!=0)
  {
    // turned about the centre of the rectangle; the y axis of svg points down, so the angle is minus
    snprintf(Rotate, sizeof(Rotate), " transform=\"rotate(%.4f %.2f %.2f)\"", -BitMap->Angle*180/Pi,
             X0 + scaleX * (BitMap->LB.X+BitMap->RT.X)/2, Y0 + scaleY * (BitMap->LB.Y+BitMap->RT.Y)/2);
  }
  sprintfAdd(Buff,1024+(int)A64.length,"<image x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" preserveAspectRatio=\"none\"%s%s xlink:href=\"data:image/bmp;base64,%s\"/>\n",
             (int)round(X0 + scaleX * BitMap->LB.X),(int)round(Y0 + scaleY * BitMap->RT.Y),
             (int)round(scaleX * (BitMap->RT.X-BitMap->LB.X)),(int)round(-scaleY * (BitMap->RT.Y - BitMap->LB.Y)),
             BitMap->Pixelated ? " style=\"image-rendering:pixelated\"" : "", Rotate,
             A64.str);
  FreeDynamicString(&A64);
  return 0;
}

int AddTextToSVG(DynamicString *Buff, TextType *Text, int X0, int Y0, double scaleX, double scaleY, double scaleV)
{
  if(!ColorVisible(Text->Color)) return 0;
  sprintfAdd(Buff, "<text x=\"%d\" y=\"%d\" ",(int)round(X0 + scaleX * Text->Point.X),(int)round(Y0 + scaleY * Text->Point.Y));
  if((int)round(Text->Angle*180/Pi) != 0 ) sprintfAdd(Buff, "transform=\"rotate(%d %d,%d)\" ",(int)round(-Text->Angle*180/Pi),(int)round(X0 + scaleX * Text->Point.X),(int)round(Y0 + scaleY * Text->Point.Y));
  sprintfAdd(Buff, "opacity=\"%.5f\" ",ColorOpacity(Text->Color));
  sprintfAdd(Buff, "font-family=\"%s\" ",Text->Font);
  sprintfAdd(Buff,"style=\"stroke:none; fill:#%02x%02x%02x; font-size: %dpx; ",ColorRed(Text->Color),ColorGreen(Text->Color),ColorBlue(Text->Color),(int)round(scaleV * Text->Size));
  switch(Text->Alignment)
  {
    case 'l':
      sprintfAdd(Buff,"text-anchor: start\" ");
    break;
    case 'r':
      sprintfAdd(Buff,"text-anchor: end\" ");
    break;
    case 'c':
      sprintfAdd(Buff,"text-anchor: middle\"\n");
    break;
    default:
      printf("ERROR: unrecognized allignment token for text in graph to svg\n");
      sprintfAdd(Buff,"\"\n");
  }
  sprintfAdd(Buff, ">\n%s\n</text>\n",Text->String);
  return 0;
}


int AddLayersToSVG(DynamicString *Buff, GraphicsType *G, int YMax)
{
  int i;
  
  int X0,Y0;
  double const scaleX =  SVG_PIXELS_PER_MM * G->CanvasSize.X ;
  double const scaleY = -SVG_PIXELS_PER_MM * G->CanvasSize.Y ;
  double const scaleV = sqrt(Abs(scaleX * scaleY));
  
  X0 = round(SVG_PIXELS_PER_MM * G->CanvasOrigin.X);
  Y0 = YMax - round(SVG_PIXELS_PER_MM * G->CanvasOrigin.Y);
  
  for(i=0;i<G->NBitMap;i++)
  {
    AddBitMapToSVG(Buff, &(G->BitMap[i]), X0, Y0, scaleX, scaleY, scaleV);
  }
  for(i=0;i<G->NLine;i++)
  {
    AddLineToSVG(Buff, &(G->Line[i]), X0, Y0, scaleX, scaleY, scaleV);
  }
  for(i=0;i<G->NCurve;i++)
  {
    AddCurveToSVG(Buff, &(G->Curve[i]), X0, Y0, scaleX, scaleY, scaleV);
  }
  for(i=0;i<G->NArc;i++)
  {
    AddArcToSVG(Buff, &(G->Arc[i]), X0, Y0, scaleX, scaleY, scaleV);
  }
  for(i=0;i<G->NText;i++)
  {
    AddTextToSVG(Buff, &(G->Text[i]), X0, Y0, scaleX, scaleY, scaleV);
  }
  for(i=1;i<G->NLayer;i++)
  {
    AddLayersToSVG(Buff, G->Layer[i], YMax);
  }
  return 0;
}

// G is a graphics type
// *PS is a string that will contain the ps text

int GraphicsToSVG(GraphicsType *G, DynamicString *SVG, bool RelativeSize)
{
  int XMax,YMax,X0,Y0;
  
  X0 = round(POSTSCIPT_PIXELS_PER_MM * G->CanvasOrigin.X);
  Y0 = round(POSTSCIPT_PIXELS_PER_MM * G->CanvasOrigin.Y);
  XMax = round(POSTSCIPT_PIXELS_PER_MM * G->CanvasSize.X ) + round(POSTSCIPT_PIXELS_PER_MM *  G->CanvasOrigin.X );
  YMax = round(POSTSCIPT_PIXELS_PER_MM * G->CanvasSize.Y ) + round(POSTSCIPT_PIXELS_PER_MM *  G->CanvasOrigin.Y );
  
  *SVG = InitDynamicString();
  if(RelativeSize)
  {
    AddRelativeSizePrologToSVG(SVG, X0, Y0, XMax, YMax );
  }
  else
  {
    AddPrologToSVG(SVG, X0, Y0, XMax, YMax ) ;
  }
  AddLayersToSVG(SVG, G, YMax);
  AddCloseToSVG(SVG);
  return 0;
}




int GraphicsToSVG(GraphicsType *G, char** SVG)
{
  DynamicString Buff;
  GraphicsToSVG(G, &Buff, false);
  *SVG = Buff.str;
  return 0;
}

/*  creates a string that contains an HTML page description to view SVG files
 *
 */
const char* SVGHTMLViewer =
  "<!DOCTYPE html>\n"
  "<html lang=\"en\">\n"
  "<head>\n"
  "    <meta charset=\"UTF-8\">\n"
  "    <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n"
  "    <title>Quanty SVG Viewer</title>\n"
  "    <style>\n"
  "        body, html { margin: 0; padding: 0; width: 100%; height: 100%; display: flex; flex-direction: column; align-items: center; background-color: #f0f0f0; }\n"
  "        #controlPanel { display: flex; align-items: center; margin: 20px 0; }\n"
  "        #svgContainer { flex-grow: 1; width: 100%; display: flex; justify-content: center; align-items: center; overflow: hidden; background-color: #f0f0f0; }\n"
  "        #svgContainer svg { display: block; max-width: 100%; max-height: calc(100% - 60px); }\n"
  "    </style>\n"
  "</head>\n"
  "<body>\n"
  "    <div id=\"controlPanel\">\n"
  "        <button id=\"blackBackground\">Black</button>\n"
  "        <button id=\"greyBackground\">Grey</button>\n"
  "        <button id=\"whiteBackground\">White</button> &nbsp;&nbsp;\n"
  "        <button id=\"toggleAspectRatio\">Toggle Aspect Ratio</button> &nbsp;&nbsp;\n"
  "        <input type=\"file\" id=\"fileInput\" accept=\".svg\">\n"
  "    </div>\n"
  "    <div id=\"svgContainer\"></div>\n"
  "    <script>\n"
  "        let keepAspectRatio = true;\n"
  "        document.getElementById('blackBackground').addEventListener('click', function() { document.getElementById('svgContainer').style.backgroundColor = 'black'; });\n"
  "        document.getElementById('greyBackground').addEventListener('click', function() { document.getElementById('svgContainer').style.backgroundColor = '#f0f0f0'; });\n"
  "        document.getElementById('whiteBackground').addEventListener('click', function() { document.getElementById('svgContainer').style.backgroundColor = 'white'; });\n"
  "        document.getElementById('toggleAspectRatio').addEventListener('click', function() { keepAspectRatio = !keepAspectRatio; const svgElement = document.getElementById('svgContainer').querySelector('svg'); if (svgElement) { updateSvgAttributes(svgElement); } });\n"
  "        document.getElementById('fileInput').addEventListener('click', function(event) { fileInput.value = null; });\n"
  "        document.getElementById('fileInput').addEventListener('change', function(event) { const file = event.target.files[0]; if (file && file.type === 'image/svg+xml') { const reader = new FileReader(); reader.onload = function(e) { const svgContent = e.target.result; document.getElementById('svgContainer').innerHTML = svgContent; const svgElement = document.getElementById('svgContainer').querySelector('svg'); if (svgElement) { updateSvgAttributes(svgElement); } }; reader.readAsText(file); } else { alert('Please select a valid SVG file.'); } });\n"
  "        function updateSvgAttributes(svgElement) { if (keepAspectRatio) { svgElement.style.width = '100%'; svgElement.style.height = '100%'; svgElement.setAttribute('preserveAspectRatio', 'xMidYMid meet'); } else { svgElement.style.width = '100vw'; svgElement.style.height = '100%'; svgElement.setAttribute('preserveAspectRatio', 'none'); } }\n"
  "    </script>\n"
  "</body>\n"
  "</html>\n";


const char* SVGSequencyHTMLViewer =
  "<!DOCTYPE html>\n"
  "<html lang=\"en\">\n"
  "<head>\n"
  "  <meta charset=\"UTF-8\">\n"
  "  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n"
  "  <title>Quanty SVG Sequence Player</title>\n"
  "  <style>\n"
  "      body { margin: 0; display: flex; flex-direction: column; height: 100vh; overflow: hidden;}\n"
  "      #controls { flex: 0 0 auto; display: flex; gap: 10px; padding: 10px; box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1); background-color: #f8f8f8; z-index: 1; }\n"
  "      #svgContainer { flex: 1 1 auto; display: flex; align-items: center; justify-content: center; overflow: hidden; }\n"
  "      #svgContainer svg { max-width: 100%; max-height: 100%; }\n"
  "      #fileNameDisplay { width: 200px; }\n"
  "  </style>\n"
  "</head>\n"
  "<body>\n"
  "  <div id=\"controls\"> <input type=\"file\" id=\"fileInput\" accept=\".svg\" multiple> <button id=\"playButton\">Play</button> <input type=\"text\" id=\"fileNameDisplay\" readonly> </div>\n"
  "  <div id=\"svgContainer\"></div>\n"
  "  <script>\n"
  "      document.getElementById('fileInput').addEventListener('change', handleFileSelect, false);\n"
  "      document.getElementById('playButton').addEventListener('click', playSequence, false);\n"
  "      let svgFiles = [];\n"
  "      let intervalId;\n"
  "      function handleFileSelect(event) { const files = event.target.files; svgFiles = []; for (let i = 0; i < files.length; i++) { svgFiles.push(files[i]); } }\n"
  "      function playSequence() {\n"
  "        if (svgFiles.length === 0) { alert('Please select some SVG files first.'); return; }\n"
  "        const svgContainer = document.getElementById('svgContainer');\n"
  "        const fileNameDisplay = document.getElementById('fileNameDisplay');\n"
  "        let index = 0;\n"
  "        if (intervalId) { clearInterval(intervalId); }\n"
  "        intervalId = setInterval(() => { const file = svgFiles[index]; fileNameDisplay.value = file.name; const reader = new FileReader(); reader.onload = function (e) { svgContainer.innerHTML = e.target.result; }; reader.readAsText(file); index++; if (index >= svgFiles.length) { index = 0; } }, 1000);\n"
  "      }\n"
  "  </script>\n"
  "</body>\n"
  "</html>\n";



