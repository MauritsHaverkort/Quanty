int LuaNewGraphics(lua_State * L)
{
  GraphicsType *Graphics;
  uint32_t len = (uint32_t) lua_gettop(L);
  Graphics = (GraphicsType*) lua_newuserdata(L, sizeof(GraphicsType));
  luaL_setmetatable(L, "Graphics_Type");
  InitGraphics(Graphics);
  if(len==1)
  {
    struct
    {
      double CanvasSizeX;
      double CanvasSizeY;
      double CanvasOriginX;
      double CanvasOriginY;
      uint32_t NLayer;
    } Options;
    LuaOptionsType LuaNewGraphicsOptions[]={
      {"Canvas.Size.X","double",&(Options.CanvasSizeX)},
      {"Canvas.Size.Y","double",&(Options.CanvasSizeY)},
      {"Canvas.Origin.X","double",&(Options.CanvasOriginX)},
      {"Canvas.Origin.Y","double",&(Options.CanvasOriginY)},
      {"NLayer","uint32_t",&(Options.NLayer)},
      {NULL, NULL, NULL}};
    if(LuaTestOptionsList(L, -2))
    {
      Options.CanvasSizeX=Graphics->CanvasSize.X;
      Options.CanvasSizeY=Graphics->CanvasSize.Y;
      Options.CanvasOriginX=Graphics->CanvasOrigin.X;
      Options.CanvasOriginY=Graphics->CanvasOrigin.Y;
      LuaCheckOptions(L, -2, LuaNewGraphicsOptions);
      Graphics->CanvasSize.X=Options.CanvasSizeX;
      Graphics->CanvasSize.Y=Options.CanvasSizeY;
      Graphics->CanvasOrigin.X=Options.CanvasOriginX;
      Graphics->CanvasOrigin.Y=Options.CanvasOriginY;
    }
    else
    {
      return luaL_error(L, "Error in LuaNewGraphics Graphics.New called with 1 argument that has the wrong format for options\n");
    }
  }
  else if (len!=0)
  {
    return luaL_error(L, "Error in LuaNewGraphics Graphics.New should be called with 0 or 1 argument, not with %d\n",len);
  }
  return 1;
}

int LuaAddToGraphicsMethod(lua_State *L)
{
  lua_pushvalue(L, lua_upvalueindex(1));
  lua_insert(L, 1);
  LuaAddToGraphics(L);
  return 0;
}

int LuaAddToGraphicsFunction(lua_State *L)
{
  GraphicsType *GraphicsOld = (GraphicsType *) luaL_checkudata(L, 1, "Graphics_Type");
  GraphicsType *Graphics = (GraphicsType*) lua_newuserdata(L, sizeof(GraphicsType));
  luaL_setmetatable(L, "Graphics_Type");
  CopyGraphics(GraphicsOld, Graphics);
  lua_replace(L, 1);
  LuaAddToGraphics(L);
  lua_pushvalue(L, 1);
  return 1;
}

/*   On stack at position 1 we have a graphics object
 *   On position -1 we can have options
 *   On position 2 to 2+m we have a definition of a new object that should be added to the graphics object
 *   We add the object on position 2 to 2+m and remove these elements from the stack and then go back to the algorithm untill there is nothing on the stack to add
 */
int LuaAddToGraphics(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  int i,j;
  uint32_t val;
  GraphicsType *Graphics = (GraphicsType *) luaL_checkudata(L, 1, "Graphics_Type");
  struct
  {
    uint32_t RGB;
    uint8_t  alpha;
    uint32_t FillRGB;
    uint8_t  Fillalpha;
    double   LineWidth;
    double   FontSize;
    const char* Alignment;
    double   Angle;
    const char* Font;
    bool     Close;
  } Options;
  char defaultAlignment = 'l';
  Options.RGB=0xFF000000;
  Options.FillRGB=0xFF000000;
  Options.alpha=0xFF;
  Options.Fillalpha=0xFF;
  Options.LineWidth=0.001;
  Options.FontSize=0.05;
  Options.Alignment=&defaultAlignment;
  Options.Angle=0;
  Options.Font=TIMES_NEW_ROMAN;
  Options.Close=false;
  LuaOptionsType LuaAddGraphicsOptions[]={
    {"alpha","uint8_t",&(Options.alpha)},
    {"Transparancy","uint8_t",&(Options.alpha)},
    {"Fill alpha","uint8_t",&(Options.Fillalpha)},
    {"Fill Transparancy","uint8_t",&(Options.Fillalpha)},
    {"Color","uint32_t",&(Options.RGB)},
    {"Fill Color","uint32_t",&(Options.FillRGB)},
    {"Close","bool",&(Options.Close)},
    {"Line Width","double",&(Options.LineWidth)},
    {"Size","double",&(Options.FontSize)},
    {"Alignment","char*",&(Options.Alignment)},
    {"Angle","double",&(Options.Angle)},
    {"Font","char*",&(Options.Font)},
    {NULL, NULL, NULL}};
  if(LuaTestOptionsList(L, -1))
  {
    len=len-1;
    LuaCheckOptions(L, -1, LuaAddGraphicsOptions);
    Options.RGB     = (Options.RGB & 0xFFFFFF) + (Options.alpha<<24);
    Options.FillRGB = (Options.FillRGB & 0xFFFFFF) + (Options.Fillalpha<<24);
  }
  while (len>1)
  {
    switch (lua_type(L,2))
    {
      case LUA_TNIL:
      {
        return luaL_error(L, "Error: Add nil to graphics\n");
      }
      case LUA_TBOOLEAN:
      {
        return luaL_error(L, "Error: add boolean to graphics\n");
      }
      case LUA_TLIGHTUSERDATA:
      {
        return luaL_error(L, "Error: add light userdata to graphics\n");
      }
      case LUA_TNUMBER:
      {
        return luaL_error(L, "Error: add number to graphics\n");
      }
      case LUA_TSTRING:
      {
        if(len<4)
        {
          luaL_error(L, "Error: add text to graphics expects 4 inputs: Graphics, Text, p.X, p.Y\n");
        }
        TextType text;
        text.String = (char*)lua_tostring(L, 2);
        text.Point.X = (double)luaL_checknumber(L,3);
        text.Point.Y = (double)luaL_checknumber(L,4);
        text.Color[3] = (uint8_t)(Options.RGB >> 24);
        text.Color[0] = (uint8_t)(Options.RGB >> 16);
        text.Color[1] = (uint8_t)(Options.RGB >> 8);
        text.Color[2] = (uint8_t)Options.RGB;
        text.Angle    = Options.Angle;
        text.Alignment = Options.Alignment[0];
        text.Size = Options.FontSize;
        text.Font = Options.Font;
        AddText(Graphics, &text);
        lua_remove(L,4);
        lua_remove(L,3);
        lua_remove(L,2);
        len=len-3;
        break;
      }
      case LUA_TTABLE:
      {
        // add line, curve, or arc depending on format of table
        uint32_t Tablelen = (uint32_t)lua_objlen(L, 2);
        lua_rawgeti(L, 2, 1);
        switch (lua_type(L,-1))
        {
          case LUA_TSTRING:
          {
            if( strcmp((char*)lua_tostring(L, -1), "Arc")==0 )
            {
              ArcType Arc;
              Arc.Color[3] = (uint8_t)(Options.RGB >> 24);
              Arc.Color[0] = (uint8_t)(Options.RGB >> 16);
              Arc.Color[1] = (uint8_t)(Options.RGB >> 8);
              Arc.Color[2] = (uint8_t)Options.RGB;
              Arc.Width        = Options.LineWidth;
              Arc.CloseAndFill = Options.Close;
              Arc.FillColor[3] = (uint8_t)(Options.FillRGB >> 24);
              Arc.FillColor[0] = (uint8_t)(Options.FillRGB >> 16);
              Arc.FillColor[1] = (uint8_t)(Options.FillRGB >> 8);
              Arc.FillColor[2] = (uint8_t)Options.FillRGB;
              Arc.ThetaR=0;
              Arc.ThetaL=2*Pi;
              lua_pop(L, 1);
              if( !( (Tablelen==3) || (Tablelen==5) ) )
              {
                return luaL_error(L, "Error: Input should either be \"Arc\" or not present.\n(For a line or Bezier curve, just give a list with 2 or 6 points per element)\n");
              }
              // Format {"Arc", {px,py}, R} or {"Arc", {px,py}, R, ThetaMin, ThetaMax}
              lua_rawgeti(L, 2, 2);
              if(lua_type(L,-1) != LUA_TTABLE) return luaL_error(L, "Error: First element of an arc should be a point defining its center given as a list of two numbers\n");
              if(lua_objlen(L, -1) != 2) return luaL_error(L, "Error: First element of an arc should be a point defining its center given as a list of two numbers\n");
              lua_rawgeti(L, -1, 1); //px
              Arc.Point.X=(double)luaL_checknumber(L, -1);
              lua_pop(L, 1);
              lua_rawgeti(L, -1, 2); // py
              Arc.Point.Y=(double)luaL_checknumber(L, -1);
              lua_pop(L, 2);
              lua_rawgeti(L, 2, 3); // R
              Arc.R=(double)luaL_checknumber(L, -1);
              lua_pop(L, 1);
              if(Tablelen==5)
              {
                lua_rawgeti(L, 2, 4); // thetamin
                Arc.ThetaR=(double)luaL_checknumber(L, -1);
                lua_pop(L, 1);
                while(Arc.ThetaR > TwoPi )
                {
                  Arc.ThetaR = Arc.ThetaR - TwoPi;
                }
                if(Arc.ThetaR < 0 )
                {
                  Arc.ThetaR = 0;
                }
                lua_rawgeti(L, 2, 5); // thetamax
                Arc.ThetaL=(double)luaL_checknumber(L, -1);
                lua_pop(L, 1);
                while(Arc.ThetaL > TwoPi )
                {
                  Arc.ThetaL = Arc.ThetaL - TwoPi;
                }
                if(Arc.ThetaL < 0 )
                {
                  Arc.ThetaL = 0;
                }
                if(Arc.ThetaL <= Arc.ThetaR )
                {
                  Arc.ThetaL = Arc.ThetaL + TwoPi;
                }
              }
              AddArc(Graphics, &Arc);
              lua_remove(L,2);
              len--;
              break;
            }
            else if ( (strcmp((char*)lua_tostring(L, -1), "BitMap")==0 ) || (strcmp((char*)lua_tostring(L, -1), "Bitmap")==0 ) )
            {
              BitMapType BitMap;
              lua_pop(L, 1);
              if( !( (Tablelen==4) || (Tablelen==5) ) )
              {
                return luaL_error(L, "Error: Input for BitMap should either be {\"BitMap\", {px,py}, {px,py}, Array_uint8_t*4} or {\"BitMap\", {px,py}, {px,py}, Array_uint8_t*4, {Transparant color uint8_t*4}} or {\"BitMap\", {px,py}, {px,py}, ColorFunction, Matrix} \n");
              }
              // Format {"BitMap", {px,py}, {px,py}, Array_uint8_t*4}
              //     or {"BitMap", {px,py}, {px,py}, Array_uint8_t*4, {Transparent_Color_uint8_t*4}}
              //     or {"BitMap", {px,py}, {px,py}, Colorfunction, Matrix}
              lua_rawgeti(L, 2, 2);
              if(lua_type(L,-1) != LUA_TTABLE) return luaL_error(L, "Error: First element of a BitMap should be a point defining its left bottom point given as a list of two numbers\n");
              if(lua_objlen(L, -1) != 2) return luaL_error(L, "Error: First element of a BitMap should be a point defining its left bottom point given as a list of two numbers\n");
              lua_rawgeti(L, -1, 1); //px
              BitMap.LB.X=(double)luaL_checknumber(L, -1);
              lua_pop(L, 1);
              lua_rawgeti(L, -1, 2); // py
              BitMap.LB.Y=(double)luaL_checknumber(L, -1);
              lua_pop(L, 2);
              lua_rawgeti(L, 2, 3);
              if(lua_type(L,-1) != LUA_TTABLE) return luaL_error(L, "Error: Second element of a BitMap should be a point defining its right top point given as a list of two numbers\n");
              if(lua_objlen(L, -1) != 2) return luaL_error(L, "Error: Second element of a BitMap should be a point defining its right top point given as a list of two numbers\n");
              lua_rawgeti(L, -1, 1); //px
              BitMap.RT.X=(double)luaL_checknumber(L, -1);
              lua_pop(L, 1);
              lua_rawgeti(L, -1, 2); // py
              BitMap.RT.Y=(double)luaL_checknumber(L, -1);
              lua_pop(L, 2);
              lua_rawgeti(L, 2, 4);
              if(lua_type(L,-1) == LUA_TTABLE)
              {
                BitMap.NpointsY = (int)lua_objlen(L, -1);
                lua_rawgeti(L, -1, 1);
                if(lua_type(L,-1) != LUA_TTABLE) return luaL_error(L, "Error: Bitmap should be a matrix, i.e. a table of tables\n");
                BitMap.NpointsX = (int)lua_objlen(L, -1);
                BitMap.Pixels = (uint32_t*)malloc(BitMap.NpointsX * BitMap.NpointsY * sizeof(uint32_t));
                lua_pop(L, 1);
                for(j=0;j<BitMap.NpointsY ;j++)
                {
                  lua_rawgeti(L, -1, j+1);
                  if(BitMap.NpointsX != lua_objlen(L, -1)) return luaL_error(L, "Error: at %d: Bitmap should be a matrix, i.e. a table of tables with same dimensions\n",j);
                  for(i=0;i<BitMap.NpointsX ;i++)
                  {
                    lua_rawgeti(L, -1, i+1);
                    if(lua_type(L,-1) == LUA_TTABLE)
                    {
                      switch (lua_objlen(L, -1))
                      {
                        case 1: // BW
                        {
                          lua_rawgeti(L, -1, 1);
                          val = (uint32_t)luaL_checkinteger(L,-1);
                          if(val<0)
                          {
                            val=0;
                          }
                          if(val>255)
                          {
                            val=255;
                          }
                          BitMap.Pixels[i+j*BitMap.NpointsX] = 0xFF000000 + (val) + (val<<8) + (val<<16);
                          lua_pop(L, 1);
                          break;
                        }
                        case 3: //RGB
                        {
                          
                          BitMap.Pixels[i+j*BitMap.NpointsX] = 0xFF;

                          lua_rawgeti(L, -1, 1);
                          val = (uint32_t)luaL_checkinteger(L,-1);
                          lua_pop(L, 1);
                          if(val<0)
                          {
                            val=0;
                          }
                          if(val>255)
                          {
                            val=255;
                          }
                          BitMap.Pixels[i+j*BitMap.NpointsX] = ( BitMap.Pixels[i+j*BitMap.NpointsX]<<8 ) + val;

                          lua_rawgeti(L, -1, 2);
                          val = (uint32_t)luaL_checkinteger(L,-1);
                          lua_pop(L, 1);
                          if(val<0)
                          {
                            val=0;
                          }
                          if(val>255)
                          {
                            val=255;
                          }

                          BitMap.Pixels[i+j*BitMap.NpointsX] = ( BitMap.Pixels[i+j*BitMap.NpointsX]<<8 ) + val;
                          
                          lua_rawgeti(L, -1, 3);
                          val = (uint32_t)luaL_checkinteger(L,-1);
                          lua_pop(L, 1);
                          if(val<0)
                          {
                            val=0;
                          }
                          if(val>255)
                          {
                            val=255;
                          }
                          BitMap.Pixels[i+j*BitMap.NpointsX] = ( BitMap.Pixels[i+j*BitMap.NpointsX]<<8 ) + val;
                          break;
                        }
                        case 4: //RGBalpha
                        {
                          lua_rawgeti(L, -1, 4);
                          val = (uint32_t)luaL_checkinteger(L,-1);
                          lua_pop(L, 1);
                          if(val<0)
                          {
                            val=0;
                          }
                          if(val>255)
                          {
                            val=255;
                          }
                          BitMap.Pixels[i+j*BitMap.NpointsX] = val;
                          
                          lua_rawgeti(L, -1, 1);
                          val = (uint32_t)luaL_checkinteger(L,-1);
                          lua_pop(L, 1);
                          if(val<0)
                          {
                            val=0;
                          }
                          if(val>255)
                          {
                            val=255;
                          }
                          BitMap.Pixels[i+j*BitMap.NpointsX] = ( BitMap.Pixels[i+j*BitMap.NpointsX]<<8 ) + val;
                          
                          lua_rawgeti(L, -1, 2);
                          val = (uint32_t)luaL_checkinteger(L,-1);
                          lua_pop(L, 1);
                          if(val<0)
                          {
                            val=0;
                          }
                          if(val>255)
                          {
                            val=255;
                          }
                          BitMap.Pixels[i+j*BitMap.NpointsX] = ( BitMap.Pixels[i+j*BitMap.NpointsX]<<8 ) + val;
                          
                          lua_rawgeti(L, -1, 3);
                          val = (uint32_t)luaL_checkinteger(L,-1);
                          lua_pop(L, 1);
                          if(val<0)
                          {
                            val=0;
                          }
                          if(val>255)
                          {
                            val=255;
                          }
                          BitMap.Pixels[i+j*BitMap.NpointsX] = ( BitMap.Pixels[i+j*BitMap.NpointsX]<<8 ) + val;

                          break;
                        }
                        default:
                        {
                          return luaL_error(L, "ERROR; A pixel in a bitmap should have length 1 (BW), 3 (RGB) or 4 (RGBalpha)\n");
                          break;
                        }
                      }
                    }
                    else
                    {
                      val = (uint32_t)luaL_checkinteger(L,-1);
                      if(val<0)
                      {
                        val=0;
                      }
                      if(val>255)
                      {
                        val=255;
                      }
                      BitMap.Pixels[i+j*BitMap.NpointsX] = 0xFF000000 + (val) + (val<<8) + (val<<16);
                    }
                    lua_pop(L, 1);
                  }
                  lua_pop(L, 1);
                }
                lua_pop(L, 1);
              }
              else
              { // length = 5 and position 4 is not a talbe, collorfun matrix input option
                
              }
              AddBitMapInherentMemory(Graphics, &BitMap);
              lua_remove(L,2);
              len--;
              break;
            }
            else
            {
              return luaL_error(L, "Error: Input %s should either be \"Arc\", \"BitMap\" or not present.\n(For a line or Bezier curve, just give a list with 2 or 6 points per element)\n",(char*)lua_tostring(L, -1));
            }
            break;
          }
          case LUA_TTABLE:
          {
            switch (lua_objlen(L, -1))
            {
              case 2:
              {
                lua_pop(L, 1);
                LineType Line;
                Line.Color[3] = (uint8_t)(Options.RGB >> 24);
                Line.Color[0] = (uint8_t)(Options.RGB >> 16);
                Line.Color[1] = (uint8_t)(Options.RGB >> 8);
                Line.Color[2] = (uint8_t)Options.RGB;
                Line.FillColor[3] = (uint8_t)(Options.FillRGB >> 24);
                Line.FillColor[0] = (uint8_t)(Options.FillRGB >> 16);
                Line.FillColor[1] = (uint8_t)(Options.FillRGB >> 8);
                Line.FillColor[2] = (uint8_t)Options.FillRGB;
                Line.Width=Options.LineWidth;
                Line.CloseAndFill=Options.Close;
                Line.NPoints=Tablelen;
                Line.Point=(PointType*)malloc(Tablelen*sizeof(PointType));
                for(i=0;i<Tablelen;i++)
                {
                  lua_rawgeti(L, 2, i+1);
                  if ( lua_objlen(L, -1)!=2 ) return luaL_error(L, "ERROR: at point position %d.\nDefinition of a line should be a table of points with length 2",i+1);
                  lua_rawgeti(L, -1, 1); // X
                  Line.Point[i].X=(double)luaL_checknumber(L, -1);
                  lua_pop(L, 1);
                  lua_rawgeti(L, -1, 2); // Y
                  Line.Point[i].Y=(double)luaL_checknumber(L, -1);
                  lua_pop(L, 1);
                  lua_pop(L, 1);
                }
                AddLineInherentLineMemory(Graphics, &Line);
                lua_remove(L,2);
                len--;
                break;
              }
              case 6:
              {
                lua_pop(L, 1);
                CurveType Curve;
                Curve.Color[3] = (uint8_t)(Options.RGB >> 24);
                Curve.Color[0] = (uint8_t)(Options.RGB >> 16);
                Curve.Color[1] = (uint8_t)(Options.RGB >> 8);
                Curve.Color[2] = (uint8_t)Options.RGB;
                Curve.FillColor[3] = (uint8_t)(Options.FillRGB >> 24);
                Curve.FillColor[0] = (uint8_t)(Options.FillRGB >> 16);
                Curve.FillColor[1] = (uint8_t)(Options.FillRGB >> 8);
                Curve.FillColor[2] = (uint8_t)Options.FillRGB;
                Curve.Width=Options.LineWidth;
                Curve.CloseAndFill=Options.Close;
                Curve.NPoints=Tablelen;
                Curve.VL=(PointType*)malloc(Tablelen*sizeof(PointType));
                Curve.Point=(PointType*)malloc(Tablelen*sizeof(PointType));
                Curve.VR=(PointType*)malloc(Tablelen*sizeof(PointType));
                for(i=0;i<Tablelen;i++)
                {
                  lua_rawgeti(L, 2, i+1);
                  if ( lua_objlen(L, -1)!=6 ) return luaL_error(L, "ERROR: at point position %d.\nDefinition of a curve should be a table of 3 points with length 6",i+1);
                  lua_rawgeti(L, -1, 1); // X
                  Curve.VL[i].X=(double)luaL_checknumber(L, -1);
                  lua_pop(L, 1);
                  lua_rawgeti(L, -1, 2); // Y
                  Curve.VL[i].Y=(double)luaL_checknumber(L, -1);
                  lua_pop(L, 1);
                  
                  lua_rawgeti(L, 2, i+1);
                  lua_rawgeti(L, -1, 3); // X
                  Curve.Point[i].X=(double)luaL_checknumber(L, -1);
                  lua_pop(L, 1);
                  lua_rawgeti(L, -1, 4); // Y
                  Curve.Point[i].Y=(double)luaL_checknumber(L, -1);
                  lua_pop(L, 1);
                  
                  lua_rawgeti(L, 2, i+1);
                  lua_rawgeti(L, -1, 5); // X
                  Curve.VR[i].X=(double)luaL_checknumber(L, -1);
                  lua_pop(L, 1);
                  lua_rawgeti(L, -1, 6); // Y
                  Curve.VR[i].Y=(double)luaL_checknumber(L, -1);
                  lua_pop(L, 1);
                  
                  lua_pop(L, 1);
                }
                AddCurveInherentMemory(Graphics, &Curve);
                lua_remove(L,2);
                len--;
                break;
              }
              default:
              {
                lua_pop(L, 1);
                return luaL_error(L, "Error: elmenet should be either a list of two coordinates (x,y for a line) or a list of 6 coordinates VL, P, VR for a Bezier curve\n");
              }
            }
            break;
          }
          default:
          {
            return luaL_error(L, "Error: elmenet should be either {\"Arc\",R,P.x,P.y} or a line or a bezier curve\n");
          }
        }
        break;
      }
      case LUA_TUSERDATA:
      {
        if (luaL_testudata(L, 2, "Graphics_Type") != NULL)
        {
          GraphicsType *NewLayer = (GraphicsType *) luaL_checkudata(L, 2, "Graphics_Type");
          AddLayer(Graphics,NewLayer);
          lua_remove(L,2);
          len--;
          break;
        }
        if (luaL_testudata(L, 2, "BitMap_Type") != NULL)
        {
          BitMapType *BitMap = (BitMapType *) luaL_checkudata(L, 2, "BitMap_Type");
          AddBitMap(Graphics, BitMap);
          lua_remove(L,2);
          len--;
          break;
        }
        return luaL_error(L, "Error: add incorrect user data to graph\n");
      }
      default:
      {
        return luaL_error(L, "Error: add unrecognized object to graphics\n");
      }
    }
  }

  return 0;
}

int LuaPlotGraphics(lua_State * L)
{
  //First element on stack is either a
  // * Spectra Type
  // * Interpolation function
  // * Array of values {{x1,y1},{x2,y2},...,{xn,yn}} with x real and y real or complex
  // * Response function type
  // or a two dimensional array of these
  
  // second element can be a list of options
  int i;
  GraphicsType* Graphics;
  const char* ScaleFormat = "%.0f";
  const char* EmptyString = "";
  SpectraTypeArray *a,aa;
  aa.NSpectra=0;
  aa.Spec=NULL;
  aa.Name=NULL;
  a=&aa;
  uint32_t len = (uint32_t) lua_gettop(L);
  Graphics = (GraphicsType*) lua_newuserdata(L, sizeof(GraphicsType));
  luaL_setmetatable(L, "Graphics_Type");
  InitGraphics(Graphics);
  struct
  {
    double   ScaleOffDiagonal;
    double   ShiftOffDiagonal;
    double   WatterfallShift;
    
    int      Nrow;
    int      Ncolumn;
    int      PlotPoints;
    
    FramePlotOptions Frame;
    LinePlotOptionsArray Style;
    
    bool     PlotIm;
    bool     PlotRe;
    bool     Negate;
  } Options;
  Options.Frame.PlotRange[0].X=-INFINITY;
  Options.Frame.PlotRange[0].Y=-INFINITY;
  Options.Frame.PlotRange[1].X=INFINITY;
  Options.Frame.PlotRange[1].Y=INFINITY;
  Options.Frame.Margins[0].X=0.08;
  Options.Frame.Margins[0].Y=0.08*sqrt(2.0);
  Options.Frame.Margins[1].X=0.01;
  Options.Frame.Margins[1].Y=0.01*sqrt(2.0);
  Options.Frame.LineWidth=0.001;
  Options.Frame.LineColor[0]=0;
  Options.Frame.LineColor[1]=0;
  Options.Frame.LineColor[2]=0;
  Options.Frame.LineColor[3]=0xFF;
  Options.Frame.FontSize=0.04;
  Options.Frame.Font=TIMES_NEW_ROMAN;
  Options.Frame.XLabel=(char*)EmptyString;
  Options.Frame.YLabel=(char*)EmptyString;
  Options.Frame.FirstXTick=-INFINITY;
  Options.Frame.dXTick=1.0;
  Options.Frame.NMinorXTick=9;
  Options.Frame.FirstYTick=-INFINITY;
  Options.Frame.dYTick=1.0;
  Options.Frame.NMinorYTick=9;
  Options.Frame.PlotXScaleFormat=(char*)ScaleFormat;
  Options.Frame.PlotYScaleFormat=(char*)ScaleFormat;
  Options.Frame.XAxis=0;
  Options.Frame.YAxis=0;
  Options.Style.NAllocLine=3;
  Options.Style.Line = (LinePlotOptions*)malloc(sizeof(LinePlotOptions)*Options.Style.NAllocLine);
  Options.Style.Line[0].LineWidth=0.002;
  Options.Style.Line[0].LineColor[0]=0;
  Options.Style.Line[0].LineColor[1]=0;
  Options.Style.Line[0].LineColor[2]=0;
  Options.Style.Line[0].LineColor[3]=0xFF;
  Options.Style.Line[0].PointSize=0;
  Options.Style.Line[0].PointColor[0]=0;
  Options.Style.Line[0].PointColor[1]=0;
  Options.Style.Line[0].PointColor[2]=0;
  Options.Style.Line[0].PointColor[3]=0xFF;
  Options.Style.Line[1].LineWidth=0.002;
  Options.Style.Line[1].LineColor[0]=0xFF;
  Options.Style.Line[1].LineColor[1]=0;
  Options.Style.Line[1].LineColor[2]=0;
  Options.Style.Line[1].LineColor[3]=0xFF;
  Options.Style.Line[1].PointSize=0;
  Options.Style.Line[1].PointColor[0]=0xFF;
  Options.Style.Line[1].PointColor[1]=0;
  Options.Style.Line[1].PointColor[2]=0;
  Options.Style.Line[1].PointColor[3]=0xFF;
  Options.Style.Line[2].LineWidth=0.002;
  Options.Style.Line[2].LineColor[0]=0;
  Options.Style.Line[2].LineColor[1]=0;
  Options.Style.Line[2].LineColor[2]=0xFF;
  Options.Style.Line[2].LineColor[3]=0xFF;
  Options.Style.Line[2].PointSize=0;
  Options.Style.Line[2].PointColor[0]=0;
  Options.Style.Line[2].PointColor[1]=0;
  Options.Style.Line[2].PointColor[2]=0xFF;
  Options.Style.Line[2].PointColor[3]=0xFF;
  
  Options.PlotIm=true;
  Options.PlotRe=true;
  Options.Negate=false;
  Options.Nrow=1;
  Options.Ncolumn=1;
  Options.PlotPoints=600;
  Options.ScaleOffDiagonal=1;
  Options.ShiftOffDiagonal=0;
  Options.WatterfallShift=0;
  LuaOptionsType LuaGraphicsOptions[]={
    {"Frame","FramePlotOptions",&(Options.Frame)},
    {"Axis", "FramePlotOptions",&(Options.Frame)},
    {"LineStyle", "LinePlotOptionsArray",&(Options.Style)},
    {"plotIm","bool",&(Options.PlotIm)},
    {"plotRe","bool",&(Options.PlotRe)},
    {"Negate","bool",&(Options.Negate)},
    {"Nrow","int",&(Options.Nrow)},
    {"Ncolumn","int",&(Options.Ncolumn)},
    {"PlotPoints","int",&(Options.PlotPoints)},
    {"ScaleOffDiagonal","double",&(Options.ScaleOffDiagonal)},
    {"ShiftOffDiagonal","double",&(Options.ShiftOffDiagonal)},
    {"WatterfallShift", "double",&(Options.WatterfallShift)},
    {NULL, NULL, NULL}};
  if(LuaTestOptionsList(L, len))
  {
    LuaCheckOptions(L, len, LuaGraphicsOptions);
    len=len-1;
  }
  // select the type of graph to plot
  if(lua_type(L,1)==LUA_TTABLE)
  {
    // test if the table has a metatable related to a response function
    // if this is the case jump to the processing of response functions
    if (lua_getmetatable(L, 1)) // Get the table's metatable and push it on the stack if excists
    {
      lua_getglobal(L, "ResponseFunctionMeta");
      if (lua_rawequal(L, -1, -2))
      {
        lua_pop(L, 2);
        goto LabelResponseFunctionPlot;
      }
      lua_pop(L, 2);
    }
    // test if we can turn the table to a SpectraTypeArray and then jump to the corresponding section
    // we want all spectra to be in the same grid.
    // Define Emin, Emax, NE and dE
    double GammaG, GammaL;
    lua_getfield(L,1,"GammaG");
    if (lua_isnil(L,-1))
    {
      GammaG=0;
    }
    else
    {
      GammaG=(double)luaL_checknumber(L, -1);
    }
    lua_pop(L,1);
    lua_getfield(L,1,"GammaL");
    if (lua_isnil(L,-1))
    {
      GammaL=0;
    }
    else
    {
      GammaL=(double)luaL_checknumber(L, -1);
    }
    lua_pop(L,1);
    double Emin=Options.Frame.PlotRange[0].X;
    double Emax=Options.Frame.PlotRange[1].X;
    if(Emin==-INFINITY) Emin=-1;
    if(Emax== INFINITY) Emax= 1;
    int NE=Options.PlotPoints;
    a->NSpectra = (uint32_t)lua_objlen(L, 1);
    a->Spec = (SpectraType*)malloc(a->NSpectra*sizeof(SpectraType));
    a->Name=NULL;
    for(i=0;i<a->NSpectra;i++)
    {
      // read element i of the table and see if we can convert it to a spectrum
      lua_rawgeti(L, 1, i+1);
      // a list of if then else that goes throught the options
      if(lua_type(L,-1)==LUA_TTABLE)
      {
        // test if the table has a metatable related to a response function
        if (lua_getmetatable(L, -1)) // Get the table's metatable and push it on the stack if excists
        {
          lua_getglobal(L, "ResponseFunctionMeta");
          if (lua_rawequal(L, -1, -2))
          {
            lua_pop(L, 2);
            // we found a response function as table. Convert it to spectra
            ResponsefunctionType *A;
            SpectraTypeArray Spec;
            A = luaL_checkResponsefunction(L,-1);
            ResponseFunctionToSpectra(A, &Spec, 0, Emin, Emax, GammaG, GammaL, NE);
            if(Options.WatterfallShift!=0)
            {
              for(int j=0;j<Spec.NSpectra;j++)
              {
                for(int k=0;k<Spec.Spec[j].N;k++)
                {
                  Spec.Spec[j].Re[k]+=i*Options.WatterfallShift;
                }
              }
            }
            free(Spec.Name);
            a->Spec[i]=Spec.Spec[0];
            FreeResponsefunctionIncludingOriginalPointer(A);
          }
          else
          {
            lua_pop(L, 2);
            return luaL_error(L, "Found a table to plot, but on position %d there is an object I can not convert to a spectrum\n",i+1);
          }
        }
        else // we found a table without meta data. We can not convert this to a spectra
        {
          return luaL_error(L, "Found a table to plot, but on position %d there is an object I can not convert to a spectrum\n",i+1);
        }
      }
      else if (lua_type(L,-1)==LUA_TUSERDATA)
      {
        if (luaL_testudata(L, -1, "Spectra_Type") != NULL)
        {
          //We found a spectra copy it
          SpectraTypeArray *Spec;
          Spec = (SpectraTypeArray *) luaL_checkudata(L, -1, "Spectra_Type");
          if(Spec->NSpectra!=1) return luaL_error(L, "A table of spectra can only contain spectra objects with 1 spectrum. At position %d we found %d spectra\n",i+1,Spec->NSpectra);
          SpectrumAllocate(Spec->Spec[0], &(a->Spec[i]));
          CopySpectra(Spec->Spec[0], &(a->Spec[i]));
        }
        else if (luaL_testudata(L, -1, "GaussianRadialFunction_Type") != NULL)
        {
          // TODO: we found a Gaussian radial function, convert it to a spectra
          return luaL_error(L, "TODO: Plot of Gaussian Radial function\n");
        }
        else if (luaL_testudata(L, -1, "InterpolatingFunction_Type") != NULL)
        {
          // TODO: we found an interpolating function, convert it to a spectra
          return luaL_error(L, "TODO plot of Interpolating Function\n");
        }
        else if (luaL_testudata(L, -1, "ResponseFunction_Type") != NULL)
        {
          // We found a response function as user type, convert it
          ResponsefunctionType *A;
          SpectraTypeArray Spec;
          A = luaL_checkResponsefunction(L,-1);
          ResponseFunctionToSpectra(A, &Spec, 0, Emin, Emax, GammaG, GammaL, NE);
          free(Spec.Name);
          a->Spec[i]=Spec.Spec[0];
        }
      }
      else
      {
        return luaL_error(L, "Found a table to plot, but on position %d there is an object I can not convert to a spectrum\n",i+1);
      }
      lua_pop(L, 1);
    }
    goto labelSpectraTypeArrayPlot;
  }
  else if (lua_type(L,1)==LUA_TUSERDATA)
  {
    if (luaL_testudata(L, 1, "Spectra_Type") != NULL)
    {
      a = (SpectraTypeArray *) luaL_checkudata(L, 1, "Spectra_Type");
      labelSpectraTypeArrayPlot:
      if(Options.Frame.PlotRange[0].X==-INFINITY) Options.Frame.PlotRange[0].X = a->Spec[0].Emin;
      if(Options.Frame.PlotRange[1].X== INFINITY) Options.Frame.PlotRange[1].X = a->Spec[0].Emax;
      if(Options.Frame.PlotRange[0].Y==-INFINITY) Options.Frame.PlotRange[0].Y = -1;
      if(Options.Frame.PlotRange[1].Y== INFINITY) Options.Frame.PlotRange[1].Y =  1;
      if(Options.Frame.FirstXTick==-INFINITY) Options.Frame.FirstXTick = Options.Frame.PlotRange[0].X;
      if(Options.Frame.FirstYTick==-INFINITY) Options.Frame.FirstYTick = Options.Frame.PlotRange[0].Y;
      
      if((Options.Nrow==1)&&(Options.Ncolumn==1))
      {
        // make a single graph with possible multiple spectra
        LinePlotOptions *LPORe = (LinePlotOptions*) malloc(sizeof(LinePlotOptions) * a->NSpectra);
        LinePlotOptions *LPOIm = (LinePlotOptions*) malloc(sizeof(LinePlotOptions) * a->NSpectra);
        if(a->NSpectra==1)
        {
          if((Options.PlotIm==true) && (Options.PlotRe==false))
          {
            LPORe[0] = Options.Style.Line[0];
            LPOIm[0] = Options.Style.Line[0];
            LPORe[0].LineWidth=0;
            LPORe[0].LineColor[3]=0;
            LPORe[0].PointSize=0;
            LPORe[0].PointColor[3]=0;
          }
          else if((Options.PlotIm==false) && (Options.PlotRe==true))
          {
            LPORe[0] = Options.Style.Line[0];
            LPOIm[0] = Options.Style.Line[0];
            LPOIm[0].LineWidth=0;
            LPOIm[0].LineColor[3]=0;
            LPOIm[0].PointSize=0;
            LPOIm[0].PointColor[3]=0;
          }
          else if((Options.PlotIm==true) && (Options.PlotRe==true))
          {
            LPORe[0] = Options.Style.Line[1];
            LPOIm[0] = Options.Style.Line[2];
          }
          else
          {
            return luaL_error(L, "Plot a spectrum but both PlotRe and PlotIm are set to false so there is nothing to do here\n");
          }
        }
        else
        { // plot multiple spectra at once, set the styles for Re and Im
          if((Options.PlotIm==true) && (Options.PlotRe==false))
          {
            for(i=0;i<a->NSpectra;i++)
            {
              LPORe[i] = Options.Style.Line[(i%(Options.Style.NAllocLine))];
              LPOIm[i] = Options.Style.Line[(i%(Options.Style.NAllocLine))];
              LPORe[i].LineWidth=0;
              LPORe[i].LineColor[3]=0;
              LPORe[i].PointSize=0;
              LPORe[i].PointColor[3]=0;
            }
          }
          else if((Options.PlotIm==false) && (Options.PlotRe==true))
          {
            for(i=0;i<a->NSpectra;i++)
            {
              LPORe[i] = Options.Style.Line[(i%(Options.Style.NAllocLine))];
              LPOIm[i] = Options.Style.Line[(i%(Options.Style.NAllocLine))];
              LPOIm[i].LineWidth=0;
              LPOIm[i].LineColor[3]=0;
              LPOIm[i].PointSize=0;
              LPOIm[i].PointColor[3]=0;
            }
          }
          else if((Options.PlotIm==true) && (Options.PlotRe==true))
          {
            for(i=0;i<a->NSpectra;i++)
            {
              LPORe[i] = Options.Style.Line[(  (2*i)%(Options.Style.NAllocLine))];
              LPOIm[i] = Options.Style.Line[((2*i+1)%(Options.Style.NAllocLine))];
            }
          }
          else
          {
            return luaL_error(L, "Plot a spectrum but both PlotRe and PlotIm are set to false so there is nothing to do here\n");
          }
        }
        PlotSpectraTypeArrayToGraphics(Graphics, a, LPORe, LPOIm, &(Options.Frame), Options.Negate, Options.WatterfallShift);
        free(LPORe);
        free(LPOIm);
      }
      else
      {
        // plot a matrix of spectra
        if ( (Options.Nrow * Options.Ncolumn) != a->NSpectra ) return luaL_error(L, "plotting a spectrum array with %d rows and %d columns, but the data array has %d elements\n",Options.Nrow,Options.Ncolumn,a->NSpectra);
        LinePlotOptions LPORe, LPOIm;
        if((Options.PlotIm==true) && (Options.PlotRe==false))
        {
          LPORe = Options.Style.Line[0];
          LPOIm = Options.Style.Line[0];
          LPORe.LineWidth=0;
          LPORe.LineColor[3]=0;
          LPORe.PointSize=0;
          LPORe.PointColor[3]=0;
        }
        else if((Options.PlotIm==false) && (Options.PlotRe==true))
        {
          LPORe = Options.Style.Line[0];
          LPOIm = Options.Style.Line[0];
          LPOIm.LineWidth=0;
          LPOIm.LineColor[3]=0;
          LPOIm.PointSize=0;
          LPOIm.PointColor[3]=0;
        }
        else if((Options.PlotIm==true) && (Options.PlotRe==true))
        {
          LPORe = Options.Style.Line[1];
          LPOIm = Options.Style.Line[2];
        }
        else
        {
          return luaL_error(L, "Plot a spectrum but both PlotRe and PlotIm are set to false so there is nothing to do here\n");
        }
        PlotSpectraTypeArrayToGraphicsArray(Graphics, a, Options.Ncolumn, Options.Nrow, &LPORe, &LPOIm, &(Options.Frame),  Options.Negate, Options.ScaleOffDiagonal, Options.ShiftOffDiagonal);
      }
    }
    else if (luaL_testudata(L, 1, "GaussianRadialFunction_Type") != NULL)
    {
      // TODO:
      return luaL_error(L, "TODO: Plot of Gaussian Radial function\n");
    }
    else if (luaL_testudata(L, 1, "InterpolatingFunction_Type") != NULL)
    {
      // TODO:
      return luaL_error(L, "TODO plot of Interpolating Function\n");
    }
    else if (luaL_testudata(L, 1, "ResponseFunction_Type") != NULL)
    {
      LabelResponseFunctionPlot:
      SpectraTypeArray Spec;
      ResponsefunctionType *A;
      double GammaL=0;
      double GammaG=0;
      double Emin,Emax;
      int NE=0;
      A = luaL_checkResponsefunction(L,1);
      if(len>1) GammaL=(double)luaL_checknumber(L, 2);
      if(len>2) GammaG=(double)luaL_checknumber(L, 3);
      if(len>3) NE=(int)luaL_checkinteger(L, 3);
      if(len>4) return luaL_error(L, "For plotting a response function I expect the first argument to be a response function and the 2nd, 3rd and 4th are optional the Lorentzian and Gaussian broadingings and the number of points in the spectrum.\n In this case I got %d argumens\n",len);
      
      if(Options.Frame.PlotRange[0].X==-INFINITY) Options.Frame.PlotRange[0].X = -1;
      if(Options.Frame.PlotRange[1].X== INFINITY) Options.Frame.PlotRange[1].X =  1;
      if(Options.Frame.PlotRange[0].Y==-INFINITY) Options.Frame.PlotRange[0].Y = -1;
      if(Options.Frame.PlotRange[1].Y== INFINITY) Options.Frame.PlotRange[1].Y =  1;
      if(Options.Frame.FirstXTick==-INFINITY) Options.Frame.FirstXTick = Options.Frame.PlotRange[0].X;
      if(Options.Frame.FirstYTick==-INFINITY) Options.Frame.FirstYTick = Options.Frame.PlotRange[0].Y;
      
      if(Options.Nrow!=1) return luaL_error(L, "Invalid option Nrow when plotting response functions\n");
      if(Options.Ncolumn!=1) return luaL_error(L, "Invalid option Ncolumn when plotting response functions\n");
          
      Emin=Options.Frame.PlotRange[0].X;
      Emax=Options.Frame.PlotRange[1].X;
      if(NE==0)
      {
        if(GammaL==0)
        {
          if(GammaG==0)
          {
            NE = Options.PlotPoints;
          }
          else
          {
            NE = 1+10*abs((Emax-Emin) / GammaG);
          }
        }
        else
        {
          NE = 1+10*abs((Emax-Emin) / GammaL);
        }
      }
      ResponseFunctionToSpectra(A, &Spec, 0, Emin, Emax, GammaG, GammaL, NE);
      
      if(A->Block==false)
      {
        // make a single graph with possible multiple spectra
        LinePlotOptions *LPORe = (LinePlotOptions*) malloc(sizeof(LinePlotOptions) * Spec.NSpectra);
        LinePlotOptions *LPOIm = (LinePlotOptions*) malloc(sizeof(LinePlotOptions) * Spec.NSpectra);
        if((Options.PlotIm==true) && (Options.PlotRe==false))
        {
          LPORe[0] = Options.Style.Line[0];
          LPOIm[0] = Options.Style.Line[0];
          LPORe[0].LineWidth=0;
          LPORe[0].LineColor[3]=0;
          LPORe[0].PointSize=0;
          LPORe[0].PointColor[3]=0;
        }
        else if((Options.PlotIm==false) && (Options.PlotRe==true))
        {
          LPORe[0] = Options.Style.Line[0];
          LPOIm[0] = Options.Style.Line[0];
          LPOIm[0].LineWidth=0;
          LPOIm[0].LineColor[3]=0;
          LPOIm[0].PointSize=0;
          LPOIm[0].PointColor[3]=0;
        }
        else if((Options.PlotIm==true) && (Options.PlotRe==true))
        {
          LPORe[0] = Options.Style.Line[1];
          LPOIm[0] = Options.Style.Line[2];
        }
        else
        {
          return luaL_error(L, "Plot a spectrum but both PlotRe and PlotIm are set to false so there is nothing to do here\n");
        }
        PlotSpectraTypeArrayToGraphics(Graphics, &Spec, LPORe, LPOIm, &(Options.Frame), Options.Negate, Options.WatterfallShift);
        free(LPORe);
        free(LPOIm);
      }
      else
      {
        switch (A->type|32)
        {
          case 'l':
          {
            Options.Nrow    = ((BlockListOfPoles*)(A->data))->BlockSize;
            Options.Ncolumn = ((BlockListOfPoles*)(A->data))->BlockSize;
          }
          break;
          case 't':
          {
            Options.Nrow    = ((BlockTriDiagonalMatrixType*)(A->data))->BlockSize;
            Options.Ncolumn = ((BlockTriDiagonalMatrixType*)(A->data))->BlockSize;
          }
          break;
          case 'a':
          {
            Options.Nrow    = ((BlockAndersonMatrixType*)(A->data))->BlockSize;
            Options.Ncolumn = ((BlockAndersonMatrixType*)(A->data))->BlockSize;         }
          break;
          case 'n':
          {
            Options.Nrow    = ((BlockNaturalImpurityOrbitalMatrixType*)(A->data))->BlockSize;
            Options.Ncolumn = ((BlockNaturalImpurityOrbitalMatrixType*)(A->data))->BlockSize;          }
          break;
        }
        // plot a matrix of spectra
        if ( (Options.Nrow * Options.Ncolumn) != Spec.NSpectra ) return luaL_error(L, "plotting a spectrum array with %d rows and %d columns, but the data array has %d elements\n",Options.Nrow,Options.Ncolumn,Spec.NSpectra);
        LinePlotOptions LPORe, LPOIm;
        if((Options.PlotIm==true) && (Options.PlotRe==false))
        {
          LPORe = Options.Style.Line[0];
          LPOIm = Options.Style.Line[0];
          LPORe.LineWidth=0;
          LPORe.LineColor[3]=0;
          LPORe.PointSize=0;
          LPORe.PointColor[3]=0;
        }
        else if((Options.PlotIm==false) && (Options.PlotRe==true))
        {
          LPORe = Options.Style.Line[0];
          LPOIm = Options.Style.Line[0];
          LPOIm.LineWidth=0;
          LPOIm.LineColor[3]=0;
          LPOIm.PointSize=0;
          LPOIm.PointColor[3]=0;
        }
        else if((Options.PlotIm==true) && (Options.PlotRe==true))
        {
          LPORe = Options.Style.Line[1];
          LPOIm = Options.Style.Line[2];
        }
        else
        {
          return luaL_error(L, "Plot a spectrum but both PlotRe and PlotIm are set to false so there is nothing to do here\n");
        }
        PlotSpectraTypeArrayToGraphicsArray(Graphics, &Spec, Options.Nrow, Options.Ncolumn, &LPORe, &LPOIm, &(Options.Frame),  Options.Negate, Options.ScaleOffDiagonal, Options.ShiftOffDiagonal);
      }
      FreeResponsefunction(A);
      FreeSpectraTypeArray(&Spec);
    }
    else if (luaL_testudata(L, 1, "Matrix_Type") != NULL)
    {
      return luaL_error(L, "Use the function MatrixPlot to plot matrices instead of the function Plot\n");
    }
    else
    {
      return luaL_error(L, "called Graphics.Plot with unkown userdata type for first variable\n");
    }
  }
  else
  {
    return luaL_error(L, "called Graphics.Plot with unkown first variable\n");
  }
  free(Options.Style.Line);
  FreeSpectraTypeArray(&aa);
  return 1;
}

int LuaGraphicsToPS(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len!=1) return luaL_error(L, "Error: Graphics.ToPS should be called with exactly one argument, not with %d\n",len);
  GraphicsType *Graphics = (GraphicsType *) luaL_checkudata(L, 1, "Graphics_Type");
  DynamicString PS;
  GraphicsToPS(Graphics, &PS);
  lua_pushstring(L,PS.str);
  FreeDynamicString(&PS);
  return 1;
}

int LuaGraphicsToSVG(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  bool RelativeSize=false;
  LuaOptionsType GraphicsToSVGOptions[]={
    {"RelativeSize","bool",&(RelativeSize)},
    {NULL, NULL, NULL}};
  if(LuaTestOptionsList(L, len))
  {
    LuaCheckOptions(L, len, GraphicsToSVGOptions);
    len=len-1;
  }
  if(len!=1) return luaL_error(L, "Error: Graphics.ToSVG should be called with exactly one argument, not with %d\n",len);
  GraphicsType *Graphics = (GraphicsType *) luaL_checkudata(L, 1, "Graphics_Type");
  DynamicString SVG;
  GraphicsToSVG(Graphics, &SVG,RelativeSize);
  lua_pushstring(L,SVG.str);
  FreeDynamicString(&SVG);
  return 1;
}


// G1+G2 results in a new graphics with 2 layers containing G1 and G2.
int LuaGraphicsAdd(lua_State * L)
{
  GraphicsType * GraphicsA = (GraphicsType *) luaL_checkudata(L, 1, "Graphics_Type");
  GraphicsType * GraphicsB = (GraphicsType *) luaL_checkudata(L, 2, "Graphics_Type");
  
  GraphicsType * GraphicsAplusB = (GraphicsType*) lua_newuserdata(L, sizeof(GraphicsType));
  luaL_setmetatable(L, "Graphics_Type");
  
  CopyGraphics(GraphicsA, GraphicsAplusB);
  AddLayer(GraphicsAplusB, GraphicsB);
  
  return 1;
}

// returns a table with elements Layers, Line, Curve, Arc, Text and BitMap and values corresponding to the number of elements present in the graph
int LuaGraphicsLen(lua_State * L)
{
  //TODO:
  return 0;
}

// transforms the graphics object to a string
int LuaGraphicsToString(lua_State * L)
{
  lua_pushfstring(L, "Graphics");
  return 1;
}



int LuaGraphicsIndex(lua_State * L)
{
  GraphicsType *Graphics = (GraphicsType *) luaL_checkudata(L, 1, "Graphics_Type");
  const char *str = lua_tostring(L, 2);
  if (strcmp(str, "SizeX") == 0)
  {
    lua_pushnumber(L, Graphics->CanvasSize.X);
    return 1;
  }
  if (strcmp(str, "SizeY") == 0)
  {
    lua_pushnumber(L, Graphics->CanvasSize.Y);
    return 1;
  }
  if (strcmp(str, "OriginX") == 0)
  {
    lua_pushnumber(L, Graphics->CanvasOrigin.X);
    return 1;
  }
  if (strcmp(str, "OriginY") == 0)
  {
    lua_pushnumber(L, Graphics->CanvasOrigin.Y);
    return 1;
  }
  if (strcmp(str, "NLayer") == 0)
  {
    lua_pushinteger(L, Graphics->NLayer);
    return 1;
  }
  if (strcmp(str, "NLine") == 0)
  {
    lua_pushinteger(L, Graphics->NLine);
    return 1;
  }
  if (strcmp(str, "NCurve") == 0)
  {
    lua_pushinteger(L, Graphics->NCurve);
    return 1;
  }
  if (strcmp(str, "NArc") == 0)
  {
    lua_pushinteger(L, Graphics->NArc);
    return 1;
  }
  if (strcmp(str, "NText") == 0)
  {
    lua_pushinteger(L, Graphics->NText);
    return 1;
  }
  if (strcmp(str, "NBitMap") == 0)
  {
    lua_pushinteger(L, Graphics->NBitMap);
    return 1;
  }
  if (strcmp(str, "Add") == 0)
  {
    lua_pushvalue(L, 1);
    lua_pushcclosure(L, LuaAddToGraphicsMethod,1);
    return 1;
  }
  return luaL_error(L, "Error index %s not found in Graphics Object\n",str);
}

int LuaGraphicsNewIndex(lua_State * L)
{
  GraphicsType *Graphics = (GraphicsType *) luaL_checkudata(L, 1, "Graphics_Type");
  const char *str = lua_tostring(L, 2);
  double a;
  if(strcmp(str, "SizeX")==0)
  {
    a = (double)luaL_checknumber(L, 3);
    Graphics->CanvasSize.X = a;
    return 0;
  }
  if(strcmp(str, "SizeY")==0)
  {
    a = (double)luaL_checknumber(L, 3);
    Graphics->CanvasSize.Y = a;
    return 0;
  }
  if(strcmp(str, "OriginX")==0)
  {
    a = (double)luaL_checknumber(L, 3);
    Graphics->CanvasOrigin.X = a;
    return 0;
  }
  if(strcmp(str, "OriginY")==0)
  {
    a = (double)luaL_checknumber(L, 3);
    Graphics->CanvasOrigin.Y = a;
    return 0;
  }
  return luaL_error(L, "Error index %s not found in Graphics Object\n",str);
}

// frees graphics
int LuaGraphicsDestroy(lua_State * L)
{
  GraphicsType *Graphics = (GraphicsType *) luaL_checkudata(L, 1, "Graphics_Type");
  FreeGraphics(Graphics);
  return 0;
}


int LuaNewBitMap(lua_State * L)
{
  int i;
  BitMapType *BitMap;
  uint32_t fillcolor=0xFFFFFFFF;
  uint32_t len = (uint32_t) lua_gettop(L);
  BitMap = (BitMapType*) lua_newuserdata(L, sizeof(BitMapType));
  luaL_setmetatable(L, "BitMap_Type");
  InitBitMap(BitMap);
  if(len==3)
  {
    struct
    {
      double CanvasSizeX;
      double CanvasSizeY;
      double CanvasOriginX;
      double CanvasOriginY;
      int    NpointsX;
      int    NpointsY;
      uint32_t color;
      uint8_t  alpha;
    } Options;
    LuaOptionsType LuaNewGraphicsOptions[]={
      {"Canvas.Size.X","double",&(Options.CanvasSizeX)},
      {"Canvas.Size.Y","double",&(Options.CanvasSizeY)},
      {"Canvas.Origin.X","double",&(Options.CanvasOriginX)},
      {"Canvas.Origin.Y","double",&(Options.CanvasOriginY)},
      {"Color","uint32_t",&(Options.color)},
      {"alpha","uint8_t",&(Options.alpha)},
      {NULL, NULL, NULL}};
    if(LuaTestOptionsList(L, -2))
    {
      Options.CanvasSizeX=BitMap->RT.X-BitMap->LB.X;
      Options.CanvasSizeY=BitMap->RT.Y-BitMap->LB.Y;
      Options.CanvasOriginX=BitMap->LB.X;
      Options.CanvasOriginY=BitMap->LB.Y;
      Options.color=0xFFFFFF;
      Options.alpha=0xFF;
      LuaCheckOptions(L, -2, LuaNewGraphicsOptions);
      fillcolor = ( (Options.color & 0xFFFFFF) + (Options.alpha<<24) );
      BitMap->RT.X=Options.CanvasSizeX;
      BitMap->RT.Y=Options.CanvasSizeY;
      BitMap->LB.X=Options.CanvasOriginX;
      BitMap->LB.Y=Options.CanvasOriginY;
      BitMap->RT.X+=BitMap->LB.X;
      BitMap->RT.Y+=BitMap->LB.Y;
    }
    else
    {
      return luaL_error(L, "Error in LuaNewBitMap BitMap.New called with last argument that has the wrong format for options\n");
    }
  }
  if ( (len==2) || (len=3) )
  {
    BitMap->NpointsX=(uint32_t)luaL_checkinteger(L, 1);
    BitMap->NpointsY=(uint32_t)luaL_checkinteger(L, 2);
    BitMap->Pixels = (uint32_t*)malloc(BitMap->NpointsX*BitMap->NpointsY*sizeof(uint32_t));
    for (i=0;i<BitMap->NpointsX*BitMap->NpointsY;i++)
    {
      BitMap->Pixels[i]=fillcolor;
    }
  }
  else
  {
    return luaL_error(L, "Error in LuaNewBitMap BitMap.New should be called with 2 or 3 argument, not with %d\n",len);
  }
  return 1;
}


int LuaBitMapFromMatrix(lua_State * L)
{
  //TODO:
  return 0;
}


int LuaBitMapOverlay(lua_State * L)
{
  //TODO:
  return 0;
}
int LuaBitMapAdd(lua_State * L)
{
  //TODO:
  return 0;
}

int LuaBitMapLen(lua_State * L)
{
  //TODO:
  return 0;
}

int LuaBitMapToString(lua_State * L)
{
  lua_pushfstring(L, "BitMap");
  return 1;
}

int LuaBitMapSetPixel(lua_State *L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  BitMapType *BitMap = (BitMapType *) luaL_checkudata(L, 1, "BitMap_Type");
  int i = (int) luaL_checkinteger(L, 2);
  if(i>BitMap->NpointsX) return luaL_error(L, "Error: can not index column %d from %d available\n",i, BitMap->NpointsX);
  int j = (int) luaL_checkinteger(L, 3);
  if(j>BitMap->NpointsY) return luaL_error(L, "Error: can not index row %d from %d available\n",j, BitMap->NpointsY);
  uint32_t rgb = (uint32_t) luaL_checkinteger(L, 4);
  BitMap->Pixels[(i-1) + (j-1)*BitMap->NpointsX] = ( (BitMap->Pixels[(i-1) + (j-1) * BitMap->NpointsX]) & 0xFF000000 ) + ( rgb & 0xFFFFFF );
  if(len==5)
  {
    uint8_t alpha = (uint8_t) luaL_checkinteger(L, 5);
    BitMap->Pixels[(i-1) + (j-1)*BitMap->NpointsX] = ( (BitMap->Pixels[(i-1) + (j-1) * BitMap->NpointsX]) & 0xFFFFFF ) + ( alpha << 24 );
  }
  return 0;
}

int LuaBitMapToFile(lua_State * L)
{
  BitMapType *BitMap = (BitMapType *) luaL_checkudata(L, 1, "BitMap_Type");
  const char *filename = lua_tostring(L, 2);
  BitMapToFile(BitMap, filename);
  return 0;
}

int LuaBitMapDestroy(lua_State * L)
{
  BitMapType *BitMap = (BitMapType *) luaL_checkudata(L, 1, "BitMap_Type");
  FreeBitMap(BitMap);
  return 0;
}


/*
 *
 */
static int LuaPlotDOS(lua_State *L)
{
  //TODO:
  
  return 0;
}

/*  The poles a[i] with weights b[i] (b NULL for weight 1) as a spectral function on the NG energies
 *  E0 + i dE: a Lorentzian of FWHM GammaL per pole, or for GammaL 0 the weight divided by dE and
 *  split linearly over the two energies around the pole. Poles outside the grid add nothing.
 */
static void BandPolesToGrid(double const *a, double const *b, int N, double E0, double dE, int NG, double GammaL, double *col)
{
  int i, j, i0;
  double x, f, bb;
  for(i=0;i<NG;i++) col[i]=0;
  for(j=0;j<N;j++)
  {
    bb = (b==NULL) ? (1) : (b[j]);
    if(GammaL>0)
    {
      for(i=0;i<NG;i++) col[i] += bb * (GammaL/2) / (Pi * (Sqr(GammaL/2) + Sqr(E0 + i*dE - a[j])));
    }
    else
    {
      x  = (a[j] - E0) / dE;
      i0 = (int)floor(x);
      f  = x - i0;
      if(i0>=0   && i0<NG)   col[i0]   += bb * (1-f) / dE;
      if(i0+1>=0 && i0+1<NG) col[i0+1] += bb * f / dE;
    }
  }
}

/*  The spectral function -Im R / Pi of the scalar response function R on the NG energies E0 + i dE with
 *  a Lorentzian of FWHM GammaL. A list of poles directly, the other types by ResponseFunctionToSpectra.
 */
static void BandResponseFunctionToGrid(lua_State *L, ResponsefunctionType const *R, double E0, double dE, int NG, double GammaL, double *col)
{
  if(((R->type|32)=='l') && (R->Block==false))
  {
    ListOfPoles const *P = (ListOfPoles*)R->data;
    BandPolesToGrid(P->a, P->b, P->N, E0, dE, NG, GammaL, col);
    return;
  }
  SpectraTypeArray S;
  // one extra energy on both sides, as a spectrum without broadening collects the poles outside on its end points
  ResponseFunctionToSpectra(R, &S, 0, E0-dE, E0+NG*dE, 0, GammaL, NG+1);
  if(S.NSpectra!=1)
  {
    luaL_error(L, "Error PlotBands: the spectral functions per k-point should be scalar response functions, not of dimension %d\n", (int)sqrt((double)S.NSpectra));
  }
  for(int i=0;i<NG;i++) col[i] = -S.Spec[0].Im[i+1] / Pi;
  FreeSpectraTypeArray(&S);
  free(S.Spec);
}

/*  Check that positions 2 to len on the stack are tables with Nktot response functions. If Ymin or Ymax
 *  is not set (infinite) it is determined from the poles of the spectral functions on position 2 with a
 *  weight above 1E-4, widened by Gamma and rounded outward. That needs lists of poles.
 */
static void LuaBandSpectraCheckAndRange(lua_State *L, int len, int Nktot, double *Ymin, double *Ymax, double Gamma)
{
  double lo=INFINITY, hi=-INFINITY;
  for(int i=2;i<=len;i++)
  {
    if(!lua_istable(L, i) || (int)lua_objlen(L, i)!=Nktot)
    {
      luaL_error(L, "Error PlotBands: input at position %d should be a matrix, or a table with one response function per k-point (%d)\n", i, Nktot);
    }
  }
  if(!(*Ymin==-INFINITY || *Ymax==INFINITY)) return;
  for(int ik=0;ik<Nktot;ik++)
  {
    lua_rawgeti(L, 2, ik+1);
    ResponsefunctionType *R = luaL_checkResponsefunction(L, -1);
    if(!(((R->type|32)=='l') && (R->Block==false)))
    {
      luaL_error(L, "Error PlotBands: the energy range is found automatically for lists of poles only, set it with {\"Frame\",{{\"Ymin\",...},{\"Ymax\",...}}}\n");
    }
    ListOfPoles const *P = (ListOfPoles*)R->data;
    for(uint32_t j=0;j<P->N;j++)
    {
      if(P->b[j] < 1E-4) continue;
      if(P->a[j] < lo) lo = P->a[j];
      if(P->a[j] > hi) hi = P->a[j];
    }
    FreeResponsefunctionIncludingOriginalPointer(R);
    lua_pop(L, 1);
  }
  if(lo>hi) luaL_error(L, "Error PlotBands: the spectral functions have no weight, set the energy range with {\"Frame\",{{\"Ymin\",...},{\"Ymax\",...}}}\n");
  if(*Ymin==-INFINITY) *Ymin = floor(lo - Gamma);
  if(*Ymax== INFINITY) *Ymax = ceil(hi + Gamma);
}

/*  The band spectral functions on the pixel rows of BandPlot, drawn by ConvertBandSpectraToBitMap.
 *  On stack 2..len either tables with one response function per k-point (BandEnergies NULL), or the
 *  band energies (BandEnergies) and band characters as matrices, taken as poles of weight 1 and of
 *  weight the character. Position 2 is the total, 3..len the 0, 1 or 3 characters.
 *  Every pixel row holds the spectral function at the energy of its centre. It is broadened by a
 *  Lorentzian of FWHM GammaL, and then by a Gaussian of FWHM GammaG as a sum over the pixel rows,
 *  including the rows within 4 sigma outside the plot. The opacity is Saturation times the total over
 *  the peak height of a pole of weight 1, so a band of weight 1 is drawn in full color.
 */
static void LuaBandSpectraToBitMap(lua_State *L, int len, CompactMatrixType const *BandEnergies, BandstructurePathType *kList, int Nktot, BitMapType *BandPlot, double Emin, double Emax,
                                   double GammaL, double GammaG, double Saturation, uint8_t *LineColor, uint8_t *SecondLineColor, double ZeroDistance)
{
  int const NSet = len-1;
  int const NY = BandPlot->NpointsY;
  int i, j, ip, s, ik, NExt, NG;
  double dE, sigma, E0, sum, Aref, a;
  CompactMatrixType *BandCharacters[3];
  if(!(Emax>Emin)) luaL_error(L, "Error PlotBands: Ymax needs to be larger than Ymin\n");
  if(GammaL<0 || GammaG<0) luaL_error(L, "Error PlotBands: GammaL and GammaG can not be negative\n");
  if(BandEnergies!=NULL)
  {
    for(s=1;s<NSet;s++)
    {
      BandCharacters[s-1]=(CompactMatrixType*) luaL_checkudata(L, 2+s, "Matrix_Type");
      if(BandCharacters[s-1]->Ni!=BandEnergies->Ni || BandCharacters[s-1]->Nj!=BandEnergies->Nj)
      {
        luaL_error(L, "Error PlotBands: the character matrix on position %d does not have the dimensions of the band energies\n", 2+s);
      }
    }
  }
  dE    = (Emax-Emin)/NY;
  sigma = GammaG / (2*sqrt(2*log(2.0)));
  NExt  = (sigma>0) ? ((int)ceil(4*sigma/dE)) : (0);
  NG    = NY + 2*NExt;
  E0    = Emin + (0.5-NExt)*dE;
  double *kern = (double*)malloc((2*NExt+1)*sizeof(double));
  double *col  = (double*)malloc(NG*sizeof(double));
  double *A    = (double*)malloc(NSet*Nktot*NY*sizeof(double));
  sum=0;
  for(j=0;j<=2*NExt;j++)
  {
    kern[j] = (NExt==0) ? (1) : (exp(-Sqr((j-NExt)*dE)/(2*Sqr(sigma))));
    sum += kern[j];
  }
  for(j=0;j<=2*NExt;j++) kern[j]/=sum;
  for(s=0;s<NSet;s++)
  {
    for(ik=0;ik<Nktot;ik++)
    {
      if(BandEnergies!=NULL)
      {
        BandPolesToGrid(&(BandEnergies->VRe[ik*BandEnergies->Nj]), (s==0) ? (NULL) : (&(BandCharacters[s-1]->VRe[ik*BandEnergies->Nj])), BandEnergies->Nj, E0, dE, NG, GammaL, col);
      }
      else
      {
        lua_rawgeti(L, 2+s, ik+1);
        ResponsefunctionType *R = luaL_checkResponsefunction(L, -1);
        BandResponseFunctionToGrid(L, R, E0, dE, NG, GammaL, col);
        FreeResponsefunctionIncludingOriginalPointer(R);
        lua_pop(L, 1);
      }
      // row ip (from the top) holds energy index i = NY-1-ip (from the bottom)
      for(ip=0;ip<NY;ip++)
      {
        i = NY-1-ip;
        sum = 0;
        for(j=0;j<=2*NExt;j++) sum += kern[j] * col[i+j];
        A[(s*Nktot+ik)*NY+ip] = sum;
      }
    }
  }
  // the peak height of a pole of weight 1 on a pixel centre
  a = E0 + NExt*dE;
  BandPolesToGrid(&a, NULL, 1, E0, dE, NG, GammaL, col);
  Aref = 0;
  for(j=0;j<=2*NExt;j++) Aref += kern[j] * col[j];
  ConvertBandSpectraToBitMap(A, NSet, Nktot, kList, BandPlot, LineColor, SecondLineColor, Saturation/Aref, ZeroDistance);
  free(kern);
  free(col);
  free(A);
}

/* On stack we have
 *     1) a path in k-space of the form {{Label, k-point}, N1 , {Label, k-point} , N2, {Label, k-point}}
 *     2) a matrix with band energies
 *     3) a matrix with band characters used as Black-White grey-scale or as band-thickness
 *     3,4)     same as 3 (used as Black-White grey-scale and band-thickness)
 *     3,4,5)   same as 3 (used as RGB band-color)
 *     3,4,5,6) same as 3 (used as RGB band-color and band-thickness)
 *     -1) possible options
 *
 *  In place of the matrices, tables with one response function per k-point can be given (as returned by
 *  TightBinding.CalculateBands with a self energy): 2) the band spectral function and 3) (3,4,5) the
 *  characters. They are drawn broadened by a Lorentzian of FWHM GammaL and a Gaussian of FWHM GammaG,
 *  see LuaBandSpectraToBitMap. With GammaL or GammaG the matrices are drawn the same way, as poles of
 *  weight 1.
 *
 *     we return
 *     1) a plot of the band-structure
 */
static int LuaPlotBands(lua_State *L)
{
  // define options
  const char* ScaleFormat = "%.0f";
  const char* EmptyString = "";
  struct
  {
    int      NPixel[2];
    double   EmptyPathWidth;
    uint32_t LineColor;
    uint8_t  LineTransparancy;
    uint32_t SecondLineColor;
    uint8_t  SecondLineTransparancy;
    double   LineWidth;
    FramePlotOptions Frame;
    bool     Interpolate;
    double   GammaL;
    double   GammaG;
    double   Saturation;
  } Options;
  Options.NPixel[0]=1000;
  Options.NPixel[1]=1000;
  Options.EmptyPathWidth = 0.05;
  Options.LineColor = 0xFF000080;
  Options.LineTransparancy = 0xFF;
  Options.SecondLineColor = 0xFFFF0000;
  Options.SecondLineTransparancy = 0xFF;
  Options.LineWidth = 0.002;
  Options.Frame.PlotRange[0].X=0;
  Options.Frame.PlotRange[0].Y=-INFINITY;
  Options.Frame.PlotRange[1].X=0;
  Options.Frame.PlotRange[1].Y=INFINITY;
  Options.Frame.Margins[0].X=0.08;
  Options.Frame.Margins[0].Y=0.08*sqrt(2.0);
  Options.Frame.Margins[1].X=0.01;
  Options.Frame.Margins[1].Y=0.01*sqrt(2.0);
  Options.Frame.LineWidth=0.001;
  Options.Frame.LineColor[0]=0;
  Options.Frame.LineColor[1]=0;
  Options.Frame.LineColor[2]=0;
  Options.Frame.LineColor[3]=0xFF;
  Options.Frame.FontSize=0.04;
  Options.Frame.Font=TIMES_NEW_ROMAN;
  Options.Frame.XLabel=(char*)EmptyString;
  Options.Frame.YLabel=(char*)EmptyString;
  Options.Frame.FirstXTick=0;
  Options.Frame.dXTick=1.0;
  Options.Frame.NMinorXTick=9;
  Options.Frame.FirstYTick=-INFINITY;
  Options.Frame.dYTick=1.0;
  Options.Frame.NMinorYTick=9;
  Options.Frame.PlotXScaleFormat=(char*)ScaleFormat;
  Options.Frame.PlotYScaleFormat=(char*)ScaleFormat;
  Options.Frame.XAxis=-INFINITY;
  Options.Frame.YAxis=-INFINITY;
  Options.Interpolate = true;
  Options.GammaL = 0;
  Options.GammaG = 0;
  Options.Saturation = 1;
  LuaOptionsType LuaPlotBandOptions[]={
    {"NPixelX","int",&(Options.NPixel[0])},
    {"NPixelY","int",&(Options.NPixel[1])},
    {"Color","uint32_t",&(Options.LineColor)},
    {"Transparancy","uint8_t",&(Options.LineTransparancy)},
    {"SecondColor","uint32_t",&(Options.SecondLineColor)},
    {"SecondTransparancy","uint8_t",&(Options.SecondLineTransparancy)},
    {"LineWidth","double",&(Options.LineWidth)},
    {"EmptyPathWidth","double",&(Options.EmptyPathWidth)},
    {"Frame","FramePlotOptions",&(Options.Frame)},
    {"Interpolate","bool",&(Options.Interpolate)},
    {"GammaL","double",&(Options.GammaL)},
    {"GammaG","double",&(Options.GammaG)},
    {"Saturation","double",&(Options.Saturation)},
    {NULL, NULL, NULL}};
  // check (and read) the input
  int len = (int)lua_gettop(L);
  if (LuaTestOptionsList(L,-1))
  {
    LuaCheckOptions(L,-1,LuaPlotBandOptions);
    Options.LineColor = Options.LineColor & 0x00FFFFFF + (((uint32_t)Options.LineTransparancy)<<24);
    len--;
  }
  if(len<2)
  {
    return luaL_error(L, "Error PlotBands: called with %d arguments (excluding options), minimal 2 expected (k-path, and Band Energies).\n",len);
  }
  if(len>6)
  {
    return luaL_error(L, "Error PlotBands: called with %d arguments (excluding options), maximal 6 expected (k-path, Band Energies, (0 to 4 times) a matrix representing the band character related to color and line width ).\n",len);
  }
  
  // check (and read) the k-list
  BandstructurePathType kList;
  if(!(lua_istable(L, 1)))
  {
    return luaL_error(L, "Error LuaPlotBands: input at position 2 should be a path in k-space of the form {{Label, k-point}, N1 , {Label, k-point} , N2, {Label, k-point}}).\n");
  }
  int lenklist = (int)lua_objlen(L, 1);
  int Nktot = (lenklist+1)/2;
  if(lenklist&&1==0)
  {
    return luaL_error(L, "Error LuaPlotBands: input at position 2 should be a path in k-space of the form {{Label, k-point}, N1 , {Label, k-point} , N2, {Label, k-point}}).\nThe length should be odd and not %d\n",lenklist);
  }
  kList.NQ = (lenklist+1)/2;
  kList.Q  = (BandstructurePointType*)malloc(kList.NQ*sizeof(BandstructurePointType));
  kList.Units = 0; // not used
  for(int i=1;i<=lenklist;i+=2)
  {
    // test k-point
    lua_rawgeti(L, 1, i); // put the table {Label, k-vec} on the stack
    if(!(lua_istable(L, -1)))
    {
      lua_pop(L,1);
      return luaL_error(L, "Wrong format for k-list at position %d\n",i);
    }
    if(!(((int)lua_objlen(L, -1))==2))
    {
      lua_pop(L,1);
      return luaL_error(L, "Wrong format for k-list at position %d\n",i);
    }
    lua_rawgeti(L,-1,1); // put the string Label on the stack
    if(!(lua_type(L, -1)==LUA_TSTRING))
    {
      lua_pop(L,1);
      lua_pop(L,1);
      return luaL_error(L, "Wrong format for k-list at position %d\n",i);
    }
    kList.Q[(i-1)/2].label = lua_tolstring(L, -1, NULL)[0];
    lua_pop(L,1); // remove the string from the stack
    lua_rawgeti(L,-1,2); // put the k-point on the stack
    if(!(lua_istable(L, -1)))
    {
      lua_pop(L,1);
      return luaL_error(L, "Wrong format for k-list at position %d\n",i);
    }
    if(!(((int)lua_objlen(L, -1))==3))
    {
      lua_pop(L,1);
      return luaL_error(L, "Wrong format for k-list at position %d\n",i);
    }
    for(int j=0;j<3;j++)
    {
      lua_rawgeti(L, -1, j+1);
      kList.Q[(i-1)/2].q[j] = (double)luaL_checknumber(L, -1);
      lua_pop(L,1);
    }
    lua_pop(L,1); // remove the k-point from the stack
    lua_pop(L,1); // remove the table {Lable, k-vec} from the stack
    // test integer N
    if(i==lenklist)
    {
      kList.Q[(i-1)/2].Nq = 0;
    }
    else
    {
      lua_rawgeti(L, 1, i+1); // put the integer Ni on the stack
      if (!lua_isnumber(L, -1))
      {
        return luaL_error(L, "Wrong format for k-list at position %d\nExpected an integer\n",i+1);
      }
      kList.Q[(i-1)/2].Nq = (int) lua_tointeger(L, -1);
      if(kList.Q[(i-1)/2].Nq<0)
      {
        return luaL_error(L, "Number of k-points in the path can not be negative at position %d\n",i+1);
      }
      Nktot += kList.Q[(i-1)/2].Nq;
      lua_pop(L,1);
    }
  }
  // finished reading the k-list
  // read the matrix with eigen energies, or a table of spectral functions per k-point
  CompactMatrixType *BandEnergies = (CompactMatrixType*) luaL_testudata(L, 2, "Matrix_Type");
  bool Spectral = (BandEnergies==NULL) || (Options.GammaL>0) || (Options.GammaG>0);
  if(Spectral && !(len==2 || len==3 || len==5))
  {
    return luaL_error(L, "Error PlotBands: spectral functions (or GammaL, GammaG) take 0, 1 or 3 characters, not %d (a band-thickness is not defined)\n", len-2);
  }
  if(BandEnergies==NULL)
  {
    LuaBandSpectraCheckAndRange(L, len, Nktot, &(Options.Frame.PlotRange[0].Y), &(Options.Frame.PlotRange[1].Y), Options.GammaL+Options.GammaG);
  }
  // if range not set determine it automatically
  if(Options.Frame.PlotRange[0].Y == -INFINITY)
  {
    Options.Frame.PlotRange[0].Y = BandEnergies->VRe[0];
    for(int i=1;i<BandEnergies->Ni*BandEnergies->Nj;i++)
    {
      if(BandEnergies->VRe[i] < Options.Frame.PlotRange[0].Y) Options.Frame.PlotRange[0].Y=BandEnergies->VRe[i];
    }
    Options.Frame.PlotRange[0].Y = floor(Options.Frame.PlotRange[0].Y);
  }
  if(Options.Frame.PlotRange[1].Y == INFINITY)
  {
    Options.Frame.PlotRange[1].Y = BandEnergies->VRe[0];
    for(int i=1;i<BandEnergies->Ni*BandEnergies->Nj;i++)
    {
      if(BandEnergies->VRe[i] > Options.Frame.PlotRange[1].Y) Options.Frame.PlotRange[1].Y=BandEnergies->VRe[i];
    }
    Options.Frame.PlotRange[1].Y = ceil(Options.Frame.PlotRange[1].Y);
  }
  if(Options.Frame.FirstYTick==-INFINITY) Options.Frame.FirstYTick = Options.Frame.PlotRange[0].Y;
  // set color
  uint8_t LineColor[4] = {(uint8_t)(Options.LineColor>>16),(uint8_t)(Options.LineColor>>8),(uint8_t)(Options.LineColor),(uint8_t)(Options.LineColor>>24)};
  uint8_t SecondLineColor[4] = {(uint8_t)(Options.SecondLineColor>>16),(uint8_t)(Options.SecondLineColor>>8),(uint8_t)(Options.SecondLineColor),(uint8_t)(Options.SecondLineColor>>24)};
  // allocate the bitmap
  BitMapType BandPlot;
  BandPlot.LB.X=Options.Frame.Margins[0].X;
  BandPlot.LB.Y=Options.Frame.Margins[0].Y;
  BandPlot.RT.X=1-Options.Frame.Margins[1].X;
  BandPlot.RT.Y=1-Options.Frame.Margins[1].Y;
  BandPlot.NpointsX = Options.NPixel[0];
  BandPlot.NpointsY = Options.NPixel[1];
  BandPlot.Pixels = (uint32_t*)malloc(BandPlot.NpointsX*BandPlot.NpointsY*sizeof(uint32_t));
  if(Spectral)
  {
    LuaBandSpectraToBitMap(L, len, BandEnergies, &kList, Nktot, &BandPlot, Options.Frame.PlotRange[0].Y, Options.Frame.PlotRange[1].Y, Options.GammaL, Options.GammaG, Options.Saturation, LineColor, SecondLineColor, Options.EmptyPathWidth);
    GraphicsType* Graphics;
    Graphics = (GraphicsType*) lua_newuserdata(L, sizeof(GraphicsType));
    luaL_setmetatable(L, "Graphics_Type");
    InitGraphics(Graphics);
    AddBitMapInherentMemory(Graphics, &BandPlot);
    PlotAllButDataToBandGraphics(Graphics, &(Options.Frame), &kList, Options.EmptyPathWidth);
    return 1;
  }
  // read character
  CompactMatrixType *BandCharacters[4];
  for(int i=3; i<=len;i++)
  {
    BandCharacters[i-3]=(CompactMatrixType*) luaL_checkudata(L, i, "Matrix_Type");
  }
  
  
  // convert to bitmap depending on number of character and if interpolation is true or not
  if(Options.Interpolate==true)
  {
    switch (len-2)
    {
      case 0:
        // no characters defined
        ConvertBandEnergiesToBitMapInterpolate(BandEnergies, &kList, &BandPlot,Options.Frame.PlotRange[0].Y, Options.Frame.PlotRange[1].Y, Options.LineWidth, LineColor, Options.EmptyPathWidth);
      break;
      case 1:
        // one character used for monochrome color
        ConvertBandEnergiesToBitMapBWInterpolate(BandEnergies, BandCharacters[0], &kList, &BandPlot, Options.Frame.PlotRange[0].Y, Options.Frame.PlotRange[1].Y, Options.LineWidth, LineColor, SecondLineColor, Options.EmptyPathWidth);
      break;
      case 2:
        // two characgters used for monochrome color and thickness
        //TODO:
        ConvertBandEnergiesToBitMapBWInterpolate(BandEnergies, BandCharacters[0], &kList, &BandPlot, Options.Frame.PlotRange[0].Y, Options.Frame.PlotRange[1].Y, Options.LineWidth, LineColor, SecondLineColor, Options.EmptyPathWidth);
      break;
      case 3:
        // three characters used for RGB color
        ConvertBandEnergiesToBitMapRGBInterpolate(BandEnergies, BandCharacters[0], BandCharacters[1], BandCharacters[2], &kList, &BandPlot, Options.Frame.PlotRange[0].Y, Options.Frame.PlotRange[1].Y, Options.LineWidth, Options.EmptyPathWidth);
      break;
      case 4:
        // four characters used for RGB color and thickness
        //TODO:
        ConvertBandEnergiesToBitMapRGBInterpolate(BandEnergies, BandCharacters[0], BandCharacters[1], BandCharacters[2], &kList, &BandPlot, Options.Frame.PlotRange[0].Y, Options.Frame.PlotRange[1].Y, Options.LineWidth, Options.EmptyPathWidth);
      break;
    }
  }
  else
  { // no interpollation
    switch (len-2)
    {
      case 0:
        // no characters defined
        ConvertBandEnergiesToBitMap(BandEnergies, &kList, &BandPlot,Options.Frame.PlotRange[0].Y, Options.Frame.PlotRange[1].Y, Options.LineWidth, LineColor, Options.EmptyPathWidth);
      break;
      case 1:
        // one character used for monochrome color
        ConvertBandEnergiesToBitMapBW(BandEnergies, BandCharacters[0], &kList, &BandPlot, Options.Frame.PlotRange[0].Y, Options.Frame.PlotRange[1].Y, Options.LineWidth, LineColor, SecondLineColor, Options.EmptyPathWidth);
      break;
      case 2:
        // two characgters used for monochrome color and thickness
        //TODO:
        ConvertBandEnergiesToBitMapBW(BandEnergies, BandCharacters[0], &kList, &BandPlot, Options.Frame.PlotRange[0].Y, Options.Frame.PlotRange[1].Y, Options.LineWidth, LineColor, SecondLineColor, Options.EmptyPathWidth);
      break;
      case 3:
        // three characters used for RGB color
        ConvertBandEnergiesToBitMapRGB(BandEnergies, BandCharacters[0], BandCharacters[1], BandCharacters[2], &kList, &BandPlot, Options.Frame.PlotRange[0].Y, Options.Frame.PlotRange[1].Y, Options.LineWidth, Options.EmptyPathWidth);
      break;
      case 4:
        // four characters used for RGB color and thickness
        //TODO:
        ConvertBandEnergiesToBitMapRGB(BandEnergies, BandCharacters[0], BandCharacters[1], BandCharacters[2], &kList, &BandPlot, Options.Frame.PlotRange[0].Y, Options.Frame.PlotRange[1].Y, Options.LineWidth, Options.EmptyPathWidth);
      break;
    }
  }
  // add bitmap to a grahics object and add axis etc
  GraphicsType* Graphics;
  Graphics = (GraphicsType*) lua_newuserdata(L, sizeof(GraphicsType));
  luaL_setmetatable(L, "Graphics_Type");
  InitGraphics(Graphics);
  AddBitMapInherentMemory(Graphics, &BandPlot);
  PlotAllButDataToBandGraphics(Graphics, &(Options.Frame), &kList, Options.EmptyPathWidth);
  return 1;
}

/*
 *
 */
static int LuaPlotBandsAndDOS(lua_State *L)
{
  //TODO:
  
  return 0;
}

/*
 *
 */
static int LuaPlotFermi(lua_State *L)
{
  //TODO:
  
  return 0;
}

/*
 *
 */
static int LuaSVGHTMLViewer(lua_State *L)
{
  lua_pushstring(L, SVGHTMLViewer);
  lua_pushstring(L, SVGSequencyHTMLViewer);
  return 2;
}






static const struct luaL_Reg Graphics_func[] = {
  {"New", LuaNewGraphics},
  {"Add", LuaAddToGraphicsFunction},
  {"Plot", LuaPlotGraphics},
  {"ToPS", LuaGraphicsToPS},
  {"ToSVG", LuaGraphicsToSVG},
  {"PlotDOS", LuaPlotDOS},
  {"PlotBands", LuaPlotBands},
  {"PlotBandsAndDOS", LuaPlotBandsAndDOS},
  {"PlotFermi",LuaPlotFermi},
  {"SVGHTMLViewer",LuaSVGHTMLViewer},
  {NULL, NULL}
};

static const struct luaL_Reg Graphics_meta[] = {
  {"__add", LuaGraphicsAdd},
  {"__len", LuaGraphicsLen},
  {"__tostring", LuaGraphicsToString},
  {"__index", LuaGraphicsIndex},
  {"__newindex", LuaGraphicsNewIndex},
  {"__gc", LuaGraphicsDestroy},
  {NULL, NULL}
};

static const struct luaL_Reg BitMap_func[] = {
  {"New", LuaNewBitMap},
  {"Overlay", LuaBitMapOverlay},
  {"FromMatrix", LuaBitMapFromMatrix},
  {"Save",LuaBitMapToFile},
  {"Pixel", LuaBitMapSetPixel},
  {NULL, NULL}
};

static const struct luaL_Reg BitMap_meta[] = {
  {"__add", LuaBitMapAdd},
  {"__len", LuaBitMapLen},
  {"__tostring", LuaBitMapToString},
  {"__gc", LuaBitMapDestroy},
  {NULL, NULL}
};
