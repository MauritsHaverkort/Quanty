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
      if(len!=4)
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
      return 0;
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
            return 0;
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
                      case 1:
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
                      case 3:
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
                      case 4:
                      {
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
                        BitMap.Pixels[i+j*BitMap.NpointsX] = val;
                        
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
            return 0;
          }
          else
          {
            return luaL_error(L, "Error: Input %s should either be \"Arc\", \"BitMap\" or not present.\n(For a line or Bezier curve, just give a list with 2 or 6 points per element)\n",(char*)lua_tostring(L, -1));
          }
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
              return 0;
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
              return 0;
            }
            default:
            {
              lua_pop(L, 1);
              return luaL_error(L, "Error: elmenet should be either a list of two coordinates (x,y for a line) or a list of 6 coordinates VL, P, VR for a Bezier curve\n");
            }
          }
        }
        default:
        {
          return luaL_error(L, "Error: elmenet should be either {\"Arc\",R,P.x,P.y} or a line or a bezier curve\n");
        }
      }
      return 0;
    }
    case LUA_TUSERDATA:
    {
      if (luaL_testudata(L, len, "Graphics_Type") != NULL)
      {
        GraphicsType *NewLayer = (GraphicsType *) luaL_checkudata(L, len, "Graphics_Type");
        AddLayer(Graphics,NewLayer);
        return 0;
      }
      if (luaL_testudata(L, len, "BitMap_Type") != NULL)
      {
        BitMapType *BitMap = (BitMapType *) luaL_checkudata(L, len, "BitMap_Type");
        AddBitMap(Graphics, BitMap);
        return 0;
      }
      return luaL_error(L, "Error: add incorrect user data to graph\n");
    }
    default:
    {
      return luaL_error(L, "Error: add unrecognized object to graphics\n");
    }
  }
  return 0;
}

int LuaPlotGraphics(lua_State * L)
{
  //TODO: make scientific plot of data
  return 0;
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
  if(len!=1) return luaL_error(L, "Error: Graphics.ToSVG should be called with exactly one argument, not with %d\n",len);
  GraphicsType *Graphics = (GraphicsType *) luaL_checkudata(L, 1, "Graphics_Type");
  DynamicString SVG;
  GraphicsToSVG(Graphics, &SVG);
  lua_pushstring(L,SVG.str);
  FreeDynamicString(&SVG);
  return 1;
}


// G1+G2 results in a new graphics with 2 layers containing G1 and G2.
int LuaGraphicsAdd(lua_State * L)
{
  //TODO: G1+G2
  return 0;
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
  
  
  
  
  unsigned char *str;
  FILE *fp;
  int i,j,k;
  uint32_t filesize = 54 + 4 * BitMap->NpointsX * BitMap->NpointsY;
  str = (unsigned char *)calloc(filesize,sizeof(char));
  str[0]='B';
  str[1]='M';
  str[2] = (unsigned char)(filesize);
  str[3] = (unsigned char)(filesize >> 8);
  str[4] = (unsigned char)(filesize >> 16);
  str[5] = (unsigned char)(filesize >> 24);
  str[10]=54;
  str[14]=40;
  str[18] = (unsigned char)(BitMap->NpointsX);
  str[19] = (unsigned char)(BitMap->NpointsX >> 8);
  str[20] = (unsigned char)(BitMap->NpointsX >> 16);
  str[21] = (unsigned char)(BitMap->NpointsX >> 24);
  str[22] = (unsigned char)(BitMap->NpointsY);
  str[23] = (unsigned char)(BitMap->NpointsY >> 8);
  str[24] = (unsigned char)(BitMap->NpointsY >> 16);
  str[25] = (unsigned char)(BitMap->NpointsY >> 24);
  str[26]=1;
  str[28]=32;
  k=54;
  for(i=BitMap->NpointsY-1;i>=0;i--)
  {
    for(j=0;j<BitMap->NpointsX;j++)
    {
      str[k]=(unsigned char)( BitMap->Pixels[j+i*BitMap->NpointsX] );
      k++;
      str[k]=(unsigned char)( BitMap->Pixels[j+i*BitMap->NpointsX]>>8 );
      k++;
      str[k]=(unsigned char)( BitMap->Pixels[j+i*BitMap->NpointsX]>>16 );
      k++;
      str[k]=(unsigned char)( BitMap->Pixels[j+i*BitMap->NpointsX]>>24 );
      k++;
    }
  }
  fp = fopen( filename , "wb" );
  fwrite(str , 1 , filesize , fp );
  fclose(fp);
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

/*
 *
 */
static int LuaPlotBands(lua_State *L)
{
  //TODO:
  
  return 0;
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






static const struct luaL_Reg Graphics_func[] = {
  {"New", LuaNewGraphics},
  {"Add", LuaAddToGraphics},
  {"Plot", LuaPlotGraphics},
  {"ToPS", LuaGraphicsToPS},
  {"ToSVG", LuaGraphicsToSVG},
  {"PlotDOS", LuaPlotDOS},
  {"PlotBands", LuaPlotBands},
  {"PlotBandsAndDOS", LuaPlotBandsAndDOS},
  {"PlotFermi",LuaPlotFermi},
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
