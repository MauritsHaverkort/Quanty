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
      {"Canvas.Origin.X","double",&(Options.CanvasSizeX)},
      {"Canvas.Origin.Y","double",&(Options.CanvasSizeY)},
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
  GraphicsType *Graphics = (GraphicsType *) luaL_checkudata(L, 1, "Graphics_Type");
  struct
  {
    uint32_t RGBalpha;
    uint32_t FillRGBalpha;
    double   LineWidth;
    double   FontSize;
    char     Alignment;
    double   Angle;
    const char* Font;
    bool     Close;
  } Options;
  Options.RGBalpha=0x000000FF;
  Options.FillRGBalpha=0x00000000;
  Options.LineWidth=0.001;
  Options.FontSize=0.05;
  Options.Alignment='l';
  Options.Angle=0;
  Options.Font=TIMES_NEW_ROMAN;
  Options.Close=false;
  LuaOptionsType LuaAddGraphicsOptions[]={
    {"Color","uint32_t",&(Options.RGBalpha)},
    {"Fill Color","uint32_t",&(Options.FillRGBalpha)},
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
      text.Color[0] = (uint8_t)(Options.RGBalpha >> 24);
      text.Color[1] = (uint8_t)(Options.RGBalpha >> 16);
      text.Color[2] = (uint8_t)(Options.RGBalpha >> 8);
      text.Color[3] = (uint8_t)Options.RGBalpha;
      text.Angle    = Options.Angle;
      text.Alignment = Options.Alignment;
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
          lua_pop(L, 1);
          
        }
        case LUA_TTABLE:
        {
          switch (lua_objlen(L, -1))
          {
            case 2:
            {
              lua_pop(L, 1);
              
            }
            case 4:
            {
              lua_pop(L, 1);
              
            }
            default:
            {
              lua_pop(L, 1);
              
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
        GraphicsType *NewLayer = (GraphicsType *) luaL_checkudata(L, 1, "Graphics_Type");
        AddLayer(Graphics,NewLayer);
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
  
  return 0;
}

int LuaGraphicsToPS(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  GraphicsType *Graphics = (GraphicsType *) luaL_checkudata(L, 1, "Graphics_Type");
  DynamicString PS;
  GraphicsToPS(Graphics, &PS);
  lua_pushstring(L,PS.str);
  FreeDynamicString(&PS);
  return 1;
}

int LuaGraphicsToSVG(lua_State * L)
{
  
  return 0;
}


// G1+G2 results in a new graphics with 2 layers containing G1 and G2.
int LuaGraphicsAdd(lua_State * L)
{
  
  return 0;
}

// returns a table with elements Layers, Line, Curve, Arc, Text and BitMap and values corresponding to the number of elements present in the graph
int LuaGraphicsLen(lua_State * L)
{
  
  return 0;
}

// transforms the graphics object to a string
int LuaGraphicsToString(lua_State * L)
{
  lua_pushfstring(L, "Graphics");
  return 1;
}

//returns layer N
int LuaGraphicsIndex(lua_State * L)
{
  
  return 0;
}

// frees graphics
int LuaGraphicsDestroy(lua_State * L)
{
  
  return 0;
}

static const struct luaL_Reg Graphics_func[] = {
  {"New", LuaNewGraphics},
  {"Add", LuaAddToGraphics},
  {"Plot", LuaPlotGraphics},
  {"ToPS", LuaGraphicsToPS},
  {"ToSVG", LuaGraphicsToSVG},
  {NULL, NULL}
};

static const struct luaL_Reg Graphics_meta[] = {
  {"__add", LuaGraphicsAdd},
  {"__len", LuaGraphicsLen},
  {"__tostring", LuaGraphicsToString},
  {"__index", LuaGraphicsIndex},
  {"__gc", LuaGraphicsDestroy},
  {NULL, NULL}
};
