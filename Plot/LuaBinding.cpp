/*  Graphics objects as lua tables, used by Graphics.ToTable, Graphics.New(Table) and io.put
 *    { CanvasSize = {x, y}, CanvasOrigin = {x, y},
 *      Layer  = { Graphics, ... },   -- layers 1 and up in this same form; layer 0 is the object itself
 *      Line   = { {Color=, Width=, Point=, CloseAndFill=, FillColor=}, ... },
 *      Curve  = { {Color=, Width=, Point=, VL=, VR=, CloseAndFill=, FillColor=}, ... },
 *      Arc    = { {Color=, Width=, R=, Point=, ThetaL=, ThetaR=, CloseAndFill=, FillColor=}, ... },
 *      Text   = { {Color=, Point=, Angle=, Alignment=, Size=, String=, Font=}, ... },
 *      BitMap = { {LB=, RT=, NpointsX=, NpointsY=, Pixels=, Pixelated=, Angle=}, ... } }
 *  with the fields of the structures in Plot/StructDef.cpp. A color is the number 0xAARRGGBB,
 *  AA the transparency (00 opaque, FF fully transparent), a point {x, y} and a list of points
 *  {{x, y}, ...}. Alignment is a string of one character. Pixels is a string in ASCIIHexDecode
 *  format: eight hexadecimal digits per pixel (its uint32_t value 0xAARRGGBB), NpointsX * NpointsY
 *  pixels in the order of the structure, whitespace ignored, closed by '>'.
 *  Writing leaves out empty lists. Reading takes every field as optional, with the defaults of
 *  InitGraphics, InitLine, etc., and refuses keys it does not know.
 */

static void LuaPushGraphicsPoint(lua_State *L, PointType const P)
{
  lua_createtable(L, 2, 0);
  lua_pushnumber(L, P.X);
  lua_rawseti(L, -2, 1);
  lua_pushnumber(L, P.Y);
  lua_rawseti(L, -2, 2);
}

static void LuaPushGraphicsPoints(lua_State *L, int const N, PointType const * const P)
{
  int i;
  lua_createtable(L, N, 0);
  for(i=0;i<N;i++)
  {
    LuaPushGraphicsPoint(L, P[i]);
    lua_rawseti(L, -2, i+1);
  }
}

static void LuaPushGraphicsColor(lua_State *L, uint32_t const C)
{
  lua_pushinteger(L, (lua_Integer)C);
}

void LuaPushTable_Graphics(lua_State *L, GraphicsType const * const G)
{
  int i, j;
  lua_createtable(L, 0, 8);
  LuaPushGraphicsPoint(L, G->CanvasSize);
  lua_setfield(L, -2, "CanvasSize");
  LuaPushGraphicsPoint(L, G->CanvasOrigin);
  lua_setfield(L, -2, "CanvasOrigin");
  if(G->NLayer>1)
  {
    lua_createtable(L, G->NLayer-1, 0);
    for(i=1;i<G->NLayer;i++)
    {
      LuaPushTable_Graphics(L, G->Layer[i]);
      lua_rawseti(L, -2, i);
    }
    lua_setfield(L, -2, "Layer");
  }
  if(G->NLine>0)
  {
    lua_createtable(L, G->NLine, 0);
    for(i=0;i<G->NLine;i++)
    {
      LineType const * const E = &(G->Line[i]);
      lua_createtable(L, 0, 5);
      LuaPushGraphicsColor(L, E->Color);                   lua_setfield(L, -2, "Color");
      lua_pushnumber(L, E->Width);                         lua_setfield(L, -2, "Width");
      LuaPushGraphicsPoints(L, E->NPoints, E->Point);      lua_setfield(L, -2, "Point");
      lua_pushboolean(L, E->CloseAndFill);                 lua_setfield(L, -2, "CloseAndFill");
      LuaPushGraphicsColor(L, E->FillColor);               lua_setfield(L, -2, "FillColor");
      lua_rawseti(L, -2, i+1);
    }
    lua_setfield(L, -2, "Line");
  }
  if(G->NCurve>0)
  {
    lua_createtable(L, G->NCurve, 0);
    for(i=0;i<G->NCurve;i++)
    {
      CurveType const * const E = &(G->Curve[i]);
      lua_createtable(L, 0, 7);
      LuaPushGraphicsColor(L, E->Color);                   lua_setfield(L, -2, "Color");
      lua_pushnumber(L, E->Width);                         lua_setfield(L, -2, "Width");
      LuaPushGraphicsPoints(L, E->NPoints, E->Point);      lua_setfield(L, -2, "Point");
      LuaPushGraphicsPoints(L, E->NPoints, E->VL);         lua_setfield(L, -2, "VL");
      LuaPushGraphicsPoints(L, E->NPoints, E->VR);         lua_setfield(L, -2, "VR");
      lua_pushboolean(L, E->CloseAndFill);                 lua_setfield(L, -2, "CloseAndFill");
      LuaPushGraphicsColor(L, E->FillColor);               lua_setfield(L, -2, "FillColor");
      lua_rawseti(L, -2, i+1);
    }
    lua_setfield(L, -2, "Curve");
  }
  if(G->NArc>0)
  {
    lua_createtable(L, G->NArc, 0);
    for(i=0;i<G->NArc;i++)
    {
      ArcType const * const E = &(G->Arc[i]);
      lua_createtable(L, 0, 8);
      LuaPushGraphicsColor(L, E->Color);                   lua_setfield(L, -2, "Color");
      lua_pushnumber(L, E->Width);                         lua_setfield(L, -2, "Width");
      lua_pushnumber(L, E->R);                             lua_setfield(L, -2, "R");
      LuaPushGraphicsPoint(L, E->Point);                   lua_setfield(L, -2, "Point");
      lua_pushnumber(L, E->ThetaL);                        lua_setfield(L, -2, "ThetaL");
      lua_pushnumber(L, E->ThetaR);                        lua_setfield(L, -2, "ThetaR");
      lua_pushboolean(L, E->CloseAndFill);                 lua_setfield(L, -2, "CloseAndFill");
      LuaPushGraphicsColor(L, E->FillColor);               lua_setfield(L, -2, "FillColor");
      lua_rawseti(L, -2, i+1);
    }
    lua_setfield(L, -2, "Arc");
  }
  if(G->NText>0)
  {
    lua_createtable(L, G->NText, 0);
    for(i=0;i<G->NText;i++)
    {
      TextType const * const E = &(G->Text[i]);
      lua_createtable(L, 0, 7);
      LuaPushGraphicsColor(L, E->Color);                   lua_setfield(L, -2, "Color");
      LuaPushGraphicsPoint(L, E->Point);                   lua_setfield(L, -2, "Point");
      lua_pushnumber(L, E->Angle);                         lua_setfield(L, -2, "Angle");
      lua_pushlstring(L, &(E->Alignment), 1);              lua_setfield(L, -2, "Alignment");
      lua_pushnumber(L, E->Size);                          lua_setfield(L, -2, "Size");
      lua_pushstring(L, (E->String==NULL)?"":E->String);   lua_setfield(L, -2, "String");
      if(E->Font!=NULL)
      {
        lua_pushstring(L, E->Font);                        lua_setfield(L, -2, "Font");
      }
      lua_rawseti(L, -2, i+1);
    }
    lua_setfield(L, -2, "Text");
  }
  if(G->NBitMap>0)
  {
    lua_createtable(L, G->NBitMap, 0);
    for(i=0;i<G->NBitMap;i++)
    {
      BitMapType const * const E = &(G->BitMap[i]);
      size_t const N = (size_t)E->NpointsX * (size_t)E->NpointsY;
      // eight digits per pixel, a line break after every 64 pixels, and the closing '>'
      char *Hex = (char*)malloc((8*N + N/64 + 2) * sizeof(char));
      size_t pos = 0;
      if(Hex==NULL) luaL_error(L, "Error in Graphics.ToTable: out of memory for a bitmap of %d by %d pixels\n", E->NpointsX, E->NpointsY);
      for(j=0;j<(int)N;j++)
      {
        snprintf(Hex+pos, 9, "%08X", E->Pixels[j]);
        pos += 8;
        if((j%64)==63) Hex[pos++] = '\n';
      }
      Hex[pos++] = '>';
      lua_createtable(L, 0, 5);
      LuaPushGraphicsPoint(L, E->LB);                      lua_setfield(L, -2, "LB");
      LuaPushGraphicsPoint(L, E->RT);                      lua_setfield(L, -2, "RT");
      lua_pushinteger(L, E->NpointsX);                     lua_setfield(L, -2, "NpointsX");
      lua_pushinteger(L, E->NpointsY);                     lua_setfield(L, -2, "NpointsY");
      lua_pushlstring(L, Hex, pos);                        lua_setfield(L, -2, "Pixels");
      if(E->Pixelated) { lua_pushboolean(L, true);         lua_setfield(L, -2, "Pixelated"); }
      if(E->Angle!=0)  { lua_pushnumber(L, E->Angle);      lua_setfield(L, -2, "Angle"); }
      free(Hex);
      lua_rawseti(L, -2, i+1);
    }
    lua_setfield(L, -2, "BitMap");
  }
}

/*  every key of the table at narg must be one of Keys, What names the table in the message
 */
static void LuaCheckGraphicsKeys(lua_State *L, int narg, const char * const Keys[], int NKeys, const char *What)
{
  int i;
  narg = lua_absindex(L, narg);
  if(lua_type(L, narg)!=LUA_TTABLE) luaL_error(L, "Error in Graphics.New(Table): %s should be a table\n", What);
  lua_pushnil(L);
  while(lua_next(L, narg)!=0)
  {
    bool known=false;
    if(lua_type(L, -2)==LUA_TSTRING)
    {
      for(i=0;i<NKeys;i++) if(strcmp(lua_tostring(L, -2), Keys[i])==0) known=true;
    }
    if(!known)
    {
      if(lua_type(L, -2)==LUA_TSTRING) luaL_error(L, "Error in Graphics.New(Table): unknown field %s in %s\n", lua_tostring(L, -2), What);
      luaL_error(L, "Error in Graphics.New(Table): a key of type %s in %s, expected named fields\n", luaL_typename(L, -2), What);
    }
    lua_pop(L, 1);
  }
}

static PointType LuaReadGraphicsPoint(lua_State *L, int narg, const char *What)
{
  PointType P;
  narg = lua_absindex(L, narg);
  if((lua_type(L, narg)!=LUA_TTABLE)||(lua_objlen(L, narg)!=2)) luaL_error(L, "Error in Graphics.New(Table): %s should be a point {x, y}\n", What);
  lua_rawgeti(L, narg, 1);
  if(lua_type(L, -1)!=LUA_TNUMBER) luaL_error(L, "Error in Graphics.New(Table): %s should be a point {x, y} of two numbers\n", What);
  P.X = lua_tonumber(L, -1);
  lua_pop(L, 1);
  lua_rawgeti(L, narg, 2);
  if(lua_type(L, -1)!=LUA_TNUMBER) luaL_error(L, "Error in Graphics.New(Table): %s should be a point {x, y} of two numbers\n", What);
  P.Y = lua_tonumber(L, -1);
  lua_pop(L, 1);
  return P;
}

/*  a list of points, allocated here, NULL for an empty list
 */
static PointType* LuaReadGraphicsPoints(lua_State *L, int narg, int *N, const char *What)
{
  int i;
  PointType *P;
  narg = lua_absindex(L, narg);
  if(lua_type(L, narg)!=LUA_TTABLE) luaL_error(L, "Error in Graphics.New(Table): %s should be a list of points {{x, y}, ...}\n", What);
  *N = (int)lua_objlen(L, narg);
  if(*N==0) return NULL;
  P = (PointType*)malloc((*N)*sizeof(PointType));
  if(P==NULL) luaL_error(L, "Error in Graphics.New(Table): out of memory for %d points\n", *N);
  for(i=0;i<*N;i++)
  {
    lua_rawgeti(L, narg, i+1);
    P[i] = LuaReadGraphicsPoint(L, -1, What);
    lua_pop(L, 1);
  }
  return P;
}

static void LuaReadGraphicsColor(lua_State *L, int narg, uint32_t *C, const char *What)
{
  double c;
  narg = lua_absindex(L, narg);
  if(lua_type(L, narg)!=LUA_TNUMBER) luaL_error(L, "Error in Graphics.New(Table): %s should be a color 0xAARRGGBB\n", What);
  c = lua_tonumber(L, narg);
  if((c!=floor(c))||(c<0)||(c>4294967295.0)) luaL_error(L, "Error in Graphics.New(Table): %s should be a color 0xAARRGGBB, a whole number from 0 to 0xFFFFFFFF, found %f\n", What, c);
  *C = (uint32_t)c;
}

static double LuaReadGraphicsNumber(lua_State *L, int narg, const char *What)
{
  if(lua_type(L, narg)!=LUA_TNUMBER) luaL_error(L, "Error in Graphics.New(Table): %s should be a number\n", What);
  return lua_tonumber(L, narg);
}

static bool LuaReadGraphicsBool(lua_State *L, int narg, const char *What)
{
  if(lua_type(L, narg)!=LUA_TBOOLEAN) luaL_error(L, "Error in Graphics.New(Table): %s should be true or false\n", What);
  return (bool)lua_toboolean(L, narg);
}

/*  The fields of an element: lua_getfield pushes the value, the macro reads it if present and pops it
 */
#define GRAPHICS_FIELD(narg, Key, Read) { lua_getfield(L, narg, Key); if(!lua_isnil(L, -1)) { Read; } lua_pop(L, 1); }

/*  fills G, initialised by InitGraphics, from the table at narg
 */
void LuaReadTable_Graphics(lua_State *L, int narg, GraphicsType *G)
{
  const char * const GraphicsKeys[] = {"CanvasSize", "CanvasOrigin", "Layer", "Line", "Curve", "Arc", "Text", "BitMap"};
  const char * const LineKeys[]     = {"Color", "Width", "Point", "CloseAndFill", "FillColor"};
  const char * const CurveKeys[]    = {"Color", "Width", "Point", "VL", "VR", "CloseAndFill", "FillColor"};
  const char * const ArcKeys[]      = {"Color", "Width", "R", "Point", "ThetaL", "ThetaR", "CloseAndFill", "FillColor"};
  const char * const TextKeys[]     = {"Color", "Point", "Angle", "Alignment", "Size", "String", "Font"};
  const char * const BitMapKeys[]   = {"LB", "RT", "NpointsX", "NpointsY", "Pixels", "Pixelated", "Angle"};
  int i, n, e;
  narg = lua_absindex(L, narg);
  LuaCheckGraphicsKeys(L, narg, GraphicsKeys, 8, "the graphics");
  GRAPHICS_FIELD(narg, "CanvasSize",   G->CanvasSize   = LuaReadGraphicsPoint(L, -1, "CanvasSize"));
  GRAPHICS_FIELD(narg, "CanvasOrigin", G->CanvasOrigin = LuaReadGraphicsPoint(L, -1, "CanvasOrigin"));
  // the elements: each list is read element by element, and each element is added to G as soon as
  // it is complete, so that G is a valid object whenever an error stops the reading
  lua_getfield(L, narg, "Line");
  if(!lua_isnil(L, -1))
  {
    if(lua_type(L, -1)!=LUA_TTABLE) luaL_error(L, "Error in Graphics.New(Table): Line should be a list of lines\n");
    n = (int)lua_objlen(L, -1);
    for(i=0;i<n;i++)
    {
      LineType E;
      lua_rawgeti(L, -1, i+1);
      e = lua_gettop(L);
      LuaCheckGraphicsKeys(L, e, LineKeys, 5, "a line");
      InitLine(&E);
      GRAPHICS_FIELD(e, "Color",        LuaReadGraphicsColor(L, -1, &E.Color, "the Color of a line"));
      GRAPHICS_FIELD(e, "Width",        E.Width = LuaReadGraphicsNumber(L, -1, "the Width of a line"));
      GRAPHICS_FIELD(e, "CloseAndFill", E.CloseAndFill = LuaReadGraphicsBool(L, -1, "CloseAndFill of a line"));
      GRAPHICS_FIELD(e, "FillColor",    LuaReadGraphicsColor(L, -1, &E.FillColor, "the FillColor of a line"));
      GRAPHICS_FIELD(e, "Point",        E.Point = LuaReadGraphicsPoints(L, -1, &(E.NPoints), "the Point of a line"));
      if(AddLineInherentLineMemory(G, &E)) luaL_error(L, "Error in Graphics.New(Table): out of memory\n");
      lua_pop(L, 1);
    }
  }
  lua_pop(L, 1);
  lua_getfield(L, narg, "Curve");
  if(!lua_isnil(L, -1))
  {
    if(lua_type(L, -1)!=LUA_TTABLE) luaL_error(L, "Error in Graphics.New(Table): Curve should be a list of curves\n");
    n = (int)lua_objlen(L, -1);
    for(i=0;i<n;i++)
    {
      CurveType E;
      int NVL=0, NVR=0;
      PointType *VL=NULL, *VR=NULL;
      lua_rawgeti(L, -1, i+1);
      e = lua_gettop(L);
      LuaCheckGraphicsKeys(L, e, CurveKeys, 7, "a curve");
      InitCurve(&E);
      GRAPHICS_FIELD(e, "Color",        LuaReadGraphicsColor(L, -1, &E.Color, "the Color of a curve"));
      GRAPHICS_FIELD(e, "Width",        E.Width = LuaReadGraphicsNumber(L, -1, "the Width of a curve"));
      GRAPHICS_FIELD(e, "CloseAndFill", E.CloseAndFill = LuaReadGraphicsBool(L, -1, "CloseAndFill of a curve"));
      GRAPHICS_FIELD(e, "FillColor",    LuaReadGraphicsColor(L, -1, &E.FillColor, "the FillColor of a curve"));
      GRAPHICS_FIELD(e, "Point",        E.Point = LuaReadGraphicsPoints(L, -1, &(E.NPoints), "the Point of a curve"));
      GRAPHICS_FIELD(e, "VL",           VL = LuaReadGraphicsPoints(L, -1, &NVL, "the VL of a curve"));
      GRAPHICS_FIELD(e, "VR",           VR = LuaReadGraphicsPoints(L, -1, &NVR, "the VR of a curve"));
      E.VL = VL;
      E.VR = VR;
      if((NVL!=E.NPoints)||(NVR!=E.NPoints))
      {
        int const NP = E.NPoints;
        FreeCurve(&E);
        luaL_error(L, "Error in Graphics.New(Table): a curve has %d points, but %d VL and %d VR, the three lists should have the same length\n", NP, NVL, NVR);
      }
      if(AddCurveInherentMemory(G, &E)) luaL_error(L, "Error in Graphics.New(Table): out of memory\n");
      lua_pop(L, 1);
    }
  }
  lua_pop(L, 1);
  lua_getfield(L, narg, "Arc");
  if(!lua_isnil(L, -1))
  {
    if(lua_type(L, -1)!=LUA_TTABLE) luaL_error(L, "Error in Graphics.New(Table): Arc should be a list of arcs\n");
    n = (int)lua_objlen(L, -1);
    for(i=0;i<n;i++)
    {
      ArcType E;
      lua_rawgeti(L, -1, i+1);
      e = lua_gettop(L);
      LuaCheckGraphicsKeys(L, e, ArcKeys, 8, "an arc");
      InitArc(&E);
      GRAPHICS_FIELD(e, "Color",        LuaReadGraphicsColor(L, -1, &E.Color, "the Color of an arc"));
      GRAPHICS_FIELD(e, "Width",        E.Width = LuaReadGraphicsNumber(L, -1, "the Width of an arc"));
      GRAPHICS_FIELD(e, "R",            E.R = LuaReadGraphicsNumber(L, -1, "the R of an arc"));
      GRAPHICS_FIELD(e, "Point",        E.Point = LuaReadGraphicsPoint(L, -1, "the Point of an arc"));
      GRAPHICS_FIELD(e, "ThetaL",       E.ThetaL = LuaReadGraphicsNumber(L, -1, "the ThetaL of an arc"));
      GRAPHICS_FIELD(e, "ThetaR",       E.ThetaR = LuaReadGraphicsNumber(L, -1, "the ThetaR of an arc"));
      GRAPHICS_FIELD(e, "CloseAndFill", E.CloseAndFill = LuaReadGraphicsBool(L, -1, "CloseAndFill of an arc"));
      GRAPHICS_FIELD(e, "FillColor",    LuaReadGraphicsColor(L, -1, &E.FillColor, "the FillColor of an arc"));
      if(AddArc(G, &E)) luaL_error(L, "Error in Graphics.New(Table): out of memory\n");
      lua_pop(L, 1);
    }
  }
  lua_pop(L, 1);
  lua_getfield(L, narg, "Text");
  if(!lua_isnil(L, -1))
  {
    if(lua_type(L, -1)!=LUA_TTABLE) luaL_error(L, "Error in Graphics.New(Table): Text should be a list of texts\n");
    n = (int)lua_objlen(L, -1);
    for(i=0;i<n;i++)
    {
      TextType E;
      lua_rawgeti(L, -1, i+1);
      e = lua_gettop(L);
      LuaCheckGraphicsKeys(L, e, TextKeys, 7, "a text");
      InitText(&E);
      GRAPHICS_FIELD(e, "Color",        LuaReadGraphicsColor(L, -1, &E.Color, "the Color of a text"));
      GRAPHICS_FIELD(e, "Point",        E.Point = LuaReadGraphicsPoint(L, -1, "the Point of a text"));
      GRAPHICS_FIELD(e, "Angle",        E.Angle = LuaReadGraphicsNumber(L, -1, "the Angle of a text"));
      GRAPHICS_FIELD(e, "Size",         E.Size = LuaReadGraphicsNumber(L, -1, "the Size of a text"));
      lua_getfield(L, e, "Alignment");
      if(!lua_isnil(L, -1))
      {
        size_t len;
        const char *a;
        if(lua_type(L, -1)!=LUA_TSTRING) luaL_error(L, "Error in Graphics.New(Table): the Alignment of a text should be a string of one character\n");
        a = lua_tolstring(L, -1, &len);
        if(len!=1) luaL_error(L, "Error in Graphics.New(Table): the Alignment of a text should be a string of one character, found \"%s\"\n", a);
        E.Alignment = a[0];
      }
      lua_pop(L, 1);
      // AddText copies the string and keeps the font through InternFont, so both may point into
      // lua strings; they stay on the stack until the text is added
      lua_getfield(L, e, "String");
      if(lua_isnil(L, -1))
      {
        E.String = (char*)"";
      }
      else
      {
        if(lua_type(L, -1)!=LUA_TSTRING) luaL_error(L, "Error in Graphics.New(Table): the String of a text should be a string\n");
        E.String = (char*)lua_tostring(L, -1);
      }
      lua_getfield(L, e, "Font");
      if(!lua_isnil(L, -1))
      {
        if(lua_type(L, -1)!=LUA_TSTRING) luaL_error(L, "Error in Graphics.New(Table): the Font of a text should be a string\n");
        E.Font = lua_tostring(L, -1);
      }
      if(AddText(G, &E)) luaL_error(L, "Error in Graphics.New(Table): out of memory\n");
      lua_pop(L, 3);
    }
  }
  lua_pop(L, 1);
  lua_getfield(L, narg, "BitMap");
  if(!lua_isnil(L, -1))
  {
    if(lua_type(L, -1)!=LUA_TTABLE) luaL_error(L, "Error in Graphics.New(Table): BitMap should be a list of bitmaps\n");
    n = (int)lua_objlen(L, -1);
    for(i=0;i<n;i++)
    {
      BitMapType E;
      double Nx=0, Ny=0;
      size_t len=0, k, N, j=0;
      const char *Hex="";
      lua_rawgeti(L, -1, i+1);
      e = lua_gettop(L);
      LuaCheckGraphicsKeys(L, e, BitMapKeys, 7, "a bitmap");
      InitBitMap(&E);
      GRAPHICS_FIELD(e, "LB",       E.LB = LuaReadGraphicsPoint(L, -1, "the LB of a bitmap"));
      GRAPHICS_FIELD(e, "RT",       E.RT = LuaReadGraphicsPoint(L, -1, "the RT of a bitmap"));
      GRAPHICS_FIELD(e, "NpointsX", Nx = LuaReadGraphicsNumber(L, -1, "the NpointsX of a bitmap"));
      GRAPHICS_FIELD(e, "NpointsY", Ny = LuaReadGraphicsNumber(L, -1, "the NpointsY of a bitmap"));
      GRAPHICS_FIELD(e, "Pixelated", E.Pixelated = LuaReadGraphicsBool(L, -1, "Pixelated of a bitmap"));
      GRAPHICS_FIELD(e, "Angle",    E.Angle = LuaReadGraphicsNumber(L, -1, "the Angle of a bitmap"));
      if((Nx!=floor(Nx))||(Nx<0)||(Ny!=floor(Ny))||(Ny<0)) luaL_error(L, "Error in Graphics.New(Table): NpointsX and NpointsY of a bitmap should be whole numbers of 0 or more\n");
      E.NpointsX = (int)Nx;
      E.NpointsY = (int)Ny;
      N = (size_t)E.NpointsX * (size_t)E.NpointsY;
      lua_getfield(L, e, "Pixels");
      if(!lua_isnil(L, -1))
      {
        if(lua_type(L, -1)!=LUA_TSTRING) luaL_error(L, "Error in Graphics.New(Table): the Pixels of a bitmap should be a string of hexadecimal digits\n");
        Hex = lua_tolstring(L, -1, &len);
      }
      // ASCIIHexDecode: hexadecimal digits, whitespace ignored, '>' ends the data
      E.Pixels = (N==0)?NULL:(uint32_t*)malloc(N*sizeof(uint32_t));
      if((N!=0)&&(E.Pixels==NULL)) luaL_error(L, "Error in Graphics.New(Table): out of memory for a bitmap of %d by %d pixels\n", E.NpointsX, E.NpointsY);
      for(k=0;k<N;k++) E.Pixels[k]=0;
      for(k=0;(k<len)&&(Hex[k]!='>');k++)
      {
        char const c = Hex[k];
        uint32_t d;
        if((c==' ')||(c=='\n')||(c=='\r')||(c=='\t')||(c=='\f')) continue;
        if((c>='0')&&(c<='9'))      d = (uint32_t)(c-'0');
        else if((c>='A')&&(c<='F')) d = (uint32_t)(c-'A'+10);
        else if((c>='a')&&(c<='f')) d = (uint32_t)(c-'a'+10);
        else
        {
          FreeBitMap(&E);
          return (void)luaL_error(L, "Error in Graphics.New(Table): the Pixels of a bitmap contain '%c', which is not a hexadecimal digit\n", c);
        }
        if(j>=8*N)
        {
          FreeBitMap(&E);
          return (void)luaL_error(L, "Error in Graphics.New(Table): the Pixels of a bitmap of %d by %d contain more than %f hexadecimal digits\n", E.NpointsX, E.NpointsY, (double)(8*N));
        }
        E.Pixels[j/8] = (E.Pixels[j/8]<<4) | d;
        j++;
      }
      if(j!=8*N)
      {
        FreeBitMap(&E);
        luaL_error(L, "Error in Graphics.New(Table): the Pixels of a bitmap of %d by %d should have %f hexadecimal digits (8 per pixel), found %f\n", E.NpointsX, E.NpointsY, (double)(8*N), (double)j);
      }
      if(AddBitMapInherentMemory(G, &E)) luaL_error(L, "Error in Graphics.New(Table): out of memory\n");
      lua_pop(L, 2);
    }
  }
  lua_pop(L, 1);
  lua_getfield(L, narg, "Layer");
  if(!lua_isnil(L, -1))
  {
    if(lua_type(L, -1)!=LUA_TTABLE) luaL_error(L, "Error in Graphics.New(Table): Layer should be a list of graphics\n");
    n = (int)lua_objlen(L, -1);
    for(i=0;i<n;i++)
    {
      GraphicsType *Layer = NewGraphics();
      lua_rawgeti(L, -1, i+1);
      // added before it is read, so that it is freed with G when an error stops the reading
      if(AddLayerInherentMemory(G, Layer)) luaL_error(L, "Error in Graphics.New(Table): out of memory\n");
      LuaReadTable_Graphics(L, -1, Layer);
      lua_pop(L, 1);
    }
  }
  lua_pop(L, 1);
}
#undef GRAPHICS_FIELD

/*  true when the table at narg has named fields, i.e. is of the form of Graphics.ToTable and not a
 *  list of options {{"Canvas.Size.X", x}, ...}
 */
static bool LuaQGraphicsTable(lua_State *L, int narg)
{
  bool named=false;
  narg = lua_absindex(L, narg);
  if(lua_type(L, narg)!=LUA_TTABLE) return false;
  lua_pushnil(L);
  while(lua_next(L, narg)!=0)
  {
    if(lua_type(L, -2)==LUA_TSTRING) named=true;
    lua_pop(L, 1);
  }
  return named;
}

int LuaGraphicsToTable(lua_State * L)
{
  if(lua_gettop(L)!=1) return luaL_error(L, "Error Graphics.ToTable expects one argument, a graphics object\n");
  GraphicsType *Graphics = (GraphicsType *) luaL_checkudata(L, 1, "Graphics_Type");
  LuaPushTable_Graphics(L, Graphics);
  return 1;
}

int LuaNewGraphics(lua_State * L)
{
  GraphicsType *Graphics;
  uint32_t len = (uint32_t) lua_gettop(L);
  Graphics = (GraphicsType*) lua_newuserdata(L, sizeof(GraphicsType));
  luaL_setmetatable(L, "Graphics_Type");
  InitGraphics(Graphics);
  if((len==1)&&LuaQGraphicsTable(L, 1))
  {
    // Graphics.New(Table) with a table of the form of Graphics.ToTable
    LuaReadTable_Graphics(L, 1, Graphics);
    return 1;
  }
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

/*  G.Rotate(Angle [, Options]) turns the graphics G, or one of its layers with all layers in it, by
 *  Angle (radians, counterclockwise) about a point. Everything keeps its size and shape in mm; the
 *  canvas of each layer becomes the bounding box of the rotated canvas. Options:
 *    {"Layer", n} or {"Layer", {n1, n2, ...}}  the layer to rotate, as for Graphics.Add; 0 (the
 *                                               default) is the whole graphics
 *    {"Center", {x, y}}                         the point to rotate about, in fractions of that
 *                                               layer (default {0.5, 0.5}, its middle)
 *  On the stack: 1 the graphics, 2 the angle, 3 the options (if given).
 */
static GraphicsType *LuaSelectLayer(lua_State *L, GraphicsType *G, CompactVectorType *Layer, const char *Caller);

int LuaRotateGraphics(lua_State *L)
{
  int len = lua_gettop(L);
  GraphicsType *Graphics = (GraphicsType *) luaL_checkudata(L, 1, "Graphics_Type");
  struct
  {
    CompactVectorType Layer;
    CompactVectorType Center;
  } Options;
  Options.Layer.Ni=0;   Options.Layer.QComplex=0;   Options.Layer.VRe=NULL;
  Options.Center.Ni=0;  Options.Center.QComplex=0;  Options.Center.VRe=NULL;
  LuaOptionsType LuaRotateGraphicsOptions[]={
    {"Layer","Vector",&(Options.Layer)},
    {"Center","Vector",&(Options.Center)},
    {NULL, NULL, NULL}};
  if((len==3)&&LuaTestOptionsList(L, 3))
  {
    LuaCheckOptions(L, 3, LuaRotateGraphicsOptions);
    len--;
  }
  if(len!=2) return luaL_error(L, "Error in Graphics.Rotate: expected a graphics, an angle and possibly options, not %d arguments\n", lua_gettop(L));
  if(lua_type(L, 2)!=LUA_TNUMBER) return luaL_error(L, "Error in Graphics.Rotate: the angle should be a number (radians, counterclockwise)\n");
  double const Angle = lua_tonumber(L, 2);
  if(!isfinite(Angle)) return luaL_error(L, "Error in Graphics.Rotate: the angle should be a finite number\n");
  PointType Center = {0.5, 0.5};
  if(Options.Center.Ni!=0)
  {
    bool const ok = (Options.Center.Ni==2) && (!Options.Center.QComplex);
    if(ok) { Center.X = Options.Center.VRe[0]; Center.Y = Options.Center.VRe[1]; }
    free(Options.Center.VRe);
    if(!ok)
    {
      free(Options.Layer.VRe);
      return luaL_error(L, "Error in Graphics.Rotate: option Center should be a point {x, y} of two real numbers, in fractions of the layer\n");
    }
  }
  Graphics = LuaSelectLayer(L, Graphics, &(Options.Layer), "Graphics.Rotate");
  // the centre in mm, on the page
  Center.X = Graphics->CanvasOrigin.X + Center.X * Graphics->CanvasSize.X;
  Center.Y = Graphics->CanvasOrigin.Y + Center.Y * Graphics->CanvasSize.Y;
  if(RotateGraphics(Graphics, Angle, Center)) return luaL_error(L, "Error in Graphics.Rotate: a layer has a canvas of size 0 (or less) and can not be rotated\n");
  return 0;
}

// G.Rotate(Angle, Options) rotates G itself
int LuaRotateGraphicsMethod(lua_State *L)
{
  lua_pushvalue(L, lua_upvalueindex(1));
  lua_insert(L, 1);
  LuaRotateGraphics(L);
  return 0;
}

// Graphics.Rotate(G, Angle, Options) returns a rotated copy of G
int LuaRotateGraphicsFunction(lua_State *L)
{
  GraphicsType *GraphicsOld = (GraphicsType *) luaL_checkudata(L, 1, "Graphics_Type");
  GraphicsType *Graphics = (GraphicsType*) lua_newuserdata(L, sizeof(GraphicsType));
  luaL_setmetatable(L, "Graphics_Type");
  CopyGraphics(GraphicsOld, Graphics);
  lua_replace(L, 1);
  LuaRotateGraphics(L);
  lua_pushvalue(L, 1);
  return 1;
}

// A Lua colour function, kept in the registry under the key Ref, that returns red, green and
// blue between 0 and 1. After the first failing call Status is not LUA_OK, the error message
// is on the stack and further calls do nothing.
typedef struct LuaColorFunctionContext
{
  lua_State *L;
  int        Ref;
  int        Status;
  const char *Caller; // the name of the Lua function, for error messages
} LuaColorFunctionContext;

// calls the colour function with the NArg numbers in Arg
static void LuaCallColorFunction(LuaColorFunctionContext * const C, int NArg, double const *Arg, double RGB[3])
{
  if(C->Status!=LUA_OK) return;
  lua_rawgeti(C->L, LUA_REGISTRYINDEX, C->Ref);
  for(int a=0;a<NArg;a++) lua_pushnumber(C->L, Arg[a]);
  C->Status = lua_pcall(C->L, NArg, 3, 0);
  if(C->Status!=LUA_OK) return;
  for(int c=0;c<3;c++)
  {
    int isnum;
    RGB[c] = (double)lua_tonumberx(C->L, c-3, &isnum);
    if(!isnum)
    {
      lua_pop(C->L, 3);
      lua_pushfstring(C->L, "%s: the ColorFunction should return three numbers, red, green and blue between 0 and 1", C->Caller);
      C->Status = LUA_ERRRUN;
      return;
    }
  }
  lua_pop(C->L, 3);
}

// The colour of a colour map made by Graphics.Plot with a ColorFunction f(Re, Im)
static void LuaColorFunction(double Re, double Im, void *Context, double RGB[3])
{
  double const Arg[2] = {Re, Im};
  LuaCallColorFunction((LuaColorFunctionContext*)Context, 2, Arg, RGB);
}

// One pixel of a bitmap for Graphics.Add, at narg: grey, {grey}, {r,g,b} or {r,g,b,t},
// each from 0 to 255 and clipped to that range, t the transparency (0 opaque, 255 fully
// transparent, 0 if left out). Returns 0xAARRGGBB with AA = t.
static uint32_t LuaReadBitMapPixel(lua_State *L, int narg)
{
  lua_Integer c[4] = {0, 0, 0, 0}; // red, green, blue, transparency; signed, so that negative values are clipped to 0
  narg = lua_absindex(L, narg);
  if(lua_type(L, narg)==LUA_TTABLE)
  {
    int const n = (int)lua_objlen(L, narg);
    if(!((n==1)||(n==3)||(n==4))) return (uint32_t)luaL_error(L, "ERROR; A pixel in a bitmap should have length 1 (BW), 3 (RGB) or 4 (RGBalpha)\n");
    for(int k=0;k<n;k++)
    {
      lua_rawgeti(L, narg, k+1);
      c[k] = luaL_checkinteger(L, -1);
      lua_pop(L, 1);
    }
    if(n==1) c[1] = c[2] = c[0];
  }
  else
  {
    c[0] = c[1] = c[2] = luaL_checkinteger(L, narg);
  }
  uint32_t Pixel = 0;
  int const Order[4] = {3, 0, 1, 2}; // transparency, red, green, blue
  for(int k=0;k<4;k++)
  {
    lua_Integer v = c[Order[k]];
    if(v<0) v=0;
    if(v>255) v=255;
    Pixel = (Pixel<<8) + (uint32_t)v;
  }
  return Pixel;
}

// The option Layer of Graphics.Add and Graphics.Rotate: n is layer n of G (0 G itself, 1 the first
// layer added to it, -1 the last one), a list {n1, n2, ...} layer n2 of layer n1 and so on, and an
// empty list G itself. Frees the vector Layer, also before an error.
static GraphicsType *LuaSelectLayer(lua_State *L, GraphicsType *G, CompactVectorType *Layer, const char *Caller)
{
  if(Layer->Ni==0) return G;
  for(int i=0;i<(int)Layer->Ni;i++)
  {
    double v = Layer->QComplex ? Layer->V[i].Re : Layer->VRe[i];
    if((Layer->QComplex)||(v!=floor(v)))
    {
      free(Layer->VRe); // VRe and V share their memory
      Layer->Ni=0;
      luaL_error(L, "Error in %s: option Layer should be a whole number or a list of whole numbers\n", Caller);
    }
    // the range is checked on the double, a cast of a double out of the range of int is undefined
    double n = (v<0) ? v + G->NLayer : v;
    if(!((n>=0)&&(n<G->NLayer)))
    {
      char Number[32]; // luaL_error knows no %g
      snprintf(Number, sizeof(Number), "%g", v);
      free(Layer->VRe);
      Layer->Ni=0;
      luaL_error(L, "Error in %s: option Layer, entry %d: there is no layer %s, the layers are 0 (the graphics itself) to %d, or -1 (the last one) to -%d\n",
                 Caller, i+1, Number, G->NLayer-1, G->NLayer);
    }
    G = G->Layer[(int)n];
  }
  free(Layer->VRe);
  Layer->Ni=0;
  return G;
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
  GraphicsType *Graphics = (GraphicsType *) luaL_checkudata(L, 1, "Graphics_Type");
  struct
  {
    uint32_t RGB;       // 0xAARRGGBB, AA the transparency
    uint8_t  alpha;     // an opacity, 255 opaque (kept for old scripts)
    uint8_t  Transparency; // 0 opaque, 255 fully transparent
    uint32_t FillRGB;
    uint8_t  Fillalpha;
    uint8_t  FillTransparency;
    double   LineWidth;
    double   FontSize;
    const char* Alignment;
    double   Angle;
    const char* Font;
    bool     Close;
    uint8_t  Pixelated; // 2 if not given, else 0 or 1 as LuaCheckOptions writes a bool
    CompactVectorType Layer; // the layer to add to, Ni=0 for the graphics itself
  } Options;
  static_assert(sizeof(bool)==1, "the option Pixelated of Graphics.Add is read as a bool into a uint8_t");
  char defaultAlignment = 'l';
  Options.RGB=0x00000000;     // opaque black
  Options.FillRGB=0x00000000;
  Options.alpha=0xFF;
  Options.Fillalpha=0xFF;
  Options.Transparency=0;
  Options.FillTransparency=0;
  Options.LineWidth=0.001;
  Options.FontSize=0.05;
  Options.Alignment=&defaultAlignment;
  Options.Angle=0;
  Options.Font=TIMES_NEW_ROMAN;
  Options.Close=false;
  Options.Pixelated=2;
  Options.Layer.Ni=0;
  Options.Layer.QComplex=0;
  Options.Layer.VRe=NULL;
  LuaOptionsType LuaAddGraphicsOptions[]={
    {"alpha","uint8_t",&(Options.alpha)},
    {"Transparancy","uint8_t",&(Options.Transparency)},
    {"Fill alpha","uint8_t",&(Options.Fillalpha)},
    {"Fill Transparancy","uint8_t",&(Options.FillTransparency)},
    {"Color","uint32_t",&(Options.RGB)},
    {"Fill Color","uint32_t",&(Options.FillRGB)},
    {"Close","bool",&(Options.Close)},
    {"Line Width","double",&(Options.LineWidth)},
    {"Size","double",&(Options.FontSize)},
    {"Alignment","char*",&(Options.Alignment)},
    {"Angle","double",&(Options.Angle)},
    {"Font","char*",&(Options.Font)},
    {"Pixelated","bool",&(Options.Pixelated)},
    {"Layer","Vector",&(Options.Layer)},
    {NULL, NULL, NULL}};
  if(LuaTestOptionsList(L, -1))
  {
    len=len-1;
    LuaCheckOptions(L, -1, LuaAddGraphicsOptions);
  }
  // the option Layer: n adds to layer n (0 the graphics itself, 1 the first layer added to it,
  // -1 the last one, drawn on top), a list {n1, n2, ...} to layer n2 of layer n1 and so on
  Graphics = LuaSelectLayer(L, Graphics, &(Options.Layer), "Graphics.Add");
  // the transparency AA of the colour, the option Transparancy and the opacity alpha together
  Options.RGB     = ColorCombine(Options.RGB,     Options.Transparency,     Options.alpha);
  Options.FillRGB = ColorCombine(Options.FillRGB, Options.FillTransparency, Options.Fillalpha);
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
        text.Color = Options.RGB;
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
              Arc.Color = Options.RGB;
              Arc.Width        = Options.LineWidth;
              Arc.CloseAndFill = Options.Close;
              Arc.FillColor = Options.FillRGB;
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
              // Format {"BitMap", {px,py}, {px,py}, Pixels}
              //     or {"BitMap", {px,py}, {px,py}, Pixels, TransparentColor}
              //     or {"BitMap", {px,py}, {px,py}, ColorFunction, Matrix1, Matrix2, ...}
              // Pixels is a table of rows, the first row at the top, of pixels grey, {grey}, {r,g,b} or {r,g,b,alpha}, 0 to 255.
              // The pixels with the red, green and blue of TransparentColor (given as a pixel) become transparent,
              // or get the alpha of TransparentColor if it has one.
              // ColorFunction(Re M1[i][j], Im M1[i][j], Re M2[i][j], Im M2[i][j], ...) gives the red, green and blue,
              // 0 to 1, of the pixel in row i and column j; these bitmaps are drawn with hard pixels.
              const char * const Formats = "{\"BitMap\", {px,py}, {px,py}, Pixels}, {\"BitMap\", {px,py}, {px,py}, Pixels, TransparentColor} or {\"BitMap\", {px,py}, {px,py}, ColorFunction, Matrix, ...}";
              BitMapType BitMap;
              InitBitMap(&BitMap);
              lua_pop(L, 1);
              if(Tablelen<4) return luaL_error(L, "Error: Input for BitMap should be %s\n", Formats);
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
                if(Tablelen>5) return luaL_error(L, "Error: Input for BitMap should be %s\n", Formats);
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
                    BitMap.Pixels[i+j*BitMap.NpointsX] = LuaReadBitMapPixel(L, -1);
                    lua_pop(L, 1);
                  }
                  lua_pop(L, 1);
                }
                lua_pop(L, 1);
                if(Tablelen==5)
                {
                  lua_rawgeti(L, 2, 5);
                  uint32_t const Key   = LuaReadBitMapPixel(L, -1);
                  // the transparency of TransparentColor if it has one, else fully transparent
                  uint32_t const Transparency = ((lua_type(L,-1)==LUA_TTABLE)&&(lua_objlen(L,-1)==4)) ? (Key & 0xFF000000) : 0xFF000000;
                  lua_pop(L, 1);
                  for(int k=0;k<BitMap.NpointsX*BitMap.NpointsY;k++)
                  {
                    if(((BitMap.Pixels[k] ^ Key) & 0x00FFFFFF) == 0) BitMap.Pixels[k] = (BitMap.Pixels[k] & 0x00FFFFFF) | Transparency;
                  }
                }
              }
              else if(lua_type(L,-1) == LUA_TFUNCTION)
              {
                int const NM = (int)Tablelen - 4;
                if(NM<1) return luaL_error(L, "Error: a BitMap with a ColorFunction needs at least one matrix; input for BitMap should be %s\n", Formats);
                LuaColorFunctionContext Context;
                Context.L = L;
                Context.Status = LUA_OK;
                Context.Caller = "Graphics.Add";
                Context.Ref = luaL_ref(L, LUA_REGISTRYINDEX); // pops the function
                CompactMatrixType **M = (CompactMatrixType**)malloc(NM*sizeof(CompactMatrixType*));
                for(int m=0;m<NM;m++)
                {
                  lua_rawgeti(L, 2, 5+m);
                  M[m] = luaL_checkCompactMatrix(L, lua_gettop(L));
                  lua_pop(L, 1); // a Matrix object stays alive, the table at position 2 holds it
                  if((M[m]->Ni!=M[0]->Ni)||(M[m]->Nj!=M[0]->Nj)) return luaL_error(L, "Error: the matrices of a BitMap should all have the same size, matrix 1 is %d by %d and matrix %d is %d by %d\n", (int)M[0]->Ni, (int)M[0]->Nj, m+1, (int)M[m]->Ni, (int)M[m]->Nj);
                }
                if((M[0]->Ni==0)||(M[0]->Nj==0)) return luaL_error(L, "Error: the matrices of a BitMap should not be empty\n");
                BitMap.NpointsY = (int)M[0]->Ni;
                BitMap.NpointsX = (int)M[0]->Nj;
                BitMap.Pixels = (uint32_t*)malloc((size_t)BitMap.NpointsX * BitMap.NpointsY * sizeof(uint32_t));
                BitMap.Pixelated = true;
                double *Arg = (double*)malloc(2*NM*sizeof(double));
                for(j=0;(j<BitMap.NpointsY)&&(Context.Status==LUA_OK);j++)
                {
                  for(i=0;i<BitMap.NpointsX;i++)
                  {
                    for(int m=0;m<NM;m++)
                    {
                      size_t const k = (size_t)j * M[m]->Nj + i;
                      if(M[m]->QComplex) { Arg[2*m] = M[m]->V[k].Re; Arg[2*m+1] = M[m]->V[k].Im; }
                      else               { Arg[2*m] = M[m]->VRe[k];  Arg[2*m+1] = 0; }
                    }
                    double RGB[3] = {1,1,1};
                    LuaCallColorFunction(&Context, 2*NM, Arg, RGB);
                    BitMap.Pixels[i+j*BitMap.NpointsX] = RGBToPixel(RGB);
                  }
                }
                free(Arg);
                for(int m=0;m<NM;m++)
                {
                  if(!M[m]->OnLuaStack)
                  {
                    FreeCompactMatrix(M[m]);
                    free(M[m]);
                  }
                }
                free(M);
                luaL_unref(L, LUA_REGISTRYINDEX, Context.Ref);
                if(Context.Status!=LUA_OK)
                {
                  // the error message of the ColorFunction is on the top of the stack
                  free(BitMap.Pixels);
                  return lua_error(L);
                }
              }
              else
              {
                return luaL_error(L, "Error: the fourth element of a BitMap should be a table of pixels or a ColorFunction; input for BitMap should be %s\n", Formats);
              }
              if(Options.Pixelated!=2) BitMap.Pixelated = (Options.Pixelated==1);
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
                Line.Color = Options.RGB;
                Line.FillColor = Options.FillRGB;
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
                Curve.Color = Options.RGB;
                Curve.FillColor = Options.FillRGB;
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

// Graphics.Plot draws spectrum Start, Start+Step, Start+2*Step, ... up to Stop
// (counting from 1 and including Stop) out of N spectra. Returns how many these are.
static int NumberOfSelectedSpectra(int N, int Start, int Stop, int Step)
{
  if(Stop>N) Stop=N;
  if(Stop<Start) return 0;
  return (Stop-Start)/Step+1;
}

// Graphics.Plot: is the value at narg a list of points {{x1,y1}, {x2,y2}, ...}? A list of curves
// has Spectra objects, response functions or lists of points as elements; a response function
// in table form has a metatable.
static bool LuaIsPointList(lua_State *L, int narg)
{
  bool Is = false;
  narg = lua_absindex(L, narg);
  if(lua_type(L, narg)!=LUA_TTABLE) return false;
  if(lua_getmetatable(L, narg))
  {
    lua_pop(L, 1);
    return false;
  }
  lua_rawgeti(L, narg, 1);
  if((lua_type(L, -1)==LUA_TTABLE)&&(lua_objlen(L, -1)==2))
  {
    lua_rawgeti(L, -1, 1);
    Is = (lua_type(L, -1)==LUA_TNUMBER);
    lua_pop(L, 1);
  }
  lua_pop(L, 1);
  return Is;
}

// Graphics.Plot: reads the list of points {{x1,y1}, {x2,y2}, ...} at narg, with y a real or a
// complex number, into the spectrum S (the y values in Re and Im, Emin and Emax the smallest and
// the largest x) and its energies *X. The points may lie on any mesh and are drawn in their order.
// Curve is the number of the curve, for error messages.
static void LuaReadPointList(lua_State *L, int narg, SpectraType *S, double **X, int Curve)
{
  narg = lua_absindex(L, narg);
  int const N = (int)lua_objlen(L, narg);
  memset(S, 0, sizeof(SpectraType));
  S->ResonantEnergy = NAN;
  S->Emin =  INFINITY;
  S->Emax = -INFINITY;
  S->N  = N-1;
  S->Re = (double*)malloc(N*sizeof(double));
  S->Im = (double*)malloc(N*sizeof(double));
  *X    = (double*)malloc(N*sizeof(double));
  for(int i=0;i<N;i++)
  {
    lua_rawgeti(L, narg, i+1);
    if((lua_type(L, -1)!=LUA_TTABLE)||(lua_objlen(L, -1)!=2)) luaL_error(L, "Graphics.Plot: point %d of curve %d should be {x, y}\n", i+1, Curve);
    lua_rawgeti(L, -1, 1);
    if(lua_type(L, -1)!=LUA_TNUMBER) luaL_error(L, "Graphics.Plot: the x of point %d of curve %d should be a number\n", i+1, Curve);
    (*X)[i] = (double)lua_tonumber(L, -1);
    lua_pop(L, 1);
    lua_rawgeti(L, -1, 2);
    if(lua_type(L, -1)==LUA_TNUMBER)
    {
      S->Re[i] = (double)lua_tonumber(L, -1);
      S->Im[i] = 0;
    }
    else if(luaL_testudata(L, -1, "Complex_Type")!=NULL)
    {
      Complex16 const * const z = (Complex16*)lua_touserdata(L, -1);
      S->Re[i] = z->Re;
      S->Im[i] = z->Im;
    }
    else luaL_error(L, "Graphics.Plot: the y of point %d of curve %d should be a real or a complex number\n", i+1, Curve);
    lua_pop(L, 2);
    if((*X)[i]<S->Emin) S->Emin = (*X)[i];
    if((*X)[i]>S->Emax) S->Emax = (*X)[i];
  }
}

// The colour of a colour map made by Graphics.Plot without a ColorFunction:
// white to RGB as Saturation times the value (Im, or Re) goes from 0 to 1
typedef struct BuiltInColorContext
{
  double RGB[3];
  double Saturation;
  bool   UseRe;
} BuiltInColorContext;

static void BuiltInColorFunction(double Re, double Im, void *Context, double RGB[3])
{
  BuiltInColorContext const * const C = (BuiltInColorContext*)Context;
  double v = C->Saturation * (C->UseRe ? Re : Im);
  if(!(v>0)) v=0;
  if(v>1) v=1;
  for(int c=0;c<3;c++) RGB[c] = 1 - v * (1 - C->RGB[c]);
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
  SpectraTypeArray *a,aa,View;
  aa.NSpectra=0;
  aa.Spec=NULL;
  aa.Name=NULL;
  a=&aa;
  View.Spec=NULL; // the selected spectra, sharing their data with the spectra they are taken from
  double *YPos=NULL; // the heights of the spectra in a waterfall plot or a colour map
  double **PointX=NULL; // for curves given as lists of points their x values, NULL for spectra on an even mesh
  int NPointX=0;
  bool AutoYRange=false; // the y range was chosen from the data, the y ticks then lie on multiples of dYTick
  uint32_t len = (uint32_t) lua_gettop(L);
  Graphics = (GraphicsType*) lua_newuserdata(L, sizeof(GraphicsType));
  luaL_setmetatable(L, "Graphics_Type");
  InitGraphics(Graphics);
  struct
  {
    double   ScaleOffDiagonal;
    double   ShiftOffDiagonal;
    double   Scale;
    double   WaterfallOffset;
    double   WaterfallShift;
    
    int      SpectraStart;
    int      SpectraStop;
    int      SpectraStep;
    int      Nrow;
    int      Ncolumn;
    int      PlotPoints;
    
    FramePlotOptions Frame;
    LinePlotOptionsArray Style;
    int      LineStyleCycle; // a list of spectra cycles through styles 1 to LineStyleCycle, 0 for all styles
    
    bool     PlotIm;
    bool     PlotRe;
    bool     Negate;
    bool     UseResonantEnergy;
    
    bool     ColorMap;
    bool     Interpolate;
    int      NPixelX;
    int      NPixelY;
    int      ColorFunction; // registry key of a Lua function, LUA_NOREF if not given
    uint32_t Color;
    double   Saturation;
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
  Options.Frame.LineColor = 0x00000000;
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
  Options.Frame.XTicks=true;
  Options.Frame.YTicks=true;
  Options.Style.NAllocLine=3;
  Options.Style.Line = (LinePlotOptions*)malloc(sizeof(LinePlotOptions)*Options.Style.NAllocLine);
  Options.Style.Line[0].LineWidth=0.002;
  Options.Style.Line[0].LineColor = 0x00000000;
  Options.Style.Line[0].PointSize=0;
  Options.Style.Line[0].PointColor = 0x00000000;
  Options.Style.Line[1].LineWidth=0.002;
  Options.Style.Line[1].LineColor = 0x00FF0000;
  Options.Style.Line[1].PointSize=0;
  Options.Style.Line[1].PointColor = 0x00FF0000;
  Options.Style.Line[2].LineWidth=0.002;
  Options.Style.Line[2].LineColor = 0x000000FF;
  Options.Style.Line[2].PointSize=0;
  Options.Style.Line[2].PointColor = 0x000000FF;
  Options.LineStyleCycle=0;
  
  Options.PlotIm=true;
  Options.PlotRe=true;
  Options.Negate=false;
  Options.UseResonantEnergy=false;
  Options.ColorMap=false;
  Options.Interpolate=true;
  Options.NPixelX=600;
  Options.NPixelY=600;
  Options.ColorFunction=LUA_NOREF;
  Options.Color=0x00000000; // opaque black
  Options.Saturation=1;
  Options.Nrow=1;
  Options.Ncolumn=1;
  Options.PlotPoints=600;
  Options.ScaleOffDiagonal=1;
  Options.ShiftOffDiagonal=0;
  Options.Scale=1;
  Options.WaterfallOffset=0;
  Options.WaterfallShift=0;
  Options.SpectraStart=1;
  Options.SpectraStop=INT_MAX;
  Options.SpectraStep=1;
  LuaOptionsType LuaGraphicsOptions[]={
    {"Frame","FramePlotOptions",&(Options.Frame)},
    {"Axis", "FramePlotOptions",&(Options.Frame)},
    {"LineStyle", "LinePlotOptionsArray",&(Options.Style)},
    {"LineStyleCycle","int",&(Options.LineStyleCycle)},
    {"plotIm","bool",&(Options.PlotIm)},
    {"plotRe","bool",&(Options.PlotRe)},
    {"Negate","bool",&(Options.Negate)},
    {"Nrow","int",&(Options.Nrow)},
    {"Ncolumn","int",&(Options.Ncolumn)},
    {"PlotPoints","int",&(Options.PlotPoints)},
    {"ScaleOffDiagonal","double",&(Options.ScaleOffDiagonal)},
    {"ShiftOffDiagonal","double",&(Options.ShiftOffDiagonal)},
    {"Scale","double",&(Options.Scale)},
    {"WaterfallOffset","double",&(Options.WaterfallOffset)},
    {"WaterfallShift","double",&(Options.WaterfallShift)},
    {"WatterfallShift","double",&(Options.WaterfallShift)}, // old spelling, kept for existing scripts
    {"SpectraStart","int",&(Options.SpectraStart)},
    {"SpectraStop","int",&(Options.SpectraStop)},
    {"SpectraStep","int",&(Options.SpectraStep)},
    {"UseResonantEnergy","bool",&(Options.UseResonantEnergy)},
    {"ColorMap","bool",&(Options.ColorMap)},
    {"Interpolate","bool",&(Options.Interpolate)},
    {"NPixelX","int",&(Options.NPixelX)},
    {"NPixelY","int",&(Options.NPixelY)},
    {"ColorFunction","function",&(Options.ColorFunction)},
    {"Color","uint32_t",&(Options.Color)},
    {"Saturation","double",&(Options.Saturation)},
    {NULL, NULL, NULL}};
  if(LuaTestOptionsList(L, len))
  {
    LuaCheckOptions(L, len, LuaGraphicsOptions);
    len=len-1;
  }
  if((Options.SpectraStart<1)||(Options.SpectraStep<1)) return luaL_error(L, "Graphics.Plot: SpectraStart and SpectraStep should be 1 or larger, not %d and %d\n",Options.SpectraStart,Options.SpectraStep);
  // the styles are always allocated (at least 3, the single spectrum needs styles 2 and 3), the cycle may be shorter
  if(Options.LineStyleCycle==0) Options.LineStyleCycle = Options.Style.NAllocLine;
  if((Options.LineStyleCycle<1)||(Options.LineStyleCycle>Options.Style.NAllocLine)) return luaL_error(L, "Graphics.Plot: LineStyleCycle should be between 1 and the number of styles, %d, not %d\n",Options.Style.NAllocLine,Options.LineStyleCycle);
  if(Options.UseResonantEnergy)
  {
    // spectrum i is drawn at WaterfallOffset + its resonant energy, which only the spectra in a Spectra object have
    if(luaL_testudata(L, 1, "Spectra_Type") == NULL) return luaL_error(L, "Graphics.Plot: UseResonantEnergy needs a Spectra object as first argument, as made by CreateResonantSpectra\n");
    if(Options.WaterfallShift!=0) return luaL_error(L, "Graphics.Plot: UseResonantEnergy places each spectrum at WaterfallOffset + its resonant energy, it can not be combined with WaterfallShift\n");
    if((Options.Nrow!=1)||(Options.Ncolumn!=1)) return luaL_error(L, "Graphics.Plot: UseResonantEnergy makes a waterfall plot, it can not be combined with Nrow and Ncolumn\n");
  }
  if(Options.ColorMap)
  {
    if((Options.Nrow!=1)||(Options.Ncolumn!=1)) return luaL_error(L, "Graphics.Plot: ColorMap draws all spectra in one graph, it can not be combined with Nrow and Ncolumn\n");
    if((Options.NPixelX<1)||(Options.NPixelY<1)) return luaL_error(L, "Graphics.Plot: NPixelX and NPixelY should be 1 or larger, not %d and %d\n",Options.NPixelX,Options.NPixelY);
  }
  bool const SelectSpectra = (Options.SpectraStart!=1)||(Options.SpectraStep!=1)||(Options.SpectraStop!=INT_MAX);
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
    // a single list of points {{x1,y1}, ...} is plotted as a list with one curve
    bool const OnePointList = LuaIsPointList(L, 1);
    int const NTable = OnePointList ? 1 : (int)lua_objlen(L, 1);
    a->NSpectra = NumberOfSelectedSpectra(NTable, Options.SpectraStart, Options.SpectraStop, Options.SpectraStep);
    if(a->NSpectra==0) return luaL_error(L, "Graphics.Plot: no spectra to plot, the table has %d elements and SpectraStart is %d, SpectraStop %d\n",NTable,Options.SpectraStart,(Options.SpectraStop<NTable)?Options.SpectraStop:NTable);
    a->Spec = (SpectraType*)malloc(a->NSpectra*sizeof(SpectraType));
    a->Name=NULL;
    for(int k=0;k<a->NSpectra;k++)
    {
      // read element i of the table and see if we can convert it to a spectrum
      // only the selected elements are read, spectrum k of the plot is element i of the table
      i = Options.SpectraStart - 1 + k * Options.SpectraStep;
      if(OnePointList) lua_pushvalue(L, 1);
      else             lua_rawgeti(L, 1, i+1);
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
            free(Spec.Name);
            a->Spec[k]=Spec.Spec[0];
            FreeResponsefunctionIncludingOriginalPointer(A);
          }
          else
          {
            lua_pop(L, 2);
            return luaL_error(L, "Found a table to plot, but on position %d there is an object I can not convert to a spectrum\n",i+1);
          }
        }
        else if(LuaIsPointList(L, -1))
        {
          // a list of points on any mesh
          if(PointX==NULL)
          {
            NPointX = (int)a->NSpectra;
            PointX = (double**)calloc(NPointX, sizeof(double*));
          }
          LuaReadPointList(L, -1, &(a->Spec[k]), &(PointX[k]), i+1);
        }
        else // we found a table without meta data. We can not convert this to a spectra
        {
          return luaL_error(L, "Found a table to plot, but on position %d there is an object I can not convert to a spectrum or a list of points {{x1,y1}, {x2,y2}, ...}\n",i+1);
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
          SpectrumAllocate(Spec->Spec[0], &(a->Spec[k]));
          CopySpectra(Spec->Spec[0], &(a->Spec[k]));
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
          a->Spec[k]=Spec.Spec[0];
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
      labelSelectSpectra:
      if(SelectSpectra)
      {
        // copy the headers of the selected spectra, not their data
        View.NSpectra = NumberOfSelectedSpectra(a->NSpectra, Options.SpectraStart, Options.SpectraStop, Options.SpectraStep);
        if(View.NSpectra==0) return luaL_error(L, "Graphics.Plot: no spectra to plot, there are %d spectra and SpectraStart is %d, SpectraStop %d\n",a->NSpectra,Options.SpectraStart,(Options.SpectraStop<(int)a->NSpectra)?Options.SpectraStop:(int)a->NSpectra);
        View.Spec = (SpectraType*)malloc(View.NSpectra*sizeof(SpectraType));
        for(int k=0;k<View.NSpectra;k++) View.Spec[k] = a->Spec[Options.SpectraStart - 1 + k * Options.SpectraStep];
        View.Name = NULL;
        a = &View;
      }
      labelSpectraTypeArrayPlot:
      if((PointX!=NULL)&&Options.ColorMap) return luaL_error(L, "Graphics.Plot: a colour map needs spectra on an even energy mesh, not lists of points\n");
      if((PointX!=NULL)&&((Options.Nrow!=1)||(Options.Ncolumn!=1))) return luaL_error(L, "Graphics.Plot: lists of points can not be plotted with Nrow and Ncolumn\n");
      // the heights of the spectra in a waterfall plot or a colour map
      YPos = (double*)malloc(a->NSpectra*sizeof(double));
      if(Options.UseResonantEnergy)
      {
        for(i=0;i<a->NSpectra;i++)
        {
          if(isnan(a->Spec[i].ResonantEnergy)) return luaL_error(L, "Graphics.Plot: UseResonantEnergy is set, but plotted spectrum %d has no resonant energy\n",i+1);
        }
        WaterfallPositions(a, Options.WaterfallOffset, 0, true, YPos);
      }
      else if(Options.ColorMap && (Options.WaterfallShift==0))
      {
        // without a shift the spectra of a colour map lie at their number in the input, counted from 1
        WaterfallPositions(a, Options.WaterfallOffset + Options.SpectraStart, Options.SpectraStep, false, YPos);
      }
      else
      {
        WaterfallPositions(a, Options.WaterfallOffset, Options.WaterfallShift, false, YPos);
      }
      if(Options.ColorMap)
      {
        int const NS = (int)a->NSpectra;
        if(NS<2) return luaL_error(L, "Graphics.Plot: a colour map needs at least two spectra, there is %d\n", NS);
        bool const Increasing = YPos[1] > YPos[0];
        for(i=1;i<NS;i++)
        {
          if( Increasing ? !(YPos[i]>YPos[i-1]) : !(YPos[i]<YPos[i-1]) ) return luaL_error(L, "Graphics.Plot: the spectra of a colour map should lie at strictly increasing or decreasing heights, plotted spectrum %d lies at %f and spectrum %d at %f\n", i, YPos[i-1], i+1, YPos[i]);
        }
        // by default the map shows all spectra, the outer ones with half the distance to their neighbour
        double const Lowest     = Increasing ? YPos[0]         : YPos[NS-1];
        double const Highest    = Increasing ? YPos[NS-1]      : YPos[0];
        double const LowestGap  = Increasing ? YPos[1]-YPos[0] : YPos[NS-2]-YPos[NS-1];
        double const HighestGap = Increasing ? YPos[NS-1]-YPos[NS-2] : YPos[0]-YPos[1];
        if(Options.Frame.PlotRange[0].Y==-INFINITY) { Options.Frame.PlotRange[0].Y = Lowest  - LowestGap/2;  AutoYRange=true; }
        if(Options.Frame.PlotRange[1].Y== INFINITY) { Options.Frame.PlotRange[1].Y = Highest + HighestGap/2; AutoYRange=true; }
      }
      if(Options.Frame.PlotRange[0].X==-INFINITY) Options.Frame.PlotRange[0].X = a->Spec[0].Emin;
      if(Options.Frame.PlotRange[1].X== INFINITY) Options.Frame.PlotRange[1].X = a->Spec[0].Emax;
      if( (!Options.ColorMap) && ( (Options.WaterfallShift!=0) || (Options.WaterfallOffset!=0) || Options.UseResonantEnergy )
          && ( (Options.Frame.PlotRange[0].Y==-INFINITY) || (Options.Frame.PlotRange[1].Y==INFINITY) ) )
      {
        // the y range of a waterfall plot covers the base lines and all plotted points inside the x range,
        // with a margin of 5% of that range
        double const YScale = (Options.Negate?-1:1)*Options.Scale;
        double Low=INFINITY, High=-INFINITY;
        for(i=0;i<a->NSpectra;i++)
        {
          SpectraType const * const S = &(a->Spec[i]);
          if(YPos[i]<Low)  Low =YPos[i];
          if(YPos[i]>High) High=YPos[i];
          for(int k=0;k<=(int)S->N;k++)
          {
            double const x = ((PointX!=NULL)&&(PointX[i]!=NULL)) ? PointX[i][k] : (S->N==0) ? S->Emin : ( k * S->Emax + (S->N-k) * S->Emin ) / S->N;
            if((x<Options.Frame.PlotRange[0].X)||(x>Options.Frame.PlotRange[1].X)) continue;
            double const v[2] = {YScale*S->Re[k] + YPos[i], YScale*S->Im[k] + YPos[i]};
            bool   const Plotted[2] = {Options.PlotRe, Options.PlotIm};
            for(int p=0;p<2;p++)
            {
              if(!Plotted[p]) continue;
              if(v[p]<Low)  Low =v[p];
              if(v[p]>High) High=v[p];
            }
          }
        }
        double Margin = 0.05*(High-Low);
        if(!(Margin>0)) Margin = 0.5;
        if(Options.Frame.PlotRange[0].Y==-INFINITY) Options.Frame.PlotRange[0].Y = Low  - Margin;
        if(Options.Frame.PlotRange[1].Y== INFINITY) Options.Frame.PlotRange[1].Y = High + Margin;
        AutoYRange=true;
      }
      if(Options.Frame.PlotRange[0].Y==-INFINITY) Options.Frame.PlotRange[0].Y = -1;
      if(Options.Frame.PlotRange[1].Y== INFINITY) Options.Frame.PlotRange[1].Y =  1;
      if(Options.Frame.FirstXTick==-INFINITY) Options.Frame.FirstXTick = Options.Frame.PlotRange[0].X;
      if(Options.Frame.FirstYTick==-INFINITY)
      {
        if(AutoYRange && (Options.Frame.dYTick>0)) Options.Frame.FirstYTick = ceil(Options.Frame.PlotRange[0].Y/Options.Frame.dYTick) * Options.Frame.dYTick;
        else                                       Options.Frame.FirstYTick = Options.Frame.PlotRange[0].Y;
      }
      
      if(Options.ColorMap)
      {
        double const YScale = (Options.Negate?-1:1)*Options.Scale;
        if(Options.ColorFunction==LUA_NOREF)
        {
          BuiltInColorContext Context;
          Context.RGB[0] = (double)ColorRed(Options.Color)   / 255; // Color is 0xAARRGGBB
          Context.RGB[1] = (double)ColorGreen(Options.Color) / 255;
          Context.RGB[2] = (double)ColorBlue(Options.Color)  / 255;
          Context.Saturation = Options.Saturation;
          Context.UseRe = (Options.PlotIm==false) && (Options.PlotRe==true);
          if(PlotSpectraTypeArrayToColorMap(Graphics, a, YPos, &(Options.Frame), YScale, Options.NPixelX, Options.NPixelY, Options.Interpolate, BuiltInColorFunction, &Context)) return luaL_error(L, "Graphics.Plot: making the colour map failed\n");
        }
        else
        {
          LuaColorFunctionContext Context;
          Context.L = L;
          Context.Ref = Options.ColorFunction;
          Context.Status = LUA_OK;
          Context.Caller = "Graphics.Plot";
          if(PlotSpectraTypeArrayToColorMap(Graphics, a, YPos, &(Options.Frame), YScale, Options.NPixelX, Options.NPixelY, Options.Interpolate, LuaColorFunction, &Context)) return luaL_error(L, "Graphics.Plot: making the colour map failed\n");
          if(Context.Status!=LUA_OK)
          {
            // the error message of the ColorFunction is on the top of the stack; the pixels belong to Graphics
            luaL_unref(L, LUA_REGISTRYINDEX, Options.ColorFunction);
            free(Options.Style.Line);
            free(View.Spec);
            free(YPos);
            FreeSpectraTypeArray(&aa);
            return lua_error(L);
          }
        }
      }
      else if((Options.Nrow==1)&&(Options.Ncolumn==1))
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
            LPORe[0].LineColor = ColorWithTransparency(LPORe[0].LineColor, 0xFF);
            LPORe[0].PointSize=0;
            LPORe[0].PointColor = ColorWithTransparency(LPORe[0].PointColor, 0xFF);
          }
          else if((Options.PlotIm==false) && (Options.PlotRe==true))
          {
            LPORe[0] = Options.Style.Line[0];
            LPOIm[0] = Options.Style.Line[0];
            LPOIm[0].LineWidth=0;
            LPOIm[0].LineColor = ColorWithTransparency(LPOIm[0].LineColor, 0xFF);
            LPOIm[0].PointSize=0;
            LPOIm[0].PointColor = ColorWithTransparency(LPOIm[0].PointColor, 0xFF);
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
              LPORe[i] = Options.Style.Line[(i%(Options.LineStyleCycle))];
              LPOIm[i] = Options.Style.Line[(i%(Options.LineStyleCycle))];
              LPORe[i].LineWidth=0;
              LPORe[i].LineColor = ColorWithTransparency(LPORe[i].LineColor, 0xFF);
              LPORe[i].PointSize=0;
              LPORe[i].PointColor = ColorWithTransparency(LPORe[i].PointColor, 0xFF);
            }
          }
          else if((Options.PlotIm==false) && (Options.PlotRe==true))
          {
            for(i=0;i<a->NSpectra;i++)
            {
              LPORe[i] = Options.Style.Line[(i%(Options.LineStyleCycle))];
              LPOIm[i] = Options.Style.Line[(i%(Options.LineStyleCycle))];
              LPOIm[i].LineWidth=0;
              LPOIm[i].LineColor = ColorWithTransparency(LPOIm[i].LineColor, 0xFF);
              LPOIm[i].PointSize=0;
              LPOIm[i].PointColor = ColorWithTransparency(LPOIm[i].PointColor, 0xFF);
            }
          }
          else if((Options.PlotIm==true) && (Options.PlotRe==true))
          {
            for(i=0;i<a->NSpectra;i++)
            {
              LPORe[i] = Options.Style.Line[(  (2*i)%(Options.LineStyleCycle))];
              LPOIm[i] = Options.Style.Line[((2*i+1)%(Options.LineStyleCycle))];
            }
          }
          else
          {
            return luaL_error(L, "Plot a spectrum but both PlotRe and PlotIm are set to false so there is nothing to do here\n");
          }
        }
        PlotSpectraTypeArrayToGraphics(Graphics, a, LPORe, LPOIm, &(Options.Frame), (Options.Negate?-1:1)*Options.Scale, YPos, PointX);
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
          LPORe.LineColor = ColorWithTransparency(LPORe.LineColor, 0xFF);
          LPORe.PointSize=0;
          LPORe.PointColor = ColorWithTransparency(LPORe.PointColor, 0xFF);
        }
        else if((Options.PlotIm==false) && (Options.PlotRe==true))
        {
          LPORe = Options.Style.Line[0];
          LPOIm = Options.Style.Line[0];
          LPOIm.LineWidth=0;
          LPOIm.LineColor = ColorWithTransparency(LPOIm.LineColor, 0xFF);
          LPOIm.PointSize=0;
          LPOIm.PointColor = ColorWithTransparency(LPOIm.PointColor, 0xFF);
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
        PlotSpectraTypeArrayToGraphicsArray(Graphics, a, Options.Ncolumn, Options.Nrow, &LPORe, &LPOIm, &(Options.Frame),  (Options.Negate?-1:1)*Options.Scale, Options.ScaleOffDiagonal, Options.ShiftOffDiagonal);
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
      if(len>3) NE=(int)luaL_checkinteger(L, 4);
      if(len>4) return luaL_error(L, "For plotting a response function I expect the first argument to be a response function and the 2nd, 3rd and 4th are optional the Lorentzian and Gaussian broadingings and the number of points in the spectrum.\n In this case I got %d argumens\n",len);
      
      // the energies at which the response function is evaluated; the y range and the ticks
      // follow below, after the jump to the plot of a list of spectra
      if(Options.Frame.PlotRange[0].X==-INFINITY) Options.Frame.PlotRange[0].X = -1;
      if(Options.Frame.PlotRange[1].X== INFINITY) Options.Frame.PlotRange[1].X =  1;
      
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
      
      if(SelectSpectra || Options.ColorMap)
      {
        // plot the (selected) elements as a list of spectra in one graph; aa is freed at the end
        FreeResponsefunction(A);
        aa = Spec;
        a = &aa;
        goto labelSelectSpectra;
      }
      if(Options.Frame.PlotRange[0].Y==-INFINITY) Options.Frame.PlotRange[0].Y = -1;
      if(Options.Frame.PlotRange[1].Y== INFINITY) Options.Frame.PlotRange[1].Y =  1;
      if(Options.Frame.FirstXTick==-INFINITY) Options.Frame.FirstXTick = Options.Frame.PlotRange[0].X;
      if(Options.Frame.FirstYTick==-INFINITY) Options.Frame.FirstYTick = Options.Frame.PlotRange[0].Y;
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
          LPORe[0].LineColor = ColorWithTransparency(LPORe[0].LineColor, 0xFF);
          LPORe[0].PointSize=0;
          LPORe[0].PointColor = ColorWithTransparency(LPORe[0].PointColor, 0xFF);
        }
        else if((Options.PlotIm==false) && (Options.PlotRe==true))
        {
          LPORe[0] = Options.Style.Line[0];
          LPOIm[0] = Options.Style.Line[0];
          LPOIm[0].LineWidth=0;
          LPOIm[0].LineColor = ColorWithTransparency(LPOIm[0].LineColor, 0xFF);
          LPOIm[0].PointSize=0;
          LPOIm[0].PointColor = ColorWithTransparency(LPOIm[0].PointColor, 0xFF);
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
        PlotSpectraTypeArrayToGraphics(Graphics, &Spec, LPORe, LPOIm, &(Options.Frame), (Options.Negate?-1:1)*Options.Scale, &(Options.WaterfallOffset), NULL); // one spectrum
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
          LPORe.LineColor = ColorWithTransparency(LPORe.LineColor, 0xFF);
          LPORe.PointSize=0;
          LPORe.PointColor = ColorWithTransparency(LPORe.PointColor, 0xFF);
        }
        else if((Options.PlotIm==false) && (Options.PlotRe==true))
        {
          LPORe = Options.Style.Line[0];
          LPOIm = Options.Style.Line[0];
          LPOIm.LineWidth=0;
          LPOIm.LineColor = ColorWithTransparency(LPOIm.LineColor, 0xFF);
          LPOIm.PointSize=0;
          LPOIm.PointColor = ColorWithTransparency(LPOIm.PointColor, 0xFF);
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
        PlotSpectraTypeArrayToGraphicsArray(Graphics, &Spec, Options.Nrow, Options.Ncolumn, &LPORe, &LPOIm, &(Options.Frame),  (Options.Negate?-1:1)*Options.Scale, Options.ScaleOffDiagonal, Options.ShiftOffDiagonal);
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
  free(View.Spec);
  free(YPos);
  for(i=0;i<NPointX;i++) free(PointX[i]);
  free(PointX);
  luaL_unref(L, LUA_REGISTRYINDEX, Options.ColorFunction);
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
  if (strcmp(str, "Rotate") == 0)
  {
    lua_pushvalue(L, 1);
    lua_pushcclosure(L, LuaRotateGraphicsMethod,1);
    return 1;
  }
  return luaL_error(L, "Error index %s not found in Graphics Object\n",str);
}

// The layers of a graphics have their own origin and size in the absolute coordinates of the
// page. Moving a graphics moves all layers in it by the same distance, and scaling it scales
// them about its origin O by the factors f, so that its layers stay where they are on it.
static void MoveLayers(GraphicsType *G, double dX, double dY)
{
  for(int i=1;i<G->NLayer;i++)
  {
    G->Layer[i]->CanvasOrigin.X += dX;
    G->Layer[i]->CanvasOrigin.Y += dY;
    MoveLayers(G->Layer[i], dX, dY);
  }
}

static void ScaleLayers(GraphicsType *G, PointType O, double fX, double fY)
{
  for(int i=1;i<G->NLayer;i++)
  {
    G->Layer[i]->CanvasOrigin.X = O.X + fX * (G->Layer[i]->CanvasOrigin.X - O.X);
    G->Layer[i]->CanvasOrigin.Y = O.Y + fY * (G->Layer[i]->CanvasOrigin.Y - O.Y);
    G->Layer[i]->CanvasSize.X  *= fX;
    G->Layer[i]->CanvasSize.Y  *= fY;
    ScaleLayers(G->Layer[i], O, fX, fY);
  }
}

// setting SizeX, SizeY, OriginX or OriginY resizes or moves the graphics with all its layers;
// a size of 0 has no scale factor, then only the graphics itself gets the new size
int LuaGraphicsNewIndex(lua_State * L)
{
  GraphicsType *Graphics = (GraphicsType *) luaL_checkudata(L, 1, "Graphics_Type");
  const char *str = lua_tostring(L, 2);
  double a;
  if(strcmp(str, "SizeX")==0)
  {
    a = (double)luaL_checknumber(L, 3);
    if(Graphics->CanvasSize.X!=0) ScaleLayers(Graphics, Graphics->CanvasOrigin, a/Graphics->CanvasSize.X, 1);
    Graphics->CanvasSize.X = a;
    return 0;
  }
  if(strcmp(str, "SizeY")==0)
  {
    a = (double)luaL_checknumber(L, 3);
    if(Graphics->CanvasSize.Y!=0) ScaleLayers(Graphics, Graphics->CanvasOrigin, 1, a/Graphics->CanvasSize.Y);
    Graphics->CanvasSize.Y = a;
    return 0;
  }
  if(strcmp(str, "OriginX")==0)
  {
    a = (double)luaL_checknumber(L, 3);
    MoveLayers(Graphics, a-Graphics->CanvasOrigin.X, 0);
    Graphics->CanvasOrigin.X = a;
    return 0;
  }
  if(strcmp(str, "OriginY")==0)
  {
    a = (double)luaL_checknumber(L, 3);
    MoveLayers(Graphics, 0, a-Graphics->CanvasOrigin.Y);
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
  uint32_t fillcolor=0x00FFFFFF; // opaque white
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
      // the colour 0xAARRGGBB with the opacity alpha (255 opaque) on top of its transparency AA
      fillcolor = ColorCombine(Options.color, 0, Options.alpha);
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
  // the colour 0xAARRGGBB, AA the transparency, and with a fifth argument an opacity alpha (255 opaque) on top
  BitMap->Pixels[(i-1) + (j-1)*BitMap->NpointsX] = rgb;
  if(len==5)
  {
    uint8_t alpha = (uint8_t) luaL_checkinteger(L, 5);
    BitMap->Pixels[(i-1) + (j-1)*BitMap->NpointsX] = ColorCombine(rgb, 0, alpha);
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
                                   double GammaL, double GammaG, double Saturation, uint32_t LineColor, uint32_t SecondLineColor, double ZeroDistance)
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
  Options.LineColor = 0x00000080;       // 0xAARRGGBB, AA the transparency: opaque dark blue
  Options.LineTransparancy = 0;
  Options.SecondLineColor = 0x00FF0000; // opaque red
  Options.SecondLineTransparancy = 0;
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
  Options.Frame.LineColor = 0x00000000;
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
  Options.Frame.XTicks=true;
  Options.Frame.YTicks=true;
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
    len--;
  }
  // Transparancy (0 opaque, 255 fully transparent) on top of the transparency AA of the color
  Options.LineColor       = ColorCombine(Options.LineColor,       Options.LineTransparancy,       255);
  Options.SecondLineColor = ColorCombine(Options.SecondLineColor, Options.SecondLineTransparancy, 255);
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
  uint32_t const LineColor       = Options.LineColor;       // 0xAARRGGBB, AA the transparency
  uint32_t const SecondLineColor = Options.SecondLineColor;
  // allocate the bitmap
  BitMapType BandPlot;
  InitBitMap(&BandPlot);
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
  {"ToTable", LuaGraphicsToTable},
  {"Add", LuaAddToGraphicsFunction},
  {"Rotate", LuaRotateGraphicsFunction},
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
