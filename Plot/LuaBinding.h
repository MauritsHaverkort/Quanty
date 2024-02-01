int LuaNewGraphics(lua_State * L);
int LuaAddToGraphics(lua_State * L);
int LuaPlotGraphics(lua_State * L);
int LuaGraphicsToPS(lua_State * L);
int LuaGraphicsToSVG(lua_State * L);
int LuaGraphicsAdd(lua_State * L);
int LuaGraphicsLen(lua_State * L);
int LuaGraphicsToString(lua_State * L);
int LuaGraphicsIndex(lua_State * L);
int LuaGraphicsDestroy(lua_State * L);
int LuaNewBitMap(lua_State * L);
int LuaBitMapFromMatrix(lua_State * L);
int LuaBitMapOverlay(lua_State * L);
int LuaBitMapAdd(lua_State * L);
int LuaBitMapLen(lua_State * L);
int LuaBitMapToString(lua_State * L);
int LuaBitMapDestroy(lua_State * L);

static int LuaPlotDOS(lua_State *L);
static int LuaPlotBands(lua_State *L);
static int LuaPlotBandsAndDOS(lua_State *L);
static int LuaPlotFermi(lua_State *L);


