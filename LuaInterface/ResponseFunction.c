/*  Response functions are stored as matrix representations and reside fully on the lua stack as tables
 *
 *  There are different types of response functions
 *    Tri
 *    And
 *    ListOfPoles
 *    NaturalImpurityOrbital
 *  each type can either be single or matrix represented
 *
 *  In order for the meta methods to function one needs to set the meta table "ResponseFunction_Type"
 *  There is no type-check, (i.e. setting the meta table is optional) for the functions acting on these tables.
 *
 *  
 *
 */


/*
 *
 */
static int LuaChangeTypeResponseFunction(lua_State *L)
{
  printf("LuaChangeTypeResponseFunction\n");
  return 0;
}

/*
 *
 */
static int LuaReducePolesResponseFunction(lua_State *L)
{
  printf("LuaReducePoles\n");
  return 0;
}

/*
 *
 */
static int LuaPlaceOnMeshResponseFunction(lua_State *L)
{
  printf("LuaPlaceOnMeshResponseFunction");
  return 0;
}

/*
 *
 */
static int LuaResponseFunctionUnm(lua_State * L)
{
  printf("LuaResponseFunctionUnm\n");
  return 0;
}

/*
 *
 */
static int LuaResponseFunctionAdd(lua_State * L)
{
  printf("LuaResponseFunctionAdd\n");
  
  lua_pushnumber(L, 1.23456);
  return 1;
}

/*
 *
 */
static int LuaResponseFunctionSub(lua_State * L)
{
  printf("LuaResponseFunctionSub\n");
  return 0;
}

/*
 *
 */
static int LuaResponseFunctionMul(lua_State * L)
{
  printf("LuaResponseFunctionMul\n");
  return 0;
}

/*
 *
 */
static int LuaResponseFunctionDiv(lua_State * L)
{
  printf("LuaResponseFunctionDiv\n");
  return 0;
}

/*
 *
 */
static int LuaResponseFunctionPow(lua_State * L)
{
  printf("LuaResponseFunctionPow\n");
  return 0;
}

/*
 *
 */
static int LuaResponseFunctionConcat(lua_State * L)
{
  printf("LuaResponseFunctionConcat\n");
  return 0;
}

/*
 *
 */
static int LuaResponseFunctionEq(lua_State * L)
{
  printf("LuaResponseFunctionEq\n");
  return 0;
}

/*
 *
 */
static int LuaResponseFunctionLt(lua_State * L)
{
  printf("LuaResponseFunctionLt\n");
  return 0;
}

/*
 *
 */
static int LuaResponseFunctionLe(lua_State * L)
{
  printf("LuaResponseFunctionLe\n");
  return 0;
}

/*
 *
 */
static int LuaResponseFunctionIndex(lua_State * L)
{
  printf("LuaResponseFunctionIndex\n");
  return 0;
}

/*
 *
 */
static int LuaResponseFunctionNewIndex(lua_State * L)
{
  printf("LuaResponseFunctionNewIndex\n");
  return 0;
}

/*
 *
 */
static int LuaResponseFunctionMode(lua_State * L)
{
  printf("LuaResponseFunctionMode\n");
  return 0;
}

/*
 *
 */
static int LuaResponseFunctionCall(lua_State * L)
{
  printf("LuaResponseFunctionCall\n");
  return 0;
}

static const struct luaL_Reg ResponseFunction_func[] = {
  {"ChangeType", LuaChangeTypeResponseFunction},
  {"ReducePoles", LuaReducePolesResponseFunction},
  {"PlaceOnMesh", LuaPlaceOnMeshResponseFunction},
  {NULL, NULL}
};

static const struct luaL_Reg ResponseFunction_meta[] = {
  {"__unm", LuaResponseFunctionUnm},
  {"__add", LuaResponseFunctionAdd},
  {"__sub", LuaResponseFunctionSub},
  {"__mul", LuaResponseFunctionMul},
  {"__div", LuaResponseFunctionDiv},
  {"__pow", LuaResponseFunctionPow},
  {"__concat", LuaResponseFunctionConcat},
  {"__eq", LuaResponseFunctionEq},
  {"__lt", LuaResponseFunctionLt},
  {"__le", LuaResponseFunctionLe},
  {"__index", LuaResponseFunctionIndex},
  {"__newindex", LuaResponseFunctionNewIndex},
  {"__mode", LuaResponseFunctionMode},
  {"__call", LuaResponseFunctionCall},
  {NULL, NULL}
};