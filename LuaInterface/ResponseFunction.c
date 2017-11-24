/*  Response functions are stored either as matrix representations and reside fully on the lua stack as tables
 *                                or as BlockMatrixType
 *
 *  There are different types of response functions
 *    Tri
 *    And
 *    ListOfPoles
 *    NaturalImpurityOrbital
 *  each type can either be single or matrix represented
 *
 *  In order for the meta methods to function for table representations one needs to set the meta table "ResponseFunctionMeta" in lua
 *  There is no type-check, (i.e. setting the meta table is optional) for the functions acting on these tables.
 *
 *
 *
 */


/*   change between l,t,a, and n type
 *   change between representation as table or usertype (table capital, user data small)
 *   parameters on lua stack (responsefunction, typedef)
 *   typedef being a string of the form:
 *   Tri, And, ListOfPoles, NaturalImpurityOrbital, Userdata, Table
 */
static int LuaChangeTypeResponseFunction(lua_State *L)
{
  BlockMatrixType *Min, *Mout;
  char const * convert;
  char type;
  double T=0;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len!=2) return luaL_error(L, "Error in ChangeType expected two arguments, not %d\n",len);
  Min = luaL_checkBlockGreensfunction(L,1);
  convert = luaL_checkstring(L, 2);
  if(strcmp(convert, "Tri")==0)
  {
    type='T'|(32&Min->type);  // t = T & 32, i.e. keep capital or not from input, but change letter to t/T
  }
  else if(strcmp(convert, "And")==0)
  {
    type='A'|(32&Min->type);
  }
  else if(strcmp(convert, "ListOfPoles")==0)
  {
    type='L'|(32&Min->type);
  }
  else if(strcmp(convert, "NaturalImpurityOrbital")==0)
  {
    type='N'|(32&Min->type);
  }
  else if(strcmp(convert, "Userdata")==0)
  {
    type=(32)&Min->type; // to lower case
  }
  else if(strcmp(convert, "Table")==0)
  {
    type=(255-32)&Min->type; // to upper case
  }
  else
  {
    return luaL_error(L, "Error in ChangeType, second argument needs to be a string containgin either: \"And\", \"Tri\", \"ListOfPoles\", \"NaturalImpurityOrbital\", \"Userdata\", \"Table\" and not %s\n",convert);
  }
  Mout = ConvertBlockMatrix(Min, type, T);
  LuaPush_BlockGreensFunction(L, Mout);
  FreeBlockMatrixFromLua(Mout);
  FreeBlockMatrixFromLua(Min);
  return 1;
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
  char temp;
  printf("LuaResponseFunctionAdd\n");

  for(temp='A';temp<='Z';temp++)
  {
    printf("original %c, lowercased %c\n",temp,(temp|32));
  }

  for(temp='a';temp<='z';temp++)
  {
    printf("original %c, lowercased %c\n",temp,(temp|32));
  }



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
