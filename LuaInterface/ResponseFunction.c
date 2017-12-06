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

static int LuaResponseFunctionToTable(lua_State *L)
{
  BlockMatrixType *M;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len!=1) return luaL_error(L, "Error in ToTable expected a single argument, not %d\n",len);

  M = luaL_checkBlockGreensfunction(L,1);
  switch(M->type)
  {
    case 'L':
      LuaPushTable_BlockListOfPoles(L, M->data);
      return 1;
    case 'A':
      LuaPushTable_BlockAndersonMatrixType(L, M->data);
      return 1;
    case 'T':
      LuaPushTable_BlockTriDiagonalMatrixType(L, M->data);
      return 1;
    case 'N':
      LuaPushTable_BlockNaturalImpurityOrbital(L, M->data);
      return 1;
    default:
      printf("Error in ResponseFunctionToTable: unrecognized type %c",M->type);
      return 0;
  }

  return 0;
}

/*
 *
 */
static int LuaReducePolesResponseFunction(lua_State *L)
{
  BlockMatrixType *M, *res;
  int Nred;

  uint32_t len = (uint32_t)lua_gettop(L);
  if(len!=2) return luaL_error(L, "Error in ReducePoles: %d argument(s), 2 expected (ListOfPoles, desired number of poles).\n", len);

  M = luaL_checkBlockGreensfunction(L,1);
  if(!lua_isnumber(L, 2)) return luaL_error(L, "Error in ReducePoles: 2nd argument should be the desired number of poles.\n");
  Nred = (int)lua_tointeger(L, 2);
  if((M->type|32)!='l') return luaL_error(L, "Error in ReducePoles: The input Green's function is of %c type, it should be a ListOfPoles", M->type);

  res = (BlockMatrixType *) lua_newuserdata(L, sizeof(BlockMatrixType));
  luaL_setmetatable(L, "ResponseFunction_Type");
  res->NSuggestedGridPoints = 0;
  res->type = 'L';
  res->data = malloc(sizeof(BlockListOfPoles));
  CopyBlockListOfPoles(*(BlockListOfPoles*)(M->data),res->data);

  ReduceBlockListOfPoles((BlockListOfPoles*)(res->data), Nred, ((BlockListOfPoles*)(res->data))->mu);
  LuaPush_BlockGreensFunction(L, res);

  return 1;
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
  BlockMatrixType *Ma, *Mb, *res;

  uint32_t len = (uint32_t)lua_gettop(L);
  if(len!=2) return luaL_error(L, "Error in ResponseFunctionAdd: %d argument(s), 2 ListOfPoles expected.\n", len);

  Ma = ConvertBlockMatrix(luaL_checkBlockGreensfunction(L,1), 'l', 0);
  Mb = ConvertBlockMatrix(luaL_checkBlockGreensfunction(L,2), 'l', 0);
  AddBlockListOfPoles((BlockListOfPoles*)(Ma->data), *(BlockListOfPoles*)(Mb->data));

  res = (BlockMatrixType *) lua_newuserdata(L, sizeof(BlockMatrixType));
  luaL_setmetatable(L, "ResponseFunction_Type");
  res->NSuggestedGridPoints = 0;
  res->type = 'L';
  res->data = malloc(sizeof(BlockListOfPoles));
  CopyBlockListOfPoles(*(BlockListOfPoles*)(Ma->data),res->data);
  LuaPush_BlockGreensFunction(L, res);

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
  BlockMatrixType *res;
  double omega, gamma;
  CompactMatrixType out;
  if (lua_gettop(L)!=3)
  {
    return luaL_error(L, "Call to a response-function should contain two elements, (omega, gamma) with omega the energy at which the response function is evaluated and gamma an imaginary ofset from the real axis (FWMH lorenzian broadening)\n");
  }
  res = luaL_checkBlockGreensfunction(L,1);
  omega = (double)luaL_checknumber(L, 2);
  gamma = (double)luaL_checknumber(L, 3);
  EvaluateBlockMatrixType(res, &out, omega, gamma);
  LuaPushTable_CompactMatrixType(L, &out);
  FreeBlockMatrixFromLua(res);
  return 1;
}

static const struct luaL_Reg ResponseFunction_func[] = {
  {"ChangeType", LuaChangeTypeResponseFunction},
  {"ReducePoles", LuaReducePolesResponseFunction},
  {"PlaceOnMesh", LuaPlaceOnMeshResponseFunction},
  {"ToTable", LuaResponseFunctionToTable},
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
