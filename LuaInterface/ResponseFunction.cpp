/*  Response functions are stored either as matrix representations and reside fully on the lua stack as tables
 *                                or as ResponsefunctionType
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
  ResponsefunctionType *Min, *Mout;
  char const * convert;
  char type;
  double T=sqrt(DBL_MIN);
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len!=2) return luaL_error(L, "Error in ChangeType expected 2 arguments, not %d\n",len);
  Min = luaL_checkResponsefunction(L,1);
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
  else if((strcmp(convert, "NaturalImpurityOrbital")==0)||(strcmp(convert, "Nat")==0))
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
  if(Min->Block==true)
  {
    Mout = ConvertBlockMatrix(Min, type, T);
  }
  else
  {
    Mout = ConvertResponsefunction(Min, type, T);
  }
  LuaPush_ResponseFunction(L, Mout);
  FreeResponsefunction(Mout);
  FreeResponsefunction(Min);
  return 1;
}

static int LuaResponseFunctionToUserdata(lua_State *L)
{
  ResponsefunctionType *M;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len!=1) return luaL_error(L, "Error in ToUserdata expected a single argument, not %d\n",len);
  M = luaL_checkResponsefunction(L,1);
  if((M->type|32)==M->type) //lowercase, i.e. only pointer stored
  {
    LuaPush_ResponseFunction(L, ResponsefunctionCopy(M));
  }
  else // upercase read a table from the stack
  {
    M->type+=32;
    LuaPush_ResponseFunction(L, M);
  }
  return 1;
}


static int LuaResponseFunctionToTable(lua_State *L)
{
  ResponsefunctionType *M;
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len!=1) return luaL_error(L, "Error in ToTable expected a single argument, not %d\n",len);
  M = luaL_checkResponsefunction(L,1);
  switch(M->type)
  {
    case 'L':
    case 'l':
    {
      if(M->Block==true)
      {
        LuaPushTable_BlockListOfPoles(L, (BlockListOfPoles*)M->data);
      }
      else
      {
        LuaPushTable_ListOfPoles(L, (ListOfPoles*)(M->data));
      }
      break;
    }
    case 'A':
    case 'a':
    {
      if(M->Block==true)
      {
        LuaPushTable_BlockAndersonMatrixType(L, M->data);
      }
      else
      {
        LuaPushTable_AndersonMatrixType(L, M->data);
      }
      break;
    }
    case 'T':
    case 't':
    {
      if(M->Block==true)
      {
        LuaPushTable_BlockTriDiagonalMatrixType(L, M->data);
      }
      else
      {
        LuaPushTable_TriDiagonalMatrixType(L, M->data);
      }
      break;
    }
    case 'N':
    case 'n':
    {
      if(M->Block==true)
      {
        LuaPushTable_BlockNaturalImpurityOrbital(L, M->data);
      }
      else
      {
        LuaPushTable_NaturalImpurityOrbital(L, M->data);
      }
      break;
    }
    default:
    {
      printf("Error in ResponseFunctionToTable: unrecognized type %c",M->type);
      return 0;
    }
  }
  FreeResponsefunction(M);
  return 1;
}

/*
 *
 */
static int LuaReducePolesResponseFunction(lua_State *L)
{
  ResponsefunctionType *M;
  int Nred;

  uint32_t len = (uint32_t)lua_gettop(L);
  if(len!=2) return luaL_error(L, "Error in ReducePoles: %d argument(s), 2 expected (ListOfPoles, desired number of poles).\n", len);

  M = luaL_checkResponsefunction(L,1);
  if(!lua_isnumber(L, 2)) return luaL_error(L, "Error in ReducePoles: 2nd argument should be the desired number of poles.\n");
  Nred = (int)lua_tointeger(L, 2);
  if((M->type|32)!='l') return luaL_error(L, "Error in ReducePoles: The input Green's function is of %c type, it should be a ListOfPoles", M->type);

  ReduceBlockListOfPoles((BlockListOfPoles*)(M->data), Nred, ((BlockListOfPoles*)(M->data))->mu);
  LuaPush_ResponseFunction(L, M);

  return 1;
}

/*
 *
 */
static int LuaPlaceOnMeshResponseFunction(lua_State *L)
{
  ResponsefunctionType *M;
  uint32_t NE;
  double *mesh;

  uint32_t len = (uint32_t)lua_gettop(L);
  if(len!=2) return luaL_error(L, "Error in PlaceOnMesh: %d argument(s), 2 expected (ListOfPoles, desired list of grid points).\n", len);
  if(lua_type(L,2)!=LUA_TTABLE) return luaL_error(L, "Error in PlaceOnMesh: the 2nd argument should be a list of energy grid points.\n");

  M = luaL_checkResponsefunction(L,1);
  if((M->type|32)!='l')
  {
    return luaL_error(L, "Error in PlaceOnMesh: The first argument should be of type 'L'/'l' (ListOfPoles).\n");
  }
  mesh = luaL_checkarray_double(L, 2, &NE);
  BlockListOfPolesPlaceOnMesh((BlockListOfPoles*)M->data, mesh, NE);
  LuaPush_ResponseFunction(L, M);

  return 1;
}

/*
 *
 */
static int LuaResponseFunctionEigensystem(lua_State *L)
{
  ResponsefunctionType *M;
  double *val;
  CompactMatrixType fun;
  M = luaL_checkResponsefunction(L,-1);
  
  switch(M->type)
  {
    case 'L':
    case 'l':
    {
      if(M->Block==true)
      {
        val = (double*)malloc(((BlockListOfPoles*)M->data)->N*((BlockListOfPoles*)M->data)->BlockSize*sizeof(double));
        BlockListOfPolesEigenSystem((BlockListOfPoles*)M->data, val, &fun);
      }
      else
      {
        val = (double*)malloc(((ListOfPoles*)M->data)->N*sizeof(double));
        ListOfPolesToEigensystemCompact((ListOfPoles*)(M->data), val, &fun);
      }
      break;
    }
    case 'A':
    case 'a':
    {
      if(M->Block==true)
      {
        val = (double*)malloc(((BlockAndersonMatrixType*)M->data)->N*((BlockAndersonMatrixType*)M->data)->BlockSize*sizeof(double));
        BlockAndersonMatrixEigenSystem((BlockAndersonMatrixType*)M->data, val, &fun);
      }
      else
      {
        val = (double*)malloc(((AndersonMatrixType*)M->data)->N*sizeof(double));
        AndersonMatrixToEigenSystemCompact((AndersonMatrixType*)M->data, val, &fun);
      }
      break;
    }
    case 'T':
    case 't':
    {
      if(M->Block==true)
      {
        val = (double*)malloc(((BlockTriDiagonalMatrixType*)M->data)->N*((BlockTriDiagonalMatrixType*)M->data)->BlockSize*sizeof(double));
        BlockTriDiagonalMatrixEigenSystem((BlockTriDiagonalMatrixType*)M->data, val, &fun);
      }
      else
      {
        val = (double*)malloc(((TriDiagonalMatrixType*)M->data)->N*sizeof(double));
        TriDiagonalMatrixToEigensystemCompact((TriDiagonalMatrixType*)M->data, val, &fun);
      }
      break;
    }
    case 'N':
    case 'n':
    {
      if(M->Block==true)
      {
        val = (double*)malloc((((BlockNaturalImpurityOrbitalMatrixType*)M->data)->val.N+ ((BlockNaturalImpurityOrbitalMatrixType*)M->data)->con.N )*((BlockNaturalImpurityOrbitalMatrixType*)M->data)->BlockSize*sizeof(double));
        BlockNaturalImpurityOrbitalMatrixEigenSystem((BlockNaturalImpurityOrbitalMatrixType*)M->data, val, &fun);
      }
      else
      {
        val = (double*)malloc((((BADoubleTriDiagonalMatrixType*)M->data)->val.N+((BADoubleTriDiagonalMatrixType*)M->data)->con.N+1)*sizeof(double));
        BADoubleTriDiagonalMatrixTypeToEigenSystemCompact((BADoubleTriDiagonalMatrixType*)M->data, val, &fun);
      }
      break;
    }
    default:
    {
      printf("Error in ResponseFunctionToTable: unrecognized type %c",M->type);
      return 0;
    }
  }
  LuaPushArray_double(L, val, fun.Ni);
  LuaPushTable_CompactMatrixType(L, &fun);
  FreeResponsefunction(M);
  FreeCompactMatrix(&fun);
  free(val);
  return 2;
}


/*
 *
 */
static int LuaResponseFunctionUnm(lua_State * L)
{
  ResponsefunctionType *res, *M;

  M = ConvertBlockMatrix(luaL_checkResponsefunction(L,1), 'l', 0);
  res = (ResponsefunctionType *) lua_newuserdata(L, sizeof(ResponsefunctionType));
  luaL_setmetatable(L, "ResponseFunction_Type");
  res->NSuggestedGridPoints = 0;
  res->type = 'L';
  res->data = (BlockListOfPoles*)malloc(sizeof(BlockListOfPoles));
  CopyBlockListOfPoles(*(BlockListOfPoles*)(M->data),res->data);
  ScaleBlockListOfPoles(res->data, -1);
  LuaPush_ResponseFunction(L, res);

  return 1;
}

/*
 *
 */
static int LuaResponseFunctionAdd(lua_State * L)
{
  ResponsefunctionType *Ma, *Mb, *res;

  uint32_t len = (uint32_t)lua_gettop(L);
  if(len!=2) return luaL_error(L, "Error in ResponseFunctionAdd: %d argument(s), 2 expected.\n", len);

  Ma = ConvertBlockMatrix(luaL_checkResponsefunction(L,1), 'l', 0);
  Mb = ConvertBlockMatrix(luaL_checkResponsefunction(L,2), 'l', 0);
  AddBlockListOfPoles((BlockListOfPoles*)(Ma->data), *(BlockListOfPoles*)(Mb->data));

  res = (ResponsefunctionType *) lua_newuserdata(L, sizeof(ResponsefunctionType));
  luaL_setmetatable(L, "ResponseFunction_Type");
  res->NSuggestedGridPoints = 0;
  res->type = 'L';
  res->data = (BlockListOfPoles*)malloc(sizeof(BlockListOfPoles));
  CopyBlockListOfPoles(*(BlockListOfPoles*)(Ma->data),res->data);
  LuaPush_ResponseFunction(L, res);

  FreeBlockMatrixFromLua(Ma);
  FreeBlockMatrixFromLua(Mb);

  return 1;
}

/*
 *
 */
static int LuaResponseFunctionSub(lua_State * L)
{
  ResponsefunctionType *Ma, *Mb, *res;

  uint32_t len = (uint32_t)lua_gettop(L);
  if(len!=2) return luaL_error(L, "Error in ResponseFunctionAdd: %d argument(s), 2 expected.\n", len);

  Ma = ConvertBlockMatrix(luaL_checkResponsefunction(L,1), 'l', 0);
  Mb = ConvertBlockMatrix(luaL_checkResponsefunction(L,2), 'l', 0);
  SubtractBlockListOfPoles((BlockListOfPoles*)(Ma->data), *(BlockListOfPoles*)(Mb->data));

  res = (ResponsefunctionType *) lua_newuserdata(L, sizeof(ResponsefunctionType));
  luaL_setmetatable(L, "ResponseFunction_Type");
  res->NSuggestedGridPoints = 0;
  res->type = 'L';
  res->data = (BlockListOfPoles*)malloc(sizeof(BlockListOfPoles));
  CopyBlockListOfPoles(*(BlockListOfPoles*)(Ma->data),res->data);
  LuaPush_ResponseFunction(L, res);

  FreeBlockMatrixFromLua(Ma);
  FreeBlockMatrixFromLua(Mb);

  return 1;
}

/*
 *
 */
static int LuaResponseFunctionMul(lua_State * L)
{
  //TODO: all
  printf("LuaResponseFunctionMul\n");
  return 0;
}

/*
 *
 */
static int LuaResponseFunctionDiv(lua_State * L)
{
  //TODO: all
  printf("LuaResponseFunctionDiv\n");
  return 0;
}

/*
 *
 */
static int LuaResponseFunctionPow(lua_State * L)
{
  //TODO: all
  printf("LuaResponseFunctionPow\n");
  return 0;
}

/*
 *
 */
static int LuaResponseFunctionConcat(lua_State * L)
{
  //TODO: all
  printf("LuaResponseFunctionConcat\n");
  return 0;
}

/*
 *
 */
static int LuaResponseFunctionEq(lua_State * L)
{
  //TODO: all
  printf("LuaResponseFunctionEq\n");
  return 0;
}

/*
 *
 */
static int LuaResponseFunctionLt(lua_State * L)
{
  //TODO: all
  printf("LuaResponseFunctionLt\n");
  return 0;
}

/*
 *
 */
static int LuaResponseFunctionLe(lua_State * L)
{
  //TODO: all
  printf("LuaResponseFunctionLe\n");
  return 0;
}

/*
 *
 */
static int LuaResponseFunctionIndex(lua_State * L)
{
  //TODO: all
  printf("LuaResponseFunctionIndex\n");
  luaMWH_StackDump(L);
  return 0;
}

/*
 *
 */
static int LuaResponseFunctionNewIndex(lua_State * L)
{
  //TODO: all
  printf("LuaResponseFunctionNewIndex\n");
  return 0;
}

/*
 *
 */
static int LuaResponseFunctionMode(lua_State * L)
{
  //TODO: all
  printf("LuaResponseFunctionMode\n");
  return 0;
}

/*
 *
 */
static int LuaResponseFunctionCall(lua_State * L)
{
  ResponsefunctionType *res;
  double omega, gamma;
  if (lua_gettop(L)!=3)
  {
    return luaL_error(L, "Call to a response-function should contain two elements, (omega, gamma) with omega the energy at which the response function is evaluated and gamma an imaginary ofset from the real axis (FWMH lorenzian broadening)\n");
  }
  res = luaL_checkResponsefunction(L,1);
  omega = (double)luaL_checknumber(L, 2);
  gamma = (double)luaL_checknumber(L, 3);
  if(res->Block==true)
  {
    CompactMatrixType out;
    EvaluateBlockResponsefunctionType(res, &out, omega, gamma);
    LuaPushTable_CompactMatrixType(L, &out);
    FreeCompactMatrix(&out);  // added by Martin - otherwise heap is lost.
  }
  else
  {
    Complex16 *out;
    out = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
    *out = EvaluateResponsefunctionType(res, omega, gamma);
    luaL_setmetatable(L, "Complex_Type");
  }
  FreeResponsefunctionIncludingOriginalPointer(res);
  return 1;
}

/*
 *
 */
static int LuaResponseFunctionToString(lua_State * L)
{
  lua_pushfstring(L, "ResponseFunction in userdata format use ToTable() in order to get a table form");
  return 1;
}


/*
 *
 */
static int LuaResponseFunctionToTightBinding(lua_State * L)
{
  TightBindingDef *TB;
  ResponsefunctionType *M;
  TB = (TightBindingDef*) lua_newuserdata(L, sizeof(TightBindingDef));
  luaL_setmetatable(L, "TightBinding_Type");
  M = luaL_checkResponsefunction(L,1);
  switch(M->type)
  {
    case 'L':
    case 'l':
    {
      if(M->Block==true)
      {
        BlockMatrixToTightBinding((BlockListOfPoles *)(M->data),TB);
      }
      else
      {
        //TODO: listofpoles to TB
      }
      break;
    }
    case 'A':
    case 'a':
    {
      if(M->Block==true)
      {
        BlockMatrixToTightBinding((BlockAndersonMatrixType *)(M->data),TB);
      }
      else
      {
        //TODO: andersonmatrixToTightbinding
      }
      break;
    }
    case 'T':
    case 't':
    {
      if(M->Block==true)
      {
        BlockMatrixToTightBinding((BlockTriDiagonalMatrixType *)(M->data),TB);
      }
      else
      {
        //TODO: tridiagonal matrix to TB
      }
      break;
    }
    case 'N':
    case 'n':
    {
      if(M->Block==true)
      {
        BlockMatrixToTightBinding((BlockNaturalImpurityOrbitalMatrixType *)(M->data),TB);
      }
      else
      {
        //TODO: NatrualImpuritiOrbitalMatrixType to TB
      }
      break;
    }
    default:
    {
      return luaL_error(L, "Error in Responsefunction.ToTightbinding() unrecognized type for conversion (%s)\n",M->type);
      return 0;
    }
  }
  FreeResponsefunction(M);
  return 1;
}

/*
 *
 */
static int LuaResponseFunctionDestroy(lua_State * L)
{
  ResponsefunctionType *a = (ResponsefunctionType *) luaL_checkudata(L, 1, "ResponseFunction_Type");
  // convert type to captial letter, i.e. claim controle over the allocation
  a->type^=32;
  FreeResponsefunction(a);  
  return 0;
}

static const struct luaL_Reg ResponseFunction_func[] = {
  {"ChangeType", LuaChangeTypeResponseFunction},
  {"ReducePoles", LuaReducePolesResponseFunction},
  {"PlaceOnMesh", LuaPlaceOnMeshResponseFunction},
  {"ToTable", LuaResponseFunctionToTable},
  {"ToUserdata", LuaResponseFunctionToUserdata},
  {"ToTightBinding", LuaResponseFunctionToTightBinding},
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
  {"__mode", LuaResponseFunctionMode},
  {"__call", LuaResponseFunctionCall},
  {"__tostring", LuaResponseFunctionToString},
  {"__gc", LuaResponseFunctionDestroy},
  {NULL, NULL}
};

static const struct luaL_Reg ResponseFunctionTable_meta[] = {
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
  {"__mode", LuaResponseFunctionMode},
  {"__call", LuaResponseFunctionCall},
  {NULL, NULL}
};


