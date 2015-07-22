static int LuaWavefunctionAdd(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  if(lua_type(L,len)==LUA_TTABLE)
  {
    lua_newtable(L);
    lua_pushnil(L);
    while(lua_next(L, len) != 0)
    {
      LuaWavefunctionAdd(L);
      lua_pushvalue(L, -3);
      lua_insert(L, -2);
      lua_settable(L, -5);
      lua_pop(L, 1);
    }
    return 1;
  }
  if(lua_type(L,1)==LUA_TTABLE)
  {
    if(len!=2) return luaL_error(L, "Error in __add found 3 arguments, wheras a + b should only have two arguments\n");
    lua_insert(L,1);
    return LuaWavefunctionAdd(L);
  }
  WaveFunctionType *l = (WaveFunctionType *) luaL_checkudata(L, 1, "Wavefunction_Type");
  WaveFunctionType *r = (WaveFunctionType *) luaL_checkudata(L, len, "Wavefunction_Type");
  WaveFunctionType *res = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
  luaL_setmetatable(L, "Wavefunction_Type");
  WaveFunctionCopy(*l, res);
  WaveFunctionAddUnionBasis(res, *r);
  return 1;
}


static int LuaWavefunctionSubMin(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  if(lua_type(L,len)==LUA_TTABLE)
  {
    lua_newtable(L);
    lua_pushnil(L);
    while(lua_next(L, len) != 0)
    {
      LuaWavefunctionSubMin(L);
      lua_pushvalue(L, -3);
      lua_insert(L, -2);
      lua_settable(L, -5);
      lua_pop(L, 1);
    }
    return 1;
  }
  WaveFunctionType *l = (WaveFunctionType *) luaL_checkudata(L, 1, "Wavefunction_Type");
  WaveFunctionType *r = (WaveFunctionType *) luaL_checkudata(L, len, "Wavefunction_Type");
  WaveFunctionType *res = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
  luaL_setmetatable(L, "Wavefunction_Type");
  WaveFunctionCopy(*r, res);
  WaveFunctionAddMultipliedUnionBasis(res, *l, -1.0);
  return 1;
}

static int LuaWavefunctionSub(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  if(lua_type(L,len)==LUA_TTABLE)
  {
    lua_newtable(L);
    lua_pushnil(L);
    while(lua_next(L, len) != 0)
    {
      LuaWavefunctionSub(L);
      lua_pushvalue(L, -3);
      lua_insert(L, -2);
      lua_settable(L, -5);
      lua_pop(L, 1);
    }
    return 1;
  }
  if(lua_type(L,1)==LUA_TTABLE)
  {
    if(len!=2) return luaL_error(L, "Error in __sub found 3 arguments, wheras a - b should only have two arguments\n");
    lua_insert(L,1);
    return LuaWavefunctionSubMin(L);
  }
  WaveFunctionType *l = (WaveFunctionType *) luaL_checkudata(L, 1, "Wavefunction_Type");
  WaveFunctionType *r = (WaveFunctionType *) luaL_checkudata(L, len, "Wavefunction_Type");
  WaveFunctionType *res = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
  luaL_setmetatable(L, "Wavefunction_Type");
  WaveFunctionCopy(*l, res);
  WaveFunctionAddMultipliedUnionBasis(res, *r, -1.0);
  return 1;
}

/*  multiplies two objects, whereby one is a wavefunction
 *  a (complex) number times wave-function scales the wave-function - multiplication commutes
 *  a wavefunction times wave-function is implemented as: psiA * psiB == Conjugate(psiA) . psiB - cummutation leads to complex conjugation
 *  a wavefunction times operator is implemented as psi * O == ConjugateTranspose(O) * psi
 *  an operator times wavefunction is implemented as psi * O
 *  a table times wavefunction is implemented as: {a,b,c} * psi == {a * psi, b * psi, c * psi}
 *  a wavefunction times table is implemented as: psi * {a,b,c} == {psi * a, psi * b, psi * c}
 */
static int LuaWavefunctionMul(lua_State * L)
{
  
  uint32_t len = (uint32_t)lua_gettop(L);
  if(lua_type(L,len)==LUA_TTABLE)
  {
    /* psi * table */
    lua_newtable(L);
    lua_pushnil(L);                /* push a nil to the stack - needed for lua_next */
    while(lua_next(L, len) != 0)   /* push the next index and value to the end of the stack +1 on stack (last element replaced)*/
    {
      LuaWavefunctionMul(L);       /* multimply the value with the wavefunction (first times last of the stack) +1 on stack*/
      lua_pushvalue(L, -3);        /* push the index to the end of the stack +1 on stack*/
      lua_insert(L, -2);           /* swap index and value */
      lua_settable(L, -5);         /* set new table value -2 on stack */
      lua_pop(L, 1);               /* remove one element from the stack -1 on stack */
    }
    return 1;
  }
  if(lua_type(L,1)==LUA_TTABLE)
  {
    /* psi * table */
    lua_newtable(L);
    lua_pushnil(L);                /* push a nil to the stack - needed for lua_next */
    while(lua_next(L, 1) != 0)     /* push the next index and value to the end of the stack +1 on stack (last element replaced)*/
    {
      lua_insert(L, 1);            /* set value from table to first element on stack */
      lua_pushvalue(L, len+1);     /* push the wavefunctoin to the end of the stack +1 on stack*/
      LuaWavefunctionMul(L);       /* multimply the value with the wavefunction (first times last of the stack) +1 on stack*/
      lua_pushvalue(L, -3);        /* push the index to the end of the stack +1 on stack*/
      lua_insert(L, -2);           /* swap index and value */
      lua_settable(L, -5);         /* set new table value -2 on stack */
      lua_pop(L, 1);               /* remove one element from the stack -1 on stack */
      lua_remove(L, 1);            /* remove first element from the stack -1 on stack */
    }
    return 1;
  }
  if (lua_isnumber(L, 1))
  {
    double b = (double) luaL_checknumber(L, 1);
    WaveFunctionType *psi = (WaveFunctionType*) luaL_checkudata(L, len, "Wavefunction_Type");
    WaveFunctionType *psiNew;
    psiNew = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
    WaveFunctionCopy(*psi, psiNew);
    WaveFunctionMultiply(psiNew, b, 0);
    luaL_setmetatable(L, "Wavefunction_Type");
    return 1;
  }
  if (lua_isnumber(L, len))
  {
    double b = (double) luaL_checknumber(L, len);
    WaveFunctionType *psi = (WaveFunctionType*) luaL_checkudata(L, 1, "Wavefunction_Type");
    WaveFunctionType *psiNew;
    psiNew = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
    WaveFunctionCopy(*psi, psiNew);
    WaveFunctionMultiply(psiNew, b, 0);
    luaL_setmetatable(L, "Wavefunction_Type");
    return 1;
  }
  if (luaL_testudata(L, 1, "Complex_Type") != NULL)
  {
    Complex16 *b = (Complex16 *) luaL_checkudata(L, 1, "Complex_Type");
    WaveFunctionType *psi = (WaveFunctionType*) luaL_checkudata(L, len, "Wavefunction_Type");
    WaveFunctionType *psiNew;
    psiNew = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
    WaveFunctionCopy(*psi, psiNew);
    WaveFunctionMultiply(psiNew, b->Re, b->Im);
    luaL_setmetatable(L, "Wavefunction_Type");
    return 1;
  }
  if (luaL_testudata(L, len, "Complex_Type") != NULL)
  {
    Complex16 *b = (Complex16 *) luaL_checkudata(L, len, "Complex_Type");
    WaveFunctionType *psi = (WaveFunctionType*) luaL_checkudata(L, 1, "Wavefunction_Type");
    WaveFunctionType *psiNew;
    psiNew = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
    WaveFunctionCopy(*psi, psiNew);
    WaveFunctionMultiply(psiNew, b->Re, b->Im);
    luaL_setmetatable(L, "Wavefunction_Type");
    return 1;
  }
  if (luaL_testudata(L, 1, "Operator_Type") != NULL)
  {
    /* O | psi> */
    WaveFunctionType *psi = (WaveFunctionType*) luaL_checkudata(L, len, "Wavefunction_Type");
    OperatorType *O = (OperatorType *) luaL_checkudata(L,1,"Operator_Type");
    WaveFunctionType *Opsi;
    Opsi = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
    OperatorPsi(*O, *psi, Opsi);
    luaL_setmetatable(L, "Wavefunction_Type");
    return 1;
  }
  if (luaL_testudata(L, len, "Operator_Type") != NULL)
  {
    /* <psi | O */
    WaveFunctionType *psi = (WaveFunctionType*) luaL_checkudata(L, 1, "Wavefunction_Type");
    OperatorType *O = (OperatorType *) luaL_checkudata(L,len,"Operator_Type");
    OperatorType OT;
    WaveFunctionType *Opsi;
    OperatorConjugateTransposeNormalOrder(*O, &OT, &QDetNotOrderedFcFaBaBc);
    Opsi = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
    luaL_setmetatable(L, "Wavefunction_Type");
    OperatorPsi(OT, *psi, Opsi);
    FreeOperator(&OT);
    return 1;
  }
  {
    /* <psi | psi> */
    WaveFunctionType *l = (WaveFunctionType *) luaL_checkudata(L, 1, "Wavefunction_Type");
    WaveFunctionType *r = (WaveFunctionType *) luaL_checkudata(L, len, "Wavefunction_Type");
    double DotRe, DotIm;
    WaveFunctionDotProduct(*l, *r, &DotRe, &DotIm);
    if(Abs(DotIm)<10*(DBL_EPSILON*Abs(DotRe)+DBL_MIN))
    {
      lua_pushnumber(L, DotRe);
    }
    else
    {
      Complex16 *res;
      res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
      res->Re = DotRe;
      res->Im = DotIm;
      luaL_setmetatable(L, "Complex_Type");
    }
  }
  return 1;
}



static int LuaWavefunctionPrint(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len==0)
  {
    PrintWaveFunctionReadable(*(WaveFunctionType *) currentudata);
  }
  else if(lua_istable(L,-1))
  {
    printf("Options for print recognized but not implemented\n");
    PrintWaveFunctionFullReadable(*(WaveFunctionType *) currentudata);
  }
  else
  {
    return luaL_error(L, "Error in .Print(arg) unrecognized arguments in Print()\n");
  }
  return 0;
}

static int LuaWavefunctionChop(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  double border;
  switch (len)
  {
    case 0:
      border=10*DBL_EPSILON;
      break;
    case 1:
      border=luaL_checknumber(L, 1);
      break;
    default:
      return luaL_error(L, "Errorr, call to chop with %d arguments, expected 0 or 1\n",len);
  }
  WaveFunctionTruncateBasis((WaveFunctionType *) currentudata, border);
  return 0;
}

static int LuaWavefunctionConjugate(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len==0)
  {
    WaveFunctionConjugate(*(WaveFunctionType *)currentudata);
  }
  else
  {
    return luaL_error(L, "Error in .Conjugate(arg) unrecognized arguments in Conjugate()\n");
  }
  return 0;
}

static int LuaWavefunctionNorm(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len==0)
  {
    WaveFunctionNormalize((WaveFunctionType *) currentudata);
  }
  else
  {
    return luaL_error(L, "Error in .Norm(arg) unrecognized arguments in Norm()\n");
  }
  return 0;
}

static int LuaWavefunctionRandomize(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  if(len==0)
  {
    WaveFunctionRandomize((WaveFunctionType *) currentudata);
  }
  else
  {
    return luaL_error(L, "Error in .Randomize(arg) unrecognized arguments in Randomize()\n");
  }
  return 0;
}

static int LuaWavefunctionIndex(lua_State * L)
{
  WaveFunctionType *a = (WaveFunctionType *) luaL_checkudata(L, 1, "Wavefunction_Type");
  const char *str = lua_tostring(L, 2);
  currentudata = a;
  if (strcmp(str, "Print") == 0)
  {
    lua_pushcfunction(L, LuaWavefunctionPrint);
    return 1;
  }
  if (strcmp(str, "Chop") == 0)
  {
    lua_pushcfunction(L, LuaWavefunctionChop);
    return 1;
  }
  if ((strcmp(str, "Conjugate") == 0)||(strcmp(str, "C")==0))
  {
    lua_pushcfunction(L, LuaWavefunctionConjugate);
    return 1;
  }
  if (strcmp(str, "Norm") == 0)
  {
    lua_pushcfunction(L, LuaWavefunctionNorm);
    return 1;
  }
  if (strcmp(str, "Randomize") == 0)
  {
    lua_pushcfunction(L, LuaWavefunctionRandomize);
    return 1;
  }
  if (strcmp(str, "Name") == 0)
  {
    lua_pushstring(L, a->Name);
    return 1;
  }
  if (strcmp(str, "NF") == 0)
  {
    lua_pushinteger(L, a->NFermion);
    return 1;
  }
  if (strcmp(str, "NB") == 0)
  {
    lua_pushinteger(L, a->NBoson);
    return 1;
  }
  if (strcmp(str, "Complex") == 0)
  {
    lua_pushboolean(L, a->QComplex);
    return 1;
  }
  if (strcmp(str, "N") == 0)
  {
    lua_pushinteger(L, a->N);
    return 1;
  }
  if (strcmp(str, "NBitsKey") == 0)
  {
    lua_pushinteger(L, a->NBitsKey);
    return 1;
  }
  return luaL_error(L, "Error operation \".%s()\" not defined on wavefunction\n",str);
}


static int LuaWavefunctionNewIndex(lua_State * L)
{
  WaveFunctionType *a = (WaveFunctionType *) luaL_checkudata(L, 1, "Wavefunction_Type");
  const char *str = lua_tostring(L, 2);
  // NBitsKey
  if(strcmp(str, "Name")==0)
  {
    WavefunctionSetName(a, lua_tostring(L, 3));
  }
  else
  {
    return luaL_error(L, "Error index %s not found in Wavefunction Object\n",str);
  }
  return 0;
}

static int LuaWavefunctionToString(lua_State * L)
{
  WaveFunctionType *res = (WaveFunctionType *) luaL_checkudata(L, 1, "Wavefunction_Type");
  char *Buffer;
  StringPrintWaveFunctionReadable(*res, &Buffer);
  lua_pushstring(L, Buffer);
  return 1;
}

static int LuaWavefunctionDestroy(lua_State * L)
{
  WaveFunctionType *psi = (WaveFunctionType*) luaL_checkudata(L, 1, "Wavefunction_Type");
  WaveFunctionFree(psi);
  return 0;
}

static const struct luaL_Reg Wavefunction_meta[] = {
  {"__add", LuaWavefunctionAdd},
  {"__sub", LuaWavefunctionSub},
  {"__mul", LuaWavefunctionMul},
  {"__index", LuaWavefunctionIndex},
  {"__newindex", LuaWavefunctionNewIndex},
  {"__tostring", LuaWavefunctionToString},
  {"__gc", LuaWavefunctionDestroy},
  {NULL, NULL}
};
