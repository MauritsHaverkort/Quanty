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



/* divide a wavefunction by a number or complex number
 *
 */
static int LuaWavefunctionDiv(lua_State * L)
{
  if (lua_isnumber(L, 1))
  {
    return luaL_error(L, "Error, number devided by an operator\n");
  }
  if(lua_type(L,1)==LUA_TTABLE)
  {
    return luaL_error(L, "Error, table devided by an operator\n");
  }
  if (lua_isnumber(L, 2))
  {
    double b = (double) luaL_checknumber(L, 2);
    WaveFunctionType *psi = (WaveFunctionType*) luaL_checkudata(L, 1, "Wavefunction_Type");
    WaveFunctionType *psiNew;
    psiNew = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
    WaveFunctionCopy(*psi, psiNew);
    WaveFunctionMultiply(psiNew, 1/b, 0);
    luaL_setmetatable(L, "Wavefunction_Type");
    return 1;
  }
  if (luaL_testudata(L, 2, "Complex_Type") != NULL)
  {
    Complex16 *b = (Complex16 *) luaL_checkudata(L, 2, "Complex_Type");
    WaveFunctionType *psi = (WaveFunctionType*) luaL_checkudata(L, 1, "Wavefunction_Type");
    WaveFunctionType *psiNew;
    psiNew = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
    WaveFunctionCopy(*psi, psiNew);
    WaveFunctionMultiply(psiNew, b->Re/(Sqr(b->Re)+Sqr(b->Im)), -b->Im/(Sqr(b->Re)+Sqr(b->Im)));
    luaL_setmetatable(L, "Wavefunction_Type");
    return 1;
  }
  return luaL_error(L, "Error, wavefunction devided by an object of unrecognized type\n");
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
    if(QWaveFunctionSameBasis(l, r))
    {
      WaveFunctionDotProductSameBasis(*l, *r, &DotRe, &DotIm);
    }
    else
    {
      WaveFunctionDotProduct(*l, *r, &DotRe, &DotIm);
    }
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

static int LuaWavefunctionUnm(lua_State * L)
{
  WaveFunctionType *res, *a;
  a = (WaveFunctionType *) luaL_checkudata(L, 1, "Wavefunction_Type");
  res = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
  WaveFunctionCopy(*a, res);
  WaveFunctionMultiply(res, -1, 0);
  luaL_setmetatable(L, "Wavefunction_Type");
  return 1;
}

static int LuaWavefunctionPrint(lua_State * L)
{
  struct
  {
    char *stream;
    Bool Full;
    Bool Binary;
  } Options;
  uint32_t len = (uint32_t)lua_gettop(L);
  LuaOptionsType LuaWavefunctionPrintOptions[]={
    {"file", "char*", &(Options.stream)},
    {"Full", "bool", &(Options.Full)},
    {"Binary", "bool", &(Options.Binary)},
    {NULL, NULL, NULL}};
  if(!((len==0)||(len==1)))
  {
    return luaL_error(L, "Error in .Print(arg) unrecognized arguments in Print()\n");
  }
  if(len==0)
  {
    PrintWaveFunctionReadable(*(WaveFunctionType *) currentudata);
    return 0;
  }
  Options.Full=False;
  Options.Binary=False;
  Options.stream=NULL;
  LuaCheckOptions(L, -1, LuaWavefunctionPrintOptions);
  if((Options.Binary==True)&&(Options.Full==True)) return luaL_error(L, "Error in .Print(arg) Binary and full text format are exclusive\n");
  if((Options.Binary==True)&&(Options.stream==NULL)) return luaL_error(L, "Error in .Print(arg) Binary mode only allowed for printing to file\n");
  if(Options.stream==NULL)
  {
    if(Options.Full==True)
    {
      PrintWaveFunctionFullReadable(*(WaveFunctionType *) currentudata);
    }
    else
    {
      PrintWaveFunctionReadable(*(WaveFunctionType *) currentudata);
    }
  }
  else
  {
    if(Options.Full==True)
    {
      FilePrintWaveFunctionFullReadable(*(WaveFunctionType *) currentudata, Options.stream, "w");
    }
    else if(Options.Binary==True)
    {
      FilePrintWaveFunctionPlainBinary(*(WaveFunctionType *) currentudata, Options.stream, "w");
    }
    else
    {
      FilePrintWaveFunctionReadable(*(WaveFunctionType *) currentudata, Options.stream, "w");
    }
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

/*  pushes a new wavefunction on the lua stack
 *  On the stack should be a wavefunction and an integer
 *  The wave-function returned is determinant found at the position of the integer
 */
static int LuaWavefunctionDet(lua_State *L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  WaveFunctionType *a = (WaveFunctionType *)currentudata;
  WaveFunctionType *psi;
  uint32_t idet;
  if(len!=1) return luaL_error(L, "Error in .Det(arg) arg should be an integer between 1 and .N (%d)\n",a->N);
  idet = (uint32_t) luaL_checkinteger(L, 1);
  if((idet>a->N)||(idet==0))return luaL_error(L, "Error in .Det(arg = %d) arg should be an integer between 1 and .N (%d)\n",idet,a->N);
  psi = (WaveFunctionType*) lua_newuserdata(L, sizeof(WaveFunctionType));
  luaL_setmetatable(L, "Wavefunction_Type");
  psi->NFermion=a->NFermion;
  psi->NBoson=a->NBoson;
  psi->QComplex=a->QComplex;
  psi->NBitsKey=a->NBitsKey;
  if(WaveFunctionInit(psi))
  {
    lua_gc(L, LUA_GCCOLLECT, 0);
    if(WaveFunctionInit(psi))
    {
      return luaL_error(L, "Error in .Det(%d) WaveFunctionInit failed out of memory?",idet);
    }
  }
  snprintf(psi->Name,MAXStrLength,"%s.Det(%d)",a->Name,idet);
  if(psi->QComplex)
  {
    ComplexWaveFunctionAddElement(psi, &(a->Determinant[BlokElementDet(idet, a->NWord, 0)]), a->ValueRe[BlokElement(idet)], a->ValueIm[BlokElement(idet)]);
  }
  else
  {
    RealWaveFunctionAddElement(psi, &(a->Determinant[BlokElementDet(idet, a->NWord, 0)]), a->ValueRe[BlokElement(idet)]);
  }
  return 1;
}


/*  pushes the prefactor of a determinant to the stack
 *  On the stack should be a wavefunction and an integer or a string representing the determinant
 *  The (complex) number returned is the prefactor of the given determinant
 */
static int LuaWavefunctionPrefactorDet(lua_State *L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  WaveFunctionType *a = (WaveFunctionType *)currentudata;
  uint32_t idet;
  if(len!=1) return luaL_error(L, "Error in .PrefactorDet(arg) arg should be an integer between 1 and .N (%d) or a string representing the determinant\n",a->N);
  //TODO: read a string representing the determinant
  
  idet = (uint32_t) luaL_checkinteger(L, 1);
  if((idet>a->N)||(idet==0))return luaL_error(L, "Error in .PrefactorDet(arg = %d) arg should be an integer between 1 and .N (%d) or a string representing the determinant\n",idet,a->N);
  if(a->QComplex)
  {
    Complex16 *res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
    res->Re = a->ValueRe[BlokElement(idet)];
    res->Im = a->ValueIm[BlokElement(idet)];
    luaL_setmetatable(L, "Complex_Type");
  }
  else
  {
    lua_pushnumber(L,a->ValueRe[BlokElement(idet)]);
  }
  return 1;
}

/*  On stack:
 *    a wavefunction
 *    a table of integers representing the determinants included in the expansion
 *  Returns:
 *    a table of wavefunctions
 */
static int LuaWavefunctionExpand(lua_State * L)
{
  WaveFunctionType *psi,*psilist;
  uint16_t *index;
  uint32_t Ndet;
  if(lua_gettop(L)!=2) luaL_error(L, "Error in Expand(psi,index) two arguments expected, recieved (%d)\n",lua_gettop(L));
  psi = (WaveFunctionType *) luaL_checkudata(L, 1, "Wavefunction_Type");
  index = luaL_checkarray_unit16_t(L, 2, &Ndet);
  if(Ndet>24) luaL_error(L, "Error in Expand(psi,index). Expanding a wavefunction in determinants yields 2^N new functions for N=%d this would not fit in your memory\n",psi->N);
  psilist = malloc((1<<Ndet)*sizeof(WaveFunctionType));
  if(psilist==NULL)
  {
    luaL_error(L, "Error in Expand(psi,index). Out off memory when expanding a function with %d determinants for % determinants into %d seperate functions\n",psi->N,Ndet,(1<<Ndet));
  }
  if(WaveFunctionExpand(psi, index, Ndet, psilist))
  {
    luaL_error(L, "Error in Expand(psi,index). Out off memory when expanding a function with %d determinants for % determinants into %d seperate functions\n",psi->N,Ndet,(1<<Ndet));
  }
  LuaPushTable_WaveFunctionType(L, psilist, (1<<Ndet));
  free(psilist);
  free(index);
  return 1;
}


static int LuaWavefunctionIndex(lua_State * L)
{
  WaveFunctionType *a = (WaveFunctionType *) luaL_checkudata(L, 1, "Wavefunction_Type");
  const char *str = lua_tostring(L, 2);
  currentudata = a;
  if (strcmp(str, "Chop") == 0)
  {
    lua_pushcfunction(L, LuaWavefunctionChop);
    return 1;
  }
  if (strcmp(str, "Complex") == 0)
  {
    lua_pushboolean(L, a->QComplex);
    return 1;
  }
  if ((strcmp(str, "Conjugate") == 0)||(strcmp(str, "C")==0))
  {
    lua_pushcfunction(L, LuaWavefunctionConjugate);
    return 1;
  }
  if (strcmp(str, "Det") == 0)
  {
    lua_pushcfunction(L, LuaWavefunctionDet);
    return 1;
  }
  if (strcmp(str, "PrefactorDet") == 0)
  {
    lua_pushcfunction(L, LuaWavefunctionPrefactorDet);
    return 1;
  }
  if (strcmp(str, "Print") == 0)
  {
    lua_pushcfunction(L, LuaWavefunctionPrint);
    return 1;
  }
  if (strcmp(str, "Randomize") == 0)
  {
    lua_pushcfunction(L, LuaWavefunctionRandomize);
    return 1;
  }
  if (strcmp(str, "N") == 0)
  {
    lua_pushinteger(L, a->N);
    return 1;
  }
  if (strcmp(str, "Name") == 0)
  {
    lua_pushstring(L, a->Name);
    return 1;
  }
  if (strcmp(str, "NB") == 0)
  {
    lua_pushinteger(L, a->NBoson);
    return 1;
  }
  if (strcmp(str, "NBitsKey") == 0)
  {
    lua_pushinteger(L, a->NBitsKey);
    return 1;
  }
  if (strcmp(str, "NF") == 0)
  {
    lua_pushinteger(L, a->NFermion);
    return 1;
  }
  if (strcmp(str, "Normalize") == 0)
  {
    lua_pushcfunction(L, LuaWavefunctionNorm);
    return 1;
  }
  return luaL_error(L, "Error operation \".%s()\" not defined on wavefunction\n",str);
}


static int LuaWavefunctionNewIndex(lua_State * L)
{
  WaveFunctionType *a = (WaveFunctionType *) luaL_checkudata(L, 1, "Wavefunction_Type");
  const char *str = lua_tostring(L, 2);
  if (strcmp(str, "Complex") == 0)
  {
    return luaL_error(L, "Errorr, you can not change the value Complex of a wavefunction %s.\nIn order to change you have to modify the wavefunction itselves, the value of Complex will follow automatically.",a->Name);
  }
  if (strcmp(str, "Det") == 0)
  {
    return luaL_error(L, "Errorr, you can not change the value of a determiant via indexing for wavefunction %s",a->Name);
  }
  if (strcmp(str, "PrefactorDet") == 0)
  {
    return luaL_error(L, "Errorr, you can not change the value of a determiant via indexing for wavefunction %s",a->Name);
  }
  if (strcmp(str, "N") == 0)
  {
    return luaL_error(L, "Errorr, you can not change the value representing the number of determinants in %s\nIf you add or remove determinants from the wavefunction the number will follow automatically.",a->Name);
  }
  if(strcmp(str, "Name")==0)
  {
    WavefunctionSetName(a, lua_tostring(L, 3));
    return 0;
  }
  if (strcmp(str, "NB") == 0)
  {
    uint32_t NBnew;
    NBnew = (uint32_t)luaL_checknumber(L, 3);
    if(NBnew+a->NFermion>65535) return luaL_error(L, "Maximum number of fermions + bosons allowed is 2^16 = 65535\nTrying to change the number of bosons from %d to %d with %d fermions\n",a->NBoson,NBnew,a->NFermion);
    WaveFunctionChangeNFandNB(a, a->NFermion, NBnew);
    //TODO: allow bosons - remove this line when debuged
    return luaL_error(L, "Current version does not suport Bosons\n");
    return 0;
  }
  if (strcmp(str, "NBitsKey") == 0)
  {
    uint32_t Nbitskey;
    Nbitskey = (uint32_t)luaL_checknumber(L, 3);
    if(WaveFunctionResetKeylist(a, Nbitskey))
    {
      return luaL_error(L, "Internal error in WaveFunctionResetKeylist\n");
    }
    return 0;
  }
  if (strcmp(str, "NF") == 0)
  {
    uint32_t NFnew;
    NFnew = (uint32_t)luaL_checknumber(L, 3);
    if(NFnew+a->NBoson>65535) return luaL_error(L, "Maximum number of fermions + bosons allowed is 2^16 = 65535\nTrying to change the number of fermions from %d to %d with %d bosons\n",a->NFermion,NFnew,a->NBoson);
    WaveFunctionChangeNFandNB(a, NFnew, a->NBoson);
    return 0;
  }
  return luaL_error(L, "Error index %s not found in Wavefunction Object\n",str);
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
  {"__div", LuaWavefunctionDiv},
  {"__unm", LuaWavefunctionUnm},
  {"__index", LuaWavefunctionIndex},
  {"__newindex", LuaWavefunctionNewIndex},
  {"__tostring", LuaWavefunctionToString},
  {"__gc", LuaWavefunctionDestroy},
  {NULL, NULL}
};
