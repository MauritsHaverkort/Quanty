
static int LuaNewOperator(lua_State * L)
{
  if(lua_gettop(L)==0)
  {
    lua_pushfstring (L, "ERROR: LuaNewOperator called with 0 arguments, expect at least one\n");
    lua_error(L);
    return 1;
  }
  if(lua_type(L, 1)!=LUA_TSTRING)
  {
    lua_pushfstring (L, "ERROR: first argument in LuaNewOperator should be a string (use \"Oppname\" as first argument)\n");
    lua_error(L);
    return 1;
  }
  char const *OT =lua_tostring(L, 1);
  if ((strcmp(OT, "Coulomb") == 0)||(strcmp(OT, "U") == 0))
  {
    lua_remove(L, 1);
    return LuaOperatorNewCoulomb(L);
  }
  if (strcmp(OT, "Sx") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewSx(L);
  }
  if (strcmp(OT, "Sy") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewSy(L);
  }
  if (strcmp(OT, "Sz") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewSz(L);
  }
  if (strcmp(OT, "Splus") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewSplus(L);
  }
  if (strcmp(OT, "Smin") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewSmin(L);
  }
  if ((strcmp(OT, "Ssqr") == 0)||(strcmp(OT, "S^2") == 0)||(strcmp(OT, "S2") == 0))
  {
    lua_remove(L, 1);
    return LuaOperatorNewSsqr(L);
  }
  if (strcmp(OT, "Lx") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewLx(L);
  }
  if (strcmp(OT, "Ly") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewLy(L);
  }
  if (strcmp(OT, "Lz") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewLz(L);
  }
  if (strcmp(OT, "Lplus") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewLplus(L);
  }
  if (strcmp(OT, "Lmin") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewLmin(L);
  }
  if ((strcmp(OT, "Lsqr") == 0)||(strcmp(OT, "L^2") == 0)||(strcmp(OT, "L2") == 0))
  {
    lua_remove(L, 1);
    return LuaOperatorNewLsqr(L);
  }
  if (strcmp(OT, "Jx") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewJx(L);
  }
  if (strcmp(OT, "Jy") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewJy(L);
  }
  if (strcmp(OT, "Jz") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewJz(L);
  }
  if (strcmp(OT, "Jplus") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewJplus(L);
  }
  if (strcmp(OT, "Jmin") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewJmin(L);
  }
  if ((strcmp(OT, "Jsqr") == 0)||(strcmp(OT, "J^2") == 0)||(strcmp(OT, "J2") == 0))
  {
    lua_remove(L, 1);
    return LuaOperatorNewJsqr(L);
  }
  if (strcmp(OT, "ldots") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewldots(L);
  }
  if (strcmp(OT, "CF") == 0)
  {
    lua_remove(L, 1);
    return LuaOperatorNewCF(L);
  }
  if ((strcmp(OT, "Create") == 0)||(strcmp(OT, "Cr") == 0))
  {
    lua_remove(L, 1);
    return LuaOperatorNewCreate(L);
  }
  if ((strcmp(OT, "Annihilate") == 0)||(strcmp(OT, "An") == 0))
  {
    lua_remove(L, 1);
    return LuaOperatorNewAnnihilate(L);
  }
  if ((strcmp(OT, "Number") == 0)||(strcmp(OT, "CrAn") == 0)||(strcmp(OT, "n") == 0))
  {
    lua_remove(L, 1);
    return LuaOperatorNewNumber(L);
  }
  
  lua_pushfstring (L, "ERROR: first argument in LuaNewOperator (%s) is not recognized as a known operatortype\n",OT);
  lua_error(L);
  return 1;
}

static int LuaNewWaveFunction(lua_State * L)
{
  Complex16 temp;
  uint32_t QComplex=0,NDet,i,detleng;
  uint8_t *det;
	WaveFunctionType *psi;
  uint32_t len = (uint32_t)lua_gettop(L);
	psi = (WaveFunctionType*) lua_newuserdata(L, sizeof(WaveFunctionType));
  luaL_setmetatable(L, "ssp.wavefunction");
  psi->NFermion=(uint32_t)luaL_checkunsigned(L, 1);
  psi->NBoson  =(uint32_t)luaL_checkunsigned(L, 2);
  psi->QComplex=1;
  if(len==3)
  {
    psi->NBitsKey=16;
  }
  else if(len==4)
  {
    psi->NBitsKey  =(uint32_t)luaL_checkunsigned(L, 4);
  }
  else
  {
    lua_pushfstring (L, "Error in LuaWavefunctionNew, expected 3 or 4 arguments, not %d\n",len);
    lua_error(L);
    return 1;
  }
  WaveFunctionInit(psi);
  luaL_checktype(L, 3, LUA_TTABLE);
  NDet = (uint32_t)lua_objlen(L, 3);
  for(i=1;i<=NDet;i++)
  {
    lua_rawgeti(L, 3, i);
    luaL_checktype(L, -1, LUA_TTABLE);
    if(lua_objlen(L, -1)!=2)
    {
      lua_pushfstring (L, "Error in LuaWavefunctionNew, determinant number %d is not of the form {Determinant,prefactor}\n",i);
      lua_error(L);
      return 1;
    }
    else
    {
      lua_rawgeti(L, -1, 1);
      det=luaL_checkarray_unit8_t(L, -1, &detleng);
      if(detleng!=psi->NWord)
      {
        lua_pushfstring (L, "Error in LuaWavefunctionNew, determinant number %d should be a table of %d words\n",i,psi->NWord);
        lua_error(L);
        return 1;
      }
      lua_pop(L, 1);
      lua_rawgeti(L, -1, 2);
      if(lua_isnumber(L, -1))
      {
        temp.Re=(double)luaL_checknumber(L, -1);
        temp.Im=0;
      }
      else
      {
        QComplex=1;
        temp=*(Complex16 *)luaL_checkudata(L, -1, "Lua.complex");
      }
      lua_pop(L, 1);
      CheckDetPhysical(det, psi->NFermion, psi->NBoson);
      ComplexWaveFunctionAddElement(psi, det, temp.Re, temp.Im);
      free(det);
    }
    lua_pop(L, 1);
  }
  if(QComplex==0) MakeRealWaveFunctionFromComplex(*psi, psi);
  return 1;
}

/*  copy or clone the last element on the lua stack and place it on the stack
 *  works recursively on tables
 *  allocates new memory for user types
 */
static int Luacopy(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  switch (lua_type(L,len))
  {
    case LUA_TNIL:
      lua_pushnil(L);
      break;
    case LUA_TBOOLEAN:
      lua_pushboolean(L,lua_toboolean(L,len));
      break;
    case LUA_TNUMBER:
      lua_pushnumber(L,lua_tonumber(L,len));
      break;
    case LUA_TSTRING:
      lua_pushlstring(L,lua_tostring(L,len),lua_strlen(L,len));
      break;
    case LUA_TLIGHTUSERDATA:
      lua_pushlightuserdata(L,(void*)lua_touserdata(L,len));
      break;
    case
      LUA_TUSERDATA:
      if (luaL_testudata(L, len, "Lua.complex") != NULL)
      {
        Complex16 *a = (Complex16 *) luaL_checkudata(L, len, "Lua.complex");
        Complex16 *res;
        res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
        res->Re = a->Re;
        res->Im = a->Im;
        luaL_setmetatable(L, "Lua.complex");
        return 1;
      }
      if (luaL_testudata(L, len, "ssp.wavefunction") != NULL)
      {
        WaveFunctionType *psi = (WaveFunctionType*) luaL_checkudata(L, len, "ssp.wavefunction");
        WaveFunctionType *psiNew;
        psiNew = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
        WaveFunctionCopy(*psi, psiNew);
        luaL_setmetatable(L, "ssp.wavefunction");
        return 1;
      }
      if (luaL_testudata(L, len, "ssp.operator") != NULL)
      {
        OperatorType *O = (OperatorType*) luaL_checkudata(L, len, "ssp.operator");
        OperatorType *ONew;
        ONew = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
        OperatorCopy(ONew, *O);
        luaL_setmetatable(L, "ssp.operator");
        return 1;
      }
      lua_pushfstring (L, "Error: copying or cloning an object of unrecognised type\n");
      lua_error(L);
      break;
    case LUA_TTABLE:
      lua_newtable(L);
      lua_pushnil(L);
      while(lua_next(L, len) != 0)
      {
        Luacopy(L);
        lua_pushvalue(L, -3);
        lua_insert(L, -2);
        lua_settable(L, -5);
        lua_pop(L, 1);
      }
      break;
    default:
      lua_pushfstring (L, "Error: copying or cloning an object of unrecognised type\n");
      lua_error(L);
      break;
  }
  return 1;
}

/*  copy or clone the last element on the lua stack and place it on the stack
 *  works recursively on tables
 *  allocates new memory for user types
 *  complex conjugates user types that can be conjugated
 */
static int Luaconjugate(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  switch (lua_type(L,len))
  {
    case LUA_TNIL:
      lua_pushnil(L);
      break;
    case LUA_TBOOLEAN:
      lua_pushboolean(L,lua_toboolean(L,len));
      break;
    case LUA_TNUMBER:
      lua_pushnumber(L,lua_tonumber(L,len));
      break;
    case LUA_TSTRING:
      lua_pushlstring(L,lua_tostring(L,len),lua_strlen(L,len));
      break;
    case LUA_TLIGHTUSERDATA:
      lua_pushlightuserdata(L,(void*)lua_touserdata(L,len));
      break;
    case
    LUA_TUSERDATA:
      if (luaL_testudata(L, len, "Lua.complex") != NULL)
      {
        Complex16 *a = (Complex16 *) luaL_checkudata(L, len, "Lua.complex");
        Complex16 *res;
        res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
        res->Re = a->Re;
        res->Im = -a->Im;
        luaL_setmetatable(L, "Lua.complex");
        return 1;
      }
      if (luaL_testudata(L, len, "ssp.wavefunction") != NULL)
      {
        WaveFunctionType *psi = (WaveFunctionType*) luaL_checkudata(L, len, "ssp.wavefunction");
        WaveFunctionType *psiNew;
        psiNew = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
        WaveFunctionCopy(*psi, psiNew);
        WaveFunctionConjugate(*psiNew);
        luaL_setmetatable(L, "ssp.wavefunction");
        return 1;
      }
      if (luaL_testudata(L, len, "ssp.operator") != NULL)
      {
        OperatorType *O = (OperatorType*) luaL_checkudata(L, len, "ssp.operator");
        OperatorType *ONew;
        ONew = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
        OperatorCopy(ONew, *O);
        OperatorComplexConjugate(ONew);
        luaL_setmetatable(L, "ssp.operator");
        return 1;
      }
      lua_pushfstring (L, "Error: copying or cloning an object of unrecognised type\n");
      lua_error(L);
      break;
    case LUA_TTABLE:
      lua_newtable(L);
      lua_pushnil(L);
      while(lua_next(L, len) != 0)
      {
        Luacopy(L);
        lua_pushvalue(L, -3);
        lua_insert(L, -2);
        lua_settable(L, -5);
        lua_pop(L, 1);
      }
      break;
    default:
      lua_pushfstring (L, "Error: copying or cloning an object of unrecognised type\n");
      lua_error(L);
      break;
  }
  return 1;
}

static int Luaconjugatetranspose(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  switch (lua_type(L,len))
  {
    case LUA_TNIL:
      lua_pushnil(L);
      break;
    case LUA_TBOOLEAN:
      lua_pushboolean(L,lua_toboolean(L,len));
      break;
    case LUA_TNUMBER:
      lua_pushnumber(L,lua_tonumber(L,len));
      break;
    case LUA_TSTRING:
      lua_pushlstring(L,lua_tostring(L,len),lua_strlen(L,len));
      break;
    case LUA_TLIGHTUSERDATA:
      lua_pushlightuserdata(L,(void*)lua_touserdata(L,len));
      break;
    case
    LUA_TUSERDATA:
      if (luaL_testudata(L, len, "Lua.complex") != NULL)
      {
        Complex16 *a = (Complex16 *) luaL_checkudata(L, len, "Lua.complex");
        Complex16 *res;
        res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
        res->Re = a->Re;
        res->Im = -a->Im;
        luaL_setmetatable(L, "Lua.complex");
        return 1;
      }
      if (luaL_testudata(L, len, "ssp.wavefunction") != NULL)
      {
        WaveFunctionType *psi = (WaveFunctionType*) luaL_checkudata(L, len, "ssp.wavefunction");
        WaveFunctionType *psiNew;
        psiNew = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
        WaveFunctionCopy(*psi, psiNew);
        WaveFunctionConjugate(*psiNew);
        luaL_setmetatable(L, "ssp.wavefunction");
        return 1;
      }
      if (luaL_testudata(L, len, "ssp.operator") != NULL)
      {
        OperatorType *O = (OperatorType*) luaL_checkudata(L, len, "ssp.operator");
        OperatorType *ONew;
        ONew = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
        OperatorConjugateTranspose(*O, ONew);
        luaL_setmetatable(L, "ssp.operator");
        return 1;
      }
      lua_pushfstring (L, "Error: copying or cloning an object of unrecognised type\n");
      lua_error(L);
      break;
    case LUA_TTABLE:
      lua_newtable(L);
      lua_pushnil(L);
      while(lua_next(L, len) != 0)
      {
        Luacopy(L);
        lua_pushvalue(L, -3);
        lua_insert(L, -2);
        lua_settable(L, -5);
        lua_pop(L, 1);
      }
      break;
    default:
      lua_pushfstring (L, "Error: copying or cloning an object of unrecognised type\n");
      lua_error(L);
      break;
  }
  return 1;
}

static int Luatranspose(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  switch (lua_type(L,len))
  {
    case LUA_TNIL:
      lua_pushnil(L);
      break;
    case LUA_TBOOLEAN:
      lua_pushboolean(L,lua_toboolean(L,len));
      break;
    case LUA_TNUMBER:
      lua_pushnumber(L,lua_tonumber(L,len));
      break;
    case LUA_TSTRING:
      lua_pushlstring(L,lua_tostring(L,len),lua_strlen(L,len));
      break;
    case LUA_TLIGHTUSERDATA:
      lua_pushlightuserdata(L,(void*)lua_touserdata(L,len));
      break;
    case
    LUA_TUSERDATA:
      if (luaL_testudata(L, len, "Lua.complex") != NULL)
      {
        Complex16 *a = (Complex16 *) luaL_checkudata(L, len, "Lua.complex");
        Complex16 *res;
        res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
        res->Re = a->Re;
        res->Im = a->Im;
        luaL_setmetatable(L, "Lua.complex");
        return 1;
      }
      if (luaL_testudata(L, len, "ssp.wavefunction") != NULL)
      {
        WaveFunctionType *psi = (WaveFunctionType*) luaL_checkudata(L, len, "ssp.wavefunction");
        WaveFunctionType *psiNew;
        psiNew = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
        WaveFunctionCopy(*psi, psiNew);
        luaL_setmetatable(L, "ssp.wavefunction");
        return 1;
      }
      if (luaL_testudata(L, len, "ssp.operator") != NULL)
      {
        OperatorType *O = (OperatorType*) luaL_checkudata(L, len, "ssp.operator");
        OperatorType *ONew;
        ONew = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
        OperatorTranspose(*O, ONew);
        luaL_setmetatable(L, "ssp.operator");
        return 1;
      }
      lua_pushfstring (L, "Error: copying or cloning an object of unrecognised type\n");
      lua_error(L);
      break;
    case LUA_TTABLE:
      lua_newtable(L);
      lua_pushnil(L);
      while(lua_next(L, len) != 0)
      {
        Luacopy(L);
        lua_pushvalue(L, -3);
        lua_insert(L, -2);
        lua_settable(L, -5);
        lua_pop(L, 1);
      }
      break;
    default:
      lua_pushfstring (L, "Error: copying or cloning an object of unrecognised type\n");
      lua_error(L);
      break;
  }
  return 1;
}

static const struct luaL_Reg ssp_functions[] = {
  {"NewOperator", LuaNewOperator},
  {"NewWavefunction", LuaNewWaveFunction},
  {"Copy", Luacopy},
  {"Clone", Luacopy},
  {"Conjugate", Luaconjugate},
  {"Conjugatetranspose", Luaconjugatetranspose},
  {"Transpose", Luatranspose},
  {"BraKet", LuaBraKet},           /* function in file BraKet.c */
  {"Eval", LuaEval},               /* function in file eval.c */
  {NULL, NULL}
};


int luaopen_ssp(lua_State * L)
{
  luaL_newmetatable(L, "ssp.operator");
  luaL_setfuncs(L, Operator_meta, 0);
  
  
  luaL_newmetatable(L, "ssp.wavefunction");
  luaL_setfuncs(L, Wavefunction_meta, 0);
  
  luaL_newlib(L, ssp_functions);
  lua_setglobal(L, "ssp");
  
  
  return 1;
}