static int LuaWavefunctionNew(lua_State *L)
{
  Complex16 temp;
  uint32_t QComplex=0,NDet,i,detleng;
  uint8_t *det;
	WaveFunctionType *psi;
  uint32_t len = (uint32_t)lua_gettop(L);
	psi = (WaveFunctionType*) lua_newuserdata(L, sizeof(WaveFunctionType));
  luaL_setmetatable(L, "Lua.wavefunction");
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

static int LuaWavefunctionDestroy(lua_State * L)
{
  WaveFunctionType *psi = (WaveFunctionType*) luaL_checkudata(L, 1, "Lua.wavefunction");
  WaveFunctionFree(psi);
  return 0;
}


static int LuaWavefunctionPrintReadable(lua_State * L)
{
  PrintWaveFunctionReadable(*(WaveFunctionType *) currentudata);
  return 0;
}

/***********************/
/*   lua interfacing   */
/***********************/


static int LuaWavefunctionIndex(lua_State * L)
{
  WaveFunctionType *a = (WaveFunctionType *) luaL_checkudata(L, 1, "Lua.wavefunction");
  const char *str = lua_tostring(L, 2);
  currentudata = a;
  if (strcmp(str, "Print") == 0)
  {
    lua_pushcfunction(L, LuaWavefunctionPrintReadable);
    return 1;
  }
  return luaL_error(L, "Error wavefunction variable not found");;
}

static const struct luaL_Reg Wavefunction_meta[] = {
  {"__index", LuaWavefunctionIndex},
  {"__gc", LuaWavefunctionDestroy},
  {NULL, NULL}
};

static const struct luaL_Reg Wavefunction_func[] = {
  {"New", LuaWavefunctionNew},
  {NULL, NULL}
};

int luaopen_Wavefunction(lua_State * L)
{
  luaL_newmetatable(L, "Lua.wavefunction");
  luaL_setfuncs(L, Wavefunction_meta, 0);
  luaL_newlib(L, Wavefunction_func);
  lua_setglobal(L, "Wavefunction");
  return 1;
}