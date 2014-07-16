static int LuaWavefunctionAdd(lua_State * L)
{
  uint32_t len = (uint32_t)lua_gettop(L);
  switch (lua_type(L,len))
  {
    case LUA_TNIL:
      return luaL_error(L, "Error cannot add a Nil to an operator\n");
      break;
    case LUA_TBOOLEAN:
      return luaL_error(L, "Error cannot add a true or false to an operator\n");
      break;
    case LUA_TNUMBER:
      return luaL_error(L, "Error cannot add a number to an operator\n");
      break;
    case LUA_TSTRING:
      return luaL_error(L, "Error cannot add a string to an operator\n");
      break;
    case LUA_TLIGHTUSERDATA:
      return luaL_error(L, "Error cannot add a lightuserdata object to an operator\n");
      break;
    case
    LUA_TUSERDATA:
      if (luaL_testudata(L, len, "Lua.complex") != NULL)
      {
        return luaL_error(L, "Error cannot add a complex number to an operator\n");
      }
      if (luaL_testudata(L, len, "ssp.operator") != NULL)
      {
        return luaL_error(L, "Error cannot add a wave function to an operator\n");
      }
      if (luaL_testudata(L, len, "ssp.wavefunction") != NULL)
      {
        WaveFunctionType *l = (WaveFunctionType *) luaL_checkudata(L, 1, "ssp.wavefunction");
        WaveFunctionType *r = (WaveFunctionType *) luaL_checkudata(L, len, "ssp.wavefunction");
        //lua_settop(L, len-2);
        WaveFunctionType *res = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
        luaL_setmetatable(L, "ssp.wavefunction");
        WaveFunctionCopy(*l, res);
        WaveFunctionAddUnionBasis(res, *r);
        return 1;
      }
      return luaL_error(L, "Error cannot add a userdata object of unrecognized type to an operator\n");
      break;
    case LUA_TTABLE:
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
      break;
    default:
      return luaL_error(L, "Error cannot add unrecognized type to operator\n");
      break;
  }
}


static int LuaWavefunctionSub(lua_State * L)
{
  return 0;
}

static int LuaWavefunctionMul(lua_State * L)
{
  return 0;
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
    PrintWaveFunctionReadable(*(WaveFunctionType *) currentudata);
  }
  else
  {
    return luaL_error(L, "Error in .Print(arg) unrecognized arguments in Print()\n");
  }
  return 0;
}

static int LuaWavefunctionIndex(lua_State * L)
{
  WaveFunctionType *a = (WaveFunctionType *) luaL_checkudata(L, 1, "ssp.wavefunction");
  const char *str = lua_tostring(L, 2);
  currentudata = a;
  if (strcmp(str, "Print") == 0)
  {
    lua_pushcfunction(L, LuaWavefunctionPrint);
    return 1;
  }
  return luaL_error(L, "Error wavefunction variable not found");
}

static int LuaWavefunctionToString(lua_State * L)
{
  WaveFunctionType *res = (WaveFunctionType *) luaL_checkudata(L, 1, "ssp.wavefunction");
  lua_pushstring(L, res->Name);
  return 1;
}

static int LuaWavefunctionDestroy(lua_State * L)
{
  WaveFunctionType *psi = (WaveFunctionType*) luaL_checkudata(L, 1, "ssp.wavefunction");
  WaveFunctionFree(psi);
  return 0;
}

static const struct luaL_Reg Wavefunction_meta[] = {
  {"__add", LuaWavefunctionAdd},
  {"__sub", LuaWavefunctionSub},
  {"__mul", LuaWavefunctionMul},
  {"__index", LuaWavefunctionIndex},
  {"__tostring", LuaWavefunctionToString},
  {"__gc", LuaWavefunctionDestroy},
  {NULL, NULL}
};
