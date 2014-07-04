static int LuaOperatorNewSx(lua_State * L)
{
  OperatorType *op;
  int ml;
  int NFermion = luaL_checknumber(L, 1);
  int NShell = luaL_checknumber(L, 2);
  uint16_t LindexUp[7], LindexDn[7];
  for (ml = 0; ml < NFermion / 2; ml++)
  {
    LindexDn[ml] = 2 * ml;
    LindexUp[ml] = 2 * ml + 1;
  }
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorSx(NFermion, NShell, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Lua.operator");
  return 1;
}

static int LuaOperatorNewSy(lua_State * L)
{
  OperatorType *op;
  int ml;
  int NFermion = luaL_checknumber(L, 1);
  int NShell = luaL_checknumber(L, 2);
  
  uint16_t LindexUp[7], LindexDn[7];
  for (ml = 0; ml < NFermion / 2; ml++)
  {
    LindexDn[ml] = 2 * ml;
    LindexUp[ml] = 2 * ml + 1;
  }
  
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorSy(NFermion, NShell, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Lua.operator");
  return 1;
}

static int LuaOperatorNewSz(lua_State * L)
{
  OperatorType *op;
  int NFermion = luaL_checknumber(L, 1);
  int NShell = luaL_checknumber(L, 2);
  int ml;
  uint16_t LindexUp[7], LindexDn[7];
  for (ml = 0; ml < NFermion / 2; ml++)
  {
    LindexDn[ml] = 2 * ml;
    LindexUp[ml] = 2 * ml + 1;
  }
  
  op = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  CreateOperatorSz(NFermion, NShell, LindexUp, LindexDn, op);
  luaL_setmetatable(L, "Lua.operator");
  return 1;
}

static int LuaOperatorMul(lua_State * L)
{
  OperatorType *l = (OperatorType *) luaL_checkudata(L, 1, "Lua.operator");
  OperatorType *res = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  luaL_setmetatable(L, "Lua.operator");
  
  if (lua_isnumber(L, 2))
  {
    double r = luaL_checknumber(L, 2);
    OperatorCopy(res, *l);
    OperatorScale(res, r, 0);
  }
  else
    if (luaL_testudata(L, 2, "Lua.complex") != NULL)
    {
      Complex16 *r = (Complex16 *)luaL_checkudata(L, 2, "Lua.complex");
      OperatorCopy(res, *l);
      OperatorScale(res, r->Re, r->Im);
    } else
    {
      OperatorType *r = (OperatorType *) luaL_checkudata(L, 2, "Lua.operator");
      InitOperatorForMult(res, 2, l, r);
      //OperatorMultiply2(res, *l, *r, 'N', 'N');
      OperatorMultiplyNormalOrder(res, *l, *r, &QDetNotOrderedFcFaBcBa);
    }
  return 1;
}

static int LuaOperatorAdd(lua_State * L)
{
  OperatorType *l = (OperatorType *) luaL_checkudata(L, 1, "Lua.operator");
  OperatorType *r = (OperatorType *) luaL_checkudata(L, 2, "Lua.operator");
  lua_settop(L, 0);
  OperatorType *res = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  luaL_setmetatable(L, "Lua.operator");
  OperatorCopy(res, *l);
  OperatorAdd(res, *r, 1., 0);
  return 1;
}

static int LuaOperatorSub(lua_State * L)
{
  OperatorType *l = (OperatorType *) luaL_checkudata(L, 1, "Lua.operator");
  OperatorType *r = (OperatorType *) luaL_checkudata(L, 2, "Lua.operator");
  OperatorType *res = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
  luaL_setmetatable(L, "Lua.operator");
  OperatorCopy(res, *l);
  OperatorAdd(res, *r, -1., 0);
  return 1;
}


static int LuaOperatorToString(lua_State * L)
{
  OperatorType *res = (OperatorType *) luaL_checkudata(L, 1, "Lua.operator");
  lua_pushstring(L, res->Name);
  return 1;
}

static int LuaOperatorPrintReadable(lua_State * L)
{
  PrintOperatorReadable(*((OperatorType *) currentudata));
  return 0;
}

static int LuaOperatorChop(lua_State * L)
{
  //int chop = luaL_checknumber(L, 1);
  OperatorChop((OperatorType *) currentudata, 10.*DBL_MIN);
  return 0;
}

static int LuaOperatorIndex(lua_State * L)
{
  OperatorType *a = (OperatorType *) luaL_checkudata(L, 1, "Lua.operator");
  const char *str = lua_tostring(L, 2);
  currentudata = a;
  
  if (strcmp(str, "PrintReadable") == 0)
  {
    lua_pushcfunction(L, LuaOperatorPrintReadable);
    return 1;
  }
  else if (strcmp(str, "Chop") == 0)
  {
    lua_pushcfunction(L, LuaOperatorChop);
    return 1;
  }
  return luaL_error(L, "Error operator variable not found");;
}

static int LuaOperatorDestroy(lua_State * L)
{
  OperatorType *a = (OperatorType *) luaL_checkudata(L, 1, "Lua.operator");
  FreeOperator(a);
  return 0;
}


static const struct luaL_Reg Operator_meta[] = {
  {"__index", LuaOperatorIndex},
  {"__add", LuaOperatorAdd},
  {"__sub", LuaOperatorSub},
  {"__mul", LuaOperatorMul},
  {"__tostring", LuaOperatorToString},
  {"__gc", LuaOperatorDestroy},
  {NULL, NULL}
};


static const struct luaL_Reg Operator_func[] = {
  {"NewSx", LuaOperatorNewSx},
  {"NewSy", LuaOperatorNewSy},
  {"NewSz", LuaOperatorNewSz},
  {NULL, NULL}
};

int luaopen_Operator(lua_State * L)
{
  luaL_newmetatable(L, "Lua.operator");
  luaL_setfuncs(L, Operator_meta, 0);
  luaL_newlib(L, Operator_func);
  lua_setglobal(L, "Operator");
  return 1;
}