static int LuaComplexNew(lua_State * L)
{
  Complex16 *res;
  double re = luaL_checknumber(L, 1);
  double im = luaL_checknumber(L, 2);
  res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
  res->Re = re;
  res->Im = im;
  luaL_setmetatable(L, "Lua.complex");
  return 1;
}

static int LuaComplexI(lua_State * L)
{
  Complex16 *res;
  res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
  res->Re = 0.;
  res->Im = 1.;
  luaL_setmetatable(L, "Lua.complex");
  return 1;
}


static int LuaComplexIndex(lua_State * L)
{
  Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Lua.complex");
  const char *str = lua_tostring(L, 2);
  
  switch (str[0])
  {
    case 'r':
    case 'R':
      lua_pushnumber(L, a->Re);
      return 1;
      
    case 'i':
    case 'I':
      lua_pushnumber(L, a->Im);
      return 1;
  }
  return luaL_error(L, "Error Complex variable not found");;
}



static int LuaComplexNewIndex(lua_State * L)
{
  Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Lua.complex");
  const char *str = lua_tostring(L, 2);
  double v = lua_tonumber(L, 3);
  
  switch (str[0])
  {
    case 'r':
    case 'R':
      a->Re = v;
      return 0;
      
    case 'i':
    case 'I':
      a->Im = v;
      return 0;
  }
  return luaL_error(L, "Error index %s not found in Complex Object\n",str);
}



static int LuaComplexAdd(lua_State * L)
{
  Complex16 *res;
  res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
  if(lua_isnumber(L, 1))
  {
    double a = luaL_checknumber(L, 1);
    Complex16 *b = (Complex16 *) luaL_checkudata(L, 2, "Lua.complex");
    res->Re = a + b->Re;
    res->Im = b->Im;
  }
  else if(lua_isnumber(L,2))
  {
    Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Lua.complex");
    double b = luaL_checknumber(L, 2);
    res->Re = a->Re + b;
    res->Im = a->Im ;
  }
  else
  {
    Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Lua.complex");
    Complex16 *b = (Complex16 *) luaL_checkudata(L, 2, "Lua.complex");
    res->Re = a->Re + b->Re;
    res->Im = a->Im + b->Im;
  }
  luaL_setmetatable(L, "Lua.complex");
  return 1;
}


static int LuaComplexMul(lua_State * L)
{
  if(lua_isnumber(L, 1))
  {
    Complex16 *res;
    res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
    double a = luaL_checknumber(L, 1);
    Complex16 *b = (Complex16 *) luaL_checkudata(L, 2, "Lua.complex");
    res->Re = a * b->Re;
    res->Im = a * b->Im;
  }
  else if(lua_isnumber(L,2))
  {
    Complex16 *res;
    res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
    Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Lua.complex");
    double b = luaL_checknumber(L, 2);
    res->Re = a->Re * b;
    res->Im = a->Im * b;
  }
  else if(luaL_testudata(L, 2,"ssp.operator")!=NULL)
  {
    OperatorType *res = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
    OperatorType *l = (OperatorType *) luaL_checkudata(L, 2, "ssp.operator");
    Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Lua.complex");
    OperatorCopy(res, *l);
    OperatorScale(res, a->Re, a->Im);
    luaL_setmetatable(L, "ssp.operator");
    return 1;
  }
  else
  {
    Complex16 *res;
    res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
    Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Lua.complex");
    Complex16 *b = (Complex16 *) luaL_checkudata(L, 2, "Lua.complex");
    res->Re = a->Re * b->Re - a->Im * b->Im;
    res->Im = a->Im * b->Re + a->Re * b->Im;
  }
  luaL_setmetatable(L, "Lua.complex");
  return 1;
}




static int LuaComplexSub(lua_State * L)
{
  Complex16 *res;
  res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
  if(lua_isnumber(L, 1))
  {
    double a = luaL_checknumber(L, 1);
    Complex16 *b = (Complex16 *) luaL_checkudata(L, 2, "Lua.complex");
    res->Re = a - b->Re;
    res->Im = -b->Im;
  }
  else if(lua_isnumber(L,2))
  {
    Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Lua.complex");
    double b = luaL_checknumber(L, 2);
    res->Re = a->Re - b;
    res->Im = a->Im ;
  }
  else
  {
    Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Lua.complex");
    Complex16 *b = (Complex16 *) luaL_checkudata(L, 2, "Lua.complex");
    res->Re = a->Re - b->Re;
    res->Im = a->Im - b->Im;
  }
  luaL_setmetatable(L, "Lua.complex");
  return 1;
}



static int LuaComplexExp(lua_State * L)
{
  Complex16 *res;
  Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Lua.complex");
  res = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
  res->Re = exp(a->Re) * cos(a->Im);
  res->Im = exp(a->Re) * sin(a->Im);
  luaL_setmetatable(L, "Lua.complex");
  return 1;
}



static int LuaComplexRe(lua_State * L)
{
  if(lua_isnumber(L, 1))
  {
    double a = luaL_checknumber(L, 1);
    lua_pushnumber(L, a);
  }
  else if (luaL_testudata(L, 1, "Lua.complex") != NULL)
  {
    Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Lua.complex");
    lua_pushnumber(L, a->Re);
  }
  else
  {
    return luaL_error(L, "Can only take complex part of (complex) numbers");
  }
  return 1;
}

static int LuaComplexIm(lua_State * L)
{
  if(lua_isnumber(L, 1))
  {
    lua_pushnumber(L, 0);
  }
  else if (luaL_testudata(L, 1, "Lua.complex") != NULL)
  {
    Complex16 *a = (Complex16 *) luaL_checkudata(L, 1, "Lua.complex");
    lua_pushnumber(L, a->Im);
  }
  else
  {
    return luaL_error(L, "Can only take complex part of (complex) numbers");
  }
  return 1;
}


static int LuaComplexToString(lua_State * L)
{
  Complex16 *v = (Complex16 *) luaL_checkudata(L, 1, "Lua.complex");
  lua_pushfstring(L, "(%f + %f I)", v->Re, v->Im);
  return 1;
}


static const struct luaL_Reg Complex_meta[] = {
  {"__index", LuaComplexIndex},
  {"__add", LuaComplexAdd},
  {"__sub", LuaComplexSub},
  {"__mul", LuaComplexMul},
  {"__newindex", LuaComplexNewIndex},
  {"__tostring", LuaComplexToString},
  {NULL, NULL}
};


static const struct luaL_Reg Complex_func[] = {
  {"new", LuaComplexNew},
  {"I", LuaComplexI},
  {"exp", LuaComplexExp},
  {"Re", LuaComplexRe},
  {"re", LuaComplexRe},
  {"Im", LuaComplexIm},
  {"im", LuaComplexIm},
  {NULL, NULL}
};



int luaopen_Complex(lua_State * L)
{
  luaL_newmetatable(L, "Lua.complex");
  luaL_setfuncs(L, Complex_meta, 0);
  luaL_newlib(L, Complex_func);
  lua_setglobal(L, "Complex");
  return 1;
}
