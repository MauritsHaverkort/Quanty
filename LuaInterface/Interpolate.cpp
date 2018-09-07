
int qsortcomparedoubles(const void *a, const void *b)
{
  if(*((double*)a)>*((double*)b)) return -1;
  return 1;
}


/*
 *
 */
static int LuaSplineInterpolate(lua_State * L)
{
  uint32_t i;
  InterpolatingFunctionType *fun;
  double yp1,ypn,*xy;
  uint32_t len = (uint32_t)lua_gettop(L);;
  if(len<2)
  {
    return luaL_error(L, "Error in SplineInterpolate expected at least two arguments, not (%d)",len);
  }
  luaL_checktype(L, 1, LUA_TTABLE);
  luaL_checktype(L, 2, LUA_TTABLE);
  if(lua_objlen(L, 1)!=lua_objlen(L, 2))
  {
    return luaL_error(L, "Error, the length of the first table with x values (%d) needs to be the same as the length of the second table with y values (%d)",lua_objlen(L, 1),lua_objlen(L, 2));
  }
  // we need to sort the input list acording to x
  xy=(double*)malloc(2*lua_objlen(L, 1)*sizeof(double));
  for(i=0;i<lua_objlen(L, 1);i++)
  {
    lua_rawgeti(L, 1, i+1);
    xy[2*i]=(double)luaL_checknumber(L, -1);
    lua_pop(L, 1);
  }
  for(i=0;i<lua_objlen(L, 1);i++)
  {
    lua_rawgeti(L, 2, i+1);
    xy[2*i+1]=(double)luaL_checknumber(L, -1);
    lua_pop(L, 1);
  }
  qsort(xy, lua_objlen(L, 1), 2*sizeof(double), &qsortcomparedoubles);
  // add the values to an Interpolating function type
  fun = (InterpolatingFunctionType*) lua_newuserdata(L, sizeof(InterpolatingFunctionType));
  luaL_setmetatable(L, "InterpolatingFunction_Type");
  InitInterpolatingFunctionAsCubicSpline(fun, (uint32_t)lua_objlen(L, 1));
  for(i=0;i<((CubicSplineType*)(fun->f))->N;i++)
  {
    ((CubicSplineType*)(fun->f))->x[i]=xy[2*i];
    ((CubicSplineType*)(fun->f))->y[i]=xy[2*i+1];
  }
  free(xy);
  if(len==2)
  {
    CreateCubicSpline((CubicSplineType*)(fun->f), NAN, NAN);
    return 1;
  }
  if(len==4)
  {
    yp1 = (double)luaL_checknumber(L, 3);
    ypn = (double)luaL_checknumber(L, 4);
    CreateCubicSpline((CubicSplineType*)(fun->f), yp1, ypn);
    return 1;
    
  }
  return luaL_error(L, "expected 2 or 4 arguments and not (%d)",len);
}

/*
 *
 */
static int LuaInterpolateUnm(lua_State * L)
{
  return 0;
}
/*
 *
 */
static int LuaInterpolateAdd(lua_State * L)
{
  return 0;
}
/*
 *
 */
static int LuaInterpolateSub(lua_State * L)
{
  return 0;
}
/*
 *
 */
static int LuaInterpolateMul(lua_State * L)
{
  return 0;
}
/*
 *
 */
static int LuaInterpolateDiv(lua_State * L)
{
  return 0;
}
/*
 *
 */
static int LuaInterpolatePow(lua_State * L)
{
  return 0;
}
/*
 *
 */
static int LuaInterpolateConcat(lua_State * L)
{
  return 0;
}
/*
 *
 */
static int LuaInterpolateEq(lua_State * L)
{
  return 0;
}
/*
 *
 */
static int LuaInterpolateLt(lua_State * L)
{
  return 0;
}
/*
 *
 */
static int LuaInterpolateLe(lua_State * L)
{
  return 0;
}
/*
 *
 */
static int LuaInterpolateMode(lua_State * L)
{
  return 0;
}
/*
 *
 */
static int LuaInterpolateCall(lua_State * L)
{
  InterpolatingFunctionType *fun = (InterpolatingFunctionType *) luaL_checkudata(L, 1, "InterpolatingFunction_Type");
  double x=luaL_checknumber(L, 2);
  lua_pushnumber(L, EvaluateInterpolatingFunction(fun, x));
  return 1;
}
/*
 *
 */
static int LuaInterpolateToString(lua_State * L)
{
  lua_pushstring(L, "Interpolating Function");
  return 1;
}
/*
 *
 */
static int LuaInterpolateDestroy(lua_State * L)
{
  InterpolatingFunctionType *fun = (InterpolatingFunctionType *) luaL_checkudata(L, 1, "InterpolatingFunction_Type");
  FreeInterpolatingFunction(fun);
  return 0;
}

static const struct luaL_Reg InterpolatingFunction_func[] = {
  {"Spline", LuaSplineInterpolate},
  {NULL, NULL}
};

static const struct luaL_Reg InterpolatingFunction_meta[] = {
  {"__unm", LuaInterpolateUnm},
  {"__add", LuaInterpolateAdd},
  {"__sub", LuaInterpolateSub},
  {"__mul", LuaInterpolateMul},
  {"__div", LuaInterpolateDiv},
  {"__pow", LuaInterpolatePow},
  {"__concat", LuaInterpolateConcat},
  {"__eq", LuaInterpolateEq},
  {"__lt", LuaInterpolateLt},
  {"__le", LuaInterpolateLe},
  {"__mode", LuaInterpolateMode},
  {"__call", LuaInterpolateCall},
  {"__tostring", LuaInterpolateToString},
  {"__gc", LuaInterpolateDestroy},
  {NULL, NULL}
};
