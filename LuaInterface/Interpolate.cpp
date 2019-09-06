int qsortcomparedoubles(const void *a, const void *b)
{
  if(*((double*)a)<*((double*)b)) return -1;
  return 1;
}


/*
 *
 */
static int LuaSplineInterpolate(lua_State * L)
{
  uint32_t i;
  Interpolation::CubicSpline *fun;
  double yp0,ypn,*xy;
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
  void* memory = lua_newuserdata(L, sizeof(Interpolation::CubicSpline));
  fun = new(memory)Interpolation::CubicSpline;
  luaL_setmetatable(L, "InterpolatingFunction_Type");

  std::vector<double> u = std::vector<double>(lua_objlen(L, 1));
  std::vector<double> v = std::vector<double>(lua_objlen(L, 1));
  for(i=0;i<lua_objlen(L, 1);i++)
  {
    u[i]=xy[2*i];
    v[i]=xy[2*i+1];
  }
  free(xy);
  if(len==2)
  {
    fun->Create(u,v);
    return 1;
  }
  if(len==4)
  {
    yp0 = (double)luaL_checknumber(L, 3);
    ypn = (double)luaL_checknumber(L, 4);
 //   printf("Error in LuaSplineInterpolate: Hermite boundary not yet implemented\n");
   // return 0;
    fun->Create(u,v,yp0,ypn); // TODO implement that whatever it does
    return 1;
    
  }
  return luaL_error(L, "expected 2 or 4 arguments and not (%d)",len);
}

static int LuaCubicHermiteSplineInterpolate(lua_State * L)
{
	uint32_t i;
	Interpolation::CubicHermiteSpline *fun;
	double *xym;
	uint32_t len = (uint32_t)lua_gettop(L);;
	if(len<2)
	{
	return luaL_error(L, "Error in SplineInterpolate expected at least 2 arguments, not (%d)",len);
	}
	if (len == 2)
	{
		luaL_checktype(L, 1, LUA_TTABLE);
		luaL_checktype(L, 2, LUA_TTABLE);
		if(lua_objlen(L, 1)!=lua_objlen(L, 2))
		{
		return luaL_error(L, "Error, the length of the first table with x values (%d) needs to be the same as the length of the second table with y values (%d)",lua_objlen(L, 1),lua_objlen(L, 2));
		}
		// we need to sort the input list acording to x
		xym=(double*)malloc(2*lua_objlen(L, 1)*sizeof(double));
		for(i=0;i<lua_objlen(L, 1);i++)
		{
		lua_rawgeti(L, 1, i+1);
		xym[2*i]=(double)luaL_checknumber(L, -1);
		lua_pop(L, 1);
		}
		for(i=0;i<lua_objlen(L, 1);i++)
		{
		lua_rawgeti(L, 2, i+1);
		xym[2*i+1]=(double)luaL_checknumber(L, -1);
		lua_pop(L, 1);
		}
		qsort(xym, lua_objlen(L, 1), 2*sizeof(double), &qsortcomparedoubles);
		void* memory = lua_newuserdata(L, sizeof(Interpolation::CubicHermiteSpline));
		fun = new(memory)Interpolation::CubicHermiteSpline;
		luaL_setmetatable(L, "InterpolatingFunction_Type");

		std::vector<double> u = std::vector<double>(lua_objlen(L, 1));
		std::vector<double> v = std::vector<double>(lua_objlen(L, 1));
		for(i=0;i<lua_objlen(L, 1);i++)
		{
		u[i]=xym[2*i];
		v[i]=xym[2*i+1];
		}
		free(xym);
		fun->Create(u,v);
	}
	else
	{
		luaL_checktype(L, 1, LUA_TTABLE);
		luaL_checktype(L, 2, LUA_TTABLE);
		luaL_checktype(L, 3, LUA_TTABLE);
		if(lua_objlen(L, 1)!=lua_objlen(L, 2))
		{
			return luaL_error(L, "Error, the length of the first table with x values (%d) needs to be the same as the length of the second table with y values (%d)",lua_objlen(L, 1),lua_objlen(L, 2));
		}
		if(lua_objlen(L, 1)!=lua_objlen(L, 3))
		{
			return luaL_error(L, "Error, the length of the first table with x values (%d) needs to be the same as the length of the third table with m values (%d)",lua_objlen(L, 1),lua_objlen(L, 3));
		}
		// we need to sort the input list acording to x
		xym=(double*)malloc(3*lua_objlen(L, 1)*sizeof(double));
		for(i=0;i<lua_objlen(L, 1);i++)
		{
			lua_rawgeti(L, 1, i+1);
			xym[3*i]=(double)luaL_checknumber(L, -1);
			lua_pop(L, 1);
		}
		for(i=0;i<lua_objlen(L, 1);i++)
		{
			lua_rawgeti(L, 2, i+1);
			xym[3*i+1]=(double)luaL_checknumber(L, -1);
			lua_pop(L, 1);
		}
		for(i=0;i<lua_objlen(L, 1);i++)
		{
			lua_rawgeti(L, 3, i+1);
			xym[3*i+2]=(double)luaL_checknumber(L, -1);
			lua_pop(L, 1);
		}
		qsort(xym, lua_objlen(L, 1), 3*sizeof(double), &qsortcomparedoubles);
		// add the values to an Interpolating function type
		void* memory = lua_newuserdata(L, sizeof(Interpolation::CubicHermiteSpline));
		fun = new(memory)Interpolation::CubicHermiteSpline;
		luaL_setmetatable(L, "InterpolatingFunction_Type");

		std::vector<double> u = std::vector<double>(lua_objlen(L, 1));
		std::vector<double> v = std::vector<double>(lua_objlen(L, 1));
		std::vector<double> w = std::vector<double>(lua_objlen(L, 1));
		for(i=0;i<lua_objlen(L, 1);i++)
		{
		u[i]=xym[3*i];
		v[i]=xym[3*i+1];
		w[i]=xym[3*i+2];
		}
		free(xym);
		fun->Create(u,v,w);
	}
	return 1;
}
/*
 *
 */
static int LuaInterpolateUnm(lua_State * L)
{
	void* memory;
	Interpolation::InterpolatingFunction *f = (Interpolation::InterpolatingFunction *) luaL_checkudata(L, 1, "InterpolatingFunction_Type");
	double n;
	if (f->type == "CubicSpline")
	{
		memory = lua_newuserdata(L, sizeof(Interpolation::CubicSpline));
		luaL_setmetatable(L, "InterpolatingFunction_Type");
		Interpolation::CubicSpline *fplusg = new(memory)Interpolation::CubicSpline;
		Interpolation::UnaryMinus(*((Interpolation::CubicSpline *)f),*fplusg);
	}
	else if (f->type == "CubicHermiteSpline")
	{
		memory = lua_newuserdata(L, sizeof(Interpolation::CubicHermiteSpline));
		luaL_setmetatable(L, "InterpolatingFunction_Type");
		Interpolation::CubicHermiteSpline *fplusg = new(memory)Interpolation::CubicHermiteSpline;
		Interpolation::UnaryMinus(*((Interpolation::CubicHermiteSpline *)f),*fplusg);
	}
	else
	{
		printf("Error in InterpolateAdd: unrecognized type of Interpolation\n");
		return 0;
	}

  return 1;
}
/*
 *
 */


static int LuaInterpolateAdd(lua_State * L)
{
	void* memory;
    Interpolation::InterpolatingFunction *f = (Interpolation::InterpolatingFunction *) luaL_checkudata(L, 1, "InterpolatingFunction_Type");
    Interpolation::InterpolatingFunction *g = (Interpolation::InterpolatingFunction *) luaL_checkudata(L, 2, "InterpolatingFunction_Type");
    if (f->type == "CubicSpline")
    {
        memory = lua_newuserdata(L, sizeof(Interpolation::CubicSpline));
        luaL_setmetatable(L, "InterpolatingFunction_Type");
        Interpolation::CubicSpline *fplusg = new(memory)Interpolation::CubicSpline;
        Interpolation::Add(*f,*g,*fplusg);
    }
    else if (f->type == "CubicHermiteSpline")
    {
        memory = lua_newuserdata(L, sizeof(Interpolation::CubicHermiteSpline));
        luaL_setmetatable(L, "InterpolatingFunction_Type");
        Interpolation::CubicHermiteSpline *fplusg = new(memory)Interpolation::CubicHermiteSpline;
        Interpolation::Add(*f,*g,*fplusg);
    }
    else
    {
        printf("Error in InterpolateAdd: unrecognized type of Interpolation\n");
        return 0;
    }
  return 1;
}
/*
 *
 */
static int LuaInterpolateSub(lua_State * L)
{
	void* memory;
    Interpolation::InterpolatingFunction *f = (Interpolation::InterpolatingFunction *) luaL_checkudata(L, 1, "InterpolatingFunction_Type");
    Interpolation::InterpolatingFunction *g = (Interpolation::InterpolatingFunction *) luaL_checkudata(L, 2, "InterpolatingFunction_Type");
    if (f->type == "CubicSpline")
    {
        memory = lua_newuserdata(L, sizeof(Interpolation::CubicSpline));
        luaL_setmetatable(L, "InterpolatingFunction_Type");
        Interpolation::CubicSpline *fplusg = new(memory)Interpolation::CubicSpline;
        Interpolation::Sub(*f,*g,*fplusg);
    }
    else if (f->type == "CubicHermiteSpline")
    {
        memory = lua_newuserdata(L, sizeof(Interpolation::CubicHermiteSpline));
        luaL_setmetatable(L, "InterpolatingFunction_Type");
        Interpolation::CubicHermiteSpline *fplusg = new(memory)Interpolation::CubicHermiteSpline;
        Interpolation::Sub(*f,*g,*fplusg);
    }
    else
    {
        printf("Error in InterpolateAdd: unrecognized type of Interpolation\n");
        return 0;
    }
    
  return 1;
}
/*
 *
 */
static int LuaInterpolateMul(lua_State * L)
{
	void* memory;
    Interpolation::InterpolatingFunction *f = (Interpolation::InterpolatingFunction *) luaL_checkudata(L, 1, "InterpolatingFunction_Type");
    Interpolation::InterpolatingFunction *g = (Interpolation::InterpolatingFunction *) luaL_checkudata(L, 2, "InterpolatingFunction_Type");
    if (f->type == "CubicSpline")
    {
        memory = lua_newuserdata(L, sizeof(Interpolation::CubicSpline));
        luaL_setmetatable(L, "InterpolatingFunction_Type");
        Interpolation::CubicSpline *fplusg = new(memory)Interpolation::CubicSpline;
        Interpolation::Mul(*f,*g,*fplusg);
    }
    else if (f->type == "CubicHermiteSpline")
    {
        memory = lua_newuserdata(L, sizeof(Interpolation::CubicHermiteSpline));
        luaL_setmetatable(L, "InterpolatingFunction_Type");
        Interpolation::CubicHermiteSpline *fplusg = new(memory)Interpolation::CubicHermiteSpline;
        Interpolation::Mul(*f,*g,*fplusg);
    }
    else
    {
        printf("Error in InterpolateAdd: unrecognized type of Interpolation\n");
        return 0;
    }
  return 1;
}
/*
 *
 */
static int LuaInterpolateDiv(lua_State * L)
{
	void* memory;
    Interpolation::InterpolatingFunction *f = (Interpolation::InterpolatingFunction *) luaL_checkudata(L, 1, "InterpolatingFunction_Type");
    Interpolation::InterpolatingFunction *g = (Interpolation::InterpolatingFunction *) luaL_checkudata(L, 2, "InterpolatingFunction_Type");
    if (f->type == "CubicSpline")
    {
        memory = lua_newuserdata(L, sizeof(Interpolation::CubicSpline));
        luaL_setmetatable(L, "InterpolatingFunction_Type");
        Interpolation::CubicSpline *fplusg = new(memory)Interpolation::CubicSpline;
        Interpolation::Div(*f,*g,*fplusg);
    }
    else if (f->type == "CubicHermiteSpline")
    {
        memory = lua_newuserdata(L, sizeof(Interpolation::CubicHermiteSpline));
        luaL_setmetatable(L, "InterpolatingFunction_Type");
        Interpolation::CubicHermiteSpline *fplusg = new(memory)Interpolation::CubicHermiteSpline;
        Interpolation::Div(*f,*g,*fplusg);
    }
    else
    {
        printf("Error in InterpolateAdd: unrecognized type of Interpolation\n");
        return 0;
    }
    
  return 1;
}
/*
 *
 */
static int LuaInterpolatePow(lua_State * L)
{
	void* memory;
    Interpolation::InterpolatingFunction *f = (Interpolation::InterpolatingFunction *) luaL_checkudata(L, 1, "InterpolatingFunction_Type");
    double n;
    lua_check(L,2,n);
    if (f->type == "CubicSpline")
    {
        memory = lua_newuserdata(L, sizeof(Interpolation::CubicSpline));
        luaL_setmetatable(L, "InterpolatingFunction_Type");
        Interpolation::CubicSpline *fplusg = new(memory)Interpolation::CubicSpline;
        Interpolation::Pow(*f,n,*fplusg);
    }
    else if (f->type == "CubicHermiteSpline")
    {
        memory = lua_newuserdata(L, sizeof(Interpolation::CubicHermiteSpline));
        luaL_setmetatable(L, "InterpolatingFunction_Type");
        Interpolation::CubicHermiteSpline *fplusg = new(memory)Interpolation::CubicHermiteSpline;
        Interpolation::Pow(*f,n,*fplusg);
    }
    else
    {
        printf("Error in InterpolateAdd: unrecognized type of Interpolation\n");
        return 0;
    }
    
  return 1;
}
/*
 *
 */
static int LuaInterpolateConcat(lua_State * L)
{
	void* memory;
    Interpolation::InterpolatingFunction *f = (Interpolation::InterpolatingFunction *) luaL_checkudata(L, 1, "InterpolatingFunction_Type");
    Interpolation::InterpolatingFunction *g = (Interpolation::InterpolatingFunction *) luaL_checkudata(L, 2, "InterpolatingFunction_Type");
    if (f->type == "CubicSpline")
    {
        memory = lua_newuserdata(L, sizeof(Interpolation::CubicSpline));
        luaL_setmetatable(L, "InterpolatingFunction_Type");
        Interpolation::CubicSpline *fplusg = new(memory)Interpolation::CubicSpline;
        Interpolation::Concat(*f,*g,*fplusg);
    }
    else if (f->type == "CubicHermiteSpline")
    {
        memory = lua_newuserdata(L, sizeof(Interpolation::CubicHermiteSpline));
        luaL_setmetatable(L, "InterpolatingFunction_Type");
        Interpolation::CubicHermiteSpline *fplusg = new(memory)Interpolation::CubicHermiteSpline;
        Interpolation::Concat(*f,*g,*fplusg);
    }
    else
    {
        printf("Error in InterpolateAdd: unrecognized type of Interpolation\n");
        return 0;
    }
  return 1;
}
/*
 *
 */
static int LuaInterpolateEq(lua_State * L)
{
    Interpolation::InterpolatingFunction *f = (Interpolation::InterpolatingFunction *) luaL_checkudata(L, 1, "InterpolatingFunction_Type");
    Interpolation::InterpolatingFunction *g = (Interpolation::InterpolatingFunction *) luaL_checkudata(L, 2, "InterpolatingFunction_Type");

    bool equal = false;
    if (f->type == "CubicSpline" && f->type == g->type)
    {
    	equal = (*((Interpolation::CubicSpline*)f) == (*((Interpolation::CubicSpline*)g)));
    }
    else if (f->type == "CubicHermiteSpline" && f->type == g->type)
    {
    	equal = (*((Interpolation::CubicHermiteSpline*)f) == (*((Interpolation::CubicHermiteSpline*)g)));
    }
    lua_push(L,equal);
  return 1;
}
/*
 *
 */
static int LuaInterpolateLt(lua_State * L)
{
	cout << endl << "Warning: Lesser comparison not yet implemented for interpolating functions. Whatever that's supposed to be." << endl;
  return 0;
}
/*
 *
 */
static int LuaInterpolateLe(lua_State * L)
{
	cout << endl << "Warning: Lesser or equal comparison not yet implemented for interpolating functions. Whatever that's supposed to be." << endl;
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
	uint32_t len = (uint32_t)lua_gettop(L);
	Interpolation::InterpolatingFunction *fun = (Interpolation::InterpolatingFunction *) luaL_checkudata(L, 1, "InterpolatingFunction_Type");
	if (lua_istable(L,2))
	{
		vector<double> x = luaL_checkarray_double(L,2);
		if (len == 3)
		{
			uint32_t k = luaL_checkunsigned(L,3);
			vector<double> y = (*fun)(x,k);
			LuaPushArray_double(L, &y[0],(int)y.size());
		}
		else
		{
			vector<double> y = (*fun)(x);
			LuaPushArray_double(L, &y[0],(int)y.size());
		}
	}
	else
	{
		double x=luaL_checknumber(L, 2);
		if (len == 3)
		{
			uint32_t k = luaL_checkunsigned(L,3);
			lua_pushnumber(L, (*fun)(x,k));
		}
		else
		{
			lua_pushnumber(L, (*fun)(x));
		}
	}

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
  Interpolation::InterpolatingFunction *fun = (Interpolation::InterpolatingFunction *) luaL_checkudata(L, 1, "InterpolatingFunction_Type");
  fun->~InterpolatingFunction();
  return 0;
}

static const struct luaL_Reg InterpolatingFunction_func[] = {
  {"Spline", LuaSplineInterpolate},
  {"CubicHermiteSpline", LuaCubicHermiteSplineInterpolate},
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
