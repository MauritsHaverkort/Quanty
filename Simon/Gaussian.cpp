 
/*
 *
 */
static int LuaGaussianRadialUnm(lua_State * L)
{
	cout << endl << "Warning: Unary minus not yet implemented for Gaussian functions." << endl;
  return 0;
}
/*
 *
 */
static int LuaGaussianRadialAdd(lua_State * L)
{
	cout << endl << "Warning: Addition not yet implemented for Gaussian functions." << endl;
  return 0;
}
/*
 *
 */
static int LuaGaussianRadialSub(lua_State * L)
{
	cout << endl << "Warning: Subtraction not yet implemented for Gaussian functions." << endl;
  return 0;
}
/*
 *
 */
static int LuaGaussianRadialMul(lua_State * L)
{
	cout << endl << "Warning: Multiplication not yet implemented for Gaussian functions." << endl;
  return 0;
}
/*
 *
 */
static int LuaGaussianRadialDiv(lua_State * L)
{
	cout << endl << "Warning: Division not yet implemented for Gaussian functions." << endl;
  return 0;
}
/*
 *
 */
static int LuaGaussianRadialPow(lua_State * L)
{
	cout << endl << "Warning: Power not yet implemented for Gaussian functions." << endl;
  return 0;
}
/*
 *
 */
static int LuaGaussianRadialConcat(lua_State * L)
{
	cout << endl << "Warning: Concatination not yet implemented for Gaussian functions." << endl;
  return 0;
}
/*
 *
 */
static int LuaGaussianRadialEq(lua_State * L)
{
	cout << endl << "Warning: Equality comparison not yet implemented for Gaussian functions." << endl;
  return 0;
}
/*
 *
 */
static int LuaGaussianRadialLt(lua_State * L)
{
	cout << endl << "Warning: Lesser comparison not yet implemented for Gaussian functions. Whatever that's supposed to be." << endl;
  return 0;
}
/*
 *
 */
static int LuaGaussianRadialLe(lua_State * L)
{
	cout << endl << "Warning: Lesser or equal comparison not yet implemented for Gaussian functions. Whatever that's supposed to be." << endl;
  return 0;
}
/*
 *
 */
static int LuaGaussianRadialMode(lua_State * L)
{
  return 0;
}
/*
 *
 */
static int LuaGaussianRadialCall(lua_State * L)
{
	uint32_t len = (uint32_t) lua_gettop(L);
	GaussianRadial *fun = (GaussianRadial *) luaL_checkudata(L, 1, "GaussianRadialFunction_Type");
	if (lua_istable(L, 2))
	{
		vector<double> x = luaL_checkarray_double(L, 2);
		vector<double> y = (*fun)(x);
		LuaPushArray_double(L, &y[0], (int) y.size());
	}
	else
	{
		double x = luaL_checknumber(L, 2);
		lua_pushnumber(L, (*fun)(x));
	}

	return 1;
}
/*
 *
 */
static int LuaGaussianRadialToString(lua_State * L)
{
  lua_pushstring(L, "Gaussian Radial Function");
  return 1;
}
/*
 *
 */
static int LuaGaussianRadialDestroy(lua_State * L)
{
  GaussianRadial *fun = (GaussianRadial *) luaL_checkudata(L, 1, "GaussianRadialFunction_Type");
  fun->~GaussianRadial();
  return 0;
}

void lua_push(lua_State *L, const GaussianRadial &pushkin)
{
	void *mem = lua_newuserdata(L, sizeof(GaussianRadial));
	GaussianRadial * G = new(mem)GaussianRadial;
	*G = pushkin;
	luaL_setmetatable(L, "GaussianRadialFunction_Type");
}
//
//void lua_push(lua_State *L, const GaussianRadial * const pushkin)
//{
//	GaussianRadial * G = (GaussianRadial *) lua_newuserdata(L, sizeof(GaussianRadial));
//	*G = *pushkin;
//	luaL_setmetatable(L, "GaussianRadialFunction_Type");
//}

void lua_push(lua_State *L, GaussianRadial * const pushkin)
{
	GaussianRadial * G = (GaussianRadial *) lua_newuserdata(L, sizeof(GaussianRadial));
	G = pushkin;
	luaL_setmetatable(L, "GaussianRadialFunction_Type");
}

void lua_check(lua_State *L, int arg, GaussianRadial &in)
{
	in = *(GaussianRadial *) luaL_checkudata(L, arg, "GaussianRadialFunction_Type");
}

void lua_check(lua_State *L, int arg, GaussianRadial * in)
{
	in = (GaussianRadial *) luaL_checkudata(L, arg, "GaussianRadialFunction_Type");
}

static int LuaNewGaussianRadial(lua_State *L)
{
	uint l;
	vector<double> amplitudes, widths;
	lua_check(L, 1, l);
	lua_check(L, 2, amplitudes);
	lua_check(L, 3, widths);
	lua_push(L, GaussianRadial(l, amplitudes, widths));
	return 1;
}

vector<double> GaussianRadialToInterpolationYValsTimesR(const GaussianRadial &g, const vector<double> &grid)
{
	vector<double> yVals(grid.size());
	for (uint i = 0; i < (uint) grid.size(); ++i)
	{
		yVals[i] = g(grid[i]) * grid[i];
	}
	return yVals;
}

static int LuaGaussianRadialToCubicSplineTimesR(lua_State *L)
{
	vector<double> grid;
	uint32_t len = (uint32_t) lua_gettop(L);
	if (len < 1 || len > 2)
	{
		return luaL_error(L, "Error: GaussianRadialToCubicSpline called with %d arguments. 1 or 2 expected, GaussianRadialFunction and optional grid.", len);
	}
	GaussianRadial g;
	lua_check(L, 1, g);
	if (len == 2)
	{
		lua_check(L, 2, grid);
	}
	else
	{
		for (uint i = 0; i < 200; ++i)
		{
			grid.push_back((double) (i*i*i) * 2.5 /1000000.);
		}
	}
	void *memory = lua_newuserdata(L, sizeof(Interpolation::CubicSpline));
	Interpolation::CubicSpline *fun = new (memory) Interpolation::CubicSpline;
	luaL_setmetatable(L, "InterpolatingFunction_Type");
	const vector<double> yVals = GaussianRadialToInterpolationYValsTimesR(g, grid);
	fun->Create(grid, yVals);
	return 1;
}

static int LuaGaussianRadialIndex(lua_State *L)
{
	GaussianRadial g;
	lua_check(L, 1, g);
	string ind;
	lua_check(L, 2, ind);
	if (ind.compare("l") == 0)
	{
		lua_push(L, g.getL());
		return 1;
	}
	else if (ind.compare("amplitudes") == 0)
	{
		lua_push(L, g.getAmplitudes());
		return 1;
	}
	else if (ind.compare("widths") == 0)
	{
		lua_push(L, g.getWidths());
		return 1;
	}
	else
	{
		return luaL_error(L, "Error: GaussianRadialFunctions have no property or method with name \"%s\".", ind.c_str());
	}
}

static int LuaGaussianRadialNewIndex(lua_State *L)
{
	GaussianRadial *g;
	lua_check(L, 1, g);
	string ind;
	lua_check(L, 2, ind);
//	if (ind.compare("l") == 0)
//	{
//		uint l;
//		lua_check(L, 3, l);
//		g->setL(l);
//		return 0;
//	}
//	else if (ind.compare("amplitudes") == 0)
//	{
//		vector<double> amplitudes;
//		lua_check(L, 3, amplitudes);
//		g->setAmplitudes(amplitudes);
//		return 0;
//	}
//	else if (ind.compare("widths") == 0)
//	{
//		vector<double> widths;
//		lua_check(L, 3, widths);
//		g->setWidths(widths);
//		return 0;
//	}
//	else
	{
		return luaL_error(L, "Error: GaussianRadialFunctions have no property of method with name \"%s\".", ind.c_str());
	}
}


static const struct luaL_Reg GaussianRadialFunction_func[] = {
  {"NewGaussianRadialFunction", LuaNewGaussianRadial},
	{"ToInterpolationTimesR", LuaGaussianRadialToCubicSplineTimesR},
  {NULL, NULL}
};

static const struct luaL_Reg GaussianRadialFunction_meta[] = {
  {"__unm", LuaGaussianRadialUnm},
  {"__add", LuaGaussianRadialAdd},
  {"__sub", LuaGaussianRadialSub},
  {"__mul", LuaGaussianRadialMul},
  {"__div", LuaGaussianRadialDiv},
  {"__pow", LuaGaussianRadialPow},
  {"__concat", LuaGaussianRadialConcat},
  {"__eq", LuaGaussianRadialEq},
  {"__lt", LuaGaussianRadialLt},
  {"__le", LuaGaussianRadialLe},
  {"__mode", LuaGaussianRadialMode},
  {"__call", LuaGaussianRadialCall},
  {"__tostring", LuaGaussianRadialToString},
  {"__gc", LuaGaussianRadialDestroy},
	{"__index", LuaGaussianRadialIndex},
	{"__newindex", LuaGaussianRadialNewIndex},
  {NULL, NULL}
};
