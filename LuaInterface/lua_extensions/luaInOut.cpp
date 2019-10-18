//Begin of lua_push, which pushes an arbitrary object onto the stack.

void lua_push(lua_State * L, const double pushkin)
{
	lua_pushnumber(L, pushkin);
}

void lua_push(lua_State * L, const double *pushkin)
{
	lua_pushnumber(L, *pushkin);
}

void lua_push(lua_State * L, const lua_Integer pushkin)
{
	lua_pushinteger(L,pushkin);
}

void lua_push(lua_State * L, const int pushkin)
{
	lua_push(L, (lua_Integer)pushkin);
}

void lua_push(lua_State * L, const uint pushkin)
{
	lua_push(L, (lua_Integer)pushkin);
}

void lua_push(lua_State * L, const uint16_t pushkin)
{
	lua_push(L, (lua_Integer)pushkin);
}

void lua_push(lua_State * L, const Complex16 &pushkin)
{
	if (pushkin.Im == 0)
	{
		lua_pushnumber(L, pushkin.Re);
	} else
	{
		Complex16 *out = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
		out->Re = pushkin.Re;
		out->Im = pushkin.Im;
		luaL_setmetatable(L, "Complex_Type");
	}
}

void lua_push(lua_State * L, const char * pushkin)
{
	lua_pushstring(L, pushkin);
}

void lua_push(lua_State * L, const string * const pushkin)
{
	lua_push(L, pushkin->c_str());
}

void lua_push(lua_State * L, const string pushkin)
{
	lua_push(L, &pushkin);
}

void lua_push(lua_State *L, const OperatorType * const pushkin)
{
	OperatorType * O = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
	*O = *pushkin;
	luaL_setmetatable(L, "Operator_Type");
}

void lua_push(lua_State *L, const WaveFunctionType * const pushkin)
{
	WaveFunctionType * O = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
	*O = *pushkin;
	luaL_setmetatable(L, "Wavefunction_Type");
}

void lua_push(lua_State *L, const WaveFunctionType pushkin)
{
	WaveFunctionType * O = (WaveFunctionType *) lua_newuserdata(L, sizeof(WaveFunctionType));
	*O = pushkin;
	luaL_setmetatable(L, "Wavefunction_Type");
}

void lua_push(lua_State *L, const OperatorType &pushkin)
{
	lua_push(L, &pushkin);
}

void lua_push(lua_State *L, const TightBindingDef * const pushkin)
{
	TightBindingDef * TB = (TightBindingDef *) lua_newuserdata(L, sizeof(TightBindingDef));
	*TB = *pushkin;
	luaL_setmetatable(L, "TightBinding_Type");
}

void lua_push(lua_State *L, const TightBindingDef &pushkin)
{
	lua_push(L, &pushkin);
}

void lua_push(lua_State *L, const AtomType * const pushkin)
{
	LuaPushAtom(L, pushkin);
}

void lua_push(lua_State *L, const AtomType &pushkin)
{
	LuaPushAtom(L, &pushkin);
}

template<typename T>
void lua_push(lua_State * L, const vector<T> &pushkin)
{
	lua_createtable(L, (int) pushkin.size(), 0);
	for (uint32_t i = 0; i < pushkin.size(); ++i)
	{
		lua_push(L, pushkin[i]);
		lua_rawseti(L, -2, i + 1);
	}
}

template<typename T>
void lua_push(lua_State * L, const T *pushkin, uint32_t n)
{
	lua_createtable(L, (int) n, 0);
	for (uint32_t i = 0; i < n; ++i)
	{
		lua_push(L, pushkin[i]);
		lua_rawseti(L, -2, i + 1);
	}
}

void lua_push(lua_State * L, const CompactMatrixType * const pushkin)
{
	if (pushkin->OnLuaStack)
	{
		CompactMatrixType *B = (CompactMatrixType *) lua_newuserdata(L, sizeof(CompactMatrixType));
		luaL_setmetatable(L, "Matrix_Type");
		*B = *pushkin;
	} else
	{
		LuaPushTable_CompactMatrixType(L, pushkin);
	}
}

void lua_push(lua_State * L, const CompactMatrixType &pushkin)
{
	lua_push(L, &pushkin);
}

void lua_push(lua_State * L, const MatrixType * const pushkin)
{
	LuaPushTable_MatrixType(L, pushkin);
}

void lua_push(lua_State * L, const MatrixType &pushkin)
{
	LuaPushTable_MatrixType(L, &pushkin);
}

template<typename keyType, typename valType>
void lua_push(lua_State * L, const map<keyType, valType> &pushkin)
{
	lua_newtable(L);
	int top = (int) lua_gettop(L);
	typename map<keyType, valType>::const_iterator iter;
	for (iter = pushkin.begin(); iter != pushkin.end(); ++iter)
	{
		lua_push(L, iter->first);
		lua_push(L, iter->second);
		lua_settable(L, top);
	}
}

void lua_push(lua_State * L, const bool pushkin)
{
	if (pushkin)
		lua_pushboolean(L, 1);
	else
		lua_pushboolean(L, 0);
}


//Begin of lua_check, which sets the value of the third element to the object on the stack.

void lua_check(lua_State * L, int arg, string &in)
{
	in = (string) luaL_checkstring(L, arg);
}

void lua_check(lua_State * L, int arg, double &in)
{
	int isnum;
	lua_Number infTest = lua_tonumberx (L, arg, &isnum);
	if (infTest == HUGE_VAL) in = HUGE_VAL;
	else in = (double) luaL_checknumber(L, arg);
}

void lua_check(lua_State * L, int arg, int &in)
{
	in = (int) luaL_checkinteger(L, arg);
}

void lua_check(lua_State * L, int arg, uint16_t &in)
{
	in = (uint16_t) luaL_checkunsigned(L, arg);
}

void lua_check(lua_State * L, int arg, uint32_t &in)
{
	in = (uint32_t) luaL_checkunsigned(L, arg);
}

void lua_check(lua_State *L, int arg, bool &in)
{
	if (!lua_isboolean(L,arg))
		luaL_error(L,"Error: Expected a boolean at position %d",arg);
	in = (bool) lua_toboolean(L,arg);
}

void lua_check(lua_State *L, int arg, CompactMatrixType &in)
{
	in = *luaL_checkCompactMatrix(L, arg);
}

void lua_check(lua_State *L, int arg, CompactMatrixType* &in)
{
	in = luaL_checkCompactMatrix(L, arg);
}

void lua_check(lua_State *L, int arg, WaveFunctionType &in)
{
	in = *(WaveFunctionType *) luaL_checkudata(L, arg, "Wavefunction_Type");
}

void lua_check(lua_State *L, int arg, WaveFunctionType* &in)
{
	in = (WaveFunctionType *) luaL_checkudata(L, arg, "Wavefunction_Type");
}

void lua_check(lua_State *L, int arg, OperatorType &in)
{
	in = *(OperatorType *)luaL_checkudata(L, arg, "Operator_Type");
}

void lua_check(lua_State *L, int arg, OperatorType* &in)
{
	in = (OperatorType *)luaL_checkudata(L, arg, "Operator_Type");
}

void lua_check(lua_State *L, int arg, TightBindingDef &in)
{
	in = *(TightBindingDef *) luaL_checkudata(L, arg, "TightBinding_Type");
}

void lua_check(lua_State *L, int arg, TightBindingDef* &in)
{
	in = (TightBindingDef *) luaL_checkudata(L, arg, "TightBinding_Type");
}

void lua_check(lua_State *L, int arg, Complex16 &in)
{
	if (lua_isnumber(L, arg))
	{
		in.Re = (double) luaL_checknumber(L, -1);
		in.Im = 0;
	} else
	{
		in = *(Complex16 *) luaL_checkudata(L, -1, "Complex_Type");
	}
}

void lua_check(lua_State *L, int arg, AtomType &in)
{
	LuaReadAtom(L, arg, &in);
}

void lua_check(lua_State *L, int arg, AtomType* &in)
{
	LuaReadAtom(L, arg, in);
}

template<typename T>
void lua_check(lua_State * L, int arg, vector<T> &in)
{
	T tmp;
	if (lua_istable(L, arg))
	{
		luaL_checktype(L, arg, LUA_TTABLE);
		uint32_t len = (uint32_t) lua_objlen(L, arg);
		in = vector < T > (len);
		for (int i = 0; i < len; i++)
		{
			lua_rawgeti(L, arg, i + 1);
			lua_check(L, -1, tmp);
			in[i] = tmp;
			lua_pop(L, 1);
		}
	}
  else
	{
		in = vector < T > (1);
		lua_check(L, arg, tmp);
		in[0] = tmp;
	}
}

// Format is { {key1, val1}, {key2, val2} } or, in case of a single value, {key, val}
template<typename keyType, typename valType>
void lua_check(lua_State * L, int arg, map<keyType, valType> &in)
{
	keyType tmpKey;
	valType tmpVal;
	luaL_checktype(L, arg, LUA_TTABLE);
	lua_rawgeti(L, arg, 1);
	bool isSinglePair = !((bool) lua_istable(L, -1));
	lua_pop(L, 1);
	if (isSinglePair)
	{
		int len = (int) lua_objlen(L, arg);
		if (len != 0)
		{
			if (len != 2)
				luaL_error(L, "Error: Map entries need to be 2-element tables in the format {key,value}. Entry has %d elements.\n", len);
			lua_rawgeti(L, arg, 1);
			lua_check(L, -1, tmpKey);
			lua_pop(L, 1);
			lua_rawgeti(L, arg, 2);
			lua_check(L, -1, tmpVal);
			lua_pop(L, 1);
			in.insert( { tmpKey, tmpVal });
		}
	} else
	{
		int len = (int) lua_objlen(L, arg);
		for (int i = 1; i < len + 1; ++i)
		{
			lua_rawgeti(L, arg, i);
			if ((int) lua_objlen(L, -1) != 0)
			{
				if ((int) lua_objlen(L, -1) != 2)
					luaL_error(L, "Error: Map entries need to be 2-element tables in the format {key,value}. Entry %d has %d elements.\n", i, (int) lua_objlen(L, -1));
				lua_rawgeti(L, -1, 1);
				lua_check(L, -1, tmpKey);
				lua_pop(L, 1);
				lua_rawgeti(L, -1, 2);
				lua_check(L, -1, tmpVal);
				lua_pop(L, 1);
				in.insert( { tmpKey, tmpVal });
			}
			lua_pop(L, 1);
		}
	}
}

void lua_check(lua_State *L, int arg, Interpolation::InterpolatingFunction* &in)
{
	in = (Interpolation::InterpolatingFunction*) luaL_checkudata(L, arg, "InterpolatingFunction_Type");
}

//Begin of luaL_isType, which checks if the type of an object on the stack is the same as that of the example given as reference, returning a bool without throwing an error.

bool luaL_isudata(lua_State *L, int ud, const char *tname)
{
	return (luaL_testudata(L, ud, tname) != NULL);
}

bool luaL_isType(lua_State *L, int arg, const OperatorType &example)
{
	return luaL_isudata(L, arg, "Operator_Type");
}

bool luaL_isType(lua_State *L, int arg, const WaveFunctionType &example)
{
	return luaL_isudata(L, arg, "Wavefunction_Type");
}

bool luaL_isType(lua_State *L, int arg, const TightBindingDef &example)
{
	return luaL_isudata(L, arg, "TightBinding_Type");
}

bool luaL_isType(lua_State *L, int arg, const double &example)
{
	int isnum;
	lua_Number infTest = lua_tonumberx (L, arg, &isnum);
	if (infTest == HUGE_VAL) return true;
	else return (bool) lua_isnumber(L, arg);
}

bool luaL_isType(lua_State *L, int arg, const Complex16 &example)
{
	return (luaL_isudata(L, arg, "Complex_Type") || (bool) lua_isnumber(L, arg));
}

bool luaL_isType(lua_State *L, int arg, const CompactMatrixType &example)
{
	if (luaL_isudata(L, arg, "Matrix_Type"))
		return true;
	if (lua_istable(L, arg))
	{
		uint len = (uint) lua_objlen(L, arg);
		if (len == 0)
			return false;

		lua_rawgeti(L, arg, 1);
		uint len_row = (uint) lua_objlen(L, -1);
		lua_pop(L, 1);

		for (uint i = 1; i <= len; ++i)
		{
			lua_rawgeti(L, arg, 1);
			if (!luaL_isType(L, -1, vector<double>(), len_row) && !luaL_isType(L, -1, vector<Complex16>(), len_row))
			{
				lua_pop(L, 1);
				return false;
			}
			lua_pop(L, 1);
		}
		return true;
	}
	return false;
}

bool luaL_isType(lua_State *L, int arg, const string &example)
{
	return (bool) lua_isstring(L, arg);
}

template<typename T>
bool luaL_isType(lua_State *L, int arg, const vector<T> &example)
{
	if (lua_istable(L, arg))
	{
		uint len = (uint) lua_objlen(L, arg);
		for (uint i = 1; i <= len; ++i)
		{
			lua_rawgeti(L, arg, i);
			if (!luaL_isType(L, -1, T()))
			{
				lua_pop(L, 1);
				return false;
			}
			lua_pop(L, 1);
		}
		return true;
	} else
	{
		return luaL_isType(L, arg, T());
	}
}

template<typename T>
bool luaL_isType(lua_State *L, int arg, const vector<T> &example, const uint vecsize)
{
	if (vecsize != (uint) lua_objlen(L, arg))
		return false;
	else
		return luaL_isType(L, arg, example);
}

template<typename keyType, typename valType>
bool luaL_isType(lua_State *L, int arg, const map<keyType, valType> &example)
{
	if (!lua_istable(L, arg))
		return false;
	lua_rawgeti(L, arg, 1);
	bool isSinglePair = !((bool) lua_istable(L, -1));
	lua_pop(L, 1);
	if (isSinglePair)
	{
		lua_rawgeti(L, arg, 1);
		if (!luaL_isType(L, -1, keyType()))
		{
			lua_pop(L, 1);
			return false;
		}
		lua_pop(L, 1);
		lua_rawgeti(L, arg, 2);
		if (!luaL_isType(L, -1, valType()))
		{
			lua_pop(L, 1);
			return false;
		}
		lua_pop(L, 1);
		return true;
	} else
	{
		uint len = (uint) lua_objlen(L, arg);
		for (uint i = 1; i <= len; ++i)
		{
			lua_rawgeti(L, arg, i);

			lua_rawgeti(L, -1, 1);
			if (!luaL_isType(L, -1, keyType()))
			{
				lua_pop(L, 1);
				return false;
			}
			lua_pop(L, 1);
			lua_rawgeti(L, -1, 2);
			if (!luaL_isType(L, -1, valType()))
			{
				lua_pop(L, 1);
				return false;
			}
			lua_pop(L, 1);

			lua_pop(L, 1);
		}
		return true;
	}
}

template<typename keyType, typename valType>
bool luaL_isType(lua_State *L, int arg, const map<keyType, vector<valType> > &example, const uint vecsize)
{
	if (!lua_istable(L, arg))
		return false;
	lua_rawgeti(L, arg, 1);
	bool isSinglePair = !((bool) lua_istable(L, -1));
	lua_pop(L, 1);
	if (isSinglePair)
	{
		lua_rawgeti(L, arg, 2);
		if (vecsize != (uint) lua_objlen(L, -1))
		{
			lua_pop(L, 1);
			return false;
		}
		lua_pop(L, 1);
	} else
	{
		uint len = (uint) lua_objlen(L, arg);
		for (uint i = 1; i <= len; ++i)
		{
			lua_rawgeti(L, arg, i);

			lua_rawgeti(L, -1, 2);
			if (vecsize != (uint) lua_objlen(L, -1))
			{
				lua_pop(L, 1);
				return false;
			}
			lua_pop(L, 1);

			lua_pop(L, 1);
		}
	}
	return luaL_isType(L, arg, example);
}

void lua_check(lua_State *L, int arg, BSpline::TBSpline* &fun)
{
	fun = (BSpline::TBSpline*) luaL_checkudata(L, -1, "InterpolatingFunction_Type");
}
