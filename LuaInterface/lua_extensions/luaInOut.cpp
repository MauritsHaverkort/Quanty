//Begin of lua_push, which pushes an arbitrary object onto the stack.

void lua_push(lua_State * L, const double pushkin)
{
	lua_pushnumber(L, pushkin);
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

void lua_push(lua_State * L, const char * pushkin)
{
	lua_pushstring(L, pushkin);
}

void lua_push(lua_State * L, const string pushkin)
{
	lua_push(L, pushkin.c_str());
}

void lua_push(lua_State *L, const OperatorType &pushkin)
{
	OperatorType * O = (OperatorType *) lua_newuserdata(L, sizeof(OperatorType));
	InitOperatorCopySettings(O, pushkin);
	O->QComplex = pushkin.QComplex;
	O->l = pushkin.l;
	luaL_setmetatable(L, "Operator_Type");
}

void lua_push(lua_State *L, const OperatorType* &pushkin)
{
	lua_push(L, *pushkin);
}

void lua_push(lua_State *L, const TightBindingDef &pushkin)
{
	TightBindingDef * TB = (TightBindingDef *) lua_newuserdata(L, sizeof(TightBindingDef));
	TightBindingCopy(pushkin,*TB);
	luaL_setmetatable(L, "TightBinding_Type");
}

void lua_push(lua_State *L, const TightBindingDef* &pushkin)
{
	lua_push(L, *pushkin);
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

void lua_push(lua_State * L, const CompactMatrixType * const pushkin)
{
	LuaPushTable_CompactMatrixType(L, pushkin);
}

void lua_push(lua_State * L, const CompactMatrixType &pushkin)
{
	LuaPushTable_CompactMatrixType(L, &pushkin);
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
	in = (double) luaL_checknumber(L, arg);
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
	in = *(OperatorType *)luaL_checkudata(L, -1, "Operator_Type");
}

void lua_check(lua_State *L, int arg, OperatorType* &in)
{
	in = (OperatorType *)luaL_checkudata(L, -1, "Operator_Type");
}

void lua_check(lua_State *L, int arg, TightBindingDef &in)
{
	in = *(TightBindingDef *) luaL_checkudata(L, arg, "TightBinding_Type");
}

void lua_check(lua_State *L, int arg, TightBindingDef* &in)
{
	in = (TightBindingDef *) luaL_checkudata(L, arg, "TightBinding_Type");
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
	} else
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
	return (bool) lua_isnumber(L, arg);
}

bool luaL_isType(lua_State *L, int arg, const Complex16 &example)
{
	return luaL_isudata(L, arg, "Complex_Type");
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
			if (!luaL_isType(L, -1, vector<double>(), len) && !luaL_isType(L, -1, vector<Complex16>(), len))
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
