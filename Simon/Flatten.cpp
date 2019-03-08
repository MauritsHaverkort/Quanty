static int LuaFlatten(lua_State * L)
{
	uint32_t len = (uint32_t) lua_gettop(L);
	if (len < 1 || len > 2)
		return luaL_error(L, "Error in Flatten: One or two elements expected, table and optional maximal flatten recursion depth.");
	uint maxRecDepth;
	if (len == 2)
	{
		lua_check(L, 2, maxRecDepth);
		lua_pop(L, 1);
	} else
	{
		maxRecDepth = 1;
	}
	if (lua_istable(L, 1))
	{
		for (uint recDepth = 1; recDepth <= maxRecDepth; ++recDepth)
		{
			lua_newtable(L);
			uint32_t size = lua_objlen(L, recDepth);
			uint ind = 1;
			for (uint i = 1; i <= size; ++i)
			{
				lua_rawgeti(L, recDepth, i);
				if (lua_istable(L, -1))
				{
					uint32_t size2 = lua_objlen(L, -1);
					for (uint j = 1; j <= size2; ++j)
					{
						lua_rawgeti(L, -1, j);
						lua_rawseti(L, recDepth + 1, ind++);
					}
					lua_pop(L, 1);
				} else
				{
					lua_rawseti(L, recDepth + 1, ind++);
				}
			}
		}
	}
	return 1;
}
