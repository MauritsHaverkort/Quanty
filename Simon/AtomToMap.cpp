void lua_push_AtomMap(lua_State * L, const AtomType &A)
{
	lua_newtable(L);

	lua_push(L, "Name");
	lua_push(L, A.Name);
	lua_settable(L, -3);

	lua_push(L, "Position");
	lua_push(L, A.Position, (uint32_t) 3);
	lua_settable(L, -3);

	lua_push(L, "Wyckoff");
	lua_push(L, A.Wyckoff);
	lua_settable(L, -3);

	lua_push(L, "Shells");
	lua_newtable(L);
	for (uint i = 0; i < (uint) A.Nshell; ++i)
	{
		lua_push(L, i + 1);
		lua_newtable(L);

		lua_push(L, "Name");
		lua_push(L, A.ShellName[i]);
		lua_settable(L, -3);

		lua_push(L, "States");
		lua_newtable(L);
		for (uint j = 0; j < (uint) A.shelldegeneracy[i]; ++j)
		{
			lua_push(L, j + 1);
			lua_push(L, A.OrbitalName[i][j]);
			lua_settable(L, -3);
		}
		lua_settable(L, -3);
		lua_settable(L, -3);
	}
	lua_settable(L, -3);
}
