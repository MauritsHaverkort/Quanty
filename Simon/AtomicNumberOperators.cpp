
map<string, OperatorType> AtomicNumberOperators(const vector<string> &orbitals, const map<string, vector<string> > &groupings)
{
	map < string, OperatorType > numberOperators;
	map<string, vector<uint16_t> > allIndices = createAtmoicIndicesDict(orbitals, groupings);
	int NF = NumberOfFermionicStates(orbitals);
	typedef map<string, vector<uint16_t> >::iterator iterType;
	for (iterType iter = allIndices.begin(); iter != allIndices.end(); ++iter)
	{
		OperatorType numOp;
		vector<uint16_t> indices = iter-> second;
		int NOpp = (int) indices.size();
		double ones[NOpp], zeroes[NOpp];
		for (int i = 0; i < NOpp; ++i)
		{
			ones[i] = 1.;
			zeroes[i] = 0.;
		}
		CreateOperatorCreateAnnihilateList(NF, &indices[0], &indices[0], ones, zeroes, NOpp, &numOp);
		string name = "N_" + iter->first;
			strcpy(numOp.Name, name.c_str());
		numberOperators.insert( { iter->first, numOp });
	}
	return numberOperators;
}




static int LuaAtomicNumberOperators(lua_State * L)
{
	uint32_t len = (uint32_t) lua_gettop(L);
	if (len < 1 || len > 3)
		return luaL_error(L, "NumberOperator called with %d arguments. 1 or 2 expected, list of orbitals, and optional list of groupings.\n", len);
	vector < string > orbitals;
	lua_check(L, 1, orbitals);
	map<string, vector<string> > groupings;
	if (len == 2)
		lua_check(L, 2, groupings);
	groupings.insert( { "all", orbitals });
	groupings.insert( { "All", orbitals });
	lua_push(L, AtomicNumberOperators(orbitals, groupings));
	return 1;
}
