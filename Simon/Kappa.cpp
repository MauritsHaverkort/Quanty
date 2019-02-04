 
bool sorterByKappa(const string orb1, const string orb2)
{
	int kappa1 = kappa(orb1);
	int kappa2 = kappa(orb2);
	if (abs(kappa1) == abs(kappa2))
		return kappa1 < kappa2;
	else
		return abs(kappa1) < abs(kappa2);
}

bool hasNonRelOrb(const vector<string> orbs)
{
	for (int i = 0; i < orbs.size(); ++i)
	{
		if (!isRelOrbital(orbs[i]))
			return true;
	}
	return false;
}

vector<string> orbsSortedByKappa(const vector<string> in)
{
	vector < string > out;
	if (hasNonRelOrb(in))
		out = nonRelToRelOrbs(in);
	else
		out = in;
	sort(out.begin(), out.end(), sorterByKappa);
	return out;
}

static int LuaSortByKappa(lua_State * L)
{
	uint32_t len = (uint32_t) lua_gettop(L);
	if (len != 1)
	{
		return luaL_error(L, "SortByKappa called with %d arguments. Only list of orbitals expected.\n", len);
	}
	vector < string > orbitalsIn;
	lua_check(L, 1, orbitalsIn);
	lua_push(L, orbsSortedByKappa(orbitalsIn));
	return 1;
}
