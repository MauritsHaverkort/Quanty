/*
 * If the operator is only ever to be multiplied with wavefunctions that have some empty orbitals, any annihilation term
 * acting on them can be deleted to save time. This is very specific. I need it.
 */
void PrepareOperatorForEmptyOrbitals(OperatorType *O, const vector<uint16_t> &ind, uint minLength)
{
	uint16_t *oadd;
	for (uint length = minLength; length <= O->MaxLength; ++length)
	{
		for (uint i = 0; i < (O->l[length].N); ++i)
		{
			oadd = &(O->l[length].o[i * length]);
			for (uint j = 0; j < length; ++j)
			{
				if (isElementOf(oadd[j], ind))
				{
					if (O->l[length].QComplex)
					{
						O->l[length].V[i].Re = 0;
						O->l[length].V[i].Im = 0;
					}
					else
					{
						O->l[length].VRe[i] = 0;
					}
					break;
				}
			}
		}
	}
	OperatorChop(O, 0);
}

/*
 * If the operator is only ever to be multiplied with wavefunctions that have some empty orbitals, any annihilation term
 * acting on them can be deleted to save time. This is very specific. I need it.
 */
static int LuaPrepareOperatorForEmptyOrbitals(lua_State *L)
{
	OperatorType in, *out;
	vector < uint16_t > inds;
	uint minLength;
	uint32_t len = (uint32_t) lua_gettop(L);
	if (len < 2 || len > 3) return luaL_error(L, "Error: PrepareOperatorForEmptyOrbitals called with %d arguments, 2 or 3 expected, Operator, list of indices and optional minimal length!", len);
	if (!luaL_isType(L, 1, in)) return luaL_error(L, "Error in PrepareOperatorForEmptyOrbitals: First argument should be an operator!");
	if (!luaL_isType(L, 2, inds)) return luaL_error(L, "Error in PrepareOperatorForEmptyOrbitals: Second argument should be an list of indices!");
	lua_check(L, 1, in);
	lua_check(L, 2, inds);
	if (len == 3)
		lua_check(L, 3, minLength);
	else
		minLength = 1;
	out = (OperatorType*) lua_newuserdata(L, sizeof(OperatorType));
	OperatorCopy(out, in);
	PrepareOperatorForEmptyOrbitals(out, inds, minLength);
  luaL_setmetatable(L, "Operator_Type");
	return 1;
}
