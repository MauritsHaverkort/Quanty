/* Some makeshift functions */


void debug(string debugMessage = "")
{
	cout << endl << endl << "DEBUG:" << debugMessage << endl << endl;
	fflush(stdout);
}

void debug(char * debugMessage)
{
	printf("\n\nDEBUG:");
	printf("%s", debugMessage);
	printf("\n\n");
}

//void debug(int i)
//{
//	debug(std::to_string(i));
//}
//
//void debug(double d)
//{
//	debug(std::to_string(d));
//}



void PartialMeanFieldEigensystem()
{
	//TODO
}

void EigensystemSimon()
{
	//TODO: ruft partialMeanFieldEigensystem auf, wenn die Restrictions dies diktieren.
}

static int LuaPartialMeanFieldEigensystem(lua_State *L)
{
	//TODO
	return 0;
}

//initialises and modifies a list of indices that can, without any approximation, be treated on a mean field level.
void mandatoryMeanFieldList(uint16_t **indices, int * const indSetSize, const DeterminantRestrictionType det, int NFermion, int NBoson) //TODO: rewrite this in C++, using vectors
{
	uint32_t NWord = NBoson + Min(NFermion, 1 + (NFermion - 1) / 8);
	*indices = (uint16_t*) malloc(det.N * NWord * sizeof(uint16_t));
	int i;
	char * detStr;
	for (i = 0; i < det.N; ++i)
	{
		if (det.Min[i] == det.Max[i])
		{
			detStr = DetToString(NFermion, NBoson, det.Determinant[i]);
			debug(detStr);
//		TODO: add every index with 1 to a temporary vector.
//		TODO: compare number of entries to number of max and min electrons. If equal, add those indices to the return vector that are not already in it.
		}
	}
}



static int LuaSimontest(lua_State * L)
{
	DeterminantRestrictionType testDet = *(DeterminantRestrictionType*) luaL_checkarray_DeterminantRestrictionType(L, 1);
	uint16_t * indList;
	int listLength;

//	char * detStr = DetToString(testDet.Determinant[1]);

	mandatoryMeanFieldList(&indList, &listLength, testDet, 3, 0);
	return 0;
}

// accepts a list of Wave Functions as well as the Hamiltonian and a list of other operators. The Hamiltonian is special here because it can be used to calculate the Deviance. It is optional.
static int LuaPrintExpectationValues(lua_State * L)
{
	WaveFunctionType *psiList;
	OperatorType *hamil, *opList;
	uint32_t numberOfWavefunctions, numberOfOperators;
	uint32_t len = (uint32_t) lua_gettop(L);
	if (len < 2 || len > 3)
		return luaL_error(L, "PrintExpectationValues called with %d arguments. 2 or 3 expected, List of Wavefunctions, optional Hamiltonian, and list of other operators.\n", len);
	if (lua_istable(L, 1))
	{
		psiList = luaL_checkarray_WaveFunctionType(L, 1, &numberOfWavefunctions);
	} else
	{
		numberOfWavefunctions = 1;
		psiList = (WaveFunctionType *) luaL_checkudata(L, 1, "Wavefunction_Type");
	}
	if (len == 2)
		hamil = (OperatorType*) luaL_checkudata(L, 2, "Operator_Type");
	if (lua_istable(L, -1))
	{
		opList = luaL_checkarray_OperatorType(L, -1, &numberOfOperators);
	} else
	{
		numberOfOperators = 1;
		opList = (OperatorType *) luaL_checkudata(L, -1, "Operator_Type");
	}

//TODO: weiter, falls BraKet nicht mehr komplett auf Lua Ebene stattfindet

	return 0;
}







static const struct luaL_Reg Simon_meta[] = {
  {NULL, NULL}
};

static const struct luaL_Reg Simon_func[] = {
  {"test", LuaSimontest},
	{"PartialMeanFieldEigensystem",LuaPartialMeanFieldEigensystem},
	{"PrintExpectationValues",LuaPrintExpectationValues},
  {NULL, NULL}
};
