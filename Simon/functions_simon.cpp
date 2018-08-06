/* Some makeshift functions */

void debug(string debugMessage = "")
{
	cout << endl << endl << "DEBUG:" << debugMessage << endl << endl;
	fflush (stdout);
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
vector<uint16_t> mandatoryMeanFieldList(const DeterminantRestrictionType det, int NFermion, int NBoson) //TODO: rewrite this in C++, using vectors
{
	uint32_t NWord = NBoson + Min(NFermion, 1 + (NFermion - 1) / 8);
	char * detStr;
	for (int i = 0; i < det.N; ++i)
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

	mandatoryMeanFieldList(testDet, 3, 0);
	return 0;
}





static const struct luaL_Reg Simon_meta[] = {
  {NULL, NULL}
};

static const struct luaL_Reg Simon_func[] = {
  {"test", LuaSimontest},
	{"PartialMeanFieldEigensystem",LuaPartialMeanFieldEigensystem},
  {NULL, NULL}
};
