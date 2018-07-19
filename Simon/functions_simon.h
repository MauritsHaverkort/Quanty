void debug(char* debugMessage);
void PartialMeanFieldEigensystem();
void EigensystemSimon();
static int LuaPartialMeanFieldEigensystem(lua_State *L);
void mandatoryMeanFieldList(uint16_t **indices, int * const indSetSize, const DeterminantRestrictionType det, int NFermion, int NBoson) ;
static int LuaSimontest(lua_State * L);
static int LuaPrintExpectationValues(lua_State * L);
