
static int Luacopy();
static int LuaOperatorExp();
static int RegisterEnergyUnits();
static void luaMWH_StackDump ();

uint32_t LuaGet();
uint32_t LuaPut();

uint32_t LuaGetASCII();
uint32_t LuaPutASCII();


CompactMatrixType* luaL_checkCompactMatrix();
void LuaPushTable_CompactMatrixType();

static int LuaChangeTypeResponseFunction();

/**********************************************************/
/**********    lua_extensions/TableToArray.c     **********/
/**********************************************************/

uint8_t* luaL_checkarray_unit8_t();
uint16_t* luaL_checkarray_unit16_t();
uint32_t* luaL_checkarray_unit32_t();
double* luaL_checkarray_double();
Complex16* luaL_checkarray_complex16();
CompactVectorType luaL_checkVector();
CompactMatrixType* luaL_checkCompactMatrix();
MatrixType* luaL_checkMatrix();
TriDiagonalMatrixType* luaL_checkTriDiagonalMatrixType();
BlockTriDiagonalMatrixType* luaL_checkBlockTriDiagonalMatrixType();
AndersonMatrixType* luaL_checkAndersonMatrixType();
BlockAndersonMatrixType* luaL_checkBlockAndersonMatrixType();
BADoubleTriDiagonalMatrixType * luaL_checkNatrualImpurityOrbitalMatrixType();
BlockNaturalImpurityOrbitalMatrixType * luaL_checkBlockNatrualImpurityOrbitalMatrixType();
ListOfPoles * luaL_checkListOfPoles();
BlockListOfPoles * luaL_checkBlockListOfPoles();
BlockListOfPoles * luaL_checkBlockListOfPolesVariableGrid();
ResponsefunctionType * luaL_checkResponsefunction();
ResponsefunctionType * luaL_checkResponsefunction();
int lua_isBlockGreensfunction();
int lua_isResponsefunction();
WaveFunctionType* luaL_checkarray_WaveFunctionType();
WaveFunctionType** luaL_checkarray_PointerWaveFunctionType();
OperatorType* luaL_checkarray_OperatorType();
void luaL_setarray_DeterminatRestrictionType();
DeterminantRestrictionType* luaL_checkarray_DeterminantRestrictionType();
HarmonicsExpansion* luaL_checkarray_HarmonicsExpansion();


static int LuaMatrixConjugate();
static int LuaMatrixConjugateTranspose();
static int LuaMatrixTranspose();
