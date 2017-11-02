
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
