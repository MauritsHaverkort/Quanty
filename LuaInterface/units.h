static int PhysicalConstantsFailed (lua_State *L);
static int AllPhysicalConstantsToString (lua_State *L);
static int PhysicalConstantsToString (lua_State *L);
static int RegisterPhysicalConstants(lua_State *L);
static int EnergyUnitsNewIndex (lua_State *L);
static int EnergyUnitsIndex (lua_State *L);
static int AllEnergyUnitsToString (lua_State *L);
int RegisterEnergyUnits(lua_State *L, char* EnUnits);
int luaopen_units(lua_State * L);
