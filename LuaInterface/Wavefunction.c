static int LuaWavefunctionNewFromOperator(lua_State *L)
{
	WaveFunctionType *psi;
	psi = (WaveFunctionType*) lua_newuserdata(L, sizeof(WaveFunctionType));
	// create
  OperatorType *op = (OperatorType *)luaL_checkudata(L, 1, "Lua.operator");
	CreateRealRandomWaveFunction(psi, op->NFermion, op->NFermion);
  
	//uint32_t WaveFunctionInit(WaveFunctionType * Psi)
  
	luaL_setmetatable(L, "Lua.wavefunction");
  return 0;
}

static int LuaWavefunctionDestroy(lua_State * L)
{
  WaveFunctionType *psi = (WaveFunctionType*) luaL_checkudata(L, 1, "Lua.wavefunction");
  WaveFunctionFree(psi);
  return 0;
}


static const struct luaL_Reg Wavefunction_meta[] = {
  {"__gc", LuaWavefunctionDestroy},
  {NULL, NULL}
};

static const struct luaL_Reg Wavefunction_func[] = {
  {"NewFromOperator", LuaWavefunctionNewFromOperator},
  {NULL, NULL}
};

int luaopen_Wavefunction(lua_State * L)
{
  luaL_newmetatable(L, "Lua.wavefunction");
  luaL_setfuncs(L, Operator_meta, 0);
  
  luaL_newlib(L, Wavefunction_func);
  lua_setglobal(L, "Wavefunction");
  
  return 1;
}