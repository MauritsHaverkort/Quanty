#include "functions_yi.h"

/* Some makeshift functions */

static int LuaFileReadWavefunction(lua_State * L)
{
  struct
  {
    uint32_t NBitsKey;
    Bool Binary;
  } Options;
  uint32_t len = (uint32_t)lua_gettop(L);
  LuaOptionsType LuaFileReadWavefunctionOptions[]={
    {"NBitsKey","uint32_t",&(Options.NBitsKey)},
    {"Binary",  "bool",  &(Options.Binary)},
    {NULL, NULL, NULL}};
  Options.NBitsKey = 16;
  Options.Binary = False;
  
  if (!((len==1)||(len==2)))
    return luaL_error(L, "Error in FileReadWaveFunction(arg).\n The arguments should at least contain the file name.");
  
  const char *flname = luaL_checkstring(L, 1);
  if (len==2) LuaCheckOptions(L, -1, LuaFileReadWavefunctionOptions);
  lua_settop(L, 0);
  
  WaveFunctionType *psi = NULL;
  psi = (WaveFunctionType*) lua_newuserdata(L, sizeof(WaveFunctionType));
  psi->NBitsKey = Options.NBitsKey;
  if (Options.Binary==True)
    FileReadWaveFunctionPlainBinary(psi, flname);
  else
    FileReadWaveFunctionPlain(psi, flname);
  
  luaL_setmetatable(L, "Wavefunction_Type");
  
  return 1;
}

static const struct luaL_Reg Yi_meta[] = {
  {NULL, NULL}
};

static const struct luaL_Reg Yi_func[] = {
  {"ReadWavefunction", LuaFileReadWavefunction},
  {"readwavefunction", LuaFileReadWavefunction},
  {NULL, NULL}
};
