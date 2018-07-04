

#include "../Headers/header.h"

#include "../Core/Core.cpp"
#include "../IO/IO.cpp"

#include "../Model/StandardOperators.cpp"
#include "../Hubbard_Model/CellFunctions.cpp"
#include "../Tight_Binding/Tight_Binding_lib.cpp"
#include "../Tight_Binding/TightBindingToG.cpp"
#include "../Tight_Binding/ImpurityUnfolding.cpp"
#include "../Tight_Binding/StuttgartNMTO.cpp"
#include "../Tight_Binding/DresdenFPLO.cpp"
#include "../DMFT/SingleSiteFunctions.cpp"

#include "lua/includeLua.c"
static void *currentudata = NULL;
#include "lua_extensions/tools.cpp"
#include "lua_extensions/TableToArray.cpp"
#include "lua_extensions/get_put.cpp"
#include "units.cpp"
#include "Complex.cpp"
#include "Operator.cpp"
#include "Wavefunction.cpp"
#include "Spectra.cpp"
#include "Matrix.cpp"
#include "TightBinding.cpp"
#include "ResponseFunction.cpp"
#include "Plot.cpp"
#include "BraKet.cpp"
#include "eval.cpp"
#include "../Yi/functions_yi.cpp"
#include "../Martin/functions_martin.cpp"
#include "../Simon/functions_simon.cpp"
#include "ssp.cpp"



int main(int argc, char *argv[])
{
  // PrintDebugInfo prints additional information to the standard output
  // information printed depends on settings of GlobalVerbosity
  // undef in release version
  // #define PrintDebugInfo
  #undef PrintDebugInfo
  #ifdef PrintDebugInfo
    printf("\n\nCompiled with Debug info on\n\n\n");
  #endif

  TimeStampAll();

  int result;
  if (argc != 2)
  {
    printf("Usage: %s luafile\n", argv[0]);
    return 1;
  }

  lua_State *L = luaL_newstate();
  luaL_openlibs(L);
  luaopen_units(L);
  luaopen_Complex(L);
  luaopen_ssp(L);

  luaconstants(L);

  lua_gc(L, LUA_GCSETPAUSE, 100); //do not weight between garbage cycles
  //lua_gc(L, LUA_GCGEN, 0); //use generation mode, i.e. visit new data more often
  //lua_gc(L, LUA_GCSETSTEPMUL, 500);

  result = luaL_loadfile(L, argv[1]);
  if (result)
  {
    printf("Error while loading the script:\n%s\n\n", lua_tostring(L, -1));
    lua_close(L);
    return 1;
  }
  result = lua_pcall(L, 0, 0, 0);
  if (result)
  {
    fprintf(stderr, "Error while executing the script:\n%s\n\n", lua_tostring(L, -1));
    lua_close(L);
    return 1;
  }
  lua_close(L);
  return 0;
}
