
// Include the header files for all functions including the definitions of constants, global variables and structures
#include "../Headers/header.h"

// Include the main libraries
#include "../Core/Core.cpp"
#include "../IO/IO.cpp"
#include "../Tight_Binding/Tight_Binding.cpp"
#include "../AtomicOperators/AtomicAux.cpp"

// Include files that should be reoredered into the main libraries
#include "../Model/StandardOperators.cpp"
#include "../Hubbard_Model/CellFunctions.cpp"
#include "../DMFT/SingleSiteFunctions.cpp"
#include "../Yi/functions_yi.cpp"
#include "../Martin/functions_martin.cpp"
#include "../Simon/Simon.cpp"



// Include the Lua interface
#include "lua/includeLua.c"
#include "lua_extensions/tools.cpp"
#include "lua_extensions/luaInOut.cpp"
#include "lua_extensions/TableToArray.cpp"
#include "lua_extensions/get_put.cpp"
#include "units.cpp"
#include "Complex.cpp"
#include "Interpolate.cpp" // Interpolation using Polymorphism
#include "Operator.cpp"
#include "Wavefunction.cpp"
#include "Spectra.cpp"
#include "Matrix.cpp"
#include "TightBinding.cpp"
#include "ResponseFunction.cpp"
#include "Plot.cpp"
#include "PostScript.cpp"
#include "BraKet.cpp"
#include "eval.cpp"
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

  int status,result;
  
  lua_State *L = luaL_newstate();
  luaL_openlibs(L);
  luaopen_units(L);
  luaopen_Complex(L);
  luaopen_ssp(L);

  luaconstants(L);

  lua_gc(L, LUA_GCSETPAUSE, 100); //do not weight between garbage cycles
  //lua_gc(L, LUA_GCGEN, 0); //use generation mode, i.e. visit new data more often
  //lua_gc(L, LUA_GCSETSTEPMUL, 500);
  
  if (argc < 2)
  {
    // run Quanty interactively
    dotty(L);
    status = lua_pcall(L, 2, 1, 0);
    result = lua_toboolean(L, -1);  /* get result */
    finalreport(L, status);
  }
  else
  {
    result = luaL_loadfile(L, argv[1]);
    if (result)
    {
      printf("Error while loading the script:\n%s\n\n", lua_tostring(L, -1));
      lua_close(L);
      return 1;
    }
    if (argc > 2)
    {
      lua_push(L,&(argv[2]),argc-2);
      lua_setglobal(L,"args");
    }
    result = lua_pcall(L, 0, 0, 0);
    if (result)
    {
      fprintf(stderr, "Error while executing the script:\n%s\n\n", lua_tostring(L, -1));
      lua_close(L);
      return 1;
    }
  }
  lua_close(L);
  return 0;
}
