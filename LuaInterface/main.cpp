
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

// Include Plot functions
#include "../Plot/Plot.cpp"
#include "../Plot/LuaBinding.cpp"

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
#include "Tensor.cpp"
#include "TightBinding.cpp"
#include "ResponseFunction.cpp"
//#include "Plot.cpp"
//#include "PostScript.cpp"
#include "BraKet.cpp"
#include "eval.cpp"
#include "ssp.cpp"



/**
 * The main function of the Lua interface to Quanty.
 * If Quanty is called without an input file you get an interactive prompt
 * If Quanty is called with an input file the file is run as a script
 */
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
  
  // Initiallize Mersenne Twister random generator
  // note that the standard rand functions use system libraries and are
  // not perfectly homogeneous.
  // lua replaces integers by doubles. For small values this is fine, for larger
  // ranges there are massive problems with random numbers, i.e. the distribution is not flat or homogeneous
  
  srand((unsigned int)time(NULL)); rand(); rand(); rand();
  unsigned long init[4];
  init[0] = (unsigned long) rand();
  init[1] = (unsigned long) rand();
  init[2] = (unsigned long) rand();
  init[3] = (unsigned long) rand();
  mt_genrand_int32();
  mt_init_by_array(init, 4);
  mt_genrand_int32();
  
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
    // the directory of the script, such that a script can find files next to
    // itself independent of the directory from which Quanty was called
    char * ScriptDirectory = MWHFileDirectory(argv[1]);
    if (ScriptDirectory != NULL)
    {
      lua_pushstring(L,ScriptDirectory);
      lua_setglobal(L,"ScriptDirectory");
      free(ScriptDirectory);
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
