#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <unistd.h>
#include <stdarg.h>
#include <sys/types.h>
#include <sys/stat.h>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_num_procs() 0
#define omp_get_max_threads() 0
#endif

#include "../Headers/header.h"
#include "../Core/TimeStamp.cpp"
#include "../Core/Hash.cpp"
#include "../Core/BasicMath.cpp"
#include "../IO/Print.cpp"
#include "../IO/Read.cpp"
#include "../IO/BitMap.cpp"
#include "../Core/Error.cpp"
#include "../Core/WaveFunctionArithmetic.cpp"
#include "../Core/WaveFunctionPerturbationArithmetic.cpp"
#include "../Core/OperatorArithmetic.cpp"
#include "../Core/OperatorWaveFunction.cpp"
#include "../Core/OperatorWaveFunctionPerturbation.cpp"
#include "../Core/BlockWaveFunctionArithmetic.cpp"
#include "../Core/Lanczos.cpp"
#include "../Core/BlockLanczos.cpp"
#include "../Core/MeanField.cpp"
#include "../Core/Spectral.cpp"
#include "../Core/ResonantSpectra.cpp"
#include "../Core/DensityMatrix.cpp"
#include "../Core/Random.cpp"
#include "../Core/MatrixOperator.cpp"
#include "../Core/ImpurityLanczos.cpp"
#include "../Model/StandardOperators.cpp"
#include "../Hubbard_Model/CellFunctions.cpp"
#include "../Tight_Binding/Tight_Binding_lib.cpp"
#include "../Tight_Binding/TightBindingToG.cpp"
#include "../Tight_Binding/ImpurityUnfolding.cpp"
#include "../Tight_Binding/StuttgartNMTO.cpp"
#include "../Tight_Binding/DresdenFPLO.cpp"
#include "../DMFT/SingleSiteFunctions.cpp"
#include "../IO/StuttgartNMTO.cpp"
#include "../IO/DresdenFPLO.cpp"

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
