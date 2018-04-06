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
#include <complex.h>
#include <sys/types.h>
#include <sys/stat.h>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_num_procs() 0
#define omp_get_max_threads() 0
#endif

#include "../Core/StructDef.c"
#include "../Headers/header.h"
#include "../Core/TimeStamp.c"
#include "../Core/Hash.c"
#include "../Core/BasicMath.c"
#include "../IO/Print.c"
#include "../IO/Read.c"
#include "../IO/BitMap.c"
#include "../Core/Error.c"
#include "../Core/WaveFunctionArithmetic.c"
#include "../Core/WaveFunctionPerturbationArithmetic.c"
#include "../Core/OperatorArithmetic.c"
#include "../Core/OperatorWaveFunction.c"
#include "../Core/OperatorWaveFunctionPerturbation.c"
#include "../Core/BlockWaveFunctionArithmetic.c"
#include "../Core/Lanczos.c"
#include "../Core/BlockLanczos.c"
#include "../Core/MeanField.c"
#include "../Core/Spectral.c"
#include "../Core/ResonantSpectra.c"
#include "../Core/DensityMatrix.c"
#include "../Core/Random.c"
#include "../Core/MatrixOperator.c"
#include "../Core/ImpurityLanczos.c"
#include "../Model/StandardOperators.c"
#include "../Hubbard_Model/CellFunctions.c"
#include "../Tight_Binding/Tight_Binding_lib.c"
#include "../Tight_Binding/TightBindingToG.c"
#include "../Tight_Binding/ImpurityUnfolding.c"
#include "../Tight_Binding/StuttgartNMTO.c"
#include "../Tight_Binding/DresdenFPLO.c"
#include "../DMFT/SingleSiteFunctions.c"
#include "../IO/StuttgartNMTO.c"
#include "../IO/DresdenFPLO.c"

#include "Documentation.c"
#include "lua/includeLua.h"
static void *currentudata = NULL;
#include "lua_extensions/tools.c"
#include "lua_extensions/TableToArray.c"
#include "lua_extensions/get_put.c"
#include "units.c"
#include "Complex.c"
#include "Operator.c"
#include "Wavefunction.c"
#include "Spectra.c"
#include "Matrix.c"
#include "TightBinding.c"
#include "ResponseFunction.c"
#include "Plot.c"
#include "BraKet.c"
#include "eval.c"
#include "../Yi/functions_yi.c"
#include "ssp.c"

#include "main.h"


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
