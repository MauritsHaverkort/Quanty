#include <stdio.h>
#include <iostream>
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
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_num_procs() 0
#define omp_get_max_threads() 0
#endif

uint32_t GlobalTimeOn=0xF;
/*      1111 1111 1111 1111 - 0xFFFF - print all timings
 *      0000 0000 0000 0001 - 0x0001 - main specific
 *      0000 0000 0000 0010 - 0x0002 - main specific
 *      0000 0000 0000 0100 - 0x0004 - main specific
 *      0000 0000 0000 1000 - 0x0008 - main specific
 *      0000 0000 0001 0000 - 0x0010 -
 *      0000 0000 0010 0000 - 0x0020 -
 *      0000 0000 0100 0000 - 0x0040 -
 *      0000 0000 1000 0000 - 0x0080 -
 */

uint32_t GlobalVerbosity=0xF;
/*      1111 1111 1111 1111 - 0xFFFF - print all statments on what the code does
 *      0000 0000 0000 0001 - 0x0001 - main specific
 *      0000 0000 0000 0010 - 0x0002 - main specific
 *      0000 0000 0000 0100 - 0x0004 - main specific
 *      0000 0000 0000 1000 - 0x0008 - main specific
 *      0000 0000 0001 0000 - 0x0010 - Lanczos function start
 *      0000 0000 0010 0000 - 0x0020 - Lanczos function converge info
 *      0000 0000 0100 0000 - 0x0040 -
 *      0000 0000 1000 0000 - 0x0080 -
 *
 *      1111 1111 1111 1111 - 0xFFFF0000 - print all debug statements
 *      0000 0000 0000 0001 - 0x00010000 - print start of each subroutine
 *      0000 0000 0000 0010 - 0x00020000 -
 *      0000 0000 0000 0100 - 0x00040000 -
 *      0000 0000 0000 1000 - 0x00080000 -
 *      0000 0000 0001 0000 - 0x00100000 -
 *      0000 0000 0010 0000 - 0x00200000 -
 *      0000 0000 0100 0000 - 0x00400000 -
 *      0000 0000 1000 0000 - 0x00800000 -
 *
 *      0000 0001 0000 0000 - 0x01000000 - used for debuging of functions during development
 */

#define PrintDebugInfo

using std::cout;
using std::endl;
using std::string;
using std::vector;


#include "../Core/StructDef.cpp"
#include "../Core/BasicMath/BasicMath_Constants.cpp"
#include "../Core/TimeStamp.cpp"

#include "../Core/Core.h"
#include "../IO/IO.h"
#include "../Tight_Binding/Tight_Binding.h"

#include "../DMFT/SingleSiteFunctions.h"
#include "../Hubbard_Model/CellFunctions.h"
#include "../Ligand_Field_Theory/OperatorDeffinitions.h"
#include "../Model/StandardOperators.h"
#include "../Tools/TemperatureAvverage.h"

#include "../LuaInterface/lua/includeLua.h"
#include "../LuaInterface/BraKet.h"
#include "../LuaInterface/ComplexHeader.h"
#include "../LuaInterface/Matrix.h"
#include "../LuaInterface/Operator.h"
#include "../LuaInterface/Plot.h"
#include "../LuaInterface/ResponseFunction.h"
#include "../LuaInterface/Spectra.h"
#include "../LuaInterface/TightBinding.h"
#include "../LuaInterface/Wavefunction.h"
#include "../LuaInterface/eval.h"
#include "../LuaInterface/lua_extensions/TableToArray.h"
#include "../LuaInterface/lua_extensions/get_put.h"
#include "../LuaInterface/lua_extensions/tools.h"
// #include "../LuaInterface/main.h"
#include "../LuaInterface/ssp.h"
#include "../LuaInterface/units.h"


#include "../Martin/SlaterIntegrals.h"
#include "../Simon/functions_simon.h"
#include "../Yi/functions_yi.h"

#include "blas.h"
#include "lapack.h"
