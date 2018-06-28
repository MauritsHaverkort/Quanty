
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

#include "../Core/StructDef.c"

#include "../Core/BasicMath/BasicMath_Constants.c"
#include "../Core/BasicMath.h"
#include "../Core/BasicMath/BasicMath_AndersonMatrixVector.h"
#include "../Core/BasicMath/BasicMath_BADoubleTriDiagonalMatrixVector.h"
#include "../Core/BasicMath/BasicMath_BlokAndersonMatrixVector.h"
#include "../Core/BasicMath/BasicMath_BlokDoubleTriDiagonalMatrixVector.h"
#include "../Core/BasicMath/BasicMath_BlokTriDiagonalMatrixVector.h"
#include "../Core/BasicMath/BasicMath_Combinatorial.h"
#include "../Core/BasicMath/BasicMath_CompactMatrixVector.h"
#include "../Core/BasicMath/BasicMath_Complex.h"
#include "../Core/BasicMath/BasicMath_DoubleTriDiagonalMatrixVector.h"
#include "../Core/BasicMath/BasicMath_GreenFunctionArithmetic.h"
#include "../Core/BasicMath/BasicMath_ListOfPoles.h"
#include "../Core/BasicMath/BasicMath_MatrixTypeConversion.h"
#include "../Core/BasicMath/BasicMath_MatrixVector.h"
#include "../Core/BasicMath/BasicMath_SparseMatrixVector.h"
#include "../Core/BasicMath/BasicMath_StandardFunctions.h"
#include "../Core/BasicMath/BasicMath_TriDiagonalMatrixVector.h"
#include "../Core/BasicMath/BasicMath_VectorLinearAlgabra.h"

#include "../Core/BasicMathBlock/BasicMathBlock.h"
#include "../Core/BasicMathBlock/BasicMath_BlockAndersonMatrix.h"
#include "../Core/BasicMathBlock/BasicMath_BlockListOfPoles.h"
#include "../Core/BasicMathBlock/BasicMath_BlockMatrixToOperator.h"
#include "../Core/BasicMathBlock/BasicMath_BlockNaturalImpurityOrbital.h"
#include "../Core/BasicMathBlock/BasicMath_BlockTridiagonalMatrix.h"
#include "../Core/BasicMathBlock/BasicMath_BlockTypeConversion.h"

#include "../Core/BlockLanczos.h"
#include "../Core/BlockWaveFunctionArithmetic.h"
#include "../Core/ComplexRealWaveFunctionArithmetic.h"
#include "../Core/ComplexWaveFunctionArithmetic.h"
#include "../Core/DensityMatrix.h"
#include "../Core/Error.h"
#include "../Core/Hash.h"
#include "../Core/ImpurityLanczos.h"
#include "../Core/Lanczos.h"
#include "../Core/MatrixOperator.h"
#include "../Core/MeanField.h"
#include "../Core/OperatorArithmetic.h"
#include "../Core/OperatorWaveFunction.h"
#include "../Core/OperatorWaveFunctionPerturbation.h"
#include "../Core/Random.h"
#include "../Core/RealComplexWaveFunctionArithmetic.h"
#include "../Core/RealWaveFunctionArithmetic.h"
#include "../Core/ResonantSpectra.h"
#include "../Core/Spectral.h"
#include "../Core/TimeStamp.h"
#include "../Core/WaveFunctionArithmetic.h"
#include "../Core/WaveFunctionPerturbationArithmetic.h"

#include "../DMFT/SingleSiteFunctions.h"

#include "../Hubbard_Model/CellFunctions.h"

#include "../IO/BitMap.h"
#include "../IO/DresdenFPLO.h"
#include "../IO/FilePrint.h"
#include "../IO/FileRead.h"
#include "../IO/Print.h"
#include "../IO/PrintMathematicaExpression.h"
#include "../IO/Read.h"
#include "../IO/StreamPrint.h"
#include "../IO/StreamRead.h"
#include "../IO/StringPrint.h"
#include "../IO/StuttgartNMTO.h"

#include "../Ligand_Field_Theory/OperatorDeffinitions.h"

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

#include "../Model/AndersonImpurity.h"
#include "../Model/CrystalField.h"
#include "../Model/LigandField.h"
#include "../Model/OneDBand.h"
#include "../Model/StandardOperators.h"

#include "../Tight_Binding/DresdenFPLO.h"
#include "../Tight_Binding/ImpurityUnfolding.h"
#include "../Tight_Binding/StuttgartNMTO.h"
#include "../Tight_Binding/TightBindingToG.h"
#include "../Tight_Binding/Tight_Binding.h"
#include "../Tight_Binding/Tight_Binding_lib.h"
#include "../Tools/TemperatureAvverage.h"

#include "blas.h"
#include "lapack.h"
