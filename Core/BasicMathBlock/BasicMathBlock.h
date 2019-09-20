#include "BasicMath_BlockAndersonMatrix.h"
#include "BasicMath_BlockListOfPoles.h"
#include "BasicMath_BlockMatrixToOperator.h"
#include "BasicMath_BlockNaturalImpurityOrbital.h"
#include "BasicMath_BlockTridiagonalMatrix.h"
#include "BasicMath_BlockTypeConversion.h"

uint32_t FreeBlockMatrix(ResponsefunctionType *M);
uint32_t FreeBlockMatrixFromLua(ResponsefunctionType *M);
ResponsefunctionType *ConvertBlockMatrix(ResponsefunctionType *M, const char type, double T, double SingularValue);
