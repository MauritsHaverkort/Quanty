#include "BasicMath_BlockAndersonMatrix.h"
#include "BasicMath_BlockListOfPoles.h"
#include "BasicMath_BlockMatrixToOperator.h"
#include "BasicMath_BlockNaturalImpurityOrbital.h"
#include "BasicMath_BlockTridiagonalMatrix.h"
#include "BasicMath_BlockTypeConversion.h"

uint32_t FreeBlockMatrix(ResponsefunctionType *M);
uint32_t FreeBlockMatrixFromLua(ResponsefunctionType *M);
int ConvertBlockMatrix(const ResponsefunctionType *const Min, ResponsefunctionType *Mout, const char type, double T, double SingularValue);
int BlockResponseFunctionRotate(const ResponsefunctionType *const G, const CompactMatrixType *const U, ResponsefunctionType * const GRot);
