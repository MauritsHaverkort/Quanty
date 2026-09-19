#include "BasicMath_BlockAndersonMatrix.h"
#include "BasicMath_BlockListOfPoles.h"
#include "BasicMath_BlockMatrixToOperator.h"
#include "BasicMath_BlockNaturalImpurityOrbital.h"
#include "BasicMath_BlockTridiagonalMatrix.h"
#include "BasicMath_BlockDoubleTridiagnonalMatrix.h"
#include "BasicMath_BlockTypeConversion.h"

uint32_t FreeBlockMatrix(ResponsefunctionType *M);
uint32_t FreeBlockMatrixFromLua(ResponsefunctionType *M);
int ConvertBlockMatrix(const ResponsefunctionType *const Min, ResponsefunctionType *Mout, const char type, const char * const Method, const double T, const double SingularValue, const double MinimalPoleWeigth, const double MinimalPoleDistance, bool deflate, bool relativeDistance, bool relativeWeigth);
int ConvertBlockMatrix(const ResponsefunctionType *const Min, ResponsefunctionType *Mout, const char type, const char * const Method, const double T, const double SingularValue, const double MinimalPoleWeigth, const double MinimalPoleDistance, bool deflate, bool relativeDistance, bool relativeWeigth, uint32_t const NTriMaxVal, uint32_t const NTriMaxCon);
int BlockResponseFunctionRotate(const ResponsefunctionType *const G, const CompactMatrixType *const U, ResponsefunctionType * const GRot, double const T, double const SingularValue, double const MinimalPoleWeigth, double const MinimalPoleDistance);
