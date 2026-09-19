uint32_t SumGreensfunctionsTriDiagonalRepresentationNonOrtho(TriDiagonalMatrixType const A, TriDiagonalMatrixType const B, TriDiagonalMatrixType * const M);
uint32_t SumGreensfunctionsTriDiagonalRepresentation(TriDiagonalMatrixType const A, TriDiagonalMatrixType const B, TriDiagonalMatrixType * const M);
uint32_t SymmetrizeAndersonMatrix(AndersonMatrixType * const M);
uint32_t SymmetrizeBADoubleTriDiagonalMatrix(BADoubleTriDiagonalMatrixType * const M);
uint32_t SymmetrizeTridiagonalMatrix(TriDiagonalMatrixType * const M);
void ReduceSqrtListOfPolesNumberFermi(ListOfPoles * M, uint32_t N, double const Ef);
void ReduceSqrtListOfPolesNumber(ListOfPoles * M, uint32_t N);
void ReduceSqrtListOfPolesMoment(ListOfPoles * M, uint32_t N);
uint32_t ReducePolesAndersonMatrixNumberFermi(AndersonMatrixType * const A, uint32_t N, double const Ef);
uint32_t ReducePolesAndersonMatrixNumber(AndersonMatrixType * const A, uint32_t N);
uint32_t ReducePolesAndersonMatrixMoment(AndersonMatrixType * const A, uint32_t N);
uint32_t ReducePolesAndersonMatrix(AndersonMatrixType * const A, double tMin, double Ef);
uint32_t ReducePolesTriDiagonalMatrix(TriDiagonalMatrixType * const A, uint32_t const NpolesMax, double const Ef);
uint32_t DifferenceGreensFunctionTriDiagonalRepresentation(TriDiagonalMatrixType const A, TriDiagonalMatrixType const B,TriDiagonalMatrixType * const M);
uint32_t SelfEnergyTriDiagonalRepresentation(TriDiagonalMatrixType const G0, TriDiagonalMatrixType const G,TriDiagonalMatrixType * const Sigma);
uint32_t SelfEnergyListOfPolesFromAndersonRepresentation(AndersonMatrixType const G0, AndersonMatrixType const G,ListOfPoles * const Sigma);
uint32_t SelfEnergyListOfPolesFromTriDiagonalRepresentation(TriDiagonalMatrixType const G0,TriDiagonalMatrixType const G, ListOfPoles * const Sigma);

int   CalculateGBathFromGandSigma(ResponsefunctionType  * const G, ResponsefunctionType  * const Sigma, ResponsefunctionType  * const GBath, double Epsilon, CompactVectorType *const Grid, uint32_t Nv, uint32_t Nc, ResponsefunctionType  * const Gv, ResponsefunctionType  * const Gc, char const * const Method, char const * const Type);


uint32_t GBathAndTriAndRepresentation(AndersonMatrixType const G0, double const mu, TriDiagonalMatrixType const Sigma,AndersonMatrixType * const G);
uint32_t GBathTriDiagonalRepresentation(TriDiagonalMatrixType const G0, TriDiagonalMatrixType const Sigma,TriDiagonalMatrixType * const GBath, uint32_t NTriDiag);
uint32_t GBathAndListAndRepresentation(AndersonMatrixType const G0, double const mu, ListOfPoles const Sigma,AndersonMatrixType * const G);
uint32_t GBathAndListAndRepresentationOnGrid(AndersonMatrixType const G0, double const mu, ListOfPoles const Sigma, double *Grid, uint32_t NGridPoints, AndersonMatrixType * const G);
uint32_t BlockGBathAndListAndRepresentation(BlockAndersonMatrixType const * const G0, double const mu, BlockListOfPoles const * const Sigma, BlockAndersonMatrixType * const G, double Epsilon);
uint32_t BlockGBathAndListAndRepresentationOnGrid(BlockAndersonMatrixType const * const G0, double const mu, BlockListOfPoles const * const Sigma, double const * const Grid, uint32_t NGridPoints, BlockAndersonMatrixType * const G);
uint32_t GBathAndTriAndRepresentationParallel(AndersonMatrixType const G0, TriDiagonalMatrixType const Sigma,AndersonMatrixType * const G);
uint32_t OldGBathAndTriAndRepresentation(AndersonMatrixType const G0, TriDiagonalMatrixType const Sigma,AndersonMatrixType * const G);
uint32_t OldGBathNewTriDiagonalRepresentation(TriDiagonalMatrixType const Gcl, TriDiagonalMatrixType const Sigma,TriDiagonalMatrixType * const GBathNew);
uint32_t VeryOldGBathNewTriDiagonalRepresentation(TriDiagonalMatrixType const Glocal, TriDiagonalMatrixType const Sigma,TriDiagonalMatrixType * const GBathNew);
int GListListListRepresentation(ListOfPoles * G0, ListOfPoles * Sigma, ListOfPoles *G);
int GListListListRepresentationOnGrid(ListOfPoles * G0, ListOfPoles * Sigma, double *Grid, uint32_t NGridPoints,  ListOfPoles *G);
int BlockGListListListRepresentation(BlockListOfPoles * G0, BlockListOfPoles * Sigma, BlockListOfPoles *G);
uint32_t TriDiagonalMatrixDysonScalar(TriDiagonalMatrixType const G, TriDiagonalMatrixType * const Gnew,MatrixType const Rho, uint32_t const NRho, TriDiagonalMatrixType * const Gp);
int BlockGListListListRepresentationOnGrid(BlockListOfPoles * G0, BlockListOfPoles * Sigma, double *Grid, uint32_t NGridPoints,  BlockListOfPoles *G);
uint32_t TridiagonalizeImpurityCoupledToBathDyson(TriDiagonalMatrixType const O, MatrixType const *V,MatrixType const H, TriDiagonalMatrixType const *g, uint32_t N,TriDiagonalMatrixType * const gnew);
uint32_t TridiagonalizeImpurityCoupledToBathDysonSecondOrderA(TriDiagonalMatrixType const O, MatrixType const *V1,MatrixType const *V2, MatrixType const H1,MatrixType const H2, TriDiagonalMatrixType const *g,uint32_t N, TriDiagonalMatrixType * const gnew);
uint32_t TridDiagonalMatrixDysonToMuchSimplified(TriDiagonalMatrixType const G1, MatrixType const *M12,TriDiagonalMatrixType const *g2, MatrixType const *Emu,int32_t const Ni, TriDiagonalMatrixType * const Gnew);
uint32_t GreenFunctionConvoluteTriDiagonalRepresentation(TriDiagonalMatrixType * const fg,TriDiagonalMatrixType const f, TriDiagonalMatrixType const g,uint32_t const NTriDiag);
Complex16 EvaluateTriDiagonalMatrix(TriDiagonalMatrixType const * A, const double omega, const double Gamma);
Complex16 EvaluateTriDiagonalMatrix(void const * A, const double omega, const double Gamma);
Complex16 EvaluateListOfPoles(ListOfPoles const * A, const double omega, const double Gamma);
Complex16 EvaluateListOfPoles(void const * A, const double omega, const double Gamma);
Complex16 EvaluateAndersonMatrix(AndersonMatrixType const * A, const double omega, const double Gamma);
Complex16 EvaluateAndersonMatrix(void const * A, const double omega, const double Gamma);
Complex16 EvaluateBADoubleTriDiagonalMatrix(BADoubleTriDiagonalMatrixType const * A, const double omega, const double Gamma);
Complex16 EvaluateBADoubleTriDiagonalMatrix(void const * A, const double omega, const double Gamma);
Complex16 EvaluateResponsefunctionType(ResponsefunctionType const *A, const double omega, const double Gamma);
/*  NTriMaxNoLimit asks a conversion not to cut the chain it produces. Zero is a chain of
 *  a[0] alone and so a value a caller can legitimately ask for, which is why the "not set"
 *  marker has to sit at the top of the range rather than at the bottom.
 */
#define NTriMaxNoLimit UINT32_MAX
uint32_t NTriMaxOr(uint32_t const NTriMax, uint32_t const N);
int ConvertResponsefunction(const ResponsefunctionType *const Min, ResponsefunctionType *Mout, const char type, const char * const Method, const double T, const double SingularValue, const double MinimalPoleWeigth, const double MinimalPoleDistance, bool deflate, bool relativeDistance, bool relativeWeigth);
int ConvertResponsefunction(const ResponsefunctionType *const Min, ResponsefunctionType *Mout, const char type, const char * const Method, const double T, const double SingularValue, const double MinimalPoleWeigth, const double MinimalPoleDistance, bool deflate, bool relativeDistance, bool relativeWeigth, uint32_t const NTriMaxVal, uint32_t const NTriMaxCon);
int ResponsefunctionChop(ResponsefunctionType *M, double epsilon, bool Deflate);
int ResponsefunctionTruncate(ResponsefunctionType *M, uint32_t const NTriMaxVal, uint32_t const NTriMaxCon);
int ResponseFunctionInvertEnergy(ResponsefunctionType * const Gin, ResponsefunctionType * const Gout);
int ResponseFunctionTranspose(ResponsefunctionType * const Gin, ResponsefunctionType * const Gout);
int ResponseFunctionDownfold(ResponsefunctionType * const Gin, CompactMatrixType const * const R, ResponsefunctionType * const Gout, CompactMatrixType * const S, double const T, double const SingularValue, double const MinimalPoleWeigth, double const MinimalPoleDistance);
int ResponsefunctionCopy(ResponsefunctionType *Min, ResponsefunctionType *Mout);
uint32_t FreeResponsefunction(ResponsefunctionType *M);
uint32_t FreeResponsefunctionIncludingOriginalPointer(ResponsefunctionType *M);
