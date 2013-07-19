/**********************************************************/
/*****************       blas/lapack      *****************/
/**********************************************************/
int dstev();
int dstevd();
int dsyev();
int dsyevd();
int zheev();
int zheevd();
int dgesvd();
int zgesvd();
int dsygv();
int zhegv();
int dsygvd();
int zhegvd();


/**********************************************************/
/*****************       IO/Print.c       *****************/
/**********************************************************/

uint32_t PrintSpectra();
uint32_t PrintWaveFunctionReadable();
uint32_t PrintWaveFunctionFullReadable();
uint32_t PrintOperatorReadable();
uint32_t PrintOperatorFullReadable();
uint32_t PrintMatrix();
uint32_t PrintSparseMatrix();
uint32_t PrintTriDiagonalMatrix();
uint32_t PrintBlokTriDiagonalMatrix();
uint32_t PrintDoubleTriDiagonalMatrix();
uint32_t PrintBlokDoubleTriDiagonalMatrix();
uint32_t PrintAndersonMatrix();
uint32_t PrintBlokAndersonMatrix();
uint32_t PrintOperatorFullReadable(OperatorType O);

uint32_t FilePrintSpectra();
uint32_t FilePrintWaveFunctionReadable();
uint32_t FilePrintWaveFunctionFullReadable();
uint32_t FilePrintOperatorReadable();
uint32_t FilePrintOperatorFullReadable();
uint32_t FilePrintMatrix();
uint32_t FilePrintSparseMatrix();
uint32_t FilePrintPlainTriDiagonalMatrix();
uint32_t FilePrintBlokTriDiagonalMatrix();
uint32_t FilePrintDoubleTriDiagonalMatrix();
uint32_t FilePrintBlokDoubleTriDiagonalMatrix();
uint32_t FilePrintAndersonMatrix();
uint32_t FilePrintBlokAndersonMatrix();
uint32_t FilePrintOperatorFullReadable(OperatorType O, const char *const flname, const char *const mode);

uint32_t StreamPrintSpectra();
uint32_t StreamPrintWaveFunctionReadable();
uint32_t StreamPrintWaveFunctionFullReadable();
uint32_t StreamPrintOperatorReadable();
uint32_t StreamPrintOperatorFullReadable();
uint32_t StreamPrintMatrix();
uint32_t StreamPrintSparseMatrix();
uint32_t StreamPrintTriDiagonalMatrix();
uint32_t StreamPrintBlokTriDiagonalMatrix();
uint32_t StreamPrintDoubleTriDiagonalMatrix();
uint32_t StreamPrintBlokDoubleTriDiagonalMatrix();
uint32_t StreamPrintAndersonMatrix();
uint32_t StreamPrintBlokAndersonMatrix();
uint32_t StreamPrintOperatorFullReadable(OperatorType O, FILE *stream);

uint32_t PrintSparseMatrixMathematica();
uint32_t StreamPrintSparseMatrixMathematica();
uint32_t FilePrintSparseMatrixMathematica();

/**********************************************************/
/**********       BasisMath/MatrixVector.c       **********/
/**********************************************************/
void     FreeMatrix();
uint32_t InitMatrix();
uint32_t InitIdentityMatrix();
uint32_t MakeRealMatrixComplex();
void     MatrixChop();
uint32_t MatrixTranspose();
void     MatrixComplexConjugate();
uint32_t MatrixConjugateTranspose();
double   MatrixInverse();
uint32_t MatrixAdd();
uint32_t MatrixScalarMultiply();
uint32_t MatrixTransform();
uint32_t MatrixTransformPart();
uint32_t MatrixMultiply();
uint32_t MatrixPart();
uint32_t MatrixCopy();
uint32_t RealMatrixVector();
uint32_t RealTransposeMatrixVector();
uint32_t ComplexMatrixVector();
uint32_t ComplexConjugateTransposeMatrixVector();
uint32_t RealVectorMatrixVector();
uint32_t ComplexVectorMatrixVector();
void     MatrixOrthonormalize();

/**********************************************************/
/*******       BasisMath/SparseMatrixVector.c       *******/
/**********************************************************/
void     FreeSparseMatrix();
uint32_t InitSparseMatrix();
uint32_t InitSparseIdentityMatrix();
uint32_t MakeRealSparseMatrixComplex();
uint32_t SparseMatrixAddElement();
void     SparseMatrixExtractElement();
void     MinimalAlocSpareseMatrix();
void     RemoveZeroSpareseMatrix();
void     SparseMatrixChop();
void     SparseMatrixTranspose();
void     SparseMatrixComplexConjugate();
void     SparseMatrixConjugateTranspose();
double   SparseMatrixInverse();
uint32_t SparseMatrixAdd();
uint32_t SparseMatrixScalarMultiply();
uint32_t SparseMatrixTransform();
uint32_t SparseMatrixDenseTransform();
uint32_t DenseMatrixSparseTransform();
uint32_t SparseMatrixTransformPart();
uint32_t SparseMatrixDenseTransformPart();
uint32_t DenseMatrixSparseTransformPart();
uint32_t SparseMatrixMultiply();
uint32_t SparseMatrixDenseMultiply();
uint32_t DenseMatrixSparseMultiply();
uint32_t SparseMatrixPart();
uint32_t SparseMatrixCopy();
uint32_t RealSparseMatrixVector();
uint32_t RealTransposeSparseMatrixVector();
void     ComplexSparseMatrixVector();
void     ComplexConjugateTransposeSparseMatrixVector();
uint32_t RealVectorSparseMatrixVector();
void     ComplexVectorSparseMatrixVector();
void     SparseMatrixOrthonormalize();


/**********************************************************/
/****       BasisMath/TriDiagonalMatrixVector.c        ****/
/**********************************************************/
void     FreeTriDiagonalMatrix();
uint32_t InitTriDiagonalMatrix();
uint32_t InitIdentityTriDiagonalMatrix();
void     TriDiagonalMatrixChop();
double   TriDiagonalMatrixInverse();
void     TriDiagonalMatrixAdd();
void     TriDiagonalMatrixScalarMultiply();
uint32_t TriDiagonalMatrixPart();
uint32_t TriDiagonalMatrixCopy();
void     RealTriDiagonalMatrixVector();
void     ComplexTriDiagonalMatrixVector();
void     RealVectorTriDiagonalMatrixVector();
void     ComplexVectorTriDiagonalMatrixVector();


/**********************************************************/
/**      BasisMath/DoubleTriDiagonalMatrixVector.c       **/
/**********************************************************/
void     FreeDoubleTriDiagonalMatrix();
uint32_t InitDoubleTriDiagonalMatrix();
uint32_t InitIdentityDoubleTriDiagonalMatrix();
void     DoubleTriDiagonalMatrixChop();
double   DoubleTriDiagonalMatrixInverse();
void     DoubleTriDiagonalMatrixAdd();
void     DoubleTriDiagonalMatrixScalarMultiply();
uint32_t DoubleTriDiagonalMatrixPart();
uint32_t DoubleTriDiagonalMatrixCopy();
uint32_t RealDoubleTriDiagonalMatrixVector();
uint32_t ComplexDoubleTriDiagonalMatrixVector();
uint32_t RealVectorDoubleTriDiagonalMatrixVector();
uint32_t ComplexVectorDoubleTriDiagonalMatrixVector();


/**********************************************************/
/******       BasisMath/AndersonMatrixVector.c       ******/
/**********************************************************/
void     FreeAndersonMatrix();
uint32_t InitAndersonMatrix();
uint32_t InitIdentityAndersonMatrix();
void     AndersonMatrixChop();
double   AndersonMatrixInverse();
void     AndersonMatrixAdd();
void     AndersonMatrixScalarMultiply();
uint32_t AndersonmatrixPart();
uint32_t AndersonmatrixCopy();
void     RealAndersonMatrixVector();
void     ComplexAndersonMatrixVector();
void     RealVectorAndersonMatrixVector();
void     ComplexVectorAndersonMatrixVector();


/**********************************************************/
/**       BasisMath/BlokTriDiagonalMatrixVector.c        **/
/**********************************************************/
void     FreeBlokTriDiagonalMatrix();
void     FreeFrameBlokTriDiagonalMatrix();
uint32_t InitBlokTriDiagonalMatrix();
uint32_t InitFrameBlokTriDiagonalMatrix();
uint32_t InitIdentityBlokTriDiagonalMatrix();
void     BlokTriDiagonalMatrixChop();
void     BlokTriDiagonalMatrixAdd();
void     BlokTriDiagonalMatrixScalarMultiply();
uint32_t BlokTriDiagonalMatrixPart();
uint32_t BlokTriDiagonalMatrixCopy();
void     BlokTriDiagonalMatrixTranspose();
void     BlokTriDiagonalMatrixComplexConjugate();
void     BlokTriDiagonalMatrixConjugateTranspose();
uint32_t RealBlokTriDiagonalMatrixVector();
uint32_t ComplexBlokTriDiagonalMatrixVector();
uint32_t RealVectorBlokTriDiagonalMatrixVector();
uint32_t ComplexVectorBlokTriDiagonalMatrixVector();


/**********************************************************/
/****      BasisMath/BlokAndersonMatrixVector.c        ****/
/**********************************************************/
void     FreeBlokAndersonMatrix();
void     FreeFrameBlokAndersonMatrix();
uint32_t InitBlokAndersonMatrix();
uint32_t InitFrameBlokAndersonMatrix();
uint32_t InitIdentityBlokAndersonMatrix();
void     BlokAndersonMatrixChop();
void     BlokAndersonMatrixAdd();
void     BlokAndersonMatrixScalarMultiply();
uint32_t BlokAndersonMatrixPart();
uint32_t BlokAndersonMatrixCopy();
void     BlokAndersonMatrixTranspose();
void     BlokAndersonMatrixComplexConjugate();
void     BlokAndersonMatrixConjugateTranspose();
uint32_t RealBlokAndersonMatrixVector();
uint32_t ComplexBlokAndersonMatrixVector();
uint32_t RealVectorBlokAndersonMatrixVector();
uint32_t ComplexVectorBlokAndersonMatrixVector();


/**********************************************************/
/**    BasisMath/BlokDoubleTriDiagonalMatrixVector.c     **/
/**********************************************************/
void     FreeBlokDoubleTriDiagonalMatrix();
void     FreeFrameBlokDoubleTriDiagonalMatrix();
uint32_t InitBlokDoubleTriDiagonalMatrix();
uint32_t InitFrameBlokDoubleTriDiagonalMatrix();
uint32_t InitIdentityBlokDoubleTriDiagonalMatrix();
uint32_t BlokDoubleTriDiagonalMatrixChop();
uint32_t BlokDoubleTriDiagonalMatrixAdd();
uint32_t BlokDoubleTriDiagonalMatrixScalarMultiply();
uint32_t BlokDoubleTriDiagonalMatrixPart();
uint32_t BlokDoubleTriDiagonalMatrixCopy();
uint32_t BlokDoubleTriDiagonalMatrixTranspose();
uint32_t BlokDoubleTriDiagonalMatrixComplexConjugate();
uint32_t BlokDoubleTriDiagonalMatrixConjugateTranspose();
uint32_t RealBlokDoubleTriDiagonalMatrixVector();
uint32_t ComplexBlokDoubleTriDiagonalMatrixVector();
uint32_t RealVectorBlokDoubleTriDiagonalMatrixVector();
uint32_t ComplexVectorBlokDoubleTriDiagonalMatrixVector();



/**********************************************************/
/***********     BasisMath/Combinatorial.c      ***********/
/**********************************************************/
double Factorial(int a);
uint32_t BinomialCoefficient(uint32_t n, uint32_t k);
uint32_t PerturbationOptions(uint32_t NV, uint32_t order);
void PerturbationOrderListToIndex(uint32_t *const POList, uint32_t *const POIndex, uint32_t const NV, uint32_t const order);
void PerturbationOrderIndexToList(uint32_t const POIndex, uint32_t *const POList, uint32_t const NV, uint32_t const order);
double ThreeJSymbol(int j1, int m1, int j2, int m2, int j3, int m3);
double ThreeJSymbolHalf(int j1, int m1, int j2, int m2, int j3, int m3);
double SlaterCoefficientC(int k, int l1, int m1, int l2, int m2);


/**********************************************************/
/********    BasisMath/MatrixTypeConversion.c      ********/
/**********************************************************/
/********    AndersonMatrixTo                      ********/
uint32_t AndersonMatrixToSparseMatrix();
uint32_t AndersonMatrixToTriDiagonalMatrix();
uint32_t AndersonMatrixToEigensystem();
uint32_t MatrixToDoubleTriDiagonalMatrix();
/********    MatrixTo                              ********/
uint32_t BlokAndersonMatrixToSparseMatrix();
uint32_t MatrixToEigensystem();
uint32_t BlokAndersonMatrixToBlokTriDiagonalMatrix();
uint32_t BlokAndersonMatrixToBlokDoubleTriDiagonalMatrix();
/********    TriDiagonalMatrixTo                   ********/
uint32_t TriDiagonalMatrixToSparseMatrix();
uint32_t TriDiagonalMatrixToEigensystem();
uint32_t TriDiagonalMatrixToAndersonMatrix();
uint32_t TriDiagonalMatrixToDoubleTriDiagonalMatrix();
/********    BlokTriDiagonalMatrixTo               ********/
uint32_t BlokTriDiagonalMatrixToSparseMatrix();
uint32_t BlokTriDiagonalMatrixToEigensystem();
uint32_t BlokTriDiagonalMatrixToBlokAndersonMatrix();
uint32_t BlokTriDiagonalMatrixToBlokDoubleTriDiagonalMatrix();
/********    DoubleTriDiagonalMatrixTo             ********/
uint32_t DoubleTriDiagonalMatrixToSparseMatrix();
uint32_t DoubleTriDiagonalMatrixToTriDiagonalMatrix();
uint32_t DoubleTriDiagonalMatrixToEigensystem();
uint32_t DoubleTriDiagonalMatrixToAndersonMatrix();
/********    BlokDoubleTriDiagonalMatrixTo         ********/
uint32_t BlokDoubleTriDiagonalMatrixToSparseMatrix();
uint32_t BlokDoubleTriDiagonalMatrixToEigensystem();
uint32_t BlokDoubleTriDiagonalMatrixToBlokTriDiagonalMatrix();
uint32_t BlokDoubleTriDiagonalMatrixToBlokAndersonMatrix();
/********    SparseMatrixTo                        ********/
uint32_t SparseMatrixToMatrix();
uint32_t SparseMatrixToTriDiagonalMatrix();
uint32_t SparseMatrixToEigensystem();
uint32_t SparseMatrixToAndersonMatrix();
uint32_t SparseMatrixToBlokTriDiagonalMatrix();
uint32_t SparseMatrixToBlokAndersonMatrix();
uint32_t SparseMatrixToDoubleTriDiagonalMatrix();
uint32_t SparseMatrixToBlokDoubleTriDiagonalMatrix();
/********    MatrixTo                              ********/
uint32_t MatrixToSparseMatrix();
uint32_t MatrixToTriDiagonalMatrix();          uint32_t TriDiagonalize();
uint32_t MatrixToEigensystem();                uint32_t Eigensystem();
uint32_t MatrixToGeneralizedEigensystem();     uint32_t GeneralizedEigensystem();
uint32_t MatrixToAndersonMatrix();
uint32_t MatrixToBlokTriDiagonalMatrix();      uint32_t BlokTriDiagonalize();
uint32_t MatrixToBlokAndersonMatrix();
uint32_t MatrixToDoubleTriDiagonalMatrix();
uint32_t MatrixToBlokDoubleTriDiagonalMatrix();
/********    ListOfPolesTo                              ********/
uint32_t ListOfPolesToBADoubleTridiagonalMatrix();


/**********************************************************/
/********    BasisMath/MatrixTypeConversion.c      ********/
/**********************************************************/

uint32_t TriDiagonalMatrixToEigensystem();
uint32_t MatrixToSparseMatrix();
uint32_t TriDiagonalize();
uint32_t MatrixToTriDiagonalMatrix();
uint32_t Eigensystem();
uint32_t MatrixToEigensystem();
uint32_t MatrixToAndersonMatrix();
uint32_t BlokTriDiagonalize();
uint32_t MatrixToBlokTriDiagonalMatrix();
uint32_t MatrixToBlokAndersonMatrix();
uint32_t MatrixToDoubleTriDiagonalMatrixNoOverlap();
uint32_t MatrixToBlokDoubleTriDiagonalMatrixNoOverlap();






/**********************************************************/
/********     Core/WaveFunctionArithmetic.c        ********/
/**********************************************************/

uint32_t WaveFunctionInit(WaveFunctionType *Psi);
void WaveFunctionFree(WaveFunctionType *const psi);
uint32_t PsiHermitianOperatorPsi(const OperatorType O, const WaveFunctionType psi, double *VRe);
uint32_t PsiHermitianOperatorSqrPsi(const OperatorType O, const WaveFunctionType psi, double *VRe);

/**********************************************************/
/*****  Core/WaveFunctionPerturbationArithmetic.c     *****/
/**********************************************************/

uint32_t WaveFunctionPerturbationInit(WaveFunctionPerturbationType *Psi);
uint32_t WaveFunctionPerturbationInitAll(WaveFunctionPerturbationType *Psi, uint32_t NFermion, uint32_t NBoson, uint32_t QComplex, uint32_t NBitsKey);
uint32_t QWaveFunctionPerturbationMinimumBloksAlocated(WaveFunctionPerturbationType psi);
uint32_t MakeComplexWaveFunctionPerturbationFromReal(WaveFunctionPerturbationType *Cpsi, WaveFunctionPerturbationType psi);
uint32_t MakeRealWaveFunctionPerturbationFromComplex(WaveFunctionPerturbationType Cpsi, WaveFunctionPerturbationType *psi);
void WaveFunctionPerturbationFreeFramework(WaveFunctionPerturbationType *const psi);
void WaveFunctionPerturbationFree(WaveFunctionPerturbationType *const psi);
void WaveFunctionPerturbationFreeValues(WaveFunctionPerturbationType *const psi);
uint32_t WaveFunctionPerturbationAddElement(WaveFunctionPerturbationType *const psi, const uint8_t *const det, const uint32_t order, const uint32_t iV, const double aRe, const double aIm);
uint32_t WaveFunctionPerturbationAddElementConserveBasis(WaveFunctionPerturbationType *const psi, const uint8_t *const det, const uint32_t order, const uint32_t iV, const double aRe, const double aIm);
void WaveFunctionPerturbationPrefactorDet(const WaveFunctionPerturbationType psi, const uint8_t *const det, const uint32_t order, const uint32_t iV, double *aRe, double *aIm);
uint32_t WaveFunctionPerturbationAddSameBasis(WaveFunctionPerturbationType * const a, const WaveFunctionPerturbationType b);
uint32_t WaveFunctionPerturbationAddMultipliedSameBasis(WaveFunctionPerturbationType * const a, const WaveFunctionPerturbationType b, const double alpha);
uint32_t WaveFunctionPerturbationAddComplexMultipliedSameBasis(WaveFunctionPerturbationType * const a, const WaveFunctionPerturbationType b, const double alphaRe, const double alphaIm);
uint32_t WaveFunctionPerturbationAddConserveBasis(WaveFunctionPerturbationType * const a, const WaveFunctionPerturbationType b);
uint32_t WaveFunctionPerturbationAddMultipliedConserveBasis(WaveFunctionPerturbationType * const a, const WaveFunctionPerturbationType b, const double alpha);
uint32_t WaveFunctionPerturbationAddComplexMultipliedConserveBasis(WaveFunctionPerturbationType * const a, const WaveFunctionPerturbationType b, const double alphaRe, const double alphaIm);
uint32_t WaveFunctionPerturbationAddUnionBasis(WaveFunctionPerturbationType *const a, const WaveFunctionPerturbationType b);
uint32_t WaveFunctionPerturbationAddMultipliedUnionBasis(WaveFunctionPerturbationType *const a, const WaveFunctionPerturbationType b, const double alpha);
uint32_t WaveFunctionPerturbationAddComplexMultipliedUnionBasis(WaveFunctionPerturbationType *const a, const WaveFunctionPerturbationType b, const double alphaRe, const double alphaIm);
uint32_t WaveFunctionPerturbationMultiply(WaveFunctionPerturbationType a, const double bRe, const double bIm);
void WaveFunctionPerturbationCalculateNorm(const WaveFunctionPerturbationType a, double *Norm);
void WaveFunctionPerturbationNormalize(WaveFunctionPerturbationType *const a);
void WaveFunctionPerturbationDotProductSameBasis(const WaveFunctionPerturbationType a, const WaveFunctionPerturbationType b, double *DotRe, double *DotIm);
void WaveFunctionPerturbationDotProduct(const WaveFunctionPerturbationType a, const WaveFunctionPerturbationType b, double *DotRe, double *DotIm);
void WaveFunctionPerturbationRandomize(WaveFunctionPerturbationType *const a);
void WaveFunctionPerturbationTruncateBasis(WaveFunctionPerturbationType *const a, const double epsilon);
void WaveFunctionPerturbationOrthogonalizeConserveBasis(WaveFunctionPerturbationType *const psi, const WaveFunctionPerturbationType *const OtherPsi, const uint32_t Npsi);
uint32_t WaveFunctionPerturbationOrthogonalize(WaveFunctionPerturbationType *const psi, const WaveFunctionPerturbationType *const OtherPsi, const uint32_t Npsi);
uint32_t CreateStartWaveFunctionPerturbationOnsiteEnergy(const OperatorPerturbationType O, WaveFunctionPerturbationType *psi, uint32_t order, uint32_t nFermion, uint32_t nBoson, uint32_t NBitsKey);
uint32_t WaveFunctionPerturbationAddFermion(WaveFunctionPerturbationType const psiIn, WaveFunctionPerturbationType *const psiOut, uint32_t FermionPos, uint32_t occupied);
uint32_t WaveFunctionPerturbationRemoveFermion(WaveFunctionPerturbationType const psiIn, WaveFunctionPerturbationType *const psiOut, uint32_t FermionPos);


/**********************************************************/
/********   Core/RealWaveFunctionArithmetic.c      ********/
/**********************************************************/

void     RealWaveFunctionFree();
void     RealWaveFunctionFreeValues();
uint32_t RealWaveFunctionAddElement();
void     RealWaveFunctionAddElementConserveBasis();
double   RealWaveFunctionPrefactorDet();
uint32_t RealWaveFunctionAddElementOMPMiniFlush();
void     RealWaveFunctionAddSameBasis();
void     RealWaveFunctionAddMultipliedSameBasis();
uint32_t RealWaveFunctionAddUnionBasis();
uint32_t RealWaveFunctionAddMultipliedUnionBasis();
void     RealWaveFunctionAddMultipliedConserveBasis();
void     RealWaveFunctionNormalize();
void     RealWaveFunctionMultiply();
double   RealWaveFunctionCalculateNorm();
double   RealWaveFunctionDotProductSameBasis();
double   RealWaveFunctionDotProduct();
void     RealWaveFunctionRandomize();
void     RealWaveFunctionTruncateBasis();
void     RealWaveFunctionOrthogonalizeConserveBasis();
uint32_t RealWaveFunctionOrthogonalize();

/**********************************************************/
/********   Core/ComplexWaveFunctionArithmetic.c   ********/
/**********************************************************/

void     ComplexWaveFunctionFree();
void     ComplexWaveFunctionFreeValues();
uint32_t ComplexWaveFunctionAddElement();
void     ComplexWaveFunctionAddElementConserveBasis();
void     ComplexWaveFunctionPrefactorDet();
uint32_t ComplexWaveFunctionAddElementOMPMiniFlush();
void     ComplexWaveFunctionAddSameBasis();
void     ComplexWaveFunctionAddMultipliedSameBasis();
void     ComplexwaveFunctionAddComplexMultipliedSamebasis();
uint32_t ComplexWaveFunctionAddUnionBasis();
uint32_t ComplexWaveFunctionAddMultipliedUnionbasis();
uint32_t ComplexWaveFunctionAddComplexMultipliedUnionBasis();
void     ComplexWaveFunctionAddMultipliedSameBasis();
void     ComplexWaveFunctionNormalize();
void     ComplexWaveFunctionMultiply();
void     ComplexWaveFunctionComplexMultiply();
double   ComplexWaveFunctionCalculateNorm();
void     ComplexWaveFunctionDotProductSameBasis();
void     ComplexWaveFunctionDotProduct();
void     ComplexWaveFunctionRandomize();
void     ComplexWaveFunctionTruncatebasis();
uint32_t ComplexWaveFunctionOrthogonalize();
uint32_t MakeComplexWaveFunctionFromReal();
uint32_t MakeRealWaveFunctionFromComplex();



/**********************************************************/
/*************       Core/Lanczos.c          **************/
/**********************************************************/

uint32_t LanczosGroundStateRR(const OperatorType H, WaveFunctionType *const psi, double *epsilon, const double zero, double shiftin, WaveFunctionType *OtherPsi, uint32_t NOtherPsi);
uint32_t LanczosGroundStateRC(const OperatorType H, WaveFunctionType *const psi, double *epsilon, const double zero, double shiftin, WaveFunctionType *OtherPsi, uint32_t NOtherPsi);
uint32_t LanczosGroundStateCR(const OperatorType H, WaveFunctionType *const psi, double *epsilon, const double zero, double shiftin, WaveFunctionType *OtherPsi, uint32_t NOtherPsi);
uint32_t LanczosGroundStateCC(const OperatorType H, WaveFunctionType *const psi, double *epsilon, const double zero, double shiftin, WaveFunctionType *OtherPsi, uint32_t NOtherPsi);
uint32_t LanczosGroundStateMR(const OperatorType H, WaveFunctionType *const psi, double *epsilon, const double zero, double shiftin, WaveFunctionType *OtherPsi, uint32_t NOtherPsi);
uint32_t LanczosGroundStateMC(const OperatorType H, WaveFunctionType *const psi, double *epsilon, const double zero, double shiftin, WaveFunctionType *OtherPsi, uint32_t NOtherPsi);
uint32_t LanczosTriDiagonalizeRR(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType *const M, double Gamma);
uint32_t LanczosTriDiagonalizeCR(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType *const M, double Gamma);
uint32_t LanczosTriDiagonalizeMR(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType *const M, double Gamma);
uint32_t LanczosTriDiagonalizeRC(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType *const M, double Gamma);
uint32_t LanczosTriDiagonalizeCC(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType *const M, double Gamma);
uint32_t LanczosTriDiagonalizeMC(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType *const M, double Gamma);
uint32_t LanczosTriDiagonalizePlusDensityRR(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType *const M, MatrixType *const Rho, double Gamma);
uint32_t LanczosTriDiagonalizePlusDensityCR(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType *const M, MatrixType *const Rho, double Gamma);
uint32_t LanczosTriDiagonalizePlusDensityMR(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType *const M, MatrixType *const Rho, double Gamma);
uint32_t LanczosTriDiagonalizePlusDensityRC(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType *const M, MatrixType *const Rho, double Gamma);
uint32_t LanczosTriDiagonalizePlusDensityCC(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType *const M, MatrixType *const Rho, double Gamma);
uint32_t LanczosTriDiagonalizePlusDensityMC(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType *const M, MatrixType *const Rho, double Gamma);


uint32_t LanczosGroundState(const OperatorType H, WaveFunctionType *const psi, double *epsilon, const double zero, double shiftin, WaveFunctionType *OtherPsi, uint32_t NOtherPsi);
uint32_t LanczosTriDiagonalize(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType *const M, double Gamma);
uint32_t LanczosTriDiagonalizePlusDensity(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType *const M, MatrixType *const Rho, double Gamma);



/**********************************************************/
/********       Core/OperatorArithmetic.c          ********/
/**********************************************************/

uint32_t FreeOperator();
uint32_t InitOperator();
uint32_t MakeRealOperatorFixedLengthComplex();
uint32_t InitOperatorFixedLength();
uint32_t InitOperatorCopySettings();
uint32_t RealOperatorAddLader();
uint32_t ComplexOperatorAddLader();
uint32_t RealOperatorAddLaderNormalOrder();
uint32_t ComplexOperatorAddLaderNormalOrder();
uint32_t OperatorAddLader(OperatorType *const O, const uint16_t *const oadd, const uint32_t length, const double VRe, const double VIm);
uint32_t OperatorAddLaderNormalOrder(OperatorType *const O, const uint16_t *const oadd, const uint32_t length, const double VRe, const double VIm, uint32_t (*QDetNotOrdered)(uint16_t,uint16_t,uint32_t));
uint32_t OperatorRemoveSmallValues();
uint32_t QDetNotOrderedFaFcBaBc(const uint16_t A, const uint16_t B, const uint32_t NFermion);
uint32_t QDetNotOrderedFcFaBcBa(const uint16_t A, const uint16_t B, const uint32_t NFermion);
uint32_t QDetNotOrderedFcFaBaBc(const uint16_t A, const uint16_t B, const uint32_t NFermion);
uint32_t QDetNotOrderedFaFcBcBa(const uint16_t A, const uint16_t B, const uint32_t NFermion);
uint32_t NormalOrderOperator();
uint32_t OperatorMultiply();
uint32_t OperatorScale();
uint32_t OperatorAdd();
uint32_t MappedOperatorAddAtLevel();
uint32_t MappedOperatorAdd();

/**********************************************************/
/********       Core/OperatorWaveFunction.c        ********/
/**********************************************************/

double   OperatorCr(uint8_t *const Determinant, const uint16_t i, const uint32_t NFermion);
double   OperatorAn(uint8_t *const Determinant, const uint16_t i, const uint32_t NFermion);
double   OperatorLader(uint8_t *const Determinant, const uint16_t *const o, const uint32_t lenght, const uint32_t NFermion);
uint32_t OperatorPsiConserveBasisCC();
uint32_t OperatorPsiConserveBasisCR();
uint32_t OperatorPsiConserveBasisRC();
uint32_t OperatorPsiConserveBasisRR();
uint32_t OperatorPsiConserveBasisMC();
uint32_t OperatorPsiConserveBasisMR();
uint32_t OperatorPsiConserveBasisPreAllocatedCC();
uint32_t OperatorPsiConserveBasisPreAllocatedCR();
uint32_t OperatorPsiConserveBasisPreAllocatedRC();
uint32_t OperatorPsiConserveBasisPreAllocatedRR();
uint32_t OperatorPsiConserveBasisPreAllocatedMC();
uint32_t OperatorPsiConserveBasisPreAllocatedMR();
uint32_t OperatorPsiCC();
uint32_t OperatorPsiCR();
uint32_t OperatorPsiRC();
uint32_t OperatorPsiRR();
uint32_t OperatorPsiMC();
uint32_t OperatorPsiMR();
uint32_t OperatorPsiSerialCC();
uint32_t OperatorPsiSerialCR();
uint32_t OperatorPsiSerialRC();
uint32_t OperatorPsiSerialRR();
uint32_t OperatorPsiSerialMC();
uint32_t OperatorPsiSerialMR();
uint32_t PsiOperatorPsiCC();
uint32_t PsiOperatorPsiCR();
uint32_t PsiOperatorPsiRC();
uint32_t PsiOperatorPsiRR();
uint32_t PsiOperatorPsiMC();
uint32_t PsiOperatorPsiMR();
uint32_t OperatorPsiConserveBasis();
uint32_t OperatorPsiConserveBasisPreAllocated();
uint32_t OperatorPsi();
uint32_t OperatorPsiSerial();
uint32_t PsiOperatorPsi(const OperatorType O, const WaveFunctionType psi, double *VRe, double *VIm);

/**********************************************************/
/*****    Core/OperatorWaveFunctionPerturbation.c     *****/
/**********************************************************/

uint32_t OperatorPsiPerturbationConserveBasis(const OperatorPerturbationType O, const WaveFunctionPerturbationType psi, WaveFunctionPerturbationType *const Opsi);
uint32_t OperatorPsiPerturbationConserveBasisPreAllocated(const OperatorPerturbationType O, const WaveFunctionPerturbationType psi, WaveFunctionPerturbationType *const Opsi);
uint32_t OperatorPsiPerturbation(const OperatorPerturbationType O, const WaveFunctionPerturbationType psi, WaveFunctionPerturbationType *const Opsi);
uint32_t PsiOperatorPsiPerturbation(const OperatorPerturbationType O, const WaveFunctionPerturbationType psi, double *VRe, double *VeIm);
uint32_t PsiHermitianOperatorPsiPerturbation(const OperatorPerturbationType O, const WaveFunctionPerturbationType psi, double *VRe);
uint32_t PsiOperatorSqrPsiPerturbation(const OperatorPerturbationType O, const WaveFunctionPerturbationType psi, double *VRe, double *VeIm);
uint32_t PsiHermitianOperatorSqrPsiPerturbation(const OperatorPerturbationType O, const WaveFunctionPerturbationType psi, double *VRe);
uint32_t PsiPrimeOperatorPsiPerturbationMatrix(const OperatorPerturbationType O, const WaveFunctionPerturbationType *psiPrime, WaveFunctionPerturbationType *const psi, MatrixType *const V);
uint32_t PsiOperatorPsiPerturbationMatrix(const OperatorPerturbationType O, const WaveFunctionPerturbationType *psi, MatrixType *const V);
uint32_t PsiPrimeOperatorPsiPerturbation(const OperatorPerturbationType O, const WaveFunctionPerturbationType psiPrime, WaveFunctionPerturbationType const psi, double *VRe, double *VIm);



/**********************************************************/
/**************        Core/Spectra.c        **************/
/**********************************************************/
void CalculateMoment(SpectraType *const G);


/**********************************************************/
/**************     Core/DensityMatrix.c     **************/
/**********************************************************/
uint32_t CalculateDensity(double *const Rho, const WaveFunctionType psi);


/**********************************************************/
/********      Hubbard Model/CellFunctions.c      *********/
/**********************************************************/
uint32_t AddAditionalInformationToSuperCell(CrystalStructureSuperCellType * A);
uint32_t CopyCrystalStructure(CrystalStructureType A, CrystalStructureType *B);


/**********************************************************/
/********    Tight Binding/Tight_Binding_lib.c    *********/
/**********************************************************/
int SetNFermionIndexTightBindingDef(TightBindingDef *A);
uint32_t CreateUSperCell(OperatorType U, OperatorType *USuperCell, CrystalStructureSuperCellType SuperCell);
int CreateTightBindingSuperCell();



/**********************************************************/
/********    DMFT.c    *********/
/**********************************************************/
uint32_t DMFTSpec();









