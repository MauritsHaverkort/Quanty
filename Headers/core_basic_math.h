/**********************************************************/
/******       BasisMath/BasicMatch.c       ******/
/**********************************************************/
uint32_t PowerInt(uint32_t x, uint32_t a);
double PowerDoubleInt(double x, uint32_t a);
void Cross(const double x[], const double y[], double z[]);
double ParallelalPipedVolume(const double x[], const double y[], const double z[]);
void Swap(double *const a, double *const b);
double CoVariance(double x[], double y[], int N);
double Variance(double x[], int N);
double Avverage(double x[], int N);
void LinearInterpolate(double x[], double y[], int N, double *a, double *b);
void QuadraticInterpolate(double x[], double y[], int N, double *a, double *b, double *c);
double InterpolateRoot(double const *x, double const *y, uint32_t N);
double InterpolateRootOfset(double const *x, double const *y, uint32_t N, double const ofset);
void InterpolateRoot2D(double const *x1, double const *x2, double const *y1, double const *y2, uint32_t N, double *const x1next, double *const x2next);
void InterpolateRoot2DOfset(double const *x1, double const *x2, double const *y1, double const *y2, uint32_t N, double const ofsety1, double const ofsety2, double *const x1next, double *const x2next);
double FindSelfconsistentnumber(double const *valin, double const *valout, uint32_t N);
void ReciprocalLatticeVector(const double r[3][3], double g[3][3], uint32_t const Units);


/**********************************************************/
/******       BasisMath/AndersonMatrixVector.c       ******/
/**********************************************************/
void   FreeAndersonMatrix(AndersonMatrixType * const M);
uint32_t InitAndersonMatrix(AndersonMatrixType * const M);
uint32_t InitIdentityAndersonMatrix(AndersonMatrixType * const M);
void   AndersonMatrixChop(AndersonMatrixType * const M, double const epsilon);
double AndersonMatrixInverse(AndersonMatrixType const M, MatrixType * const MI);
void   AndersonMatrixAdd(AndersonMatrixType * const M, AndersonMatrixType B, double alphaRe);
void   AndersonMatrixScalarMultiply(AndersonMatrixType * const M, double alphaRe);
uint32_t AndersonmatrixPart(AndersonMatrixType * const M, AndersonMatrixType A, const uint32_t iStart, const uint32_t iEnd);
uint32_t AndersonmatrixCopy(AndersonMatrixType * const M, AndersonMatrixType A);
uint32_t AndersonMatrixOrder(AndersonMatrixType * const M, AndersonMatrixType A);
uint32_t AndersonMatrixOrderAndAbs(AndersonMatrixType * const M, AndersonMatrixType A);
void   RealAndersonMatrixVector(AndersonMatrixType const M, const double *const x, double *const y);
void   ComplexAndersonMatrixVector(AndersonMatrixType const M, const double *const xRe, const double *const xIm, double *const yRe, double *const yIm);
void   RealVectorAndersonMatrixVector(AndersonMatrixType const M, const double *const x, double *const y);
void   ComplexVectorAndersonMatrixVector(AndersonMatrixType const M, const double *const xRe, const double *const xIm, double *const yRe, double *const yIm);
uint32_t AndersonMatrixCalculateChemicalPotential(AndersonMatrixType const M, double const occupation, double *const mu);

/**********************************************************/
/****       BasisMath/TriDiagonalMatrixVector.c        ****/
/**********************************************************/
void   FreeTriDiagonalMatrix(TriDiagonalMatrixType * const M);
uint32_t InitTriDiagonalMatrix(TriDiagonalMatrixType * const M);
uint32_t InitIdentityTriDiagonalMatrix(TriDiagonalMatrixType * const M);
void   TriDiagonalMatrixChop(TriDiagonalMatrixType * const M, double const epsilon);
double TriDiagonalMatrixInverse(TriDiagonalMatrixType const M, MatrixType * const MI);
void   TriDiagonalMatrixAdd(TriDiagonalMatrixType * const M, TriDiagonalMatrixType B, double alphaRe);
void   TriDiagonalMatrixScalarMultiply(TriDiagonalMatrixType * const M, double alphaRe);
uint32_t TriDiagonalMatrixPart(TriDiagonalMatrixType * const M, TriDiagonalMatrixType A, const uint32_t iStart, const uint32_t iEnd);
uint32_t TriDiagonalMatrixCopy(TriDiagonalMatrixType * const M, TriDiagonalMatrixType A);
void   RealTriDiagonalMatrixVector(TriDiagonalMatrixType const M, const double *const x, double *const y);
void   ComplexTriDiagonalMatrixVector(TriDiagonalMatrixType const M, const double *const xRe, const double *const xIm, double *const yRe, double *const yIm);
void   RealVectorTriDiagonalMatrixVector(TriDiagonalMatrixType const M, const double *const x, double *const y);
void   ComplexVectorTriDiagonalMatrixVector(TriDiagonalMatrixType const M, const double *const xRe, const double *const xIm, double *const yRe, double *const yIm);
uint32_t SumGreensfunctionsTriDiagonalRepresentation(TriDiagonalMatrixType const A, TriDiagonalMatrixType const B, TriDiagonalMatrixType * const M);

/**********************************************************/
/********         BasisMath/ListOfPoles.c          ********/
/**********************************************************/
uint32_t InitListOfPoles(ListOfPoles * M);
void   FreeListOfPoles(ListOfPoles * M);
void   ScaleListOfPoles(ListOfPoles * const A, double const scale);
uint32_t AddListOfPoles(ListOfPoles * const A, ListOfPoles B);
uint32_t MixListOfPoles(ListOfPoles * const A, ListOfPoles const B, double const alpha);
uint32_t SymmetrizeListOfPoles(ListOfPoles * const A);
void   OrderListOfPoles(ListOfPoles * const M);
void   ReduceListOfPolesWeigth(ListOfPoles * M, double const bmin, double const Ef);
void   ReduceListOfPolesNumberFermi(ListOfPoles * M, uint32_t N, double Ef);
void   ReduceListOfPolesNumber(ListOfPoles * M, uint32_t N);
void   ReduceListOfPolesMoment(ListOfPoles * M, uint32_t N);
void   ReduceListOfPolesSpacing(ListOfPoles * M, double da);
uint32_t SubtractListOfPoles(ListOfPoles * const A, const ListOfPoles B);
uint32_t GaussianBroadenListOfPoles(ListOfPoles * const M, double const G, uint32_t const Nmul);
uint32_t ListOfPolesCopy(ListOfPoles * M, ListOfPoles B);
uint32_t ListOfPolesCopyVal(ListOfPoles * M, ListOfPoles B, double Ef, double Epsilon);
uint32_t ListOfPolesCopyCon(ListOfPoles * M, ListOfPoles B, double Ef, double Epsilon);
uint32_t ListOfPolesCopyValOccupation(ListOfPoles * M, ListOfPoles B, double *Ef, double occupation);
uint32_t ListOfPolesCopyConOccupation(ListOfPoles * M, ListOfPoles B, double *Ef, double occupation);
uint32_t ListOfPolesDetermineOccupation(ListOfPoles M, double const Ef, double *const occupation);
uint32_t ListOfPoleCopyChoped(ListOfPoles M, ListOfPoles * A, double Epsilon);
uint32_t ShiftListOfPolesToMakeChemPotZeroMetal(ListOfPoles * const M, double const occupation);
uint32_t ShiftListOfPolesToMakeChemPotZeroInsulator(ListOfPoles * const M, double const occupation);
uint32_t ShiftListOfPolesToMakeChemPotZero(ListOfPoles * const M, double const occupation, double const epsilon);
uint32_t ShiftListOfPolesToMakeFirstMomentZero(ListOfPoles * const M);
uint32_t ListOfPolesCalculateChemicalPotential(ListOfPoles const M, double const occupation, double *const mu);

/**********************************************************/
/**********       BasisMath/MatrixVector.c       **********/
/**********************************************************/
void   FreeMatrix(MatrixType * const M);
uint32_t InitMatrix(MatrixType * const M);
uint32_t InitIdentityMatrix(MatrixType * const M);
uint32_t MakeRealMatrixComplex(MatrixType * const M);
uint32_t MakeComplexMatrixReal(MatrixType * const M);
void   MatrixChop(MatrixType * const M, double const epsilon);
uint32_t MatrixTranspose(MatrixType * const M);
void   MatrixComplexConjugate(MatrixType * const M);
uint32_t MatrixConjugateTranspose(MatrixType * M);
double MatrixInverse(MatrixType * M);
uint32_t MatrixAdd(MatrixType * const M, const MatrixType B, double const alphaRe, double const alphaIm);
uint32_t MatrixScalarMultiply(MatrixType * const M, double const alphaRe, double const alphaIm);
uint32_t MatrixTransform(MatrixType * const M, MatrixType const A, MatrixType const T);
uint32_t MatrixTransformPart(MatrixType * const M, MatrixType const A, MatrixType const T, const uint32_t MNiMin, const uint32_t MNiMax, const uint32_t MNjMin, const uint32_t MNjMax);
uint32_t MatrixMultiply(MatrixType * const M, MatrixType const A, MatrixType const B);
uint32_t MatrixPart(MatrixType * const M, MatrixType const A, const uint32_t iStart, const uint32_t iEnd, const uint32_t jStart, const uint32_t jEnd);
uint32_t MatrixCopy(MatrixType * const M, MatrixType const A);
uint32_t RealMatrixVector(MatrixType const M, const double *const a, double *const b);
uint32_t RealTransposeMatrixVector(MatrixType const M, const double *const a, double *const b);
uint32_t ComplexMatrixVector(MatrixType const M, const double *const aRe, const double *const aIm, double *const bRe, double *const bIm);
uint32_t ComplexConjugateTransposeMatrixVector(MatrixType const M, const double *const aRe, const double *const aIm, double *const bRe, double *const bIm);
uint32_t RealVectorMatrixVector(MatrixType const M, const double *const a, double *const b);
uint32_t ComplexVectorMatrixVector(MatrixType const M, const double *const aRe, const double *const aIm, double *const bRe, double *const bIm);
void   MatrixOrthonormalize(MatrixType * const M);



/**********************************************************/
/**********    BasisMath/CompactMatrixVector.c   **********/
/**********************************************************/
void   FreeCompactMatrix(CompactMatrixType * const M);
uint32_t InitCompactMatrix(CompactMatrixType * const M);
uint32_t InitCompactIdentityMatrix(CompactMatrixType * const M);
uint32_t MakeCompactRealMatrixComplex(CompactMatrixType * const M);
uint32_t CompactMatrixChop(CompactMatrixType * const M, double const epsilon);
uint32_t CompactMatrixConjugate(CompactMatrixType * const M);
uint32_t CompactMatrixTranspose(CompactMatrixType * const M);
uint32_t CompactMatrixConjugateTranspose(CompactMatrixType * const M);
double CompactMatrixInverse(CompactMatrixType * const M);

/**********************************************************/
/*******       BasisMath/SparseMatrixVector.c       *******/
/**********************************************************/
void   FreeSparseMatrix(SparseMatrixType * const M);
uint32_t InitSparseMatrix(SparseMatrixType * const M);
uint32_t InitSparseIdentityMatrix(SparseMatrixType * const M);
uint32_t MakeRealSparseMatrixComplex(SparseMatrixType * const M);
uint32_t SparseMatrixAddElement(SparseMatrixType * const M, uint32_t i, uint32_t j, double VRe, double VIm);
void   SparseMatrixExtractElement(SparseMatrixType const M, uint32_t const i, uint32_t const j, double *VRe, double *VIm);
void   MinimalAlocSpareseMatrix(SparseMatrixType * const M);
void   RemoveZeroSpareseMatrix(SparseMatrixType * const M, double epsilon);
void   SparseMatrixChop(SparseMatrixType * const M, double const epsilon);
void   SparseMatrixTranspose(SparseMatrixType * const M);
void   SparseMatrixComplexConjugate(SparseMatrixType * const M);
void   SparseMatrixConjugateTranspose(SparseMatrixType * const M);
double SparseMatrixInverse(SparseMatrixType * const M, double epsilon);
uint32_t SparseMatrixAdd(SparseMatrixType * const M, const SparseMatrixType B, double const alphaRe, double const alphaIm);
uint32_t SparseMatrixScalarMultiply(SparseMatrixType * const M, double const alphaRe, double const alphaIm);
uint32_t SparseMatrixTransform(SparseMatrixType * const M, SparseMatrixType const A, SparseMatrixType const T, double epsilon);
uint32_t SparseMatrixDenseTransform(SparseMatrixType * const M, SparseMatrixType const A, MatrixType const T, double epsilon);
uint32_t DenseMatrixSparseTransform(SparseMatrixType * const M, MatrixType const A, SparseMatrixType const T, double epsilon);
uint32_t SparseMatrixTransformPart(SparseMatrixType * const M, SparseMatrixType const A, SparseMatrixType const T, double epsilon, const uint32_t MNiMin, const uint32_t MNiMax, const uint32_t MNjMin, const uint32_t MNjMax);
uint32_t SparseMatrixDenseTransformPart(SparseMatrixType * const M, SparseMatrixType const A, MatrixType const T, double epsilon, const uint32_t MNiMin, const uint32_t MNiMax, const uint32_t MNjMin, const uint32_t MNjMax);
uint32_t DenseMatrixSparseTransformPart(SparseMatrixType * const M, MatrixType const A, SparseMatrixType const T, double epsilon, const uint32_t MNiMin, const uint32_t MNiMax, const uint32_t MNjMin, const uint32_t MNjMax);
uint32_t SparseMatrixMultiply(SparseMatrixType * const M, SparseMatrixType const A, SparseMatrixType const B);
uint32_t SparseMatrixDenseMultiply(SparseMatrixType * const M, SparseMatrixType const A, MatrixType const B);
uint32_t DenseMatrixSparseMultiply(SparseMatrixType * const M, MatrixType const A, SparseMatrixType const B);
uint32_t SparseMatrixPart(SparseMatrixType * const M, SparseMatrixType const A, const uint32_t iStart, const uint32_t iEnd, const uint32_t jStart, const uint32_t jEnd);
uint32_t SparseMatrixCopy(SparseMatrixType * const M, SparseMatrixType const A);
uint32_t RealSparseMatrixVector(SparseMatrixType const M, const double *const a, double *const b);
uint32_t RealTransposeSparseMatrixVector(SparseMatrixType const M, const double *const a, double *const b);
void   ComplexSparseMatrixVector(SparseMatrixType const M, const double *const aRe, const double *const aIm, double *const bRe, double *const bIm);
void   ComplexConjugateTransposeSparseMatrixVector(SparseMatrixType const M, const double *const aRe, const double *const aIm, double *const bRe, double *const bIm);
uint32_t RealVectorSparseMatrixVector(SparseMatrixType const M, const double *const a, double *const b);
void   ComplexVectorSparseMatrixVector(SparseMatrixType const M, const double *const aRe, const double *const aIm, double *const bRe, double *const bIm);
void   SparseMatrixOrthonormalize(MatrixType * const M, double epsilon);




/**********************************************************/
/**      BasisMath/DoubleTriDiagonalMatrixVector.c       **/
/**********************************************************/
void   FreeDoubleTriDiagonalMatrix(DoubleTriDiagonalMatrixType * const M);
uint32_t InitDoubleTriDiagonalMatrix(DoubleTriDiagonalMatrixType * const M);
uint32_t InitIdentityDoubleTriDiagonalMatrix(DoubleTriDiagonalMatrixType * const M);
void   DoubleTriDiagonalMatrixChop(DoubleTriDiagonalMatrixType * const M, double const epsilon);
double DoubleTriDiagonalMatrixInverse(DoubleTriDiagonalMatrixType const M, MatrixType * const MI);
void   DoubleTriDiagonalMatrixAdd(DoubleTriDiagonalMatrixType * const M, DoubleTriDiagonalMatrixType B, double alphaRe);
void   DoubleTriDiagonalMatrixScalarMultiply(DoubleTriDiagonalMatrixType * const M, double alphaRe);
uint32_t DoubleTriDiagonalMatrixPart(DoubleTriDiagonalMatrixType * const M, DoubleTriDiagonalMatrixType A, const uint32_t iEndval, const uint32_t iEndcon);
uint32_t DoubleTriDiagonalMatrixCopy(DoubleTriDiagonalMatrixType * const M, DoubleTriDiagonalMatrixType A);
uint32_t RealDoubleTriDiagonalMatrixVector(DoubleTriDiagonalMatrixType const M, const double *const x, double *const y);
uint32_t ComplexDoubleTriDiagonalMatrixVector(DoubleTriDiagonalMatrixType const M, const double *const xRe, const double *const xIm, double *const yRe, double *const yIm);
uint32_t RealVectorDoubleTriDiagonalMatrixVector(DoubleTriDiagonalMatrixType const M, const double *const x, double *const y);
uint32_t ComplexVectorDoubleTriDiagonalMatrixVector(DoubleTriDiagonalMatrixType const M, const double *const xRe, const double *const xIm, double *const yRe, double *const yIm);



/**********************************************************/
/**       BasisMath/BlokTriDiagonalMatrixVector.c        **/
/**********************************************************/
void   FreeBlokTriDiagonalMatrix(BlokTriDiagonalMatrixType * const M);
void   FreeFrameBlokTriDiagonalMatrix(BlokTriDiagonalMatrixType * const M);
uint32_t InitBlokTriDiagonalMatrix(BlokTriDiagonalMatrixType * const M, uint32_t bloksize, uint32_t QComplex);
uint32_t InitFrameBlokTriDiagonalMatrix(BlokTriDiagonalMatrixType * const M);
uint32_t InitIdentityBlokTriDiagonalMatrix(BlokTriDiagonalMatrixType * const M, uint32_t bloksize, uint32_t QComplex);
void   BlokTriDiagonalMatrixChop(BlokTriDiagonalMatrixType * const M, double const epsilon);
void   BlokTriDiagonalMatrixAdd(BlokTriDiagonalMatrixType * const M, BlokTriDiagonalMatrixType B, double const alphaRe);
void   BlokTriDiagonalMatrixScalarMultiply(BlokTriDiagonalMatrixType * const M, double alphaRe);
uint32_t BlokTriDiagonalMatrixPart(BlokTriDiagonalMatrixType * const M, BlokTriDiagonalMatrixType A, uint32_t iStart, uint32_t iEnd);
uint32_t BlokTriDiagonalMatrixCopy(BlokTriDiagonalMatrixType * const M, BlokTriDiagonalMatrixType A);
void   BlokTriDiagonalMatrixTranspose(BlokTriDiagonalMatrixType * const M);
void   BlokTriDiagonalMatrixComplexConjugate(BlokTriDiagonalMatrixType * const M);
void   BlokTriDiagonalMatrixConjugateTranspose(BlokTriDiagonalMatrixType * const M);
uint32_t RealBlokTriDiagonalMatrixVector(BlokTriDiagonalMatrixType const M, const double *const x, double *const y);
uint32_t ComplexBlokTriDiagonalMatrixVector(BlokTriDiagonalMatrixType const M, const double *const xRe, const double *const xIm, double *const yRe, double *const yIm);
uint32_t RealVectorBlokTriDiagonalMatrixVector(BlokTriDiagonalMatrixType const M, const double *const x, double *const y);
uint32_t ComplexVectorBlokTriDiagonalMatrixVector(BlokTriDiagonalMatrixType const M, const double *const xRe, const double *const xIm, double *const yRe, double *const yIm);


/**********************************************************/
/**         BasisMath_BlockTriDiagonalMatrix.c           **/
/**********************************************************/
uint32_t InitBlockTridiagonalMatrix(BlockTriDiagonalMatrixType *const A);
uint32_t FreeBlockTridiagonalMatrix(BlockTriDiagonalMatrixType *A);
uint32_t FreeBlockMatrixFromLua(ResponsefunctionType *M);


/**********************************************************/
/****      BasisMath/BlokAndersonMatrixVector.c        ****/
/**********************************************************/
void   FreeBlokAndersonMatrix(BlokAndersonMatrixType * const M);
void   FreeFrameBlokAndersonMatrix(BlokAndersonMatrixType * const M);
uint32_t InitBlokAndersonMatrix(BlokAndersonMatrixType * const M, uint32_t bloksize, uint32_t QComplex);
uint32_t InitFrameBlokAndersonMatrix(BlokAndersonMatrixType * const M);
uint32_t InitIdentityBlokAndersonMatrix(BlokAndersonMatrixType * const M, uint32_t bloksize, uint32_t QComplex);
void   BlokAndersonMatrixChop(BlokAndersonMatrixType * const M, double const epsilon);
void   BlokAndersonMatrixAdd(BlokAndersonMatrixType * const M, BlokAndersonMatrixType B, double alphaRe);
void   BlokAndersonMatrixScalarMultiply(BlokAndersonMatrixType * const M, double alphaRe);
uint32_t BlokAndersonMatrixPart(BlokAndersonMatrixType * const M, BlokAndersonMatrixType A, const uint32_t iStart, const uint32_t iEnd);
uint32_t BlokAndersonMatrixCopy(BlokAndersonMatrixType * const M, BlokAndersonMatrixType A);
void   BlokAndersonMatrixTranspose(BlokAndersonMatrixType * const M);
void   BlokAndersonMatrixComplexConjugate(BlokAndersonMatrixType * const M);
void   BlokAndersonMatrixConjugateTranspose(BlokAndersonMatrixType * const M);
uint32_t RealBlokAndersonMatrixVector(BlokAndersonMatrixType const M, const double *const x, double *const y);
uint32_t ComplexBlokAndersonMatrixVector(BlokAndersonMatrixType const M, const double *const xRe, const double *const xIm, double *const yRe, double *const yIm);
uint32_t RealVectorBlokAndersonMatrixVector(BlokAndersonMatrixType const M, const double *const x, double *const y);
uint32_t ComplexVectorBlokAndersonMatrixVector(BlokAndersonMatrixType const M, const double *const xRe,const double *const xIm, double *const yRe, double *const yIm);


/**********************************************************/
/**    BasisMath/BlokDoubleTriDiagonalMatrixVector.c     **/
/**********************************************************/
void   FreeBlokDoubleTriDiagonalMatrix(BlokDoubleTriDiagonalMatrixType * const M);
void   FreeFrameBlokDoubleTriDiagonalMatrix(BlokDoubleTriDiagonalMatrixType * const M);
uint32_t InitBlokDoubleTriDiagonalMatrix(BlokDoubleTriDiagonalMatrixType * const M, uint32_t bloksize, uint32_t QComplex);
uint32_t InitFrameBlokDoubleTriDiagonalMatrix(BlokDoubleTriDiagonalMatrixType * const M);
uint32_t InitIdentityBlokDoubleTriDiagonalMatrix(BlokDoubleTriDiagonalMatrixType * const M, uint32_t bloksize, uint32_t QComplex);
uint32_t BlokDoubleTriDiagonalMatrixChop(BlokDoubleTriDiagonalMatrixType * const M, double const epsilon);
uint32_t BlokDoubleTriDiagonalMatrixAdd(BlokDoubleTriDiagonalMatrixType * const M, BlokDoubleTriDiagonalMatrixType B, double const alphaRe);
uint32_t BlokDoubleTriDiagonalMatrixScalarMultiply(BlokDoubleTriDiagonalMatrixType * const M, double alphaRe);
uint32_t BlokDoubleTriDiagonalMatrixPart(BlokDoubleTriDiagonalMatrixType * const M, BlokDoubleTriDiagonalMatrixType A, uint32_t iEndval, uint32_t iEndcon);
uint32_t BlokDoubleTriDiagonalMatrixCopy(BlokDoubleTriDiagonalMatrixType * const M, BlokDoubleTriDiagonalMatrixType A);
uint32_t BlokDoubleTriDiagonalMatrixTranspose(BlokDoubleTriDiagonalMatrixType * const M);
uint32_t BlokDoubleTriDiagonalMatrixComplexConjugate(BlokDoubleTriDiagonalMatrixType * const M);
uint32_t BlokDoubleTriDiagonalMatrixConjugateTranspose(BlokDoubleTriDiagonalMatrixType * const M);
uint32_t RealBlokDoubleTriDiagonalMatrixVector(BlokDoubleTriDiagonalMatrixType const M, double *const x, double *const y);
uint32_t ComplexBlokDoubleTriDiagonalMatrixVector(BlokDoubleTriDiagonalMatrixType const M, double *const xRe, double *const xIm, double *const yRe, double *const yIm);
uint32_t RealVectorBlokDoubleTriDiagonalMatrixVector(BlokDoubleTriDiagonalMatrixType const M, const double *const x, double *const y);
uint32_t ComplexVectorBlokDoubleTriDiagonalMatrixVector(BlokDoubleTriDiagonalMatrixType const M, const double *const xRe, const double *const xIm, double *const yRe, double *const yIm);



/**********************************************************/
/********    BasisMath/MatrixTypeConversion.c      ********/
/**********************************************************/
/********    AndersonMatrixTo                      ********/
uint32_t AndersonMatrixToMatrix(AndersonMatrixType const M, MatrixType * const Mfull);
uint32_t AndersonMatrixToTriDiagonalMatrixNonOrthogonal(AndersonMatrixType const M, TriDiagonalMatrixType * const H, uint32_t NTriDiag);
uint32_t AndersonMatrixToTriDiagonalMatrix(AndersonMatrixType const M, TriDiagonalMatrixType * const H, uint32_t NTriDiag);
uint32_t AndersonMatrixToDoubleTriDiagonalMatrixNonOrthogonal(AndersonMatrixType const M, DoubleTriDiagonalMatrixType * const H, double const Ef, double const epsilon);
uint32_t AndersonMatrixToDoubleTriDiagonalMatrix(AndersonMatrixType const M, DoubleTriDiagonalMatrixType * const H, double const Ef, double const epsilon);
uint32_t AndersonMatrixInverseToListOfPoles(AndersonMatrixType const M, ListOfPoles * const H);
uint32_t AndersonMatrixToEigensystem(AndersonMatrixType const M, MatrixType * const T, double *val);
uint32_t AndersonMatrixToListOfPoles(AndersonMatrixType const M, ListOfPoles * const A);
uint32_t AndersonMatrixToBADoubleTriDiagonalMatrix(AndersonMatrixType const M, BADoubleTriDiagonalMatrixType * const A, uint32_t NTriDiagVal, uint32_t NTriDiagCon, double Ef, double Epsilon);
uint32_t AndersonMatrixToBADoubleTriDiagonalMatrixOccupation(AndersonMatrixType const M, double const onsiteshift, BADoubleTriDiagonalMatrixType * const A, double *const occupation);
/********    MatrixTo                              ********/
//uint32_t BlokAndersonMatrixToSparseMatrix();
uint32_t MatrixToEigensystem(MatrixType const M, MatrixType * const T, double *val);
//uint32_t MatrixToDoubleTriDiagonalMatrix();
//uint32_t BlokAndersonMatrixToBlokTriDiagonalMatrix();
//uint32_t BlokAndersonMatrixToBlokDoubleTriDiagonalMatrix();
/********    TriDiagonalMatrixTo                   ********/
//uint32_t TriDiagonalMatrixToSparseMatrix();
uint32_t TriDiagonalMatrixToEigensystem(TriDiagonalMatrixType const M, MatrixType * const T, double *val);
uint32_t TriDiagonalMatrixToAndersonMatrix(TriDiagonalMatrixType const M, AndersonMatrixType * const A);
uint32_t TriDiagonalMatrixToDoubleTriDiagonalMatrix(TriDiagonalMatrixType const M, DoubleTriDiagonalMatrixType * const A, double Ef, double epsilon);
/********    BlokTriDiagonalMatrixTo               ********/
//uint32_t BlokTriDiagonalMatrixToSparseMatrix();
//uint32_t BlokTriDiagonalMatrixToEigensystem();
//uint32_t BlokTriDiagonalMatrixToBlokAndersonMatrix();
//uint32_t BlokTriDiagonalMatrixToBlokDoubleTriDiagonalMatrix();
/********    DoubleTriDiagonalMatrixTo             ********/
//uint32_t DoubleTriDiagonalMatrixToSparseMatrix();
uint32_t DoubleTriDiagonalMatrixToTriDiagonalMatrix(DoubleTriDiagonalMatrixType const M, TriDiagonalMatrixType * const H, uint32_t NTriDiag);
//uint32_t DoubleTriDiagonalMatrixToEigensystem();
//uint32_t DoubleTriDiagonalMatrixToAndersonMatrix();
/********    BlokDoubleTriDiagonalMatrixTo         ********/
//uint32_t BlokDoubleTriDiagonalMatrixToSparseMatrix();
//uint32_t BlokDoubleTriDiagonalMatrixToEigensystem();
//uint32_t BlokDoubleTriDiagonalMatrixToBlokTriDiagonalMatrix();
//uint32_t BlokDoubleTriDiagonalMatrixToBlokAndersonMatrix();
/********    SparseMatrixTo                        ********/
//uint32_t SparseMatrixToMatrix();
uint32_t SparseMatrixToTriDiagonalMatrix(SparseMatrixType M, MatrixType * const T, TriDiagonalMatrixType * const H, uint32_t NTriDiag);
//uint32_t SparseMatrixToEigensystem();
//uint32_t SparseMatrixToAndersonMatrix();
//uint32_t SparseMatrixToBlokTriDiagonalMatrix();
//uint32_t SparseMatrixToBlokAndersonMatrix();
//uint32_t SparseMatrixToDoubleTriDiagonalMatrix();
//uint32_t SparseMatrixToBlokDoubleTriDiagonalMatrix();
/********    MatrixTo                              ********/
uint32_t MatrixToSparseMatrix(MatrixType const M, SparseMatrixType * const U, double const epsilon);
uint32_t MatrixToTriDiagonalMatrix(MatrixType const M, MatrixType * const T, TriDiagonalMatrixType * const H);
uint32_t TriDiagonalize(MatrixType const M, MatrixType * const T, TriDiagonalMatrixType * const H);
uint32_t MatrixToEigensystem(MatrixType const M, MatrixType * const T, double *val);
uint32_t Eigensystem(MatrixType const M, MatrixType * const T, double *val);
uint32_t MatrixToGeneralizedEigensystem(MatrixType const M, MatrixType const S, MatrixType * const T, MatrixType * const Sroot, double *val);
uint32_t GeneralizedEigensystem(MatrixType const M, MatrixType const S, MatrixType * const T, MatrixType * const Sroot,  double *val);
uint32_t MatrixToAndersonMatrix(MatrixType const M, MatrixType * const T, AndersonMatrixType * const H);
uint32_t MatrixToBlokTriDiagonalMatrix(MatrixType const M, MatrixType * const T, BlokTriDiagonalMatrixType * const H, const uint32_t NBlok);
uint32_t BlokTriDiagonalize(MatrixType const M, MatrixType * const T, BlokTriDiagonalMatrixType * const H, const uint32_t BlokDim);
uint32_t MatrixToBlokAndersonMatrix(MatrixType const M, MatrixType * const T, BlokAndersonMatrixType * const H, const uint32_t BlokDim);
//uint32_t MatrixToDoubleTriDiagonalMatrix();
//uint32_t MatrixToBlokDoubleTriDiagonalMatrix();
/********    ListOfPolesTo                              ********/
uint32_t ListOfPolesToBADoubleTridiagonalMatrix(const ListOfPoles A, BADoubleTriDiagonalMatrixType * const M, uint32_t NTriDiagval, uint32_t NTriDiagcon, double mu, double Epsilon);
uint32_t ListOfPolesToBADoubleTridiagonalMatrixOccupation(const ListOfPoles A, BADoubleTriDiagonalMatrixType * const M, uint32_t NTriDiagval, uint32_t NTriDiagcon, double mu, double occupation);

/**********************************************************/
/********    BasisMath/MatrixTypeConversion.c      ********/
/**********************************************************/

uint32_t TriDiagonalMatrixToEigensystem(TriDiagonalMatrixType const M, MatrixType * const T, double *val);
uint32_t MatrixToSparseMatrix(MatrixType const M, SparseMatrixType * const U, double const epsilon);
uint32_t TriDiagonalize(MatrixType const M, MatrixType * const T, TriDiagonalMatrixType * const H);
uint32_t MatrixToTriDiagonalMatrix(MatrixType const M, MatrixType * const T, TriDiagonalMatrixType * const H);
uint32_t Eigensystem(MatrixType const M, MatrixType * const T, double *val);
uint32_t MatrixToEigensystem(MatrixType const M, MatrixType * const T, double *val);
uint32_t MatrixToAndersonMatrix(MatrixType const M, MatrixType * const T, AndersonMatrixType * const H);
uint32_t BlokTriDiagonalize(MatrixType const M, MatrixType * const T, BlokTriDiagonalMatrixType * const H, const uint32_t BlokDim);
uint32_t MatrixToBlokTriDiagonalMatrix(MatrixType const M, MatrixType * const T, BlokTriDiagonalMatrixType * const H, const uint32_t NBlok);
uint32_t MatrixToBlokAndersonMatrix(MatrixType const M, MatrixType * const T, BlokAndersonMatrixType * const H, const uint32_t BlokDim);
uint32_t MatrixToDoubleTriDiagonalMatrixNoOverlap(MatrixType const M, MatrixType * const Tval, MatrixType * const Tcon, DoubleTriDiagonalMatrixType * const H, double const Ef);
uint32_t MatrixToBlokDoubleTriDiagonalMatrixNoOverlap(MatrixType const M, MatrixType * const Tval, MatrixType * const Tcon, BlokDoubleTriDiagonalMatrixType * const H, const uint32_t BlokDim, double const Ef);
uint32_t TriDiagonalMatrixToSecondDoubleTriDiagonalMatrix(TriDiagonalMatrixType const M, SecondDoubleTriDiagonalMatrixType * const A, double Ef, double epsilon);
uint32_t MapSecondoubleTridiagonalMatrixToBA(SecondDoubleTriDiagonalMatrixType const A,  BADoubleTriDiagonalMatrixType * const M);


/**********************************************************/
/**    BasisMathBlock/BasicMath_BlockTypeConversion.c    **/
/**********************************************************/

uint32_t EvaluateBlockListOfPoles(BlockListOfPoles const *A, CompactMatrixType * Gw, const double omega, const double Gamma);
uint32_t EvaluateBlockAndersonMatrix(BlockListOfPoles const *A, CompactMatrixType * Gw, const double omega, const double Gamma);
uint32_t EvaluateBlockTriDiagonalMatrix(BlockTriDiagonalMatrixType const * A, CompactMatrixType * Gw, const double omega, const double Gamma);
uint32_t EvaluateBlockNaturalImpurityOrbital(BlockNaturalImpurityOrbitalMatrixType const *A, CompactMatrixType * Gw, const double omega, const double Gamma);
uint32_t BlockListOfPolesToSpectra(BlockListOfPoles const A, SpectraType *G, double Emin, double Emax, uint32_t N);
uint32_t BlockAndersonMatrixToSpectra(BlockAndersonMatrixType const A, SpectraType *G, double Emin, double Emax, uint32_t N);
uint32_t BlockTriDiagonalMatrixToSpectra(const double Energy, BlockTriDiagonalMatrixType const * A, SpectraType *const G, double const Emin, double const Emax, uint32_t const N, double const Gamma);
uint32_t BlockNaturalImpurityOrbitalToSpectra(BlockNaturalImpurityOrbitalMatrixType const * A, SpectraType *const G, double const Emin, double const Emax, uint32_t const N, double const Gamma);
uint32_t BlockListOfPolesInversedToCompactMatrix(BlockListOfPoles const * const A, CompactMatrixType * const B);
uint32_t TightBindingHamiltonianToBlockListOfPoles(TightBindingDef const A, int Nk[3], BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE);
uint32_t TightBindingHamiltonianMultiSiteToBlockListOfPoles(TightBindingDef const A, int Nk[3], uint32_t const Nvec, int const *vec, BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE);
uint32_t TightBindingHamiltonianPlusSelfEnergyToBlockListOfPoles(TightBindingDef const A, int Nk[3], BlockListOfPoles const Sigma, BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE);
uint32_t TightBindingHamiltonianPlusSelfEnergyMultiSiteToBlockListOfPoles(TightBindingDef const A, int Nk[3], uint32_t const Nvec, int const *vec, BlockListOfPoles const Sigma, BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE);
uint32_t TightBindingHamiltonianAddSpinPlusSelfEnergyToBlockListOfPoles(TightBindingDef const A, int Nk[3], BlockListOfPoles const Sigma, BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE);
uint32_t TightBindingHamiltonianAddSpinPlusSelfEnergyMultiSiteToBlockListOfPoles(TightBindingDef const A, int Nk[3], uint32_t const Nvec, int const *vec, BlockListOfPoles const Sigma, BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE);
uint32_t DiagonalMatrixBlockTridiagonalize(double const *A, CompactMatrixType const *psi, CompactMatrixType *const T, BlockTriDiagonalMatrixType *const ABlock, double NBlock, double SingularValue);
uint32_t BlockListOfPolesToBlockTridiagonalMatrix(BlockListOfPoles const A, BlockTriDiagonalMatrixType *const B, int const NBlockMax, double const SingularValue);
uint32_t BlockListOfPolesToBlockTridiagonalMatrixReduceDimension(BlockListOfPoles const A, BlockTriDiagonalMatrixType *B, CompactMatrixType const *AT,  int const NBlockMax, double const SingularValue);
uint32_t BlockListOfPolesToBlockNaturalImpurityOrbitalMatrix(BlockListOfPoles const *A, BlockNaturalImpurityOrbitalMatrixType * const B, double const kbT, uint32_t const NBlockMaxVal, uint32_t const NBlockMaxCon);
uint32_t BlockListOfPolesToBlockNaturalImpurityOrbitalMatrixReduceDimension(BlockListOfPoles const *A, BlockNaturalImpurityOrbitalMatrixType * const B, CompactMatrixType const *AT, double const kbT, uint32_t const NBlockMaxVal, uint32_t const NBlockMaxCon);
uint32_t BlockListOfPolesToBlockAndersonMatrix(BlockListOfPoles const A, BlockAndersonMatrixType *const B, int const NBlockMax);
uint32_t BlockListOfPolesToBlockAndersonMatrixReduceDimension(BlockListOfPoles const A, BlockAndersonMatrixType *B, CompactMatrixType const *AT, int const NBlockMax);
uint32_t BlockTriDiagonalMatrixEigenValues();
uint32_t BlockTriDiagonalMatrixEigenSystem();
uint32_t BlockTriDiagonalMatrixToBlockListOfPoles(BlockTriDiagonalMatrixType const * A, BlockListOfPoles * const G);
uint32_t BlockTriDiagonalMatrixToBlockListOfPolesVariableGrid(BlockTriDiagonalMatrixType const * A, BlockListOfPoles * const G, const double * const Grid, const uint32_t NGridPoints);
uint32_t BlockTriDiagonalMatrixToBlockListOfPolesFixedGrid(BlockTriDiagonalMatrixType const * A, BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE);
uint32_t BlockTriDiagonalMatrixToBlockAndersonMatrix(BlockTriDiagonalMatrixType const ATri, BlockAndersonMatrixType * const B);
uint32_t BlockTriDiagonalMatrixToBlockNaturalImpurityMatrix(BlockTriDiagonalMatrixType const ATri, BlockNaturalImpurityOrbitalMatrixType * const B);
uint32_t BlockAndersonMatrixToBlockListOfPoles(BlockAndersonMatrixType const AAnd, BlockListOfPoles * const B);
uint32_t BlockAndersonMatrixToBlockListOfPolesVariableGrid(BlockAndersonMatrixType const AAnd, BlockListOfPoles * const B, const double * const Grid, const uint32_t NGridPoints);
uint32_t BlockAndersonMatrixToBlockTridiagonalMatrix(BlockAndersonMatrixType const AAnd, BlockTriDiagonalMatrixType * const B);
uint32_t BlockAndersonMatrixToBlockNaturalImpurityOrbitalMatrix(BlockAndersonMatrixType const AAnd, BlockNaturalImpurityOrbitalMatrixType * const B);
uint32_t BlockNaturalImpurityOrbitalMatrixToBlockAndersonMatrix(BlockNaturalImpurityOrbitalMatrixType const AAnd, BlockAndersonMatrixType * const B);
uint32_t BlockNaturalImpurityOrbitalMatrixToBlockListOfPoles(BlockNaturalImpurityOrbitalMatrixType const AAnd, BlockListOfPoles * const B);
uint32_t BlockNaturalImpurityOrbitalMatrixToBlockListOfPolesVariableGrid(BlockNaturalImpurityOrbitalMatrixType const AAnd, BlockListOfPoles * const B, const double * const Grid, const uint32_t NGridPoints);
uint32_t BlockNaturalImpurityOrbitalMatrixToBlockTridiagonalMatrix(BlockNaturalImpurityOrbitalMatrixType const AAnd, BlockTriDiagonalMatrixType * const B);
ResponsefunctionType *ConvertBlockMatrix(ResponsefunctionType *M, const char type, double T);
uint32_t CopyBlockNaturalImpurityOrbitalMatrix(BlockNaturalImpurityOrbitalMatrixType A, BlockNaturalImpurityOrbitalMatrixType *B);












/**********************************************************/
/***********     BasisMath/Combinatorial.c      ***********/
/**********************************************************/
double Factorial(int a);
double Factorial2(int a);
uint32_t BinomialCoefficient(uint32_t n, uint32_t k);
uint32_t PerturbationOptions(uint32_t NV, uint32_t order);
void   PerturbationOrderListToIndex(uint32_t * const POList, uint32_t * const POIndex, uint32_t const NV, uint32_t const order);
void   PerturbationOrderIndexToList(uint32_t const POIndex, uint32_t * const POList, uint32_t const NV, uint32_t const order);
double ThreeJSymbol(int j1, int m1, int j2, int m2, int j3, int m3);
double ThreeJSymbolHalf(int j1, int m1, int j2, int m2, int j3, int m3);
double SlaterCoefficientC(int k, int l1, int m1, int l2, int m2);



/**********************************************************/
/***********        BasisMath/Complex           ***********/
/**********************************************************/

Complex16 FaddeevaW(Complex16 z);
double FaddeevaW_Im(double x);
double erfcx(const double x);
