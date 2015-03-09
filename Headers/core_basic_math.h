/**********************************************************/
/******       BasisMath/BasicMatch.c       ******/
/**********************************************************/
uint32_t PowerInt();
double PowerDoubleInt();
void Cross();
double ParallelalPipedVolume();
void Swap();
double CoVariance();
double Variance();
double Avverage();
void LinearInterpolate();
void QuadraticInterpolate();
double InterpolateRoot();
double InterpolateRootOfset();
void InterpolateRoot2D();
void InterpolateRoot2DOfset();
double FindSelfconsistentnumber();
void ResiprocalLatticeVector();


/**********************************************************/
/******       BasisMath/AndersonMatrixVector.c       ******/
/**********************************************************/
void   FreeAndersonMatrix();
uint32_t InitAndersonMatrix();
uint32_t InitIdentityAndersonMatrix();
void   AndersonMatrixChop();
double AndersonMatrixInverse();
void   AndersonMatrixAdd();
void   AndersonMatrixScalarMultiply();
uint32_t AndersonmatrixPart();
uint32_t AndersonmatrixCopy();
uint32_t AndersonMatrixOrder();
uint32_t AndersonMatrixOrderAndAbs();
void   RealAndersonMatrixVector();
void   ComplexAndersonMatrixVector();
void   RealVectorAndersonMatrixVector();
void   ComplexVectorAndersonMatrixVector();
uint32_t AndersonMatrixCalculateChemicalPotential();

/**********************************************************/
/****       BasisMath/TriDiagonalMatrixVector.c        ****/
/**********************************************************/
void   FreeTriDiagonalMatrix();
uint32_t InitTriDiagonalMatrix();
uint32_t InitIdentityTriDiagonalMatrix();
void   TriDiagonalMatrixChop();
double TriDiagonalMatrixInverse();
void   TriDiagonalMatrixAdd();
void   TriDiagonalMatrixScalarMultiply();
uint32_t TriDiagonalMatrixPart();
uint32_t TriDiagonalMatrixCopy();
void   RealTriDiagonalMatrixVector();
void   ComplexTriDiagonalMatrixVector();
void   RealVectorTriDiagonalMatrixVector();
void   ComplexVectorTriDiagonalMatrixVector();
uint32_t SumGreensfunctionsTriDiagonalRepresentation();

/**********************************************************/
/********         BasisMath/ListOfPoles.c          ********/
/**********************************************************/
uint32_t InitListOfPoles();
void   FreeListOfPoles();
void   ScaleListOfPoles();
uint32_t AddListOfPoles();
uint32_t MixListOfPoles();
uint32_t SymmetrizeListOfPoles();
void   OrderListOfPoles();
void   ReduceListOfPolesWeigth();
void   ReduceListOfPolesNumberFermi();
void   ReduceListOfPolesNumber();
void   ReduceListOfPolesMoment();
void   ReduceListOfPolesSpacing();
uint32_t SubtractListOfPoles();
uint32_t GaussianBroadenListOfPoles();
uint32_t ListOfPolesCopy();
uint32_t ListOfPolesCopyVal();
uint32_t ListOfPolesCopyCon();
uint32_t ListOfPolesCopyValOccupation();
uint32_t ListOfPolesCopyConOccupation();
uint32_t ListOfPolesDetermineOccupation();
uint32_t ListOfPoleCopyChoped();
uint32_t ShiftListOfPolesToMakeChemPotZeroMetal();
uint32_t ShiftListOfPolesToMakeChemPotZeroInsulator();
uint32_t ShiftListOfPolesToMakeChemPotZero();
uint32_t ShiftListOfPolesToMakeFirstMomentZero();
uint32_t ListOfPolesCalculateChemicalPotential();

/**********************************************************/
/**********       BasisMath/MatrixVector.c       **********/
/**********************************************************/
void   FreeMatrix();
uint32_t InitMatrix();
uint32_t InitIdentityMatrix();
uint32_t MakeRealMatrixComplex();
uint32_t MakeComplexMatrixReal();
void   MatrixChop();
uint32_t MatrixTranspose();
void   MatrixComplexConjugate();
uint32_t MatrixConjugateTranspose();
double MatrixInverse();
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
void   MatrixOrthonormalize();



/**********************************************************/
/**********    BasisMath/CompactMatrixVector.c   **********/
/**********************************************************/
void   FreeCompactMatrix();
uint32_t InitCompactMatrix();
uint32_t InitCompactIdentityMatrix();
uint32_t MakeCompactRealMatrixComplex();
uint32_t CompactMatrixChop();
uint32_t CompactMatrixConjugate();
uint32_t CompactMatrixTranspose();
uint32_t CompactMatrixConjugateTranspose();
double CompactMatrixInverse();

/**********************************************************/
/*******       BasisMath/SparseMatrixVector.c       *******/
/**********************************************************/
void   FreeSparseMatrix();
uint32_t InitSparseMatrix();
uint32_t InitSparseIdentityMatrix();
uint32_t MakeRealSparseMatrixComplex();
uint32_t SparseMatrixAddElement();
void   SparseMatrixExtractElement();
void   MinimalAlocSpareseMatrix();
void   RemoveZeroSpareseMatrix();
void   SparseMatrixChop();
void   SparseMatrixTranspose();
void   SparseMatrixComplexConjugate();
void   SparseMatrixConjugateTranspose();
double SparseMatrixInverse();
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
void   ComplexSparseMatrixVector();
void   ComplexConjugateTransposeSparseMatrixVector();
uint32_t RealVectorSparseMatrixVector();
void   ComplexVectorSparseMatrixVector();
void   SparseMatrixOrthonormalize();




/**********************************************************/
/**      BasisMath/DoubleTriDiagonalMatrixVector.c       **/
/**********************************************************/
void   FreeDoubleTriDiagonalMatrix();
uint32_t InitDoubleTriDiagonalMatrix();
uint32_t InitIdentityDoubleTriDiagonalMatrix();
void   DoubleTriDiagonalMatrixChop();
double DoubleTriDiagonalMatrixInverse();
void   DoubleTriDiagonalMatrixAdd();
void   DoubleTriDiagonalMatrixScalarMultiply();
uint32_t DoubleTriDiagonalMatrixPart();
uint32_t DoubleTriDiagonalMatrixCopy();
uint32_t RealDoubleTriDiagonalMatrixVector();
uint32_t ComplexDoubleTriDiagonalMatrixVector();
uint32_t RealVectorDoubleTriDiagonalMatrixVector();
uint32_t ComplexVectorDoubleTriDiagonalMatrixVector();



/**********************************************************/
/**       BasisMath/BlokTriDiagonalMatrixVector.c        **/
/**********************************************************/
void   FreeBlokTriDiagonalMatrix();
void   FreeFrameBlokTriDiagonalMatrix();
uint32_t InitBlokTriDiagonalMatrix();
uint32_t InitFrameBlokTriDiagonalMatrix();
uint32_t InitIdentityBlokTriDiagonalMatrix();
void   BlokTriDiagonalMatrixChop();
void   BlokTriDiagonalMatrixAdd();
void   BlokTriDiagonalMatrixScalarMultiply();
uint32_t BlokTriDiagonalMatrixPart();
uint32_t BlokTriDiagonalMatrixCopy();
void   BlokTriDiagonalMatrixTranspose();
void   BlokTriDiagonalMatrixComplexConjugate();
void   BlokTriDiagonalMatrixConjugateTranspose();
uint32_t RealBlokTriDiagonalMatrixVector();
uint32_t ComplexBlokTriDiagonalMatrixVector();
uint32_t RealVectorBlokTriDiagonalMatrixVector();
uint32_t ComplexVectorBlokTriDiagonalMatrixVector();


/**********************************************************/
/**         BasisMath_BlockTriDiagonalMatrix.c           **/
/**********************************************************/
uint32_t InitBlockTridiagonalMatrix();
uint32_t FreeBlockTridiagonalMatrix();




/**********************************************************/
/****      BasisMath/BlokAndersonMatrixVector.c        ****/
/**********************************************************/
void   FreeBlokAndersonMatrix();
void   FreeFrameBlokAndersonMatrix();
uint32_t InitBlokAndersonMatrix();
uint32_t InitFrameBlokAndersonMatrix();
uint32_t InitIdentityBlokAndersonMatrix();
void   BlokAndersonMatrixChop();
void   BlokAndersonMatrixAdd();
void   BlokAndersonMatrixScalarMultiply();
uint32_t BlokAndersonMatrixPart();
uint32_t BlokAndersonMatrixCopy();
void   BlokAndersonMatrixTranspose();
void   BlokAndersonMatrixComplexConjugate();
void   BlokAndersonMatrixConjugateTranspose();
uint32_t RealBlokAndersonMatrixVector();
uint32_t ComplexBlokAndersonMatrixVector();
uint32_t RealVectorBlokAndersonMatrixVector();
uint32_t ComplexVectorBlokAndersonMatrixVector();


/**********************************************************/
/**    BasisMath/BlokDoubleTriDiagonalMatrixVector.c     **/
/**********************************************************/
void   FreeBlokDoubleTriDiagonalMatrix();
void   FreeFrameBlokDoubleTriDiagonalMatrix();
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
/********    BasisMath/MatrixTypeConversion.c      ********/
/**********************************************************/
/********    AndersonMatrixTo                      ********/
uint32_t AndersonMatrixToMatrix();
uint32_t AndersonMatrixToTriDiagonalMatrixNonOrthogonal();
uint32_t AndersonMatrixToTriDiagonalMatrix();
uint32_t AndersonMatrixToDoubleTriDiagonalMatrixNonOrthogonal();
uint32_t AndersonMatrixToDoubleTriDiagonalMatrix();
uint32_t AndersonMatrixInverseToListOfPoles();
uint32_t AndersonMatrixToEigensystem();
uint32_t AndersonMatrixToListOfPoles();
uint32_t AndersonMatrixToBADoubleTriDiagonalMatrix();
/********    MatrixTo                              ********/
uint32_t BlokAndersonMatrixToSparseMatrix();
uint32_t MatrixToEigensystem();
uint32_t MatrixToDoubleTriDiagonalMatrix();
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
uint32_t MatrixToTriDiagonalMatrix();
uint32_t TriDiagonalize();
uint32_t MatrixToEigensystem();
uint32_t Eigensystem();
uint32_t MatrixToGeneralizedEigensystem();
uint32_t GeneralizedEigensystem();
uint32_t MatrixToAndersonMatrix();
uint32_t MatrixToBlokTriDiagonalMatrix();
uint32_t BlokTriDiagonalize();
uint32_t MatrixToBlokAndersonMatrix();
uint32_t MatrixToDoubleTriDiagonalMatrix();
uint32_t MatrixToBlokDoubleTriDiagonalMatrix();
/********    ListOfPolesTo                              ********/
uint32_t ListOfPolesToBADoubleTridiagonalMatrix();
uint32_t ListOfPolesToBADoubleTridiagonalMatrixOccupation();

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
uint32_t TriDiagonalMatrixToSecondDoubleTriDiagonalMatrix();
uint32_t MapSecondoubleTridiagonalMatrixToBA();


/**********************************************************/
/***********     BasisMath/Combinatorial.c      ***********/
/**********************************************************/
double Factorial();
double Factorial2();
uint32_t BinomialCoefficient();
uint32_t PerturbationOptions();
void   PerturbationOrderListToIndex();
void   PerturbationOrderIndexToList();
double ThreeJSymbol();
double ThreeJSymbolHalf();
double SlaterCoefficientC();
