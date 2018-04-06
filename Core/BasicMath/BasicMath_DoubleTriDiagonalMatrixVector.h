/* This file was automatically generated.  Do not edit! */
uint32_t ComplexVectorDoubleTriDiagonalMatrixVector(DoubleTriDiagonalMatrixType const M,const double *const xRe,const double *const xIm,double *const yRe,double *const yIm);
uint32_t RealVectorDoubleTriDiagonalMatrixVector(DoubleTriDiagonalMatrixType const M,const double *const x,double *const y);
uint32_t ComplexDoubleTriDiagonalMatrixVector(DoubleTriDiagonalMatrixType const M,const double *const xRe,const double *const xIm,double *const yRe,double *const yIm);
uint32_t RealDoubleTriDiagonalMatrixVector(DoubleTriDiagonalMatrixType const M,const double *const x,double *const y);
uint32_t DoubleTriDiagonalMatrixCopy(DoubleTriDiagonalMatrixType *const M,DoubleTriDiagonalMatrixType A);
uint32_t DoubleTriDiagonalMatrixPart(DoubleTriDiagonalMatrixType *const M,DoubleTriDiagonalMatrixType A,const uint32_t iEndval,const uint32_t iEndcon);
void DoubleTriDiagonalMatrixScalarMultiply(DoubleTriDiagonalMatrixType *const M,double alphaRe);
void DoubleTriDiagonalMatrixAdd(DoubleTriDiagonalMatrixType *const M,DoubleTriDiagonalMatrixType B,double alphaRe);
double MatrixInverse(MatrixType *M);
uint32_t InitMatrix(MatrixType *const M);
double DoubleTriDiagonalMatrixInverse(DoubleTriDiagonalMatrixType const M,MatrixType *const MI);
void DoubleTriDiagonalMatrixChop(DoubleTriDiagonalMatrixType *const M,double const epsilon);
uint32_t InitIdentityTriDiagonalMatrix(TriDiagonalMatrixType *const M);
uint32_t InitIdentityDoubleTriDiagonalMatrix(DoubleTriDiagonalMatrixType *const M);
uint32_t InitTriDiagonalMatrix(TriDiagonalMatrixType *const M);
uint32_t InitDoubleTriDiagonalMatrix(DoubleTriDiagonalMatrixType *const M);
void FreeTriDiagonalMatrix(TriDiagonalMatrixType *const M);
void FreeDoubleTriDiagonalMatrix(DoubleTriDiagonalMatrixType *const M);
