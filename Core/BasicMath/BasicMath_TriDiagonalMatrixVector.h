
uint32_t TriDiagonalMatrixToGreensFunction(TriDiagonalMatrixType *const SpecTri,double wRes,double energy,double gamma,MatrixType *const invMat);
uint32_t TriDiagonalMatrixBProduct(TriDiagonalMatrixType *SpecTri,uint32_t blow,uint32_t bhigh,double *bProduct);
void ComplexVectorTriDiagonalMatrixVector(TriDiagonalMatrixType const M,const double *const xRe,const double *const xIm,double *const yRe,double *const yIm);
void RealVectorTriDiagonalMatrixVector(TriDiagonalMatrixType const M,const double *const x,double *const y);
void ComplexTriDiagonalMatrixVector(TriDiagonalMatrixType const M,const double *const xRe,const double *const xIm,double *const yRe,double *const yIm);
void RealTriDiagonalMatrixVector(TriDiagonalMatrixType const M,const double *const x,double *const y);
uint32_t TriDiagonalMatrixCopy(TriDiagonalMatrixType *const M,TriDiagonalMatrixType A);
uint32_t TriDiagonalMatrixPart(TriDiagonalMatrixType *const M,TriDiagonalMatrixType A,const uint32_t iStart,const uint32_t iEnd);
void TriDiagonalMatrixScalarMultiply(TriDiagonalMatrixType *const M,double alphaRe);
void TriDiagonalMatrixAdd(TriDiagonalMatrixType *const M,TriDiagonalMatrixType B,double alphaRe);
double MatrixInverse(MatrixType *M);
uint32_t InitMatrix(MatrixType *const M);
double TriDiagonalMatrixInverse(TriDiagonalMatrixType const M,MatrixType *const MI);
void ShiftTriDiagonalMatrix(TriDiagonalMatrixType *const M,double const delta);
void TriDiagonalMatrixChop(TriDiagonalMatrixType *const M,double const epsilon);
uint32_t InitIdentityTriDiagonalMatrix(TriDiagonalMatrixType *const M);
uint32_t InitTriDiagonalMatrix(TriDiagonalMatrixType *const M);
void FreeTriDiagonalMatrix(TriDiagonalMatrixType *const M);
