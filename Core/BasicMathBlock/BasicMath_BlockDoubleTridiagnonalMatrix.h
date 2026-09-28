uint32_t InitBlockDoubleTriDiagonalMatrix(BlockDoubleTriDiagonalMatrixType *const A);
uint32_t FreeBlockDoubleTriDiagonalMatrix(BlockDoubleTriDiagonalMatrixType *const A);
uint32_t SetBlockDoubleTriDiagonalMatrixQComplex(BlockDoubleTriDiagonalMatrixType *const A);
bool BlockDoubleTriDiagonalMatrixB0IsZero(BlockDoubleTriDiagonalMatrixType const *const A);
uint32_t BlockDoubleTriDiagonalMatrixChop(BlockDoubleTriDiagonalMatrixType *const A, double const Epsilon);
uint32_t BlockDoubleTriDiagonalMatrixChopReduceDimension(BlockDoubleTriDiagonalMatrixType *const A, double const Epsilon);
uint32_t BlockDoubleTriDiagonalMatrixTruncate(BlockDoubleTriDiagonalMatrixType *const A, uint32_t const Nval, uint32_t const Ncon);
uint32_t CopyBlockDoubleTriDiagonalMatrix(BlockDoubleTriDiagonalMatrixType A, BlockDoubleTriDiagonalMatrixType * B);


