uint32_t InitBlockDoubleTriDiagonalMatrix(BlockDoubleTriDiagonalMatrixType *const A);
uint32_t FreeBlockDoubleTriDiagonalMatrix(BlockDoubleTriDiagonalMatrixType *const A);
uint32_t BlockDoubleTriDiagonalMatrixChop(BlockDoubleTriDiagonalMatrixType *const A, double const Epsilon);
uint32_t BlockDoubleTriDiagonalMatrixChopReduceDimension(BlockDoubleTriDiagonalMatrixType *const A, double const Epsilon);
uint32_t CopyBlockDoubleTriDiagonalMatrix(BlockDoubleTriDiagonalMatrixType A, BlockDoubleTriDiagonalMatrixType * B);


