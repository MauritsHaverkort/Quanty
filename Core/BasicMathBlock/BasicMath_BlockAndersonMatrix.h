uint32_t InitBlockAndersonMatrix(BlockAndersonMatrixType *const A);
uint32_t FreeBlockAndersonMatrix(BlockAndersonMatrixType *A);
uint32_t ComplexBlockAndersonMatrixToReal(BlockAndersonMatrixType *const A);
uint32_t RealBlockAndersonMatrixToComplex(BlockAndersonMatrixType *const A);
uint32_t BlockAndersonMatrixChop(BlockAndersonMatrixType *const A, double const Epsilon);
uint32_t BlockAndersonMatrixChopReduceDimension(BlockAndersonMatrixType *const A, double const Epsilon);
uint32_t CopyBlockAndersonMatrix(BlockAndersonMatrixType A, BlockAndersonMatrixType *B);
int BlockAndersonMatrixInverseToListOfPoles(BlockAndersonMatrixType const * const A, BlockListOfPoles * const B, double epsilon);
int BlockAndersonMatrixFromListOfPolesAsSelfEnergy(BlockListOfPoles const * const Sigma, BlockAndersonMatrixType * const G);
