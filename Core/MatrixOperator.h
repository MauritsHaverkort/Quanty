int DoubleTriDiagonalMatrixToImpurityHamiltonian(OperatorType * O, DoubleTriDiagonalMatrixType * M, bool AddSpin, uint32_t FermionShift, uint32_t NF);



uint32_t TriDiagonalMatrixToAndersonImpurityOperator(OperatorType * O, TriDiagonalMatrixType Mup, TriDiagonalMatrixType Mdn, double U, double mu);
uint32_t AndersonMatrixToAndersonImpurityOperator(OperatorType * O, AndersonMatrixType Mup, AndersonMatrixType Mdn,double U, double mu);
uint32_t SecondDoubleTriDiagonalMatrixToAndersonImpurityOperator(OperatorType * O,SecondDoubleTriDiagonalMatrixType Mup,SecondDoubleTriDiagonalMatrixType Mdn, double U,double mu);
uint32_t DoubleTriDiagonalMatrixToAndersonImpurityOperator(OperatorType * O, DoubleTriDiagonalMatrixType Mup,DoubleTriDiagonalMatrixType Mdn, double U, double mu);
uint32_t BADoubleTriDiagonalMatrixToAndersonImpurityOperator(OperatorType * O, BADoubleTriDiagonalMatrixType Mup, BADoubleTriDiagonalMatrixType Mdn, double U, double mu);
uint32_t OperatorToDeterminantMatrix(const OperatorType * const O, const WaveFunctionType * const psi, CompactMatrixType *M);
uint32_t OperatorToDeterminantMatrix2(const OperatorType * const O, const WaveFunctionType * const psi_i, const WaveFunctionType * const psi_j, CompactMatrixType *M);
uint32_t OperatorToDeterminantMatrixWithPrefactorsInPsi(const OperatorType * const O, const WaveFunctionType * const psi, CompactMatrixType *M);
