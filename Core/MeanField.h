uint32_t CompactMatrixToOneParticleOperator(const CompactMatrixType * const Omat, OperatorType * const O);
uint32_t OperatorToOneElectronMatrix(const OperatorType * const O, MatrixType * const Oone);
uint32_t OperatorCrAnLength2ToOneElectronCompactMatrix(const OperatorType * const O, CompactMatrixType * const Oone);
uint32_t MeanFieldOperatorOneElectronMatrix(const OperatorType * const O, CompactMatrixType * const MO, const CompactMatrixType Rho);
uint32_t MeanFieldOperator(const OperatorType * const O, OperatorType * const MO, const CompactMatrixType Rho);
uint32_t CalculateDensityMatrixFromOneElectronEigensystem(const CompactMatrixType * const fun, double *val, CompactMatrixType *Rho, const double kbT, const double Nelec, double *Etot);
uint32_t MeanFieldGroundstate(const OperatorType *const O, CompactMatrixType * const Rho, const double kbT, double * Etot);
