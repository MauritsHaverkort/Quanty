uint32_t StreamReadSpectraSettings(SpectraSettings * O, FILE * stream);
uint32_t FileReadSpectraSettings(SpectraSettings * O, const char *const flname);
uint32_t StreamPrintSpectraSettings(SpectraSettings O, FILE * stream);
uint32_t PrintSpectraSettings(SpectraSettings O);
uint32_t DetermineImpurityShiftFromOccupationAndersonMatrix2(AndersonMatrixType M, double const occupation, double *const shift);
uint32_t XASBathBATrunc(SpectraSettings const O, BADoubleTriDiagonalMatrixType const BathBA,BADoubleTriDiagonalMatrixType * const BathBATrunc);
uint32_t SpectraVPerturbation(SpectraSettings const O,BADoubleTriDiagonalMatrixType const GBathBAa1g,BADoubleTriDiagonalMatrixType const GBathBAb1g,BADoubleTriDiagonalMatrixType const GBathBAb2g,BADoubleTriDiagonalMatrixType const GBathBAeg,OperatorType * const VPerturbation,TriDiagonalMatrixType * const gPerturbation);
uint32_t get2p3dAtomicInfo(uint32_t Z, uint32_t dCount, double *Fdd, double *Fpd, double *Gpd, double *zeta2p, double *zeta3d);
int LigandField(uint32_t Z, uint32_t dCount);
int ImpurityModel3d(uint32_t Z, uint32_t dCount);
