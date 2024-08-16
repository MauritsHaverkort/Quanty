uint32_t ResponseFunctionEstimateEminAndEmaxPole(const ResponsefunctionType *const M, double *EminPole, double *EmaxPole);
uint32_t ResponseFunctionToSpectra(const ResponsefunctionType *M, SpectraTypeArray *A, double GroundStateEnergy, double Emin, double Emax, double GammaG, double GammaL, uint32_t NE);
uint32_t BlockResponseFunctionToSpectra(const ResponsefunctionType *M, SpectraTypeArray *A, double GroundStateEnergy, double Emin, double Emax, double GammaG, double GammaL, uint32_t NE);
uint32_t NoneBlockResponseFunctionToSpectra(const ResponsefunctionType *M, SpectraType *A, double GroundStateEnergy, double Emin, double Emax, double GammaG, double GammaL, uint32_t NE);

