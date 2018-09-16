uint32_t StreamReadModelSettings(ModelSettings * O, FILE * stream);
int32_t InitModelSettings(ModelSettings * O);
uint32_t FileReadModelSettings(ModelSettings * O, const char *const flname);
uint32_t StreamPrintModelSettings(ModelSettings O, FILE * stream);
uint32_t PrintModelSettings(ModelSettings O, const char *const flname);
uint32_t InitNaturalImpurityWavefunctionRandom(WaveFunctionType * const psi, uint32_t NBitsKey, uint32_t const NFermions, uint32_t NFermionShift, uint32_t NImpOrb, uint32_t NBathV,uint32_t minNIm,uint32_t maxNIm);
uint32_t ImpurityHFromModelBath(OperatorType *const Himp, int geometry, int shape, double width, int NPoles,int NBathSitesAnd, int NBathSitesV, int NBathSitesC, double mu, double kBT, int blocksize, int spin, double *eshiftsImpurity, double *eshifts, double *hopping, int NFermion, int NFermionShift);
uint32_t ImpurityHFromModelBathPerturbation(OperatorType *const Himp, OperatorType * const VPert, TriDiagonalMatrixType * const gPert, int NCellV, int geometry, int shape, double width, int NPoles,int NBathSitesAnd, int NBathSitesV, int NBathSitesC, double mu, double kBT, int blocksize, int spin, double *eshiftsImpurity, double *eshifts, double *hopping, int NFermion, int NFermionShift);
uint32_t ImpurityHFromTB(OperatorType *const Himp, ModelSettings O);
uint32_t GetIntraShellAtomicInteractions(uint32_t Z, uint32_t n, uint32_t ell, uint32_t vCount, double Fvv[MAXNSHELLS][MAXCOULOMBVV],uint32_t iShell);
uint32_t GetSpinOrbitEnergy(uint32_t Z, uint32_t n, uint32_t ell, uint32_t vCount, double *zeta);
uint32_t GetInterShellAtomicCoulombInteractions(uint32_t Z, uint32_t nc, uint32_t lc, uint32_t nv, uint32_t lv, uint32_t vCount, double Fcv[MAXNSHELLS][MAXNSHELLS][MAXCOULOMBCV],uint32_t iShellc,uint32_t iShellv);
uint32_t GetAtomicInteractions(ModelSettings *O);
 uint32_t ScaleAtomicInteractions(ModelSettings *O);
uint32_t RunModelCalculation(ModelSettings MSet);
