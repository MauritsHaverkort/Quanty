/* This file was automatically generated.  Do not edit! */
uint32_t NiImpurity(SpectraSettings SpecSet);
void NiLigandField();
uint32_t LNOInitXASWavefunctionRandom(WaveFunctionType *const psi,uint32_t const NFermions,uint32_t const NBath);
uint32_t InitXASWavefunctionRandom2(WaveFunctionType *const psi,uint32_t const NFermions,uint32_t const ncore,uint32_t const nval,uint32_t const nlig);
int twiddle(int *x,int *y,int *z,int *p);
void inittwiddle(int m,int n,int *p);
uint32_t InitXASWavefunctionRandom(WaveFunctionType *const psi,uint32_t const NFermions,uint32_t const ncore,uint32_t const nval);
uint32_t SpectraVPerturbation(SpectraSettings const O,BADoubleTriDiagonalMatrixType const GBathBAa1g,BADoubleTriDiagonalMatrixType const GBathBAb1g,BADoubleTriDiagonalMatrixType const GBathBAb2g,BADoubleTriDiagonalMatrixType const GBathBAeg,OperatorType *const VPerturbation,TriDiagonalMatrixType *const gPerturbation);
uint32_t XASBathBATrunc(SpectraSettings const O,BADoubleTriDiagonalMatrixType const BathBA,BADoubleTriDiagonalMatrixType *const BathBATrunc);
uint32_t DetermineImpurityShiftFromOccupationAndersonMatrix2(AndersonMatrixType M,double const occupation,double *const shift);
uint32_t PrintSpectraSettings(SpectraSettings O);
uint32_t StreamPrintSpectraSettings(SpectraSettings O,FILE *stream);
uint32_t FileReadSpectraSettings(SpectraSettings *O,const char *const flname);
uint32_t StreamReadSpectraSettings(SpectraSettings *O,FILE *stream);
