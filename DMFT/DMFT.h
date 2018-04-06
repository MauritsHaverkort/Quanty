/* This file was automatically generated.  Do not edit! */
uint32_t DMFTSpec(DMFTSettings const HSet,uint32_t const iOuter,double const Egrd,WaveFunctionType const psiGrd,OperatorType HAnd,OperatorType const *VPerturbation,TriDiagonalMatrixType const *gPerturbation,OperatorType const OAndn,OperatorType const OAnup,OperatorType const OCrdn,OperatorType const OCrup);
uint32_t CreateIPESupOperatorDMFT(OperatorType *O,uint32_t NFermion);
uint32_t CreateIPESdnOperatorDMFT(OperatorType *O,uint32_t NFermion);
uint32_t CreatePESupOperatorDMFT(OperatorType *O,uint32_t NFermion);
uint32_t CreatePESdnOperatorDMFT(OperatorType *O,uint32_t NFermion);
uint32_t CreateGclDMFT(DMFTSettings const O,OperatorType const HAnd,WaveFunctionType const Psi,OperatorType const *VPerturbation,TriDiagonalMatrixType const *gPerturbation,double const mu,AndersonMatrixType *const GclAnd,TimeType *const Time);
uint32_t TestIfPsiHasCorrectDimensions(BADoubleTriDiagonalMatrixType const GBathBATrunc,WaveFunctionType *const Psi);
uint32_t createGBathBATruncandVPerturbation(DMFTSettings const O,BADoubleTriDiagonalMatrixType const GBathBA,BADoubleTriDiagonalMatrixType *const GBathBATrunc,OperatorType *const VPerturbation,TriDiagonalMatrixType *const gPerturbation);
uint32_t FindGroundStateFromG0andSigma(DMFTSettings const O,AndersonMatrixType const G0,ListOfPoles const Sigma,double *const mu,double *const Nelectron,WaveFunctionType *const Psi,BADoubleTriDiagonalMatrixType *const GBathBARet,AndersonMatrixType *const GBathAndRet,TimeType *const Time,double *const onsiteshiftstart);
uint32_t initializeDMFT(const char *const flname,DMFTSettings *O,AndersonMatrixType *const G0And,ListOfPoles *const SigmaLP,WaveFunctionType *const PsiStart);
int main(int argc,const char *argv[]);
int main(int argc,const char *argv[]);
