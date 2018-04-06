/* This file was automatically generated.  Do not edit! */
uint32_t CreateMultiSitesGclDMFT(MultiSitesDMFTSettings *const O,OperatorType const HAnd,WaveFunctionType const Psi,OperatorType const *VPerturbation,TriDiagonalMatrixType const *gPerturbation,double const mu,AndersonMatrixType *const GclupAnd,AndersonMatrixType *const GcldnAnd,TimeType *const Time);
uint32_t MultiSitesBADoubleTriDiagonalMatrixToAndersonImpurityOperator(OperatorType *O,OperatorType *Oloc,BADoubleTriDiagonalMatrixType *Mup,BADoubleTriDiagonalMatrixType *Mdn,uint32_t NSite);
uint32_t DCAHamiltonian(OperatorType *O,BADoubleTriDiagonalMatrixType *Mup,BADoubleTriDiagonalMatrixType *Mdn,uint32_t NSite,double U,double mu);
uint32_t MultiSitesTestPsiDimensions(MultiSitesDMFTSettings *const O,BADoubleTriDiagonalMatrixType *const GBathupBATrunc,BADoubleTriDiagonalMatrixType *const GBathdnBATrunc,WaveFunctionType *const Psi);
uint32_t MultiSitesCreateGBathBATruncandVPerturbation(MultiSitesDMFTSettings *const O,BADoubleTriDiagonalMatrixType *const GBathupBA,BADoubleTriDiagonalMatrixType *const GBathdnBA,BADoubleTriDiagonalMatrixType *const GBathupBATrunc,BADoubleTriDiagonalMatrixType *const GBathdnBATrunc,OperatorType *const VPerturbation,TriDiagonalMatrixType *const gPerturbation);
uint32_t MultiSitesFindGroundStateFromG0andSigma(MultiSitesDMFTSettings *const O,AndersonMatrixType *const G0,ListOfPoles *const Sigmaup,ListOfPoles *const Sigmadn,double *const mu,double *const Nelectron,WaveFunctionType *const Psi,BADoubleTriDiagonalMatrixType *const GBathupBARet,BADoubleTriDiagonalMatrixType *const GBathdnBARet,AndersonMatrixType *const GBathupAndRet,AndersonMatrixType *const GBathdnAndRet,TimeType *const Time,double *const onsiteshiftstart);
double read_oneline(const char *filename);
uint32_t initializeMultiSitesDMFT(const char *const flname,MultiSitesDMFTSettings *const O,AndersonMatrixType *const G0And,ListOfPoles *const SigmaupLP,ListOfPoles *const SigmadnLP,WaveFunctionType *const PsiStart);
uint32_t FileReadMultiSitesDMFTSettings(MultiSitesDMFTSettings *O,const char *const flname);
int main();
int main(int argc,const char *argv[]);
int main(int argc,const char *argv[]);
