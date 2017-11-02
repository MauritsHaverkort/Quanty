

/**********************************************************/
/********     Core/WaveFunctionArithmetic.c        ********/
/**********************************************************/

uint32_t WaveFunctionInit(WaveFunctionType * Psi);
uint32_t WaveFunctionInitCopyBasis(WaveFunctionType * Psi, WaveFunctionType const PsiOrig);
void     WaveFunctionFree(WaveFunctionType * const psi);
void     WaveFunctionFreeValues(WaveFunctionType * const psi);
uint32_t PsiHermitianOperatorPsi(const OperatorType O, const WaveFunctionType psi, double *VRe);
uint32_t PsiHermitianOperatorSqrPsi(const OperatorType O, const WaveFunctionType psi, double *VRe);




/**********************************************************/
/*****  Core/WaveFunctionPerturbationArithmetic.c     *****/
/**********************************************************/

uint32_t WaveFunctionPerturbationInit(WaveFunctionPerturbationType * Psi);
uint32_t WaveFunctionPerturbationInitAll(WaveFunctionPerturbationType * Psi, uint32_t NFermion, uint32_t NBoson, uint32_t QComplex, uint32_t NBitsKey);
uint32_t QWaveFunctionPerturbationMinimumBloksAlocated(WaveFunctionPerturbationType psi);
uint32_t MakeComplexWaveFunctionPerturbationFromReal(WaveFunctionPerturbationType * Cpsi, WaveFunctionPerturbationType psi);
uint32_t MakeRealWaveFunctionPerturbationFromComplex(WaveFunctionPerturbationType Cpsi,  WaveFunctionPerturbationType * psi);
void   WaveFunctionPerturbationFreeFramework(WaveFunctionPerturbationType * const psi);
void   WaveFunctionPerturbationFree(WaveFunctionPerturbationType * const psi);
void   WaveFunctionPerturbationFreeValues(WaveFunctionPerturbationType * const psi);
uint32_t WaveFunctionPerturbationAddElement(WaveFunctionPerturbationType * const psi, const uint8_t * const det, const uint32_t order, const uint32_t iV, const double aRe, const double aIm);
uint32_t WaveFunctionPerturbationAddElementConserveBasis(WaveFunctionPerturbationType * const psi, const uint8_t * const det, const uint32_t order, const uint32_t iV, const double aRe, const double aIm);
void   WaveFunctionPerturbationPrefactorDet(const WaveFunctionPerturbationType psi, const uint8_t * const det, const uint32_t order, const uint32_t iV, double *aRe, double *aIm);
uint32_t WaveFunctionPerturbationAddSameBasis(WaveFunctionPerturbationType * const a, const WaveFunctionPerturbationType b);
uint32_t WaveFunctionPerturbationAddMultipliedSameBasis(WaveFunctionPerturbationType * const a, const WaveFunctionPerturbationType b, const double alpha);
uint32_t WaveFunctionPerturbationAddComplexMultipliedSameBasis(WaveFunctionPerturbationType * const a, const WaveFunctionPerturbationType b, const double alphaRe, const double alphaIm);
uint32_t WaveFunctionPerturbationAddConserveBasis(WaveFunctionPerturbationType * const a, const WaveFunctionPerturbationType b);
uint32_t WaveFunctionPerturbationAddMultipliedConserveBasis(WaveFunctionPerturbationType * const a, const WaveFunctionPerturbationType b, const double alpha);
uint32_t WaveFunctionPerturbationAddComplexMultipliedConserveBasis(WaveFunctionPerturbationType * const a, const WaveFunctionPerturbationType b, const double alphaRe, const double alphaIm);
uint32_t WaveFunctionPerturbationAddUnionBasis(WaveFunctionPerturbationType * const a, const WaveFunctionPerturbationType b);
uint32_t WaveFunctionPerturbationAddMultipliedUnionBasis(WaveFunctionPerturbationType * const a, const WaveFunctionPerturbationType b, const double alpha);
uint32_t WaveFunctionPerturbationAddComplexMultipliedUnionBasis(WaveFunctionPerturbationType * const a, const WaveFunctionPerturbationType b, const double alphaRe, const double alphaIm);
uint32_t WaveFunctionPerturbationMultiply(WaveFunctionPerturbationType a, const double bRe, const double bIm);
void   WaveFunctionPerturbationCalculateNorm(const WaveFunctionPerturbationType a, double *Norm);
void   WaveFunctionPerturbationNormalize(WaveFunctionPerturbationType * const a);
void   WaveFunctionPerturbationDotProductSameBasis(const WaveFunctionPerturbationType a, const WaveFunctionPerturbationType b, double *DotRe, double *DotIm);
void   WaveFunctionPerturbationDotProduct(const WaveFunctionPerturbationType a, const WaveFunctionPerturbationType b, double *DotRe, double *DotIm);
void   WaveFunctionPerturbationRandomize(WaveFunctionPerturbationType * const a);
void   WaveFunctionPerturbationTruncateBasis(WaveFunctionPerturbationType * const a, const double epsilon);
void   WaveFunctionPerturbationOrthogonalizeConserveBasis(WaveFunctionPerturbationType * const psi, const WaveFunctionPerturbationType * const OtherPsi, const uint32_t Npsi);
uint32_t WaveFunctionPerturbationOrthogonalize(WaveFunctionPerturbationType * const psi, const WaveFunctionPerturbationType * const OtherPsi, const uint32_t Npsi);
uint32_t CreateStartWaveFunctionPerturbationOnsiteEnergy(const OperatorPerturbationType O, WaveFunctionPerturbationType * psi, uint32_t order, uint32_t nFermion, uint32_t nBoson, uint32_t NBitsKey);
uint32_t WaveFunctionPerturbationAddFermion(WaveFunctionPerturbationType const psiIn, WaveFunctionPerturbationType * const psiOut, uint32_t FermionPos, uint32_t occupied);
uint32_t WaveFunctionPerturbationRemoveFermion(WaveFunctionPerturbationType const psiIn, WaveFunctionPerturbationType * const psiOut, uint32_t FermionPos);


/**********************************************************/
/********   Core/RealWaveFunctionArithmetic.c      ********/
/**********************************************************/

void   RealWaveFunctionFree(WaveFunctionType * const psi);
void   RealWaveFunctionFreeValues(WaveFunctionType * const psi);
uint32_t RealWaveFunctionAddElement(WaveFunctionType * const psi, const uint8_t * const det, const double a);
void   RealWaveFunctionAddElementConserveBasis(WaveFunctionType const *psi, const uint8_t * const det, const double a);
double RealWaveFunctionPrefactorDet(const WaveFunctionType psi, const uint8_t * const det);
uint32_t RealWaveFunctionAddElementOMPMiniFlush(WaveFunctionType * const psi, uint32_t * NMini, const uint8_t * const MiniDet, const double *const MiniValue);
void   RealWaveFunctionAddSameBasis(WaveFunctionType a, const WaveFunctionType b);
void   RealWaveFunctionAddMultipliedSameBasis(WaveFunctionType a, const WaveFunctionType b, const double alpha);
uint32_t RealWaveFunctionAddUnionBasis(WaveFunctionType * const a, const WaveFunctionType b);
uint32_t RealWaveFunctionAddMultipliedUnionBasis(WaveFunctionType * const a, const WaveFunctionType b, const double alpha);
void   RealWaveFunctionAddMultipliedConserveBasis(WaveFunctionType * const a, const WaveFunctionType b,const double alpha);
void   RealWaveFunctionNormalize(WaveFunctionType * const a);
void   RealWaveFunctionMultiply(WaveFunctionType a, const double b);
double RealWaveFunctionCalculateNorm(const WaveFunctionType a);
double RealWaveFunctionDotProductSameBasis(const WaveFunctionType a, const WaveFunctionType b);
double RealWaveFunctionDotProduct(const WaveFunctionType *a, const WaveFunctionType *b);
void   RealWaveFunctionRandomize(WaveFunctionType * const a);
void   RealWaveFunctionTruncateBasis(WaveFunctionType * const a, const double epsilon);
void   RealWaveFunctionOrthogonalizeConserveBasis(WaveFunctionType * const psi, const WaveFunctionType * const OtherPsi, const uint32_t Npsi);
uint32_t RealWaveFunctionOrthogonalize(WaveFunctionType * const psi, const WaveFunctionType * const OtherPsi, const uint32_t Npsi);

/**********************************************************/
/********   Core/ComplexWaveFunctionArithmetic.c   ********/
/**********************************************************/

void   ComplexWaveFunctionFree(WaveFunctionType * const psi);
void   ComplexWaveFunctionFreeValues(WaveFunctionType * const psi);
uint32_t ComplexWaveFunctionAddElement(WaveFunctionType * const psi, const uint8_t * const det, const double aRe, const double aIm);
void   ComplexWaveFunctionAddElementConserveBasis(WaveFunctionType psi, const uint8_t * const det, const double aRe, const double aIm);
void   ComplexWaveFunctionPrefactorDet(const WaveFunctionType psi, const uint8_t * const det, double *aRe, double *aIm);
uint32_t ComplexWaveFunctionAddElementOMPMiniFlush(WaveFunctionType * const psi, uint32_t * NMini, const uint8_t * const MiniDet, const double *const MiniValueRe, const double *const MiniValueIm);
void   ComplexWaveFunctionAddSameBasis(WaveFunctionType a, const WaveFunctionType b);
void   ComplexWaveFunctionAddMultipliedSameBasis(WaveFunctionType a, const WaveFunctionType b, const double alpha);
void   ComplexwaveFunctionAddComplexMultipliedSamebasis(WaveFunctionType a, const WaveFunctionType b, const double alphaRe, const double alphaIm);
uint32_t ComplexWaveFunctionAddUnionBasis(WaveFunctionType * const a, const WaveFunctionType b);
uint32_t ComplexWaveFunctionAddMultipliedUnionbasis(WaveFunctionType * const a, const WaveFunctionType b, const double alphaRe);
uint32_t ComplexWaveFunctionAddComplexMultipliedUnionBasis(WaveFunctionType * const a, const WaveFunctionType b, const double alphaRe, const double alphaIm);
void   ComplexWaveFunctionAddMultipliedSameBasis(WaveFunctionType a, const WaveFunctionType b, const double alpha);
void   ComplexWaveFunctionNormalize(WaveFunctionType * const a);
void   ComplexWaveFunctionMultiply(WaveFunctionType a, const double b);
void   ComplexWaveFunctionComplexMultiply(WaveFunctionType * a, const double bRe, const double bIm);
double ComplexWaveFunctionCalculateNorm(const WaveFunctionType a);
void   ComplexWaveFunctionDotProductSameBasis(const WaveFunctionType a, const WaveFunctionType b, double *DotRe, double *DotIm);
void   ComplexWaveFunctionDotProduct(const WaveFunctionType *a, const WaveFunctionType *b, double *DotRe, double *DotIm);
void   ComplexWaveFunctionRandomize(WaveFunctionType * const a);
void   ComplexWaveFunctionTruncatebasis(WaveFunctionType * const a, const double epsilon);
uint32_t ComplexWaveFunctionOrthogonalize(WaveFunctionType * const psi, const WaveFunctionType * const OtherPsi, const uint32_t Npsi);
uint32_t MakeComplexWaveFunctionFromReal(WaveFunctionType * Cpsi, WaveFunctionType psi);
uint32_t MakeRealWaveFunctionFromComplex(WaveFunctionType Cpsi, WaveFunctionType * psi);



/**********************************************************/
/*************       Core/Lanczos.c          **************/
/**********************************************************/

uint32_t LanczosGroundStateRR(const OperatorType H, WaveFunctionType * const psi, double *epsilon, const double zero, double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi);
uint32_t LanczosGroundStateRC(const OperatorType H, WaveFunctionType * const psi, double *epsilon, const double zero, double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi);
uint32_t LanczosGroundStateCR(const OperatorType H, WaveFunctionType * const psi, double *epsilon, const double zero, double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi);
uint32_t LanczosGroundStateCC(const OperatorType H, WaveFunctionType * const psi, double *epsilon, const double zero, double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi);
uint32_t LanczosGroundStateMR(const OperatorType H, WaveFunctionType * const psi, double *epsilon, const double zero, double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi);
uint32_t LanczosGroundStateMC(const OperatorType H, WaveFunctionType * const psi, double *epsilon, const double zero, double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi);
uint32_t LanczosGroundStateRestrictedRR(const OperatorType H, WaveFunctionType * const psi, DeterminantRestrictionType const Restrictions, double *epsilon, const double zero, double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi);
uint32_t LanczosGroundStateRestrictedRC(const OperatorType H, WaveFunctionType * const psi, DeterminantRestrictionType const Restrictions, double *epsilon, const double zero, double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi);
uint32_t LanczosGroundStateRestrictedCR(const OperatorType H, WaveFunctionType * const psi, DeterminantRestrictionType const Restrictions, double *epsilon, const double zero, double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi);
uint32_t LanczosGroundStateRestrictedCC(const OperatorType H, WaveFunctionType * const psi, DeterminantRestrictionType const Restrictions, double *epsilon, const double zero, double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi);
uint32_t LanczosGroundStateRestrictedMR(const OperatorType H, WaveFunctionType * const psi, DeterminantRestrictionType const Restrictions, double *epsilon, const double zero, double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi);
uint32_t LanczosGroundStateRestrictedMC(const OperatorType H, WaveFunctionType * const psi, DeterminantRestrictionType const Restrictions, double *epsilon, const double zero, double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi);
uint32_t LanczosTriDiagonalizeRR(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType * const M, double Gamma);
uint32_t LanczosTriDiagonalizeCR(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType * const M, double Gamma);
uint32_t LanczosTriDiagonalizeMR(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType * const M, double Gamma);
uint32_t LanczosTriDiagonalizeRC(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType * const M, double Gamma);
uint32_t LanczosTriDiagonalizeCC(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType * const M, double Gamma);
uint32_t LanczosTriDiagonalizeMC(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType * const M, double Gamma);
uint32_t LanczosTriDiagonalizePlusDensityRR(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType * const M, MatrixType * const Rho, double Gamma);
uint32_t LanczosTriDiagonalizePlusDensityCR(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType * const M, MatrixType * const Rho, double Gamma);
uint32_t LanczosTriDiagonalizePlusDensityMR(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType * const M, MatrixType * const Rho, double Gamma);
uint32_t LanczosTriDiagonalizePlusDensityRC(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType * const M, MatrixType * const Rho, double Gamma);
uint32_t LanczosTriDiagonalizePlusDensityCC(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType * const M, MatrixType * const Rho, double Gamma);
uint32_t LanczosTriDiagonalizePlusDensityMC(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType * const M, MatrixType * const Rho, double Gamma);


uint32_t LanczosGroundState(const OperatorType H, WaveFunctionType * const psi, double *epsilon, const double zero, double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi);
uint32_t LanczosGroundStateRestricted(const OperatorType H, WaveFunctionType * const psi, DeterminantRestrictionType const Restrictions, double *epsilon, const double zero, double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi);
uint32_t LanczosTriDiagonalize(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType * const M, double Gamma);
uint32_t LanczosTriDiagonalizePlusDensity(const OperatorType H, WaveFunctionType psi, double epsilon, TriDiagonalMatrixType * const M, MatrixType * const Rho, double Gamma);



/**********************************************************/
/********       Core/OperatorArithmetic.c          ********/
/**********************************************************/

uint32_t FreeOperator(OperatorType * const H);
uint32_t InitOperator(OperatorType * const O);
uint32_t InitNamedOperator(OperatorType * const O, char *name);
uint32_t MakeRealOperatorFixedLengthComplex(OperatorType * const O, const uint32_t length);
uint32_t InitOperatorFixedLength(OperatorType * const O, const uint32_t length, const uint32_t NBitsKey);
uint32_t InitOperatorCopySettings(OperatorType * const O, const OperatorType Oset);
uint32_t RealOperatorAddLader(OperatorType * const O, const uint16_t * const oadd, const uint32_t length, const double VRe);
uint32_t ComplexOperatorAddLader(OperatorType * const O, const uint16_t * const oadd, const uint32_t length, const double VRe, const double VIm);
uint32_t RealOperatorAddLaderNormalOrder(OperatorType * const O, const uint16_t * const oadd, const uint32_t length, const double VRe, uint32_t(*QDetNotOrdered) (uint16_t, uint16_t, uint32_t));
uint32_t ComplexOperatorAddLaderNormalOrder(OperatorType * const O, const uint16_t * const oadd, const uint32_t length, const double VRe, const double VIm, uint32_t(*QDetNotOrdered) (uint16_t,uint16_t,uint32_t));
uint32_t OperatorAddLader(OperatorType * const O, const uint16_t * const oadd, const uint32_t length, const double VRe,const double VIm);
uint32_t OperatorAddLaderNormalOrder(OperatorType * const O, const uint16_t * const oadd, const uint32_t length,const double VRe, const double VIm, uint32_t(*QDetNotOrdered) (uint16_t, uint16_t, uint32_t));
uint32_t OperatorRemoveSmallValues(OperatorType * const O, const double epsilon);
uint32_t QDetNotOrderedFaFcBaBc(const uint16_t A, const uint16_t B, const uint32_t NFermion);
uint32_t QDetNotOrderedFcFaBcBa(const uint16_t A, const uint16_t B, const uint32_t NFermion);
uint32_t QDetNotOrderedFcFaBaBc(const uint16_t A, const uint16_t B, const uint32_t NFermion);
uint32_t QDetNotOrderedFaFcBcBa(const uint16_t A, const uint16_t B, const uint32_t NFermion);
uint32_t NormalOrderOperator(const OperatorType O, OperatorType * const Oo, uint32_t(*QDetNotOrdered) (uint16_t, uint16_t, uint32_t));
uint32_t OperatorMultiply(OperatorType * const O, OperatorType const B, OperatorType const A);
uint32_t OperatorScale(OperatorType * const O, const double alphaRe, const double alphaIm);
uint32_t OperatorAdd(OperatorType * const OA, const OperatorType OB, const double alphaRe, const double alphaIm);
//uint32_t MappedOperatorAddAtLevel();
uint32_t MappedOperatorAdd(OperatorType * const OA, const OperatorType OB, const MatrixType U);

/**********************************************************/
/********       Core/OperatorWaveFunction.c        ********/
/**********************************************************/

double OperatorCr(uint8_t * const Determinant, const uint16_t i, const uint32_t NFermion);
double OperatorAn(uint8_t * const Determinant, const uint16_t i, const uint32_t NFermion);
double OperatorLader(uint8_t * const Determinant, const uint16_t * const o, const uint32_t lenght, const uint32_t NFermion);
uint32_t OperatorPsiConserveBasisCC(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi);
uint32_t OperatorPsiConserveBasisCR(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi);
uint32_t OperatorPsiConserveBasisRC(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi);
uint32_t OperatorPsiConserveBasisRR(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi);
uint32_t OperatorPsiConserveBasisMC(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi);
uint32_t OperatorPsiConserveBasisMR(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi);
uint32_t OperatorPsiConserveBasisPreAllocatedCC(const OperatorType *O, const WaveFunctionType *psi, WaveFunctionType * const Opsi);
uint32_t OperatorPsiConserveBasisPreAllocatedCR(const OperatorType *O, const WaveFunctionType *psi, WaveFunctionType * const Opsi);
uint32_t OperatorPsiConserveBasisPreAllocatedRC(const OperatorType *O, const WaveFunctionType *psi, WaveFunctionType * const Opsi);
uint32_t OperatorPsiConserveBasisPreAllocatedRR(const OperatorType *O, const WaveFunctionType *psi,WaveFunctionType * const Opsi);
uint32_t OperatorPsiConserveBasisPreAllocatedMC(const OperatorType *O, const WaveFunctionType *psi, WaveFunctionType * const Opsi);
uint32_t OperatorPsiConserveBasisPreAllocatedMR(const OperatorType *O, const WaveFunctionType *psi, WaveFunctionType * const Opsi);
uint32_t OperatorPsiCC(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi);
uint32_t OperatorPsiCR(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi);
uint32_t OperatorPsiRC(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi);
uint32_t OperatorPsiRR(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi);
uint32_t OperatorPsiMC(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi);
uint32_t OperatorPsiMR(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi);
//uint32_t OperatorPsiSerialCC();
//uint32_t OperatorPsiSerialCR();
//uint32_t OperatorPsiSerialRC();
//uint32_t OperatorPsiSerialRR();
//uint32_t OperatorPsiSerialMC();
//uint32_t OperatorPsiSerialMR();
uint32_t PsiOperatorPsiCC(const OperatorType O, const WaveFunctionType psi, double *VRe, double *VIm);
uint32_t PsiOperatorPsiCR(const OperatorType O, const WaveFunctionType psi, double *VRe, double *VIm);
uint32_t PsiOperatorPsiRC(const OperatorType O, const WaveFunctionType psi, double *VRe, double *VIm);
uint32_t PsiOperatorPsiRR(const OperatorType O, const WaveFunctionType psi, double *VRe);
uint32_t PsiOperatorPsiMC(const OperatorType O, const WaveFunctionType psi, double *VRe, double *VIm);
uint32_t PsiOperatorPsiMR(const OperatorType O, const WaveFunctionType psi, double *VRe, double *VIm);
uint32_t OperatorPsiConserveBasis(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi);
uint32_t OperatorPsiConserveBasisPreAllocated(const OperatorType *O, const WaveFunctionType *psi, WaveFunctionType * const Opsi);
uint32_t OperatorPsi(const OperatorType O, const WaveFunctionType psi, WaveFunctionType * const Opsi);
//uint32_t OperatorPsiSerial();
uint32_t PsiOperatorPsi(const OperatorType O, const WaveFunctionType psi, double *VRe, double *VIm);
uint32_t QDetNotRestricted(uint8_t const *det, uint32_t const NFermion, uint32_t const NBoson, uint32_t const NWord, DeterminantRestrictionType const Restrictions);



/**********************************************************/
/*****    Core/OperatorWaveFunctionPerturbation.c     *****/
/**********************************************************/

uint32_t OperatorPsiPerturbationConserveBasis(const OperatorPerturbationType O, const WaveFunctionPerturbationType psi, WaveFunctionPerturbationType * const Opsi);
uint32_t OperatorPsiPerturbationConserveBasisPreAllocated(const OperatorPerturbationType O, const WaveFunctionPerturbationType psi, WaveFunctionPerturbationType * const Opsi);
uint32_t OperatorPsiPerturbation(const OperatorPerturbationType O, const WaveFunctionPerturbationType psi, WaveFunctionPerturbationType * const Opsi);
uint32_t PsiOperatorPsiPerturbation(const OperatorPerturbationType O, const WaveFunctionPerturbationType psi, double *VRe, double *VeIm);
uint32_t PsiHermitianOperatorPsiPerturbation(const OperatorPerturbationType O, const WaveFunctionPerturbationType psi, double *VRe);
uint32_t PsiOperatorSqrPsiPerturbation(const OperatorPerturbationType O, const WaveFunctionPerturbationType psi, double *VRe, double *VIm);
uint32_t PsiHermitianOperatorSqrPsiPerturbation(const OperatorPerturbationType O, const WaveFunctionPerturbationType psi, double *VRe);
uint32_t PsiPrimeOperatorPsiPerturbationMatrix(const OperatorPerturbationType O, const WaveFunctionPerturbationType * psiPrime, WaveFunctionPerturbationType * const psi, MatrixType * const V);
uint32_t PsiOperatorPsiPerturbationMatrix(const OperatorPerturbationType O, const WaveFunctionPerturbationType * psi, MatrixType * const V);
uint32_t PsiPrimeOperatorPsiPerturbation(const OperatorPerturbationType O, const WaveFunctionPerturbationType psiPrime, WaveFunctionPerturbationType const psi, double *VRe, double *VIm);




/**********************************************************/
/**************        Core/Spectra.c        **************/
/**********************************************************/
void   CalculateMoment(SpectraType * const G);
int SpectrumAllocate(const SpectraType Settings, SpectraType * const G);
void SpectrumFree(SpectraType * const G);



/**********************************************************/
/**************     Core/DensityMatrix.c     **************/
/**********************************************************/
uint32_t CalculateDensity(double *const Rho, const WaveFunctionType psi);
uint32_t CalculateDensityMatrix(MatrixType * const Rho, const WaveFunctionType psi);
uint32_t CalculatePrefactorHistogram(double *const Hist, const WaveFunctionType psi, const double exp, const uint32_t N);
uint32_t CalculatePartialDensityCorrelationTensor(Complex16 * Rho, const WaveFunctionType psi, uint32_t const NFermionIncluded, uint16_t const * FermionIncluded);



double InterpolateRootOfset(double const *x, double const *y, uint32_t N, double const ofset);
double InterpolateRoot(double const *x, double const *y, uint32_t N);
