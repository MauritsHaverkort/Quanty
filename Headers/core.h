

/**********************************************************/
/********     Core/WaveFunctionArithmetic.c        ********/
/**********************************************************/

uint32_t WaveFunctionInit(WaveFunctionType * Psi);
void   WaveFunctionFree(WaveFunctionType * const psi);
uint32_t PsiHermitianOperatorPsi(const OperatorType O, const WaveFunctionType psi, double *VRe);
uint32_t PsiHermitianOperatorSqrPsi(const OperatorType O, const WaveFunctionType psi, double *VRe);




/**********************************************************/
/*****  Core/WaveFunctionPerturbationArithmetic.c     *****/
/**********************************************************/

uint32_t WaveFunctionPerturbationInit(WaveFunctionPerturbationType * Psi);
uint32_t WaveFunctionPerturbationInitAll(WaveFunctionPerturbationType * Psi, uint32_t NFermion, uint32_t NBoson,
                                         uint32_t QComplex, uint32_t NBitsKey);
uint32_t QWaveFunctionPerturbationMinimumBloksAlocated(WaveFunctionPerturbationType psi);
uint32_t MakeComplexWaveFunctionPerturbationFromReal(WaveFunctionPerturbationType * Cpsi,
                                                     WaveFunctionPerturbationType psi);
uint32_t MakeRealWaveFunctionPerturbationFromComplex(WaveFunctionPerturbationType Cpsi,
                                                     WaveFunctionPerturbationType * psi);
void   WaveFunctionPerturbationFreeFramework(WaveFunctionPerturbationType * const psi);
void   WaveFunctionPerturbationFree(WaveFunctionPerturbationType * const psi);
void   WaveFunctionPerturbationFreeValues(WaveFunctionPerturbationType * const psi);
uint32_t WaveFunctionPerturbationAddElement(WaveFunctionPerturbationType * const psi, const uint8_t * const det,
                                            const uint32_t order, const uint32_t iV, const double aRe,
                                            const double aIm);
uint32_t WaveFunctionPerturbationAddElementConserveBasis(WaveFunctionPerturbationType * const psi,
                                                         const uint8_t * const det, const uint32_t order,
                                                         const uint32_t iV, const double aRe, const double aIm);
void   WaveFunctionPerturbationPrefactorDet(const WaveFunctionPerturbationType psi, const uint8_t * const det,
                                            const uint32_t order, const uint32_t iV, double *aRe, double *aIm);
uint32_t WaveFunctionPerturbationAddSameBasis(WaveFunctionPerturbationType * const a,
                                              const WaveFunctionPerturbationType b);
uint32_t WaveFunctionPerturbationAddMultipliedSameBasis(WaveFunctionPerturbationType * const a,
                                                        const WaveFunctionPerturbationType b, const double alpha);
uint32_t WaveFunctionPerturbationAddComplexMultipliedSameBasis(WaveFunctionPerturbationType * const a,
                                                               const WaveFunctionPerturbationType b,
                                                               const double alphaRe, const double alphaIm);
uint32_t WaveFunctionPerturbationAddConserveBasis(WaveFunctionPerturbationType * const a,
                                                  const WaveFunctionPerturbationType b);
uint32_t WaveFunctionPerturbationAddMultipliedConserveBasis(WaveFunctionPerturbationType * const a,
                                                            const WaveFunctionPerturbationType b, const double alpha);
uint32_t WaveFunctionPerturbationAddComplexMultipliedConserveBasis(WaveFunctionPerturbationType * const a,
                                                                   const WaveFunctionPerturbationType b,
                                                                   const double alphaRe, const double alphaIm);
uint32_t WaveFunctionPerturbationAddUnionBasis(WaveFunctionPerturbationType * const a,
                                               const WaveFunctionPerturbationType b);
uint32_t WaveFunctionPerturbationAddMultipliedUnionBasis(WaveFunctionPerturbationType * const a,
                                                         const WaveFunctionPerturbationType b, const double alpha);
uint32_t WaveFunctionPerturbationAddComplexMultipliedUnionBasis(WaveFunctionPerturbationType * const a,
                                                                const WaveFunctionPerturbationType b,
                                                                const double alphaRe, const double alphaIm);
uint32_t WaveFunctionPerturbationMultiply(WaveFunctionPerturbationType a, const double bRe, const double bIm);
void   WaveFunctionPerturbationCalculateNorm(const WaveFunctionPerturbationType a, double *Norm);
void   WaveFunctionPerturbationNormalize(WaveFunctionPerturbationType * const a);
void   WaveFunctionPerturbationDotProductSameBasis(const WaveFunctionPerturbationType a,
                                                   const WaveFunctionPerturbationType b, double *DotRe, double *DotIm);
void   WaveFunctionPerturbationDotProduct(const WaveFunctionPerturbationType a, const WaveFunctionPerturbationType b,
                                          double *DotRe, double *DotIm);
void   WaveFunctionPerturbationRandomize(WaveFunctionPerturbationType * const a);
void   WaveFunctionPerturbationTruncateBasis(WaveFunctionPerturbationType * const a, const double epsilon);
void   WaveFunctionPerturbationOrthogonalizeConserveBasis(WaveFunctionPerturbationType * const psi,
                                                          const WaveFunctionPerturbationType * const OtherPsi,
                                                          const uint32_t Npsi);
uint32_t WaveFunctionPerturbationOrthogonalize(WaveFunctionPerturbationType * const psi,
                                               const WaveFunctionPerturbationType * const OtherPsi,
                                               const uint32_t Npsi);
uint32_t CreateStartWaveFunctionPerturbationOnsiteEnergy(const OperatorPerturbationType O,
                                                         WaveFunctionPerturbationType * psi, uint32_t order,
                                                         uint32_t nFermion, uint32_t nBoson, uint32_t NBitsKey);
uint32_t WaveFunctionPerturbationAddFermion(WaveFunctionPerturbationType const psiIn,
                                            WaveFunctionPerturbationType * const psiOut, uint32_t FermionPos,
                                            uint32_t occupied);
uint32_t WaveFunctionPerturbationRemoveFermion(WaveFunctionPerturbationType const psiIn,
                                               WaveFunctionPerturbationType * const psiOut, uint32_t FermionPos);


/**********************************************************/
/********   Core/RealWaveFunctionArithmetic.c      ********/
/**********************************************************/

void   RealWaveFunctionFree();
void   RealWaveFunctionFreeValues();
uint32_t RealWaveFunctionAddElement();
void   RealWaveFunctionAddElementConserveBasis();
double RealWaveFunctionPrefactorDet();
uint32_t RealWaveFunctionAddElementOMPMiniFlush();
void   RealWaveFunctionAddSameBasis();
void   RealWaveFunctionAddMultipliedSameBasis();
uint32_t RealWaveFunctionAddUnionBasis();
uint32_t RealWaveFunctionAddMultipliedUnionBasis();
void   RealWaveFunctionAddMultipliedConserveBasis();
void   RealWaveFunctionNormalize();
void   RealWaveFunctionMultiply();
double RealWaveFunctionCalculateNorm();
double RealWaveFunctionDotProductSameBasis();
double RealWaveFunctionDotProduct();
void   RealWaveFunctionRandomize();
void   RealWaveFunctionTruncateBasis();
void   RealWaveFunctionOrthogonalizeConserveBasis();
uint32_t RealWaveFunctionOrthogonalize();

/**********************************************************/
/********   Core/ComplexWaveFunctionArithmetic.c   ********/
/**********************************************************/

void   ComplexWaveFunctionFree();
void   ComplexWaveFunctionFreeValues();
uint32_t ComplexWaveFunctionAddElement();
void   ComplexWaveFunctionAddElementConserveBasis();
void   ComplexWaveFunctionPrefactorDet();
uint32_t ComplexWaveFunctionAddElementOMPMiniFlush();
void   ComplexWaveFunctionAddSameBasis();
void   ComplexWaveFunctionAddMultipliedSameBasis();
void   ComplexwaveFunctionAddComplexMultipliedSamebasis();
uint32_t ComplexWaveFunctionAddUnionBasis();
uint32_t ComplexWaveFunctionAddMultipliedUnionbasis();
uint32_t ComplexWaveFunctionAddComplexMultipliedUnionBasis();
void   ComplexWaveFunctionAddMultipliedSameBasis();
void   ComplexWaveFunctionNormalize();
void   ComplexWaveFunctionMultiply();
void   ComplexWaveFunctionComplexMultiply();
double ComplexWaveFunctionCalculateNorm();
void   ComplexWaveFunctionDotProductSameBasis();
void   ComplexWaveFunctionDotProduct();
void   ComplexWaveFunctionRandomize();
void   ComplexWaveFunctionTruncatebasis();
uint32_t ComplexWaveFunctionOrthogonalize();
uint32_t MakeComplexWaveFunctionFromReal();
uint32_t MakeRealWaveFunctionFromComplex();



/**********************************************************/
/*************       Core/Lanczos.c          **************/
/**********************************************************/

uint32_t LanczosGroundStateRR(const OperatorType H, WaveFunctionType * const psi, double *epsilon, const double zero,
                              double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi);
uint32_t LanczosGroundStateRC(const OperatorType H, WaveFunctionType * const psi, double *epsilon, const double zero,
                              double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi);
uint32_t LanczosGroundStateCR(const OperatorType H, WaveFunctionType * const psi, double *epsilon, const double zero,
                              double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi);
uint32_t LanczosGroundStateCC(const OperatorType H, WaveFunctionType * const psi, double *epsilon, const double zero,
                              double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi);
uint32_t LanczosGroundStateMR(const OperatorType H, WaveFunctionType * const psi, double *epsilon, const double zero,
                              double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi);
uint32_t LanczosGroundStateMC(const OperatorType H, WaveFunctionType * const psi, double *epsilon, const double zero,
                              double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi);
uint32_t LanczosTriDiagonalizeRR(const OperatorType H, WaveFunctionType psi, double epsilon,
                                 TriDiagonalMatrixType * const M, double Gamma);
uint32_t LanczosTriDiagonalizeCR(const OperatorType H, WaveFunctionType psi, double epsilon,
                                 TriDiagonalMatrixType * const M, double Gamma);
uint32_t LanczosTriDiagonalizeMR(const OperatorType H, WaveFunctionType psi, double epsilon,
                                 TriDiagonalMatrixType * const M, double Gamma);
uint32_t LanczosTriDiagonalizeRC(const OperatorType H, WaveFunctionType psi, double epsilon,
                                 TriDiagonalMatrixType * const M, double Gamma);
uint32_t LanczosTriDiagonalizeCC(const OperatorType H, WaveFunctionType psi, double epsilon,
                                 TriDiagonalMatrixType * const M, double Gamma);
uint32_t LanczosTriDiagonalizeMC(const OperatorType H, WaveFunctionType psi, double epsilon,
                                 TriDiagonalMatrixType * const M, double Gamma);
uint32_t LanczosTriDiagonalizePlusDensityRR(const OperatorType H, WaveFunctionType psi, double epsilon,
                                            TriDiagonalMatrixType * const M, MatrixType * const Rho, double Gamma);
uint32_t LanczosTriDiagonalizePlusDensityCR(const OperatorType H, WaveFunctionType psi, double epsilon,
                                            TriDiagonalMatrixType * const M, MatrixType * const Rho, double Gamma);
uint32_t LanczosTriDiagonalizePlusDensityMR(const OperatorType H, WaveFunctionType psi, double epsilon,
                                            TriDiagonalMatrixType * const M, MatrixType * const Rho, double Gamma);
uint32_t LanczosTriDiagonalizePlusDensityRC(const OperatorType H, WaveFunctionType psi, double epsilon,
                                            TriDiagonalMatrixType * const M, MatrixType * const Rho, double Gamma);
uint32_t LanczosTriDiagonalizePlusDensityCC(const OperatorType H, WaveFunctionType psi, double epsilon,
                                            TriDiagonalMatrixType * const M, MatrixType * const Rho, double Gamma);
uint32_t LanczosTriDiagonalizePlusDensityMC(const OperatorType H, WaveFunctionType psi, double epsilon,
                                            TriDiagonalMatrixType * const M, MatrixType * const Rho, double Gamma);


uint32_t LanczosGroundState(const OperatorType H, WaveFunctionType * const psi, double *epsilon, const double zero,
                            double shiftin, WaveFunctionType * OtherPsi, uint32_t NOtherPsi);
uint32_t LanczosTriDiagonalize(const OperatorType H, WaveFunctionType psi, double epsilon,
                               TriDiagonalMatrixType * const M, double Gamma);
uint32_t LanczosTriDiagonalizePlusDensity(const OperatorType H, WaveFunctionType psi, double epsilon,
                                          TriDiagonalMatrixType * const M, MatrixType * const Rho, double Gamma);



/**********************************************************/
/********       Core/OperatorArithmetic.c          ********/
/**********************************************************/

uint32_t FreeOperator();
uint32_t InitOperator();
uint32_t MakeRealOperatorFixedLengthComplex();
uint32_t InitOperatorFixedLength();
uint32_t InitOperatorCopySettings();
uint32_t RealOperatorAddLader();
uint32_t ComplexOperatorAddLader();
uint32_t RealOperatorAddLaderNormalOrder();
uint32_t ComplexOperatorAddLaderNormalOrder();
uint32_t OperatorAddLader(OperatorType * const O, const uint16_t * const oadd, const uint32_t length, const double VRe,
                          const double VIm);
uint32_t OperatorAddLaderNormalOrder(OperatorType * const O, const uint16_t * const oadd, const uint32_t length,
                                     const double VRe, const double VIm, uint32_t(*QDetNotOrdered) (uint16_t, uint16_t,
                                                                                                    uint32_t));
uint32_t OperatorRemoveSmallValues();
uint32_t QDetNotOrderedFaFcBaBc(const uint16_t A, const uint16_t B, const uint32_t NFermion);
uint32_t QDetNotOrderedFcFaBcBa(const uint16_t A, const uint16_t B, const uint32_t NFermion);
uint32_t QDetNotOrderedFcFaBaBc(const uint16_t A, const uint16_t B, const uint32_t NFermion);
uint32_t QDetNotOrderedFaFcBcBa(const uint16_t A, const uint16_t B, const uint32_t NFermion);
uint32_t NormalOrderOperator();
uint32_t OperatorMultiply();
uint32_t OperatorScale();
uint32_t OperatorAdd();
uint32_t MappedOperatorAddAtLevel();
uint32_t MappedOperatorAdd();

/**********************************************************/
/********       Core/OperatorWaveFunction.c        ********/
/**********************************************************/

double OperatorCr(uint8_t * const Determinant, const uint16_t i, const uint32_t NFermion);
double OperatorAn(uint8_t * const Determinant, const uint16_t i, const uint32_t NFermion);
double OperatorLader(uint8_t * const Determinant, const uint16_t * const o, const uint32_t lenght,
                     const uint32_t NFermion);
uint32_t OperatorPsiConserveBasisCC();
uint32_t OperatorPsiConserveBasisCR();
uint32_t OperatorPsiConserveBasisRC();
uint32_t OperatorPsiConserveBasisRR();
uint32_t OperatorPsiConserveBasisMC();
uint32_t OperatorPsiConserveBasisMR();
uint32_t OperatorPsiConserveBasisPreAllocatedCC();
uint32_t OperatorPsiConserveBasisPreAllocatedCR();
uint32_t OperatorPsiConserveBasisPreAllocatedRC();
uint32_t OperatorPsiConserveBasisPreAllocatedRR();
uint32_t OperatorPsiConserveBasisPreAllocatedMC();
uint32_t OperatorPsiConserveBasisPreAllocatedMR();
uint32_t OperatorPsiCC();
uint32_t OperatorPsiCR();
uint32_t OperatorPsiRC();
uint32_t OperatorPsiRR();
uint32_t OperatorPsiMC();
uint32_t OperatorPsiMR();
uint32_t OperatorPsiSerialCC();
uint32_t OperatorPsiSerialCR();
uint32_t OperatorPsiSerialRC();
uint32_t OperatorPsiSerialRR();
uint32_t OperatorPsiSerialMC();
uint32_t OperatorPsiSerialMR();
uint32_t PsiOperatorPsiCC();
uint32_t PsiOperatorPsiCR();
uint32_t PsiOperatorPsiRC();
uint32_t PsiOperatorPsiRR();
uint32_t PsiOperatorPsiMC();
uint32_t PsiOperatorPsiMR();
uint32_t OperatorPsiConserveBasis();
uint32_t OperatorPsiConserveBasisPreAllocated();
uint32_t OperatorPsi();
uint32_t OperatorPsiSerial();
uint32_t PsiOperatorPsi();

/**********************************************************/
/*****    Core/OperatorWaveFunctionPerturbation.c     *****/
/**********************************************************/

uint32_t OperatorPsiPerturbationConserveBasis();
uint32_t OperatorPsiPerturbationConserveBasisPreAllocated();
uint32_t OperatorPsiPerturbation();
uint32_t PsiOperatorPsiPerturbation();
uint32_t PsiHermitianOperatorPsiPerturbation();
uint32_t PsiOperatorSqrPsiPerturbation();
uint32_t PsiHermitianOperatorSqrPsiPerturbation();
uint32_t PsiPrimeOperatorPsiPerturbationMatrix();
uint32_t PsiOperatorPsiPerturbationMatrix();
uint32_t PsiPrimeOperatorPsiPerturbation();



/**********************************************************/
/**************        Core/Spectra.c        **************/
/**********************************************************/
void   CalculateMoment();
int SpectrumAllocate();
void SpectrumFree();



/**********************************************************/
/**************     Core/DensityMatrix.c     **************/
/**********************************************************/
uint32_t CalculateDensity();
uint32_t CalculatePrefactorHistogram();
uint32_t CalculatePartialDensityCorrelationTensor();



double InterpolateRootOfset();
double InterpolateRoot();
