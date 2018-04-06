/* This file was automatically generated.  Do not edit! */
void AllocFailed();
uint32_t CalculatePrefactorHistogram(double *const Hist,const WaveFunctionType psi,const double exp,const uint32_t N);
uint32_t CalculateOccupationHistogram(double *const Hist,const WaveFunctionType psi,const uint16_t *fermionlist,const uint32_t N);
uint32_t CalculateDensity(double *const Rho,const WaveFunctionType psi);
uint32_t CalculatePartialDensityMatrix(MatrixType *const Rho,const WaveFunctionType psi,const uint32_t Nparticles,const uint16_t *const jindex);
uint32_t CalculateSpinCorrelationMatrix(MatrixType *const Rho,const WaveFunctionType psi);
uint32_t FreeOperator(OperatorType *const H);
uint32_t CalculatePartialDensityCorrelationTensor(Complex16 *Rho,const WaveFunctionType psi,uint32_t const NFermionIncluded,uint16_t const *FermionIncluded);
uint32_t CalculateDensityCorrelationMatrix(MatrixType *const Rho,const WaveFunctionType psi);
uint32_t CalculateDensityMatrixWithOffDiagonals(MatrixType *const Rho,uint32_t const NPsi,const WaveFunctionType *const psi);
uint32_t PsiOperatorPsiMatrix(const OperatorType O,const WaveFunctionType *psi,MatrixType *const V);
uint32_t CalculatePartialDensityMatrixWithOffDiagonals(MatrixType *const Rho,uint32_t const NPsi,const WaveFunctionType *const psi,uint32_t const NFermionIncluded,uint16_t const *FermionIncluded);
uint32_t PsiOperatorPsiRR(const OperatorType O,const WaveFunctionType psi,double *VRe);
uint32_t PsiOperatorPsiRC(const OperatorType O,const WaveFunctionType psi,double *VRe,double *VIm);
uint32_t CalculateDensityMatrix(MatrixType *const Rho,const WaveFunctionType psi);
