uint32_t EvaluateBlockListOfPoles(BlockListOfPoles const *A, CompactMatrixType * Gw, const double omega, const double Gamma);
uint32_t EvaluateBlockListOfPoles(void const *A, CompactMatrixType * Gw, const double omega, const double Gamma);
uint32_t EvaluateBlockAndersonMatrix(BlockAndersonMatrixType const *A, CompactMatrixType * Gw, const double omega, const double Gamma);
uint32_t EvaluateBlockAndersonMatrix(void const *A, CompactMatrixType * Gw, const double omega, const double Gamma);
uint32_t EvaluateBlockTriDiagonalMatrix(BlockTriDiagonalMatrixType const * A, CompactMatrixType * Gw, const double omega, const double Gamma);
uint32_t EvaluateBlockTriDiagonalMatrix(void const * A, CompactMatrixType * Gw, const double omega, const double Gamma);
uint32_t EvaluateBlockDoubleTriDiagonalMatrix(BlockDoubleTriDiagonalMatrixType const * A, CompactMatrixType * Gw, const double omega, const double Gamma);
uint32_t EvaluateBlockDoubleTriDiagonalMatrix(void const * A, CompactMatrixType * Gw, const double omega, const double Gamma);
uint32_t EvaluateBlockNaturalImpurityOrbital(BlockNaturalImpurityOrbitalMatrixType const *A, CompactMatrixType * Gw, const double omega, const double Gamma);
uint32_t EvaluateBlockNaturalImpurityOrbital(void const *A, CompactMatrixType * Gw, const double omega, const double Gamma);
uint32_t EvaluateBlockResponsefunctionType(ResponsefunctionType const *A, CompactMatrixType * Gw, const double omega, const double Gamma);

uint32_t BlockListOfPolesToSpectra(double GroundStateEnergy, BlockListOfPoles const *A, SpectraType *G, double Emin, double Emax, uint32_t N);
uint32_t BlockListOfPolesToSpectraWithBroadening(double GroundStateEnergy, BlockListOfPoles const *A, SpectraType *G, double Emin, double Emax, uint32_t N, double Gamma);
uint32_t BlockListOfPolesEigenValues(const BlockListOfPoles * const A, double * val);
uint32_t BlockListOfPolesEigenSystem(const BlockListOfPoles * const A, double * val, CompactMatrixType * fun);
uint32_t BlockListOfPolesToBlockTridiagonalMatrix(BlockListOfPoles const A, BlockTriDiagonalMatrixType *const B, int const NBlockMax, double const SingularValue);
uint32_t BlockListOfPolesToBlockTridiagonalMatrixReduceDimension(BlockListOfPoles const A, BlockTriDiagonalMatrixType *B, CompactMatrixType const *AT,  int const NBlockMax, double const SingularValue);
uint32_t BlockListOfPolesToBlockDoubleTriDiagonalMatrix(BlockListOfPoles const *A, BlockDoubleTriDiagonalMatrixType * const B, double const kbT, double const SingularValue, double const MinimalPoleDistance);
uint32_t BlockListOfPolesToBlockNaturalImpurityOrbitalMatrix(BlockListOfPoles const *A, BlockNaturalImpurityOrbitalMatrixType * const B, double const kbT, uint32_t const NBlockMaxVal, uint32_t const NBlockMaxCon);
uint32_t BlockListOfPolesToBlockNaturalImpurityOrbitalMatrixReduceDimension(BlockListOfPoles const *A, BlockNaturalImpurityOrbitalMatrixType * const B, CompactMatrixType const *AT, double const kbT, uint32_t const NBlockMaxVal, uint32_t const NBlockMaxCon);
uint32_t BlockListOfPolesToBlockAndersonMatrix(BlockListOfPoles const A, BlockAndersonMatrixType *const B, int const NBlockMax, double const SingularValue);
uint32_t BlockListOfPolesToBlockAndersonMatrixReduceDimension(BlockListOfPoles const A, BlockAndersonMatrixType *B, CompactMatrixType const *AT, int const NBlockMax);
uint32_t BlockListOfPolesToDensityMatrix(BlockListOfPoles const * const A, CompactMatrixType * const rho, double const mu, double const kbT);

uint32_t BlockAndersonMatrixToSpectra(const double Energy, const BlockAndersonMatrixType *const A, SpectraType *G, double Emin, double Emax, uint32_t N);
uint32_t BlockAndersonMatrixToSpectraWithBroadening(const double Energy, const BlockAndersonMatrixType *const A, SpectraType *G, double Emin, double Emax, uint32_t N, double const Gamma);
uint32_t BlockAndersonMatrixEigenValues(const BlockAndersonMatrixType * const A, double * val);
uint32_t BlockAndersonMatrixEigenSystem(const BlockAndersonMatrixType * const A, double * val, CompactMatrixType * fun);
uint32_t BlockAndersonMatrixToBlockListOfPoles(BlockAndersonMatrixType const * A, BlockListOfPoles * const G, double const MinimalPoleDistance, bool relativeMinimumPoleDistance);
uint32_t BlockAndersonMatrixToBlockListOfPolesVariableGrid(BlockAndersonMatrixType const AAnd, BlockListOfPoles * const B, const double * const Grid, const uint32_t NGridPoints);
uint32_t BlockAndersonMatrixToBlockTridiagonalMatrix(BlockAndersonMatrixType const * A, BlockTriDiagonalMatrixType * const B, const double SingularValue, double const MinimalPoleDistance, bool const relativeDistance);
uint32_t BlockAndersonMatrixInverseToBlockListOfPoles(BlockAndersonMatrixType const * A, BlockListOfPoles * const Delta);
uint32_t BlockAndersonMatrixToBlockDoubleTriDiagonalMatrix(BlockAndersonMatrixType const * A, BlockDoubleTriDiagonalMatrixType * const B, const double kbT, double const SingularValue, double const MinimalPoleDistance);
uint32_t BlockAndersonMatrixToBlockNaturalImpurityOrbitalMatrix(BlockAndersonMatrixType const * A, BlockNaturalImpurityOrbitalMatrixType * const B, const double kbT, double const MinimalPoleDistance, bool const relativeDistance);

uint32_t BlockTriDiagonalMatrixToSpectra(const double Energy, BlockTriDiagonalMatrixType const * A, SpectraType *const G, double const Emin, double const Emax, uint32_t const N);
uint32_t BlockTriDiagonalMatrixToSpectraWithBroadening(const double Energy, BlockTriDiagonalMatrixType const * A, SpectraType *const G, double const Emin, double const Emax, uint32_t const N, double const Gamma);
uint32_t BlockTriDiagonalMatrixEigenValues(const BlockTriDiagonalMatrixType * const A, double * val);
uint32_t BlockTriDiagonalMatrixEigenSystem(const BlockTriDiagonalMatrixType * const A, double * val, CompactMatrixType * fun);
uint32_t BlockTriDiagonalMatrixEigenSystemWithDeflation(const BlockTriDiagonalMatrixType * const A, double * val, CompactMatrixType * fun, const uint32_t * const BlockSize);
uint32_t BlockTriDiagonalMatrixToBlockListOfPoles(BlockTriDiagonalMatrixType const * A, BlockListOfPoles * const G, double const MinPoleDistance, bool const relativeDistance);
uint32_t BlockTriDiagonalMatrixToBlockListOfPolesVariableGrid(BlockTriDiagonalMatrixType const * A, BlockListOfPoles * const G, const double * const Grid, const uint32_t NGridPoints);
uint32_t BlockTriDiagonalMatrixToBlockListOfPolesFixedGrid(BlockTriDiagonalMatrixType const * A, BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE);
uint32_t BlockTriDiagonalMatrixToBlockAndersonMatrix(BlockTriDiagonalMatrixType const ATri, BlockAndersonMatrixType * const B);
uint32_t BlockTriDiagonalMatrixToBlockDoubleTriDiagonalMatrix(BlockTriDiagonalMatrixType const ATri, BlockDoubleTriDiagonalMatrixType * const B, const double kBT, double const SingularValue, double const MinimalPoleDistance);
uint32_t BlockTriDiagonalMatrixToBlockNaturalImpurityMatrix(BlockTriDiagonalMatrixType const ATri, BlockNaturalImpurityOrbitalMatrixType * const B, const double kBT, double const MinimalPoleDistance, bool const relativeDistance);

uint32_t BlockDoubleTriDiagonalMatrixToSpectra(const double Energy, BlockDoubleTriDiagonalMatrixType const * A, SpectraType *const G, double const Emin, double const Emax, uint32_t const N);
uint32_t BlockDoubleTriDiagonalMatrixToSpectraWithBroadening(const double Energy, BlockDoubleTriDiagonalMatrixType const * A, SpectraType *const G, double const Emin, double const Emax, uint32_t const N, double const Gamma);
uint32_t BlockDoubleTriDiagonalMatrixToBlockTriDiagonalMatrix(BlockDoubleTriDiagonalMatrixType const * M, BlockTriDiagonalMatrixType * const H, double const SingularValue, double const MinimalPoleDistance, bool const relativeDistance);
uint32_t BlockDoubleTriDiagonalMatrixToBlockAndersonMatrix(BlockDoubleTriDiagonalMatrixType const * M, BlockAndersonMatrixType * const H, double const MinimalPoleDistance, bool const relativeDistance);
uint32_t BlockDoubleTriDiagonalMatrixToBlockNaturalImpurityMatrix(BlockDoubleTriDiagonalMatrixType const M, BlockNaturalImpurityOrbitalMatrixType * const H);
uint32_t BlockDoubleTriDiagonalMatrixToBlockListOfPoles(BlockDoubleTriDiagonalMatrixType const M, BlockListOfPoles * const H);
uint32_t BlockDoubleTriDiagonalMatrixToCompactMatrix(BlockDoubleTriDiagonalMatrixType const M, CompactMatrixType * const H);

uint32_t BlockNaturalImpurityOrbitalMatrixToSpectra(BlockNaturalImpurityOrbitalMatrixType const * A, SpectraType *const G, double const Emin, double const Emax, uint32_t const N);
uint32_t BlockNaturalImpurityOrbitalMatrixToSpectraWithBroadening(BlockNaturalImpurityOrbitalMatrixType const * A, SpectraType *const G, double const Emin, double const Emax, uint32_t const N, double const Gamma);
uint32_t BlockNaturalImpurityOrbitalMatrixEigenValues(const BlockNaturalImpurityOrbitalMatrixType * const A, double * val);
uint32_t BlockNaturalImpurityOrbitalMatrixEigenSystem(const BlockNaturalImpurityOrbitalMatrixType * const A, double * val, CompactMatrixType * fun);
uint32_t BlockNaturalImpurityOrbitalMatrixToBlockAndersonMatrix(BlockNaturalImpurityOrbitalMatrixType const * A, BlockAndersonMatrixType * const B);
uint32_t BlockNaturalImpurityOrbitalMatrixToBlockListOfPoles(BlockNaturalImpurityOrbitalMatrixType const * A, BlockListOfPoles * const B, double const MinimalPoleDistance, bool const relativeDistance);
uint32_t BlockNaturalImpurityOrbitalMatrixToBlockListOfPolesVariableGrid(BlockNaturalImpurityOrbitalMatrixType const AAnd, BlockListOfPoles * const B, const double * const Grid, const uint32_t NGridPoints);
uint32_t BlockNaturalImpurityOrbitalMatrixToBlockTridiagonalMatrix(BlockNaturalImpurityOrbitalMatrixType const * A, BlockTriDiagonalMatrixType * const B, const double SingularValue);
uint32_t BlockNaturalImpurityOrbitalMatrixToBlockDoubleTridiagonalMatrix(BlockNaturalImpurityOrbitalMatrixType const * A, BlockDoubleTriDiagonalMatrixType * const B, const double kBT, const double SingularValue, double const MinimalPoleDistance);
uint32_t BlockListOfPolesInversedToCompactMatrix(BlockListOfPoles const * const A, CompactMatrixType * const B);

uint32_t TightBindingHamiltonianToDensityMatrix(TightBindingDef const * const A, int Nk[3], CompactMatrixType *rho, double mu, double kbT);
uint32_t TightBindingHamiltonianToBlockListOfPoles(TightBindingDef const A, int Nk[3], BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE);
uint32_t TightBindingHamiltonianToBlockListOfPolesGrid(TightBindingDef const A, int Nk[3], BlockListOfPoles * const G, double const * const Egrid, uint32_t const NEgrid);
uint32_t TightBindingHamiltonianMultiSiteToBlockListOfPoles(TightBindingDef const A, int Nk[3], uint32_t const Nvec, int const *vec, BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE);
uint32_t TightBindingHamiltonianPlusSelfEnergyToBlockListOfPoles(TightBindingDef const A, int Nk[3], BlockListOfPoles const Sigma, BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE);
uint32_t TightBindingHamiltonianPlusSelfEnergyToBlockListOfPolesGrid(TightBindingDef const A, int Nk[3], BlockListOfPoles const Sigma, BlockListOfPoles * const G, double const * const Egrid, uint32_t const NEgrid);
uint32_t TightBindingHamiltonianPlusSelfEnergyMultiSiteToBlockListOfPoles(TightBindingDef const A, int Nk[3], uint32_t const Nvec, int const *vec, BlockListOfPoles const Sigma, BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE);
uint32_t TightBindingHamiltonianAddSpinPlusSelfEnergyToBlockListOfPoles(TightBindingDef const A, int Nk[3], BlockListOfPoles const Sigma, BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE);
uint32_t TightBindingHamiltonianAddSpinPlusSelfEnergyMultiSiteToBlockListOfPoles(TightBindingDef const A, int Nk[3], uint32_t const Nvec, int const *vec, BlockListOfPoles const Sigma, BlockListOfPoles * const G, double const Emin, double const Emax, uint32_t const NE);

uint32_t DiagonalMatrixBlockTridiagonalize(double const *A, CompactMatrixType *psi, CompactMatrixType *const T, BlockTriDiagonalMatrixType *const ABlock, double NBlock, double SingularValue);
