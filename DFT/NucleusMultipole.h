/*
 * NucleusMultipole.h
 *
 * Non-spherical nuclear charge densities, "Function3D" (DFT stage A5c, design section 3.1),
 * see NucleusMultipole.cpp.
 */

#ifndef DFT_NUCLEUSMULTIPOLE_H_
#define DFT_NUCLEUSMULTIPOLE_H_

double RealSphericalHarmonic(int L, int M, double x, double phi);
uint32_t NucleusMultipoleCount(uint32_t Lmax);
void NucleusMultipoleLM(uint32_t i, int *L, int *M);
int NucleusMultipoleProject(const double *rho, uint32_t nP, uint32_t nTheta, uint32_t nPhi, const double *xTheta,
                            const double *wTheta, uint32_t Lmax, double *rhoLM);
int NucleusMultipolePotentials(const RadialBasisType *b, uint32_t Lmax, const double *rhoLM, double *VLM);
double NucleusMultipoleMoment(const RadialBasisType *b, int L, const double *rhoLM);
void RealHarmonicMatrixElement(int kappaA, int muA2, int L, int M, int kappaB, int muB2, double *re, double *im);

#endif /* DFT_NUCLEUSMULTIPOLE_H_ */
